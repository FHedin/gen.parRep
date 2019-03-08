/**
 * \file omm_interface.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#include <cstdio>
#include <cstring>
#include <cmath>

#include <sstream>
#include <iostream>
#include <fstream>
#include <memory>

#include "logger.hpp"
#include "rand.hpp"
#include "mpi_utils.hpp"
#include "lua_interface.hpp"
#include "omm_interface.hpp"

// A few useful constants definition taken from OMM's reference/SimTKOpenMMRealType.h
// #define ANGSTROM     (1e-10) 
#define KILO         (1e3)
// #define NANO         (1e-9)
// #define PICO         (1e-12)
// #define A2NM         (ANGSTROM/NANO)
// #define NM2A         (NANO/ANGSTROM)
// #define RAD2DEG      (180.0/M_PI)
// #define CAL2JOULE    (4.184)
// #define E_CHARGE     (1.60217733e-19)

// #define AMU          (1.6605402e-27)
#define BOLTZMANN    (1.380658e-23)            /* (J/K)   */
#define AVOGADRO     (6.0221367e23)          
#define RGAS         (BOLTZMANN*AVOGADRO)      /* (J/(mol K))  */
#define BOLTZ        (RGAS/KILO)               /* (kJ/(mol K)) */
// #define FARADAY      (E_CHARGE*AVOGADRO)       /* (C/mol)      */
// #define ELECTRONVOLT (E_CHARGE*AVOGADRO/KILO)  /* (kJ/mol)   */     
// 
// #define EPSILON0     (5.72765E-4)              /* (e^2 Na/(kJ nm)) == (e^2/(kJ mol nm)) */ 
// 
// #define SPEED_OF_LIGHT   (2.9979245800E05)      /* nm/ps                */
// #define ATOMICMASS_keV   (940000.0)             /* Atomic mass in keV   */
// #define ELECTRONMASS_keV (512.0)                /* Electron mas in keV  */
// 
// #define ONE_4PI_EPS0      138.935456
// #define PRESFAC           (16.6054)             /* bar / pressure unity */
// #define ENM2DEBYE         48.0321               /* Convert electron nm to debye */
// #define DEBYE2ENM         0.02081941

/**
 * constant used when calculating temperature from kinetic energy
 * T = (2./(dofs*K_Boltzmann))*Ekin
 * where dofs in the number of degrees of freedom in the system
 * So EkinToK = 2./K_Boltzmann with K_Boltzmann in units of kJ/(mol.K)
 */
const double OMM_interface::EkinToK = 2.0/BOLTZ;

using OpenMM::Vec3;
using namespace std;

const std::map<PLATFORMS,std::string> OMM_interface::omm_platforms_names = 
{
  {PLATFORMS::REF, "Reference"},
  {PLATFORMS::CPU, "CPU"},
  {PLATFORMS::OCL, "OpenCL"},
  {PLATFORMS::CUDA,"CUDA"}
};

//------------------------------------------------------------------------
/**
 * C++ code interacting with the OpenMM framework
 */
//------------------------------------------------------------------------

/**
 * Constructor : build system from serialization instead
 * then most of the calls to other methods below are not useful anymore
 */
OMM_interface::OMM_interface(DATA& _dat,
                             const string& omm_name_version,
                             OpenMM::System& syst,
                             OpenMM::Integrator& integr,
                             OpenMM::State& state,
                             PLATFORMS platformDesired,
                             const PLATFORM_PROPERTIES& platformProperties
) : MD_interface(_dat,
                 MD_ENGINES::OPENMM,
                 omm_name_version,
                 MD_DISTANCE_UNIT::NANOMETER,
                 MD_ENERGY_UNIT::KJ_PER_MOL,
                 true
                )
/* MD_interface::MD_interface(DATA& _dat,MD_ENGINES _engine_type, string& _engine_description, MD_DISTANCE_UNIT _distance_unit, MD_ENERGY_UNIT _energy_unit, bool _engine_supports_groups_splitting */
{
  // Load all available OpenMM plugins from their default location.
  const string& plugins_dir = OpenMM::Platform::getDefaultPluginsDirectory();
  OpenMM::Platform::loadPluginsFromDirectory(plugins_dir);
  
  LOG_PRINT(LOG_INFO,"OMM plugins successfully loaded from directory %s\n",plugins_dir.c_str());
  
  // deserialize the system
  system = unique_ptr<OpenMM::System>(&syst);
  
  // deserialize the integrator
  integrator = unique_ptr<OpenMM::Integrator>(&integr);
  
  // re-seed the integrator
  if(dynamic_cast<OpenMM::LangevinIntegrator*>(integrator.get()))
  {
    lint = dynamic_cast<OpenMM::LangevinIntegrator*>(integrator.get());
    integType = INTEGRATORS::LANGEVIN;
    int32_t seed = get_int32();
    lint->setRandomNumberSeed(seed);
    LOG_PRINT(LOG_INFO,"OMM unserialised Langevin integrator re-initialised with seed %d\n",seed);
  }
  else if (dynamic_cast<OpenMM::BrownianIntegrator*>(integrator.get()))
  {
    bint = dynamic_cast<OpenMM::BrownianIntegrator*>(integrator.get());
    integType = INTEGRATORS::BROWNIAN;
    int32_t seed = get_int32();
    bint->setRandomNumberSeed(seed);
    LOG_PRINT(LOG_INFO,"OMM unserialised Brownian integrator re-initialised with seed %d\n",seed);
  }
  else
  {
    throw runtime_error("Error when re-seeding the OMM integrator : integrator type is neither OpenMM::LangevinIntegrator or OpenMM::BrownianIntegrator !\n");
  }
  
  // from system and integrator, create a platform and a context using the state
  addPlatform(platformDesired,platformProperties,state);
  
  // retireve pbc data and store it in dat
  state.getPeriodicBoxVectors(pbc[0],pbc[1],pbc[2]);
  memcpy(dat.pbc.data(),pbc.data(),sizeof(pbc));
  
  // fill dat structure with unserialised data
  dat.natom    = system->getNumParticles();
  dat.timestep = integrator->getStepSize();

  // retrieve initial positions and velocities
  pos = vector<Vec3>(dat.natom);
  vel = vector<Vec3>(dat.natom);
  
  /*
   * calculate the number of degrees of freedom of the system
   * it is 3 times the number of particles (not virtual ones for which mass=0)
   * minus the number of constraints
   * and minus 3 if a center of mass removal is enabled
   */
  dofs = 0.0;
  
  for(int32_t i=0; i<system->getNumParticles(); i++)
  {
    if(system->getParticleMass(i) > 1e-8)
      dofs += 3;
  }
  
  dofs -= system->getNumConstraints();
  
  for(int32_t i=0; i<system->getNumForces(); i++)
  {
    OpenMM::Force& f = system->getForce(i);
    if(dynamic_cast<OpenMM::CMMotionRemover*>(&f))
    {
      dofs -=3;
      break;
    }
  }
  
}

OMM_interface::~OMM_interface()
{
}

void OMM_interface::addPlatform(PLATFORMS platformDesired, const PLATFORM_PROPERTIES& properties, OpenMM::State& state)
{
  OpenMM::Platform* platform = &(OpenMM::Platform::getPlatformByName("Reference"));
  
  // try to use the best available platform
  if(platformDesired == AUTO)
  {
    // count number of platforms
    int nplatforms = OpenMM::Platform::getNumPlatforms() ;
    fprintf(stdout,"Number of OpenMM platforms detected : %d\n",nplatforms);
    
    // and select platform based on the speed attribute
    double best_speed=-1.0;
    int best_platform=-1;
    for(int index=0; index<nplatforms; index++)
    {
      OpenMM::Platform& lplatform = OpenMM::Platform::getPlatform(index);
      double lspeed = lplatform.getSpeed();
      fprintf(stdout," Platform[%d] is : %s | speed is %lf \n",index,lplatform.getName().c_str(),lspeed);
      best_platform = (lspeed > best_speed) ?  index:best_platform;
      best_speed    = (lspeed > best_speed) ? lspeed:best_speed;
    }
        
    platform = &(OpenMM::Platform::getPlatform(best_platform));
    fprintf(stdout,"Will use Platform '%s', which is apparently the fastest\n",platform->getName().c_str());
    
  }
  // or if the user specified one, try to initialise it directly
  else
  {
    const string wantPlatform = omm_platforms_names.at(platformDesired);
    platform = &(OpenMM::Platform::getPlatformByName(wantPlatform));
    fprintf(stdout,"Will use Platform '%s', provided by user.\n",wantPlatform.c_str());
    LOG_PRINT(LOG_INFO,"Will use Platform '%s', provided by user.\n",wantPlatform.c_str());
  }
        
  // create the context from all the objects
  context = unique_ptr<OpenMM::Context> (new OpenMM::Context(*system, *integrator, *platform, properties) );
    
  // also register initial state
  context->setState(state);
    
  platformName = context->getPlatform().getName();
  
  LOG_PRINT(LOG_INFO,"OMM context created for platform %s\n",platformName.c_str());
}

void OMM_interface::doNsteps(uint32_t numSteps)
{
  integrator->step((int32_t)numSteps);
}

void OMM_interface::minimise(double tol, uint32_t maxsteps, ENERGIES& ener)
{
  // minimisation with L-BGFS
  OpenMM::LocalEnergyMinimizer::minimize(*context,tol,(int32_t)maxsteps);
  
  // get back data from context and extract energies
  const int infoMask = OpenMM::State::Energy;
  const OpenMM::State state = context->getState(infoMask,true);
  
  ener.epot() = state.getPotentialEnergy();
  ener.ekin() = state.getKineticEnergy();
  ener.etot() = ener.epot() + ener.ekin();
}

uint32_t OMM_interface::runForPhysicalTime(chrono::milliseconds t)
{
  uint32_t steps = 0;
  
  chrono::nanoseconds lt = chrono::duration_cast<std::chrono::nanoseconds>(t);
  
  if(first_call_of_runForPhysicalTime)
  {
    auto start  = chrono::high_resolution_clock::now();
    integrator->step(1);
    auto end    = chrono::high_resolution_clock::now();
    
    cost_for_1step_dyna = chrono::duration_cast<std::chrono::nanoseconds>(end-start);
    int64_t nsteps = (int64_t) lt.count() / (int64_t) cost_for_1step_dyna.count();
    
    if(nsteps != 0)
      integrator->step((int32_t)nsteps - 1);
    
    steps = (uint32_t) nsteps;
    
    first_call_of_runForPhysicalTime = false;
  }
  else
  {
    int64_t nsteps = (int64_t) lt.count() / (int64_t) cost_for_1step_dyna.count();
    
    if(nsteps != 0)
      integrator->step((int32_t)nsteps);
    
    steps = (uint32_t) nsteps;
  }

  return steps;
}

void OMM_interface::minimiseWithCopy(double tol, uint32_t maxsteps, ENERGIES& ener)
{
  // backup coordinates and velocities
  int infoMask = OpenMM::State::Positions | OpenMM::State::Velocities;
  
  OpenMM::State state = context->getState(infoMask,true);
  
  const vector<Vec3> lpos(state.getPositions());
  const vector<Vec3> lvel(state.getVelocities());
  
  //do minimisation as before and save energies
  minimise(tol,maxsteps,ener);
  
  // restore coordinates
  context->setPositions(lpos);
  context->setVelocities(lvel);
  
  context->applyConstraints(integrator->getConstraintTolerance());
  context->applyVelocityConstraints(integrator->getConstraintTolerance());
}

void OMM_interface::minimiseWithCopy(double tol, uint32_t maxsteps, ENERGIES& ener,
                                      vector<XYZ>& pos_out, vector<XYZ>& vels_out)
{
  // backup coordinates and velocities
  int infoMask = OpenMM::State::Positions | OpenMM::State::Velocities;
  
  OpenMM::State state = context->getState(infoMask,true);
  
  const vector<Vec3> lpos(state.getPositions());
  const vector<Vec3> lvel(state.getVelocities());

  //do minimisation as before and save energies
  minimise(tol,maxsteps,ener);

  // will return minimised coordinates and velocities
  state = context->getState(infoMask,true);
  
  const vector<Vec3> min_pos(state.getPositions());
  const vector<Vec3> min_vel(state.getVelocities());
  
  pos_out.resize(min_pos.size());
  vels_out.resize(min_vel.size());
  
  memcpy(&pos_out[0],&min_pos[0],min_pos.size()*sizeof(Vec3));
  memcpy(&vels_out[0],&min_vel[0],min_vel.size()*sizeof(Vec3));
  
  // restore coordinates
  context->setPositions(lpos);
  context->setVelocities(lvel);
  
  context->applyConstraints(integrator->getConstraintTolerance());
  context->applyVelocityConstraints(integrator->getConstraintTolerance());
  
}

void OMM_interface::setSimClockTime(double t)
{
  context->setTime(t);
}

void OMM_interface::setCrdsVels(ATOM at[])
{

  for (uint32_t n=0; n<dat.natom; n++)
  {
    memcpy(&pos[n],&at[n].crds,sizeof(Vec3));
    memcpy(&vel[n],&at[n].vels,sizeof(Vec3));
  }
  
  memcpy(pbc.data(),dat.pbc.data(),sizeof(pbc));
  
  context->setPeriodicBoxVectors(pbc[0],pbc[1],pbc[2]);
  
  context->setPositions(pos);
  context->setVelocities(vel);
  
  context->applyConstraints(integrator->getConstraintTolerance());
  context->applyVelocityConstraints(integrator->getConstraintTolerance());
  
}

void OMM_interface::randomiseVelocities()
{
  // set velocities to initial temperature
  int32_t velSeed = get_int32();
  context->setVelocitiesToTemperature(dat.T,velSeed);
  LOG_PRINT(LOG_INFO,"Velocities randomised to target T = %.3lf K with seed %d\n",dat.T,velSeed);
  
  double newT = 0.0;
  this->getSimData(nullptr,nullptr,&newT,nullptr);
  LOG_PRINT(LOG_INFO,"Real temperature after randomisation is : T = %.3lf K\n",newT);
}

void OMM_interface::getSimData(double* timeInPs, ENERGIES* energies, double* currentTemperature,
                               ATOM atoms[], int32_t groups)
{
  const bool wantTime        = !(timeInPs==nullptr);
  const bool wantEnergy      = !(energies==nullptr);
  const bool wantTemperature = !(currentTemperature==nullptr);
  const bool wantCrdVels     = !(atoms==nullptr);
  
  int32_t infoMask = 0;
  
  /**
   * OpenMM allows retrieving data for only a part of the System if it was splitted in several forceGroups.
   * In this case it requires a specific mask for deciding which groups to access.
   * 
   * Let us assume the Forces were divided in two groups: 0 and 1 (OMM supports groups 0 to 31)
   * 
   * If the result of the operation ((mask & (1<<i)) != 0) is true, the group i will have its data retrieved.
   *
   * For group 0, a valid mask is 0x1 i.e. 1: indeed (1<<0) = 2^0 = 1 , so 1 & 1 = 1 so group 0 is included.
   *                                          However as (1<<1) = 2^1 = 2  then 1 & 2 = 0 so group 1 is not included
   * 
   * Similarly for obtaining data for group 1 only : 
   * We use mask 0x2 (i.e. 2) because : 2 & (1<<1) = 2&2 = 2 , but group 0 not included because 2 & (1<<0) = 2 & 1 = 0
   * 
   * A mask of 0xFFFFFFFF makes sure to retrieve data for all groups
   * 
   * We assume that groups already contains a valid mask if it is not equal to -1
   */
  int32_t groupsMask = (groups==-1) ? 0xFFFFFFFF : groups;
  
  if(wantCrdVels)
    infoMask |= (OpenMM::State::Positions | OpenMM::State::Velocities);
  
  if(wantEnergy || wantTemperature)
    infoMask |= OpenMM::State::Energy;
  
  OpenMM::State state = context->getState(infoMask,true,groupsMask);
  
  if(wantEnergy || wantTemperature)
    T = EkinToK*state.getKineticEnergy()/dofs;
  
  if(wantTime)
    *timeInPs = state.getTime();
  
  if(wantCrdVels)
  {
    const vector<Vec3>& lpos = state.getPositions();
    const vector<Vec3>& lvel = state.getVelocities();
    
    state.getPeriodicBoxVectors(pbc[0],pbc[1],pbc[2]);
      
    // copy to C++ struct
    for (uint32_t n=0; n < dat.natom; n++)
    {
      memcpy(&atoms[n].crds,&lpos[n],sizeof(Vec3));
      memcpy(&atoms[n].vels,&lvel[n],sizeof(Vec3));
    }
    
    memcpy(dat.pbc.data(),pbc.data(),sizeof(pbc));
  }
          
  /* If energy has been requested, obtain it from state */
  if(wantEnergy)
  {
    energies->epot() = state.getPotentialEnergy();
    energies->ekin() = state.getKineticEnergy();
    energies->etot() = energies->epot() + energies->ekin();
  }
    
  if(wantTemperature)
    *currentTemperature = T;

}

void OMM_interface::getParticlesParams(ATOM atoms[])
{
  int32_t np = system->getNumParticles();

  for(int32_t n=0; n<np; n++)
    atoms[n].mass = system->getParticleMass(n);

}

double OMM_interface::getTime()
{
  const OpenMM::State state = context->getState(0);
  time = state.getTime();
  return time;
}

double OMM_interface::getTemperature()
{
  return T;
}

void OMM_interface::backupOMMobject()
{
  string fname = "omm.restart.chk";

  MPIutils::mpi_get_unique_name(fname);
  
  ofstream ofs(fname,ofstream::out|ofstream::binary);
  
  context->createCheckpoint(ofs);
  
  ofs.close();
  
}

void OMM_interface::restoreOMMobject()
{
  string fname = "omm.restart.chk";
  
  MPIutils::mpi_get_unique_name(fname);
  
  ifstream ifs(fname,ofstream::in|ofstream::binary);
  
  context->loadCheckpoint(ifs);
  
  ifs.close();
  
}

//------------------------------------------------------------------------
/**
 * Static methods
 */
//------------------------------------------------------------------------

unique_ptr<OMM_interface> OMM_interface::initialise_fromSerialization(DATA& _dat,
                                                                      unique_ptr<luaInterface>& luaItf
                                                                     )
{
  const lua_ParVal_map& luaParams = luaItf->get_parsed_parameters_map();
  
  // get paths to OMM XML serialised files
  const string sysXMLfile        = luaParams.at("system");
  const string integratorXMLfile = luaParams.at("integrator");
  const string stateXMLfile      = luaParams.at("state");
  // get desired OpenMM platform
  const string platformName      = luaParams.at("platform");
  
  LOG_PRINT(LOG_INFO,"Initialising OpenMM interface with serialized xml files : '%s' '%s' '%s'\n",
          sysXMLfile.c_str(),integratorXMLfile.c_str(),stateXMLfile.c_str());
  
  // --------------------------------------------- 
  
  OpenMM::System* sys = nullptr;
  ifstream file(sysXMLfile,ios_base::in);
  if(file.is_open())
  {
    sys = OpenMM::XmlSerializer::deserialize<OpenMM::System>(file);
    file.close();
    LOG_PRINT(LOG_INFO,"Serialized xml system file '%s' properly read\n",sysXMLfile.c_str());
  }
  else
  {
    throw runtime_error("Error when attempting to parse the OpenMM XML file " + sysXMLfile + " !\n");
  }

  // --------------------------------------------- 
  
  OpenMM::Integrator* integ = nullptr;
  file.open(integratorXMLfile,ios_base::in);
  if(file.is_open())
  {
    integ = OpenMM::XmlSerializer::deserialize<OpenMM::Integrator>(file);
    file.close();
    LOG_PRINT(LOG_INFO,"Serialized xml integrator file '%s' properly read\n",integratorXMLfile.c_str());
  }
  else
  {
    throw runtime_error("Error when attempting to parse the OpenMM XML file " + integratorXMLfile + " !\n");
  }
  
  // --------------------------------------------- 
  OpenMM::State st;
  file.open(stateXMLfile,ios_base::in);
  if(file.is_open())
  {
    st = *OpenMM::XmlSerializer::deserialize<OpenMM::State>(file);
    file.close();
    LOG_PRINT(LOG_INFO,"Serialized xml state file '%s' properly read\n",stateXMLfile.c_str());
  }
  else
  {
    throw runtime_error("Error when attempting to parse the OpenMM XML file " + stateXMLfile + " !\n");
  }

  // --------------------------------------------- 
  PLATFORMS platformDesired = AUTO;
  if (platformName=="AUTO")
    platformDesired = AUTO;
  else
  {
    for(auto& it : omm_platforms_names)
    {
      if (platformName == it.second)
      {
        platformDesired = it.first;
        break;
      }
    }
  }

  LOG_PRINT(LOG_INFO,"OMM requested platform is : %s\n",platformName.c_str());
  
  const PLATFORM_PROPERTIES& platformProperties = luaItf->get_omm_platform_properties();
  
  // create the c++ object corresponding to a c++ OpenMM interface
  const string omm_version = "Openmm" + OpenMM::Platform::getOpenMMVersion();
  OMM_interface* omm_itf = nullptr;
  omm_itf = new OMM_interface(_dat,omm_version,*sys,*integ,st,platformDesired,platformProperties);
  
  LOG_PRINT(LOG_INFO,"OpenMM interface properly initialised from serialised files !\n");
  
  return unique_ptr<OMM_interface>(omm_itf);
}
