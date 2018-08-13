/**
 * \file omm_interface.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
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

#include "omm_interface.hpp"
#include "mpi_utils.hpp"

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
                             OpenMM::System& syst,
                             OpenMM::Integrator& integr,
                             OpenMM::State& state,
                             PLATFORMS platformDesired
                             ) : MD_interface(_dat)
{
  // Load all available OpenMM plugins from their default location.
  const string& plugins_dir = OpenMM::Platform::getDefaultPluginsDirectory();
  OpenMM::Platform::loadPluginsFromDirectory(plugins_dir);
  
  LOG_PRINT(LOG_INFO,"OMM plugins successfully loaded from directory %s\n",plugins_dir.c_str());
  
  // deserialize the system
  system = unique_ptr<OpenMM::System>(&syst);
  
  // deserialize the integrator
  integrator = unique_ptr<OpenMM::Integrator>(&integr);
  
  // deserialize an OpenMM state : coordinates velocities forces energy
  addPlatform(platformDesired,state);
  
  state.getPeriodicBoxVectors(pbc[0],pbc[1],pbc[2]);
  
  memcpy(dat.pbc.data(),pbc.data(),sizeof(pbc));
  
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
    fprintf(stderr,"Error in %s line %d ; see error log file for details\n",__FILE__,__LINE__);
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  // fill dat structure with unserialised data
  dat.natom    = system->getNumParticles();
  dat.timestep = integrator->getStepSize();

  pos = vector<Vec3>(dat.natom);
  vel = vector<Vec3>(dat.natom);
}

void OMM_interface::addPlatform(PLATFORMS platformDesired, OpenMM::State& state)
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
  context = unique_ptr<OpenMM::Context> (new OpenMM::Context(*system, *integrator, *platform) );
    
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
  OpenMM::LocalEnergyMinimizer::minimize(*context, (int32_t)tol, maxsteps);
  
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
  
}

void OMM_interface::randomiseVelocities()
{
  // set velocities to initial temperature
  int32_t velSeed = get_int32();
  context->setVelocitiesToTemperature(dat.T,velSeed);
  LOG_PRINT(LOG_INFO,"Velocities randomised for T = %.3lf K with seed %d\n",dat.T,velSeed);
}

void OMM_interface::getState(double* timeInPs, ENERGIES* energies, double* currentTemperature,
                             ATOM atoms[])
{
  const bool wantTime        = !(timeInPs==nullptr);
  const bool wantEnergy      = !(energies==nullptr);
  const bool wantTemperature = !(currentTemperature==nullptr);
  const bool wantCrdVels     = !(atoms==nullptr);
  
  int32_t infoMask = 0;
  
  if(wantCrdVels)
    infoMask |= OpenMM::State::Positions | OpenMM::State::Velocities;
  
  if(wantEnergy)
    infoMask |= OpenMM::State::Energy;
  
  OpenMM::State state = context->getState(infoMask,true);
  
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
  {
    if(integType==LANGEVIN)
    {
      *currentTemperature = lint->getTemperature();
    }
    else if(integType==BROWNIAN)
    {
      *currentTemperature = bint->getTemperature();
    }
  }

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
  if(integType==LANGEVIN)
  {
    T = lint->getTemperature();
  }
  else if(integType==BROWNIAN)
  {
    T = bint->getTemperature();
  }

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
                                                                      const string& sysXMLfile,
                                                                      const string& integratorXMLfile,
                                                                      const string& stateXMLfile,
                                                                      const string& platformName
                                                                     )
{
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
    fprintf(stderr,"Error in %s line %d ; see error log file for details\n",__FILE__,__LINE__);
    MPI_CUSTOM_ABORT_MACRO();
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
    fprintf(stderr,"Error in %s line %d ; see error log file for details\n",__FILE__,__LINE__);
    MPI_CUSTOM_ABORT_MACRO();
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
    fprintf(stderr,"Error in %s line %d ; see error log file for details\n",__FILE__,__LINE__);
    MPI_CUSTOM_ABORT_MACRO();
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
  
  // create the c++ object corresponding to a c++ OpenMM interface
  OMM_interface* omm_itf = nullptr;
  omm_itf = new OMM_interface(_dat,*sys,*integ,st,platformDesired);
  
  LOG_PRINT(LOG_INFO,"OpenMM interface properly initialised from serialised files !\n");
  
  return unique_ptr<OMM_interface>(omm_itf);
}
