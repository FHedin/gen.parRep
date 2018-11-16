/**
 * Copyright (c) 2016-2018, Florent Hédin, Tony Lelièvre, and École des Ponts - ParisTech
 * All rights reserved.
 * 
 * The 3-clause BSD license is applied to this software.
 * 
 * See LICENSE.txt
 */

#include <memory>
#include <iostream>
#include <sstream>
#include <chrono>
#include <exception>

#include "runSim.hpp"

#include "logger.hpp"

#include "md_interface.hpp"
#include "omm_interface.hpp"

#include "lua_interface.hpp"

#include "parRep.hpp"
#include "parRep_FV.hpp"
#include "parRep_FV_multiExits.hpp"

#include "mpi_utils.hpp"

#include "rand.hpp"

using namespace std;

void run_simulation(const string& inpf)
{
  /* find out MY process ID, and how many processes were started. */
  int32_t my_id     = -1;
  int32_t num_procs = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  /*
   * Each MPI rank will initialise its own random numbers generator ; see rand.hpp
   */
  init_rand();
  
  /*
   * This is how to initialize the random generators providing seeds manually (at least 5 seeds of type uint64_t encodes as a comma separated string)
   * 
   * ATTENTION Using the directly the following line would be wrong because all the random generators on all
   *           MPI ranks would be initialised with the same seeds !!! Each rank should use a modified vrsion of those seeds.
   * 
   *   const string ten_seeds = "12813764096581876017,1660537369292506360,10652021971297968194,"
   *    "1632182917282226681,2681620832946811191,11894265520984226376,"
   *    "8987030001799979889,7928564290642933391,6948006017539476486,1814061233159141676";
   * 
   *    init_rand(ten_seeds);
   * 
   */

  fprintf(stdout,"Rank %d properly initialised its mt19937 random numbers generator\n",my_id);
  fprintf(stdout,"5 random uint32_t from rank %d : %u %u %u %u %u\n",my_id,get_uint32(),
          get_uint32(),get_uint32(),get_uint32(),get_uint32());
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  // stores simulation parameters
  DATA dat;
  
  // start the chronometer in DATA now
  dat.start_time = chrono::steady_clock::now();
  
  // for storing coordinates and velocities
  unique_ptr<ATOM[]> atmList(nullptr);
  
  // pointer to a MD interface, this can contain any interface to any engine
  unique_ptr<MD_interface> md(nullptr);
  
  // initialise the code parsing the lua input file
  unique_ptr<luaInterface> luaItf(nullptr);
  try
  {
    luaItf = unique_ptr<luaInterface>( new luaInterface(inpf,dat,atmList,md) );
  }
  catch(exception& e)
  {
    fprintf(stderr,"std::exception captured on rank %d when initialising Lua interface !\n",my_id);
    fprintf(stderr,"Error message is : %s\n",e.what());
    fprintf(stderr,"Now flushing files and terminating all other MPI ranks...\n");
    LOG_FLUSH_ALL();
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  // parse the file
  try
  {
    luaItf->parse_lua_file();
  }
  catch(exception& e)
  {
    fprintf(stderr,"std::exception captured on rank %d when parsing the Lua input file !\n",my_id);
    fprintf(stderr,"Error message is : %s\n",e.what());
    fprintf(stderr,"Now flushing files and terminating all other MPI ranks...\n");
    LOG_FLUSH_ALL();
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  // get lua parsed parameters and map them to the dat structure
  const lua_ParVal_map& luaParams = luaItf->get_parsed_parameters_map();
  
  // get max run time to know when to interrupt properly the software instead of crashing on queuing systems
  dat.max_run_time_hours = chrono::duration<double, chrono::hours::period>(stod(luaParams.at("max_run_time_hours")));
  
  dat.minutes_to_stop_before_max_run_time = chrono::minutes(stoul(luaParams.at("minutes_to_stop_before_max_run_time")));
  
  // get the number of steps of simulation to perform with MD engine ; unsigned 64 bits for allowing theoretically more than 2 billions
  dat.nsteps = stoul(luaParams.at("numSteps"));
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  // First retrieve the type of MD engine from the Lua input file
  const string required_engine_name = luaParams.at("MD_Engine");
  
  MD_ENGINES md_engine_type = MD_ENGINES::UNKNOWN_ENGINE;
  if(required_engine_name == "OpenMM")
    md_engine_type = MD_ENGINES::OPENMM;
  else
    md_engine_type = MD_ENGINES::UNKNOWN_ENGINE;
  
  fprintf(stdout,    "User required to perform simulation using the '%s' MD engine.\n",required_engine_name.c_str());
  LOG_PRINT(LOG_INFO,"User required to perform simulation using the '%s' MD engine.\n",required_engine_name.c_str());
  
  // then call the appropriate class constructor depending on the type of parsed MD engine
  try
  {
    switch(md_engine_type)
    {
      case MD_ENGINES::OPENMM :
        md =  unique_ptr<MD_interface>( OMM_interface::initialise_fromSerialization(dat,luaItf) );
      break;
      
      case MD_ENGINES::UNKNOWN_ENGINE:
        throw runtime_error("The MD engine was set to the default error type 'MD_ENGINES::UNKNOWN_ENGINE'\n");
        break;
    }
  }
  catch(exception& e)
  {
    fprintf(stderr,"std::exception captured on rank %d when initialising the MD engine !\n",my_id);
    fprintf(stderr,"Error message is : %s\n",e.what());
    fprintf(stderr,"Now flushing files and terminating all other MPI ranks...\n");
    LOG_FLUSH_ALL();
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  // from this point we can allocate the c++ side storage for atom crds and vels
  atmList = unique_ptr<ATOM[]>(new ATOM[dat.natom]);
  
  // copy back initial configuration from MD engine to the ATOM[] storage ; also get initial temperature and energies
  ENERGIES iniE;
  md->getSimData(nullptr,&iniE,&(dat.T),atmList.get());
  // getParticlesParams(...) allows us to also retrieve the atomic mass
  md->getParticlesParams(atmList.get());
  
  // export back to lua dat.natom, dat.timestep, dat.T
  luaItf->set_lua_variable("natoms",dat.natom);
  luaItf->set_lua_variable("timeStep",dat.timestep);
  luaItf->set_lua_variable("temperature",dat.T);
  
  // push to Lua some initial value for energy and physical time
  luaItf->set_lua_variable("epot",iniE.epot());
  luaItf->set_lua_variable("ekin",iniE.ekin());
  luaItf->set_lua_variable("etot",iniE.etot());
  luaItf->set_lua_variable("referenceTime",0.0);
  
  bool useFlemingViot=true;
  if(luaParams.at("algorithm")=="PARREP") useFlemingViot=false;
    
  // summary of the simulation
  fprintf(stdout,"\nStarting program with %d MPI ranks\n\n",num_procs);
  
  fprintf(stdout,"This is : rank[%d]\n\n",my_id);

  MPI_Barrier(MPI_COMM_WORLD);
  
  /*
   * PARREP_SIMULATION encapsulates a pointer to a ParRepAbstract class so that it
   * can store any type of parrep simulation 
   * no need to delete it, it is done automatically
   */
  PARREP_SIMULATION simulation = nullptr;
  
  try
  {
    if(useFlemingViot)
    {
      const uint32_t exit_events = stoul(luaParams.at("allowedExitEvents"));
      
      if(exit_events > 1)
        simulation = PARREP_SIMULATION(new ParRepFV_multiExits(dat,atmList,md,luaItf,exit_events));
      else
        simulation = PARREP_SIMULATION(new ParRepFV(dat,atmList,md,luaItf));
    }
    else
    {
      simulation = PARREP_SIMULATION(new ParRep(dat,atmList,md,luaItf));
    }
    
  }
  catch(exception& e)
  {
    fprintf(stderr,"std::exception captured on rank %d when creating PARREP_SIMULATION object !\n",my_id);
    fprintf(stderr,"Error message is : %s\n",e.what());
    fprintf(stderr,"Now flushing files and terminating all other MPI ranks...\n");
    LOG_FLUSH_ALL();
    MPI_CUSTOM_ABORT_MACRO();
  }

  try
  {
    simulation->run();
  }
  catch(exception& e)
  {
    fprintf(stderr,"std::exception captured on rank %d when running simulation !\n",my_id);
    fprintf(stderr,"Error message is : %s\n",e.what());
    fprintf(stderr,"Now flushing files and terminating all other MPI ranks...\n");
    LOG_FLUSH_ALL();
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
}



