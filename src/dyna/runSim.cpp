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

#include "runSim.hpp"

#include "logger.hpp"
#include "md_interface.hpp"
#include "omm_interface.hpp"

#include "lua_interface.hpp"

#include "parRep.hpp"
#include "parRep_FV.hpp"
#include "parRep_FV_multiExits.hpp"

#include "mpi_utils.hpp"

// GCC / G++ only
#ifdef __GNUC__

#include <execinfo.h>

/* Obtain a backtrace and print it to stdout. */
static void print_trace()
{
  void *array[200];
  int size;
  char **strings;
  
  size = backtrace(array, 10);
  strings = backtrace_symbols(array, size);
  
  LOG_PRINT(LOG_ERROR,"Obtained %zd stack frames.\n", size);
  
  for (int i = 0; i < size; i++)
    LOG_PRINT(LOG_ERROR,"%s\n", strings[i]);
  
  free (strings);
}

#endif

using namespace std;

void run_simulation(const string& inpf)
{
  /* find out MY process ID, and how many processes were started. */
  int32_t my_id     = -1;
  int32_t num_procs = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
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
    fprintf(stderr,"std::exception captured when initialising Lua interface ! See error log file for more details.\n");
    LOG_PRINT(LOG_ERROR,"A std::exception has been emmited on rank %d when initialising the Lua interface.\n",my_id);
    LOG_PRINT(LOG_ERROR,"Error message is : %s\n",e.what());
    LOG_PRINT(LOG_ERROR,"Now flushing files and terminating all other MPI jobs...\n");
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
    fprintf(stderr,"std::exception captured when parsing Lua input file ! See error log file for more details.\n");
    LOG_PRINT(LOG_ERROR,"A std::exception has been emmited on rank %d when parsing the Lua input file.\n",my_id);
    LOG_PRINT(LOG_ERROR,"Error message is : %s\n",e.what());
    LOG_PRINT(LOG_ERROR,"Now flushing files and terminating all other MPI jobs...\n");
    LOG_FLUSH_ALL();
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  // get lua parsed parameters and map them to the dat structure
  const lua_ParVal_map& luaParams = luaItf->get_parsed_parameters_map();
  
  // get max run time to know when to interrupt properly the software instead of crashing on queuing systems
  dat.max_run_time_hours = chrono::duration<double, chrono::hours::period>(stod(luaParams.at("max_run_time_hours")));
  
  dat.minutes_to_stop_before_max_run_time = chrono::minutes(stoul(luaParams.at("minutes_to_stop_before_max_run_time")));
  
  // get desired OpenMM platform
  string platform =   luaParams.at("platform");
  
  // get the number of steps of simulation to perform with openmm ; unsigned 64 bits for allowing theoretically more than 2 billions
  dat.nsteps = stoul(luaParams.at("numSteps"));
  
  // get paths to OMM XML serialised files
  string integrator = luaParams.at("integrator");
  string system =     luaParams.at("system");
  string state =      luaParams.at("state");
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  // TODO here allow possibility of other MD engine when ready
  try
  {
    md =  unique_ptr<MD_interface>(
      OMM_interface::initialise_fromSerialization(dat,
                                                  system,
                                                  integrator,
                                                  state,
                                                  platform)
    );
  }
  catch(exception& e)
  {
    fprintf(stderr,"std::exception captured when initialising OpenMM simulation ! See error log file for more details.\n");
    LOG_PRINT(LOG_ERROR,"A std::exception has been emmited on rank %d.\n",my_id);
    LOG_PRINT(LOG_ERROR,"Error message is : %s\n",e.what());
    #ifdef __GNUC__
    print_trace();
    #endif
    LOG_PRINT(LOG_ERROR,"Now flushing files and terminating all other MPI jobs...\n");
    LOG_FLUSH_ALL();
    MPI_CUSTOM_ABORT_MACRO();
  }

  // dat.natom, dat.timestep, dat.T are read from OMM xml files
  // we redefine them back to lua
  luaItf->set_lua_variable("natoms",dat.natom);
  luaItf->set_lua_variable("timeStep",dat.timestep);
  luaItf->set_lua_variable("temperature",dat.T);
    
  // from this point system was unserialised and we know the number of atoms : allocate memory
  atmList = unique_ptr<ATOM[]>(new ATOM[dat.natom]);
  
  // copy back initial OpenMM configuration to the ATOM[] ; also get initial temperature (may fluctuate)
  ENERGIES iniE;
  md->getState(nullptr,&iniE,&(dat.T),atmList.get());
  md->getParticlesParams(atmList.get());
  
  // push to Lua some initial value for global variables
  luaItf->set_lua_variable("epot",iniE.epot());
  luaItf->set_lua_variable("ekin",iniE.ekin());
  luaItf->set_lua_variable("etot",iniE.etot());
  luaItf->set_lua_variable("referenceTime",0.0);
  
  bool useFlemingViot=true;
  if(luaParams.at("algorithm")=="PARREP") useFlemingViot=false;
    
  // summary of the simulation
  fprintf(stdout,"\nStarting program with %d MPI ranks\n\n",num_procs);
  
  fprintf(stdout,"This is : rank[%d]\n\n",my_id);
  
  const string& ommVersion = OMM_interface::getOMMversion();
  fprintf(stdout,"Using OpenMM toolkit version '%s'\n",ommVersion.c_str());

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
    fprintf(stderr,"std::exception captured when creating PARREP_SIMULATION object ! See error log file for more details.\n");
    LOG_PRINT(LOG_ERROR,"An std::exception has been emmited on rank %d.\n",my_id);
    LOG_PRINT(LOG_ERROR,"Error message is : %s\n",e.what());
    #ifdef __GNUC__
    print_trace();
    #endif
    LOG_PRINT(LOG_ERROR,"Now flushing files and terminating all other MPI jobs...\n");
    LOG_FLUSH_ALL();
    MPI_CUSTOM_ABORT_MACRO();
  }

  try
  {
    simulation->run();
  }
  catch(exception& e)
  {
    fprintf(stderr,"std::exception captured when running simulation ! See error log file for more details.\n");
    LOG_PRINT(LOG_ERROR,"An std::exception has been emmited on rank %d.\n",my_id);
    LOG_PRINT(LOG_ERROR,"Error message is : %s\n",e.what());
    #ifdef __GNUC__
    print_trace();
    #endif
    LOG_PRINT(LOG_ERROR,"Now flushing files and terminating all other MPI jobs...\n");
    LOG_FLUSH_ALL();
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
}



