/**
 * \file lua_interface.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#include <iostream>
#include <fstream>
#include <chrono>

#include "sol.hpp"
#include "rand.hpp"
#include "logger.hpp"

#include "lua_interface.hpp"

#include "omm_interface.hpp"

#include "mpi_utils.hpp"

// function in lsqlite3.c in charge of registering lua interface to sqlite3
extern "C"
{
  int luaopen_lsqlite3(lua_State *L);
}

using namespace std;

// some shortcults when using tuples are defined here 
typedef tuple<double,double,double> CrdsTuple;
typedef tuple<double,double,double> VelsTuple;
typedef tuple<double,double,double> COMTuple;
typedef tuple<double,double,double> ENETuple;
typedef tuple<sol::table,sol::table> SolTabl2Tuple;

luaInterface::luaInterface(const string _inputFile,
                           DATA& _dat,
                           std::unique_ptr<ATOM[]>& _at,
                           std::unique_ptr<MD_interface>& _md)
                           : inputFile(_inputFile), dat(_dat), at(_at), md(_md)
{
  /* find out MY process ID, and how many processes were started. */
  my_id     = -1;
  num_procs = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  // load all available lua libraries by default
  lua.open_libraries();
  
  // register code from the embedded luasqlite3.c
  lua.require("sqlite3",luaopen_lsqlite3);

  register_default_lua_variables();
  register_default_lua_functions();

}

void luaInterface::parse_lua_file()
{
  string content;
  
  ifstream ifs(inputFile,ifstream::in);
  
  content = string( (istreambuf_iterator<char>(ifs) ),
                    (istreambuf_iterator<char>()  ) );

  ifs.close();
  

  // let the sol Lua interface read the whole input script
  lua.script(content);

  // max running time parsed here
  pvMap["max_run_time_hours"]                  = to_string(lua["max_run_time_hours"].get_or(24.0));
  pvMap["minutes_to_stop_before_max_run_time"] = to_string(lua["minutes_to_stop_before_max_run_time"].get_or(5));
  
  cout << "Parsed max_run_time_hours (or default if missing) : " << pvMap["max_run_time_hours"] << endl;
  cout << "Parsed minutes_to_stop_before_max_run_time (or default if missing) : " << pvMap["minutes_to_stop_before_max_run_time"] << endl;
  
  // TODO if we allow MD engines other than openMM in the future, changes will come here
  // get the OpenMM platform name
  pvMap["platform"] = lua["OMMplatform"].get_or(string("AUTO"));
  
  string& p = pvMap["platform"];
  if(p=="AUTO")
  {
  }
  else if(p == "REF")
    p = OMM_interface::omm_platforms_names.at(PLATFORMS::REF);
  else if(p == "CPU")
    p = OMM_interface::omm_platforms_names.at(PLATFORMS::CPU);
  else if(p == "OCL")
    p = OMM_interface::omm_platforms_names.at(PLATFORMS::OCL);
  else if(p == "CUDA")
    p = OMM_interface::omm_platforms_names.at(PLATFORMS::CUDA);
  else
  {
    cerr << "Error when parsing OMM platform : '" << p <<"' is not a valid platform type !" << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  cout << "Parsed OMM platform (or default if missing) : " << pvMap["platform"] << endl;
  
  // read the paths to the serialised OpenMM files
  sol::optional<string> integrator_path = lua["integrator"]["xml"];
  sol::optional<string> system_path     = lua["system"]["xml"];
  sol::optional<string> state_path      = lua["state"]["xml"];
  
  // handle errors
  if(integrator_path)
  {
    pvMap["integrator"] = integrator_path.value();
    cout << "Parsed integrator file : " << pvMap["integrator"] << endl;
  }
  else
  {
    cerr << "Error when reading integrator file name !" << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  if(system_path)
  {
    pvMap["system"] = system_path.value();
    cout << "Parsed system file : " << pvMap["system"] << endl;
  }
  else
  {
    cerr << "Error when reading system file name !" << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  if(state_path)
  {
    pvMap["state"] = state_path.value();
    cout << "Parsed state file : " << pvMap["state"] << endl;
  }
  else
  {
    cerr << "Error when reading state file name !" << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  // get total number of steps for this simulation
  sol::optional<uint64_t> numSteps = lua["numSteps"];
  if(numSteps)
  {
    pvMap["numSteps"] = lua["numSteps"];
    cout << "Parsed numSteps : " << pvMap["numSteps"] << endl;
  }
  else
  {
    cerr << "Error when reading number of steps !" << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }

  pvMap["equilibrationSteps"] = to_string(lua["equilibrationSteps"].get_or(50e3));
  
  pvMap["dbName"] = lua["database"]["name"].get_or(string("run.db"));
  pvMap["dbFreq"] = to_string(lua["database"]["backupFrequency"].get_or(500.0));
  
  /*
   * Parse mandatory lua functions (states definitions)
   */
  sol::optional<ParRep_function_state_init> mandatory_state_init = lua.get<ParRep_function_state_init>("state_init");
  if(mandatory_state_init)
  {
    func_state_init = mandatory_state_init.value();
  }
  else
  {
    cerr << "Error : cannot find the mandatory 'function state_init()' in Lua script !" << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }
    
  sol::optional<ParRep_function_check_state> mandatory_check_state_left = lua.get<ParRep_function_check_state>("check_state_left");
  if(mandatory_check_state_left)
  {
    func_check_state = mandatory_check_state_left.value();
  }
  else
  {
    cerr << "Error : cannot find the mandatory 'function check_state_left()' in Lua script !" << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  /*
   * more parsing for parRep parameters
   */
  // first need to know if std parRep or FV parRep are used
  sol::optional<string> algo_type = lua["simulation"]["algorithm"];
  if(algo_type)
  {
    pvMap["algorithm"] = algo_type.value();
  }
  else
  {
    cerr << "Error when reading type of algorithm for section 'simulation' !" << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }

  if(pvMap["algorithm"] == "PARREP")
  {
    cout << "Parsing ParRep simulation parameters." << endl;
    parse_parrep_params();
  }
  else if(pvMap["algorithm"] == "PARREP_FV")
  {
    cout << "Parsing ParRep FV simulation parameters." << endl;
    parse_parrepFV_params();
    
    uint32_t num_ex_events = (uint32_t) stoul(pvMap.at("allowedExitEvents"));
    
    if(num_ex_events > 1)
    {
      sol::optional<ParRep_function_get_serialized_state> mandatory_get_serialized_state = lua.get<ParRep_function_get_serialized_state>("get_serialized_state");
      if(mandatory_get_serialized_state)
      {
        func_get_serialized_state = mandatory_get_serialized_state.value();
      }
      else
      {
        cerr << "Error : cannot find the mandatory 'function get_serialized_state()' in Lua script !" << endl;
        MPI_CUSTOM_ABORT_MACRO();
      }
      
      sol::optional<ParRep_function_put_serialized_state> mandatory_put_serialized_state = lua.get<ParRep_function_put_serialized_state>("put_serialized_state");
      if(mandatory_put_serialized_state)
      {
        func_put_serialized_state = mandatory_put_serialized_state.value();
      }
      else
      {
        cerr << "Error : cannot find the mandatory 'function put_serialized_state()' in Lua script !" << endl;
        MPI_CUSTOM_ABORT_MACRO();
      }
    }
    else
    {
      lua["get_serialized_state"] = [](){};
      func_get_serialized_state = lua.get<ParRep_function_get_serialized_state>("get_serialized_state");
      
      lua["put_serialized_state"] = [](){};
      func_put_serialized_state = lua.get<ParRep_function_put_serialized_state>("put_serialized_state");
    }
    
  }
  else
  {
    cerr << "Error when reading type of algorithm : " << pvMap["algorithm"] << " is not a valid type of algortihm ! " << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  // get database interface lua functions
  sol::table sqldb = lua.get<sol::table>("SQLiteDB");
  db_open    = sqldb.get<SQLiteDB_open>("open");
  db_close   = sqldb.get<SQLiteDB_close>("close");
  db_insert  = sqldb.get<SQLiteDB_insert>("insert_state");
  db_backup  = sqldb.get<SQLiteDB_backup>("backup_to_file");
  
}

void luaInterface::register_default_lua_functions()
{
  /*
   * a set of lua function (defined with lambda functions) for accessing data or functions of this c++ code from lua script
   */
  
  // calling this from lua terminates the software, after a call to MPI_Finalize and LOG_FLUSH_ALL
  lua.set_function("exit_from_lua",
                   []()
                   {
                     LOG_FLUSH_ALL();
                     MPI_CUSTOM_ABORT_MACRO();
                   }
  );
  

  // expose high resolution c++ timer to lua for profiling
  typedef std::chrono::high_resolution_clock hr;
  lua.set_function("hr_timer",
                   []()->hr::time_point{return hr::now();}
  );
  
  // calculates time diff in nanoseconds from two timers
  lua.set_function("hr_timediff_ns",
                   [](hr::time_point before, hr::time_point after)
                   {
                     auto t = std::chrono::duration_cast<std::chrono::nanoseconds>(after-before);
                     return t.count();
                   }
  );
  
  // calculates time diff in microseconds from two timers
  lua.set_function("hr_timediff_us",
                   [](hr::time_point before, hr::time_point after)
                   {
                     auto t = std::chrono::duration_cast<std::chrono::microseconds>(after-before);
                     return t.count();
                   }
  );
  
  // calculates time diff in milliseconds from two timers
  lua.set_function("hr_timediff_ms",
                   [](hr::time_point before, hr::time_point after)
                   {
                     auto t = std::chrono::duration_cast<std::chrono::milliseconds>(after-before);
                     return t.count();
                   }
  );
  
  // get coordinates of one atom (values, so it is a safe read-only from lua side)
  lua.set_function("get_coordinates",
                   [&](size_t idx)->CrdsTuple{return CrdsTuple(at[idx-1].x(),at[idx-1].y(),at[idx-1].z());}
  );

  // get velocities of one atom  (values, so it is a safe read-only from lua side)
  lua.set_function("get_velocities",
                   [&](size_t idx)->VelsTuple{return VelsTuple(at[idx-1].vx(),at[idx-1].vy(),at[idx-1].vz());}
  );

  // get center of mass of the system
  lua.set_function("get_COM",
                   [&]()->COMTuple
                   {
                      XYZ cmass;

                      double sumMass = 0.0;

                      for(uint32_t i=0; i<dat.natom; i++)
                      {
                        cmass += at[i].crds*at[i].mass;
                        sumMass += at[i].mass;
                      }

                      cmass /= sumMass;
                      
                      return COMTuple(cmass.x(),cmass.y(),cmass.z());
                   }
  );
  
  // returns coordinates within a table (values, so it is a safe read-only from lua side)
  lua.set_function("get_all_coordinates",
                   [&]()->sol::table
                   {
                     sol::table crds = lua.create_table();
                     crds["x"] = lua.create_table();
                     crds["y"] = lua.create_table();
                     crds["z"] = lua.create_table();
                     
                     for(uint32_t i=1; i<=dat.natom; i++)
                     {
                       crds["x"][i] = at[i-1].x();
                       crds["y"][i] = at[i-1].y();
                       crds["z"][i] = at[i-1].z();
                     }
                     
                     return crds;
                   }
  );
  
  // returns velocities within a table (values, so it is a safe read-only from lua side)
  lua.set_function("get_all_velocities",
                   [&]()->sol::table
                   {
                     sol::table vels = lua.create_table();
                     vels["x"] = lua.create_table();
                     vels["y"] = lua.create_table();
                     vels["z"] = lua.create_table();
                     
                     for(uint32_t i=1; i<=dat.natom; i++)
                     {
                       vels["x"][i] = at[i-1].vx();
                       vels["y"][i] = at[i-1].vy();
                       vels["z"][i] = at[i-1].vz();
                     }
                     
                     return vels;
                   }
  );
  
  // returns both coordinates and velocities within 2 tables (values, so it is a safe read-only from lua side)
  lua.set_function("get_all_crdvels",
                   [&]() -> SolTabl2Tuple
                   {
                     sol::table crds = lua["get_all_coordinates"]();
                     sol::table vels = lua["get_all_velocities"]();
                     return SolTabl2Tuple(crds,vels);
                   }
  );
  
  // returns periodic boundary conditions as a table of 3 lua tables (a,b,c) (read only, safe)
  lua.set_function("get_pbc",
                   [&]() -> sol::table
                   {
                     sol::table a = lua.create_table();
                     sol::table b = lua.create_table();
                     sol::table c = lua.create_table();
                    
                     a["x"] = dat.a().x();   a["y"] = dat.a().y();   a["z"] = dat.a().z();
                     b["x"] = dat.b().x();   b["y"] = dat.b().y();   b["z"] = dat.b().z();
                     c["x"] = dat.c().x();   c["y"] = dat.c().y();   c["z"] = dat.c().z();
                    
                     sol::table pbc = lua.create_table();
                     pbc["a"] = a;
                     pbc["b"] = b;
                     pbc["c"] = c;
                     return pbc;
                   }
  );

  // set periodic boundary conditions ; lua side is a table of 3 lua tables (a,b,c), each with (x,y,z) members
  lua.set_function("set_pbc",
                   [&](sol::table pbc)
                   {
                      dat.a().x() = pbc["a"]["x"];
                      dat.a().y() = pbc["a"]["y"];
                      dat.a().z() = pbc["a"]["z"];
                      
                      dat.b().x() = pbc["b"]["x"];
                      dat.b().y() = pbc["b"]["y"];
                      dat.b().z() = pbc["b"]["z"];
                      
                      dat.c().x() = pbc["c"]["x"];
                      dat.c().y() = pbc["c"]["y"];
                      dat.c().z() = pbc["c"]["z"];
                   }
  );
  
  // set coordinates from lua to c++/openMM (and make sure the omm object is properly updated)
  //  WARNING possibly slow, use rarely
  lua.set_function("set_all_coordinates",
                   [&](sol::table crds)
                   {
                     for(uint32_t i=1; i<=dat.natom; i++)
                     {
                       at[i-1].x() = crds["x"][i];
                       at[i-1].y() = crds["y"][i];
                       at[i-1].z() = crds["z"][i];
                     }
                     md->setCrdsVels(at.get());
                   }
  );
  
  // set velocities from lua to c++/openMM (and make sure the omm object is properly updated)
  //  WARNING possibly slow, use rarely
  lua.set_function("set_all_velocities",
                   [&](sol::table vels)
                   {
                     for(uint32_t i=1; i<=dat.natom; i++)
                     {
                       at[i-1].vx() = vels["x"][i];
                       at[i-1].vy() = vels["y"][i];
                       at[i-1].vz() = vels["z"][i];
                     }
                     md->setCrdsVels(at.get());
                   }
  );

  // set coordinates and velocities from lua (and make sure the omm object is properly updated)
  //  WARNING possibly slow, use rarely
  lua.set_function("set_all_crdvels",
                   [&](sol::table crds, sol::table vels)
                   {
                     for(uint32_t i=1; i<=dat.natom; i++)
                     {
                       at[i-1].x() = crds["x"][i];
                       at[i-1].y() = crds["y"][i];
                       at[i-1].z() = crds["z"][i];
                       
                       at[i-1].vx() = vels["x"][i];
                       at[i-1].vy() = vels["y"][i];
                       at[i-1].vz() = vels["z"][i];
                     }
                     md->setCrdsVels(at.get());
                   }
  );
  
  // get center of mass of a subset of the system, using indices provided by the lua script
  lua.set_function("get_COM_idxs",
                   [&](sol::table idxs)->COMTuple
                   {
                      vector<int32_t> ref_idxs(idxs.size());
                      for(size_t i=1; i<= idxs.size(); i++)
                        ref_idxs[i-1] = idxs[i];

                      XYZ cmass;
                      double sumMass = 0.0;

                      for(int32_t i : ref_idxs)
                      {
                        cmass += at[i].crds*at[i].mass;
                        sumMass  += at[i].mass;
                      }

                      cmass /= sumMass;

                      return COMTuple(cmass.x(),cmass.y(),cmass.z());
                   }
  );

  // get mass of an atom (value, read-only)
  lua.set_function("get_mass",
                   [&](size_t idx)->double{return at[idx-1].mass;}
  );
  
  // perform energy minimisation, and return the 3 energy terms
  lua.set_function("get_minimised_energy",
                   [&](double tolerance, uint32_t maxSteps) -> ENETuple
                   {
                     ENERGIES minE;
                     md->minimiseWithCopy(tolerance,maxSteps,minE);
                     return ENETuple(minE.epot(),minE.ekin(),minE.etot());
                   }
  );

  // perform energy minimisation, and return the coordinates and velocities
  lua.set_function("get_minimised_crdvels",
                   [&](double tolerance, uint32_t maxSteps) -> SolTabl2Tuple
                  {
                    ENERGIES minE;
                    vector<XYZ> coordinates;
                    vector<XYZ> velocities;

                    md->minimiseWithCopy(tolerance,maxSteps,minE,coordinates,velocities);

                    sol::table crds = lua.create_table();
                    crds["x"] = lua.create_table();
                    crds["y"] = lua.create_table();
                    crds["z"] = lua.create_table();
                    
                    sol::table vels = lua.create_table();
                    vels["x"] = lua.create_table();
                    vels["y"] = lua.create_table();
                    vels["z"] = lua.create_table();
                    
                    for(uint32_t i=1; i<=dat.natom; i++)
                    {
                      crds["x"][i] = coordinates[i-1][0];
                      crds["y"][i] = coordinates[i-1][1];
                      crds["z"][i] = coordinates[i-1][2];
                      
                      vels["x"][i] = velocities[i-1][0];
                      vels["y"][i] = velocities[i-1][1];
                      vels["z"][i] = velocities[i-1][2];
                    }

                    return SolTabl2Tuple(crds,vels);
                  }
);
  
  // get the current temperature (from kinetic energy)
  lua.set_function("get_temperature",
                   [&]()->double{ return md->getTemperature();}
  );

}

void luaInterface::register_default_lua_variables()
{
  
  lua["restart"] = true;
  
  lua["mpi_rank_id"]    = my_id;
  lua["mpi_num_ranks"]  = num_procs;
  
  lua["natoms"]   = 0;
  lua["timeStep"] = 0.0;
  lua["temperature"] = 0.0;
  
  lua.create_named_table("minimisation");
  lua["minimisation"]["Tolerance"] = 1e-6;
  lua["minimisation"]["MaxSteps"]  = 0;
  
  lua["equilibrationSteps"] = 50e3;
  
  lua.create_named_table("database");
  lua["database"]["name"] = "run.db";
  lua["database"]["backupFrequency"] = 500.0;
  
  lua["epot"] = 0.0;
  lua["ekin"] = 0.0;
  lua["etot"] = 0.0;
  lua["referenceTime"] = 0.0;
  
  lua.create_named_table("SQLiteDB");
  
  // expose dummy lambda functions ; they do nothing
  //  and are the default in case the user does not override them in the lua script
  lua["SQLiteDB"]["open"].set_function(
    []()->void
    {
      cerr << "Dummy lambda db_open\n";
    }
  );
  
  lua["SQLiteDB"]["close"].set_function(
    []()->void
    {
      cerr << "Dummy lambda db_close\n";
    }
  );
  
  lua["SQLiteDB"]["insert_state"].set_function(
    []()->void
    {
      cerr << "Dummy lambda db_insert\n";
    }
  );
  
  lua["SQLiteDB"]["backup_to_file"].set_function(
    []()->void
    {
      cerr << "Dummy lambda db_backup\n";
    }
  );
  
}

void luaInterface::parse_parrep_params()
{
  sol::table simulation_def = lua.get<sol::table>("simulation");

  pvMap["tauDecorr"] = to_string(simulation_def["tauDecorr"].get_or(10.0));
  pvMap["checkDecorr"] = to_string(simulation_def["checkDecorr"].get_or(100));
  

  pvMap["tauDephase"] = to_string(simulation_def["tauDephase"].get_or(10.0));
  pvMap["checkDephase"] = to_string(simulation_def["checkDephase"].get_or(100));
  

  pvMap["checkDynamics"] = to_string(simulation_def["checkDynamics"].get_or(100));
  
}

void luaInterface::parse_parrepFV_params()
{
 
  sol::table simulation_def = lua.get<sol::table>("simulation");
  
  pvMap["allowedExitEvents"] = to_string(simulation_def["allowedExitEvents"].get_or(1));
  
  pvMap["afterExit"] = simulation_def["afterExit"].get_or(string("wait"));
  
  pvMap["checkFV"] = to_string(simulation_def["checkFV"].get_or(100));
  pvMap["checkGR"] = to_string(simulation_def["checkGR"].get_or(1));
  pvMap["checkDynamics"] = to_string(simulation_def["checkDynamics"].get_or(100));
  
  sol::optional<sol::table> gr_functions_list = simulation_def.get<sol::table>("GRobservables");
  if(gr_functions_list)
  {
    sol::table& tab = gr_functions_list.value();
    for(size_t n=1; n<=tab.size(); n++)
    {
      const GR_function_name fname = tab[n];
      cout << "Found a GR function declaration : " << fname << endl;
      gr_funcs_names.push_back(fname);
      GR_function f = lua.get<GR_function>(fname);
      gr_funcs[fname] = f;
    }
  }
  else
  {
    cerr << "Error when trying to read simulation.GRobservables which is required for a 'parrepFV' simulation !" << endl;
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  pvMap["GRtol"] = to_string(simulation_def["GRtol"].get_or(0.01));
  
}
