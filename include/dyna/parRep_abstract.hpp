/**
 * \file parRep_abstract.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#ifndef PARREP_ABSTRACT_HPP_INCLUDED
#define PARREP_ABSTRACT_HPP_INCLUDED

#include <cstdlib>
#include <memory>

#include "global.hpp"
#include "logger.hpp"
#include "md_interface.hpp"
#include "lua_interface.hpp"

#include "mpi_utils.hpp"

/**
 * \brief All parRep derived class will inherit from this one (directly or via ParRepAbstract below)
 *        which just contains the run method. It is also used as the common type for the below defined
 *        PARREP_SIMULATION typedef.
 * 
 * \note purely virtual class but it is important to define a virtual
 *       destructor for avoiding memory leaks when using polymorphism
 */
class ParRepRunnable
{

public:

  // virtual dtor is required for avoiding memory leaks when using polymorphism in pointers (see PARREP_SIMULATION below)
  virtual ~ParRepRunnable(){};
  
  /**
   * @brief the main function running parallel replica simulations
   *
   * pure virtual method, each parRep will redefine its own
   *
   */
  virtual void run() = 0;
  
protected:

};

/**
 * PARREP_SIMULATION is a wrapped pointer to a ParRepRunnable object.
 * 
 * It allows to put any type of parrep simulation within a given variable, in order to allows polymorphism
 */
typedef std::unique_ptr<ParRepRunnable> PARREP_SIMULATION;

/**
 * @brief Abstract class representing a parRep simulation.
 * 
 * Contains methods and attributes common to all ParRep variants
 * 
 */
class ParRepAbstract : public ParRepRunnable
{

public:

  /**
   * @brief The constructor for an abstract parrep simulation, called by classes inheriting from this one.
   * 
   * It will initialise many of the references common to all parrep simulations.
   * 
   * It is also in charge of performing some MPI setup.
   * 
   * @param _dat Simulation data, see global.hpp
   * @param _at Array of coordinates and velocities
   * @param _md The MD interface
   * @param _luaItf The Lua interface
   * @param ignore_local_setup If true do not perform some of the local setup and let the derived class do it
   */
  ParRepAbstract(DATA& _dat,
                 std::unique_ptr<ATOM[]>& _at,
                 std::unique_ptr<MD_interface>& _md,
                 std::unique_ptr<luaInterface>& _luaItf,
                 bool ignore_local_setup=false
                )
                : dat(_dat),at(_at),md(_md),luaItf(_luaItf),
                params(luaItf->get_parsed_parameters_map()),
                lua_state_init(luaItf->get_function_state_init()),
                lua_check_state(luaItf->get_function_check_state()),
                lua_check_transient(luaItf->get_function_check_transient()),
                db_open(luaItf->get_db_open()),
                db_close(luaItf->get_db_close()),
                db_insert(luaItf->get_db_insert()),
                db_backup(luaItf->get_db_backup())
  {
    if(ignore_local_setup==false)
    {
      local_setup_done = true;
      MPI_setup();
    }
  }
  
  /**
  * @brief Destructor with explicit call to some methods
  * 
  */
  ~ParRepAbstract()
  {
    if(local_setup_done)
    {
      MPI_clean();
    }
  }

  /**
  * @brief Run method from ParRepRunnable is still marked as purely virtual in order to force
  *        derived classes to reimplement it.
  */
  virtual void run() override = 0;

protected:

  // references obtained from constructor
  DATA& dat;                              ///< Simulation data, see global.hpp
  std::unique_ptr<ATOM[]>& at;            ///< Array of coordinates and velocities 
  std::unique_ptr<MD_interface>& md;      ///< The interface to the MD engine
  std::unique_ptr<luaInterface>& luaItf;  ///< The Lua interface

  // references to objects used by all parrep variants
  const lua_ParVal_map& params;                                 ///< Map of Lua parameters
  const ParRep_function_state_init&      lua_state_init;        ///< Lua function for state initialization
  const ParRep_function_check_state&     lua_check_state;       ///< Lua function for state checking
  const ParRep_function_check_transient& lua_check_transient;   ///< Lua function for state checking
  
  const SQLiteDB_open&    db_open;                        ///< Lua function for database opening
  const SQLiteDB_close&   db_close;                       ///< Lua function for database closeing
  const SQLiteDB_insert&  db_insert;                      ///< Lua function for database insert
  const SQLiteDB_backup&  db_backup;                      ///< Lua function for database backup

  /*
   * some MPI parameters and variables common to all parRep algorithms
   */
  int32_t masterRank=-1;    ///< the rank of the MPI master process ; set at initialization by MPI_setup()
  bool i_am_master = false; ///< Is this rank the masterRank ? set at initialization by MPI_setup()
  
  MPI_Comm global_comm = MPI_COMM_NULL;     ///< The global communicator: by default a simple cloning of MPI_COMM_WORLD, but may be a more restricted set of ranks (may be setup in a possible override of MPI_setup)
  MPI_Group global_group = MPI_GROUP_NULL;   ///< The global group associated to the global_comm
  
  MPI_Info rma_info = MPI_INFO_NULL;        ///< This holds optmization hints communicated to MPI when creating RMA windows
  
  int32_t mpi_id_in_gcomm = -1;  ///< the id of this mpi rank in global_comm
  int32_t mpi_gcomm_size = -1;   ///< the total number of ranks in global_comm

  /*
   * some simulation parameters common to all parRep algorithms
   */
  
  /**
   * this is the reference clock time of the simulation
   *  see articles for details
   */
  double ref_clock_time = 0.;
 
  uint32_t equil_steps;           ///< Number of equilibration steps to perform
  double equil_local_time = 0.;   ///< Physical equilibation time
  ENERGIES equil_e;               ///< a structure for storing equilibrated energies

  // parameters for parallel exit phase
  double dyna_local_time = 0.;    ///< Physical parallel dynamics time
  ENERGIES dyna_e;               ///< for storing parallel dynamics energies
  uint32_t dyna_cycles_done = 0; ///<to know how many cycles of t_poll steps have been performed in the parallel exit phase

  uint32_t breakerID = std::numeric_limits<uint32_t>::max(); ///< for storing the mpi rank id of the exiting replica : defaut to the max possible uint32_t value because we will do a MPI_AllReduce minimum search
  uint32_t dynamics_check;        ///< frequency (in steps) at which to check if there was a state exit
  double t_poll;                  ///< time intervall (physical units e.g. ps) at which to check if a replica left the current state during the parallel phase

  bool     left_state = false;    ///< for checking if the system left the reference state

  // database parameters
  double   db_backup_frequency_ps;  ///< frequency for dumping the in memory db to a file
  double   last_db_backup = 0.0;    ///<  variable holding the last time at which such a backup was performed

  /**
   * @brief This function is in charge of initialising some MPI variables used by the ParRep code
   * 
   * each derived class may re-implement its own; if overriden, use the second constructor of this class and call yourself MPI_setup in derived class ctor
   * 
   * By default, just sets masterRank to 0 and creates a copy of the MPI_COMM_WORLD and of its associated MPI_group
   * 
   */
  virtual void MPI_setup()
  {
    LOG_PRINT(LOG_DEBUG,"Call of the default MPI_setup() defined in file %s at line %d \n",__FILE__,__LINE__);
    
    // first a copy of the original MPI_COMM_WORLD is done
    MPI_Comm_dup(MPI_COMM_WORLD,&global_comm);
    // extract the original group handle
    MPI_Comm_group(global_comm, &global_group);
    
    MPI_Comm_rank(global_comm, &mpi_id_in_gcomm);
    MPI_Comm_size(global_comm, &mpi_gcomm_size);
    
    masterRank = 0;
    i_am_master = (mpi_id_in_gcomm == masterRank);
    
    // set some optimization flags for MPI_Win RMA operations
    MPI_Info_create(&rma_info);
    MPI_Info_set(rma_info,"same_size","true");
    MPI_Info_set(rma_info,"same_disp_unit","true");
    
  }
  
  /**
   * @brief This function is in charge of cleaning some MPI variables used by the ParRep code
   * 
   * each derived class may re-implement its own; if overriden, call yourself MPI_clean in derived class dtor
   */
  virtual void MPI_clean()
  {
    LOG_PRINT(LOG_DEBUG,"Call of the default MPI_clean() defined in file %s at line %d \n",__FILE__,__LINE__);
    
    MPI_Group_free(&global_group);
    MPI_Comm_free(&global_comm);
    MPI_Info_free(&rma_info);
  }
  
  /**
   * @brief Performs an equilibration of the reference structure at the beginning of simulation
   */
  virtual void do_equilibration()
  {
    LOG_PRINT(LOG_INFO,"Equilibration on rank %d for %.2lf ps\n",  mpi_id_in_gcomm, equil_steps*dat.timestep);
    fprintf(stdout,    "Equilibration on rank %d for %.2lf ps\n\n",mpi_id_in_gcomm, equil_steps*dat.timestep);
    
    // equilibrate
    md->doNsteps(equil_steps);
    equil_local_time = (double) equil_steps * dat.timestep;
    
    // we also reset the md time as if there was no equilibration before
    md->setSimClockTime(0.);
    
    // get coordinates
    md->getSimData(nullptr,&equil_e,nullptr,at.get());
  }
  
  /**
   * @brief If the current is outside of any known state,
   *        propagate it until it enters a new known metstable state.
   * 
   * Only one replica does that; then it will broadcast its new reference time and the new system to all the others
   * 
   */
  virtual void do_transient_propagation(ENERGIES& tr_e)
  {
    LOG_PRINT(LOG_INFO,"Rank %d performs transient propagation...\n",mpi_id_in_gcomm);
    fprintf(stdout,    "Rank %d performs transient propagation...\n",mpi_id_in_gcomm);
    
    bool need_transient_propagation = false;
    double time_spent_in_transient_propagation = 0.0;
    while(true)
    {
      md->doNsteps(dynamics_check);
      time_spent_in_transient_propagation += (double) dynamics_check * dat.timestep;
      
      LOG_PRINT(LOG_DEBUG,"Rank %d did %lf ps of transient propagation...\n",mpi_id_in_gcomm,time_spent_in_transient_propagation);
      fprintf(stdout,     "Rank %d did %lf ps of transient propagation...\n",mpi_id_in_gcomm,time_spent_in_transient_propagation);
      
      md->getSimData(nullptr,&tr_e,nullptr,at.get());
      luaItf->set_lua_variable("epot",tr_e.epot());
      luaItf->set_lua_variable("ekin",tr_e.ekin());
      luaItf->set_lua_variable("etot",tr_e.etot());
      
      need_transient_propagation = lua_check_transient();
      
      if(!need_transient_propagation)
        break;
    }
    
    ref_clock_time += time_spent_in_transient_propagation;
    
    LOG_PRINT(LOG_INFO,"Rank %d entered a new metastable state after %lf ps of transient propagation!\n",
              mpi_id_in_gcomm,time_spent_in_transient_propagation);
    fprintf(stdout,    "Rank %d entered a new metastable state after %lf ps of transient propagation!\n",
            mpi_id_in_gcomm,time_spent_in_transient_propagation);

  }
  
  /**
   * @brief This performs the parallel dynamics stage ; executed by all ranks at the same time,
   * and stopping as soon as one of them escapes
   */
  virtual void do_dynamics()
  {
    bool needToBreak = false;
    MPI_Win break_window = MPI_WIN_NULL;
//     MPI_Request ibarrier_req = MPI_REQUEST_NULL;
    
    // NOTE there is an implicit barrier in MPI_Win_create
    MPI_Win_create(&needToBreak,sizeof(bool),sizeof(bool),
                   rma_info,global_comm,&break_window);
    
    LOG_PRINT(LOG_INFO,"Rank %d now entering parallel dynamics stage\n",  mpi_id_in_gcomm);
    fprintf(stdout,    "Rank %d now entering parallel dynamics stage\n\n",mpi_id_in_gcomm);
    
    while(true)
    {
      // do some steps
      md->doNsteps(dynamics_check);
      
//       MPI_Wait(&ibarrier_req,MPI_STATUS_IGNORE);
//       
//       if(needToBreak)
//         break;
      
      dyna_local_time += t_poll;
      
      md->getSimData(nullptr,&dyna_e,nullptr,at.get());
      luaItf->set_lua_variable("epot",dyna_e.epot());
      luaItf->set_lua_variable("ekin",dyna_e.ekin());
      luaItf->set_lua_variable("etot",dyna_e.etot());
      
      left_state = lua_check_state();
      
      /*
       * We check if this rank left the state
       */
      if(left_state)
      {
        LOG_PRINT(LOG_INFO,"Rank %d escaped the state after %.2lf ps\n",  mpi_id_in_gcomm,dyna_local_time);
        fprintf(stdout,    "Rank %d escaped the state after %.2lf ps\n\n",mpi_id_in_gcomm,dyna_local_time);
        
        // if it left the state, save the mpi_id_in_gcomm of this rank and broadcast it later 
        breakerID = (uint32_t) mpi_id_in_gcomm;
        needToBreak = true;
        
        /*
         * Using a MPI_Win to notify the others that they should stop
         */
        MPI_Win_lock_all(0,break_window);
        for(int32_t i=0; i<mpi_gcomm_size; i++)
        {
          if(i == mpi_id_in_gcomm)
            continue;
          
          MPI_Put(&needToBreak,   //const void *origin_addr
                  1,MPI_CXX_BOOL, //int origin_count, MPI_Datatype origin_datatype,
                  i,0,            // int target_rank, MPI_Aint target_disp,
                  1,MPI_CXX_BOOL, // int target_count, MPI_Datatype target_datatype,
                  break_window);  // MPI_Win win
        }
        MPI_Win_unlock_all(break_window);
      }
      
      // instead use a non blocking barrier and check after the next doNsteps in order to hide the sync cost behind computations
//       MPI_Ibarrier(global_comm,&ibarrier_req);
      
      MPI_Barrier(global_comm);
      
      if(needToBreak)
        break;
      
      LOG_PRINT(LOG_INFO,"Still within the same state after %.2lf ps...\n",  dyna_local_time);
      if(LOG_SEVERITY == LOG_DEBUG)
        fprintf(stdout,  "Still within the same state after %.2lf ps...\n\n",dyna_local_time);
      
      dyna_cycles_done++;
      
      // check running time and depending on return value perform a backup before exiting
      if(check_running_time())
        doBackupAndExit();
      
    } // loop on while
    
    /*
     * the breaker waits here, the others should arrive as soon has they are notified via needToBreak, which will at most take one more
     * md evaluation (dynamics_check)
     */
    MPI_Barrier(global_comm);
    
    MPI_Win_free(&break_window);
    
  } // end of do_dynamics
  
  /**
   * @brief Checks if the simulation have been running for dat.max_run_time_hours
   * 
   * In fact it tries to shutdown properly the simulation when getting close to max_run_time_hours
   * (see variable dat.minutes_to_stop_before_max_run_time),
   * which is useful when running on a cluster with a time limited queuing system
   * 
   * @return TRUE if simulation should stop
   */
  virtual bool check_running_time()
  {
    // get current time
    auto now = std::chrono::steady_clock::now();
    
    // check for how many minutes simulation has been running
    std::chrono::minutes minutes_elapsed = std::chrono::duration_cast<std::chrono::minutes>(now - dat.start_time);
    
    // check for how long (in minutes) simulation is allowed to tun
    std::chrono::minutes max_allowed_minutes = std::chrono::duration_cast<std::chrono::minutes>(dat.max_run_time_hours);
    max_allowed_minutes -= dat.minutes_to_stop_before_max_run_time;
    
    bool should_stop=false;
    
    int32_t elapsed = minutes_elapsed.count();
    int32_t allowed = max_allowed_minutes.count();
    
    if(elapsed < allowed)
    {
      LOG_PRINT(LOG_DEBUG,"Simulation has been running for %d minutes of an allowed maximum of %d ; continuing ...\n",
                elapsed,allowed
      );
    }
    else
    {
      LOG_PRINT(LOG_WARNING,"Simulation has been running for the maximum allowed time of %d minutes, stopping simulation immediately !!!\n",allowed);
      
      should_stop=true;
      
      if(LOG_SEVERITY>=LOG_WARNING)
      {
        fprintf(stdout,"Simulation has been running for the maximum allowed time of %d minutes, stopping simulation immediately !!!\n",allowed);
      }
          
    }
        
    return should_stop;
  }
  
  /**
   * @brief If the max allowed running time is almost reached, this is called for performing a proper file backup
   * 
   * This virtual implementation is quite limited, there may be a proper override by derived classes
   * 
   * \todo Finish implementation of the backup features here for an eventual restart one day ?
   */
  virtual void doBackupAndExit()
  {
    LOG_PRINT(LOG_INFO,"Maximum allowed run time was reached, performing data serialization.\n");
    
    //  a lot of other things to backup also ...
    db_backup();
    db_close();
    LOG_FLUSH_ALL();
    
    MPI_Barrier(global_comm);
    
    MPI_CUSTOM_ABORT_MACRO();
  }

private:
  
  bool local_setup_done = false;

  
}; // end of parrep abstract base class


#endif // PARREP_ABSTRACT_HPP_INCLUDED

