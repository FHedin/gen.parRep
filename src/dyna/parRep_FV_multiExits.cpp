/**
 * \file parRep_FV_multiExits.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#include <cstdlib>

#include <algorithm>

#include "rand.hpp"

#include "parRep_FV_multiExits.hpp"

using namespace std;

ParRepFV_multiExits::ParRepFV_multiExits(DATA& _dat,
                                         unique_ptr<ATOM[]>& _at,
                                         unique_ptr<MD_interface>& _md,
                                         unique_ptr<luaInterface>& _luaItf,
                                         const uint32_t _max_allowed_events
                                         ) : ParRepFV(_dat,_at,_md,_luaItf),
                                         allowed_exit_events(_max_allowed_events),
                                         lua_get_serialized_state(luaItf->get_function_get_serialized_state()),
                                         lua_put_serialized_state(luaItf->get_function_put_serialized_state())
{

  if(params.at("afterExit") == "clone")
    type =  multiExitType::MULTI_CLONE;
  else if (params.at("afterExit") == "wait")
    type = multiExitType::MULTI_WAIT;
  
}

/*
 * The function performing FV ParRep
 */
void ParRepFV_multiExits::run()
{
  fprintf(stdout,"\nRunning a Generalized ParRep with Gelman-Rubin statistics and Fleming-Viot particle processes, with multiple exit events (maximum: %u) allowed\n\n",
          allowed_exit_events);

  // open database of states : saved in memory for performance but regularly backed up to a file in case of crash
  db_open();
  
  /*
   * Stage 0 : Equilibrate the system for some steps
   *  and use as reference state at the beginning
   *  only masterRank does this and then broadcasting is used for other ranks
   */
  if(equil_steps>0)
  {
    // only masterRank equilibrates
    if(i_am_master)
    {
      do_equilibration();
    }
    MPI_Barrier(global_comm);
    
    MPIutils::mpi_broadcast_atom_array(dat,at.get(),masterRank,global_comm);
    
    // not required on masterRank : provide to openmm coordinates and velocites from equilibration
    if(!i_am_master)
    {
      md->setCrdsVels(at.get());
    }
      
    MPI_Bcast(&equil_e.ene[0],3,MPI_DOUBLE,masterRank,global_comm);
    
    // update lua interface
    luaItf->set_lua_variable("epot",equil_e.epot());
    luaItf->set_lua_variable("ekin",equil_e.ekin());
    luaItf->set_lua_variable("etot",equil_e.etot());
    luaItf->set_lua_variable("referenceTime",ref_clock_time);
  }
  
  // after equilibration and before starting parrep algorithm it is time to initialise the Lua code defining a state
  lua_state_init();
  
  /*
   * First be sure that the system is within a metastable state, if not perform transient propagation
   */
  bool need_transient_propagation = lua_check_transient();
  if(need_transient_propagation)
  {
    ENERGIES tr_e;
    
    if(i_am_master)
    {
      do_transient_propagation(tr_e);
    }
    MPI_Barrier(global_comm);
    
    // broadcast to the others
    MPI_Bcast(&ref_clock_time,1,MPI_DOUBLE,masterRank,global_comm);
    MPI_Bcast(&tr_e.ene[0],3,MPI_DOUBLE,masterRank,global_comm);
    MPIutils::mpi_broadcast_atom_array(dat,at.get(),masterRank,global_comm);
    
    // update lua interface
    luaItf->set_lua_variable("epot",tr_e.epot());
    luaItf->set_lua_variable("ekin",tr_e.ekin());
    luaItf->set_lua_variable("etot",tr_e.etot());
    luaItf->set_lua_variable("referenceTime",ref_clock_time);
    
    // call state_init again now that we are sure to be within a state
    lua_state_init();
  }
  
  /*
   * NOTE main loop here
   */
  do
  {
    LOG_PRINT(LOG_INFO,"New ParRepFVMultiExit loop iteration : ref_clock_time is %.2lf ps \n",ref_clock_time);
    fprintf(stdout,"\n//---------------------------------------------------------------------------------------------//\n");
    fprintf(stdout,    "New ParRepFVMultiExit loop iteration : ref_clock_time is %.2lf ps \n\n",ref_clock_time);;
    
    fv_local_time = 0.;
    fv_e = ENERGIES();
    
    /*
     * Stage 1 (FV) : equivalent of decorrelation+dephasing done together
     *  using the FV approach and the GR analysis
     */
    do_FlemingViot_procedure();
    // there was already a barrier at the end of do_FlemingViot_procedure() so from here we suppose synchronization of all the replicas
    
    // update lua interface
    luaItf->set_lua_variable("epot",fv_e.epot());
    luaItf->set_lua_variable("ekin",fv_e.ekin());
    luaItf->set_lua_variable("etot",fv_e.etot());
    luaItf->set_lua_variable("referenceTime",ref_clock_time);

    /*
     * stage 2 : Run parallel dynamics  with multiple exit events
     */
    dyna_local_time = 0.;
    dyna_e = ENERGIES();
    
    exit_cfg.clear();
    exit_cfg.resize(0);
    exit_cfg.shrink_to_fit();
      
    // each rank does its independent dynamics
    do_dynamics();
    
    if(exit_cfg.size()>0)
    {
      for(exitState_config& cfg : exit_cfg)
      {
        double esc_time = (double)(cfg.num_workers_when_exiting - 1);
        esc_time *= (cfg.time - t_poll);
        esc_time += ((double)cfg.local_mpi_id_when_exiting)*t_poll ;
        esc_time += cfg.time;
        
        LOG_PRINT(LOG_INFO,"Rank WORLD-%d (LOCAL-%d) calculated escape_time is %.2lf ps\n",cfg.global_mpi_id_when_exiting,cfg.local_mpi_id_when_exiting,esc_time);
        fprintf(stdout,    "Rank WORLD-%d (LOCAL-%d) calculated escape_time is %.2lf ps\n",cfg.global_mpi_id_when_exiting,cfg.local_mpi_id_when_exiting,esc_time);
        
        cfg.escape_time = esc_time;
        
      }
      LOG_FLUSH_ALL();
    }
      
    MPI_Barrier(global_comm);
    
    /*
     * NOTE when doing multi exit parrep should we increase the equivalent MD clock time using the min,
     *      the max or the sum of all the exit times ?
     *      this will have an impact on benchmarking and speedup factors: here we use the largest one
     *      
     * We use the MPI_MAXLOC operation to find simultaneously the max value and its mpi index
    */
    // first find the local max time
    maxloc_type max_escape_time = {numeric_limits<double>::min(),-1};
    vector<exitState_config>::iterator max_iter;
    if(exit_cfg.size() > 0)
    {
      max_iter = max_element(exit_cfg.begin(),exit_cfg.end(),
                             [](const exitState_config& a, const exitState_config& b)->bool{return (a.escape_time < b.escape_time);}
                            );
      max_escape_time.index = mpi_id_in_gcomm;
      max_escape_time.value = max_iter->escape_time;
    }
    
    // then find the max over all the ranks
    MPI_Allreduce(MPI_IN_PLACE,&max_escape_time,1,MPI_DOUBLE_INT,MPI_MAXLOC,global_comm);
    
    LOG_PRINT(LOG_INFO,"Largest escape_time is %.2lf ps, from replica %d\n",max_escape_time.value,max_escape_time.index);
    fprintf(    stdout,"Largest escape_time is %.2lf ps, from replica %d\n",max_escape_time.value,max_escape_time.index);
    

    int32_t sharing_id = max_escape_time.index;
    
    // save the state and the escape time to the database for each replica that has escaped
    if(exit_cfg.size() > 0)
    {
      for(exitState_config& cfg : exit_cfg)
      {
        // update lua interface
        luaItf->set_lua_variable("epot",cfg.e.epot());
        luaItf->set_lua_variable("ekin",cfg.e.ekin());
        luaItf->set_lua_variable("etot",cfg.e.etot());
        
        memcpy(at.get(),cfg.atomic_cfg.get(),dat.natom*sizeof(ATOM));
        memcpy(dat.pbc.data(),cfg.pbc_cfg.data(),sizeof(dat.pbc));
        
        lua_put_serialized_state(cfg.serialized_state_def);
        
        db_insert(true,fv_local_time,cfg.escape_time, // mandatory arguments
                  "exit_order",cfg.exit_rank // starting from there, optional arguments
        );
      }
    }
    
    LOG_PRINT(LOG_INFO,"Will use configuration of replica %d as reference for next iteration\n",sharing_id);
    fprintf(    stdout,"Will use configuration of replica %d as reference for next iteration\n",sharing_id);
    
    // then broadcast the data from the rank corresponding to this max_escape_time.index to all the others
    if(mpi_id_in_gcomm == sharing_id)
    {
      memcpy(at.get(),max_iter->atomic_cfg.get(),dat.natom*sizeof(ATOM));
      memcpy(dat.pbc.data(),max_iter->pbc_cfg.data(),sizeof(dat.pbc));
      dyna_local_time = max_iter->time;
      dyna_e = max_iter->e;
      
      need_transient_propagation = lua_check_transient();
      if(need_transient_propagation)
      {
        do_transient_propagation(dyna_e);
        luaItf->set_lua_variable("referenceTime",ref_clock_time);
      }
    }
    MPI_Barrier(global_comm);
    
    ref_clock_time += max_escape_time.value;

    MPI_Bcast(&dyna_local_time,1,MPI_DOUBLE,sharing_id,global_comm);
    MPI_Bcast(&ref_clock_time,1,MPI_DOUBLE,breakerID,global_comm);
    MPI_Bcast(&(dyna_e.ene[0]),3,MPI_DOUBLE,sharing_id,global_comm);
    MPIutils::mpi_broadcast_atom_array(dat,at.get(),sharing_id,global_comm);
    
    // update lua interface
    luaItf->set_lua_variable("referenceTime",ref_clock_time);
    luaItf->set_lua_variable("epot",dyna_e.epot());
    luaItf->set_lua_variable("ekin",dyna_e.ekin());
    luaItf->set_lua_variable("etot",dyna_e.etot());
    
    // frequent sqlite in-memory db backup to a file
    if( (ref_clock_time-last_db_backup) > db_backup_frequency_ps )
    {
      db_backup();
      last_db_backup = ref_clock_time;
    }
    
    /*
     * prepare a new loop iteration
     */
    
    // we reset the omm time
    md->setSimClockTime(0.);
    // and each rank uses the system of the dyna breaking one
    md->setCrdsVels(at.get());
    
    // now all ranks need to check again what is the current state in order to be ready for next iteration
    lua_state_init();
    
    MPI_Barrier(global_comm);
    
    /*
     * Ready to loop again
     */
    
  } while( ref_clock_time < ((double)dat.nsteps * dat.timestep) );
  // NOTE End of main loop here
  
  // backup in-memory database to file before exiting simulation 
  db_backup();
  db_close();

  MPI_Barrier(global_comm);

}

/*
 * ParRep functions : represent successive stages of the algorithm
 *
 * The following are exactly the same as with standard ParRep_FV so they are not overidden :
 * 
 * void ParRepFV::do_equilibration();
 * void ParRepFV::do_FlemingViot_procedure();
 */

void ParRepFV_multiExits::do_dynamics()
{
  switch(type)
  {
    case MULTI_WAIT:
      do_dynamics_multi_wait();
      break;
      
    case MULTI_CLONE:
      do_dynamics_multi_clone();
      break;
  }
}

void ParRepFV_multiExits::do_dynamics_multi_wait()
{
  LOG_PRINT(LOG_INFO,"Rank %d now entering multiple-exits (with waiting) parallel dynamics stage\n",  mpi_id_in_gcomm);
  fprintf(stdout,    "Rank %d now entering multiple-exits (with waiting) parallel dynamics stage\n\n",mpi_id_in_gcomm);
  
  /* It is necessary to define a new MPI_Comm ensemble of replicas different from "global_comm" : exiting replicas
   *  will be progressively removed from this new MPI_Comm as they exit
   */
  
  int32_t num_workers = mpi_gcomm_size;
  int32_t local_group_id = mpi_id_in_gcomm;
  
  /* construct a vector of 0 to mpi_gcomm_size-1 ids
   *   some ranks will be progressively removed from it as they produce an exit event
   */
  vector<int32_t> worker_ranks(num_workers,0);
  for(int32_t i=0; i<num_workers; i++) worker_ranks[i] = i;
  
#ifdef PARREP_DEBUG_BUILD
  LOG_PRINT(LOG_DEBUG,"Initial dump of parallel workers : \n");
  for(int32_t& i : worker_ranks)
    LOG_PRINT(LOG_DEBUG,"Worker %d\n",i);
  LOG_FLUSH_ALL();
#endif

  // create the working group and the working communicator : contain all the ranks at the beginning, and some will be removed later
  MPI_Group working_group = MPI_GROUP_NULL;
  MPI_Comm  working_comm = MPI_COMM_NULL;
  MPI_Group new_group = MPI_GROUP_NULL;
  MPI_Comm  new_comm = MPI_COMM_NULL;
  MPI_Group_incl(global_group,num_workers,worker_ranks.data(),&working_group);
  MPI_Comm_create(global_comm, working_group, &working_comm);
  
  // initially num_workers == mpi_gcomm_size and local_group_id == mpi_id_in_gcomm but check to be sure
  MPI_Comm_rank(working_comm, &local_group_id);
  assert(local_group_id==mpi_id_in_gcomm);
  
  MPI_Comm_size(working_comm, &num_workers);
  assert(num_workers==mpi_gcomm_size);

  uint32_t num_exits = 0, l_num_exits = 0;
  
  uint32_t i_exit = 0;
  
  do
  {
    // do some steps
    md->doNsteps(dynamics_check);
    dyna_local_time += t_poll;
    
    md->getSimData(nullptr,&dyna_e,nullptr,at.get());
    luaItf->set_lua_variable("epot",dyna_e.epot());
    luaItf->set_lua_variable("ekin",dyna_e.ekin());
    luaItf->set_lua_variable("etot",dyna_e.etot());
    
    left_state = lua_check_state();
    
    /*
     * We check if this rank left the current state
     */
    if(left_state)
    {
      LOG_PRINT(LOG_INFO,"Rank WORLD-%d (LOCAL-%d) escaped the state after %.2lf ps\n",mpi_id_in_gcomm,local_group_id,dyna_local_time);
      fprintf(stdout,    "Rank WORLD-%d (LOCAL-%d) escaped the state after %.2lf ps\n\n",mpi_id_in_gcomm,local_group_id,dyna_local_time);

      exitState_config cfg;
      
      cfg.global_mpi_id_when_exiting = mpi_id_in_gcomm;
      cfg.local_mpi_id_when_exiting = local_group_id;
      cfg.num_workers_when_exiting = num_workers;
      cfg.exit_rank = num_exits + 1;
      
      cfg.atomic_cfg = unique_ptr<ATOM[]>(new ATOM[dat.natom]);
      memcpy(cfg.atomic_cfg.get(),at.get(),dat.natom*sizeof(ATOM));
      
      memcpy(cfg.pbc_cfg.data(),dat.pbc.data(),sizeof(dat.pbc));
      
      cfg.time = dyna_local_time;
      cfg.e = dyna_e;
      
      cfg.serialized_state_def = lua_get_serialized_state();
      
      exit_cfg.push_back(std::move(cfg));
      
      i_exit = 1;
    }
    else
    {
      LOG_PRINT(LOG_INFO,"Still within the same state after %.2lf ps...\n",dyna_local_time);
      if(LOG_SEVERITY == LOG_DEBUG)
        fprintf(stdout,  "Still within the same state after %.2lf ps...\n\n",dyna_local_time);
    }
    
    // first check how many reps left the state during the last t_poll windows; although quite rare multiple exits should be considered
    l_num_exits = 0;
    MPI_Allreduce(&i_exit,&l_num_exits,1,MPI_UINT32_T,MPI_SUM,working_comm);
    
    if(l_num_exits > 0)
    {
      
      LOG_PRINT(LOG_INFO,"It appears that %u ranks escaped the state after %.2lf ps\n",l_num_exits,dyna_local_time);
      fprintf(stdout,    "It appears that %u ranks escaped the state after %.2lf ps\n",l_num_exits,dyna_local_time);
      
      unique_ptr<int32_t[]> exiting_ranks(new int32_t[num_workers]);
      
      // if this ranks exits it sends a -1, if not its local id
      int32_t send_flag = (i_exit==1)?-1:local_group_id;
      MPI_Allgather(&send_flag,1,MPI_INT32_T,
                    exiting_ranks.get(),1,MPI_INT32_T,
                    working_comm);
      
#ifdef PARREP_DEBUG_BUILD
      LOG_PRINT(LOG_DEBUG,"Dump of exiting_ranks after %u exits : \n",l_num_exits);
      for(int32_t i=0; i<num_workers; i++)
        LOG_PRINT(LOG_DEBUG,"%d\n",exiting_ranks[i]);
      LOG_FLUSH_ALL();
#endif

      worker_ranks = vector<int32_t>();

      for(int32_t r=0; r<num_workers; r++)
      {
        if(exiting_ranks[r]!=-1)
        {
          worker_ranks.push_back(exiting_ranks[r]);
        }
      }
      
      num_workers = (int32_t) worker_ranks.size();
      
#ifdef PARREP_DEBUG_BUILD
      LOG_PRINT(LOG_DEBUG,"Dump of parallel workers after %u exits : total size is %d\n",l_num_exits,num_workers);
      for(int32_t& i : worker_ranks)
        LOG_PRINT(LOG_DEBUG,"Worker %d\n",i);
      LOG_FLUSH_ALL();
#endif

      num_exits += l_num_exits;
      
      // if not destroy and recreate the mpi groups and communicators, this time only with the remaining workers
      MPI_Group_incl(working_group,num_workers,worker_ranks.data(),&new_group);
      MPI_Group_free(&working_group);
      working_group = new_group;

      MPI_Comm_create(working_comm, working_group, &new_comm);
      MPI_Comm_free(&working_comm);
      working_comm = new_comm;

    }

    /*
     * if the other ranks continue but this one exited, it leaves
     * this piece of the code now and will wait for the others to finish in the run method
     */
    if(left_state)
    {
      break;
    }
    
    // if at least allowed_exit_events have been observed it is time to exit this infinite while loop
    if(num_exits >= allowed_exit_events)
    {
      LOG_PRINT(LOG_INFO,"The number of exits (%u) reached the maximum allowed limit (%u), exiting parallel phase...\n",num_exits,allowed_exit_events);
      LOG_FLUSH_ALL();
      
      break;
    }
    
    if(l_num_exits > 0)
    {
      // reset local time for those who didn't exited, and run the race again
      dyna_local_time = 0.0;
      md->setSimClockTime(0.);
      
      // get new local id for this configuration: important as the parrep exit time depends on the exiting rank AND total number of ranks when an exit is observed
      MPI_Comm_rank(working_comm, &local_group_id);
      MPI_Comm_size(working_comm, &num_workers); // redundant, but to be sure
    }

    // check running time and depending on return value perform a backup before exiting
    bool need_stop = check_running_time();
    if(need_stop) doBackupAndExit();
    
    LOG_PRINT(LOG_DEBUG,"This rank WORLD-%d (LOCAL-%d) will continue his job now...\n",mpi_id_in_gcomm,local_group_id);
    
  } while(true);

}

void ParRepFV_multiExits::do_dynamics_multi_clone()
{
  LOG_PRINT(LOG_INFO,"Rank %d now entering multiple-exits (with cloning) parallel dynamics stage\n",  mpi_id_in_gcomm);
  fprintf(stdout,    "Rank %d now entering multiple-exits (with cloning) parallel dynamics stage\n\n",mpi_id_in_gcomm);

  uint32_t num_clones = 0;
  uint32_t num_exits = 0;
  
  uint32_t l_num_exits = 0;
  uint32_t i_exit = 0;
  
  bool i_need_cloning = false;
  
  MPI_Win at_window = MPI_WIN_NULL;
  MPI_Win pbc_window = MPI_WIN_NULL;
  MPI_Win branching_window = MPI_WIN_NULL;
  
  // set the at[] as a window of memory accessible by other ranks: one sided communications
  MPI_Win_create(at.get(),dat.natom*sizeof(ATOM),1,      //void *base, MPI_Aint size, int disp_unit
                 rma_info,global_comm,&at_window);  // MPI_Info info, MPI_Comm comm, MPI_Win *win
  
  // same for the pbc
  MPI_Win_create(dat.pbc.data(),sizeof(dat.pbc),1,
                 rma_info,global_comm,&pbc_window);
  
  MPI_Win_create(&i_need_cloning,sizeof(bool),sizeof(bool),
                 rma_info,global_comm,&branching_window);
  
  do
  {
    // do some steps
    md->doNsteps(dynamics_check);
    dyna_local_time += t_poll;
    
    md->getSimData(nullptr,&dyna_e,nullptr,at.get());
    luaItf->set_lua_variable("epot",dyna_e.epot());
    luaItf->set_lua_variable("ekin",dyna_e.ekin());
    luaItf->set_lua_variable("etot",dyna_e.etot());

    /*
     * We check if this rank left the reference state
     */
    left_state = lua_check_state();
    if(left_state)
    {
      LOG_PRINT(LOG_INFO,"Rank %d escaped the state after %.2lf ps\n",mpi_id_in_gcomm,dyna_local_time);
      fprintf(stdout,"Rank %d escaped the state after %.2lf ps\n\n",mpi_id_in_gcomm,dyna_local_time);
      
      i_exit = 1;
      i_need_cloning = true;
    }
    else
    {
      LOG_PRINT(LOG_INFO,"Still within the same state after %.2lf ps...\n",dyna_local_time);
      if(LOG_SEVERITY == LOG_DEBUG)
        fprintf(stdout,  "Still within the same state after %.2lf ps...\n\n",dyna_local_time);
    }
    
    // first check how many reps left the state during the last t_poll windows; although quiet rare multiple exits should be considered
    l_num_exits = 0;
    MPI_Allreduce(&i_exit,&l_num_exits,1,MPI_UINT32_T,MPI_SUM,global_comm);
    
    num_clones += l_num_exits;
    num_exits  += l_num_exits;
    
    if(l_num_exits > 0)
    {
      LOG_PRINT(LOG_INFO,"It appears that %u ranks escaped the state after %.2lf ps\n",l_num_exits,dyna_local_time);
      fprintf(stdout,    "It appears that %u ranks escaped the state after %.2lf ps\n",l_num_exits,dyna_local_time);
    }
    
    if(left_state)
    {
      exitState_config cfg;
      
      // if it left the state, save the exit config of this rank (global and local are the same with this cloning variant); the num_workers is also constant
      cfg.global_mpi_id_when_exiting = mpi_id_in_gcomm;
      cfg.local_mpi_id_when_exiting  = mpi_id_in_gcomm;
      cfg.num_workers_when_exiting   = mpi_gcomm_size;
      cfg.exit_rank                  = num_clones + 1;
      
      cfg.atomic_cfg = unique_ptr<ATOM[]>(new ATOM[dat.natom]);
      memcpy(cfg.atomic_cfg.get(),at.get(),dat.natom*sizeof(ATOM));
      
      memcpy(cfg.pbc_cfg.data(),dat.pbc.data(),sizeof(dat.pbc));
      
      cfg.time = dyna_local_time;
      cfg.e = dyna_e;
      
      cfg.serialized_state_def = lua_get_serialized_state();
      
      exit_cfg.push_back(std::move(cfg));
      
      // clone one of the replica who didn't exit and continue to run
      // first we find a proper candidate
      int32_t candidate = -1;
      do
      {
        candidate = get_int32_min_max(0,mpi_gcomm_size-1);
        
        if(candidate==mpi_id_in_gcomm)
          continue;
        
        bool candidate_also_exited = false;
        
        MPI_Win_lock(MPI_LOCK_SHARED,candidate,0,branching_window);
        MPI_Get(&candidate_also_exited,1,     // void *origin_addr, int origin_count,
                MPI_CXX_BOOL,candidate,         // MPI_Datatype origin_datatype, int target_rank,
                0,1,                          // MPI_Aint target_disp, int target_count,
                MPI_CXX_BOOL,branching_window); // MPI_Datatype target_datatype, MPI_Win win
        MPI_Win_unlock(candidate,branching_window);
        
        if(!candidate_also_exited)
          break;
        
      }while(true);
      
      LOG_PRINT(LOG_DEBUG,"Rank %d left the state, will clone from candidate %d...\n",mpi_id_in_gcomm,candidate);
      
      // put a lock on at_window and pbc_window on rank candidate so that it is not written by candidate
      //  all other ranks are allowed to read from candidate because we use MPI_LOCK_SHARED (i.e. several ranks can query the same candidate)
      MPI_Win_lock(MPI_LOCK_SHARED,candidate,0,at_window);
      MPI_Win_lock(MPI_LOCK_SHARED,candidate,0,pbc_window);
      
      MPI_Get(at.get(),dat.natom*sizeof(ATOM), // void *origin_addr, int origin_count,
              MPI_BYTE,candidate,              // MPI_Datatype origin_datatype, int target_rank,
              0,dat.natom*sizeof(ATOM),        // MPI_Aint target_disp, int target_count,
              MPI_BYTE,at_window);             // MPI_Datatype target_datatype, MPI_Win win
      
      MPI_Get(dat.pbc.data(),sizeof(dat.pbc), // void *origin_addr, int origin_count,
              MPI_BYTE,candidate,                    // MPI_Datatype origin_datatype, int target_rank,   
              0,sizeof(dat.pbc),              // MPI_Aint target_disp, int target_count,
              MPI_BYTE,pbc_window);                  // MPI_Datatype target_datatype, MPI_Win win
      
      // remove lock and proceed
      MPI_Win_unlock(candidate,at_window);
      MPI_Win_unlock(candidate,pbc_window);
      
      md->setCrdsVels(at.get());
      
      md->getSimData(nullptr,&dyna_e,nullptr,nullptr);
      luaItf->set_lua_variable("epot",dyna_e.epot());
      luaItf->set_lua_variable("ekin",dyna_e.ekin());
      luaItf->set_lua_variable("etot",dyna_e.etot());
      
      lua_state_init();
      md->setSimClockTime(0.);
      
      LOG_PRINT(LOG_DEBUG,"State has been successfully cloned from candidate %d\n",candidate);
      
      // reset exit flag
      i_exit = 0;
      i_need_cloning = false;

    }
    
    if(l_num_exits > 0)
    {
      // reset local time for those who didn't exited, and run the race again
      dyna_local_time = 0.0;
      md->setSimClockTime(0.);
    }
    
    // required
    MPI_Barrier(global_comm);
    
    // if at least allowed_exit_events have been observed it is time to exit this infinite while loop
    if(num_exits >= allowed_exit_events)
    {
      LOG_PRINT(LOG_INFO,"The number of exits (%u) reached the maximum allowed limit (%u), exiting parallel phase...\n",num_exits,allowed_exit_events);
      fprintf(stdout,    "The number of exits (%u) reached the maximum allowed limit (%u), exiting parallel phase...\n",num_exits,allowed_exit_events);
      LOG_FLUSH_ALL();
      
      break;
    }
    
    // check running time and depending on return value perform a backup before exiting
    bool need_stop = check_running_time();
    if(need_stop) doBackupAndExit();
    
  }while(true);
  
  MPI_Win_free(&at_window);
  MPI_Win_free(&pbc_window);
  MPI_Win_free(&branching_window);
  
}
