/**
 * \file parRep_FV.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#include <utility>
#include <limits>
#include <array>

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "logger.hpp"
#include "rand.hpp"

#include "parRep_FV.hpp"

using namespace std;

ParRepFV::ParRepFV(DATA& _dat,
                   unique_ptr<ATOM[]>& _at,
                   unique_ptr<MD_interface>& _md,
                   unique_ptr<luaInterface>& _luaItf,
                   bool ignore_parrepabstract_setup
                  ) : ParRepAbstract(_dat,_at,_md,_luaItf,ignore_parrepabstract_setup),
                      gr_functions(luaItf->get_gr_functions()),
                      gr_names(luaItf->get_gr_names())
{
  equil_steps    = (uint32_t) stoul(params.at("equilibrationSteps"));
  gr_check       = (uint32_t) stoul(params.at("checkGR"));
  fv_check       = (uint32_t) stoul(params.at("checkFV"));
  discard_first_N_grObs = (uint32_t) stoul(params.at("minAccumulatedObs"));
  dynamics_check = (uint32_t) stoul(params.at("checkDynamics"));
  grTol          = stod(params.at("GRtol"));
  
  // and also get database parameters
  db_backup_frequency_ps = stod(params.at("dbFreq"));
  
  // calculate polling time from input params
  t_poll = ((double)dynamics_check)*dat.timestep ;
  
  if(!ignore_parrepabstract_setup)
  {
    // assign a Fleming-Viot role to this replica, baster on MPI rank
    my_FV_role = (i_am_master) ? REF_WALKER : FV_WALKER;
  }
}

///////////////////////////////////////////////////
/*
 * The function performing FV ParRep
 */
///////////////////////////////////////////////////

void ParRepFV::run()
{
  fprintf(stdout,"\nRunning a Generalized ParRep with Gelman-Rubin statistics and Fleming-Viot particle processes.\n"
                 "Role of this replica is: %s\n",role_string.at(my_FV_role).c_str());
  
  // open database of states
  if(i_am_master)
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
    // not required on masterRank : provide to openmm coordinates and velocities from equilibration
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
  
  MPI_Barrier(global_comm);
  
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
    MPI_Bcast(tr_e.ene.data(),3,MPI_DOUBLE,masterRank,global_comm);
    MPIutils::mpi_broadcast_atom_array(dat,at.get(),masterRank,global_comm);
    
    // update lua interface
    luaItf->set_lua_variable("epot",tr_e.epot());
    luaItf->set_lua_variable("ekin",tr_e.ekin());
    luaItf->set_lua_variable("etot",tr_e.etot());
    luaItf->set_lua_variable("referenceTime",ref_clock_time);
    
    // call state_init again now that we are sure to be within a state
    lua_state_init();
  }
  
  MPI_Barrier(global_comm);
  
  /*
   * NOTE main loop here
   */
  do
  {
    LOG_PRINT(LOG_INFO,"New ParRepFV loop iteration : ref_clock_time is %.2lf ps \n",  ref_clock_time);
    fprintf(stdout,"\n//---------------------------------------------------------------------------------------------//\n");
    fprintf(stdout,    "New ParRepFV loop iteration : ref_clock_time is %.2lf ps \n\n",ref_clock_time);
    
    fv_local_time = 0.;
    fv_e = ENERGIES();
    
    /*
     * Stage 1 (FV) : equivalent of decorrelation + dephasing done together
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
     * stage 2 : Run parallel dynamics 
     */
    dyna_local_time = 0.;
    dyna_e = ENERGIES();
    breakerID = numeric_limits<uint32_t>::max();
    dyna_cycles_done = 0;
      
    // each rank does its independent dynamics
    do_dynamics();
    // there was already a barrier at the end of do_dynamics() so from here we suppose synchronization of all the replicas
    
    /*
     * find and send the id of the rank that escaped the state to all others
     * MPI_MIN corresponds to a min of the argmin
     */
    MPI_Allreduce(MPI_IN_PLACE,&breakerID,1,MPI_UINT32_T,
                  MPI_MIN,global_comm);
    
    // transient propagation if required
    if(mpi_id_in_gcomm==(int32_t)breakerID)
    {
      need_transient_propagation = lua_check_transient();
      if(need_transient_propagation)
      {
        do_transient_propagation(dyna_e);
        luaItf->set_lua_variable("referenceTime",ref_clock_time);
      }
    }
    MPI_Barrier(global_comm);
    
    /*
     * broacast data from breaking replica to the others
     */
    MPI_Bcast(&dyna_local_time,1,MPI_DOUBLE,breakerID,global_comm);
    MPI_Bcast(&dyna_cycles_done,1,MPI_UINT32_T,breakerID,global_comm);
    MPI_Bcast(&ref_clock_time,1,MPI_DOUBLE,breakerID,global_comm);
    MPI_Bcast(dyna_e.ene.data(),3,MPI_DOUBLE,breakerID,global_comm);
    
    // coordinates and velocities are broadcasted using a non blocking ibroadcast so that we can do something else in the mean time
    MPIutils::mpi_broadcast_atom_array(dat,at.get(),breakerID,global_comm);
    
    /*
     * the ref_clock_time is updated using the escape time
     * This formula comes from the article "parRep for markov chains" (Aristoff, Lelièvre, Simpson)
     */
    double escape_time =  (double)(mpi_gcomm_size-1);
    escape_time *= ((double)dyna_cycles_done)*t_poll;
    escape_time += ((double)breakerID)*t_poll ;
    escape_time += dyna_local_time;
    
    ref_clock_time += escape_time;
    
    LOG_PRINT(LOG_INFO,"Parallel dynamics stopped because rank %u escaped (after %.2lf ps).\n",  breakerID,dyna_local_time);
    fprintf(    stdout,"Parallel dynamics stopped because rank %u escaped (after %.2lf ps).\n\n",breakerID,dyna_local_time);
    
    LOG_PRINT(LOG_INFO,"ParRep_FV escape_time is %.2lf ps\n",escape_time);
    fprintf(    stdout,"ParRep_FV escape_time is %.2lf ps\n\n",escape_time);
    
    // update lua interface
    luaItf->set_lua_variable("epot",dyna_e.epot());
    luaItf->set_lua_variable("ekin",dyna_e.ekin());
    luaItf->set_lua_variable("etot",dyna_e.etot());
    luaItf->set_lua_variable("referenceTime",ref_clock_time);

    // save the state and the escape time to the database
    if(i_am_master)
    {
      db_insert(true,fv_local_time,escape_time);
    }
    
    // if necessary db backup is performed
    if(i_am_master && ((ref_clock_time-last_db_backup) > db_backup_frequency_ps) )
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
    
    // now all ranks need to re-initialize the state identification code on Lua's side
    lua_state_init();
    
    /*
     * Ready to loop again
     */
    MPI_Barrier(global_comm);
    
  }while( ref_clock_time < ((double)dat.nsteps * dat.timestep) );
  // NOTE End of main loop here
  
  if(i_am_master)
  {
    // backup database to file if required before exiting simulation 
    db_backup();
    db_close();
  }
  
  MPI_Barrier(global_comm);
      
} // end of function run

///////////////////////////////////////////////////
/*
 * Parrep functions : represent successive stages of the algorithm
 */
///////////////////////////////////////////////////

void ParRepFV::do_FlemingViot_procedure()
{
  LOG_PRINT(LOG_INFO,"Rank %d performing Fleming-Viot procedure\n",  mpi_id_in_gcomm);
  fprintf(stdout,    "Rank %d performing Fleming-Viot procedure\n\n",mpi_id_in_gcomm);
  
  /**
   *  the C++ interface for collecting Gelman-Rubin statistics
   *  For the moment only rank masterRank will manage it.
   *  The others regularly send their data (stored in gr_observations) to rank masterRank
   */
  std::unique_ptr<GelmanRubinAnalysis> gr = nullptr;
  
  // The ref walker is in charge of allocating the GelmanRubinAnalysis object, it will also check the convergence
  if(i_am_master)
  {
    gr = unique_ptr<GelmanRubinAnalysis>(new GelmanRubinAnalysis((uint32_t)mpi_gcomm_size,discard_first_N_grObs));
    // register grObservables parsed from lua file
    for(const pair<GR_function_name,GR_function>& p : gr_functions)
    {
      gr->registerObservable(Observable(p.first),grTol);
    }
  }
  
  // number of observables that the workers send to the ref walker every fv_check steps
  const uint32_t numObs = fv_check/gr_check;
  
  map<GR_function_name,double*> recvObsMap;
  if(my_FV_role == REF_WALKER)
  {
    for(const GR_function_name& n: gr_names)
      recvObsMap[n] = new double[mpi_gcomm_size*numObs];
  }
  
  // stores requests used when gathering in non blocking mode on ref walker the observables sent by the workers
  vector<MPI_Request> igather_obs_reqs(gr_names.size(),MPI_REQUEST_NULL);
  
  // if this rank is requiring branching
  bool i_need_branching = false;

  // declaration and create of sattic memory windows for one sided MPI comms
  MPI_Win at_window        = MPI_WIN_NULL;
  MPI_Win pbc_window       = MPI_WIN_NULL;
  MPI_Win branching_window = MPI_WIN_NULL;

  // set the at[] as a window of memory accessible by other ranks: one sided communications
  MPI_Win_create(at.get(),dat.natom*sizeof(ATOM),1,      //void *base, MPI_Aint size, int disp_unit
                 rma_info,global_comm,&at_window);  // MPI_Info info, MPI_Comm comm, MPI_Win *win
  // same for the pbc
  MPI_Win_create(dat.pbc.data(),sizeof(dat.pbc),1,
                 rma_info,global_comm,&pbc_window);
  // and also for the i_need_branching variable
  MPI_Win_create(&i_need_branching,sizeof(bool),sizeof(bool),
                 rma_info,global_comm,&branching_window);

  // we need a window for exchanging Gelman-Rubin data but its size is varying so we create it here using a dynamic MPI window
  map<GR_function_name,MPI_Win> gr_windows;
  for(const GR_function_name& n: gr_names)
  {
    gr_windows[n] = MPI_WIN_NULL;
    MPI_Win_create_dynamic(rma_info,
                           global_comm,
                           &gr_windows[n]);
  }
  
  /*
   * Each mpi rank uses a map for storing observations of each observable during simulation
   * 
   * For each observable with a given GR_function_name, a vector of double stores the associated observations
   * 
   * From time to time this data is sent to masterRank which owns the GelmanRubinAnalysis object, and the masterRank checks convergence
   */
  map<GR_function_name,double*> gr_observations;
  
  /*
   * In order to use MPI_Win_create_dynamic we also need a record of the MPI_Aint 
   * address at which the vector<double>() are stored, and a static window for allowing RMA access to each of them
   */
  map<GR_function_name,MPI_Aint> gr_addr;
  map<GR_function_name,MPI_Win>  gr_windows_addr;
  
  for(const GR_function_name& n: gr_names)
  {
    //gr_observations[n] = vector<double>();
    gr_observations[n] = new double[numObs];
    gr_addr[n] = 0;
    gr_windows_addr[n] = MPI_WIN_NULL;
    MPI_Win_create(&(gr_addr[n]),sizeof(MPI_Aint),sizeof(MPI_Aint),
                   rma_info,global_comm,&(gr_windows_addr[n]));
  }

  md->getSimData(nullptr,&fv_e,nullptr,at.get());
  luaItf->set_lua_variable("epot",fv_e.epot());
  luaItf->set_lua_variable("ekin",fv_e.ekin());
  luaItf->set_lua_variable("etot",fv_e.etot());
  luaItf->set_lua_variable("referenceTime",ref_clock_time);

  ////////////////////////////////////////////////////////////////////////////
  /*
   * loop until there is a converged distribution of FV walkers (i.e. ranks) within the state
   */
  uint32_t steps_since_last_check_state = 0;
  uint32_t steps_since_beginning = 0;
  uint32_t obsIndex = 0;
  
  bool converged = false;
  bool reset_loop = false;
  
  MPI_Win converged_window = MPI_WIN_NULL;
  MPI_Win reset_window = MPI_WIN_NULL;
  
  MPI_Win_create(&converged,sizeof(bool),sizeof(bool),
                 rma_info,global_comm,&converged_window);
  MPI_Win_create(&reset_loop,sizeof(bool),sizeof(bool),
                 rma_info,global_comm,&reset_window);
  
  /*
   * we will use when possible a non blocking barrier (ibarrier), and before checking the request associated to it,
   * we will try to perform independent computations in order to hide latency of the ibarrier
   */
  MPI_Request barrier_req = MPI_REQUEST_NULL;
  
  while(true)
  {

    md->doNsteps(gr_check);

    // wait until data sent in the igather (observables) arrived to the master node
    MPI_Waitall(gr_names.size(),igather_obs_reqs.data(),MPI_STATUSES_IGNORE);

    if(steps_since_last_check_state == 0)
    {
      obsIndex = 0;
    }
    
    /*
     * wait for the last ibarrier at the end of the previous loop ieration completed
     * if first iteration or if there was a loop iteration without barrier this has no effect
     */
    MPI_Wait(&barrier_req,MPI_STATUS_IGNORE);

    if(converged)
      break;
    
    steps_since_last_check_state += gr_check;
    steps_since_beginning += gr_check;
    
    //increment local time and continue
    fv_local_time += (double) gr_check * dat.timestep;

    md->getSimData(nullptr,&fv_e,nullptr,at.get());
    luaItf->set_lua_variable("epot",fv_e.epot());
    luaItf->set_lua_variable("ekin",fv_e.ekin());
    luaItf->set_lua_variable("etot",fv_e.etot());
    
    for(const pair<GR_function_name,GR_function>& p : gr_functions)
    {
      // first retrieve the type of observable (name), the pointer to the lua function (f) corresponding to it, and the vector of observations (v) corresponding to name 
      const GR_function_name& name = p.first;
      const GR_function&         f = p.second;

      gr_observations[name][obsIndex] = f();
    }
    obsIndex += 1;
    
    // do the following F-V branching procedure less often (every fv_check steps)
    //  than the GR accumulation phase (gr_check)
    if( steps_since_last_check_state%fv_check != 0 )
      continue;
    
    // check running time and depending on return value perform a backup before exiting
    if(check_running_time())
      doBackupAndExit();

    left_state = lua_check_state();
    i_need_branching = left_state;
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //  if REF_WALKER exiting: it is required to tell the FV_WORKER s that we have to reset
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    
    if(i_am_master && left_state)
    {
      LOG_PRINT(LOG_INFO,"Rank %d (Ref. Walker) exited the ref. state...\n",  mpi_id_in_gcomm);
      fprintf(stdout,    "Rank %d (Ref. Walker) exited the ref. state...\n\n",mpi_id_in_gcomm);
      
      db_insert(false,0.0,fv_local_time);
      
      // set the reset_loop flag to true for all other replicas so that they stop soon
      reset_loop = true;
      
      MPI_Win_lock_all(MPI_MODE_NOCHECK,reset_window);
      for(int32_t i=1; i<mpi_gcomm_size; i++)
      {
        MPI_Put(&reset_loop, //const void *origin_addr
                1,MPI_CXX_BOOL, //int origin_count, MPI_Datatype origin_datatype,
                i,0,  // int target_rank, MPI_Aint target_disp,
                1,MPI_CXX_BOOL, // int target_count, MPI_Datatype target_datatype,
                reset_window); // MPI_Win win
      }
      MPI_Win_unlock_all(reset_window);
      
      // perform transient propagation if required
      bool need_transient_propagation = lua_check_transient();
      if(need_transient_propagation)
      {
        do_transient_propagation(fv_e);
        luaItf->set_lua_variable("referenceTime",ref_clock_time);
        luaItf->set_lua_variable("epot",fv_e.epot());
        luaItf->set_lua_variable("ekin",fv_e.ekin());
        luaItf->set_lua_variable("etot",fv_e.etot());
      }
    }

    /*
     * we use a non-blocking barrier because we execute between the barrier and the corresponding MPI_Wait routine
     * some independent computations
     */
    MPI_Ibarrier(global_comm,&barrier_req);
    
    steps_since_last_check_state = 0;
    
    // check running time and depending on return value perform a backup before exiting
    if(check_running_time())
      doBackupAndExit();
    
    /*
     * attach the windows to the observable values
     * NOTE useless if reset_loop is triggered below
     * But this is supposed to be something really rare, so we put this code here in order to hibe the latency of the barrier.
     * if reset_loop is trigerred, this step is undone immediately.
     */
    for(const GR_function_name& n: gr_names)
    {
      MPI_Win& w = gr_windows[n];
      MPI_Aint& addr = gr_addr[n];
      MPI_Win_attach(w, gr_observations[n], obsIndex*sizeof(double));
      MPI_Get_address(gr_observations[n],&addr);
    }
    
    // we did our best to do something else, but now we have no choice and we should wait for the ibarrier completion ...
    MPI_Wait(&barrier_req,MPI_STATUS_IGNORE);
    
    // executed by all ranks when a reset is required
    if(reset_loop)
    {
      MPIutils::mpi_broadcast_atom_array(dat,at.get(),masterRank,global_comm);
      MPI_Bcast(fv_e.ene.data(),3,MPI_DOUBLE,masterRank,global_comm);
      MPI_Bcast(&ref_clock_time,1,MPI_DOUBLE,masterRank,global_comm);
      
      // before waiting for the ibcast to complete, reset other things here to hide latency 
      
      converged = false;
      reset_loop = false;
      
      steps_since_last_check_state = 0;
      steps_since_beginning = 0;
      
      // enforce requests to be null
      // NOTE already the case if everything done properly
      if(barrier_req != MPI_REQUEST_NULL)
        MPI_Request_free(&barrier_req);

      // undo the useless memory attachment
      for(const GR_function_name& n: gr_names)
      {
        MPI_Win& w = gr_windows[n];

        MPI_Win_detach(w, gr_observations[n]);
        
        gr_addr[n] = 0;
        obsIndex = 0;
      }
      
      for(size_t i=0; i<gr_names.size(); i++)
      {
        if(igather_obs_reqs[i] != MPI_REQUEST_NULL)
          MPI_Request_free(&igather_obs_reqs[i]);
      }
      
      md->setSimClockTime(0.);
      
      luaItf->set_lua_variable("epot",fv_e.epot());
      luaItf->set_lua_variable("ekin",fv_e.ekin());
      luaItf->set_lua_variable("etot",fv_e.etot());
      
      // update reference clock time
      ref_clock_time += fv_local_time;
      luaItf->set_lua_variable("referenceTime",ref_clock_time);
      
      if(!i_am_master)
      {
        md->setCrdsVels(at.get());
      }
      else /* REF WALKER reset chains of accumulated G-R Observables*/
      {
        gr->reset_all_chains();
      }
      
      // now all ranks need to check again what is the current state in order to be ready for next iteration
      lua_state_init();

      LOG_PRINT(LOG_INFO,"Reset of ref_walker and fv_walker s, starting again the FV loop...\n");
      
      continue; // go back to beginning of the do-while loop
    }
    
    /*
     *  F-V branching procedure will start here, for ranks 1 to N-1 (as rank 0 is the ref walker it is not concerned by the branching)
     *  We should randomly choose one of the other ranks, and ask for its values for :
     *    + coordinates
     *    + velocities
     *    + PBCs
     *    + GR observables accumulation history
     * 
     *  Then update local ones with those values, and run again
     */

    // if required, find a valid candidate for branching
    if( (my_FV_role == FV_WALKER) && i_need_branching)
    {
      // when a rank will branch it will send a communication request to another node (a 'candidate')
      int32_t candidate = -1;
      
      /*
       * A good candidate is :
       *  + another chain ...
       *  + ... not requiring branching
       */
      do
      {
        candidate = get_int32_min_max(0,mpi_gcomm_size-1);
        
        if( (candidate == mpi_id_in_gcomm) || (candidate == masterRank) )
          continue;
        
        bool candidate_also_branching = false;
        
        MPI_Win_lock(MPI_LOCK_SHARED,candidate,0,branching_window);
        MPI_Get(&candidate_also_branching,1, // void *origin_addr, int origin_count,
                MPI_CXX_BOOL,candidate,        // MPI_Datatype origin_datatype, int target_rank,
                0,1,                         // MPI_Aint target_disp, int target_count,
                MPI_CXX_BOOL,branching_window);// MPI_Datatype target_datatype, MPI_Win win
        MPI_Win_unlock(candidate,branching_window);
        
        if(!candidate_also_branching)
          break;
        
      }while(true);

      LOG_PRINT(LOG_DEBUG,"Rank %d left the ref state, initiating a F-V branching from candidate %d...\n",mpi_id_in_gcomm,candidate);
      
      // other ranks are allowed to read from the same candidate ( use of MPI_LOCK_SHARED )
      MPI_Win_lock(MPI_LOCK_SHARED,candidate,0,at_window);
      MPI_Win_lock(MPI_LOCK_SHARED,candidate,0,pbc_window);
      
      MPI_Get(at.get(),dat.natom*sizeof(ATOM), // void *origin_addr, int origin_count,
              MPI_BYTE,candidate,              // MPI_Datatype origin_datatype, int target_rank,
              0,dat.natom*sizeof(ATOM),        // MPI_Aint target_disp, int target_count,
              MPI_BYTE,at_window);             // MPI_Datatype target_datatype, MPI_Win win
      
      MPI_Get(dat.pbc.data(),sizeof(dat.pbc),  // void *origin_addr, int origin_count,
              MPI_BYTE,candidate,              // MPI_Datatype origin_datatype, int target_rank,   
              0,sizeof(dat.pbc),               // MPI_Aint target_disp, int target_count,
              MPI_BYTE,pbc_window);            // MPI_Datatype target_datatype, MPI_Win win

      // remove lock and proceed
      MPI_Win_unlock(candidate,at_window);
      MPI_Win_unlock(candidate,pbc_window);

      // do the same for gr data
      for(const GR_function_name& n: gr_names)
      {
        // step 1 : retrieve mem address for using RMA op on a dynamic window below
        MPI_Aint addr = 0;
        MPI_Win& gr_win_addr = gr_windows_addr[n];
        MPI_Win_lock(MPI_LOCK_SHARED,candidate,0,gr_win_addr);
        MPI_Get(&addr,1,MPI_AINT,candidate,
                0,1,MPI_AINT,gr_win_addr);
        MPI_Win_unlock(candidate,gr_win_addr);
        
        // step 2 : data exchange
        MPI_Win& gr_win = gr_windows[n];
        
        MPI_Win_lock(MPI_LOCK_SHARED,candidate,0,gr_win);
        MPI_Get(gr_observations[n],obsIndex,// void *origin_addr, int origin_count,
                MPI_DOUBLE,candidate, // MPI_Datatype origin_datatype, int target_rank,
                addr,obsIndex,         // MPI_Aint target_disp, int target_count,
                MPI_DOUBLE,gr_win);   // MPI_Datatype target_datatype, MPI_Win win
        MPI_Win_unlock(candidate,gr_win);
      }

      LOG_PRINT(LOG_DEBUG,"Rank %d properly branched data from rank %d\n",mpi_id_in_gcomm,candidate);
      LOG_PRINT(LOG_DEBUG,"Rank %d also branched GR observations from rank %d\n",mpi_id_in_gcomm,candidate);
      
      md->setCrdsVels(at.get());
      
      md->getSimData(nullptr,&fv_e,nullptr,nullptr);
      luaItf->set_lua_variable("epot",fv_e.epot());
      luaItf->set_lua_variable("ekin",fv_e.ekin());
      luaItf->set_lua_variable("etot",fv_e.etot());
      
      i_need_branching=false;
      candidate=-1;
    }

    // again this barrier can't be avoided otherwise we will detach memory too early
    // NOTE but we can still do something else before checking the barrier completion ...
    MPI_Ibarrier(global_comm,&barrier_req);
    
    // the following is independent computations to hide latency ...
    // check running time and depending on return value perform a backup before exiting
    if(check_running_time())
      doBackupAndExit();
    
    // now we check barrier completion here
    MPI_Wait(&barrier_req,MPI_STATUS_IGNORE);
    
    // detach the windows now that the ibarrier has been checked 
    for(const GR_function_name& n: gr_names)
    {
      MPI_Win& w = gr_windows[n];
      MPI_Aint& addr = gr_addr[n];
      MPI_Win_detach(w, gr_observations[n]);
      addr = 0;
    }
    
    // then each rank sends its data to the REF_WALKER
    for(size_t i=0; i<gr_names.size(); i++)
    {
      const string& n = gr_names[i];
      double* recvObs = (my_FV_role == REF_WALKER) ? recvObsMap[n] : nullptr;
      
      MPI_Igather(gr_observations[n],numObs,MPI_DOUBLE,
                  recvObs,numObs,MPI_DOUBLE,
                  masterRank,global_comm,&igather_obs_reqs[i]);
    }

    /*
     * update GR statistics : check if each observable has converged
     * if convergence notify the workers by setting the corresponding boolean flag to true
     * the workers do not wait for this, they already looped to the beginning of the while loop and will start doing md again
     * after doing gr_check steps they will wait for the ibarrier at the end of this block to be completed;
     * this way we hide the latency time required for notifying the walkers behind the md dynamics time
     */
    if(my_FV_role == REF_WALKER)
    {
      MPI_Waitall(gr_names.size(),igather_obs_reqs.data(),MPI_STATUSES_IGNORE);
      
      for(const GR_function_name& n : gr_names)
      {
        double* recvObs = recvObsMap[n]; //.get();
        for(int32_t i=0; i<mpi_gcomm_size; i++)
        {
          const size_t from = i*numObs;
          gr->addNObservations(n,recvObs+from,numObs,i);
        }
      }
      
      gr->updateStatistics();

      gr->describe();
      if(LOG_SEVERITY == LOG_DEBUG)
      {
        for(uint32_t n=0; n<(uint32_t)mpi_gcomm_size; n++)
          gr->describeChain(n);
      }

      vector<bool> convergedVec = gr->check_convergence_all();
      
      bool lconverged = convergedVec[0];
      for(size_t n=1; n<gr_names.size(); n++)
        lconverged &= convergedVec[n];
      
      if(lconverged)
      {
        converged = true;
        MPI_Win_lock_all(MPI_MODE_NOCHECK,converged_window);
        for(int32_t i=1; i<mpi_gcomm_size; i++)
        {
          MPI_Put(&converged, //const void *origin_addr
                  1,MPI_CXX_BOOL, //int origin_count, MPI_Datatype origin_datatype,
                  i,0,  // int target_rank, MPI_Aint target_disp,
                  1,MPI_CXX_BOOL, // int target_count, MPI_Datatype target_datatype,
                  converged_window); // MPI_Win win
        }
        MPI_Win_unlock_all(converged_window);
      }
      
    }

    /*
     * again barrier can't be avoided but we hide it using a non blocking variant
     * the completion check is performed at the beginning of the while loop, AFTER a call to md->doNsteps(gr_check) in order to hide letency behind a small
     * amount of computations, ie md integration of gr_check steps
     */
    MPI_Ibarrier(global_comm,&barrier_req);
    
  }// convergence while loop
  
  LOG_PRINT(LOG_INFO,"Fleming-Viot converged after %.2lf ps\n",fv_local_time);
  fprintf(    stdout,"Fleming-Viot converged after %.2lf ps\n",fv_local_time);
  
  if(my_FV_role == REF_WALKER)
  {
    LOG_PRINT(LOG_INFO,"Fleming-Viot converged after %.2lf ps : Gelman-Rubin statistics are : \n",fv_local_time);
    
    gr->describe();
    if(LOG_SEVERITY == LOG_DEBUG)
    {
      for(uint32_t n=0; n<(uint32_t)mpi_gcomm_size; n++)
        gr->describeChain(n);
    }
    
    for(const GR_function_name& n : gr_names)
      delete[] recvObsMap[n];

    recvObsMap.clear();
  }

  // update reference clock time
  ref_clock_time += fv_local_time;

  //----------------------
  
  // reset GR object before continuing
  if(i_am_master)
    gr->reset_all_chains();
  
  for(const GR_function_name& n: gr_names)
    delete[] gr_observations[n];
  
  gr_observations.clear();

  for(size_t i=0; i<gr_names.size(); i++)
  {
    if(igather_obs_reqs[i] != MPI_REQUEST_NULL)
      MPI_Request_free(&igather_obs_reqs[i]);
  }
  
  MPI_Win_free(&at_window);
  MPI_Win_free(&pbc_window);
  MPI_Win_free(&branching_window);
  MPI_Win_free(&reset_window);
  MPI_Win_free(&converged_window);
  
  for(const GR_function_name& n: gr_names)
  {
    MPI_Win_free(&(gr_windows[n]));
    MPI_Win_free(&(gr_windows_addr[n]));
  }

  // this one costs nothing 
  MPI_Barrier(global_comm);
  
} // end of do_FlemingViot_procedure

