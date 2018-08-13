/**
 * \file parRep.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "logger.hpp"

#include "mpi_utils.hpp"

#include "parRep.hpp"

using namespace std;

ParRep::ParRep(DATA& _dat,
              unique_ptr<ATOM[]>& _at,
              unique_ptr<MD_interface>& _md,
              unique_ptr<luaInterface>& _luaItf
              ) : ParRepAbstract(_dat,_at,_md,_luaItf)
{
  // get simulation parameters parsed from lua script
  equil_steps     = (uint32_t) stoul(params.at("equilibrationSteps"));
  t_corr          = stod(params.at("tauDecorr"));
  corr_check      = (uint32_t) stoul(params.at("checkDecorr"));
  t_dephase       = stod(params.at("tauDephase"));
  dephase_check   = (uint32_t) stoul(params.at("checkDephase"));
  dynamics_check  = (uint32_t) stoul(params.at("checkDynamics"));
  
  // and also get database parameters
  db_backup_frequency_ps = stod(params.at("dbFreq"));
  
  // calculate polling time from input params
  t_poll = ((double)dynamics_check)*dat.timestep ;
}

void ParRep::run()
{
  fprintf(stdout,"\nRunning a standard ParRep.\n\n");

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
    LOG_PRINT(LOG_INFO,"New ParRep loop iteration : ref_clock_time is %.2lf ps \n",ref_clock_time);
    fprintf(stdout,"\n//---------------------------------------------------------------------------------------------//\n");
    fprintf(stdout    ,"New ParRep loop iteration : ref_clock_time is %.2lf ps \n\n",ref_clock_time);
    
    /*
     * Phase 1 : Perform decorrelation on masterRank only
     */
    corr_local_time = 0.;
    corr_e = ENERGIES();
    // only masterRank decorrelates
    if(i_am_master)
    {
      do_decorrelation();
    }// end
    MPI_Barrier(global_comm);

    MPI_Bcast(&corr_local_time,1,MPI_DOUBLE,masterRank,global_comm);
    MPIutils::mpi_broadcast_atom_array(dat,at.get(),masterRank,global_comm);
    MPI_Bcast(&ref_clock_time,1,MPI_DOUBLE,masterRank,global_comm);
    
    // not required on masterRank : provide to openmm coordinates and velocities from decorrelation
    if(!i_am_master)
    {
      md->setCrdsVels(at.get());
    }

    // update lua interface
    luaItf->set_lua_variable("epot",corr_e.epot());
    luaItf->set_lua_variable("ekin",corr_e.ekin());
    luaItf->set_lua_variable("etot",corr_e.etot());
    luaItf->set_lua_variable("referenceTime",ref_clock_time);

    /*
     * Phase 2 : Perform dephasing on each rank
     */
    
    dephase_local_time = 0.;
    dephase_e = ENERGIES();
    // all ranks dephase
    do_dephasing();
    
    // update lua interface
    luaItf->set_lua_variable("epot",dephase_e.epot());
    luaItf->set_lua_variable("ekin",dephase_e.ekin());
    luaItf->set_lua_variable("etot",dephase_e.etot());

    /*
     * phase 3 : Run parallel dynamics 
     */
    dyna_local_time = 0.;
    dyna_e = ENERGIES();
    breakerID = numeric_limits<uint32_t>::max();
    dyna_cycles_done = 0;
      
    // each rank does its independent dynamics until exit
    do_dynamics();
    
    /*
     * send id of rank that escaped the state to all others;
     * in case several ranks left the state at the same time, we choose the one with the minimum breakerID and store it to breakID
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
    
    //then broadcast from the breakID rank to all others
    MPI_Bcast(&ref_clock_time,1,MPI_DOUBLE,masterRank,global_comm);
    MPI_Bcast(&dyna_local_time,1,MPI_DOUBLE,breakerID,global_comm);
    MPI_Bcast(&dyna_e.ene[0],3,MPI_DOUBLE,breakerID,global_comm);
    MPI_Bcast(&dyna_cycles_done,1,MPI_UINT32_T,breakerID,global_comm);
    
    MPIutils::mpi_broadcast_atom_array(dat,at.get(),breakerID,global_comm);
    
    /*
     * Calculating the escape time
     * This formula comes from the article "parRep for markov chains" (Aristoff, Lelièvre, Simpson)
     */
    double escape_time =  (double)(mpi_gcomm_size-1);
    escape_time *= ((double)dyna_cycles_done)*t_poll;
    escape_time += ((double)breakerID)*t_poll ;
    escape_time += dyna_local_time;
    
    // ref_clock_time already contains the tau_corr so just add the escape time
    ref_clock_time += escape_time;
    
    LOG_PRINT(LOG_INFO,"Parallel dynamics stopped because rank %u escaped (after %.2lf ps).\n",  breakerID,dyna_local_time);
    fprintf(stdout    ,"Parallel dynamics stopped because rank %u escaped (after %.2lf ps).\n\n",breakerID,dyna_local_time);
    
    LOG_PRINT(LOG_INFO,"ParRep escape time is %.2lf ps\n",  escape_time);
    fprintf(stdout,    "ParRep escape time is %.2lf ps\n\n",escape_time);
    
    // update lua interface
    luaItf->set_lua_variable("epot",dyna_e.epot());
    luaItf->set_lua_variable("ekin",dyna_e.ekin());
    luaItf->set_lua_variable("etot",dyna_e.etot());
    luaItf->set_lua_variable("referenceTime",ref_clock_time);
    
    // save the state and the escape time to the database
    if(i_am_master)
    {
      db_insert(true,t_corr,escape_time);
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
    
    // now all ranks need to check again what is the current state in order to be ready for next iteration
    lua_state_init();
    
    /*
     * Ready to loop again
     */
    
  }while( ref_clock_time < ((double)dat.nsteps * dat.timestep) );
  
  if(i_am_master)
  {
    // backup database to file if required before exiting simulation 
    db_backup();
    db_close();
  }
  
} // end of function run

// Definition of static functions declared at the top of this file

///////////////////////////////////////////////////
/*
 * Parrep functions : represent successive stages of the algorithm
 */
///////////////////////////////////////////////////

void ParRep::do_decorrelation()
{
    LOG_PRINT(LOG_INFO,"Rank 0 performing decorrelation with tau_corr = %.2lf ps\n",  t_corr);
    fprintf(stdout,    "Rank 0 performing decorrelation with tau_corr = %.2lf ps\n\n",t_corr);
    
    const double rt = (double) corr_check * dat.timestep;
    
    do
    {
      // do some steps
      md->doNsteps(corr_check);
      // time update
      corr_local_time += rt;
      
      md->getState(nullptr,&corr_e,nullptr,at.get());
      luaItf->set_lua_variable("epot",corr_e.epot());
      luaItf->set_lua_variable("ekin",corr_e.ekin());
      luaItf->set_lua_variable("etot",corr_e.etot());

      left_state = lua_check_state();

      if(left_state)
      {
        LOG_PRINT(LOG_DEBUG,"Rank %d left the state during decorrelation\n",mpi_id_in_gcomm);
        
        ref_clock_time += corr_local_time;
        luaItf->set_lua_variable("referenceTime",ref_clock_time);
        
        // save the previously left state to the list of visited ones
        db_insert(false,0.0,corr_local_time);
        
        // re-initialise local time and set the new reference energy
        corr_local_time = 0.;
        
        // if the system is not in a new state, also perform transient propagation
        if(lua_check_transient())
        {
          do_transient_propagation(corr_e);
          
          // update lua interface
          luaItf->set_lua_variable("epot",corr_e.epot());
          luaItf->set_lua_variable("ekin",corr_e.ekin());
          luaItf->set_lua_variable("etot",corr_e.etot());
          luaItf->set_lua_variable("referenceTime",ref_clock_time);
          
          // call state_init again now that we are sure to be within a state
          lua_state_init();
        }
      }
      
      // check running time and depending on return value perform a backup before exiting
      bool need_stop = check_running_time();
      if(need_stop) doBackupAndExit();
      
    } while(corr_local_time < t_corr);
    
    LOG_PRINT(LOG_INFO,"Rank 0 found a valid decorrelated state after %.2lf ps\n",corr_local_time);
    fprintf(stdout,    "Rank 0 found a valid decorrelated state after %.2lf ps\n\n",corr_local_time);

    ref_clock_time += corr_local_time;
    luaItf->set_lua_variable("referenceTime",ref_clock_time);
    
    // get coordinates
    md->getState(nullptr,&corr_e,nullptr,at.get());

}

void ParRep::do_dephasing()
{
  /* a copy of the decorrelated coordinates is made.
   *  If during dephasing stage a system leaves the state, the system is rolled-back
   *  to the initial state
   */ 
  unique_ptr<ATOM[]> rollBack(new ATOM[dat.natom]);
  array<PBC,3> rollPBC;
  memcpy(rollBack.get(),at.get(),dat.natom*sizeof(ATOM));
  memcpy(rollPBC.data(),dat.pbc.data(),sizeof(dat.pbc));
  double rollB_time = 0.0;
  
  LOG_PRINT(LOG_INFO,"Rank %d performing dephasing with tau_dephasing = %.2lf ps\n",  mpi_id_in_gcomm,t_dephase);
  fprintf(stdout,    "Rank %d performing dephasing with tau_dephasing = %.2lf ps\n\n",mpi_id_in_gcomm,t_dephase);
  
  const double rt = (double) dephase_check * dat.timestep;
  
  do
  {
    // do some steps
    md->doNsteps(dephase_check);
    dephase_local_time += rt;

    md->getState(nullptr,&dephase_e,nullptr,at.get());
    luaItf->set_lua_variable("epot",dephase_e.epot());
    luaItf->set_lua_variable("ekin",dephase_e.ekin());
    luaItf->set_lua_variable("etot",dephase_e.etot());
    
    left_state = lua_check_state();
    
    // if this rank left the initial state, reset coordinates and run again
    if(left_state)
    {
      dephase_local_time = rollB_time;
      memcpy(at.get(),rollBack.get(),dat.natom*sizeof(ATOM));
      memcpy(dat.pbc.data(),rollPBC.data(),sizeof(dat.pbc));
      md->setCrdsVels(at.get());
      md->setSimClockTime(dephase_local_time);
      LOG_PRINT(LOG_DEBUG,"Rank %d left the ref. state, resetting coordinates to a previous backup and run again.\n",mpi_id_in_gcomm);
    }
    else
    {
      rollB_time = dephase_local_time;
      memcpy(rollBack.get(),at.get(),dat.natom*sizeof(ATOM));
      memcpy(rollPBC.data(),dat.pbc.data(),sizeof(dat.pbc));
    }
    
    // check running time and depending on return value perform a backup before exiting
    bool need_stop = check_running_time();
    if(need_stop) doBackupAndExit();
    
  } while(dephase_local_time < t_dephase);
  
  LOG_PRINT(LOG_INFO,"Rank %d dephased properly (after %.2lf ps).\n",  mpi_id_in_gcomm,dephase_local_time);
  fprintf(stdout,    "Rank %d dephased properly (after %.2lf ps).\n\n",mpi_id_in_gcomm,dephase_local_time);
  
  md->getState(nullptr,&dephase_e,nullptr,at.get());

}
