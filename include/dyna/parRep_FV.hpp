/**
 * \file parRep_FV.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#ifndef PARREP_FV_HPP_INCLUDED
#define PARREP_FV_HPP_INCLUDED

#include "parRep_abstract.hpp"

#include "GelmanRubin.hpp"

/**
 * @brief This class represents a parRep simulation using a Fleming-Viot particle process in order to converge to the QSD
 * 
 * Convergence is achieved by using Gelman-Rubin statistics
 * 
 * This class inherits from the abstract class ParRepAbstract ; it is necessary to implement overriding methods declared
 * as purely virtual
 * 
 * This class can also be inherited from, for defining variants of ParRep still based on the
 * Fleming-Viot procedure: in that class the following methods marked as virtual can be overriden
 */
class ParRepFV : public ParRepAbstract
{

public:

  /**
   * @brief The constructor for the F-V parRep variant
   * 
   * @param _dat Simulation data, see global.hpp
   * @param _at Array of coordinates and velocities
   * @param _md The MD interface
   * @param _luaItf The Lua interface
   * @param ignore_parrepabstract_setup Whether the abstract class should ignore some of its initilaisation procedures, in cas the derived class takescare of 
   *                                    doing it
   * 
   */
  ParRepFV(DATA& _dat,
           std::unique_ptr<ATOM[]>& _at,
           std::unique_ptr<MD_interface>& _md,
           std::unique_ptr<luaInterface>& _luaItf,
           bool ignore_parrepabstract_setup = false
  );
  
  /**
   * @brief the main function running a Fleming-Viot ParRep simulation
   */
  virtual void run() override;
  
protected:
  
  /**
   * \enum 
   *  During the Fleming-Viot procedure, the MPI ranks are divided in 2 groups:
   *  + One reference walker: it will be the MPI masterRank, if it leaves the state the procedure is
   *    stopped (equivalent to the decorrelation phase of the original ParRep Algorithm) 
   *  + Fleming-Viot walkers: their role is to accumulate Gelman-Rubin statistics for assessing convergence
   *    (equivalent to dephasing), if they leave the state they are cloned from one of other other FV_WALKER
   */
  enum FV_ROLE
  {
    REF_WALKER, ///< the masterRank is a REF_WALKER
    FV_WALKER   ///< the N-1 other ranks have a role of FV_WALKER
  };
  
  FV_ROLE my_FV_role; ///< keeps trace of the role of this F-V replica
  
  /// set string names corresponding to the FV_ROLE 's
  const std::map<FV_ROLE,std::string> role_string = { {FV_ROLE::REF_WALKER,"REF_WALKER"} , {FV_ROLE::FV_WALKER,"FV_WALKER"} };
  
  /*
   * some simulation parameters only used by this version of ParRepFV
   *  see parRep_abstract.hpp for common parameters
   */
  
  // parameters for Fleming-Viot combined stage
  uint32_t gr_check; ///< frequency at which to accumulate G-R observables 
  uint32_t fv_check; ///< frequency at whick to check if the F-V reference walker left the state
  uint32_t discard_first_N_grObs;  ///< FV convergence will be tested when at least this amount of observations have already been accumulated
  
  double fv_local_time = 0.; ///< physical time during Fleming Viot procedure
  ENERGIES fv_e;   ///< ENERGIES of the system during Fleming Viot procedure
  
  double grTol;  ///< Tolerance for the Gelman-Rubin statistics
  
  /** a map for storing GR observables parsed from lua input file :
   *  -> the key is a string containing the name of the observable ;
   *  -> the value is a an object typedefed as a 'GR_function', see lua interface header for info ;
   * such functions take no parameters and returns a double precision number, i.e. an 'observation'
   */
  const std::map<GR_function_name,GR_function>& gr_functions;
  
  /**
   * a vector of the names used for indexing the above defined map
   */
  const std::vector<GR_function_name>& gr_names;
  
//   /**
//    *  the C++ interface for collecting Gelman-Rubin statistics
//    *  For the moment only rank masterRank will manage it.
//    *  The others regularly send their data (stored in gr_observations) to rank masterRank
//    */
//   std::unique_ptr<GelmanRubinAnalysis> gr = nullptr;
  
  /*
   * DECLARATIONS : parrep simulation split in several methods
   */

  /**
   * @brief This performs simultaneously a decorrelation and a dephasing stage
   *  It is using a Fleming viot particle process:
   *  no a priori t_corr or t_dephase are required, as Gelman-Rubin statistics
   *  are used for achieving convergence
   */
  void do_FlemingViot_procedure();

};

#endif // PARREP_FV_HPP_INCLUDED
