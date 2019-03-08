/**
 * \file parRep.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#ifndef PARREP_HPP_INCLUDED
#define PARREP_HPP_INCLUDED

#include <limits>

#include "parRep_abstract.hpp"

/**
 * \brief This class performs an "original" parrep simulation, i.e. similar in nature to the original concept from A. Voter
 */
class ParRep final : public ParRepAbstract
{

public:

  /**
   * @brief The constructor for the original parrep simulation
   * 
  * @param _dat Simulation data, see global.hpp
  * @param _at Array of coordinates and velocities
  * @param _md The MD interface
  * @param _luaItf The Lua interface
   */
  
  ParRep(DATA& _dat,
         std::unique_ptr<ATOM[]>& _at,
         std::unique_ptr<MD_interface>& _md,
         std::unique_ptr<luaInterface>& _luaItf
        );

  /**
   * @brief the main function running parralel replica simulations
   */
  virtual void run() override;
  
private:

  ///////////////////////////////////////////////////
  /*
   * some simulation parameters only used by this version of parRep
   *  see parRep_abstract.hpp for common parameters
   */
  ///////////////////////////////////////////////////
  
  // parameters for decorrelation stage
  double   t_corr; ///< arbitrary (user-defined) decorrelation time
  uint32_t corr_check; ///< frequency at which to check, during decorrelation, if the ref. walker left the state 
  double corr_local_time = 0.; ///< Physical decorrelation time
  ENERGIES corr_e; ///< a structure for storing decorrelated energies
  
  // parameters for dephasing stage
  double   t_dephase; ///< arbitrary (user-defined) dephasing time
  uint32_t dephase_check; ///< frequency at which to check, during dephasing, if a given replica left the state 
  double dephase_local_time = 0.; ///< Physical dephasing time
  ENERGIES dephase_e; ///< a structure for storing dephased energies
  
  /*
   * DECLARATIONS
   */
  
  /**
   * @brief This performs the decorrelation stage ; executed by masterRank only
   */
  void do_decorrelation();
  
  /**
   * @brief This performs the dephasing stage ; executed by all ranks in parallel
   */
  void do_dephasing();

};

#endif // PARREP_HPP_INCLUDED
