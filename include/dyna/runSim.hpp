/**
 * \file runSim.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#ifndef RUNSIM_HPP_INCLUDED
#define RUNSIM_HPP_INCLUDED

#include <string>

#include "global.hpp"

/**
 * @brief This function, called from main, is in charge of:
 * 
 *  + Parsing Input file
 *  + Setting up the desired type of parRep run_simulation
 *  + Running the simulation
 */
void run_simulation(const std::string& inpf);

#endif // RUNSIM_HPP_INCLUDED

