/**
 * \file runSim.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
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
 * 
 * @param inpf The path to the input Lua file, mandatory.
 * @param seeds_inp_file The path to the unique binary file containing random seeds to be loaded, optional.
 * @param seeds_out_file The path to the unique binary file containing random seeds to be loaded, optional.
 */
void run_simulation(const std::string& inpf,
                    const std::string& seeds_inp_file = NULLFILE,
                    const std::string& seeds_out_file = NULLFILE
                   );

#endif // RUNSIM_HPP_INCLUDED

