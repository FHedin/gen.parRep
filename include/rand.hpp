/**
 * \file rand.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#ifndef RAND_HPP_INCLUDED
#define RAND_HPP_INCLUDED

#include <cstdint>

/**
 * @brief Call this function for initialising a Mersenne twister 19937 random numbers generator
 * (one per MPI rank, independantly initialised)
 */
void init_rand();

/**
 * @brief Call this function for obtaining a uniformly distributed random number in the range [0.0,1.0]
 *
 * @return A random number uniformly distributed in the range [0.0,1.0]
 */
double get_double_0_1();

/**
 * @brief Call this for obtaining a random unsigned 32 bits integer
 * 
 * @return a random unsigned 32 bits integer between 0 and UINT32_MAX
 */
uint32_t get_uint32();

/**
 * @brief Call this for obtaining a random signed 32 bits integer
 * 
 * @return a random signed 32 bits integer between INT32_MIN and INT32_MAX
 */
int32_t  get_int32();

/**
 * @brief Returns an unsigned int within [min;max] (like a dice roll)
 * 
 * @param min lower bound of random range
 * @param max upper bound of random range
 * 
 * @return a random unsigned 32 bits integer between min and max
 */
int32_t get_int32_min_max(const int32_t min, const int32_t max);

#endif // RAND_HPP_INCLUDED
