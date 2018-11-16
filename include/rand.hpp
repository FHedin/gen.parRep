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

#include <string>

/**
 * @brief Call this function for initialising a Mersenne twister 19937 random numbers generator, using random seeds.
 * 
 * The generator is seeded by using integers taken from a std::random_device (usually by reading /dev/urandom where available, or using the hardware instruction RDRND on compatible CPUs)
 */
void init_rand();

/**
 * @brief Call this function for initialising a Mersenne twister 19937 random numbers generator, using as seeds a list of unsigned 64 bits integers (encoded as a comma separated string, at least 5 seeds required).
 * 
 * @attention it is crucial to provide a different string on each MPI rank, otherwise the generators will be correlated/identical !!!!
 * 
 * @param seeds seeds is a std::string containing n seeds, comma separated (e.g. seed = "123,456,789"). n is >= 5
 */
void init_rand(const std::string& seeds);

/**
 * @brief Call this function for obtaining a uniformly distributed random number in the range [0.0,1.0]
 *
 * @return A random number uniformly distributed in the range [0.0,1.0]
 */
double get_double_unif_0_1();

/**
 * @brief Call this function for obtaining a uniformly distributed random number in the range [-1.0,1.0]
 *
 * @return A random number uniformly distributed in the range [-1.0,1.0]
 */
double get_double_unif_m1_p1();

/**
 * @brief Call this function for obtaining a normally distributed random number (mean = 0.0 | sd = 1.0)
 *
 * @return A random number normally distributed (mean = 0.0 | sd = 1.0)
 */
double get_double_normal_0_1();

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
