/**
 * \file rand.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#include <random>
#include <array>
#include <limits>
#include <iostream>
#include <sstream>

#include <cstdio>

#include "logger.hpp"
#include "rand.hpp"
#include "global.hpp"

using namespace std;

// this is the random numbers generator
static mt19937_64 generator;

// this distribution generates a real random number distributed between 0 and 1
static uniform_real_distribution<double>  double_0_1_unif(0.0,1.0);

// this distribution generates a real random number uniformly distributed between -1.0 and 1.0
static uniform_real_distribution<double>  double_m1_p1_unif(-1.0,1.0);

// this distribution generates a real random number normally distributed (mean 0.0 sd 1.0)
static normal_distribution<double>  double_0_1_norm(0.0,1.0);

// the two following generate a random int or unsigned int between a min and a max, user defined
static uniform_int_distribution<int32_t>  int32_t_dist(numeric_limits<int32_t>::min(),
                                                       numeric_limits<int32_t>::max());

static uniform_int_distribution<uint32_t> uint32_t_dist(numeric_limits<uint32_t>::min(),
                                                        numeric_limits<uint32_t>::max());

void init_rand()
{
  random_device rd;
  LOG_PRINT(LOG_INFO,"Entropy of 'std::random_device' on this machine is : %lf\n",rd.entropy());
  
  /*
   * Each MPI ranks reads 512 uint64_t 'seeds' from random_device
   */
  array<uint64_t,512> seeds;
  for(size_t n=0; n<512; n++)
    seeds[n] = (uint64_t) rd();

  // then each ranks initialises its own mt19937 with the seed sequence
  seed_seq seq(seeds.begin(),seeds.end());
  generator.seed(seq);
  
  // finally for each generator we discard some random numbers to be sure to avoid any correlation
  generator.discard(4*generator.state_size);
  
}

void init_rand(const string& seeds_str)
{
  istringstream split_stream(seeds_str);
  string token;
  
  vector<uint64_t> seeds;
  
  /*
   * Exctract the seeds from the comma separated string ;
   * NOTE It is crucial to use a different string on each MPI rank otherwise the generators will provide the same numbers on each rank !!
   */
  while(getline(split_stream, token, ','))
  {
    seeds.push_back((uint64_t)stoul(token));
  }
  
  if(seeds.size()<5)
  {
    throw runtime_error("Error : when providing a string encoded list of seeds for the random generator, you need to provide at least 5 seeds, but it appears that you have provided only " + to_string(seeds.size()) + " seeds ! \n");
  }
  
  seed_seq seq(seeds.begin(),seeds.end());
  generator.seed(seq);
  
  // finally for each generator we discard some random numbers to be sure to avoid any correlation
  generator.discard(4*generator.state_size);
  
}

double get_double_unif_0_1()
{
  return double_0_1_unif(generator);
}

double get_double_unif_m1_p1()
{
  return double_m1_p1_unif(generator);
}

double get_double_normal_0_1()
{
  return double_0_1_norm(generator);
}

uint32_t get_uint32()
{
  return uint32_t_dist(generator);
}

int32_t  get_int32()
{
  return int32_t_dist(generator);
}

int32_t get_int32_min_max(const int32_t min, const int32_t max)
{
  return (min + (int32_t)(get_double_unif_0_1()*((max - min) + 1)) );
}
