/**
 * \file rand.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#include <random>
#include <array>
#include <limits>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <cstdio>

#include "global.hpp"
#include "logger.hpp"
#include "mpi_utils.hpp"
#include "rand.hpp"

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

void init_rand(SEEDS_IO io_type, const string& seeds_file_name)
{
  
  bool save_seeds = false;
  bool load_seeds = false;
  
  switch(io_type)
  {
    case SEEDS_IO::NONE :
      break;
      
    case SEEDS_IO::SAVE_TO_FILE :
      save_seeds = true;
      break;
      
    case SEEDS_IO::LOAD_FROM_FILE :
      load_seeds = true;
      break;
      
    default:
      throw runtime_error("Error : unsupported SEEDS_IO mode in function '" + string(__func__) + "' !!!");
      break;
  }

  /*
   * Each MPI ranks requires 512 uint32_t 'seeds' for initialising its mt19937_64 generator
   */
  array<uint32_t,512> seeds;
  
  if(!load_seeds)
  {
    random_device rd;
    LOG_PRINT(LOG_INFO,"Entropy of 'std::random_device' on this machine is : %lf\n",rd.entropy());
    
    for(size_t n=0; n<512; n++)
      seeds[n] = (uint32_t) rd();
  }
  else
  {
    // First enforce sync with a barrier
    MPI_Barrier(MPI_COMM_WORLD);
    
    int32_t myrank = -1, numranks = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    
    // see comments in the following block "if(save_seeds)..." below
    
    MPI_File myfile;
    
    MPI_File_open(MPI_COMM_WORLD, seeds_file_name.c_str(),
                  MPI_MODE_RDONLY, MPI_INFO_NULL,
                  &myfile);
    
    MPI_Offset disp = myrank*512*sizeof(uint32_t);
    
    MPI_File_set_view(myfile, disp,
                      MPI_UINT32_T, MPI_UINT32_T, 
                      "native", MPI_INFO_NULL);
    
    MPI_File_read(myfile, seeds.data(),
                  512, MPI_UINT32_T,
                  MPI_STATUS_IGNORE);
    
    MPI_File_close(&myfile);
    
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  /*
   * Save seeds to a unique binary file for reproducibility
   */
  if(save_seeds)
  {
    // First enforce sync with a barrier
    MPI_Barrier(MPI_COMM_WORLD);
    
    int32_t myrank = -1, numranks = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numranks);
    
    /* Shared file */
    MPI_File myfile;
    
    /*
     * Open the file : it will be a binary file open in write mode, with overwritting
     * 
     * int MPI_File_open(MPI_Comm comm, const char *filename,
     *                   int amode, MPI_Info info,
     *                   MPI_File *fh)
     */
    MPI_File_open(MPI_COMM_WORLD, seeds_file_name.c_str(),
                  MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL,
                  &myfile);
    
    // displacment == absolute displacement in the file where the data will be written (from beginning of file)
    MPI_Offset disp = myrank*512*sizeof(uint32_t);
    
    /*
     * Set the file view : tels at which value of disp data will be written
     * 
     * MPI_File_set_view(MPI_File fh, MPI_Offset disp,
     *                   MPI_Datatype etype, MPI_Datatype filetype,
     *                   const char *datarep, MPI_Info info)
     */
    MPI_File_set_view(myfile, disp,
                      MPI_UINT32_T, MPI_UINT32_T, 
                      "native", MPI_INFO_NULL);
    /*
     * Write buf to the file : 
     * 
     int MPI_File_write(MPI_File fh, const void *buf,
                        int count, MPI_Datatype datatype,
                        MPI_Status *status)
     */
    MPI_File_write(myfile, seeds.data(),
                   512, MPI_UINT32_T,
                   MPI_STATUS_IGNORE);
    
    /*
     * Close the file
     */
    MPI_File_close(&myfile);
    
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // then each ranks initialises its own mt19937 with the seed sequence
  seed_seq seq(seeds.begin(),seeds.end());
  generator.seed(seq);
  
  // finally for each generator we discard some random numbers to be sure to avoid any correlation
  generator.discard(4*generator.state_size);
  
  MPI_Barrier(MPI_COMM_WORLD);
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
