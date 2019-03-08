/**
 * \file mpi_utils.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#ifndef MPI_UTILS_H_INCLUDED
#define MPI_UTILS_H_INCLUDED

#include <string>

#include <mpi.h>

// Explicitly require support of the MPI 3 (or higher) standard
#if MPI_VERSION < 3

#error "An MPI implementation supporting the MPI-3 (or higher) standard is required !"

#endif

#include "global.hpp"
#include "logger.hpp"

namespace MPIutils
{

  /**
   * @brief broadcasting functions for sharing data
   *        used when initialising simulation, or during simulation
   * 
   * @param dat Simulatiuon data
   * @param at Atomic structure
   * @param mpi_root The root process of the broadcast
   * @param comm The MPI communicator
   */
  void mpi_broadcast_atom_array(DATA& dat, ATOM at[], const int32_t mpi_root, MPI_Comm& comm);
  
  /**
   * @brief broadcasting functions for sharing data (non blocking version)
   *        used when initialising simulation, or during simulation
   * 
   * @param dat Simulatiuon data
   * @param at Atomic structure
   * @param mpi_root The root process of the broadcast
   * @param comm The MPI communicator
   */
//   void mpi_ibroadcast_atom_array(DATA& dat, ATOM at[], const int32_t mpi_root, MPI_Comm& comm);

  /**
   * @brief broadcasting functions for sharing data (non blocking version)
   *        to use during simulation for hidding latency of the broadcast by performing other computations before checking the requests
   * 
   * @param dat Simulatiuon data
   * @param at Atomic structure
   * @param mpi_root The root process of the broadcast
   * @param comm The MPI communicator
   * 
   * @return Pointer to an array of requests of size 2, for checking completion of the broadcasts of ATOM[] and PBC structs
   */
//   MPI_Request* mpi_ibroadcast_atom_array_nowait(DATA& dat, ATOM at[], const int32_t mpi_root, MPI_Comm& comm);
  
  /**
   * @brief broadcasting functions for sharing data (non blocking version)
   *        to use during simulation for hidding latency of the broadcast by performing other computations before checking the requests
   * 
   * @param dat Simulatiuon data
   * @param at Atomic structure
   * @param mpi_root The root process of the broadcast
   * @param comm The MPI communicator
   * @param reqA Pointer to a MPI_Request
   * @param reqB Pointer to a MPI_Request
   */
//   void mpi_ibroadcast_atom_array_nowait(DATA& dat, ATOM at[], const int32_t mpi_root,
//                                                  MPI_Comm& comm, MPI_Request* reqA, MPI_Request* reqB);
  
  /**
   * @brief  Function that alters the file name in a unique way when using MPI
   *         (the MPI my_id is added before extension) so that different processes write to different files
   * 
   * @param orig_name Reference to the file name
   */
  void mpi_get_unique_name(std::string& orig_name);

  /**
   * @brief Use this instead of a direct call to MPI_Abort, as it will before print an error message
   * indicating to the user from which file, function, and line of code the abort was called.
   * 
   * example of use : mpi_abort_commworld_with_message(__FILE__,__func__,__LINE__)
   * or you can pass srings and line number manually
   * 
   * See the following MACRO MPI_CUSTOM_ABORT_MACRO which simplifies the use of mpi_abort_commworld_with_message
   * 
   * @param fileName     Name of the file caling this functon, usually __FILE__
   * @param functionName Name of the function caling this functon, usually __func__
   * @param lineNumber   Line number where this function was called, usually __LINE__
   */
  void mpi_abort_commworld_with_message(const char* fileName, const char* functionName, int32_t lineNumber);

} // namespace MPIutils

/**
 * @brief This macro simplifies the use of mpi_abort_commworld_with_message
 * 
 * It is equivalent to a call like : MPIutils::mpi_abort_commworld_with_message(__FILE__,__func__,__LINE__)
 * 
 * Use : MPI_CUSTOM_ABORT_MACRO();
 */
#define MPI_CUSTOM_ABORT_MACRO() MPIutils::mpi_abort_commworld_with_message(__FILE__,__func__,__LINE__)

#endif // MPI_UTILS_H_INCLUDED
