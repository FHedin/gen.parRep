/**
 * \file mpi_utils.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#include <cstdint>

#include "mpi_utils.hpp"

#include "logger.hpp"

using namespace std;

void MPIutils::mpi_broadcast_atom_array(DATA& dat, ATOM at[], const int32_t mpi_root, MPI_Comm& comm)
{
  MPI_Bcast(&(at[0]),dat.natom*sizeof(ATOM),MPI_BYTE,
            mpi_root,comm);
  
  MPI_Bcast(dat.pbc.data(),sizeof(dat.pbc),MPI_BYTE,
            mpi_root,comm);
}

// void MPIutils::mpi_ibroadcast_atom_array(DATA& dat, ATOM at[], const int32_t mpi_root, MPI_Comm& comm)
// {
//   array<MPI_Request,2> reqs = {{MPI_REQUEST_NULL,MPI_REQUEST_NULL}};
//   
//   MPI_Ibcast(&(at[0]),dat.natom*sizeof(ATOM),MPI_BYTE,
//              mpi_root,comm,&(reqs[0]));
//   
//   MPI_Ibcast(dat.pbc.data(),sizeof(dat.pbc),MPI_BYTE,
//              mpi_root,comm,&(reqs[1]));
//   
//   MPI_Waitall(2,reqs.data(),MPI_STATUSES_IGNORE);
//   
//   MPI_Request_free(&reqs[0]);
//   MPI_Request_free(&reqs[1]);
// }

// MPI_Request* MPIutils::mpi_ibroadcast_atom_array_nowait(DATA& dat, ATOM at[], const int32_t mpi_root, MPI_Comm& comm)
// {
//   MPI_Ibcast(&(at[0]),dat.natom*sizeof(ATOM),MPI_BYTE,
//              mpi_root,comm,&(reqs[0]));
//   
//   MPI_Ibcast(dat.pbc.data(),sizeof(dat.pbc),MPI_BYTE,
//              mpi_root,comm,&(reqs[1]));
//   
//   return reqs.data();
// 
// }

// void MPIutils::mpi_ibroadcast_atom_array_nowait(DATA& dat, ATOM at[], const int32_t mpi_root,
//                                                 MPI_Comm& comm, MPI_Request* reqA, MPI_Request* reqB)
// {
//   MPI_Ibcast(&(at[0]),dat.natom*sizeof(ATOM),MPI_BYTE,
//              mpi_root,comm,reqA);
//   
//   MPI_Ibcast(dat.pbc.data(),sizeof(dat.pbc),MPI_BYTE,
//              mpi_root,comm,reqB);
//   
// }

void MPIutils::mpi_get_unique_name(string& orig_name)
{
  if(orig_name == NULLFILE)
  {
    return;
  }
  
  /* find out MY process ID */
  int32_t my_id     = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  
  // store the my_id in a string and add a . after
  string id = to_string(my_id);
  id += ".";
  
  // insert string id just after the last . of the original file name (i.e. before extension)
  auto pos = orig_name.find_last_of('.');
  orig_name.insert(pos+1,id);
}

void MPIutils::mpi_abort_commworld_with_message(const char* fileName, const char* functionName, int32_t lineNumber)
{
  const string fileString = string(fileName);
  const string funcString = string(functionName);
  const string lineString = to_string(lineNumber);
  
  const string message = "MPI ABORT required from file " + fileString + " in function " +
                   funcString + " at line " + lineString +
                   " !!! Calling MPI_Abort with error code -1 now ...";
  
  fprintf(stdout,"%s\n",message.c_str());
  fprintf(stderr,"%s\n",message.c_str());

  close_logfiles();
  
  MPI_Abort(MPI_COMM_WORLD,-1);
}
