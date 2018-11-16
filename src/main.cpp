/**
 * \file main.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#include <string>

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include "global.hpp"
#include "logger.hpp"
#include "runSim.hpp"

#include "mpi_utils.hpp"

using namespace std;

// -----------------------------------------------------------------------------------------

/*
 * Some global variables initialisation here
 */

/*
 * Errors, warning, etc ... --> logging.
 * Default Level is LOG_WARNING, which means that everything which is at least
 * a warning is printed (it includes Errors also).
 * See logger.h for other possibilities.
 */
LOG_LEVELS LOG_SEVERITY = LOG_WARNING;

/*
 *  End of global variables initialisation
 */

// -----------------------------------------------------------------------------------------

//prototypes of functions in this c++ unit
void help(char *argv[]);

// -----------------------------------------------------------------------------------------
/**
 * \brief   Main entry of the program
 * \details The 2 variables \b argc and \b argv are used for extracting
 *          command line parameters.
 * \param   argc Number of arguments, at least one as \b argv contains at least the program's name.
 * \param   argv Array of character strings containing the arguments.
 * \return  On exit returns EXIT_SUCCESS, EXIT_FAILURE otherwise.
 */
int main(int argc, char* argv[])
{
  int32_t version,subversion;
  MPI_Get_version(&version,&subversion);
  
  char  mpi_lib_version_string[MPI_MAX_LIBRARY_VERSION_STRING];
  int32_t  str_len = MPI_MAX_LIBRARY_VERSION_STRING;
  MPI_Get_library_version(mpi_lib_version_string,&str_len);

  if(version < 3)
  {
    fprintf(stderr,"[ERROR] Your MPI runtime version was found to be %d.%d, however version 3.0 or newer is required !\n",version,subversion);
    
    fprintf(stderr,"[NOTE] If you have compiled this program yourself, you have a valid 3.0 (or newer) runtime available somewhere because this was already checked at compile time: therefore you are probably facing a dynamic library mismatch !! Please check the following environment variables:\n");
    
    fprintf(stderr,"\t PATH : %s\n",getenv("PATH"));
    fprintf(stderr,"\t LIBRARY_PATH : %s\n",getenv("LIBRARY_PATH"));
    fprintf(stderr,"\t LD_LIBRARY_PATH : %s\n",getenv("LD_LIBRARY_PATH"));
    
    return EXIT_FAILURE;
  }
  
  /*
   * This programm allows multiple threads per MPI rank for the MD engine, however only the thread who called
   * MPI_Init_thread here will perform MPI calls ; this implies that the MD engine SHALL NOT do any MPI call !!
   */
  int threads_request_level = MPI_THREAD_FUNNELED;
  int library_provided_threads_support;
  MPI_Init_thread(&argc, &argv, threads_request_level, &library_provided_threads_support);
  if(library_provided_threads_support != threads_request_level)
  {
    // If the MPI runtime does not support MPI_THREAD_FUNNELED (it should always be available but who knows...) : error
    fprintf(stderr,"[ERROR] Your MPI runtime apparently doesn't support MPI_Init_thread called with the following theads support : MPI_THREAD_FUNNELED! \n");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_CUSTOM_ABORT_MACRO();
  }
  
  // be sure that when exiting log files are properly flushed
  atexit(LOG_FLUSH_ALL);

  MPI_Barrier(MPI_COMM_WORLD);

  /* find out MY process ID, and how many processes were started. */
  int32_t my_id     = -1;
  int32_t num_procs = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  
  MPI_Barrier(MPI_COMM_WORLD);

  /* arguments parsing, we need at least "prog_name -i an_input_file"
   * prints some more instructions if needed
   */
  if(argc < 3)
  {
    fprintf(stdout,"[Info] No input file ! \n");
    help(argv);
    
    MPI_CUSTOM_ABORT_MACRO();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  char inpf[FILENAME_MAX] = "";

  // arguments parsing
  for (uint32_t i=1; i<(uint32_t)argc; i++)
  {
    // get name of input file
    if(!strcasecmp(argv[i],"-i"))
    {
      sprintf(inpf,"%s",argv[++i]);
    }
    /*
     * reopen stdout to unique file based on the user specified template name ;
     * expects one dot and an extension at the end (exception: redirection to /dev/null is allowed)
     */
    else if(!strcasecmp(argv[i],"-o"))
    {
      string my_name(argv[++i]);
      MPIutils::mpi_get_unique_name(my_name);
      stdout = freopen(my_name.c_str(),"w",stdout);
    }
    /*
     * reopen stderr to unique file based on the user specified template name ;
     * expects one dot and an extension at the end (exception: redirection to /dev/null is allowed)
     */
    else if(!strcasecmp(argv[i],"-e"))
    {
      string my_name(argv[++i]);
      MPIutils::mpi_get_unique_name(my_name);
      stderr = freopen(my_name.c_str(),"w",stderr);
    }
    // specify the logging level
    else if(!strcasecmp(argv[i],"-log"))
    {
      if(!strcasecmp(argv[++i],"no"))
      {
        LOG_SEVERITY = LOG_NOTHING;
      }
      else if(!strcasecmp(argv[i],"warn"))
      {
        LOG_SEVERITY = LOG_WARNING;
      }
      else if(!strcasecmp(argv[i],"info"))
      {
        LOG_SEVERITY = LOG_INFO;
      }
      else if(!strcasecmp(argv[i],"dbg"))
      {
        LOG_SEVERITY = LOG_DEBUG;
      }
      else
        fprintf(stdout,"[Warning] Unknown log level '%s' : default value used.\n\n",argv[i]);
    }
    // print help and proper exit
    else if( !strcasecmp(argv[i],"-h") || !strcasecmp(argv[i],"-help") || !strcasecmp(argv[i],"--help") )
    {
      help(argv);
      MPI_CUSTOM_ABORT_MACRO();
    }
    // error if unknown command line option
    else
    {
      fprintf(stdout,"[Error] Argument '%s' is unknown.\n",argv[i]);
      help(argv);
      MPI_CUSTOM_ABORT_MACRO();
    }
  }

  // wait for parsing of command line on each rank
  MPI_Barrier(MPI_COMM_WORLD);
  
  //prepare log files if necessary
  try
  {
    init_logfiles();
  }
  catch(exception& e)
  {
    fprintf(stderr,"std::exception captured when initializing log files ! Error message is %s \n\n",e.what());
    MPI_CUSTOM_ABORT_MACRO();
  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(my_id==0)
  {
    // Print date and some env. variables
    fprintf(stdout,"Welcome to %s ! Command line arguments succesfully parsed, now initialising...\n\n",argv[0]);
    fprintf(stdout,"Logging level is : %s : see the documentation to see which .log files are generated, and what they contain.\n\n",get_loglevel_string());
    fprintf(stdout,"DATE : %s\n",get_time());
    fprintf(stdout,"USER : %s\n",getenv("USER"));
    fprintf(stdout,"PWD : %s\n",getenv("PWD"));
  }

  // run the simulation
  const string input(inpf);
  run_simulation(input);

  // closing log files is the last thing to do as errors may occur at the end
  close_logfiles();
  
  fprintf(stdout,"End of program\n");
  
  fflush(stdout);
  fflush(stderr);

  MPI_Finalize();
  
  return EXIT_SUCCESS;
}

// -----------------------------------------------------------------------------------------
/**
 * \brief   This function simply prints a basic help message.
 *
 * \details If any of \b -h or \b -help or \b --help are provided on the command line this help message is printed.\n
 *          If no command line parameter is present this message is also printed.\n
 *          If an unknown command line parameter is present this message is also printed.
 *
 * \param   argv Simply the same array of command line parameters, from function \b #main.
 */
void help(char *argv[])
{
  fprintf(stdout,"Need at least one argument : %s -i an_input_file\n",argv[0]);
  fprintf(stdout,"optional args : -o [stdout_to_file] -e [stderr_to_file] -log [logging level, one of { no | warn | info | dbg }] \n");
  fprintf(stdout,"Example : \n %s -i input_file -o out.txt -log info \n\n",argv[0]);
  fprintf(stdout,"The default logging level is 'warn' \n");
}

/*!
 * \mainpage A C++ implementation of the Generalized Parallel Replica algorithm
 *
 * Metastability is one of the major encountered obstacle when performing long molecular dynamics simulations,
 * and many methods were developed with the aim of addressing this challenge.
 * 
 * The ``Parallel Replica''(parRep) dynamics is known for allowing to simulate very long
 * trajectories of metastable Langevin dynamics in the materials science community,
 * but it relies on assumptions that can hardly be transposed to the world of biochemical simulations.
 *
 * The later developed ``Generalized parRep'' variant solves those issues, but it was not applied to
 * significant systems of interest so far.
 *
 * This software is a new implementation of the Generalized parallel replica method, targeting
 * frequently encountered metastable biochemical systems.
 * 
 * \image html  main.png
 * \image latex main.pdf
 * 
 */
