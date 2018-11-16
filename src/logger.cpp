/**
 * \file logger.cpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2018
 */

#include <string>
#include <iostream>

#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <cstdarg>

#include "global.hpp"
#include "logger.hpp"
#include "mpi_utils.hpp"

#define DATE_BUFFER_SIZE 128
#define MSG_BUFFER_SIZE  256

// those variabes are persisting but only accessible from this file
static FILE *F_WARN , *F_INFO , *F_DEBUG ;
static time_t rawtime;
static bool filesOpened=false;

using namespace std;

void init_logfiles()
{
  F_WARN  = nullptr;
  F_INFO  = nullptr;
  F_DEBUG = nullptr;

  //--------------------------------------------------
  
  if(LOG_SEVERITY > LOG_NOTHING )
  {
    string my_name="warning.log";
    MPIutils::mpi_get_unique_name(my_name);
    F_WARN = fopen(my_name.c_str(),"wt");
    
    if(F_WARN == nullptr)
    {
      throw runtime_error("Error : impossible to open-write file " + my_name + " !\n");
    }
  }
  
  //--------------------------------------------------
  
  if(LOG_SEVERITY > LOG_WARNING )
  {
    string my_name="info.log";
    MPIutils::mpi_get_unique_name(my_name);
    F_INFO = fopen(my_name.c_str(),"wt");
    
    if(F_INFO == nullptr)
    {
      throw runtime_error("Error : impossible to open-write file " + my_name + " !\n");
    }
  }
  
  //--------------------------------------------------
  

  if(LOG_SEVERITY > LOG_INFO )
  {
    string my_name="debug.log";
    MPIutils::mpi_get_unique_name(my_name);
    F_DEBUG = fopen(my_name.c_str(),"wt");
    
    if(F_DEBUG == nullptr)
    {
      throw runtime_error("Error : impossible to open-write file " + my_name + " !\n");
    }
  }
  
  //--------------------------------------------------
  
  filesOpened=true;
}

void close_logfiles()
{
  LOG_FLUSH_ALL();

  if(filesOpened)
  {
    if(LOG_SEVERITY > LOG_NOTHING )
      fclose(F_WARN);

    if(LOG_SEVERITY > LOG_WARNING )
      fclose(F_INFO);

    if(LOG_SEVERITY > LOG_INFO )
      fclose(F_DEBUG);
    
    filesOpened=false;
  }
}

const char* get_loglevel_string()
{
  const char* str_log=NULL;

  switch(LOG_SEVERITY)
  {
  case LOG_NOTHING:
    str_log = "LOG_NOTHING";
    break;

  case LOG_WARNING:
    str_log = "LOG_WARNING";
    break;

  case LOG_INFO:
    str_log = "LOG_INFO";
    break;

  case LOG_DEBUG:
    str_log = "LOG_DEBUG";
    break;

  default:
    str_log = "UNKNOWN";
    break;
  }

  return str_log;
}

/**
 * \brief Puts in a string the current date, later used for logging.
 *
 * \details This function returns the date with the following format :
 * \li DAY/MONTH/YEAR-HH:MM:SS
 * \li For example : 17/10/2013-16:26:36
 *
 * \param date A character string filled with the date and hour
 */
static void get_time_ptr(char date[])
{
  struct tm * timeinfo;
  time_t newtime;

  time(&newtime);

  if(newtime != rawtime)
  {
    rawtime = newtime;
    timeinfo = localtime(&rawtime);
    strftime(date,DATE_BUFFER_SIZE,"%d/%b/%Y-%H:%M:%S-%s",timeinfo);
  }
}

char* get_time()
{
  static char date[DATE_BUFFER_SIZE];
  struct tm * timeinfo;
  time_t newtime;

  time(&newtime);

  timeinfo = localtime(&newtime);
  strftime(date,DATE_BUFFER_SIZE,"%d/%b/%Y-%H:%M:%S-%s",timeinfo);

  return date;
}

void LOG_PRINT(LOG_LEVELS mesg_severity, const char *fmt, ...)
{
  
  if(!filesOpened)
    return;
  
  // if the severity of the current message is at least equal to the global level we print it
  if (mesg_severity <= LOG_SEVERITY)
  {
    va_list list;

    FILE *FP=NULL;

    static char event_date[DATE_BUFFER_SIZE]="";
    static char message[MSG_BUFFER_SIZE]="";

    // get the date string
    get_time_ptr(event_date);

    // depending of the severity, chose the correct file for writting
    switch(mesg_severity)
    {
      case LOG_WARNING:
        FP = F_WARN;
        sprintf(message,"[Warning @ ");
        break;

      case LOG_INFO:
        FP = F_INFO;
        sprintf(message,"[Info @ ");
        break;

      case LOG_DEBUG:
        FP = F_DEBUG;
        sprintf(message,"[Debug @ ");
        break;

      default:
        break;
    }

    // write date contained in event_date to message
    strncat(message,event_date,DATE_BUFFER_SIZE);
    strncat(message,"]\t",3);

    // print message to file
    fprintf(FP,"%s",message);

    // prepare processing of the variable arguments list
    // fmt is the last non-optional argument
    va_start(list,fmt);

    // now forward what have to be printed to vfprintf
    vfprintf(FP,fmt,list);

    va_end(list);

  } // end of if (mesg_severity >= LOG_SEVERITY)

}

void LOG_FLUSH_ALL()
{
  if(filesOpened)
  {
    if(LOG_SEVERITY > LOG_NOTHING )
      fflush(F_WARN);

    if(LOG_SEVERITY > LOG_WARNING )
      fflush(F_INFO);

    if(LOG_SEVERITY > LOG_INFO )
      fflush(F_DEBUG);
  }
  
  //also flush stdout and stderr
  fflush(stdout);
  fflush(stderr);
}

