/**
 * \file logger.hpp
 * \copyright The 3-clause BSD license is applied to this software.
 * \author Florent Hédin
 * \author Tony Lelièvre
 * \author École des Ponts - ParisTech
 * \date 2016-2019
 */

#ifndef LOGGER_H_INCLUDED
#define LOGGER_H_INCLUDED

#include <cstdint>

/**
*   \enum       LOG_LEVELS
*
*   \brief      The available logging Levels.
*
*   \details This will control what is written
*   to logging files ; those files are named :
*   \li warning.log
*   \li info.log
*   \li debug.log
*
*   The higher the level, the more text is written to those text files, so be careful with long simulations
*   coupled to high levels of logging.
*
*
*   The different levels, members of the \b enum \b #LOG_LEVELS are :
*   \li \b #LOG_NOTHING : At this level there is no logging at all.
*   \li \b #LOG_WARNING : A Warning is a non-critical event, i.e. the event is just reported but the program
*       execution continues, for example a missing or wrong parameter for which a default value is going to be used.
*   \li \b #LOG_INFO : An Info message reports to the user a possibly useful information, an event with no consequence for the simulation.
*       For example, that the trajectory was succesfully written as planed.
*   \li \b #LOG_DEBUG : A debug message is really technical and verbose, for example dumping a variable.
*       This should be enabled for bug tracking.
*
*   The logging is progressive, enabling \b #LOG_INFO means that the higher level \b #LOG_WARNING is also enabled.
*
*   The default is set to \b #LOG_WARNING. \n
*
*/
enum LOG_LEVELS
{
    LOG_NOTHING, /*!< No log file is created : this is not recommended as no information is reported */
    LOG_WARNING, /*!< Warnings are reported to \b warning.log . This is the default. */
    LOG_INFO,    /*!< Info messages are reported to \b info.log , and also warnings to their respective file. */
    LOG_DEBUG    /*!< A large amount of debugging messages are written to \b debug.log . The previous levels are still written to their respective files. */
};

extern LOG_LEVELS LOG_SEVERITY;

/**
 * \brief   Prepares LOGGING I/O if necessary, depending of the value of \b #LOG_SEVERITY
 */
void init_logfiles();

/**
 * \brief   Before exiting the program, closes properly the logging files.
 */
void close_logfiles();

/**
 * \brief   Returns a string corresponding to the logging level.
 * \return An array of char containing the description of the logging level, for example "LOG_WARNING"
 */
const char* get_loglevel_string();

/**
 * \brief Puts in a string the current date, and returns it
 *
 * \details This function returns the date with the following format :
 * \li DAY/MONTH/YEAR-HH:MM:SS
 * \li For example : 17/10/2013-16:26:36
 *
 * \return A character string filled with the date and hour
 */
char* get_time();

/**
 * \brief \b #LOG_PRINT adds a given line of text (printf-like formatted) to a given log file
 *
 * \param mesg_severity The logging severity
 * \param fmt A printf like format string
 * \param ... parameters forwarded to vfprintf
 */
void LOG_PRINT(LOG_LEVELS mesg_severity, const char *fmt, ...);

/**
 * @brief Flush all opened files ; slow but useful when debugging if there may be a crash soon...
 */
void LOG_FLUSH_ALL();

#endif // LOGGER_H_INCLUDED
