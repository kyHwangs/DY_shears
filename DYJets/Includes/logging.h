#ifndef LOGGING_H
#define LOGGING_H

#include <iostream>

#include <boost/iostreams/filtering_stream.hpp>

/// \cond
namespace YAML
{
class Node;
}
/// \endcond

namespace util
{

/**
 * \brief Module providing support for some advanced logging features.
 *
 * This namespace defines five objects that can be used like the standard \c cout and \c cerr
 * output streams: \ref debug, \ref info, \ref warn, \ref error and \ref fatal. In addition to
 * printing to \c stderr, they provide some semantics: the output of \c debug is usually not
 * shown, while \c fatal is used for very serious errors. These are called logging
 * \ref level "levels". The amount of output printed to the screen can be controlled using the
 * \ref set_primary_level function.
 *
 * Example usage:
 *
 * ~~~{.cc}
 * logging::info << "Starting..." << std::endl;
 * logging::debug << "This is a debug message" << std::endl;
 * logging::warn  << "There is something suspicious" << std::endl;
 * logging::fatal << "Things got terribly wrong" << std::endl;
 * ~~~
 *
 * Output when the log level is \ref info :
 * ~~~
 * [INFO ] Starting...
 * [WARN ] There is something suspicious
 * [FATAL] Things got terribly wrong
 * ~~~
 *
 * \note
 * It is important to flush the streams at every line: use \c std::endl instead of \c "\n".
 *
 * Before piping to it \c cerr, the modules adds the information of which stream was used at the
 * beginning of every line. Colors can be used to distinguish between the different streams: this
 * is done automatically when the output is a terminal. This behavior can be controlled by means of
 * the \ref set_use_color function.
 *
 * It is also possible to define a \ref set_secondary_stream "secondary stream". A copy of the
 * output will be sent to this stream (except that colors will be discarded). This feature is
 * primarily meant to be used to write a log file.
 */
namespace logging
{

/// \brief Helps typing less when developing the framework.
using stream = boost::iostreams::filtering_ostream;

extern stream debug; ///< \brief Stream for debug messages.
extern stream info;  ///< \brief Stream for information messages.
extern stream warn;  ///< \brief Stream for warning messages.
extern stream error; ///< \brief Stream for reporting errors.
extern stream fatal; ///< \brief Stream for reporting fatal errors.

/// \brief Describes the log levels.
enum class level {
    debug, ///< \brief Display all messages.
    info,  ///< \brief Display all but debug messages.
    warn,  ///< \brief Display all messages except info and debug ones.
    error, ///< \brief Display messages from the error and fatal streams.
    fatal  ///< \brief Display fatal errors only.
};

/// \brief Sets whether the primary stream should have colors.
void set_use_color(bool color);

/// \brief Call this function to pipe ROOT's error messages through the \c logging streams.
void override_root_handler();

/**
 * \brief Changes the *primary* stream, ie the one usually used for terminal output.
 *
 * The stream has to remain valid until the end of the program (or another stream is set).
 *
 * \see \ref set_primary_level
 */
void set_primary_stream(std::ostream &stream = std::cerr);

/**
 * \brief Changes the log level of the *primary* stream.
 * \see \ref set_primary_stream
 */
void set_primary_level(level l);

/**
 * \brief Changes the *secondary* stream, ie the one usually used for file output.
 *
 * The stream has to remain valid until the end of the program (or another stream is set).
 *
 * \see \ref set_secondary_level
 */
void set_secondary_stream(std::ostream &stream);

/**
 * \brief Removes any secondary stream
 *
 * It is safe to call this function when no secondary stream is set.
 */
void unset_secondary_stream();

/**
 * \brief Changes the log level of the *secondary* stream.
 * \see \ref set_secondary_stream
 */
void set_secondary_level(level l);
} // namespace logging
} // namespace util

#endif // LOGGING_H
