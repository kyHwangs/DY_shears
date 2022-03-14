#include "logging.h"

#include <cstdio> // fileno, stdin
#include <fstream>
#include <iostream>
#include <memory>
#include <regex>
#include <sstream>
#include <string>

#include <unistd.h> // isatty

#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/filter/line.hpp>
#include <boost/iostreams/tee.hpp>
#include <boost/ref.hpp>

#include <TError.h>

#include <yaml-cpp/yaml.h>

#include "ansi_seq.h"

namespace util
{

namespace logging
{

/// \cond
// Set streams to non-null sinks in case the user doesn't call init()
stream debug(boost::ref(std::cerr));
stream info(boost::ref(std::cerr));
stream warn(boost::ref(std::cerr));
stream error(boost::ref(std::cerr));
stream fatal(boost::ref(std::cerr));
/// \endcond

namespace /* anonymous */
{

class logging_filter : public boost::iostreams::line_filter
{
    std::string _prefix;
    bool _color;

  public:
    explicit logging_filter(const std::string &prefix, bool color) : _prefix(prefix), _color(color)
    {
    }

  private:
    std::string do_filter(const std::string &line)
    {
        std::string composed_line = "[" + _prefix + "] " + line + ansi::clear_line();
        if (_color) {
            return composed_line;
        } else {
            // Erase any ANSI sequence
            const static std::regex ansi_seq("\x1B\\[[0-?]*[ -/]*[@-~]");
            return std::regex_replace(composed_line, ansi_seq, "");
        }
    }
};

std::string level_string(level l)
{
    using namespace ansi;

    switch (l) {
    case level::debug:
        return "DEBUG";
    case level::info:
        return setcolor(cyan) + "INFO " + reset();
    case level::warn:
        return setcolor(bright_magenta) + "WARN " + reset();
    case level::error:
        return setcolor(bright_red) + "ERROR" + reset();
    case level::fatal:
        return setcolor(bright_yellow) + setcolor(red, background) + "FATAL" + reset();
    }
    return "";
}

void push_chain(
    stream &log_stream, level stream_level, level min_level, bool color, std::ostream &out)
{
    if (stream_level >= min_level) {
        log_stream.push(logging_filter(level_string(stream_level), color));
        log_stream.push(boost::ref(out));
    } else {
        log_stream.push(boost::iostreams::null_sink());
    }
}

void root_error_handler(int level, bool abort, const char *location, const char *msg)
{
    if (level < kInfo) {
        debug << location << ": " << msg << std::endl;
    } else if (level < kWarning) {
        info << location << ": " << msg << std::endl;
    } else if (level < kError) {
        warn << location << ": " << msg << std::endl;
    } else if (level < kBreak) {
        error << location << ": " << msg << std::endl;
    } else {
        fatal << location << ": " << msg << std::endl;
    }
    if (abort) {
        fatal << "ROOT requested stopping the execution." << std::endl;
        std::abort();
    }
}

class stream_settings
{
    logging::stream &stream;
    logging::stream *secondary_stream;
    level stream_level;

  public:
    bool color = false;
    level primary_level = level::info;
    level secondary_level = level::debug;
    std::ostream *primary_ostream = &std::cerr; // Not owned
    std::ostream *secondary_ostream = nullptr;  // Not owned

    explicit stream_settings(logging::stream &stream, logging::level level)
        : stream(stream), stream_level(level)
    {
        update();
    }

    void update()
    {
        stream.reset();
        if (secondary_stream != nullptr) {
            secondary_stream->reset();
            delete secondary_stream;
            secondary_stream = nullptr;
        }

        if (secondary_ostream != nullptr) {
            // Setup secondary stream
            secondary_stream = new logging::stream;
            push_chain(*secondary_stream,
                       stream_level,
                       secondary_level,
                       false, // no color
                       *secondary_ostream);
            stream.push(tee(boost::ref(*secondary_stream)));
        }
        // Setup primary stream
        push_chain(stream, stream_level, primary_level, color, *primary_ostream);
    }

    ~stream_settings()
    {
        stream.reset();
        stream.push(boost::ref(std::cout));
        if (secondary_stream != nullptr) {
            secondary_stream->reset();
            delete secondary_stream;
            secondary_stream = nullptr;
        }
    }
};

stream_settings debug_settings(debug, level::debug);
stream_settings info_settings(info, level::info);
stream_settings warn_settings(warn, level::warn);
stream_settings error_settings(error, level::error);
stream_settings fatal_settings(fatal, level::fatal);
}

void set_use_color(bool color)
{
    for (auto s :
         {&debug_settings, &info_settings, &warn_settings, &error_settings, &fatal_settings}) {
        s->color = color;
        s->update();
    }
}

void override_root_handler() { SetErrorHandler(root_error_handler); }

void set_primary_stream(std::ostream &stream)
{
    for (auto s :
         {&debug_settings, &info_settings, &warn_settings, &error_settings, &fatal_settings}) {
        s->primary_ostream = &stream;
        s->update();
    }
}

void set_primary_level(level l)
{
    for (auto s :
         {&debug_settings, &info_settings, &warn_settings, &error_settings, &fatal_settings}) {
        s->primary_level = l;
        s->update();
    }
}

void set_secondary_stream(std::ostream &stream)
{
    for (auto s :
         {&debug_settings, &info_settings, &warn_settings, &error_settings, &fatal_settings}) {
        s->secondary_ostream = &stream;
        s->update();
    }
}

void unset_secondary_stream()
{
    debug << "Unsetting secondary stream" << std::endl;
    for (auto s :
         {&debug_settings, &info_settings, &warn_settings, &error_settings, &fatal_settings}) {
        s->secondary_ostream = nullptr;
        s->update();
    }
}

void set_secondary_level(level l)
{
    for (auto s :
         {&debug_settings, &info_settings, &warn_settings, &error_settings, &fatal_settings}) {
        s->secondary_level = l;
        s->update();
    }
}
} // namespace logging
} // namespace ansi
