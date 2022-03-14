#include "options.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <thread>

#include <TH1.h>

#include "ansi_seq.h"
#include "logging.h"
#include "timer.h"

namespace util
{

options::~options() { logging::unset_secondary_stream(); }

void options::add_defaults(const std::string &default_config_file)
{
    _config_file = default_config_file;
    po::options_description general = po::options_description("General options");
    general.add_options()("config,c",
                          po::value<std::string>()->default_value(default_config_file),
                          "Sets the configuration file to use");
    general.add_options()("verbose,v", "Print a lot of debugging information");
    general.add_options()("help,h", "Prints this help message");
    _all.add(general);
}

void options::print_usage()
{
    std::cerr << "Usage: " << _prog_name << " [options...]" << std::endl;
    std::cerr << _all << std::endl;
}

void options::parse_command_line(int argc, char **argv)
{
    _prog_name = argv[0];
    try {
        po::options_description with_egg = _all;
        with_egg.add_options()("make-paper", "Make analysis note");
        po::store(po::parse_command_line(argc, argv, with_egg), map);
        po::notify(map);
    } catch (po::error &e) {
        std::cerr << "Error: " << e.what() << "." << std::endl;
        print_usage();
        std::exit(EXIT_FAILURE);
    }
}

void options::process_help()
{
    if (map.count("help") > 0) {
        print_usage();
        std::exit(EXIT_SUCCESS);
    }
}

void options::process_config()
{
    if (map.count("config") > 0) {
        _config_file = map["config"].as<std::string>();
    }
    config = YAML::LoadFile(config_file());
}

namespace /* anonymous */
{
logging::level str_to_level(const std::string &str)
{
    if (str == "debug") {
        return logging::level::debug;
    } else if (str == "info") {
        return logging::level::info;
    } else if (str == "warn") {
        return logging::level::warn;
    } else if (str == "error") {
        return logging::level::error;
    } else if (str == "fatal") {
        return logging::level::fatal;
    } else {
        throw std::runtime_error("Invalid log level: " + str);
    }
}
} // namespace anonymous

void options::setup_logging()
{
    // Configuration file
    YAML::Node node = config["log"];
    if (node["color"]) {
        try {
            logging::set_use_color(node["color"].as<bool>());
        } catch (...) {
            std::string strval = node["color"].as<std::string>();
            if (strval == "auto") {
                logging::set_use_color(isatty(fileno(stderr)));
            } else {
                throw std::runtime_error("Invalid color mode: " + strval);
            }
        }
    }
    if (node["log level"]) {
        logging::set_primary_level(str_to_level(node["log level"].as<std::string>()));
    }
    if (node["log file"]) {
        _logfile_out = std::make_shared<std::ofstream>(node["log file"].as<std::string>());
        logging::set_secondary_stream(*_logfile_out);
    }
    if (node["log file level"]) {
        logging::set_secondary_level(str_to_level(node["log file level"].as<std::string>()));
    }
    if ((!node["override root handler"]) || node["override root handler"].as<bool>()) {
        logging::override_root_handler();
    }

    // Command line
    if (map.count("verbose") > 0) {
        logging::set_primary_level(logging::level::debug);
        logging::debug << "Debug messages are enabled." << std::endl;
    }
    if (map.count("verbose") > 0) {
        logging::debug << "Using configuration file: " << _config_file << std::endl;
    }
}

void options::process_easter_egg()
{
    if (map.count("make-paper") > 0) {
        logging::info << "Paper-making mode selected." << std::endl;

        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        logging::info << "Cutting down trees..." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(3000));

        std::string command = "find ~ -name \\*.root ";
        command += "-printf ";
        command += "\"[" + ansi::setcolor(ansi::cyan) + "INFO" + ansi::reset() +
                   " ] Cutting down: %h/%f\\n\" ";
        command += "| head -n50 ";
        std::system(command.c_str());

        logging::warn << "Further files omitted for brevity." << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(10000));

        logging::info << "All trees cut down successfully. Sending them to the paper mill..."
                      << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));

        int messages_shown = 0;

        const unsigned long steps = 314159265358UL; // pi * 10^n
        timer time(steps);
        time.start();
        for (double fraction_remaining = 1.0;
             (unsigned long)(steps * fraction_remaining) > 0; // Integer precision is finite...
             fraction_remaining *= 0.9) {
            time.update(steps - steps * fraction_remaining);
            std::this_thread::sleep_for(std::chrono::milliseconds(200));

            if (fraction_remaining < 1e-3 && messages_shown == 0) {
                logging::warn << "Oxygen level dropping dangerously." << std::endl;
                messages_shown++;
            } else if (fraction_remaining < 1e-6 && messages_shown == 1) {
                logging::warn << "Fatal oxygen level reached. Hope your mask is working."
                              << std::endl;
                messages_shown++;
            }
        }
        time.stop();

        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        logging::info << "Upload successful. Downloading results..." << std::endl;

        std::this_thread::sleep_for(std::chrono::milliseconds(4000));
        // Produce a real segfault with an interesting stack trace
        std::copy((int *)(nullptr), (int *)(nullptr) + 1, (int *)(nullptr));

        std::exit(EXIT_FAILURE);
    }
}

void options::setup_root() const { TH1::SetDefaultSumw2(); }

void options::default_init(int argc,
                           char **argv,
                           const std::string &default_config_file,
                           const std::initializer_list<po::options_description> &groups)
{
    add_defaults(default_config_file);
    for (auto g : groups) {
        _all.add(g);
    }
    parse_command_line(argc, argv);
    process_help();
    process_config();
    setup_logging();
    process_easter_egg();
    setup_root();
}

/// \cond
template <> std::string name_for<po::variables_map>(const std::string &name)
{
    std::string res = name;
    boost::replace_all(res, " ", "-");
    return res;
}

template <> bool is_present<YAML::Node>(const std::string &name, const YAML::Node &node)
{
    return node.IsMap() && node[name];
}

template <>
bool is_present<po::variables_map>(const std::string &name, const po::variables_map &varmap)
{
    return varmap.count(name) > 0;
}
/// \endcond

} // namespace util
