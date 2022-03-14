#ifndef OPTIONS_H
#define OPTIONS_H

#include <iostream>
#include <memory>

#include <boost/algorithm/string/replace.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <yaml-cpp/yaml.h>

#include "logging.h"

namespace po = boost::program_options;

/**
 * \brief Generally useful classes and functions.
 *
 * This namespace is meant to group utilities that have little physics interest.
 */
namespace util
{

/**
 * \brief Class to handle program options and startup.
 *
 * Starting up a program is as simple as:
 *
 * ~~~{.cc}
 * try {
 *     options opt;
 *     opt.default_init(argc, argv, "config.yml", {});
 * } catch (std::exception &e) {
 *     logging::fatal << e.what() << std::endl;
 *     return EXIT_FAILURE;
 * }
 * ~~~
 *
 * The lines above will parse the command line and the config file, handling common options such as
 * \c --help. Most errors are reported as exceptions.
 */
class options
{
    std::string _prog_name;
    po::options_description _all;
    std::string _config_file;
    std::shared_ptr<std::ostream> _logfile_out;

  public:
    /// \brief The parsed contents of the configuration file.
    YAML::Node config;

    /// \brief The parsed contents of the command line.
    po::variables_map map;

    /// \brief Constructor.
    explicit options() = default;

    /// \brief Deleted copy constructor.
    explicit options(const options &other) = delete;

    /// \brief Destructor.
    virtual ~options();

    /// \brief Deleted \c operator=.
    options &operator= (const options &other) = delete;

    /**
     * \brief Add default options \c --help, \c --verbose. and \c --config to the list of command
     *        line options.
     */
    void add_defaults(const std::string &default_config_file);

    /// \brief Prints usage information.
    void print_usage();

    /// \brief Parses the command line.
    void parse_command_line(int argc, char **argv);

    /// \brief Processes the \c --help options.
    void process_help();

    /// \brief Processes the \c --config options.
    void process_config();

    /// \brief Sets up the \ref logging module according to the config file and \c --verbose.
    void setup_logging();

    /// \brief Sets up ROOT options like \c SetDefaultSumw2.
    void setup_root() const;

    /**
     * \brief Default init sequence.
     *
     * \param argc The number of arguments from the command line.
     * \param argv The arguments from the command line.
     * \param default_config_file The name of the default config file.
     * \param groups Additionnal command line option groups.
     */
    void default_init(int argc,
                      char **argv,
                      const std::string &default_config_file,
                      const std::initializer_list<po::options_description> &groups);

    /// \brief Retrieves the name of the config file that was used.
    std::string config_file() const { return _config_file; }

  private:
    /// \brief \c --make-paper
    void process_easter_egg();
};

/**
 * \brief Transforms the name to match the convention for the given container.
 *
 * The default implementation returns its argument unchanged.
 *
 * \relates options
 */
template <class Container> std::string name_for(const std::string &name) { return name; }

/**
 * \brief Transforms the name to match the convention for the command line.
 *
 * This specialization turns whitespaces into dashes.
 *
 * \relates options
 */
template <> std::string name_for<po::variables_map>(const std::string &name);

/**
 * \brief Checks whether a value with `name` is present in the given `Container`.
 *
 * This function has specialization for \c YAML::Node and \c po::variables_map.
 *
 * \relates options
 */
template <class Container> bool is_present(const std::string &name, const Container &);

/**
 * \brief Sets the value of \c target according to the contents of \c container.
 *
 * If the value isn't present in the container, this function doesn't do anything. Else, it tries
 * to convert it to the right type, and to validate it using `check`. An exception is thrown in
 * case of error.
 *
 * \param container The container to get the value from.
 * \param target    The variable to store the value into.
 * \param name      The name of the variable. It will be translated using \ref name_for.
 * \param descr     A user-readable string describing the variable. Used for output.
 * \param check     A function that takes a value and returns \c true if it is valid.
 *
 * \relates options
 */
template <class Type, class Container, class Checker>
void set_value_safe(const Container &container,
                    Type &target,
                    const std::string &name,
                    const std::string &descr,
                    const Checker &check)
{
    std::string name_tr = name_for<Container>(name);
    if (is_present(name_tr, container)) {
        Type config_value = container[name_tr].template as<Type>();
        if (check(config_value)) {
            logging::debug << "Setting " << descr << " to " << config_value << std::endl;
            target = config_value;
        } else {
            throw std::runtime_error("Invalid " + descr + ": " + std::to_string(config_value));
        }
    }
}

/**
 * \brief Sets the value of \c target according to the contents of \c container.
 *
 * Same as \ref set_value_safe, but doesn't validate the value.
 *
 * \param container The container to get the value from.
 * \param target    The variable to store the value into.
 * \param name      The name of the variable. It will be translated using \ref name_for.
 * \param descr     A user-readable string describing the variable. Used for output.
 *
 * \relates options
 */
template <class Type, class Container>
void set_value_safe(const Container &container,
                    Type &target,
                    const std::string &name,
                    const std::string &descr)
{
    set_value_safe(container, target, name, descr, [](Type) { return true; });
}
} // namespace util

#endif // OPTIONS_H
