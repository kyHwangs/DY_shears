\page options Options Handling

Options can come from two different places: the command line or the configuration file. In general,
options coming from the command-line supersede those in the config file, but not all options are
available on the command-line. Options are handled by the \ref util::options class, usually passed
as a const reference named `opt`:

~~~{.cc}
void function(const util::options &opt)
{
    ...
}
~~~

The `opt` object exposes two variables: a list of all command-line options called
\ref util::options::map "map", and a variable that represents the contents of the configuration
file, \ref util::options::config "config". These variables will be called "containers" in what
follows.

The `util` namespace contains \ref util::options "some useful functions" that operate on
containers:

- `util::is_present` to check if a container has a value for the given option
- `util::set_value_safe` to set a value only if it is present in the config file

They are used as follows:

~~~{.cc}
void function(const util::options &opt)
{
    if (util::is_present("A", opt.config)) {
        util::logging::info << "Option 'A' is present in the config file." << std::endl;
    }
    if (util::is_present("A", opt.map)) {
        util::logging::info << "Option 'A' was given on the command line." << std::endl;
    }

    std::string B = "default";
    // Set B to the value given in the config file, if any
    util::set_value_safe(opt.config, B, "B", "option 'B'");
    // Set B to the value given on the command line, if any
    util::set_value_safe(opt.map, B, "B", "option 'B'");
}
~~~

The last parameter of `set_value_safe` is a name used in case the parameter value is invalid.

There is a second variant of `set_value_safe` that takes a function as the last argument. The
function is used to check validity of the value provided by the user. As this is a perfect
use-case for C++ lambda expressions, please go and learn them if you want input validation.

## Command line specifics

Command line options have to be declared before being used. This is how the program can display a
nice list of allowed options when called with `-h` or `--help`. You'll find examples of option
declarations in \ref util::options::add_defaults.

## Configuration file specifics

The configuration file has more structure than command line options: it handles lists and maps of
arbitrary objects. The key point to remember about the config file is that every value (be it a
key, a value, a map or a list) can be represented as a `YAML::Node`, and that `YAML::Node` can be
used as a container (`opt.config` itself is in fact a `YAML::Node`).

`YAML::Node` has built-in support for turning lists (in the config file) into `std::vector`s, and
maps into `std::map`s. Error messages when converting to C++ values are rather uninformative, so
it's best to avoid using nested containers. There's also a nice syntax to allow for arbitrary
conversions (i.e. turning a `YAML::Node` into a custom C++ object), using specializations of
`YAML::convert<T>::decode`. This mechanism is used to read e.g. the list of samples.
