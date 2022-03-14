\page architecture Architectural Overview

This page describes the code architecture from a high-level perspective. It is designed to give
you an idea of where to find the code handling a particular part of the analysis.

## Where to find files

The code is split in three:

- The `Main` directory contains programs (i.e. files with a `main` function).
- The `Includes` directory contains all header files (`*.h`).
- The `Sources` directory contains all source files (`*.cc`).

Files in the `Includes` and `Sources` directory may be used by several programs.

The main configuration files (`dyjets*.yml`, `dyjets-binnings.yml`) are kept at the root
of the tree (i.e. in the `DYJets` directory). Other, less important, directories exist, listed in
the summary below.

| Directory          | Contents
|:-------------------|:---------
| `Documentation`    | Markdown files used to build this documentation.
| `EfficiencyTables` | Scale factor tables
| `Includes`         | Header (`.h`) files
| `Main`             | Source (`.cc`) files used to build programs
| `RooUnfold`        | A (modified) copy of the `RooUnfold` library
| `Sources`          | Source (`.cc`) files

## Code architecture

The architecture was designed to make it easy to implement new analyzers: ideally, one would just
get the list of particles in the final state, combine them in the right way and fill histograms.
Technical stuff such as looping over events, applying scale factors or normalizing MC should be
implemented in a way that makes it easy-to-use.

The above goals are achieved using an object-oriented design, as is natural in C++. The code is
thus split into many classes, each being declared in a `.h` file and implemented in the
corresponding `.cc` file (albeit simple functions are implemented inline, in the `.h`). Every
public class or function is documented using Doxygen markup. The documentation should be sufficient
to understand what a class does and how to use it, but not how it works internally − such
implementation details should be irrelevant.

To help navigating the documentation, classes are grouped in several namespaces according to their
purpose:

- The \ref util namespace groups classes that have no physical interest: job control classes,
  option handling facilities, ...
- The \ref data namespace groups classes that describe datasets.
- The \ref physics namespace is the most interesting for physics.

In addition, some classes inherited from the DYJets code are in the global namespace. They should
eventually be moved to a more appropriate location.

## Dependencies

The code uses a number of select libraries to avoid reinventing the wheel. They are all available
in recent LCG environments. The following libraries are required:

- The [Boost C++ Libraries](http://boost.org), and in particular:
    - The Program Options library to handle command-line options;
    - The Filesystem library to handle filesystem operation such as creating directories;
    - The String Algorithms library for operations on strings.
- The well-known (among CERN physicists) [ROOT](http://root.cern.ch) toolkit (at least version `6.12`).
- [RooUnfold](http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html), used for unfolding. A
  modified copy is kept in the source tree.
- [yaml-cpp](https://github.com/jbeder/yaml-cpp) (version 0.5 or later), used by the `Higgs` code
  to parse the configuration files.

In addition, `cmake` >= 3.8 and a C++ 17 compiler are required. Again, recent LCG environments
should provide these tools.
