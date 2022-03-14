\mainpage Introduction

## Quick start

Clone the shears project to your working directory, and move to the analysis
folder:

~~~{.sh}
git clone ssh://git@gitlab.cern.ch:7999/shears/shears.git
cd shears/DYJets
~~~

The code relies on recent versions of the C++ standard and ROOT libraries. They
are distributed by the [CERN LCG](http://lcginfo.cern.ch/) environment. A
supported version can be activated using:

~~~{.sh}
source lcg-env.sh
~~~

You're free to try any other environment (CMSSW environments are known NOT to
work).

The code relies on [CMake](https://cmake.org/) for building (building
out-of-source is *not* supported). The initial compilation is done as follows:

~~~{.sh}
cmake . -DCMAKE_BUILD_TYPE=RelWithDebInfo
make -j$(nproc)
~~~

See below for more information about CMake options.

After the first compilation, the code can be rebuilt simply by running `make`
again:

~~~{.sh}
make -j$(nproc)
~~~

You're now ready to run the code. Programs relevant to the Drell-Yan analysis are located in the `Main`
directory and start with `dyjets-`. Use the `--help` option for usage information.

## Build options

There are some options that change how things are built. They can be passed to `cmake` after the dot:

~~~{.sh}
cmake3 . -DOPTION=VALUE
~~~

Here is a list of the most important options:

| Option             | Contents
|:-------------------|:---------
| `DEBUG_PRINTOUT`   | Enable output for event-by-event debugging.
| `CMAKE_BUILD_TYPE` | See the [documentation](https://cmake.org/cmake/help/latest/variable/CMAKE_BUILD_TYPE.html). `RelWithDebInfo` is the recommended setting.

## Pointers to relevant documentation

- The \ref architecture "Architectural Overview" documents the code structure
- The \ref coding-style "Coding Style" page contains some information that any coder should know
- There is also a page about \ref options "Options Handling"
- It's always useful to have pointers to the lists of [namespaces](namespaces.html) and
  [classes](annotated.html).
