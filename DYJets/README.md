This folder contains the code used for Drell-Yan analyses. In addition, it has a
stub Higgs (to four leptons) analysis. A fair amount of documentation is
automatically extracted from the source code, and can be found
[here](http://homepage.iihe.ac.be/~lmoureau/shears/ci/).

Quick start
============

Clone the shears project to your working directory, and move to the analysis
folder:

    $ git clone ssh://git@gitlab.cern.ch:7999/shears/shears.git
    $ cd shears/DYJets

The code relies on recent versions of the C++ standard and ROOT libraries. They
are distributed by the [CERN LCG](http://lcginfo.cern.ch/) environment. A
supported version can be activated using:

    $ source lcg-env.sh

You're free to try any other environment.

The code relies on [CMake](https://cmake.org/) for building (building
out-of-source is *not* supported). The initial compilation is done as follows:

    $ cmake . -DCMAKE_BUILD_TYPE=RelWithDebInfo
    $ make -j$(nproc)

See `Documentation/Introduction.md` for more information about CMake options.

After the first compilation, the code can be rebuilt simply by running `make`
again:

    $ make -j$(nproc)

Usage
=====

Coming "soon", in the meantime take a look at the latest version of the
[Doxygen documentation](http://homepage.iihe.ac.be/~lmoureau/shears/ci/).
