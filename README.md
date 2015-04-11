# simpimc

Simple path integral monte carlo in c++.

master: ![master build status](https://travis-ci.org/etano/simpimc.svg?branch=master)

dev: ![dev build status](https://travis-ci.org/etano/simpimc.svg?branch=dev)

## Installation

### Dependencies

* HDF5 (with c++ API) - http://www.hdfgroup.org/HDF5/
* Armadillo c++ - http://arma.sourceforge.net/
* Meinspline - https://github.com/etano/meinspline

These dependencies will automatically be installed by cmake. However, if you'd like to install them yourself, be sure to set the environmental variables:

    export HDF5_HOME=/hdf5/install/directory
    export ARMADILLO_HOME=/armadillo/install/directory
    export EINSPLINE_HOME=/einspline/install/directory
    
### Compiling

Build like:

    git clone https://github.com/etano/simpimc.git
    cd simpimc
    git submodule update --init --recursive
    mkdir build && cd build
    cmake ..
    make

To build the documentation:

    make doc

## Running

To run simpimc simply do the following:

    simpimc input.xml
    
Or if you wish to run in parallel:

    mpiexec -np 2 simpimc input.xml
    
Example XML inputs can be found in the inputs folder.
