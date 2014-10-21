# simpimc

Simple path integral monte carlo in c++.

## Installation

### Dependencies:

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
    git submodule init
    git submodule update
    mkdir build && cd build
    cmake ..
    make
