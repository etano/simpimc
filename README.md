# simpimc

Simple path integral monte carlo in c++.

## Installation

### Dependencies:

* HDF5 (with c++ API) - http://www.hdfgroup.org/HDF5/
* Armadillo c++ - http://arma.sourceforge.net/
* Meinspline - https://github.com/etano/meinspline

Be sure to set the environmental variables:

    export HDF5_HOME=/armadillo/install/directory
    export ARMA_HOME=/armadillo/install/directory
    export EINSPLINE_HOME=/einspline/install/directory
    
### Compiling

Build like:

    git clone https://github.com/etano/simpimc.git
    cd simpimc
    git submodule init
    git submodule update
    mkdir build && cd build
    cmake ..
    make install
