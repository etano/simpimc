# simpimc

Simple path integral monte carlo in c++.

## Installation

### Dependencies:

* Armadillo c++ - http://arma.sourceforge.net/
* Meinspline - https://github.com/etano/meinspline

Be sure to set the environmental variables:

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
