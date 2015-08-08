# Load modules
source comet_modules.sh

# Get the repository
git clone https://github.com/etano/simpimc.git
cd simpimc/
git checkout dev
git submodule update --init --recursive

# Start building
mkdir build && cd build/
cmake .. -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_C_COMPILER=mpicc
make

# Meinspline installation breaks because of pkgconfig
vim ../depends/meinspline/meinspline-master/configure
# Comment out PKG_CHECK_MODULE lines 22870 and 22871

# Resume make
make

# Setup pair action generation things
cd ../scripts/pagen
sh setup.sh
