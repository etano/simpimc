# Load modules
source titan_modules.sh

# Get repository
git clone https://github.com/etano/simpimc.git
cd simpimc
git checkout dev
git submodule update --init --recursive
# Edit CMakeLists.txt to build statically
# Change mpicc and mpicxx to cc and CC, respectively
mkdir build && cd build/
cmake ..
make
# Had to change how mkl links in the links.txt file
make
