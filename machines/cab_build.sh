# Load modules
source cab_modules.sh

# Get the repository
git clone https://github.com/etano/simpimc.git
cd simpimc/
git checkout dev
git submodule update --init --recursive\

# Start building
mkdir build && cd build/
cmake ..
make -j

# Resume make
make -j

# Setup pair action generation things
cd ../scripts/pagen
sh setup.sh
