module swap PrgEnv-pgi PrgEnv-gnu
module load cray-hdf5
module load cmake
module load python
module load python_scipy
module load python_numpy
module load python_h5py
export HDF5_HOME="/opt/cray/hdf5/1.8.13/GNU/48"
export MKL_ROOT="/opt/intel/composer_xe_2013_sp1.2.144/mkl"
