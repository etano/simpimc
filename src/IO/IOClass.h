#ifndef IO_H
#define IO_H

#include "../config.h"
#include "IO_Datatype.h"

class IOClass
{
public:
  H5::H5File file;
  inline void load(H5::H5File tmpFile) { file = tmpFile; }

  /// Gets for HDF5 traits
  // Get type
  template <class T>
  inline H5::PredType GetHDF5Datatype(T &val) { return IO::hdf5_type_traits<T>::get_type(val); }
  // Get address of first element
  template <class T>
  inline void* GetHDF5Addr(T &val) { return IO::hdf5_type_traits<T>::get_addr(val); }
  // Get size of whole object
  template <class T>
  inline size_t GetHDF5Size(T &val) { return IO::hdf5_type_traits<T>::get_size(val); }
  // Get constant pointer to shape array
  template <class T>
  inline const hsize_t* GetHDF5Shape(T &val) { return IO::hdf5_type_traits<T>::get_shape(val); }
  // Get rank of object
  template <class T>
  inline const int GetHDF5Rank(T &val) { return IO::hdf5_type_traits<T>::get_rank(val); }

  // Read
  template<class T>
  inline void Read(const std::string& dataset_name, T& data)
  {
    H5::DataSet dataset = file.openDataSet(dataset_name);
    H5::DataSpace dataspace = dataset.getSpace();
    H5::PredType datatype = GetHDF5Datatype(data);
    dataset.read(GetHDF5Addr(data), datatype, dataspace, dataspace);
  }

  // Write
  template<class T>
  inline void Write(const std::string& dataset_name, T& data)
  {
    H5::PredType datatype = GetHDF5Datatype(data);
    datatype.setOrder(H5T_ORDER_LE);
    H5::DataSpace dataspace(GetHDF5Rank(data), GetHDF5Shape(data));
    H5::DataSet dataset = file.createDataSet(dataset_name, datatype, dataspace);
    dataset.write(GetHDF5Addr(data), datatype);
  }

  // Create Group
  inline void CreateGroup(const std::string& group_name)
  {
    H5::Group group = file.createGroup(group_name);
  }

  // Create extendable dataset
  template<class T>
  inline void CreateExtendableDataSet(const std::string& prefix, const std::string& dataset_name, T& data)
  {
    // Get data information
    int data_rank = GetHDF5Rank(data);
    const hsize_t* data_shape = GetHDF5Shape(data);
    H5::PredType datatype = GetHDF5Datatype(data);

    // Create the data space with one unlimited dimension.
    int rank = data_rank + 1;
    hsize_t dims[rank];
    hsize_t maxdims[rank];
    dims[0] = 1;
    maxdims[0] = H5S_UNLIMITED;
    for (int i=1; i<rank; i++) {
      dims[i] = data_shape[i-1];
      maxdims[i] = data_shape[i-1];
    }
    H5::DataSpace mspace(rank, dims, maxdims);

    // Modify dataset creation properties, i.e. enable chunking.
    H5::DSetCreatPropList cparms;
    const hsize_t* chunk_dims = dims; // Default chunk size to 1 x shape(data)
    cparms.setChunk(rank, chunk_dims);

    // Set fill value for the dataset
    int fill_val = 0;
    cparms.setFillValue(datatype, &fill_val);

    // Create a new dataset within the file using cparms creation properties.
    std::string full_name = prefix + dataset_name;
    H5::DataSet dataset = file.createDataSet(full_name, datatype, mspace, cparms);

    // Extend the dataset.
    dataset.extend(dims);

    // Select a hyperslab.
    H5::DataSpace fspace = dataset.getSpace();
    hsize_t offset[rank];
    for (int i=0; i<rank; i++)
      offset[i] = 0;
    fspace.selectHyperslab(H5S_SELECT_SET, dims, offset);

    // Write original data into hyperslab
    dataset.write(GetHDF5Addr(data), datatype, mspace, fspace);
  }

  // Extend dataset
  template<class T>
  inline void AppendDataSet(const std::string& prefix, const std::string& dataset_name, T& data)
  {
    // Get data information
    int data_rank = GetHDF5Rank(data);
    const hsize_t* data_shape = GetHDF5Shape(data);
    H5::PredType datatype = GetHDF5Datatype(data);

    // Open data set
    std::string full_name = prefix + dataset_name;
    H5::DataSet dataset = file.openDataSet(full_name);

    // Get old dataspace properties
    H5::DataSpace fspace = dataset.getSpace();
    int rank = data_rank + 1;
    hsize_t dims_old[rank], maxdims[rank];
    fspace.getSimpleExtentDims(dims_old, maxdims);

    // Extend the dataset.
    hsize_t dims_new[rank];
    dims_new[0] = dims_old[0]+1;
    for (int i=1; i<rank; i++)
      dims_new[i] = data_shape[i-1];
    dataset.extend(dims_new);

    // Select a hyperslab.
    fspace = dataset.getSpace();
    hsize_t offset[rank], dims_orig[rank];
    offset[0] = dims_old[0];
    dims_orig[0] = 1;
    for (int i=1; i<rank; i++) {
      offset[i] = 0;
      dims_orig[i] = data_shape[i-1];
    }
    fspace.selectHyperslab(H5S_SELECT_SET, dims_orig, offset);

    // Define memory space
    H5::DataSpace mspace(rank, dims_orig);

    // Write the data to the hyperslab.
    dataset.write(GetHDF5Addr(data), datatype, mspace, fspace);
  }

};
#endif
