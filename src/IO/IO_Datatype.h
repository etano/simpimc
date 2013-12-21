#ifndef IO_DATATYPE_H
#define IO_DATATYPE_H

#include "H5Cpp.h"

namespace IO {

  // Template to retrieve traits of any HDF5 object
  template <class T>
  struct hdf5_type_traits {
    static H5::PredType get_type(T&& val);
    static inline size_t get_size(T& val);
    static inline void* get_addr(T& val);
    static inline const hsize_t* get_shape(T& val);
    static inline const int get_rank(T& val);
  };

  // Specialization of hdf5_type_traits for primitive types
#define PRIMITIVE(Type, H5PredType, H5Type) \
        template<> \
        inline H5::PredType hdf5_type_traits<Type>::get_type(Type&&) { return H5Type; } \
        template<> \
        inline size_t hdf5_type_traits<Type>::get_size(Type&) { return 1; } \
        template<> \
        inline void* hdf5_type_traits<Type>::get_addr(Type& val) { return &val; } \
        template<> \
        inline const hsize_t* hdf5_type_traits<Type>::get_shape(Type& val) { const hsize_t shape[] = {1}; return shape; } \
        template<> \
        inline const int hdf5_type_traits<Type>::get_rank(Type& val) { return 1; }
  PRIMITIVE(char, H5::IntType, H5::PredType::NATIVE_CHAR);
  PRIMITIVE(short, H5::IntType, H5::PredType::NATIVE_SHORT);
  PRIMITIVE(int, H5::IntType, H5::PredType::NATIVE_INT);
  PRIMITIVE(long, H5::IntType, H5::PredType::NATIVE_LONG);
  PRIMITIVE(signed char, H5::IntType, H5::PredType::NATIVE_SCHAR);
  PRIMITIVE(unsigned char, H5::IntType, H5::PredType::NATIVE_UCHAR);
  PRIMITIVE(unsigned short, H5::IntType, H5::PredType::NATIVE_USHORT);
  PRIMITIVE(unsigned int, H5::IntType, H5::PredType::NATIVE_UINT);
  PRIMITIVE(unsigned long, H5::IntType, H5::PredType::NATIVE_ULONG);
  PRIMITIVE(unsigned long long, H5::IntType, H5::PredType::NATIVE_ULLONG);
  PRIMITIVE(bool, H5::IntType, H5::PredType::NATIVE_HBOOL);
  PRIMITIVE(std::complex<double>, H5::FloatType, H5::PredType::NATIVE_DOUBLE);
  PRIMITIVE(std::complex<long double>, H5::FloatType, H5::PredType::NATIVE_FLOAT);
#if PRECISION==double
  PRIMITIVE(RealType, H5::FloatType, H5::PredType::NATIVE_DOUBLE);
#elif PRECISION==single
  PRIMITIVE(RealType, H5::FloatType, H5::PredType::NATIVE_FLOAT);
#endif
#undef PRIMITIVE

  // Specialization of hdf5_type_traits for armadillo types
#define ARMATYPE(Type, ElemType, H5PredType, H5Type) \
        template<> \
        inline H5::PredType hdf5_type_traits<Type>::get_type(Type&&) { return H5Type; } \
        template<> \
        inline size_t hdf5_type_traits<Type>::get_size(Type& val) { return val.size(); } \
        template<> \
        inline void* hdf5_type_traits<Type>::get_addr(Type& val) { return val.memptr(); } \
        template<> \
        inline const hsize_t* hdf5_type_traits<Type>::get_shape(Type& val) { const hsize_t shape[] = {val.n_cols,val.n_rows}; return shape; } \
        template<> \
        inline const int hdf5_type_traits<Type>::get_rank(Type& val) { return 2; }
  ARMATYPE(Imatrix, int, H5::IntType, H5::PredType::NATIVE_INT);
  ARMATYPE(Ivector, int, H5::IntType, H5::PredType::NATIVE_INT);
#if PRECISION==double
  ARMATYPE(Tmatrix, RealType, H5::FloatType, H5::PredType::NATIVE_DOUBLE);
  ARMATYPE(Tvector, RealType, H5::FloatType, H5::PredType::NATIVE_DOUBLE);
#elif PRECISION==single
  ARMATYPE(Tmatrix, RealType, H5::FloatType, H5::PredType::NATIVE_FLOAT);
  ARMATYPE(Tvector, RealType, H5::FloatType, H5::PredType::NATIVE_FLOAT);
#endif
#undef ARMATYPE

}

#endif
