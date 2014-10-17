PROJECT(SIMPIMC)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#if you don't want the full compiler output, set to "OFF"
SET(CMAKE_VERBOSE_MAKEFILE OFF)

#-------------------------------------------------------------------
# Options
#-------------------------------------------------------------------

SET(PRECISION double) #double float
SET(USE_MPI TRUE)
SET(USE_OPENMP TRUE)
SET(QMC_BUILD_STATIC FALSE)
SET(NDIM 3)
SET(COMPILER_MAKE "GNU") # INTEL, GNU, or IBM
SET(DEBUG TRUE)

#-------------------------------------------------------------------
# Compilers
#-------------------------------------------------------------------
IF(${COMPILER_MAKE} MATCHES "IBM")
  # set the compiler
  IF(USE_MPI)
    SET(CMAKE_C_COMPILER mpixlc_r)
    SET(CMAKE_CXX_COMPILER mpixlcxx_r)
  ELSE()
    SET(CMAKE_C_COMPILER xlc_r)
    SET(CMAKE_CXX_COMPILER xlc++_r)
  ENDIF()

  SET(AIX_ARCH "qp")
  SET(AIX_ARCH_FLAGS "-qarch=${AIX_ARCH} -qstrict -qhot=level=1 -qtune=qp -qsimd=auto")
  IF(USE_OPENMP)
    SET(AIX_ARCH_FLAGS "${AIX_ARCH_FLAGS} -qsmp=omp -qthreaded")
  ENDIF()

  SET(AIX_CXX_COMMON_FLAGS "-qinline=auto:level=10")
  SET(AIX_OPT_FLAGS "-qmaxmem=-1")

  SET(CMAKE_CXX_FLAGS "${AIX_ARCH_FLAGS} ${AIX_CXX_COMMON_FLAGS} ${AIX_OPT_FLAGS}")
  SET(CMAKE_C_FLAGS "${AIX_ARCH_FLAGS} ${AIX_OPT_FLAGS}")
ENDIF()

IF(${COMPILER_MAKE} MATCHES "INTEL")
  IF(USE_MPI)
    SET(CMAKE_C_COMPILER mpicc)
    SET(CMAKE_CXX_COMPILER mpic++)
  ELSE()
    SET(CMAKE_C_COMPILER icc)
    SET(CMAKE_CXX_COMPILER icpc)
  ENDIF()

  SET(INTEL_C_FLAGS "-w")
  IF(USE_OPENMP)
    SET(INTEL_C_FLAGS "${INTEL_C_FLAGS} -fopenmp")
  ENDIF(USE_OPENMP)

  SET(CMAKE_CXX_FLAGS "${INTEL_C_FLAGS}")
  SET(CMAKE_C_FLAGS "${INTEL_C_FLAGS}")
ENDIF()

IF(${COMPILER_MAKE} MATCHES "GNU")
  IF(USE_MPI)
    SET(CMAKE_C_COMPILER mpicc)
    SET(CMAKE_CXX_COMPILER mpicxx)
  ELSE()
    SET(CMAKE_C_COMPILER $ENV{CC})
    SET(CMAKE_CXX_COMPILER $ENV{CXX})
  ENDIF()

  SET(GNU_C_FLAGS "-w")
  IF(USE_OPENMP)
    SET(GNU_C_FLAGS "${GNU_C_FLAGS} -fopenmp")
  ENDIF(USE_OPENMP)

  SET(CMAKE_CXX_FLAGS "${GNU_C_FLAGS}")
  SET(CMAKE_C_FLAGS "${GNU_C_FLAGS}")
ENDIF()

IF(${COMPILER_MAKE} MATCHES "CLANG")
  IF(USE_MPI)
    SET(CMAKE_C_COMPILER mpicc)
    SET(CMAKE_CXX_COMPILER mpicxx)
  ELSE()
    SET(CMAKE_C_COMPILER clang)
    SET(CMAKE_CXX_COMPILER clang)
  ENDIF()

  SET(CLANG_C_FLAGS "-w")
  IF(USE_OPENMP)
    SET(CLANG_C_FLAGS "${CLANG_C_FLAGS} -fopenmp")
  ENDIF(USE_OPENMP)

  SET(CMAKE_CXX_FLAGS "${CLANG_C_FLAGS}")
  SET(CMAKE_C_FLAGS "${CLANG_C_FLAGS}")
ENDIF()

SET(COMMON_FLAGS "-O3 -ffast-math -DNDIM=${NDIM}")
IF(DEBUG)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -pg")
ELSE(DEBUG)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -DARMA_NO_DEBUG")
ENDIF(DEBUG)

IF(USE_MPI)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -DUSE_MPI")
ENDIF(USE_MPI)
IF(USE_OPENMP)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -DUSE_OPENMP")
ENDIF(USE_OPENMP)

SET(CMAKE_CXX_FLAGS "-std=c++11 -msse4 ${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")
SET(CMAKE_CXX_SOURCE_FILE_EXTENSIONS ,"cpp")

SET(CMAKE_C_FLAGS "-std=c99 -msse4 ${CMAKE_C_FLAGS} ${COMMON_FLAGS}")
SET(CMAKE_C_SOURCE_FILE_EXTENSIONS ,"c")

IF(QMC_BUILD_STATIC)
  # Enable static linking
  set(BUILD_STATIC_LIBS TRUE)
  set(LINK_SEARCH_END_STATIC TRUE)
  set(CMAKE_EXE_LINKER_FLAGS "-static -Wl,-Bstatic")
  set(MPI_LINK_FLAGS "-Wl,-Bstatic")
  set(LINK_SEARCH_END_STATIC TRUE)
  set(HDF5_USE_STATIC_LIBRARIES TRUE)

  # Disable shared linking
  set(BUILD_SHARED_LIBS FALSE)
  set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS)       # remove -Wl,-Bdynamic
  set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS)
  set(CMAKE_SHARED_LIBRARY_C_FLAGS)         # remove -fPIC
  set(CMAKE_SHARED_LIBRARY_CXX_FLAGS)
  set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS)    # remove -rdynamic
  set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS)
ENDIF(QMC_BUILD_STATIC)

#-------------------------------------------------------------------
# Directories
#-------------------------------------------------------------------
SET(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib" CACHE PATH "Single output directory for building all libraries")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Single output directory for building all executables")
SET(CMAKE_INSTALL_PREFIX "$ENV{HOME}" CACHE PATH "Default install to home directory")
SET(PROJECT_CMAKE "${CMAKE_SOURCE_DIR}/CMake")
SET(BUILD_DIR "${CMAKE_SOURCE_DIR}/build")

#-------------------------------------------------------------------
# Dependencies
#-------------------------------------------------------------------

# Find Lapack/Blas/Atlas/MKL/XLmass/ESSL/...
SET(LA_LIBS "")
IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "ppc") # This covers GNU and IBM Compilers on BGQ
  # Find Lapack/BLAS BGQ
  SET(LAPACK_LIBRARY /usr/local/tools/lapack/lib/liblapack.a)
  SET(BLAS_LIBRARY /usr/local/tools/blas/lib/libblas.a)
  SET(ESSL_LIBRARY /usr/local/tools/essl/5.1/lib/libesslsmpbg.a)
  SET(XLMASS_LIBRARY /opt/ibmcmp/xlmass/bg/7.3/bglib64/libmass.a
                     /opt/ibmcmp/xlmass/bg/7.3/bglib64/libmassv.a)
  SET(FORTRAN_LIBRARY /opt/ibmcmp/xlf/bg/14.1/bglib64/libxlf90_r.a
                      /opt/ibmcmp/xlf/bg/14.1/bglib64/libxl.a
                      /opt/ibmcmp/xlf/bg/14.1/bglib64/libxlfmath.a
                      /opt/ibmcmp/xlf/bg/14.1/bglib64/libxlfpmt4.a
                      /opt/ibmcmp/xlf/bg/14.1/bglib64/libxlfpmt8.a
                      /opt/ibmcmp/xlf/bg/14.1/bglib64/libxlopt.a)
  SET(LA_LIBS ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${ESSL_LIBRARY}
    ${XLMASS_LIBRARY} ${FORTRAN_LIBRARY} -ldl)
ELSEIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  SET(CMAKE_CXX_FLAGS "-framework Accelerate ${CMAKE_CXX_FLAGS}")
  SET(CMAKE_C_FLAGS "-framework Accelerate ${CMAKE_C_FLAGS}")
ELSE()
    # Find MKL
    FIND_PATH(MKL_INCLUDE_DIR mkl.h)
    FIND_PATH(MKL_LINK_DIR libmkl_core.a)
    SET(MKL_FOUND false)
    IF(MKL_INCLUDE_DIR)
      SET(MKL_FOUND true)
    ELSE()
      FIND_PATH(MKL_INCLUDE_DIR mkl.h ${MKL_HOME} ${MKL_HOME}/include $ENV{MKL_HOME} $ENV{MKL_HOME}/include)
      FIND_PATH(MKL_LINK_DIR libmkl_core.a ${MKL_HOME} ${MKL_HOME}/lib/intel64  $ENV{MKL_HOME} $ENV{MKL_HOME}/lib/intel64)
      IF(MKL_INCLUDE_DIR)
        SET(MKL_FOUND true)
      ENDIF()
    ENDIF()
    IF(MKL_FOUND)
      SET(MKL_LIBS mkl_intel_lp64 mkl_intel_thread mkl_core iomp5 pthread m mkl_scalapack_lp64 mkl_blacs_intelmpi_lp64)
      #SET(MKL_LIBS "-L$ENV{MKL_HOME};mkl_intel_lp64;mkl_sequential;mkl_core")
      LINK_DIRECTORIES(/opt/intel/lib/intel64)
      LINK_DIRECTORIES(${MKL_LINK_DIR})
      SET(LA_LIBS ${MKL_LIBS})
    ELSE()
      FIND_LIBRARY(LAPACK_LIBS lapack /usr/lib /usr/local/lib)
      FIND_LIBRARY(BLAS_LIBS blas /usr/lib /usr/local/lib)
      IF(LAPACK_LIBS AND BLAS_LIBS)
        SET(LA_LIBS ${LAPACK_LIBS} ${BLAS_LIBS})
        LINK_LIBRARIES(${LA_LIBS})
      ENDIF()
    ENDIF()
ENDIF()
MESSAGE(STATUS "LA_LIBS=${LA_LIBS}")

# Find Armadillo
FIND_PATH(ARMA_INCLUDE_DIR armadillo $ENV{ARMA_HOME}/include)
FIND_LIBRARY(ARMA_LIBS armadillo $ENV{ARMA_HOME}/lib64 $ENV{ARMA_HOME}/lib)
IF(ARMA_INCLUDE_DIR AND ARMA_LIBS)
    INCLUDE_DIRECTORIES(${ARMA_INCLUDE_DIR})
    LINK_LIBRARIES(${ARMA_LIBS})
    MESSAGE(STATUS "ARMA_INCLUDE_DIR=${ARMA_INCLUDE_DIR}")
    MESSAGE(STATUS "ARMA_LIBS=${ARMA_LIBS}")
ELSE()
    MESSAGE(FATAL_ERROR "Set ARMA_HOME")
ENDIF()

# Find HDF5
FIND_PATH(HDF5_INCLUDE_DIR hdf5.h $ENV{HDF5_HOME}/include)
FIND_LIBRARY(HDF5_LIB hdf5 $ENV{HDF5_HOME}/lib $ENV{HDF5_HOME}/lib64)
FIND_LIBRARY(HDF5_CPP_LIB hdf5_cpp $ENV{HDF5_HOME}/lib $ENV{HDF5_HOME}/lib64)
SET(HDF5_LIBS ${HDF5_LIB} ${HDF5_CPP_LIB})
IF(HDF5_INCLUDE_DIR AND HDF5_LIBS)
    INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIR})
    LINK_LIBRARIES(${HDF5_LIBS})
    MESSAGE(STATUS "HDF5_INCLUDE_DIR=${HDF5_INCLUDE_DIR}")
    MESSAGE(STATUS "HDF5_LIBS=${HDF5_LIBS}")
ELSE()
    MESSAGE(FATAL_ERROR "Set HDF5_HOME")
ENDIF()

# Find Einspline
FIND_PATH(EINSPLINE_INCLUDE_DIR einspline $ENV{EINSPLINE_HOME}/include)
FIND_LIBRARY(EINSPLINE_LIBS einspline $ENV{EINSPLINE_HOME}/lib)
IF(EINSPLINE_INCLUDE_DIR AND EINSPLINE_LIBS)
    INCLUDE_DIRECTORIES(${EINSPLINE_INCLUDE_DIR})
    LINK_LIBRARIES(${EINSPLINE_LIBS})
    MESSAGE(STATUS "EINSPLINE_INCLUDE_DIR=${EINSPLINE_INCLUDE_DIR}")
    MESSAGE(STATUS "EINSPLINE_LIBS=${EINSPLINE_LIBS}")
ELSE()
    MESSAGE(FATAL_ERROR "Set EINSPLINE_HOME")
ENDIF()

# Create target
SET (TEMPLATE_SRCS src/main.cpp
  src/AlgorithmClass.cpp
  src/BeadClass.cpp
  src/LoopClass.cpp
  src/PathClass.cpp
  src/InitPaths.cpp
  src/SimulationClass.cpp
  src/SpeciesClass.cpp
  src/WriteClass.cpp
  src/Actions/ActionClass.cpp
  src/Actions/CoulombClass.cpp
  src/Actions/KineticClass.cpp
  src/Actions/PairActionClass.cpp
  src/Actions/NodalClass.cpp
  src/Actions/DavidPairActionClass.cpp
  src/Actions/IlkkaPairActionClass.cpp
  src/Actions/TrapClass.cpp
  src/Moves/MoveClass.cpp
  src/Moves/BisectClass.cpp
  src/Moves/PermBisectClass.cpp
  src/Moves/PermBisectIterativeClass.cpp
  src/Observables/ObservableClass.cpp
  src/Observables/EnergyClass.cpp
  src/Observables/PathDumpClass.cpp
  src/Observables/TimeClass.cpp
  src/Utils/Communication/Communication.cpp
  src/Utils/IO/InputClass.cpp
  src/Utils/IO/xmlParser.cpp
  src/Utils/RNG/RNGClass.cpp)
ADD_EXECUTABLE(simpimc ${TEMPLATE_SRCS})
TARGET_LINK_LIBRARIES(simpimc ${LA_LIBS} ${ARMA_LIBS} ${HDF5_LIBS} ${EINSPLINE_LIBS})

# Install to home directory
INSTALL(TARGETS simpimc DESTINATION $ENV{HOME}/bin)
