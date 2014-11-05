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
SET(BUILD_STATIC FALSE)
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

SET(COMMON_FLAGS "-O3 -ffast-math")
IF(DEBUG)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -pg")
ENDIF(DEBUG)

IF(USE_MPI)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -DUSE_MPI")
ENDIF(USE_MPI)
IF(USE_OPENMP)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -DUSE_OPENMP")
ENDIF(USE_OPENMP)

SET(CMAKE_CXX_FLAGS "-std=c++11 -Wall ${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")
SET(CMAKE_CXX_SOURCE_FILE_EXTENSIONS ,"cpp")

SET(CMAKE_C_FLAGS "-std=c99 ${CMAKE_C_FLAGS} ${COMMON_FLAGS}")
SET(CMAKE_C_SOURCE_FILE_EXTENSIONS ,"c")

IF(BUILD_STATIC)
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
ENDIF(BUILD_STATIC)

#-------------------------------------------------------------------
# Directories
#-------------------------------------------------------------------
SET(LIBRARY_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/lib" CACHE PATH "Single output directory for building all libraries")
SET(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin" CACHE PATH "Single output directory for building all executables")
SET(CMAKE_INSTALL_PREFIX "$ENV{HOME}" CACHE PATH "Default install to home directory")
SET(PROJECT_CMAKE "${CMAKE_SOURCE_DIR}/CMake")
SET(BUILD_DIR "${CMAKE_SOURCE_DIR}/build")

#-------------------------------------------------------------------
# Scaffold
#-------------------------------------------------------------------
SET(SCAFFOLD_DIR "${CMAKE_SOURCE_DIR}/src/Utils")
SET(SCAFFOLD_DEPENDS_DIR "${CMAKE_SOURCE_DIR}/depends")
SET(SCAFFOLD_INSTALL_DIR "${CMAKE_SOURCE_DIR}")
INCLUDE("${SCAFFOLD_DIR}/CMakeLists.txt")

#-------------------------------------------------------------------
# Sources and target
#-------------------------------------------------------------------
# simpimc
SET (SRCS ${SRCS} src/main.cpp
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
  src/Observables/PairCorrelationClass.cpp
  src/Observables/PathDumpClass.cpp
  src/Observables/PermutationClass.cpp
  src/Observables/SignClass.cpp
  src/Observables/StructureFactorClass.cpp
  src/Observables/TimeClass.cpp)
ADD_EXECUTABLE(simpimc ${SRCS})
TARGET_LINK_LIBRARIES(simpimc ${LIBS})
INSTALL(TARGETS simpimc DESTINATION $ENV{HOME}/bin)

# Stats library
ExternalProject_Add(
  stats
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/scripts
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_COMMAND ${CMAKE_CXX_COMPILER} -Wall -O3 -ffast-math -shared Stats.cpp -o stats.so
  INSTALL_COMMAND ""
)

