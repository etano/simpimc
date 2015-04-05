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
SET(COMPILER_MAKE "GNU") # INTEL, GNU, or CLANG
SET(SCAFFOLD_MATRIX_LIBRARY "ARMADILLO") # ARMADILLO or EIGEN
SET(OPTIMIZE TRUE)
SET(DEBUG TRUE)

#-------------------------------------------------------------------
# Compilers
#-------------------------------------------------------------------
IF(${COMPILER_MAKE} MATCHES "INTEL")
  IF(USE_MPI)
    SET(CMAKE_C_COMPILER mpiicc)
    SET(CMAKE_CXX_COMPILER mpiicpc)
  ELSE()
    SET(CMAKE_C_COMPILER icc)
    SET(CMAKE_CXX_COMPILER icpc)
  ENDIF()
ENDIF()

IF(${COMPILER_MAKE} MATCHES "GNU")
  IF(USE_MPI)
    SET(CMAKE_C_COMPILER mpicc)
    SET(CMAKE_CXX_COMPILER mpicxx)
  ELSE()
    SET(CMAKE_C_COMPILER gcc)
    SET(CMAKE_CXX_COMPILER g++)
  ENDIF()
ENDIF()

IF(${COMPILER_MAKE} MATCHES "CLANG")
  IF(USE_MPI)
    SET(CMAKE_C_COMPILER mpicc)
    SET(CMAKE_CXX_COMPILER mpicxx)
  ELSE()
    SET(CMAKE_C_COMPILER clang)
    SET(CMAKE_CXX_COMPILER clang++)
  ENDIF()
ENDIF()

#-------------------------------------------------------------------
# Flags
#-------------------------------------------------------------------
SET(COMMON_FLAGS "-w")

IF(DEBUG)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -pg -g -Wall")
ELSE(DEBUG)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -DNDEBUG")
ENDIF(DEBUG)

IF(OPTIMIZE)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -O3 -ffast-math")
ENDIF(OPTIMIZE)

IF(USE_MPI)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -DUSE_MPI -DMPICH_IGNORE_CXX_SEEK")
ENDIF(USE_MPI)

IF(USE_OPENMP)
  SET(COMMON_FLAGS "${COMMON_FLAGS} -fopenmp -DUSE_OPENMP")
ENDIF(USE_OPENMP)

SET(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")
SET(CMAKE_CXX_SOURCE_FILE_EXTENSIONS ,"cc")

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
SET(SCAFFOLD_DIR "${CMAKE_SOURCE_DIR}/src/scaffold")
SET(SCAFFOLD_DEPENDS_DIR "${CMAKE_SOURCE_DIR}/depends")
SET(SCAFFOLD_INSTALL_DIR "${CMAKE_SOURCE_DIR}")
INCLUDE("${SCAFFOLD_DIR}/CMakeLists.txt")
SET(LIBS ${LIBS} ${SCAFFOLD_LIBS})
SET(SRCS ${SRCS} ${SCAFFOLD_SRCS})

#-------------------------------------------------------------------
# Einspline
#-------------------------------------------------------------------
INCLUDE(CMake/FindSSE.cmake)
INCLUDE(CMake/FindEinspline.cmake)
SET(LIBS ${LIBS} ${EINSPLINE_LIBS})

#-------------------------------------------------------------------
# Sources and target
#-------------------------------------------------------------------
# simpimc
SET (SRCS ${SRCS} src/main.cc
  src/algorithm_class.cc
  src/loop_class.cc
  src/path_class.cc
  src/simulation_class.cc
  src/species_class.cc
  src/write_class.cc
  src/actions/kinetic_class.cc
  src/actions/pair_action_class.cc
  src/actions/nodal_class.cc
  src/actions/free_nodal_class.cc
  src/actions/optimized_nodal_class.cc
  src/actions/bare_pair_action_class.cc
  src/actions/david_pair_action_class.cc
  src/actions/ilkka_pair_action_class.cc
  src/actions/importance_pair_action_class.cc
  src/actions/trap_class.cc
  src/moves/move_class.cc
  src/moves/bisect_class.cc
  src/moves/displace_particle_class.cc
  src/moves/perm_bisect_class.cc
  src/moves/perm_bisect_iterative_class.cc
  src/moves/shift_ref_slice_class.cc
  src/moves/vary_optimized_nodal_class.cc
  src/observables/contact_density_class.cc
  src/observables/energy_class.cc
  src/observables/importance_weight_class.cc
  src/observables/pair_correlation_class.cc
  src/observables/path_dump_class.cc
  src/observables/permutation_class.cc
  src/observables/record_optimized_nodal_class.cc
  src/observables/sign_class.cc
  src/observables/structure_factor_class.cc
  src/observables/time_class.cc)
ADD_EXECUTABLE(simpimc ${SRCS})
TARGET_LINK_LIBRARIES(simpimc ${LIBS})
INSTALL(TARGETS simpimc DESTINATION $ENV{HOME}/bin)

# Stats library
ExternalProject_Add(
  stats
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/scripts
  CONFIGURE_COMMAND ""
  BUILD_IN_SOURCE 1
  BUILD_COMMAND make
  INSTALL_COMMAND ""
)

