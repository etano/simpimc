#ifndef SIMPIMC_FRAMEWORK_CLASS_H_
#define SIMPIMC_FRAMEWORK_CLASS_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <memory>
#include <sys/time.h>
#include <atomic>

#include <scaffold/matrix/matrix.h>
#include <scaffold/algorithm/algorithm.h>
#include <scaffold/communication/communication.h>
#include <scaffold/io/io_xml.h>
#include <scaffold/io/io_hdf5.h>
#include <scaffold/rng/rng.h>

using namespace scaffold::matrix;
using namespace scaffold::algorithm;
using namespace scaffold::io;
using namespace scaffold::parallel;
using namespace scaffold::rand;

#include "simulation_class.h"

/// Class that handles parallelization. It contains the instantiation of the algorithm and random number generator
class Framework
{
private:
  Input in; ///< Input file object
  IO out; ///< Output file object

  Communicator world_comm; ///< This is the global MPIWORLD communicator.
  Communicator inter_comm; ///< This is for communication between the rank 0 procs of each walker group.
  Communicator intra_comm; ///< This is for commmunication between procs within a walker group.
public:
  /// Constructor does nothing.
  Framework(std::string in_file)
  {
    // Input
    in.Load(in_file);

    // Build MPI Model
    uint32_t procs_per_group = in.GetChild("Parallel").GetAttribute<uint32_t>("procs_per_group",1);
    uint32_t N = world_comm.NumProcs();
    assert ((N % procs_per_group) == 0);
    uint32_t n_groups = N/procs_per_group;
    uint32_t my_group = world_comm.MyProc()/procs_per_group;
    world_comm.Split(my_group, intra_comm);
    vec<int> ranks(n_groups);
    for (uint32_t group=0; group<n_groups; group++)
      ranks(group) = group*procs_per_group;
    world_comm.Subset(ranks, inter_comm);
    int n_threads = 1;
#if USE_OPENMP
    #pragma omp parallel
    {
      n_threads = omp_get_num_threads();
    }
#endif
    if (world_comm.MyProc() == 0) {
      std::cout << "Running simpimc on " << in_file << std::endl;
      std::cout <<"# Processes: "<< N << ", # Groups: "<< n_groups
           <<", # Processes/Group: "<< procs_per_group
           <<", # Threads/Process: "<< n_threads << std::endl;
    }

    // Output
    std::stringstream tmp_ss;
    tmp_ss << in.GetChild("IO").GetAttribute<std::string>("output_prefix") << "." << my_group << ".h5";
    std::string output = tmp_ss.str();
    out.Load(output);
    out.Create();

    // Write input data
    out.CreateGroup("Input");
    out.Write("Input/file_name",in_file);
    std::string in_string = in.GetString();
    out.Write("Input/contents",in_string);
  }

  /// Runs the parallel framework
  void Run()
  {
    Input rng_in = in.GetChild("RNG");
    int seed = 12345;
    if (rng_in.node.name == "null") {
#if USE_MPI
      seed = (int)time(0)*(world_comm.MyProc()+1);
#else
      seed = (int)time(0);
#endif
    } else {
#if USE_MPI
      seed = in.GetChild("RNG").GetAttribute<int>("seed",(int)time(0)*(world_comm.MyProc()+1));
#else
      seed = in.GetChild("RNG").GetAttribute<int>("seed",(int)time(0));
#endif
    }
    RNG rng(seed);
    out.CreateGroup("RNG");
    out.Write("RNG/seed",seed);

    // Simulation
    Simulation sim(in, out, rng, inter_comm.MyProc());
    sim.Run();
  }
};

#endif // SIMPIMC_FRAMEWORK_CLASS_H_
