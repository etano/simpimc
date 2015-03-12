#ifndef SIMPIMC_SIMULATION_CLASS_H_
#define SIMPIMC_SIMULATION_CLASS_H_

#include "algorithm_class.h"

class Simulation
{
public:
  // Constructor
  Simulation() {}

  // IO
  Input in;
  IO out;

  // MPI setup
  Communicator world_comm; // This is the global MPIWORLD communicator.
  Communicator inter_comm; // This is for communication between the rank 0 procs of each walker group.
  Communicator intra_comm; // This is for commmunication between procs within a walker group.
  uint procs_per_group;

  void SetupSimulation(std::string in_file);
  void Run();

};

#endif // SIMPIMC_SIMULATION_CLASS_H_
