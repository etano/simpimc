#ifndef SIMPIMC_SIMULATION_CLASS_H_
#define SIMPIMC_SIMULATION_CLASS_H_

#include "algorithm_class.h"
#include "scaffold/communication/communication.h"

using namespace scaffold::parallel;

/// Class that handles parallelization.
/// It contains the instantiation of the algorithm
/// and random number generator
class Simulation
{
private:
  Input in; ///< Input file object
  IO out; ///< Output file object

  Communicator world_comm; ///< This is the global MPIWORLD communicator.
  Communicator inter_comm; ///< This is for communication between the rank 0 procs of each walker group.
  Communicator intra_comm; ///< This is for commmunication between procs within a walker group.
public:
  /// Constructor does nothing.
  Simulation() {}

  /// Sets up the simulation.
  void SetupSimulation(std::string in_file);

  /// Runs the simulation.
  void Run();
};

#endif // SIMPIMC_SIMULATION_CLASS_H_
