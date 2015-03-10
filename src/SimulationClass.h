#ifndef SimulationClass_H
#define SimulationClass_H

#include "config.h"
#include "PathClass.h"
#include "AlgorithmClass.h"

class Simulation
{
public:
  // Constructor
  Simulation() {}

  // IO
  Input in;
  IOClass out;

  // MPI setup
  CommunicatorClass WorldComm; // This is the global MPIWORLD communicator.
  CommunicatorClass InterComm; // This is for communication between the rank 0 procs of each walker group.
  CommunicatorClass IntraComm; // This is for commmunication between procs within a walker group.
  uint procsPerGroup;

  void SetupSimulation(string inFile);
  void Run();

};

#endif
