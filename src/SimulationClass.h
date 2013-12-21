#ifndef SimulationClass_H
#define SimulationClass_H

#include <iostream>
#include "config.h"
#include "Communication/Communication.h"
#include "IO/InputClass.h"
#include "IO/IOClass.h"
#include "RNG/RNGClass.h" // RNG class
#include "PathClass.h" // paths class
#include "AlgorithmClass.h" // paths class

class Simulation
{
public:
  // Constructor
  Simulation() {}

  // IO
  Input in;
  IOClass out;
  void SetupIO(string inFile);

  // MPI setup
  CommunicatorClass WorldComm; // This is the global MPIWORLD communicator.
  CommunicatorClass InterComm; // This is for communication between the rank 0 procs of each walker group.
  CommunicatorClass IntraComm; // This is for commmunication between procs within a walker group.
  void BuildMPIModel();
  int procsPerGroup;

  // Run
  void Run();

  // Algorithm
  Algorithm algorithm;
};

#endif
