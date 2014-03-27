#ifndef SimulationClass_H
#define SimulationClass_H

#include <iostream>
#include "Utils/config.h"
#include "Utils/Communication/Communication.h"
#include "Utils/IO/InputClass.h"
#include "Utils/IO/IOClass.h"
#include "Utils/RNG/RNGClass.h" // RNG class
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

  // MPI setup
  CommunicatorClass WorldComm; // This is the global MPIWORLD communicator.
  CommunicatorClass InterComm; // This is for communication between the rank 0 procs of each walker group.
  CommunicatorClass IntraComm; // This is for commmunication between procs within a walker group.
  int procsPerGroup;

  void SetupSimulation(string inFile);
  void Run();

  // Algorithm
  Algorithm algorithm;
};

#endif
