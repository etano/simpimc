#ifndef AlgorithmClass_H
#define AlgorithmClass_H

#include <iostream>
#include "PathClass.h"
#include "LoopClass.h"
#include "EventClass.h"
#include "WriteClass.h"
#include "config.h"
#include "Actions/ActionClass.h"
#include "Actions/KineticClass.h"
#include "Actions/NodalClass.h"
#include "Actions/BarePairActionClass.h"
#include "Actions/DavidPairActionClass.h"
#include "Actions/IlkkaPairActionClass.h"
#include "Actions/TrapClass.h"
#include "Moves/MoveClass.h"
#include "Moves/BisectClass.h"
#include "Moves/DisplaceParticleClass.h"
#include "Moves/PermBisectClass.h"
#include "Moves/PermBisectIterativeClass.h"
#include "Moves/ShiftRefSliceClass.h"
#include "Observables/ObservableClass.h"
#include "Observables/EnergyClass.h"
#include "Observables/PairCorrelationClass.h"
#include "Observables/PathDumpClass.h"
#include "Observables/PermutationClass.h"
#include "Observables/SignClass.h"
#include "Observables/StructureFactorClass.h"
#include "Observables/TimeClass.h"

class Algorithm
{
public:
  // Constructor
  Algorithm(CommunicatorClass& WorldComm, CommunicatorClass& tmpInterComm, CommunicatorClass& IntraComm)
   : path(WorldComm,tmpInterComm,IntraComm), InterComm(tmpInterComm)
  {}
  void Init(Input &in, IOClass &out, RNG &rng);
  void Run();

  // Communicator
  CommunicatorClass& InterComm;

  // Algorithm Events
  std::vector<Event*> events;
  Loop mainLoop;

  // Actions
  std::vector<Action*> actions;

  // Datastructure
  Path path;
};

#endif
