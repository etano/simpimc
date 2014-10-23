#ifndef AlgorithmClass_H
#define AlgorithmClass_H

#include <iostream>
#include "PathClass.h"
#include "LoopClass.h"
#include "EventClass.h"
#include "WriteClass.h"
#include "Utils/config.h"
#include "Utils/Communication/Communication.h"
#include "Utils/IO/InputClass.h"
#include "Utils/IO/IOClass.h"
#include "Utils/RNG/RNGClass.h"
#include "Actions/ActionClass.h"
#include "Actions/CoulombClass.h"
#include "Actions/KineticClass.h"
#include "Actions/NodalClass.h"
#include "Actions/DavidPairActionClass.h"
#include "Actions/IlkkaPairActionClass.h"
#include "Actions/TrapClass.h"
#include "Moves/MoveClass.h"
#include "Moves/BisectClass.h"
#include "Moves/PermBisectClass.h"
#include "Moves/PermBisectIterativeClass.h"
#include "Observables/ObservableClass.h"
#include "Observables/EnergyClass.h"
#include "Observables/PairCorrelationClass.h"
#include "Observables/PathDumpClass.h"
#include "Observables/StructureFactorClass.h"
#include "Observables/TimeClass.h"

class Algorithm
{
public:
  // Constructor
  Algorithm(CommunicatorClass& WorldComm, CommunicatorClass& InterComm, CommunicatorClass& IntraComm)
   : path(WorldComm,InterComm,IntraComm)
  {}
  void Init(Input &in, IOClass &out, RNG &rng);
  void Run();

  // Algorithm Events
  std::vector<Event*> events;
  Loop mainLoop;

  // Actions
  std::vector<Action*> actions;

  // Datastructure
  Path path;
};

#endif
