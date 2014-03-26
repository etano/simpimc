#ifndef AlgorithmClass_H
#define AlgorithmClass_H

#include <iostream>
#include "PathClass.h"
#include "LoopClass.h"
#include "EventClass.h"
#include "WriteClass.h"
#include "Utils/config.h"
#include "Utils/IO/InputClass.h"
#include "Utils/IO/IOClass.h"
#include "Utils/RNG/RNGClass.h"
#include "Actions/ActionClass.h"
#include "Actions/CoulombClass.h"
#include "Actions/KineticClass.h"
#include "Actions/PairActionClass.h"
#include "Actions/TrapClass.h"
#include "Moves/MoveClass.h"
#include "Moves/BisectClass.h"
#include "Observables/ObservableClass.h"
#include "Observables/EnergyClass.h"
#include "Observables/TimeClass.h"

class Algorithm
{
public:
  // Constructor
  Algorithm() {}
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
