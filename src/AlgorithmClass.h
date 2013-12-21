#ifndef AlgorithmClass_H
#define AlgorithmClass_H

#include <iostream>
#include "config.h"
#include "PathClass.h"
#include "LoopClass.h"
#include "EventClass.h"
#include "WriteClass.h"
#include "IO/InputClass.h"
#include "IO/IOClass.h"
#include "RNG/RNGClass.h"
#include "Actions/ActionClass.h"
#include "Actions/KineticClass.h"
#include "Moves/MoveClass.h"
#include "Moves/BisectClass.h"
#include "Observables/ObservableClass.h"
#include "Observables/EnergyClass.h"

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

  // Datastructure
  Path path;

};

#endif
