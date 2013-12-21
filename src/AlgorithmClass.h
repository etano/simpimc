#ifndef AlgorithmClass_H
#define AlgorithmClass_H

#include <iostream>
#include "config.h"
#include "IO/InputFile.h"
#include "IO/IO.h"
#include "PathClass.h" // paths class
#include "RNG/RNG.h" // RNG class
#include "Moves/MoveClass.h" // moves class
#include "Moves/BisectClass.h" // moves class
#include "Observables/ObservableClass.h" // observables class
#include "Observables/EnergyClass.h" // observables class
#include "LoopClass.h"
#include "EventClass.h"
#include "WriteClass.h"

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
