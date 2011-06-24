#ifndef SimulationClass_H
#define SimulationClass_H

#include "StandardLibraries.h" // Standard libraries
#include "GlobalConstants.h"
#include "PathClass.h" // paths class
#include "RNGClass.h" // RNG class

// Moves
#include "Moves/MoveClass.h" // moves class

// Observables
#include "Observables/ObservableClass.h" // observables class

class Simulation
{
private:   

protected:

public:
  // Constructor
  Simulation( const int nPartIn , const int nDIn , const int nBeadIn, const double betaIn , const int fermiIn , const int halfspaceIn , const int nodeTypeIn , const int useNodeDistIn , const double LIn )
    : rng((int)time(0)) , path( nPartIn , nDIn , nBeadIn, betaIn , fermiIn , halfspaceIn , nodeTypeIn , useNodeDistIn , LIn )
  {
  }  

  // Simulation Random Number Generator
  RNG rng;  
  
  // Simulation Paths
  Path path;

  // Simulation Moves
  std::vector<Move*> moves;

  // Simulation Observables
  std::vector<Observable*> observables;
};

#endif
