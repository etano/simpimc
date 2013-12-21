#ifndef MoveClass_H
#define MoveClass_H

#include "../EventClass.h"
#include "../PathClass.h"
#include "../RNG/RNGClass.h"
#include "../IO/InputClass.h"
#include "../IO/IOClass.h"

class Move : public Event
{
private:

protected:
  // Path
  Path& path;

  // RNG
  RNG& rng;
public:
  // Constructor
  Move(Path &tmpPath, RNG &tmpRNG, Input &in, IOClass &out)
    : Event(), path(tmpPath), rng(tmpRNG)
  {
    Init(in);
  }
  void Init(Input &in);

  // Moves
  unsigned int moveSkip;
  virtual void MakeMove() {};
  inline void DoEvent() {MakeMove();}

  // Equilibration
  double perAcceptDesired;
  unsigned int nEqSweep;
  unsigned int nEqStep;
  double stepSize;
  void Equilibrate();

  // Acceptance
  unsigned int nAttempt;
  unsigned int nAccept;
  double perAccept;
  void resetCounters();
  double getPerAccept();

  // Write
  virtual void Write() {};

  // Permutation Functions
  void assignParticleLabels();
  void setPerm(const int permType, int* perm, int* iPerm, const int i, const int j, const int k);
};

#endif
