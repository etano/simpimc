#ifndef MoveClass_H
#define MoveClass_H

#include "../EventClass.h"
#include "../PathClass.h"
#include "../RNG/RNGClass.h"
#include "../IO/InputClass.h"
#include "../IO/IOClass.h"
#include "../Actions/ActionClass.h"

class Move : public Event
{
private:

protected:
  // Path
  Path& path;

  // RNG
  RNG& rng;

  vector<Action*> &actionList;
public:
  // Constructor
  Move(Path &tmpPath, RNG &tmpRNG, vector<Action*> &tmpActionList,  Input &in, IOClass &out)
    : Event(), path(tmpPath), rng(tmpRNG), actionList(tmpActionList)
  {
    Init(in);
  }
  void Init(Input &in);

  // Moves
  unsigned int moveSkip;
  virtual void MakeMove() {};
  inline void DoEvent() {MakeMove();}

  // Equilibration
  RealType perAcceptDesired;
  unsigned int nEqSweep;
  unsigned int nEqStep;
  RealType stepSize;
  void Equilibrate();

  // Acceptance
  unsigned int nAttempt;
  unsigned int nAccept;
  RealType perAccept;
  void resetCounters();
  RealType getPerAccept();

  // Write
  virtual void Write() {};

  // Permutation Functions
  void assignParticleLabels();
  void setPerm(const int permType, int* perm, int* iPerm, const int i, const int j, const int k);
};

#endif
