#ifndef MoveClass_H
#define MoveClass_H

#include "../EventClass.h"
#include "../PathClass.h"
#include "../Utils/RNG/RNGClass.h"
#include "../Utils/IO/InputClass.h"
#include "../Utils/IO/IOClass.h"
#include "../Actions/ActionClass.h"

class Move : public Event
{
private:

protected:
  Path& path;
  RNG& rng;
  IOClass& out;

  vector<Action*> &actionList;
public:
  // Constructor
  Move(Path &tmpPath, RNG &tmpRNG, vector<Action*> &tmpActionList,  Input &in, IOClass &tmpOut)
    : Event(), path(tmpPath), rng(tmpRNG), actionList(tmpActionList), out(tmpOut)
  {
    name = in.getAttribute<string>("name");
    out.CreateGroup("Moves/"+name);
    firstTime = 1;
    Reset();
  }

  // Moves
  inline void DoEvent() {
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    RealType start = time.tv_sec + (time.tv_usec / 1000000.);
    MakeMove();
    gettimeofday(&time, NULL); //END-TIME
    RealType end = time.tv_sec + (time.tv_usec / 1000000.);
    timeSpent += end - start;
  }
  virtual void Init(Input &in) {};
  virtual void MakeMove() {};

  // Acceptance
  bool firstTime;
  unsigned int nAttempt, nAccept;
  RealType getPerAccept();
  void Reset();

  // Write
  void Write();
};

#endif
