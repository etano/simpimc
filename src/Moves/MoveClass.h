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
    long totalTime = (time.tv_sec * 1000) + (time.tv_usec / 1000);
    MakeMove();
    gettimeofday(&time, NULL); //END-TIME
    totalTime = (((time.tv_sec * 1000) + (time.tv_usec / 1000)) - totalTime);
    timeSpent += totalTime;
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
