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
  Move(Path &tmpPath, RNG &tmpRNG, vector<Action*> &tmpActionList, Input &in, IOClass &tmpOut)
    : Event(), path(tmpPath), rng(tmpRNG), actionList(tmpActionList), out(tmpOut)
  {
    name = in.getAttribute<string>("name");
    type = in.getAttribute<string>("type");
    out.CreateGroup("Moves/"+name);
    out.Write("Moves/"+name+"/type",type);
    firstTime = 1;
    Reset();
  }

  string type;
  virtual void DoEvent();

  // Moves
  virtual void Init(Input &in) {};
  virtual int Attempt() {};
  virtual void Accept() {};
  virtual void Reject() {};

  // Acceptance
  bool firstTime;
  unsigned int nAttempt, nAccept;
  virtual void Reset();

  // Write
  virtual void Write();

};

#endif
