#ifndef ObservableClass_H
#define ObservableClass_H

#include "../EventClass.h"
#include "../PathClass.h"
#include "../IO/InputClass.h"
#include "../IO/IOClass.h"

class Observable : public Event
{
private:

protected:
  Path& path;
  IOClass &out;

public:
  // Constructor
  Observable(Path &tmpPath, Input &in, IOClass &tmpOut)
    : Event(), path(tmpPath), out(tmpOut)
  {
    name = in.getAttribute<string>("name");
    out.CreateGroup("Observables/"+name);
    firstTime = 1;
    nMeasure = 0;
  }

  bool firstTime;
  int nMeasure, timeSpent;

  // Functions
  inline void DoEvent() { Accumulate(); };
  virtual void Init(Input &in) {};
  virtual void Accumulate() {};
  virtual void Reset() {};
  virtual void Write() {};
};

#endif
