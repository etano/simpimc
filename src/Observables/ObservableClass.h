#ifndef ObservableClass_H
#define ObservableClass_H

#include "../EventClass.h"
#include "../PathClass.h"
#include "../Utils/IO/InputClass.h"
#include "../Utils/IO/IOClass.h"

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
  int nMeasure;

  // Functions
  inline void DoEvent() {
    struct timeval time;
    //gettimeofday(&time, NULL); // Start Time
    RealType start = time.tv_sec + (time.tv_usec / 1000000.);
    Accumulate();
    //gettimeofday(&time, NULL); //END-TIME
    RealType end = time.tv_sec + (time.tv_usec / 1000000.);
    timeSpent += end - start;
  }
  virtual void Init(Input &in) {};
  virtual void Accumulate() {};
  virtual void Reset() {};
  virtual void Write() {};
};

#endif
