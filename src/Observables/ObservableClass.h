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
  // Path
  Path& path;

  // IO
  Input &in;
  IOClass &out;

  bool FirstTime;
  unsigned int nMeasure;
public:
  // Constructor
  Observable(Path &tmpPath, Input &tmpIn, IOClass &tmpOut)
    : Event(), path(tmpPath), in(tmpIn), out(tmpOut)
  {
    Name = in.get<string>("Name");
    out.CreateGroup(Name);
    FirstTime = 1;
    nMeasure = 0;
  }

  // Functions
  inline void DoEvent() { Accumulate(); };
  virtual void Init() {};
  virtual void Accumulate() {};
  virtual void Reset() {};
  virtual void Write() {};
};

#endif
