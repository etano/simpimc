#ifndef TimeClass_H
#define TimeClass_H

#include "ObservableClass.h"
#include "../EventClass.h"

class Time : public Observable
{
private:
  std::vector< std::shared_ptr<Event> > &eventList;
  double start, end;
protected:
public:
  Time(Path &tmpPath, std::vector< std::shared_ptr<Event> > &tmpEventList, Input &in, IOClass &out)
    : eventList(tmpEventList), Observable(tmpPath, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif
