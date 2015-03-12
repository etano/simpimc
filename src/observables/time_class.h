#ifndef SIMPIMC_OBSERVABLES_TIME_CLASS_H_
#define SIMPIMC_OBSERVABLES_TIME_CLASS_H_

#include "observable_class.h"

class Time : public Observable
{
private:
  std::vector<std::shared_ptr<Event>> &event_list;
  double start, end;
protected:
public:
  Time(Path &path, std::vector<std::shared_ptr<Event>> &tmp_event_list, Input &in, IO &out)
    : event_list(tmp_event_list), Observable(path, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_TIME_CLASS_H_
