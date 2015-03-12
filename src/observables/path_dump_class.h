#ifndef SIMPIMC_OBSERVABLES_PATH_DUMP_CLASS_H_
#define SIMPIMC_OBSERVABLES_PATH_DUMP_CLASS_H_

#include "observable_class.h"

class PathDump : public Observable
{
private:
  uint n_write_calls, n_dump;
protected:
public:
  PathDump(Path &path, Input &in, IO &out)
    : Observable(path, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_PATH_DUMP_CLASS_H_
