#ifndef PathDumpClass_H
#define PathDumpClass_H

#include "ObservableClass.h"
#include "../Actions/ActionClass.h"

class PathDump : public Observable
{
private:
  int skip, nWriteCalls, nDump;
protected:
public:
  PathDump(Path &tmpPath, Input &in, IOClass &out)
    : Observable(tmpPath, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif
