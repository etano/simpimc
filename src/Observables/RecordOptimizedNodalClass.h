#ifndef RecordOptimizedNodalClass_H
#define RecordOptimizedNodalClass_H

#include "ObservableClass.h"
#include "../Actions/ActionClass.h"

class RecordOptimizedNodal : public Observable
{
private:
  vec<double> paramSetCount;
  std::shared_ptr<Action> action;
  std::vector< std::shared_ptr<Action> > &actionList;
protected:
public:
  RecordOptimizedNodal(Path &tmpPath, std::vector< std::shared_ptr<Action> > &tmpActionList, Input &in, IOClass &out)
    : actionList(tmpActionList), Observable(tmpPath, in, out)
  {
    Init(in);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif
