#ifndef ImportanceWeightClass_H
#define ImportanceWeightClass_H

#include "ObservableClass.h"
#include "../Actions/ActionClass.h"

class ImportanceWeight : public Observable
{
private:
  std::vector< std::shared_ptr<Action> > &actionList;
  vec<double> IWs;
protected:
public:
  ImportanceWeight(Path &tmpPath, std::vector< std::shared_ptr<Action> > &tmpActionList, Input &in, IOClass &out)
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
