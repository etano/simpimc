#ifndef EnergyClass_H
#define EnergyClass_H

#include "ObservableClass.h"
#include "../Actions/ActionClass.h"

class Energy : public Observable
{
private:
  vector<Action*> &actionList;
  vector<RealType> Es;
protected:
public:
  Energy(Path &tmpPath, vector<Action*> &tmpActionList, Input &in, IOClass &out)
    : actionList(tmpActionList), Observable(tmpPath, in, out)
  {
    Init();
  }

  virtual void Init();
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif
