#ifndef EnergyClass_H
#define EnergyClass_H

#include "ObservableClass.h"
#include "../Actions/ActionClass.h"

class Energy : public Observable
{
private:
  std::vector< std::shared_ptr<Action> > &actionList;
  vec<double> Es, Vs;
  bool measureV;
protected:
public:
  Energy(Path &tmpPath, std::vector< std::shared_ptr<Action> > &tmpActionList, Input &in, IOClass &out)
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
