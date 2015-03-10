#ifndef ContactDensityClass_H
#define ContactDensityClass_H

#include "ObservableClass.h"
#include "../Actions/ActionClass.h"

class ContactDensity : public Observable
{
private:
  std::vector< std::shared_ptr<Action> > actionList, &fullActionList;
  string speciesA, speciesB;
  uint iSpeciesA, iSpeciesB;
  uint ZA;
  double total;
protected:
public:
  ContactDensity(Path &tmpPath, std::vector<std::shared_ptr<Action> >& tmpActionList, Input &in, IOClass &out)
    : fullActionList(tmpActionList), Observable(tmpPath, in, out)
  {
    Init(in);
    string data_type = "scalar";
    out.Write(prefix+"/data_type",data_type);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif
