#ifndef ContactProbabilityClass_H
#define ContactProbabilityClass_H

#include "ObservableClass.h"
#include "../Actions/ActionClass.h"

class ContactProbability : public Observable
{
private:
  std::vector< std::shared_ptr<Action> > actionList, &fullActionList;
  string speciesA, speciesB;
  int iSpeciesA, iSpeciesB;
  int ZA;
  double total;
protected:
public:
  ContactProbability(Path &tmpPath, std::vector< std::shared_ptr<Action> > &tmpActionList, Input &in, IOClass &out)
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
