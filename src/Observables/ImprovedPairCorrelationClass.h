#ifndef ImprovedPairCorrelationClass_H
#define ImprovedPairCorrelationClass_H

#include "ObservableClass.h"
#include "../Actions/ActionClass.h"

class ImprovedPairCorrelation : public Observable
{
private:
  std::vector< std::shared_ptr<Action> > &actionList;
  string speciesA, speciesB;
  int iSpeciesA, iSpeciesB;
  Histogram gr;
protected:
public:
  ImprovedPairCorrelation(Path &tmpPath, std::vector< std::shared_ptr<Action> > &tmpActionList, Input &in, IOClass &out)
    : actionList(tmpActionList), Observable(tmpPath, in, out)
  {
    Init(in);
    string data_type = "histogram";
    out.Write(prefix+"/data_type",data_type);
  }

  virtual void Init(Input &in);
  virtual void Reset();
  virtual void Accumulate();
  virtual void Write();
};

#endif
