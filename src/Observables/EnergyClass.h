#ifndef EnergyClass_H
#define EnergyClass_H

#include "ObservableClass.h"

class Energy : public Observable
{
private:
  vec E, KE, VE, NE;
  double Etot, KEtot, VEtot, NEtot;
  vec dr;

  // Energy Observables
  double getKE();
  double getVE();
  double getVEint();
  double getVEext();
  double getNE();
protected:
public:
  Energy( Path& pathIn , std::string outputSuffixIn , std::string observableLabelIn , unsigned int skipIn , unsigned int blockIn )
    : Observable( pathIn , outputSuffixIn , observableLabelIn , skipIn , blockIn )
  {
    E.zeros(path.nType);
    KE.zeros(path.nType);
    VE.zeros(path.nType);
    NE.zeros(path.nType);
    dr.zeros(path.nD);
  }

  virtual void Accumulate( const int iType );
  virtual void Output();
  virtual void Print();
  virtual void Stats();
};

#endif
