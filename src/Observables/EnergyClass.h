#ifndef EnergyClass_H
#define EnergyClass_H

#include "ObservableClass.h"

class Energy : public Observable
{
private: 
  vec Etot, KEtot, VEtot, NEtot;
  double E, KE, VE, NE;

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
    Etot.zeros(path.nType);
    KEtot.zeros(path.nType);
    VEtot.zeros(path.nType);
    NEtot.zeros(path.nType);
    dr.set_size(path.nD);
  }

  virtual void Accumulate( const int iType );
  virtual void Output();
  virtual void Print();
  virtual void Stats();
};

#endif
