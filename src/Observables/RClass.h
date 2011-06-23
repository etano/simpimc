#ifndef RClass_H
#define RClass_H

#include "ObservableClass.h"

class R : public Observable
{
private: 
  vec Rtemp;
  mat Rtot;
protected:

public:
  R( Path& pathIn , std::string outputSuffixIn , std::string observableLabelIn , unsigned int skipIn , unsigned int blockIn )
    : Observable( pathIn , outputSuffixIn , observableLabelIn , skipIn , blockIn )
  { 
    Rtemp.zeros(path.nPart);
    Rtot.zeros(path.nType,path.nPart);
  }

  virtual void Accumulate( const int pType );
  virtual void Output();
  virtual void Print();
  virtual void Stats();
};

#endif
