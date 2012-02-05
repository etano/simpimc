#ifndef R2Class_H
#define R2Class_H

#include "ObservableClass.h"

class R2 : public Observable
{
private:
  mat R2block;
  vec R2tot;
protected:

public:
  R2( Path& pathIn , std::string outputSuffixIn , std::string observableLabelIn , unsigned int skipIn , unsigned int blockIn )
    : Observable( pathIn , outputSuffixIn , observableLabelIn , skipIn , blockIn )
  {
    R2block.zeros(path.nType,path.nPart);
    R2tot.zeros(path.nPart);
  }

  virtual void Accumulate( const int pType );
  virtual void Output();
  virtual void Print();
  virtual void Stats();
};

#endif
