#ifndef TrapClass_H
#define TrapClass_H

#include "ActionClass.h"

class Trap : public Action
{
private:
  RealType omega;
protected:

public:
  // Constructor
  Trap(Path &path, Input &in, IOClass &out)
    : Action(path,in,out)
  {
    Init();
  }

  // Functions
  virtual void Init();
  virtual RealType DActionDBeta();
  virtual RealType GetAction(int b0, int b1, vector<int> &particles, int level);
  virtual void Write();
};

#endif
