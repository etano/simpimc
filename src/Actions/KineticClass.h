#ifndef KineticClass_H
#define KineticClass_H

#include "ActionClass.h"

class Kinetic : public Action
{
private:
  int nImages;
protected:

public:
  // Constructor
  Kinetic(Path &path, Input &in, IOClass &out)
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
