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
  virtual void DActionDBeta();
  virtual void Action(int slice1, int slice2, const Array<int,1> &activeParticles, int level);
  virtual void Write();
};

#endif
