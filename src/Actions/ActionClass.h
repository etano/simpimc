#ifndef ActionClass_H
#define ActionClass_H

#include "../EventClass.h"
#include "../PathClass.h"
#include "../IO/InputClass.h"
#include "../IO/IOClass.h"

class Action : public Event
{
private:

protected:
  // Path
  Path& path;

  // IO
  IOClass &out;

  bool FirstTime;
public:
  // Constructor
  Action(Path &tmpPath, Input &in, IOClass &tmpOut)
    : path(tmpPath), out(tmpOut)
  {
    name = in.getAttribute<string>("name");
    out.CreateGroup("Actions/"+name);
  }

  string name;

  // Functions
  virtual void Init(Input &in) {};
  virtual RealType DActionDBeta() {};
  virtual RealType GetAction(int b0, int b1, vector<int> &particles, int level) {};
  virtual void Write() {};
  void DoEvent() {};
};

#endif
