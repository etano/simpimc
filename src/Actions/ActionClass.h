#ifndef ActionClass_H
#define ActionClass_H

#include "../EventClass.h"
#include "../PathClass.h"
#include "../Utils/IO/InputClass.h"
#include "../Utils/IO/IOClass.h"

class Action : public Event
{
private:

protected:
  // Path
  Path& path;

  // IO
  IOClass &out;

  bool FirstTime;

  // Find species offset
  void GetOffset(string species, int &iSpecies, int &offset);

public:
  // Constructor
  Action(Path &tmpPath, Input &in, IOClass &tmpOut)
    : path(tmpPath), out(tmpOut)
  {
    name = in.getAttribute<string>("name");
    type = in.getAttribute<string>("type");
    out.CreateGroup("Actions/"+name);
    out.Write("Actions/"+name+"/type",type);
  }

  string type;

  // Functions
  virtual void Init(Input &in) {};
  virtual RealType DActionDBeta() {};
  virtual RealType GetAction(int b0, int b1, vector<int> &particles, int level) {};
  virtual void Write() {};
  void DoEvent() {};
};

#endif
