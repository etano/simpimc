#ifndef ActionClass_H
#define ActionClass_H

#include "../PathClass.h"
#include "../IO/InputFile.h"
#include "../IO/IO.h"

class Action
{
private:

protected:
  // Path
  Path& path;

  // IO
  Input &in;
  IOClass &out;

  bool FirstTime;
public:
  // Constructor
  Action(Path &tmpPath, Input &tmpIn, IOClass &tmpOut)
    : path(tmpPath), in(tmpIn), out(tmpOut)
  {
    Name = in.get<string>("Name");
    out.CreateGroup(Name);
  }

  string Name;

  // Functions
  virtual void Init() {};
  virtual RealType DActionDBeta() {};
  virtual RealType Action() {};
  virtual void Write() {};
};

#endif
