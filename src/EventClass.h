#ifndef EventClass_H
#define EventClass_H

#include "config.h"       // Standard libraries

class Event
{
private:

protected:

public:
  // Constructor
  Event() {}

  string name;
  long timeSpent;
  virtual void DoEvent() {};
  virtual void Write() {};
};

#endif
