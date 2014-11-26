#ifndef EventClass_H
#define EventClass_H

#include "config.h"
#include <sys/time.h>

class Event
{
private:

protected:

public:
  // Constructor
  Event() { timeSpent = 0.; }

  string name;
  double timeSpent;
  virtual void DoEvent() {};
  virtual void Write() {};
};

#endif
