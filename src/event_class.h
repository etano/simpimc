#ifndef SIMPIMC_EVENT_CLASS_H_
#define SIMPIMC_EVENT_CLASS_H_

#include "config.h"
#include <sys/time.h>

class Event
{
private:

protected:

public:
  // Constructor
  Event() { time_spent = 0.; }

  std::string name;
  double time_spent;
  virtual void DoEvent() {};
  virtual void Write() {};
};

#endif // SIMPIMC_EVENT_CLASS_H_
