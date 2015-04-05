#ifndef SIMPIMC_EVENT_CLASS_H_
#define SIMPIMC_EVENT_CLASS_H_

#include "config.h"

/// Parent class for all moves, observables, loops, and writes
class Event
{
public:
  std::string name; ///< Name of the event
  double time_spent; ///< Amount of time spent performing event

  /// Constructor only initializes time_spent to 0
  Event() { time_spent = 0.; }

  /// Execute the event
  virtual void DoEvent() {};

  /// Write out information about the event
  virtual void Write() {};
};

#endif // SIMPIMC_EVENT_CLASS_H_
