#ifndef SIMPIMC_LOOP_CLASS_H_
#define SIMPIMC_LOOP_CLASS_H_

#include "event_class.h"

/// Class that builds loops out of events
class Loop : public Event
{
private:
  /// Reads a loop
  void ReadLoop(Input &in, std::vector<std::shared_ptr<Event>> &event_list);

  /// Find an event
  std::shared_ptr<Event> FindEvent(std::string name, std::vector<std::shared_ptr<Event>> &event_list);
protected:

public:
  /// Constructor doesn't do anything
  Loop() {};

  /// Initialize the loops
  void Init(Input &in, std::vector<std::shared_ptr<Event>> &event_list);

  /// Execute an event
  void DoEvent();

  int n_steps; ///< Number of steps in a loop
  std::vector<std::shared_ptr<Event>> events; ///< Vector of all event objects
};

#endif // SIMPIMC_LOOP_CLASS_H_
