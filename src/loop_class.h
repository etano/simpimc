#ifndef SIMPIMC_LOOP_CLASS_H_
#define SIMPIMC_LOOP_CLASS_H_

#include "event_class.h"

class Loop : public Event
{
private:

protected:

public:
  // Constructor
  Loop() {};

  void Init(Input &in, std::vector<std::shared_ptr<Event>> &event_list);
  void ReadLoop(Input &in, std::vector<std::shared_ptr<Event>> &event_list);
  std::shared_ptr<Event> FindEvent(std::string name, std::vector<std::shared_ptr<Event>> &event_list);
  void DoEvent();

  int n_steps;
  std::vector<std::shared_ptr<Event>> events;
};

#endif // SIMPIMC_LOOP_CLASS_H_
