#ifndef SIMPIMC_LOOP_CLASS_H_
#define SIMPIMC_LOOP_CLASS_H_

#include "event_class.h"

/// Class that builds loops out of events
class Loop : public Event
{
private:
  /// Find an event
  std::shared_ptr<Event> FindEvent(std::string name, std::vector<std::shared_ptr<Event>> &event_list)
  {
    std::vector<std::shared_ptr<Event>>::iterator iter;
    for (iter=event_list.begin(); iter!=event_list.end(); ++iter)
      if ((*iter)->name == name)
        return (*iter);
    return nullptr;
  }

public:
  int n_steps; ///< Number of steps in a loop
  std::vector<std::shared_ptr<Event>> events; ///< Vector of all event objects

  /// Main loop constructor
  Loop()
  {
    name = "Loop";
    n_steps = 1;
  }

  /// Constructor
  Loop(Input &in, std::vector<std::shared_ptr<Event>> &event_list)
  {
    name = "Loop";
    n_steps = in.GetAttribute<int>("n_step");
    ReadLoop(in,event_list);
  };

  void ReadLoop(Input &in, std::vector<std::shared_ptr<Event>> &event_list)
  {
    std::vector<Input> in_list = in.GetChildList();
    for (auto& in: in_list) {
      std::string type = in.GetName();
      if (type == "Move" || type == "Observable") {
        std::string name = in.GetAttribute<std::string>("name");
        std::shared_ptr<Event> event(FindEvent(name,event_list));
        if (event != nullptr)
          events.push_back(event);
        else
          std::cerr << "WARNING: Event not found, " << name << std::endl;
      } else if (type == "Write") {
        events.push_back(FindEvent("Write",event_list));
      } else if (type == "Loop") {
        events.push_back(std::make_shared<Loop>(in,event_list));
      } else if (type == "n_step") {
      } else
        std::cerr << "WARNING: Unrecognized event, " << type << std::endl;
    }
  }

  /// Execute an event
  void DoEvent()
  {
    for (int i=0; i<n_steps; i++)
      for (auto& event: events)
        event->DoEvent();
  }
};

#endif // SIMPIMC_LOOP_CLASS_H_
