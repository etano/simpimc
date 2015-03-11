#include "loop_class.h"

void Loop::Init(Input &in, std::vector<std::shared_ptr<Event>> &event_list)
{
  name = "Loop";
  n_steps = in.GetAttribute<int>("n_step");
  ReadLoop(in, event_list);
}

void Loop::ReadLoop(Input &in, std::vector<std::shared_ptr<Event>> &event_list)
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
      std::shared_ptr<Loop> new_loop(new Loop());
      new_loop->Init(in,event_list);
      events.push_back(new_loop);
    } else if (type == "n_step") {
    } else
      std::cerr << "WARNING: Unrecognized event, " << type << std::endl;
  }
}

std::shared_ptr<Event> Loop::FindEvent(std::string name, std::vector<std::shared_ptr<Event>> &event_list)
{
  std::vector<std::shared_ptr<Event>>::iterator iter;
  for (iter=event_list.begin(); iter!=event_list.end(); ++iter)
    if ((*iter)->name == name)
      return (*iter);
  return nullptr;
}

void Loop::DoEvent()
{
  for (int i=0; i<n_steps; i++) {
    for (auto& event: events) {
      //std::cerr << event->name << std::endl;
      event->DoEvent();
    }
  }
}
