#include "LoopClass.h"

void Loop::Init(Input &in, std::vector< std::shared_ptr<Event> > &eventList)
{
  name = "Loop";
  nSteps = in.getAttribute<int>("nStep");
  ReadLoop(in, eventList);
}

void Loop::ReadLoop(Input &in, std::vector< std::shared_ptr<Event> > &eventList)
{
  vector<Input> inList = in.getChildList();
  for (auto& in: inList) {
    string type = in.getName();
    if (type == "Move" || type == "Observable") {
      string name = in.getAttribute<string>("name");
      std::shared_ptr<Event> event(FindEvent(name,eventList));
      if (event != nullptr)
        events.push_back(event);
      else
        std::cerr << "WARNING: Event not found, " << name << endl;
    } else if (type == "Write") {
      events.push_back(FindEvent("Write",eventList));
    } else if (type == "Loop") {
      std::shared_ptr<Loop> newLoop(new Loop());
      newLoop->Init(in,eventList);
      events.push_back(newLoop);
    } else if (type == "nStep") {
    } else
      std::cerr << "WARNING: Unrecognized event, " << type << endl;
  }
}

std::shared_ptr<Event> Loop::FindEvent(string name, std::vector< std::shared_ptr<Event> > &eventList)
{
  std::vector< std::shared_ptr<Event> >::iterator iter;
  for (iter=eventList.begin(); iter!=eventList.end(); ++iter)
    if ((*iter)->name == name)
      return (*iter);
  return nullptr;
}

void Loop::DoEvent()
{
  for (int i=0; i<nSteps; i++) {
    for (auto& event: events) {
      //cout << event->name << endl;
      event->DoEvent();
    }
  }
}
