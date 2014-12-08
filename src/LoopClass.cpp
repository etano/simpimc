#include "LoopClass.h"

void Loop::Init(Input &in, vector<Event*> &eventList)
{
  name = "Loop";
  nSteps = in.getAttribute<int>("nStep");
  ReadLoop(in, eventList);
}

void Loop::ReadLoop(Input &in, vector<Event*> &eventList)
{
  vector<Input> inList = in.getChildList();
  for (auto& in: inList) {
    string type = in.getName();
    if (type == "Move" || type == "Observable") {
      string name = in.getAttribute<string>("name");
      Event *event = FindEvent(name,eventList);
      if (event != NULL)
        events.push_back(event);
      else
        std::cerr << "WARNING: Event not found, " << name << endl;
    } else if (type == "Write") {
      Event *event = FindEvent("Write",eventList);
      events.push_back(event);
    } else if (type == "Loop") {
      Loop *newLoop = new Loop();
      newLoop->Init(in,eventList);
      events.push_back(newLoop);
    } else if (type == "nStep") {
    } else
      std::cerr << "WARNING: Unrecognized event, " << type << endl;
  }
}

Event* Loop::FindEvent(string name, vector<Event*> &eventList)
{
  vector<Event*>::iterator iter;
  for (iter=eventList.begin(); iter!=eventList.end(); ++iter)
    if ((*iter)->name == name)
      return (*iter);
  return NULL;
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
