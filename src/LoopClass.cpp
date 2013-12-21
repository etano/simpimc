#include "LoopClass.h"

void Loop::Init(Input &In, vector<Event*> &eventList)
{
  Name = "Loop";
  nSteps = In.get<int>("nStep");
  ReadLoop(In,eventList);
}

void Loop::ReadLoop(Input &In, vector<Event*> &eventList)
{
  vector<Input*> inList;
  In.getInList(inList);
  for (int i = 0; i < inList.size(); ++i) {
    string type = inList[i]->getKey();
    if (type == "Move" || type == "Observable") {
      string name = inList[i]->getValue<std::string>();
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
      newLoop->Init(*inList[i],eventList);
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
    if ((*iter)->Name == name)
      return (*iter);
  return NULL;
}

void Loop::DoEvent()
{
  vector<Event*>::iterator iter;
  for (int i=0; i<nSteps; i++) {
    //cout << i << " ";
    for (iter=events.begin(); iter!=events.end(); ++iter) {
      //cout << (*iter)->Name << endl;
      (*iter)->DoEvent();
    }
  }
}
