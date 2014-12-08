#ifndef LoopClass_H
#define LoopClass_H

#include "EventClass.h"

class Loop : public Event
{
private:

protected:

public:
  // Constructor
  Loop() {};

  void Init(Input &In, std::vector< std::shared_ptr<Event> > &eventList);
  void ReadLoop(Input &In, std::vector< std::shared_ptr<Event> > &eventList);
  std::shared_ptr<Event> FindEvent(string name, std::vector< std::shared_ptr<Event> > &eventList);
  void DoEvent();

  int nSteps;
  vector< std::shared_ptr<Event> > events;
};

#endif
