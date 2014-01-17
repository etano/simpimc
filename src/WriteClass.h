#ifndef WritesClass_H
#define WritesClass_H

#include "EventClass.h"
#include "Utils/IO/InputClass.h"
#include "Utils/IO/IOClass.h"

class Writes : public Event
{
private:

protected:
  // Output
  IOClass &Out;

  // Event list
  vector<Event*> &events;
public:
  // Constructor
  Writes(IOClass &tmpOut, vector<Event*> &tmpEvents)
    : Event(), Out(tmpOut), events(tmpEvents)
  {
    name = "Write";
  }

  // Functions
  void DoEvent();
  virtual void Write() {};
};

#endif
