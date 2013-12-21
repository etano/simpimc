#ifndef WriteClass_H
#define WriteClass_H

#include "EventClass.h"
#include "IO/InputClass.h"
#include "IO/IOClass.h"

class Write : public Event
{
private:

protected:
  // Output
  IOClass &Out;

  // Event list
  vector<Event*> &events;
public:
  // Constructor
  Write(IOClass &tmpOut, vector<Event*> &tmpEvents)
    : Event(), Out(tmpOut), events(tmpEvents)
  {
    Name = "Write";
  }

  // Functions
  void DoEvent();
};

#endif
