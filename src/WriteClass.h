#ifndef WritesClass_H
#define WritesClass_H

#include "EventClass.h"

class Writes : public Event
{
private:
  int iBlock;
  double initial, start, end;
protected:
  // Output
  IOClass &Out;

  // Event list
  vector<Event*> &events;

  // Communicator
  CommunicatorClass &InterComm;
public:
  // Constructor
  Writes(IOClass &tmpOut, vector<Event*> &tmpEvents, CommunicatorClass &tmpInterComm)
    : Event(), Out(tmpOut), events(tmpEvents), InterComm(tmpInterComm)
  {
    name = "Write";
    iBlock = 0;
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    initial = time.tv_sec + (time.tv_usec / 1000000.);
    end = time.tv_sec + (time.tv_usec / 1000000.);
  }

  // Functions
  void DoEvent();
  virtual void Write() {};
};

#endif
