#ifndef SIMPIMC_WRITES_CLASS_H_
#define SIMPIMC_WRITES_CLASS_H_

#include "event_class.h"

class Writes : public Event
{
private:
  int block_i;
  double initial, start, end;
protected:
  // Output
  IO &out;

  // Event list
  std::vector<std::shared_ptr<Event>> &events;

  // Communicator
  Communicator &inter_comm;
public:
  // Constructor
  Writes(IO &tmp_out, std::vector<std::shared_ptr<Event>> &tmp_events, Communicator &tmp_inter_comm)
    : Event(), out(tmp_out), events(tmp_events), inter_comm(tmp_inter_comm)
  {
    name = "Write";
    block_i = 0;
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    initial = time.tv_sec + (time.tv_usec / 1000000.);
    end = time.tv_sec + (time.tv_usec / 1000000.);
  }

  // Functions
  void DoEvent();
  virtual void Write() {};
};

#endif // SIMPIMC_WRITES_CLASS_H_
