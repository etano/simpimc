#ifndef SIMPIMC_WRITES_CLASS_H_
#define SIMPIMC_WRITES_CLASS_H_

#include "event_class.h"

/// Class the controls most writing to file
class Writes : public Event
{
private:
  double initial; ///< Initial time of day
  double start; ///< Time at which block was started
  double end; ///< Time at which block plus writing has ended
  int block_i; ///< Index of current block
protected:
  const uint32_t proc_i; ///< Index of current processor
  IO &out; ///< Reference to output file
  std::vector<std::shared_ptr<Event>> &events; ///< List of all event
public:
  /// Constructor initializes block index and time of day
  Writes(IO &tmp_out, std::vector<std::shared_ptr<Event>> &tmp_events, const uint32_t t_proc_i)
    : Event(), out(tmp_out), events(tmp_events), proc_i(t_proc_i)
  {
    name = "Write";
    block_i = 0;
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    initial = time.tv_sec + (time.tv_usec / 1000000.);
    end = time.tv_sec + (time.tv_usec / 1000000.);
  }

  /// Perform the write for all other events
  void DoEvent();
};

#endif // SIMPIMC_WRITES_CLASS_H_
