#ifndef SIMPIMC_WRITES_CLASS_H_
#define SIMPIMC_WRITES_CLASS_H_

#include "event_class.h"
#include "../actions/action_class.h"

/// Class the controls most writing to file
class Writes : public Event
{
private:
  double initial; ///< Initial time of day
  double start; ///< Time at which block was started
  double end; ///< Time at which block plus writing has ended
  int block_i; ///< Index of current block
  const uint32_t proc_i; ///< Index of current processor
  IO &out; ///< Reference to output file
  std::vector<std::shared_ptr<Event>> &events; ///< List of all event
  std::vector<std::shared_ptr<Action>> &actions; ///< List of all actions
public:
  /// Constructor initializes block index and time of day
  Writes(IO &tmp_out, std::vector<std::shared_ptr<Event>> &tmp_events, std::vector<std::shared_ptr<Action>> &tmp_actions, const uint32_t t_proc_i)
    : Event(), out(tmp_out), events(tmp_events), actions(tmp_actions), proc_i(t_proc_i)
  {
    name = "Write";
    block_i = 0;
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    initial = time.tv_sec + (time.tv_usec / 1000000.);
    end = time.tv_sec + (time.tv_usec / 1000000.);
  }

  /// Perform the write for all other events
  void DoEvent()
  {
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    start = time.tv_sec + (time.tv_usec / 1000000.);

    // Write events
    for (auto &event: events)
      event->Write();

    // Write actions
    for (auto &action: actions)
      action->Write();

    gettimeofday(&time, NULL); // Current Time
    double block_time = start - end;
    double total_time = start - initial;
    std::cout << "Clone #: " << proc_i << ", Block #: " << block_i << ", Block Time: " << block_time << ", Total Time: " << total_time << std::endl;

    block_i += 1;

    gettimeofday(&time, NULL); // End Time
    end = time.tv_sec + (time.tv_usec / 1000000.);
    time_spent += end - start;
  }
};

#endif // SIMPIMC_WRITES_CLASS_H_
