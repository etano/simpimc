#include "write_class.h"

void Writes::DoEvent()
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
