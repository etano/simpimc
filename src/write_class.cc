#include "write_class.h"

void Writes::DoEvent()
{
  struct timeval time;
  gettimeofday(&time, NULL); // Start Time
  start = time.tv_sec + (time.tv_usec / 1000000.);

  std::vector<std::shared_ptr<Event>>::iterator iter;
  for (iter=events.begin(); iter!=events.end(); ++iter)
    (*iter)->Write();

  gettimeofday(&time, NULL); // Current Time
  double block_time = start - end;
  double total_time = start - initial;
  std::cout << "Clone #: " << inter_comm.MyProc() << ", Block #: " << block_i << ", Block Time: " << block_time << ", Total Time: " << total_time << std::endl;

  block_i += 1;

  gettimeofday(&time, NULL); // End Time
  end = time.tv_sec + (time.tv_usec / 1000000.);
  time_spent += end - start;
}
