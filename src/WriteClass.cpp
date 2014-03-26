#include "WriteClass.h"

void Writes::DoEvent()
{
  struct timeval time;
  gettimeofday(&time, NULL); // Start Time
  start = time.tv_sec + (time.tv_usec / 1000000.);

  vector<Event*>::iterator iter;
  for (iter=events.begin(); iter!=events.end(); ++iter)
    (*iter)->Write();

  cout << "---" << endl;
  cout << "Block #: " << iBlock << endl;
  gettimeofday(&time, NULL); // Current Time
  RealType blockTime = start - end;
  cout << "Block Time: " << blockTime << endl;
  RealType totalTime = start - initial;
  cout << "Total Time: " << totalTime << endl;
  iBlock += 1;

  gettimeofday(&time, NULL); // End Time
  end = time.tv_sec + (time.tv_usec / 1000000.);
  timeSpent += end - start;
}
