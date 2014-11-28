#include "WriteClass.h"

void Writes::DoEvent()
{
  struct timeval time;
  gettimeofday(&time, NULL); // Start Time
  start = time.tv_sec + (time.tv_usec / 1000000.);

  //vector<Event*>::iterator iter;
  //for (iter=events.begin(); iter!=events.end(); ++iter)
  //  (*iter)->Write();

  gettimeofday(&time, NULL); // Current Time
  double blockTime = start - end;
  double totalTime = start - initial;
  cout << "Clone #: " << InterComm.MyProc() << ", Block #: " << iBlock << ", Block Time: " << blockTime << ", Total Time: " << totalTime << endl;

  iBlock += 1;

  gettimeofday(&time, NULL); // End Time
  end = time.tv_sec + (time.tv_usec / 1000000.);
  timeSpent += end - start;
}
