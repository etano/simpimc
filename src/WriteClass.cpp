#include "WriteClass.h"

void Writes::DoEvent()
{
  struct timeval time;
  //gettimeofday(&time, NULL); // Start Time
  RealType start = time.tv_sec + (time.tv_usec / 1000000.);

  vector<Event*>::iterator iter;
  for (iter=events.begin(); iter!=events.end(); ++iter)
    (*iter)->Write();

  //gettimeofday(&time, NULL); //END-TIME
  RealType end = time.tv_sec + (time.tv_usec / 1000000.);
  timeSpent += end - start;
}
