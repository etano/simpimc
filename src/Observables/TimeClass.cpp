#include "TimeClass.h"

void Time::Init(Input &in)
{
  struct timeval time;
  gettimeofday(&time, NULL); //END-TIME
  end = time.tv_sec + (time.tv_usec / 1000000.);
  Reset();

}

void Time::Reset()
{
  for (int i=0; i<eventList.size(); ++i)
    eventList[i]->timeSpent = 0;

  struct timeval time;
  gettimeofday(&time, NULL); // Start Time
  start = time.tv_sec + (time.tv_usec / 1000000.);
  timeSpent += start - end;
}

void Time::Accumulate()
{
}

void Time::Write()
{
  struct timeval time;
  gettimeofday(&time, NULL); //END-TIME
  end = time.tv_sec + (time.tv_usec / 1000000.);
  RealType totalTime = end - start;
  RealType norm = 1.;
  if (totalTime != 0)
    norm /= totalTime;
  vector<RealType> eventTimes(eventList.size());
  for (int i=0; i<eventTimes.size(); ++i)
    eventTimes[i] = eventList[i]->timeSpent;

  if (firstTime) {
    firstTime = 0;
    for (int i=0; i<eventTimes.size(); ++i)
      out.CreateExtendableDataSet("/Observables/"+name+"/", eventList[i]->name, eventTimes[i]);
    out.CreateExtendableDataSet("/Observables/"+name+"/", "Block", totalTime);
  } else {
    for (int i=0; i<eventTimes.size(); ++i)
      out.AppendDataSet("/Observables/"+name+"/", eventList[i]->name, eventTimes[i]);
    out.AppendDataSet("/Observables/"+name+"/", "Block", totalTime);
  }

  Reset();
}
