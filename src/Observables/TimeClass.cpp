#include "TimeClass.h"

void Time::Init(Input &in)
{
  eventTimes.resize(eventList.size());
  Reset();
}

void Time::Reset()
{
  for (int i=0; i<eventTimes.size(); ++i) {
    eventTimes[i] = 0.;
    eventList[i]->timeSpent = 0;
  }

  struct timeval time;
  gettimeofday(&time, NULL); // Start Time
  timeSpent = (time.tv_sec * 1000) + (time.tv_usec / 1000);
}

void Time::Accumulate()
{
  for (int i=0; i<eventTimes.size(); ++i)
    eventTimes[i] += eventList[i]->timeSpent;
  struct timeval time;
  gettimeofday(&time, NULL); //END-TIME
  timeSpent = (((time.tv_sec * 1000) + (time.tv_usec / 1000)) - timeSpent);
  //cout << timeSpent << endl;
}

void Time::Write()
{
  struct timeval time;
  gettimeofday(&time, NULL); //END-TIME
  timeSpent = (((time.tv_sec * 1000) + (time.tv_usec / 1000)) - timeSpent);
  //cout << eventTimes[0] << " " << timeSpent << endl;
  RealType norm = 1.;
  //if (timeSpent != 0)
  //  norm /= timeSpent;
  for (int i=0; i<eventTimes.size(); ++i)
    eventTimes[i] *= norm;

  if (firstTime) {
    firstTime = 0;
    for (int i=0; i<eventTimes.size(); ++i)
      out.CreateExtendableDataSet("/Observables/"+name+"/", eventList[i]->name, eventTimes[i]);
  } else {
    for (int i=0; i<eventTimes.size(); ++i)
      out.AppendDataSet("/Observables/"+name+"/", eventList[i]->name, eventTimes[i]);
  }

  Reset();
}
