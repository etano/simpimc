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
  for (uint i=0; i<eventList.size(); ++i)
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
  double totalTime = end - start;
  double norm = 1.;
  if (totalTime != 0)
    norm /= totalTime;
  vector<double> eventTimes(eventList.size());
  for (uint i=0; i<eventTimes.size(); ++i)
    eventTimes[i] = eventList[i]->timeSpent;

  if (firstTime) {
    firstTime = 0;
    string data_type = "scalar";
    out.CreateGroup(prefix+"Block");
    out.Write(prefix+"Block/data_type",data_type);
    out.CreateExtendableDataSet(prefix+"Block/", "x", totalTime);
    for (uint i=0; i<eventTimes.size(); ++i) {
      out.CreateGroup(prefix+eventList[i]->name);
      out.Write(prefix+eventList[i]->name+"/data_type",data_type);
      out.CreateExtendableDataSet("/"+prefix+eventList[i]->name+"/", "x", eventTimes[i]);
    }
  } else {
    out.AppendDataSet(prefix+"Block/", "x", totalTime);
    for (uint i=0; i<eventTimes.size(); ++i)
      out.AppendDataSet("/"+prefix+eventList[i]->name+"/", "x", eventTimes[i]);
  }

  Reset();
}
