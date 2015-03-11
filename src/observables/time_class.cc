#include "time_class.h"

void Time::Init(Input &in)
{
  struct timeval time;
  gettimeofday(&time, NULL); //END-TIME
  end = time.tv_sec + (time.tv_usec / 1000000.);
  Reset();

}

void Time::Reset()
{
  for (uint i=0; i<event_list.size(); ++i)
    event_list[i]->time_spent = 0;

  struct timeval time;
  gettimeofday(&time, NULL); // Start Time
  start = time.tv_sec + (time.tv_usec / 1000000.);
  time_spent += start - end;
}

void Time::Accumulate()
{
}

void Time::Write()
{
  struct timeval time;
  gettimeofday(&time, NULL); //END-TIME
  end = time.tv_sec + (time.tv_usec / 1000000.);
  double total_time = end - start;
  double norm = 1.;
  if (total_time != 0)
    norm /= total_time;
  std::vector<double> event_times(event_list.size());
  for (uint i=0; i<event_times.size(); ++i)
    event_times[i] = event_list[i]->time_spent;

  if (first_time) {
    first_time = 0;
    std::string data_type = "scalar";
    out.CreateGroup(prefix+"Block");
    out.Write(prefix+"Block/data_type",data_type);
    out.CreateExtendableDataSet(prefix+"Block/", "x", total_time);
    for (uint i=0; i<event_times.size(); ++i) {
      out.CreateGroup(prefix+event_list[i]->name);
      out.Write(prefix+event_list[i]->name+"/data_type",data_type);
      out.CreateExtendableDataSet("/"+prefix+event_list[i]->name+"/", "x", event_times[i]);
    }
  } else {
    out.AppendDataSet(prefix+"Block/", "x", total_time);
    for (uint i=0; i<event_times.size(); ++i)
      out.AppendDataSet("/"+prefix+event_list[i]->name+"/", "x", event_times[i]);
  }

  Reset();
}
