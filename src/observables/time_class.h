#ifndef SIMPIMC_OBSERVABLES_TIME_CLASS_H_
#define SIMPIMC_OBSERVABLES_TIME_CLASS_H_

#include "observable_class.h"

/// Measure the amount of time performing each event
class Time : public Observable
{
private:
  double start; ///< Block start time
  double end; ///< Block end time
  std::vector<std::shared_ptr<Event>> &event_list; ///< Reference to vector of all events

  /// Initiate the observable
  virtual void Init(Input &in);

  /// Accumulate the observable
  virtual void Accumulate();

  /// Reset the observable's counters
  virtual void Reset();
public:
  /// Constructor calls Init
  Time(Path &path, std::vector<std::shared_ptr<Event>> &tmp_event_list, Input &in, IO &out)
    : event_list(tmp_event_list), Observable(path, in, out)
  {
    Init(in);
  }

  /// Write relevant information about an observable to the output
  virtual void Write();
};

#endif // SIMPIMC_OBSERVABLES_TIME_CLASS_H_
