#ifndef SIMPIMC_OBSERVABLES_TIME_CLASS_H_
#define SIMPIMC_OBSERVABLES_TIME_CLASS_H_

#include "observable_class.h"

/// Measure the amount of time performing each event
class Time : public Observable {
   private:
    double start;                                     ///< Block start time
    double end;                                       ///< Block end time
    std::vector<std::shared_ptr<Event>> &event_list;  ///< Reference to vector of all events

    /// Accumulate the observable
    virtual void Accumulate() {}

    /// Reset the observable's counters
    virtual void Reset() {
        for (uint32_t i = 0; i < event_list.size(); ++i)
            event_list[i]->time_spent = 0;

        struct timeval time;
        gettimeofday(&time, NULL);  // Start Time
        start = time.tv_sec + (time.tv_usec / 1000000.);
        time_spent += start - end;
    }

   public:
    /// Constructor calls Init
    Time(Path &path, std::vector<std::shared_ptr<Event>> &tmp_event_list, Input &in, IO &out)
        : event_list(tmp_event_list), Observable(path, in, out) {
        struct timeval time;
        gettimeofday(&time, NULL);  //END-TIME
        end = time.tv_sec + (time.tv_usec / 1000000.);
        Reset();
    }

    /// Write relevant information about an observable to the output
    virtual void Write() {
        struct timeval time;
        gettimeofday(&time, NULL);  //END-TIME
        end = time.tv_sec + (time.tv_usec / 1000000.);
        double total_time = end - start;
        double norm = 1.;
        if (total_time != 0)
            norm /= total_time;
        std::vector<double> event_times(event_list.size());
        for (uint32_t i = 0; i < event_times.size(); ++i)
            event_times[i] = event_list[i]->time_spent;

        if (first_time) {
            first_time = 0;
            std::string data_type = "scalar";
            out.CreateGroup(prefix + "block");
            out.Write(prefix + "block/data_type", data_type);
            out.CreateExtendableDataSet(prefix + "block/", "x", total_time);
            for (uint32_t i = 0; i < event_times.size(); ++i) {
                out.CreateGroup(prefix + event_list[i]->name);
                out.Write(prefix + event_list[i]->name + "/data_type", data_type);
                out.CreateExtendableDataSet("/" + prefix + event_list[i]->name + "/", "x", event_times[i]);
            }
        } else {
            out.AppendDataSet(prefix + "block/", "x", total_time);
            for (uint32_t i = 0; i < event_times.size(); ++i)
                out.AppendDataSet("/" + prefix + event_list[i]->name + "/", "x", event_times[i]);
        }

        Reset();
    }
};

#endif  // SIMPIMC_OBSERVABLES_TIME_CLASS_H_
