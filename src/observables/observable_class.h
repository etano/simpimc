#ifndef SIMPIMC_OBSERVABLES_OBSERVABLE_CLASS_H_
#define SIMPIMC_OBSERVABLES_OBSERVABLE_CLASS_H_

#include "../event_class.h"
#include "../actions/action_class.h"

/// Parent class for all observables
class Observable : public Event
{
protected:
  bool first_time; ///< Whether or not writing to output for the first time
  uint32_t n_measure; ///< Number of measurements
  uint32_t skip; ///< Number of calls to Accumulate to skip
  IO &out; ///< Reference to the output
  Path &path; ///< Reference to the path
  std::string data_type; ///< Datatype of observable
  std::string prefix; ///< Prefix for all output
  std::string type; ///< Type of observable

  /// A simple linear grid for histograms
  struct LinearGrid
  {
    vec<double> rs; ///< Vector of r values
    double r_min; ///< Minimum r value
    double r_max; ///< Maximum r value
    double dr; ///< Grid spacing
    double d_ir; ///< Inverse grid spacing
    uint32_t n_r; ///< Number of grid points

    /// Create the r grid
    void CreateGrid(double t_r_min, double t_r_max, uint32_t t_n_r)
    {
      r_min = t_r_min;
      r_max = t_r_max;
      n_r = t_n_r;
      dr = (r_max-r_min)/(n_r-1.);
      d_ir = 1./dr;
      rs.set_size(n_r);
      for (uint32_t i=0; i<n_r; ++i)
        rs(i) = r_min + i*dr;
    };

    /// Access an element of the grid
    inline double operator() (uint32_t i) { return rs(i); };

    /// Compute the index corresponding to an r value
    uint32_t ReverseMap(const double r) {
      uint32_t i = (uint32_t) nearbyint((r-r_min)*d_ir-0.5);
      if (i < 0)
        return 0;
      else
        return i;
    };
  };

  /// A simple histogram
  struct Histogram
  {
    vec<double> y; ///< Vector of y values
    LinearGrid x; ///< Linear grid of x values
  };

  /// Initiate the observable
  virtual void Init(Input &in) {};

  /// Accumulate the observable
  virtual void Accumulate() {};

  /// Reset the observable's counters
  virtual void Reset() {};

public:
  /// Constructor sets name, type, skip, prefix
  Observable(Path &t_path, Input &in, IO &t_out, std::string t_data_type="none")
    : Event(), path(t_path), out(t_out), data_type(t_data_type)
  {
    name = in.GetAttribute<std::string>("name");
    type = in.GetAttribute<std::string>("type");
    skip = in.GetAttribute<uint32_t>("skip",1);
    prefix = "Observables/"+name+"/";
    out.CreateGroup(prefix);
    out.Write(prefix+"/type",type);
    out.Write(prefix+"/data_type",data_type);
    first_time = 1;
    n_measure = 0;
  }

  /// Perform and time the accumulation of an observable
  inline void DoEvent() {
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    double start = time.tv_sec + (time.tv_usec / 1000000.);
    Accumulate();
    gettimeofday(&time, NULL); //END-TIME
    double end = time.tv_sec + (time.tv_usec / 1000000.);
    time_spent += end - start;
  }

  /// Write relevant information about an observable to the output
  virtual void Write() {};
};

#endif // SIMPIMC_OBSERVABLES_OBSERVABLE_CLASS_H_
