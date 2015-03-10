#ifndef ObservableClass_H
#define ObservableClass_H

#include "../EventClass.h"
#include "../PathClass.h"

class Observable : public Event
{
private:

protected:
  Path& path;
  IOClass &out;
  uint skip;

  struct LinearGrid
  {
    vec<double> rs;
    double rMin, rMax, dr, iDr;
    uint nR;

    void CreateGrid(double t_rMin, double t_rMax, uint t_nR)
    {
      rMin = t_rMin;
      rMax = t_rMax;
      nR = t_nR;
      dr = (rMax-rMin)/(nR-1.);
      iDr = 1./dr;
      rs.set_size(nR);
      for (uint i=0; i<nR; ++i)
        rs(i) = rMin + i*dr;
    };
    inline double operator() (uint i) { return rs(i); };
    uint ReverseMap(const double r) {
      uint i = (uint) nearbyint((r-rMin)*iDr-0.5);
      if (i < 0)
        return 0;
      else
        return i;
    };
  };

  struct Histogram
  {
    vec<double> y;
    LinearGrid x;
  };

  string prefix;

public:
  // Constructor
  Observable(Path &tmpPath, Input &in, IOClass &tmpOut)
    : Event(), path(tmpPath), out(tmpOut)
  {
    name = in.getAttribute<string>("name");
    type = in.getAttribute<string>("type");
    skip = in.getAttribute<uint>("skip",1);
    prefix = "Observables/"+name+"/";
    out.CreateGroup(prefix);
    out.Write(prefix+"/type",type);
    firstTime = 1;
    nMeasure = 0;
  }

  string type;
  bool firstTime;
  uint nMeasure;

  // Functions
  inline void DoEvent() {
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    double start = time.tv_sec + (time.tv_usec / 1000000.);
    Accumulate();
    gettimeofday(&time, NULL); //END-TIME
    double end = time.tv_sec + (time.tv_usec / 1000000.);
    timeSpent += end - start;
  }
  virtual void Init(Input &in) {};
  virtual void Accumulate() {};
  virtual void Reset() {};
  virtual void Write() {};
};

#endif
