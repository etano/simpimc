#ifndef ObservableClass_H
#define ObservableClass_H

#include "../EventClass.h"
#include "../PathClass.h"
#include "../Utils/IO/InputClass.h"
#include "../Utils/IO/IOClass.h"

class Observable : public Event
{
private:

protected:
  Path& path;
  IOClass &out;
  int skip;

  struct LinearGrid
  {
    Tvector rs;
    RealType rMin, rMax, dr, iDr;
    int nR;

    void CreateGrid(RealType t_rMin, RealType t_rMax, int t_nR)
    {
      rMin = t_rMin;
      rMax = t_rMax;
      nR = t_nR;
      dr = (rMax-rMin)/(nR-1.);
      iDr = 1./dr;
      rs.set_size(nR);
      for (int i=0; i<nR; ++i)
        rs(i) = rMin + i*dr;
    };
    inline RealType operator() (int i) { return rs(i); };
    int ReverseMap(RealType r) {
      int i = (int) nearbyint((r-rMin)*iDr-0.5);
      if (i < 0)
        return 0;
      else
        return i;
    };
  };

  struct Histogram
  {
    Tvector y;
    LinearGrid x;
  };

public:
  // Constructor
  Observable(Path &tmpPath, Input &in, IOClass &tmpOut)
    : Event(), path(tmpPath), out(tmpOut)
  {
    name = in.getAttribute<string>("name");
    type = in.getAttribute<string>("type");
    skip = in.getAttribute<int>("skip",1);
    out.CreateGroup("Observables/"+name);
    out.Write("Observables/"+name+"/type",type);
    firstTime = 1;
    nMeasure = 0;
  }

  string type;
  bool firstTime;
  int nMeasure;

  // Functions
  inline void DoEvent() {
    struct timeval time;
    gettimeofday(&time, NULL); // Start Time
    RealType start = time.tv_sec + (time.tv_usec / 1000000.);
    Accumulate();
    gettimeofday(&time, NULL); //END-TIME
    RealType end = time.tv_sec + (time.tv_usec / 1000000.);
    timeSpent += end - start;
  }
  virtual void Init(Input &in) {};
  virtual void Accumulate() {};
  virtual void Reset() {};
  virtual void Write() {};
};

#endif
