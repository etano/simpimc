#ifndef ActionClass_H
#define ActionClass_H

#include "../EventClass.h"
#include "../PathClass.h"

class Action : public Event
{
private:

protected:
  // Path
  Path& path;

  // IO
  IOClass &out;

public:
  // Constructor
  Action(Path &tmpPath, Input &in, IOClass &tmpOut)
    : path(tmpPath), out(tmpOut), isImportanceWeight(0)
  {
    name = in.getAttribute<string>("name");
    type = in.getAttribute<string>("type");
    out.CreateGroup("Actions/"+name);
    out.Write("Actions/"+name+"/type",type);
  }

  string type;
  std::vector<std::string> speciesList;
  bool isImportanceWeight;

  // Functions
  virtual void Init(Input &in) {};
  virtual double DActionDBeta() { return 0.; };
  virtual double GetAction(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level) { return 0.; };
  virtual vec<double> GetActionGradient(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level) { vec<double> zero_vec; zero_vec.zeros(path.nD); return zero_vec; };
  virtual double GetActionLaplacian(const uint b0, const uint b1, const vector<pair<uint,uint> >& particles, const uint level) { return 0.; };
  virtual double Potential() { return 0.; };
  virtual double ImportanceWeight() { return 0; };
  virtual void Write() {};
  virtual void Accept() {};
  virtual void Reject() {};
  void DoEvent() {};

  // FIXME: This only pertains to optimized nodes, but had to put it here for the associated move.
  virtual uint GetParamSet() {};
  virtual uint GetNumParamSets() {};
  virtual void SetParamSet(uint t_iParamSet) {};
  virtual void SetRandomParamSet() {};

};

#endif
