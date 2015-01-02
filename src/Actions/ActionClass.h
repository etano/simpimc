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
  virtual double DActionDBeta() {};
  virtual double GetAction(const int b0, const int b1, const vector< pair<int,int> >& particles, const int level) {};
  virtual double Potential() { return 0.; };
  virtual double ImportanceWeight() { return 0; };
  virtual void Write() {};
  virtual void Accept() {};
  virtual void Reject() {};
  void DoEvent() {};

  // FIXME: This only pertains to optimized nodes, but had to put it here for the associated move.
  virtual int GetParamSet() {};
  virtual int GetNumParamSets() {};
  virtual void SetParamSet(int t_iParamSet) {};
  virtual void SetRandomParamSet() {};

};

#endif
