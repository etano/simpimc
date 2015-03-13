#ifndef SIMPIMC_ACTIONS_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_ACTION_CLASS_H_

#include "../path_class.h"

class Action
{
private:

protected:
  // Path
  Path& path;

  // IO
  IO &out;

public:
  // Constructor
  Action(Path &tmp_path, Input &in, IO &tmpOut)
    : path(tmp_path), out(tmpOut), is_importance_weight(0)
  {
    name = in.GetAttribute<std::string>("name");
    type = in.GetAttribute<std::string>("type");
    out.CreateGroup("Actions/"+name);
    out.Write("Actions/"+name+"/type",type);
  }

  std::string name, type;
  std::vector<std::string> species_list;
  bool is_importance_weight;

  // Functions
  virtual void Init(Input &in) {};
  virtual double DActionDBeta() { return 0.; };
  virtual double GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level) { return 0.; };
  virtual vec<double> GetActionGradient(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level) { vec<double> zero_vec; zero_vec.zeros(path.n_d); return zero_vec; };
  virtual double GetActionLaplacian(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level) { return 0.; };
  virtual double Potential() { return 0.; };
  virtual double ImportanceWeight() { return 0; };
  virtual void Write() {};
  virtual void Accept() {};
  virtual void Reject() {};

  // FIXME: This only pertains to optimized nodes, but had to put it here for the associated move.
  virtual uint32_t GetParamSet() {};
  virtual uint32_t GetNumParamSets() {};
  virtual void SetParamSet(uint32_t t_param_set_i) {};
  virtual void SetRandomParamSet() {};

};

#endif // SIMPIMC_ACTIONS_ACTION_CLASS_H_
