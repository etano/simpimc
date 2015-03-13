#ifndef SIMPIMC_ACTIONS_DAVID_PAIR_ACTION_CLASS_H_
#define SIMPIMC_ACTIONS_DAVID_PAIR_ACTION_CLASS_H_

#include "pair_action_class.h"

class DavidPairAction : public PairAction
{
private:

public:
  // Constructor
  DavidPairAction(Path &path, Input &in, IO &out)
    : PairAction(path,in,out)
  {
    Init(in);
  }

  // Data
  NUgrid* grid;
  vec<double> taus;
  field<multi_NUBspline_1d_d*> ukj, dukj_dbeta;

  // Functions
  virtual void ReadFile(std::string file_name);

  // Pair actions
  virtual double CalcV(double r, double r_p, const uint32_t level);
  virtual double CalcVLong() { return 0.; };
  virtual double CalcU(double r, double r_p, double s, const uint32_t level);
  virtual double CalcULong(const uint32_t b0, const uint32_t b1, const uint32_t level) { return 0.; };
  virtual double CalcdUdBeta(double r, double r_p, double s, const uint32_t level);
  virtual double CalcdUdBetaLong() { return 0.; };

};

#endif // SIMPIMC_ACTIONS_DAVID_PAIR_ACTION_CLASS_H_
