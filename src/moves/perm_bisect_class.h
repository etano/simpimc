#ifndef SIMPIMC_MOVES_PERM_BISECT_CLASS_H_
#define SIMPIMC_MOVES_PERM_BISECT_CLASS_H_

#include "bisect_class.h"

/// Permuting bisection move parent class
class PermBisect : public Bisect
{
protected:
  /// Cycle structure for a permutation
  struct Cycle
  {
    double weight; ///< Sampling weight of the cycle
    double contribution; ///< Sum of all weights up to and including this cycle
    uint32_t index; ///< Index of cycle
    uint32_t type; ///< Type of cycle
    vec<uint32_t> perm; ///< Vector of permuted indices
    vec<uint32_t> i_perm; ///< Inverse vector of permuted indices
    vec<uint32_t> part; ///< Vector of particles involved in the cycle
  };

  double epsilon; ///< Tolerance when making the permutation table
  double log_epsilon; ///< Log of tolerance when making the permutation table
  uint32_t n_perm_part; ///< Number of particles involved in a permutation
  uint32_t n_perm_type; ///< Number of different types of permutations given n_perm_part
  uint32_t perm_type; ///< Type of permutation being performed
  vec<uint32_t> perm_accept; ///< Vector of permutations accepted for each permutation type
  vec<uint32_t> perm_attempt; ///< Vector of permutations attempted for each permutation type
  std::vector<Cycle*> cycles; ///< Vector of all possible cycles
  mat<double> t; ///< Table of possible permutation weights
  field<Cycle> all_cycles; ///< All possible cycles

  /// Permute the beads in the cycle
  void PermuteBeads(field<std::shared_ptr<Bead>> &b0, field<std::shared_ptr<Bead>> &b1, const Cycle &c);

  /// Assign particle labels to the affected beads
  void AssignParticleLabels();

  /// Accept the move
  virtual void Accept();

  /// Initializes the move
  virtual void Init(Input &in);

  /// Rejects the move
  virtual void Reject();

  /// Resets the relevant counters
  virtual void Reset();
public:
  /// Constructor instantiates parent class and calls Init
  PermBisect(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Bisect(path, rng, action_list, in, out)
  {
    Init(in);
  }

  /// Writes relevant information about the move to the output file
  virtual void Write();
};


#endif // SIMPIMC_MOVES_PERM_BISECT_CLASS_H_
