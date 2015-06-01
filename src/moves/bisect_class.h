#ifndef SIMPIMC_MOVES_BISECT_CLASS_H_
#define SIMPIMC_MOVES_BISECT_CLASS_H_

#include "single_species_move_class.h"
#include "../actions/free_spline_class.h"

/// Performs a bisection move
class Bisect : public SingleSpeciesMove
{
protected:
  bool adaptive; ///< Whether or not to dynamically adjust the bisection level
  double i_4_lambda_tau_n_bisect_beads; ///< 1/(4\lambda\tau*n_bisect_beads) for affected species
  double target_ratio; ///< Target accepted ratio if using adaptive bisection level
  int n_images; ///< Number of periodic images used for the bisection sampling
  int ref_accept; ///< Counter tracking the number of bisections involving the reference slice are accepted
  int ref_attempt; ///< Counter tracking the number of bisections involving the reference slice are attempted
  uint32_t bead0; ///< First time slice in bisection
  uint32_t bead1; ///< Last time slice in bisection
  uint32_t n_bisect_beads; ///< Number of beads involved in a bisection of level n_level
  uint32_t n_level; ///< Number of levels in the bisection
  std::vector<FreeSpline> rho_free_splines; ///< Holds the splined action for every time slice
  std::vector<std::shared_ptr<Bead>> affected_beads; ///< Vector of the beads affected by the current bisection

  /// Accept the move
  virtual void Accept();

  /// Attempt the move
  virtual bool Attempt();

  /// Initializes the move
  virtual void Init(Input &in);

  /// Rejects the move
  virtual void Reject();

  /// Resets the relevant counters
  virtual void Reset();

  /// Creates splined action for all time slices
  virtual void SetupSpline();
public:
  /// Constructor instantiates parent and calls Init
  Bisect(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : SingleSpeciesMove(path, rng, action_list, in, out)
  {
    Init(in);
  }

  /// Writes relevant information about the move to the output file
  virtual void Write();
};

#endif // SIMPIMC_MOVES_BISECT_CLASS_H_
