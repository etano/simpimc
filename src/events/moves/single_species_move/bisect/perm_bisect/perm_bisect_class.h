#ifndef SIMPIMC_MOVES_PERM_BISECT_CLASS_H_
#define SIMPIMC_MOVES_PERM_BISECT_CLASS_H_

#include "../bisect_class.h"

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
  void PermuteBeads(field<std::shared_ptr<Bead>> &b0, field<std::shared_ptr<Bead>> &b1, const Cycle &c)
  {
    // Set permutation type
    perm_type = c.type;

    // Execute the permutation
    uint32_t n_perm = c.part.size();
    for (uint32_t i=0; i<n_perm; i++)
      b0(i)->SetNextBead(b1(c.perm(i)));
    for (uint32_t i=0; i<n_perm; i++)
      b1(i)->SetPrevBead(b0(c.i_perm(i)));
    for (uint32_t i=0; i<n_perm; i++)
      b1(i) = b0(i)->GetNextBead(1);
  }

  /// Assign particle labels to the affected beads
  void AssignParticleLabels()
  {
    for (uint32_t p_i=0; p_i<species->GetNPart(); p_i++) {
      std::shared_ptr<Bead> b(species->GetBead(p_i,bead1-1));
      for (uint32_t b_i=species->bead_loop(bead1-1); b_i<species->GetNBead(); b_i++) {
        species->SetBead(p_i,b_i,b);
        species->GetBead(p_i,b_i)->SetP(p_i);
        b = b->GetNextBead(1);
      }
    }
  }

  /// Accept the move
  virtual void Accept()
  {
    // Accept move, so store things
    for (uint32_t p_i=0; p_i<species->GetNPart(); p_i++) { // todo: can make this more efficient by only restoring touched particles
      species->GetBead(p_i,bead1)->StorePrev();
      species->GetBead(p_i,bead1-1)->StoreNext();
    }
    if (n_perm_part > 1) // only need to reassign particle labels if actual permutation
      AssignParticleLabels();

    // Increment permutation counter
    perm_attempt(n_perm_part-1) += 1;
    perm_accept(n_perm_part-1) += 1;

    Bisect::Accept();
  }

  /// Rejects the move
  virtual void Reject()
  {
    // No need to do some things if bisection isn't attempted
    if (n_perm_part > 0) {
      // Restore things
      for (uint32_t p_i=0; p_i<species->GetNPart(); p_i++) { // TODO: can make this more efficient by only restoring touched particles
        species->GetBead(p_i,bead1)->RestorePrev();
        species->GetBead(p_i,bead1-1)->RestoreNext();
      }

      // Increment counter
      perm_attempt(n_perm_part-1) += 1;
    }

    Bisect::Reject();
  }

  /// Resets the relevant counters
  virtual void Reset()
  {
    // Reset counters
    perm_attempt.zeros();
    perm_accept.zeros();
  }

public:
  /// Constructor instantiates parent class and calls Init
  PermBisect(Path &path, RNG &rng, std::vector<std::shared_ptr<Action>> &action_list, Input &in, IO &out)
    : Bisect(path, rng, action_list, in, out)
  {
    // Read in things
    epsilon = in.GetAttribute<double>("epsilon",1.e-100);
    log_epsilon = log(epsilon);

    // Initiate permutation table
    t.zeros(species->GetNPart(),species->GetNPart());

    // Initiate acceptance ratio counters
    perm_accept.set_size(species->GetNPart());
    perm_attempt.set_size(species->GetNPart());
    Reset();
    Bisect::Reset();
    Move::Reset();
  }

  /// Writes relevant information about the move to the output file
  virtual void Write()
  {
    // Write
    if (first_time) {
      out.Write("/Moves/"+name+"/n_perm_type", n_perm_type);
      out.CreateExtendableDataSet("/Moves/"+name+"/", "perm_attempt", perm_attempt);
      out.CreateExtendableDataSet("/Moves/"+name+"/", "perm_accept", perm_accept);
    } else {
      out.AppendDataSet("/Moves/"+name+"/", "perm_attempt", perm_attempt);
      out.AppendDataSet("/Moves/"+name+"/", "perm_accept", perm_accept);
    }

    Reset();
    Bisect::Write();
  }
};

#endif // SIMPIMC_MOVES_PERM_BISECT_CLASS_H_
