#include "perm_bisect_class.h"

void PermBisect::Init(Input &in)
{
  // Read in things
  epsilon = in.GetAttribute<double>("epsilon",1.e-100);
  log_epsilon = log(epsilon);

  // Initiate permutation table
  t.zeros(n_part,n_part);

  // Initiate acceptance ratio counters
  perm_accept.set_size(n_part);
  perm_attempt.set_size(n_part);
  Reset();
}

// Accept current move
void PermBisect::Accept()
{
  // Change sign weight for fermions
  if (!(n_perm_part%2) && path.species_list[species_i]->fermi)
    path.sign *= -1;

  // Accept move, so store things
  for (uint32_t p_i=0; p_i<n_part; p_i++) { // todo: can make this more efficient by only restoring touched particles
    path(species_i,p_i,bead1)->StorePrev();
    path(species_i,p_i,bead1-1)->StoreNext();
  }
  if (n_perm_part > 1) // only need to reassign particle labels if actual permutation
    AssignParticleLabels();

  // Increment permutation counter
  perm_attempt(perm_type) += 1;
  perm_accept(perm_type) += 1;

  Bisect::Accept();
}

// Reject current move
void PermBisect::Reject()
{
  // No need to do some things if bisection isn't attempted
  if (n_perm_part > 0) {
    perm_attempt(n_perm_part-1) += 1;

    // Restore things
    for (uint32_t p_i=0; p_i<n_part; p_i++) { // TODO: can make this more efficient by only restoring touched particles
      path(species_i,p_i,bead1)->RestorePrev();
      path(species_i,p_i,bead1-1)->RestoreNext();
    }
  }

  Bisect::Reject();
}

// Reset counters
void PermBisect::Reset()
{
  // Reset counters
  perm_attempt.zeros();
  perm_accept.zeros();

  Bisect::Reset();
}

// Permute paths between b0 and b1 given cycle
void PermBisect::PermuteBeads(field<std::shared_ptr<Bead>>& b0, field<std::shared_ptr<Bead>>& b1, const Cycle &c)
{
  // Set permutation type
  perm_type = c.type;

  // Execute the permutation
  uint32_t n_perm = c.part.size();
  for (uint32_t i=0; i<n_perm; i++)
    b0(i)->next = b1(c.perm(i));
  for (uint32_t i=0; i<n_perm; i++)
    b1(i)->prev = b0(c.i_perm(i));
  for (uint32_t i=0; i<n_perm; i++)
    b1(i) = b0(i)->next;

  return;
}

// Reassign particle labels
void PermBisect::AssignParticleLabels()
{
  for (uint32_t p_i=0; p_i<n_part; p_i++) {
    std::shared_ptr<Bead> b(path(species_i,p_i,bead1-1));
    for (uint32_t b_i=path.bead_loop(bead1-1); b_i<path.n_bead; b_i++) {
      path.species_list[species_i]->bead(p_i,b_i) = b;
      path(species_i,p_i,b_i)->p = p_i;
      b = b->next;
    }
  }
}

void PermBisect::Write()
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

  Bisect::Write();
}
