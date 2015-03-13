#include "nodal_class.h"

double Nodal::DActionDBeta()
{
  return 0.;
}

double Nodal::GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
{
  // Currently old node should be fine
  if (path.mode == 0)
    return 0.;

  // Decide whether or not to check the node
  if (level > max_level || !path.species_list[species_i]->fermi)
    return 0.;
  bool check_node = false;
  for (auto& p: particles) {
    if (p.first == species_i) {
      check_node = true;
      break;
    }
  }
  if (!check_node)
    return 0.;

  // Constants
  uint32_t skip = 1<<level;
  double level_tau = skip*path.tau;

  // See if ref slice included
  bool checkAll = false;
  if (b1 < path.n_bead)
    checkAll = ((b0 <= path.ref_bead) && (b1 >= path.ref_bead));
  else
    checkAll = (path.bead_loop(b1) >= path.ref_bead);

  // Set start and end
  if (checkAll) {
    start_b = 0;
    end_b = path.n_bead-1;
  } else {
    start_b = b0;
    end_b = b1;
  }

  // Set ref beads and initial beads
  std::vector<std::shared_ptr<Bead>> ref_beads, other_beads;
  int slice_diff0 = path.bead_loop(start_b) - path.ref_bead;
  int abs_slice_diff_0 = abs(slice_diff0);
  for (uint32_t p_i=0; p_i<n_part; ++p_i) {
    ref_beads.push_back(path(species_i,p_i,path.ref_bead));
    if (abs_slice_diff_0 >= 0) {
      //other_beads.push_back(path.GetNextBead(ref_beads[p_i],abs_slice_diff_0)); // fixme: This may be the only correct form
      other_beads.push_back(path(species_i,p_i,start_b));
    } else {
      //other_beads.push_back(path.GetPrevBead(ref_beads[p_i],abs_slice_diff_0));
      other_beads.push_back(path(species_i,p_i,start_b));
    }
  }

  // Compute action
  mat<double> g(n_part,n_part);
  double tot = 0.;
  #pragma omp parallel for reduction(+:tot)
  for (uint32_t b_i=start_b; b_i<=end_b; b_i+=skip) {
    if (b_i != path.ref_bead && tot < 1.e100) {
      // Form rho_f
      int slice_diff = path.bead_loop(b_i) - path.ref_bead;
      uint32_t abs_slice_diff = abs(slice_diff);
      uint32_t inv_slice_diff = path.n_bead - abs_slice_diff;
      uint32_t min_slice_diff = std::min(abs_slice_diff, inv_slice_diff);

      for (uint32_t p_i=0; p_i<n_part; ++p_i) {
        for (uint32_t p_j=0; p_j<n_part; ++p_j) {
          vec<double> dr(path.Dr(ref_beads[p_i], other_beads[p_j]));
          g(p_i,p_j) = GetGij(dr, min_slice_diff);
        }
      }
      rho_f(b_i) = det(g);

      // Check sign
      if (rho_f(b_i) < 0.)
        tot += 1.e100;
    }

    // Move down the line
    for (uint32_t p_i=0; p_i<n_part; ++p_i)
      other_beads[p_i] = path.GetNextBead(other_beads[p_i],skip);
  }

  return tot;
}

void Nodal::Accept()
{
  //uint32_t nCheck = n_part;
  //for (uint32_t b_i=start_b; b_i<=end_b; ++b_i)
  //  rho_f_c(b_i) = rho_f(b_i);
}
