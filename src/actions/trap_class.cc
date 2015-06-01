#include "trap_class.h"

double Trap::DActionDBeta()
{
  double tot = 0.;
  #pragma omp parallel for collapse(2) reduction(+:tot)
  for (uint32_t p_i=0; p_i<path.n_part; p_i+=1) {
    for (uint32_t b_i=0; b_i<path.n_bead; b_i+=1) {
      tot += dot(path(species_i,p_i,b_i)->r, path(species_i,p_i,b_i)->r);
    }
  }

  return cofactor_b*tot;
}

double Trap::GetAction(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
{
  if (level > max_level)
    return 0.;

  bool check = false;
  for (uint32_t p=0; p<particles.size(); ++p) {
    if (particles[p].first == species_i) {
      check = true;
      break;
    }
  }
  if (!check)
    return 0.;

  uint32_t skip = 1<<level;
  double tot = 0.;
  for (uint32_t p=0; p<particles.size(); ++p) {
    uint32_t s_i = particles[p].first;
    uint32_t p_i = particles[p].second;
    for (uint32_t b_i=b0; b_i<b1; b_i+=skip) {
      vec<double> dr(path.GetR(path(s_i,p_i,b_i)));
      tot += dot(dr, dr);
    }
  }

  return cofactor_a*skip*tot;
}

void Trap::Write()
{

}
