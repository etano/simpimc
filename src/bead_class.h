#ifndef SIMPIMC_BEAD_CLASS_H_
#define SIMPIMC_BEAD_CLASS_H_

#include "config.h"

struct Bead
{
public:
  Bead();
  Bead(uint32_t n_d, uint32_t t_s, uint32_t t_p, uint32_t t_b)
    : s(t_s), p(t_p), b(t_b), r(n_d), r_c(n_d), is_ira(false), is_masha(false)
  {}

  uint32_t p, b, s;
  bool is_ira, is_masha; // head and tail (respectively)
  double n_dist, n_dist_c;
  vec<double> r, r_c;
  vec<std::complex<double>> rho_k, rho_k_c;
  std::shared_ptr<Bead> next, next_c, prev, prev_c;

  inline void StoreR() { r_c = r; };
  inline void RestoreR() { r = r_c; };
  inline void StoreRhoK() { rho_k_c = rho_k; };
  inline void RestoreRhoK() { rho_k = rho_k_c; };
  inline void StorePrev() { prev_c = prev; };
  inline void RestorePrev() { prev = prev_c; };
  inline void StoreNext() { next_c = next; };
  inline void RestoreNext() { next = next_c; };
  inline void StoreNodeDistance() { n_dist_c = n_dist; };
  inline void RestoreNodeDistance() { n_dist = n_dist_c; };
  inline void Move(const vec<double>& dr) { r += dr; };

  inline void Store()
  {
    StoreR();
    StoreRhoK();
    StorePartRecord();
    StoreNodeDistance();
  }

  inline void Restore()
  {
    RestoreR();
    RestoreRhoK();
    RestorePartRecord();
    RestoreNodeDistance();
  }

  inline void StorePartRecord()
  {
    next_c = next;
    prev_c = prev;
  }

  inline void RestorePartRecord()
  {
    next = next_c;
    prev = prev_c;
  }

  std::shared_ptr<Bead> NextB(const uint32_t n)
  {
    std::shared_ptr<Bead> bead(prev->next);
    for (uint32_t i=0; i<n; i++) bead = bead->next;
    return bead;
  }

  std::shared_ptr<Bead> NextBC(const uint32_t n)
  {
    std::shared_ptr<Bead> bead(prev->next);
    for (uint32_t i=0; i<n; i++) bead = bead->next_c;
    return bead;
  }

  std::shared_ptr<Bead> PrevB(const uint32_t n)
  {
    std::shared_ptr<Bead> bead(next->prev);
    for (uint32_t i=0; i<n; i++) bead = bead->prev;
    return bead;
  }

  std::shared_ptr<Bead> PrevBC(const uint32_t n)
  {
    std::shared_ptr<Bead> bead(next->prev);
    for (uint32_t i=0; i<n; i++) bead = bead->prev_c;
    return bead;
  }

};

#endif // SIMPIMC_BEAD_CLASS_H_
