#ifndef SIMPIMC_BEAD_CLASS_H_
#define SIMPIMC_BEAD_CLASS_H_

/// Class that contains all information about a single bead
class Bead
{
public:
  bool is_ira; ///< Whether or not bead is the head of a worm
  bool is_masha; ///< Whether or not bead is the tail of a worm
  uint32_t p; ///< Particle index
  const uint32_t b; ///< Bead index
  const uint32_t s; ///< Species index

  std::shared_ptr<Bead> next; ///< Pointer to next bead
  std::shared_ptr<Bead> next_c; ///< Copy of pointer to next bead
  std::shared_ptr<Bead> prev; ///< Pointer to previous bead
  std::shared_ptr<Bead> prev_c; ///< Copy of point to previous bead
  vec<double> r; ///< Position vector
  vec<double> r_c; ///< Copy of position vector
  vec<std::complex<double>> rho_k; ///< Single bead contribution to charge density
  vec<std::complex<double>> rho_k_c; ///< Copy of single bead contribution to charge density

  /// Constructor for new bead given spatial dimension, species, particle, and time slice
  Bead(const uint32_t n_d, const uint32_t t_s, const uint32_t t_p, const uint32_t t_b)
    : s(t_s), p(t_p), b(t_b), r(n_d), r_c(n_d), is_ira(false), is_masha(false)
  {}

  /// Sets r_c to r
  inline void StoreR() { r_c = r; };

  /// Sets r to r_c
  inline void RestoreR() { r = r_c; };

  /// Sets rho_k_c to rho_k
  inline void StoreRhoK() { rho_k_c = rho_k; };

  /// Sets rho_k to rho_k_c
  inline void RestoreRhoK() { rho_k = rho_k_c; };

  /// Sets prev_c to prev
  inline void StorePrev() { prev_c = prev; };

  /// Sets prev to prev_c
  inline void RestorePrev() { prev = prev_c; };

  /// Sets next_c to next
  inline void StoreNext() { next_c = next; };

  /// Sets next to next_c
  inline void RestoreNext() { next = next_c; };

  /// Stores all variables to copies
  inline void Store()
  {
    StoreR();
    StoreRhoK();
    StorePartRecord();
  }

  /// Restores all variables from copies
  inline void Restore()
  {
    RestoreR();
    RestoreRhoK();
    RestorePartRecord();
  }

  /// Stores next and prev to next_c and prev_c
  inline void StorePartRecord()
  {
    next_c = next;
    prev_c = prev;
  }

  /// Restores next and prev from next_c and prev_c
  inline void RestorePartRecord()
  {
    next = next_c;
    prev = prev_c;
  }

  /// Returns bead n time steps further up the path as determined by next
  inline std::shared_ptr<Bead> NextB(const uint32_t n)
  {
    return n==1 ? next : next->NextB(n-1);
  }

  /// Returns bead n time steps further up the path as determined by next_c
  inline std::shared_ptr<Bead> NextBC(const uint32_t n)
  {
    return n==1 ? next_c : next_c->NextBC(n-1);
  }

  /// Returns bead n time steps before as determined by prev
  inline std::shared_ptr<Bead> PrevB(const uint32_t n)
  {
    return n==1 ? prev : prev->PrevB(n-1);
  }

  /// Returns bead n time steps before as determined by prev_c
  inline std::shared_ptr<Bead> PrevBC(const uint32_t n)
  {
    return n==1 ? prev_c : prev_c->PrevBC(n-1);
  }

};

#endif // SIMPIMC_BEAD_CLASS_H_
