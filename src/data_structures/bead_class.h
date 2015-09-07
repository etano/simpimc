#ifndef SIMPIMC_BEAD_CLASS_H_
#define SIMPIMC_BEAD_CLASS_H_

#include "k_space_class.h"

/// Controls whether reading/writing new or old objects. This is important for the Metropolis rejection scheme.
typedef enum {OLD_MODE=0, NEW_MODE=1} ModeType;

/// Class that contains all information about a single bead
class Bead
{
private:
  ModeType& mode_; ///< Holds the current mode_.
  KSpace& ks_; ///< K space data container
  vec<double> r_; ///< Position vector
  vec<double> r_c_; ///< Copy of position vector
  uint32_t p_i_; ///< Particle index
  const uint32_t b_i_; ///< Bead index
  const uint32_t s_i_; ///< Species index
  vec<std::complex<double>> rho_k_; ///< Single bead contribution to charge density
  vec<std::complex<double>> rho_k_c_; ///< Copy of single bead contribution to charge density

  /// Returns bead n time steps further up the path as determined by next
  std::shared_ptr<Bead> NextB(const uint32_t n) { return !n ? prev_->next_ : next_->NextB(n-1); }

  /// Returns bead n time steps further up the path as determined by next_c_
  std::shared_ptr<Bead> NextBC(const uint32_t n) { return !n ? prev_c_->next_c_ : next_c_->NextBC(n-1); }

  /// Returns bead n time steps before as determined by prev
  std::shared_ptr<Bead> PrevB(const uint32_t n) { return !n ? next_->prev_ : prev_->PrevB(n-1); }

  /// Returns bead n time steps before as determined by prev_c_
  std::shared_ptr<Bead> PrevBC(const uint32_t n) { return !n ? next_c_->prev_c_ : prev_c_->PrevBC(n-1); }
public:
  /// Constructor for new bead given spatial dimension, species, particle, and time slice
  Bead(const uint32_t n_d, const uint32_t s_i, const uint32_t p_i, const uint32_t b_i, ModeType& mode, KSpace& ks)
    : s_i_(s_i), p_i_(p_i), b_i_(b_i), r_(n_d), r_c_(n_d), mode_(mode), ks_(ks)
  {}

  std::shared_ptr<Bead> next_; ///< Pointer to next bead
  std::shared_ptr<Bead> next_c_; ///< Copy of pointer to next bead
  std::shared_ptr<Bead> prev_; ///< Pointer to previous bead
  std::shared_ptr<Bead> prev_c_; ///< Copy of point to previous bead


  /// Return the bead index
  const uint32_t GetB() const { return b_i_; }

  /// Return the particle index
  const uint32_t GetP() { return p_i_; }

  /// Return the particle index
  void SetP(const uint32_t p_i) { p_i_ = p_i; }

  /// Return the species index
  const uint32_t GetS() const { return s_i_; }

  /// Sets r_c to r
  void StoreR() { r_c_ = r_; };

  /// Sets r to r_c
  void RestoreR() { r_ = r_c_; };

  /// Sets rho_k_c to rho_k
  void StoreRhoK() { rho_k_c_ = rho_k_; };

  /// Sets rho_k to rho_k_c
  void RestoreRhoK() { rho_k_ = rho_k_c_; };

  /// Sets prev_c_ to prev
  void StorePrev() { prev_c_ = prev_; };

  /// Sets prev to prev_c_
  void RestorePrev() { prev_ = prev_c_; };

  /// Sets next_c_ to next
  void StoreNext() { next_c_ = next_; };

  /// Sets next to next_c_
  void RestoreNext() { next_ = next_c_; };

  /// Stores all variables to copies
  void Store()
  {
    StoreR();
    StoreRhoK();
    StoreNext();
    StorePrev();
  }

  /// Restores all variables from copies
  void Restore()
  {
    RestoreR();
    RestoreRhoK();
    RestoreNext();
    RestorePrev();
  }

  /// Get the n_d dimensional vector defining the postion of bead b
  vec<double>& GetR() { return mode_ ? r_ : r_c_; };

  /// Get the n_d dimensional vector defining the postion of bead b
  void SetR(const vec<double>& r) { mode_ ? r_=r : r_c_=r; };

  /// Get the n_d dimensional vector defining the postion of bead b
  void SetR(const uint32_t d_i, const double x) { mode_ ? r_(d_i)=x : r_c_(d_i)=x; };

  /// Get the bead n steps after the given bead
  std::shared_ptr<Bead> GetNextBead(const uint32_t n) { return mode_ ? NextB(n) : NextBC(n); };

  /// Get the bead n steps before the given bead
  std::shared_ptr<Bead> GetPrevBead(const uint32_t n) { return mode_ ? PrevB(n) : PrevBC(n); };

  /// Set the bead after the given bead
  void SetNextBead(const std::shared_ptr<Bead>& b) { mode_ ? next_=b : next_c_=b; };

  /// Set the bead before the given bead
  void SetPrevBead(const std::shared_ptr<Bead>& b) { mode_ ? prev_=b : prev_c_=b; };

  /// Get rho_k or rho_k_c for a specific bead depending on current mode
  vec<std::complex<double>>& GetRhoK() { return mode_ ? rho_k_ : rho_k_c_; };

  /// Get rho_k or rho_k_c for a specific bead depending on asserted mode
  vec<std::complex<double>>& GetRhoK(const ModeType& mode) { return mode ? rho_k_ : rho_k_c_; };

  /// Calculate rho_k for a single bead
  void CalcRhoK()
  {
    ks_.CalcC(GetR());
    for (uint32_t k_i=0; k_i<ks_.indices.size(); ++k_i) {
      std::complex<double> factor = 1.;
      for (uint32_t d_i=0; d_i<ks_.n_d; d_i++)
        factor *= ks_.c_k(d_i)(ks_.indices[k_i](d_i));
      mode_ ? rho_k_(k_i)=factor : rho_k_c_(k_i)=factor;
    }
  }

  void InitRhoK()
  {
    rho_k_.zeros(ks_.vecs.size());
    rho_k_c_.zeros(ks_.vecs.size());
    CalcRhoK();
  }

  /// Print bead coordinates
  void Print()
  {
    std::cout << p_i_ << " " << b_i_ << " ";
    const vec<double>& r(GetR());
    for(const auto& x : r_)
      std:: cout << x << " ";
    std::cout << std::endl;
  }

};

#endif // SIMPIMC_BEAD_CLASS_H_
