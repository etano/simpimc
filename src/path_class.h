#ifndef SIMPIMC_PATH_CLASS_H_
#define SIMPIMC_PATH_CLASS_H_

#include "species_class.h"

/// Controls whether reading/writing new or old objects.
/// This is important for the Metropolis rejection scheme.
typedef enum {OLD_MODE, NEW_MODE} ModeType;

/// Main data holding object. Contains all information about the
/// path integral. Holds a vector of Species objects.
class Path
{
private:
  ModeType mode; ///< Holds the current mode.
public:
  bool approximate; ///< Whether or not to use the approximate exponential
  bool pbc; ///< Whether or not using periodic boundary conditions
  bool perm_sectors_setup; ///< Whether or not permutation sectors are setup
  double beta; ///< Inverse temperature: 1/(k_B T)
  double iL; ///< Inverse of cubic box side length
  double importance_weight; ///< Importance weight of the current configuration
  double L; ///< Cubic box side length
  double k_c; ///< Largest cutoff used for k vectors
  double tau; ///< Time step
  double vol; ///< Volume of the simulation cell
  int sign; ///< Sign of the current configuration
  uint32_t n_bead; ///< Total number of time slices
  uint32_t n_d; ///< Number of spatial dimensions
  uint32_t n_part; ///< Total number of particles
  uint32_t n_species; ///< Total number of species of particles
  uint32_t ref_bead; ///< Index of the current reference point time slice

  field<vec<std::complex<double>>> c_k; ///< Constant defined for each charge density
  field<vec<std::complex<double>>> rho_k; ///< All charge densities
  field<vec<std::complex<double>>> rho_k_c; ///< All charge densities copy
  std::map<std::vector<uint32_t>,uint32_t,CompareVec<uint32_t>> poss_perms; ///< Map defining all possible permutations
  std::map<std::vector<uint32_t>,uint32_t,CompareVec<uint32_t>>::const_iterator poss_perms_iterator; ///< Iterator over map that defines all possible permutations
  std::vector<double> mag_ks; ///< Vector holding magnitudes of k vectors
  std::vector<std::shared_ptr<Species>> species_list; ///< Vector holding pointers to all species objects
  std::vector<vec<double>> ks; ///< Vector holding k vectors
  std::vector<vec<int>> k_indices; ///< Vector holding k index vectors
  vec<uint32_t> bead_loop; ///< Helper vector for indexing the periodicity in beta
  vec<double> k_box; ///< Box defined by Brillouin zone
  vec<int> max_k_index; ///< Maximum k index used for each physical dimension

  /// Constructor does nothing
  Path() {}

  /// Initializes the path, including all species objects
  void Init(Input &in, IO &out, RNG &rng, const uint32_t proc_i);

  /// Shortcut to a specific bead given its species, particle number, and time slice
  std::shared_ptr<Bead> operator() (const uint32_t s_i, const uint32_t p_i, const uint32_t b_i) { return species_list[s_i]->bead(p_i,bead_loop(b_i)); };

  /// Sets species index by matching the species name string
  void GetSpeciesInfo(const std::string &species, uint32_t &species_i);

  /// Fast exponential
  // FIXME: Currently this does nothing.
  inline double FastExp(const double x) { return exp(x); };

  /// Set the mode to the passed mode
  inline void SetMode(ModeType t_mode) { mode = t_mode; };

  /// Returns the current mode
  inline ModeType GetMode() { return mode; };

  /// Print the current configuration to screen
  void PrintPath();

  /// Store the given vector of beads' positions
  void StoreR(std::vector<std::shared_ptr<Bead>> &affected_beads);

  /// Restore the given vector of beads' position
  void RestoreR(std::vector<std::shared_ptr<Bead>> &affected_beads);

  /// Get the n_d dimensional vector defining the postion of bead b
  inline vec<double>& GetR(const std::shared_ptr<Bead> b) { return mode==NEW_MODE ? b->r : b->r_c; };

  /// Get the bead n steps after the given bead
  inline std::shared_ptr<Bead> GetNextBead(const std::shared_ptr<Bead> b, const uint32_t n) { return mode==NEW_MODE ? b->NextB(n) : b->NextBC(n); };

  /// Get the bead n steps before the given bead
  inline std::shared_ptr<Bead> GetPrevBead(const std::shared_ptr<Bead> b, const uint32_t n) { return mode==NEW_MODE ? b->PrevB(n) : b->PrevBC(n); };

  /// Compute the vector between two vectors and put it in the box
  inline vec<double> Dr(const vec<double> &r0, const vec<double> &r1) { vec<double> dr = r0-r1; PutInBox(dr); return dr; };

  /// Compute the vector between a bead and a vector and put it in the box
  inline vec<double> Dr(const std::shared_ptr<Bead> b0, const vec<double> &r1) { return Dr(GetR(b0), r1); };

  /// Compute the vector between two beads and put it in the box
  inline vec<double> Dr(const std::shared_ptr<Bead> b0, const std::shared_ptr<Bead> b1) { return Dr(GetR(b0), GetR(b1)); };

  /// Compute the vector to the midpoint between two beads and put it in the box
  inline vec<double> RBar(const std::shared_ptr<Bead> b0, const std::shared_ptr<Bead> b1) { return GetR(b1) + 0.5*Dr(b0, b1); };

  /// Compute the distance between two objects
  template<class T>
  inline double MagDr(T &r0, T &r1) { return Mag(Dr(r0,r1)); };

  /// Enforce the periodic boundary condition on a vector
  void PutInBox(vec<double> &r);

  /// Whether or not to include a k vector
  bool Include(const vec<double> &k, const double kCut);

  /// Setup k vectors
  void SetupKs(const double kCut);

  /// Initialize rho_k
  void InitRhoK();

  /// Update rho_k for the entire path
  void UpdateRhoK();

  /// Update rho_k for specific particles between time slices b0 and b1
  void UpdateRhoK(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Update rho_k_p for specific particles between time slices b0 and b1
  void UpdateRhoKP(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level);

  /// Update rho_k_p for specific particles belonging to species s_i between time slices b0 and b1
  void UpdateRhoKP(const uint32_t b0, const uint32_t b1, const uint32_t s_i, const std::vector<uint32_t> &particles, const uint32_t level);

  /// Calculate constant for each charge density
  void CalcC(const vec<double>& r);

  /// Add tmp_rho_k to rho_k for a single particle p_i
  void AddRhoKP(field<vec<std::complex<double>>> &tmp_rho_k, const uint32_t p_i, const uint32_t b_i, const uint32_t s_i, const int pm);

  /// Calculate rho_k for a single bead
  inline void CalcRhoKP(const std::shared_ptr<Bead> b);

  /// Get rho_k or rho_k_c depending on current mode
  inline field<vec<std::complex<double>>>& GetRhoK() { return mode==NEW_MODE ? (rho_k) : (rho_k_c); };

  /// Get rho_k or rho_k_c for a specific bead depending on current mode
  inline vec<std::complex<double>>& GetRhoK(const std::shared_ptr<Bead> b) { return mode==NEW_MODE ? (b->rho_k) : (b->rho_k_c); };

  /// Set rho_k_c to rho_k for a specific time slice and species
  inline void StoreRhoK(const uint32_t b_i, const uint32_t s_i) { rho_k_c(bead_loop(b_i),s_i) = rho_k(bead_loop(b_i),s_i); };

  /// Set rho_k to rho_k_c for a specific time slice and species
  inline void RestoreRhoK(const uint32_t b_i, const uint32_t s_i) { rho_k(bead_loop(b_i),s_i) = rho_k_c(bead_loop(b_i),s_i); };

  /// Set rho_k_c to rho_k for specific beads
  void StoreRhoKP(std::vector<std::shared_ptr<Bead>> &affected_beads);

  /// Set rho_k to rho_k_c for specific beads
  void RestoreRhoKP(std::vector<std::shared_ptr<Bead>> &affected_beads);

  /// Calculate the sign of the current configuration
  int CalcSign();

  /// Get the current cycle counts
  void SetCycleCount(const uint32_t s_i, std::vector<uint32_t> &cycles);

  /// Get the current permutation sector of a given species
  uint32_t GetPermSector(const uint32_t s_i);

  /// Get the current permutation sector of a given species and set the current cycle counts
  uint32_t GetPermSector(const uint32_t s_i, std::vector<uint32_t> &cycles);

  /// Initialize permutation sectors
  void SetupPermSectors(const uint32_t n, const uint32_t sectors_max);

  /// Get dr, dr_p, and drr_p
  inline void DrDrpDrrp(const uint32_t b0, const uint32_t b1, const uint32_t s0, const uint32_t s1, const uint32_t p0, const uint32_t p1, double &r_mag, double &r_p_mag, double &r_r_p_mag)
  {
    //Dr((*this)(s1,p1,b0),(*this)(s0,p0,b0),r);
    //Dr((*this)(s1,p1,b1),(*this)(s0,p0,b1),r_p);

    vec<double> r = GetR((*this)(s1,p1,b0)) - GetR((*this)(s0,p0,b0));
    vec<double> r_p = GetR((*this)(s1,p1,b1)) - GetR((*this)(s0,p0,b1));
    for (uint32_t d_i=0; d_i<n_d; ++d_i) {
      r(d_i) -= nearbyint(r(d_i)*iL)*L;
      r_p(d_i) += nearbyint((r(d_i)-r_p(d_i))*iL)*L;
    }
    vec<double> r_r_p = r - r_p;
    for (uint32_t d_i=0; d_i<n_d; ++d_i)
      r_r_p(d_i) -= nearbyint(r_r_p(d_i)*iL)*L;
    r_mag = mag(r);
    r_p_mag = mag(r_p);
    r_r_p_mag = mag(r_r_p);

  };
};


#endif // SIMPIMC_PATH_CLASS_H_
