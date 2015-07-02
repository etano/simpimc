#ifndef SIMPIMC_PATH_CLASS_H_
#define SIMPIMC_PATH_CLASS_H_

#include "species_class.h"

/// Controls whether reading/writing new or old objects. This is important for the Metropolis rejection scheme.
typedef enum {OLD_MODE, NEW_MODE} ModeType;

/// Main data holding object. Contains all information about the path integral. Holds a vector of Species objects.
class Path
{
private:
  ModeType mode; ///< Holds the current mode.
public:
  bool approximate; ///< Whether or not to use the approximate exponential
  bool pbc; ///< Whether or not using periodic boundary conditions
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
  uint32_t proc_i; ///< Process index
  uint32_t ref_bead; ///< Index of the current reference point time slice

  field<vec<std::complex<double>>> c_k; ///< Constant defined for each charge density
  field<vec<std::complex<double>>> rho_k; ///< All charge densities
  field<vec<std::complex<double>>> rho_k_c; ///< All charge densities copy
  std::vector<std::map<std::vector<uint32_t>,uint32_t,CompareVec<uint32_t>>> poss_perms; ///< Vector of maps defining all possible permutations for each species
  std::map<std::vector<uint32_t>,uint32_t,CompareVec<uint32_t>>::const_iterator poss_perms_iterator; ///< Iterator over map that defines all possible permutations
  std::vector<bool> perm_sectors_setup; ///< Whether or not permutation sectors are setup for each species
  std::vector<double> mag_ks; ///< Vector holding magnitudes of k vectors
  std::vector<std::shared_ptr<Species>> species_list; ///< Vector holding pointers to all species objects
  std::vector<vec<double>> ks; ///< Vector holding k vectors
  std::vector<vec<int>> k_indices; ///< Vector holding k index vectors
  vec<uint32_t> bead_loop; ///< Helper vector for indexing the periodicity in beta
  vec<double> k_box; ///< Box defined by Brillouin zone
  vec<int> max_k_index; ///< Maximum k index used for each physical dimension

  /// Constructor only sets process index
  Path(uint32_t t_proc_i, Input &in, IO &out, RNG &rng)
   : proc_i(t_proc_i)
  {
    out.CreateGroup("System");
    n_d = in.GetChild("System").GetAttribute<uint32_t>("n_d");
    out.Write("System/n_d",n_d);
    n_bead = in.GetChild("System").GetAttribute<uint32_t>("n_bead");
    out.Write("System/n_bead",n_bead);
    beta = in.GetChild("System").GetAttribute<double>("beta");
    out.Write("System/beta",beta);
    pbc = in.GetChild("System").GetAttribute<bool>("pbc", true);
    out.Write("System/pbc",pbc);
    if (pbc) {
      L = in.GetChild("System").GetAttribute<double>("L");
      iL = 1./L;
      vol = pow(L,n_d);
    } else {
      L = 0.;
      iL = 0.;
      vol = 1.;
    }
    out.Write("System/L",L);

    // Approximate with fast math
    approximate = false;

    // Constants
    tau = beta/(1.*n_bead);
    out.Write("System/tau",tau);

    // Initialize species
    out.CreateGroup("System/Particles");
    std::vector<Input> species_input = in.GetChild("Particles").GetChildList("Species");
    n_species = species_input.size();
    n_part = 0;
    for (uint32_t s_i=0; s_i<n_species; s_i++) {
      species_list.push_back(std::make_shared<Species>(species_input[s_i],out,s_i,n_d,n_bead));
      n_part += species_list[s_i]->n_part;
    }
    out.Write("System/n_part",n_part);

    // Initiate bead looping
    bead_loop.set_size(2*n_bead);
    for (uint32_t b_i = 0; b_i < n_bead; b_i++) {
      bead_loop(b_i) = b_i;
      bead_loop(b_i + n_bead) = bead_loop(b_i);
    }

    // Initiate mode
    SetMode(NEW_MODE);

    // Initialize paths within each species
    for (uint32_t s_i=0; s_i<n_species; s_i++)
      species_list[s_i]->InitPaths(species_input[s_i],out,rng,proc_i,L);

    // Reset k_cut
    k_c = 0.;
    double k_cut = in.GetChild("System").GetAttribute<double>("k_cut",2.*M_PI/(pow(vol,1./n_d)));
    SetupKs(k_cut);

    // Resize permutation containers
    poss_perms.resize(n_species);
    perm_sectors_setup.resize(n_species);
    for (uint32_t i=0; i<n_species; ++i)
      perm_sectors_setup[i] = 0;

    // Initiate some global things
    sign = 1;
    importance_weight = 1;
    ref_bead = 0;
    CalcSign();
  }

  /// Shortcut to a specific bead given its species, particle number, and time slice
  std::shared_ptr<Bead> operator() (const uint32_t s_i, const uint32_t p_i, const uint32_t b_i) { return species_list[s_i]->bead(p_i,bead_loop(b_i)); };

  /// Sets species index by matching the species name string
  void GetSpeciesInfo(const std::string &species, uint32_t &species_i)
  {
    for (uint32_t s_i=0; s_i<n_species; s_i++) {
      if (species_list[s_i]->name == species) {
        species_i = s_i;
        return;
      }
    }
    std::cout << "ERROR: No species of name " << species << " !" << std::endl;
  }

  /// Set the mode to the passed mode
  inline void SetMode(ModeType t_mode) { mode = t_mode; };

  /// Returns the current mode
  inline ModeType GetMode() { return mode; };

  /// Print the current configuration to screen
  void PrintPath()
  {
    SetMode(OLD_MODE);
    for (uint32_t s_i=0; s_i<n_species; ++s_i) {
      for (uint32_t p_i=0; p_i<species_list[s_i]->n_part; ++p_i) {
        for (uint32_t b_i=0; b_i<n_bead; ++b_i) {
          std::cout << s_i << " " << p_i << " " << b_i << " ";
          vec<double> r = GetR((*this)(s_i,p_i,b_i));
          for (uint32_t d_i=0; d_i<n_d; ++d_i) {
            std:: cout << r(d_i) << " ";
          }
          std::cout << std::endl;
        }
      }
    }
    std::cout << std::endl;
  }

  /// Store the given vector of beads' positions
  void StoreR(std::vector<std::shared_ptr<Bead>> &affected_beads)
  {
    for (auto& b: affected_beads)
      b->StoreR();
  }

  /// Restore the given vector of beads' position
  void RestoreR(std::vector<std::shared_ptr<Bead>> &affected_beads)
  {
    for (auto& b: affected_beads)
      b->RestoreR();
  }

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
  inline vec<double> Dr(const vec<double> &r0, const std::shared_ptr<Bead> b1) { return Dr(r0, GetR(b1)); };

  /// Compute the vector between two beads and put it in the box
  inline vec<double> Dr(const std::shared_ptr<Bead> b0, const std::shared_ptr<Bead> b1) { return Dr(GetR(b0), GetR(b1)); };

  /// Compute the vector to the midpoint between two beads and put it in the box
  inline vec<double> RBar(const std::shared_ptr<Bead> b0, const std::shared_ptr<Bead> b1) { return GetR(b1) + 0.5*Dr(b0, b1); };

  /// Compute the distance between two objects
  template<class T>
  inline double MagDr(T &r0, T &r1) { return Mag(Dr(r0,r1)); };

  /// Enforce the periodic boundary condition on a vector
  void PutInBox(vec<double> &r)
  {
    for (uint32_t d_i=0; d_i<n_d; ++d_i)
      r(d_i) -= nearbyint(r(d_i)*iL)*L;
  }

  /// Whether or not to include a k vector
  bool Include(const vec<double> &k, const double k_cut)
  {
    double k2 = dot(k,k);
    if (k2 < k_cut*k_cut && k2 != 0.) {
      if(k(0) > 0.)
        return true;
      else if (n_d >= 2 && (k(0) == 0. && k(1) > 0.))
        return true;
      else if (n_d >= 3 && (k(0) == 0. && k(1) == 0. && k(2) > 0.))
        return true;
      else
        return false;
    } else
      return false;
  }

  /// Setup k vectors
  void SetupKs(const double k_cut)
  {
    if (k_cut <= k_c)
      return;
    else {
      ks.clear();
      k_indices.clear();
      mag_ks.clear();
      k_c = k_cut;
    }

    // Calculate k box
    k_box.set_size(n_d);
    for (uint32_t d_i=0; d_i<n_d; d_i++)
      k_box(d_i) = 2.*M_PI/L;

    // Calculate max k index based on k cutoff
    max_k_index.set_size(n_d);
    for (uint32_t d_i=0; d_i<n_d; d_i++)
      max_k_index(d_i) = (uint32_t) ceil(1.1*k_cut/k_box(d_i));

    // Set up C vector
    c_k.set_size(n_d);
    for (uint32_t d_i=0; d_i<n_d; d_i++)
      c_k(d_i).set_size(2*max_k_index(d_i)+1);

    // Generate all possible combinations and permutations of k indices
    std::vector<int> indices;
    for (uint32_t d_i=0; d_i<n_d; d_i++)
      for (int i=-max_k_index(d_i); i<=max_k_index(d_i); i++)
        indices.push_back(i);
    std::vector<std::vector<int>> tmp_k_indices;
    GenCombPermK(tmp_k_indices, indices, n_d, false, true);

    // Iterate through indices, form k vectors, and include only those that should be included
    for (auto& ik: tmp_k_indices) {
      vec<double> k(n_d);
      vec<int> ki(n_d);
      for (uint32_t d_i=0; d_i<n_d; d_i++) {
        k(d_i) = ik[d_i]*k_box(d_i);
        ki(d_i) = max_k_index(d_i) + ik[d_i];
      }
      if (Include(k, k_cut)) {
        ks.push_back(k);
        k_indices.push_back(ki);
        mag_ks.push_back(sqrt(dot(k,k)));
      }
    }

    // Initiate rho k
    InitRhoK();
  }

  /// Initialize rho_k
  void InitRhoK()
  {
    // Resize rho_k
    rho_k.set_size(n_bead,n_species);
    rho_k_c.set_size(n_bead,n_species);

    // Set to old mode
    ModeType t_mode = GetMode();
    SetMode(OLD_MODE);

    // Zero out rho_k array
    for (uint32_t b_i=0; b_i<n_bead; b_i++) {
      for (uint32_t s_i=0; s_i<n_species; s_i++) {
        rho_k(b_i,s_i).zeros(k_indices.size());
        rho_k_c(b_i,s_i).zeros(k_indices.size());
        for (uint32_t p_i=0; p_i<species_list[s_i]->n_part; ++p_i) {
          (*this)(s_i,p_i,b_i)->rho_k.zeros(k_indices.size());
          (*this)(s_i,p_i,b_i)->rho_k_c.zeros(k_indices.size());
        }
      }
    }

    // Calculate rho_k's
    for (uint32_t b_i=0; b_i<n_bead; b_i++) {
      for (uint32_t s_i=0; s_i<n_species; s_i++) {
        // Calculate rho_kp values
        for (uint32_t p_i=0; p_i<species_list[s_i]->n_part; p_i++) {
          CalcRhoKP((*this)(s_i,p_i,b_i));
          (*this)(s_i,p_i,b_i) -> RestoreRhoK();
          rho_k_c(b_i,s_i) += (*this)(s_i,p_i,b_i)->rho_k_c;
        }
        rho_k(b_i,s_i) = rho_k_c(b_i,s_i);
      }
    }

    SetMode(t_mode);
  }

  /// Update rho_k for the entire path
  void UpdateRhoK()
  {
    // Calculate rho_k's
    for (uint32_t b_i=0; b_i<n_bead; b_i++) {
      for (uint32_t s_i=0; s_i<n_species; s_i++) {
        // Zero out rho_k array
        rho_k(b_i,s_i).zeros();

        // Calculate rho_kp values
        for (uint32_t p_i=0; p_i<species_list[s_i]->n_part; p_i++) {
          rho_k(b_i,s_i) += (*this)(s_i,p_i,b_i)->rho_k;
        }
      }
    }
  }

  /// Update rho_k for specific particles between time slices b0 and b1
  void UpdateRhoK(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
  {
    ModeType t_mode = GetMode();
    uint32_t skip = 1<<level;

    // Get species numbers
    std::vector<uint32_t> species;
    for (auto& p: particles) {
      uint32_t s_i = p.first;
      uint32_t p_i = p.second;
      species.push_back(s_i);
    }

    // Reset to old copy
    for (auto& s_i: species) {
      for (uint32_t b_i=b0; b_i<b1; b_i+=skip)
        rho_k(bead_loop(b_i),s_i) = rho_k_c(bead_loop(b_i),s_i);
    }

    // Update C values of changed particles
    for (auto& p: particles) {
      uint32_t s_i = p.first;
      uint32_t p_i = p.second;
      for (uint32_t b_i=b0; b_i<b1; b_i+=skip) {
        // Add in new values
        SetMode(NEW_MODE);
        AddRhoKP(rho_k,p_i,b_i,s_i,1);

        // Subtract out old values
        SetMode(OLD_MODE);
        AddRhoKP(rho_k,p_i,b_i,s_i,-1);
      }
    }

    SetMode(t_mode);
  }

  /// Update rho_k_p for specific particles between time slices b0 and b1
  void UpdateRhoKP(const uint32_t b0, const uint32_t b1, const std::vector<std::pair<uint32_t,uint32_t>> &particles, const uint32_t level)
  {
    uint32_t skip = 1<<level;

    // Get species numbers
    std::vector<uint32_t> species, ps;
    for (auto& p: particles) {
      uint32_t s_i = p.first;
      uint32_t p_i = p.second;
      species.push_back(s_i);
      ps.push_back(p_i);
    }

    for (auto& s: species)
      UpdateRhoKP(b0, b1, s, ps, level);
  }

  /// Update rho_k_p for specific particles belonging to species s_i between time slices b0 and b1
  void UpdateRhoKP(const uint32_t b0, const uint32_t b1, const uint32_t s_i, const std::vector<uint32_t> &particles, const uint32_t level)
  {
    ModeType t_mode = GetMode();
    SetMode(NEW_MODE);
    uint32_t skip = 1<<level;

    // Reset to old copy
    for (uint32_t b_i=b0; b_i<b1; b_i+=skip)
      rho_k(bead_loop(b_i),s_i) = rho_k_c(bead_loop(b_i),s_i);

    // Update C values of changed particles
    for (auto& p_i: particles) {
      for (uint32_t b_i=b0; b_i<b1; b_i+=skip) {
        // Calculate new values
        CalcRhoKP((*this)(s_i,p_i,b_i));

        // Add in new values and subtract out old values
        rho_k(bead_loop(b_i),s_i) += (*this)(s_i,p_i,b_i)->rho_k - (*this)(s_i,p_i,b_i)->rho_k_c;
      }
    }

    SetMode(t_mode);
  }

  /// Calculate constant for each charge density
  void CalcC(const vec<double>& r)
  {
    for (uint32_t d_i=0; d_i<n_d; d_i++) {
      std::complex<double> tmp_c_k;
      double phi = r(d_i)*k_box(d_i);
      tmp_c_k = std::complex<double>(cos(phi), sin(phi));
      c_k(d_i)(max_k_index(d_i)) = 1.;
      for (uint32_t k_i=1; k_i<=max_k_index(d_i); k_i++) {
        c_k(d_i)(max_k_index(d_i)+k_i) = tmp_c_k * c_k(d_i)(max_k_index(d_i)+k_i-1);
        c_k(d_i)(max_k_index(d_i)-k_i) = conj(c_k(d_i)(max_k_index(d_i)+k_i));
      }
    }
  }

  /// Add tmp_rho_k to rho_k for a single particle p_i
  void AddRhoKP(field<vec<std::complex<double>>> &tmp_rho_k, const uint32_t p_i, const uint32_t b_i, const uint32_t s_i, const int pm)
  {
    vec<double> r = GetR((*this)(s_i,p_i,b_i));
    CalcC(r);
    for (uint32_t k_i=0; k_i<k_indices.size(); k_i++) {
      vec<int> &ik = k_indices[k_i];
      std::complex<double> factor = pm;
      for (uint32_t d_i=0; d_i<n_d; d_i++)
        factor *= c_k(d_i)(ik(d_i));
      tmp_rho_k(bead_loop(b_i),s_i)(k_i) += factor;
    }
  }

  /// Calculate rho_k for a single bead
  inline void CalcRhoKP(const std::shared_ptr<Bead> b)
  {
    vec<double> r = GetR(b);
    vec<std::complex<double>> &tmp_rho_k = GetRhoK(b);
    CalcC(r);
    for (uint32_t k_i=0; k_i<k_indices.size(); ++k_i) {
      std::complex<double> factor = 1.;
      for (uint32_t d_i=0; d_i<n_d; d_i++)
        factor *= c_k(d_i)(k_indices[k_i](d_i));
      tmp_rho_k(k_i) = factor;
    }
  }

  /// Get rho_k or rho_k_c depending on current mode
  inline field<vec<std::complex<double>>>& GetRhoK() { return mode==NEW_MODE ? (rho_k) : (rho_k_c); };

  /// Get rho_k or rho_k_c for a specific bead depending on current mode
  inline vec<std::complex<double>>& GetRhoK(const std::shared_ptr<Bead> b) { return mode==NEW_MODE ? (b->rho_k) : (b->rho_k_c); };

  /// Set rho_k_c to rho_k for a specific time slice and species
  inline void StoreRhoK(const uint32_t b_i, const uint32_t s_i) { rho_k_c(bead_loop(b_i),s_i) = rho_k(bead_loop(b_i),s_i); };

  /// Set rho_k to rho_k_c for a specific time slice and species
  inline void RestoreRhoK(const uint32_t b_i, const uint32_t s_i) { rho_k(bead_loop(b_i),s_i) = rho_k_c(bead_loop(b_i),s_i); };

  /// Set rho_k_c to rho_k for specific beads
  void StoreRhoKP(std::vector<std::shared_ptr<Bead>> &affected_beads)
  {
    if (k_c > 0)
      for (auto& b: affected_beads)
        b->StoreRhoK();
  }

  /// Set rho_k to rho_k_c for specific beads
  void RestoreRhoKP(std::vector<std::shared_ptr<Bead>> &affected_beads)
  {
    if (k_c > 0)
      for (auto& b: affected_beads)
        b->RestoreRhoK();
  }

  /// Calculate the sign of the current configuration
  int CalcSign()
  {
    sign = 1;
    for (uint32_t s_i=0; s_i<n_species; s_i++) {
      if (species_list[s_i]->fermi) {
        std::vector<uint32_t> cycles;
        SetCycleCount(s_i, cycles);
        for (auto& cycle: cycles)
          sign *= pow(-1,cycle-1);
      }
    }
  }

  /// Get the current cycle counts
  void SetCycleCount(const uint32_t s_i, std::vector<uint32_t> &cycles)
  {
    vec<uint32_t> already_counted(zeros<vec<uint32_t>>(species_list[s_i]->n_part));
    for (uint32_t p_i=0; p_i<species_list[s_i]->n_part; p_i++) {
      if (!already_counted(p_i)) {
        uint32_t cycle_length = 1;
        std::shared_ptr<Bead> b((*this)(s_i,p_i,n_bead-1));
        already_counted(p_i) = 1;
        while (GetNextBead(b,1) != (*this)(s_i,p_i,0)) {
          cycle_length++;
          b = GetNextBead(b,n_bead);
          already_counted(b->p) = 1;
        }
        cycles.push_back(cycle_length);
      }
    }
  }

  /// Get the current permutation sector of a given species
  uint32_t GetPermSector(const uint32_t s_i)
  {
    std::vector<uint32_t> cycles;
    SetCycleCount(s_i, cycles);
    return GetPermSector(s_i, cycles);
  }

  /// Get the current permutation sector of a given species and set the current cycle counts
  uint32_t GetPermSector(const uint32_t s_i, std::vector<uint32_t> &cycles)
  {
    sort(cycles.begin(),cycles.end());
    poss_perms_iterator = poss_perms[s_i].find(cycles);
    if (poss_perms_iterator == poss_perms[s_i].end()) { // TODO: Feel confident enough to remove this
      std::cout << "Broken Permutation: " << std::endl;
      for (auto& cycle: cycles)
        std::cout << cycle << " ";
      std::cout << std::endl;
      exit(1);
    }
    return poss_perms_iterator->second;
  }

  /// Initialize permutation sectors for a given species
  void SetupPermSectors(const uint32_t s_i, const uint32_t sectors_max)
  {
    int n = species_list[s_i]->n_part;
    if (!perm_sectors_setup[s_i]) {
      std::cout << "Setting up permutation sectors..." << std::endl;
      std::vector<int> a;
      a.resize(n);
      for (int i=0; i<n; i++)
        a[i] = 0;
      int k = 1;
      int y = n-1;
      std::vector<std::vector<uint32_t>> tmp_poss_perms;
      while (k != 0 && (sectors_max > poss_perms[s_i].size() || !sectors_max)) {
        int x = a[k-1] + 1;
        k -= 1;
        while (2*x <= y) {
          a[k] = x;
          y -= x;
          k += 1;
        }
        int l = k+1;
        while (x <= y && (sectors_max > poss_perms[s_i].size() || !sectors_max)) {
          a[k] = x;
          a[l] = y;
          std::vector<uint32_t> b;
          for (std::vector<int>::size_type j=0; j!=k+2; j++)
            b.push_back(a[j]);
          tmp_poss_perms.push_back(b);
          x += 1;
          y -= 1;
        }
        a[k] = x+y;
        y = x+y-1;
        std::vector<uint32_t> c;
        for (std::vector<int>::size_type j=0; j!=k+1; j++)
          c.push_back(a[j]);
        tmp_poss_perms.push_back(c);
      }

      int n_sectors = tmp_poss_perms.size();
      for (std::vector<int>::size_type j=0; j != n_sectors; j++)
        poss_perms[s_i][tmp_poss_perms[j]] = j;
      perm_sectors_setup[s_i] = 1;
    }
  }

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

  /// Get dr, dr_p, and drr_p and keep vectors
  inline void DrDrpDrrp(const uint32_t b0, const uint32_t b1, const uint32_t s0, const uint32_t s1, const uint32_t p0, const uint32_t p1, double &r_mag, double &r_p_mag, double &r_r_p_mag, vec<double> &r, vec<double> &r_p, vec<double> &r_r_p)
  {
    //Dr((*this)(s1,p1,b0),(*this)(s0,p0,b0),r);
    //Dr((*this)(s1,p1,b1),(*this)(s0,p0,b1),r_p);

    r = GetR((*this)(s1,p1,b0)) - GetR((*this)(s0,p0,b0));
    r_p = GetR((*this)(s1,p1,b1)) - GetR((*this)(s0,p0,b1));
    for (uint32_t d_i=0; d_i<n_d; ++d_i) {
      r(d_i) -= nearbyint(r(d_i)*iL)*L;
      r_p(d_i) += nearbyint((r(d_i)-r_p(d_i))*iL)*L;
    }
    r_r_p = r - r_p;
    for (uint32_t d_i=0; d_i<n_d; ++d_i)
      r_r_p(d_i) -= nearbyint(r_r_p(d_i)*iL)*L;
    r_mag = mag(r);
    r_p_mag = mag(r_p);
    r_r_p_mag = mag(r_r_p);

  };

};

#endif // SIMPIMC_PATH_CLASS_H_
