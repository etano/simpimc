#ifndef SIMPIMC_PATH_CLASS_H_
#define SIMPIMC_PATH_CLASS_H_

#include "species_class.h"


/// Main data holding object. Contains all information about the path integral. Holds a vector of Species objects.
class Path
{
private:
  ModeType mode; ///< Holds the current mode.
  double L; ///< Cubic box side length
  double iL; ///< Inverse of cubic box side length
  double vol; ///< Volume of the simulation cell
  double surface; ///<surface of the simulation cell (relevant for boundary term integration)
  bool pbc; ///< Whether or not using periodic boundary conditions
  double beta; ///< Inverse temperature: 1/(k_B T)
  double importance_weight; ///< Importance weight of current configuration
  double tau; ///< Time step
  uint32_t n_bead; ///< Total number of time slices
  uint32_t n_d; ///< Number of spatial dimensions
  uint32_t n_species; ///< Total number of species of particles
  uint32_t proc_i; ///< Process index
public:
  KSpace ks; ///< K space data container
  std::vector<std::shared_ptr<Species>> species_list; ///< Vector holding pointers to all species objects
  vec<uint32_t> bead_loop; ///< Helper vector for indexing the periodicity in beta

  /// Constructor only sets process index
  Path(uint32_t t_proc_i, Input &in, IO &out, RNG &rng)
   : proc_i(t_proc_i),
     n_d(in.GetChild("System").GetAttribute<uint32_t>("n_d")),
     n_bead(in.GetChild("System").GetAttribute<uint32_t>("n_bead")),
     beta(in.GetChild("System").GetAttribute<double>("beta")),
     pbc(in.GetChild("System").GetAttribute<bool>("pbc", true)),
     L(in.GetChild("System").GetAttribute<double>("L", 0.)),
     importance_weight(1.)
  {
    // Computed constants
    if (pbc) {
      iL = 1./L;
      vol = pow(L,n_d);
      surface = 2*n_d*pow(L,n_d-1); 
    } else {
      L = 0.;
      iL = 0.;
      vol = 1.;
      surface=1.;
    }
    tau = beta/(1.*n_bead);

    // Write out constants
    out.CreateGroup("System");
    out.Write("System/n_d",n_d);
    out.Write("System/n_bead",n_bead);
    out.Write("System/beta",beta);
    out.Write("System/pbc",pbc);
    out.Write("System/L",L);
    out.Write("System/tau",tau);

    // Initiate mode
    SetMode(NEW_MODE);

    // Setup k space
    ks.n_d = n_d;
    ks.L = L;
    ks.cutoff = 0.;
    double k_cut = in.GetChild("System").GetAttribute<double>("k_cut",2.*M_PI/L);
    ks.Setup(k_cut);

    // Initialize species
    out.CreateGroup("System/Particles");
    std::vector<Input> species_input = in.GetChild("Particles").GetChildList("Species");
    n_species = species_input.size();
    for (uint32_t s_i=0; s_i<n_species; s_i++)
      species_list.push_back(std::make_shared<Species>(species_input[s_i],out,s_i,n_d,n_bead,mode,ks,rng,proc_i));

    // Init rho_k for each species
    for (auto& species : species_list)
      species->InitRhoK();

    // Initiate bead looping
    bead_loop.set_size(2*n_bead);
    for (uint32_t b_i = 0; b_i < n_bead; b_i++) {
      bead_loop(b_i) = b_i;
      bead_loop(b_i + n_bead) = bead_loop(b_i);
    }

  }

  /// Return L
  const double GetL() { return L; }

  /// Return 1/L
  const double GetInverseL() { return iL; }

  /// Return volume
  const double GetVol() { return vol; }

  /// Return surface
  const double GetSurface() { return surface; }

  /// Sets species index by matching the species name string
  std::shared_ptr<Species> GetSpecies(const std::string &species_name)
  {
    for (const auto& species : species_list)
      if (species->GetName() == species_name)
        return species;
    std::cerr << "ERROR: No species of name " << species_name << " !" << std::endl;
  }

  /// Set the mode to the passed mode
  void SetMode(ModeType t_mode) { mode = t_mode; }

  /// Returns the current mode
  ModeType GetMode() { return mode; };

  /// Print the current configuration to screen
  void Print()
  {
    for (const auto& species : species_list)
      species->Print();
    std::cout << std::endl;
  }

  /// Compute the vector between two vectors and put it in the box
  vec<double> Dr(const vec<double> &r0, const vec<double> &r1) { vec<double> dr(r0-r1); PutInBox(dr); return std::move(dr); };

  /// Compute the vector between a bead and a vector and put it in the box
  vec<double> Dr(const std::shared_ptr<Bead> &b0, const vec<double> &r1) { return Dr(b0->GetR(), r1); };

  /// Compute the vector between two beads and put it in the box
  vec<double> Dr(const vec<double> &r0, const std::shared_ptr<Bead> &b1) { return Dr(r0, b1->GetR()); };

  /// Compute the vector between two beads and put it in the box
  vec<double> Dr(const std::shared_ptr<Bead> &b0, const std::shared_ptr<Bead> &b1) { return Dr(b0->GetR(), b1->GetR()); };

  /// Compute the vector to the midpoint between two beads and put it in the box
  vec<double> RBar(const std::shared_ptr<Bead> &b0, const std::shared_ptr<Bead> &b1) { return b1->GetR() + 0.5*Dr(b0, b1); };

  /// Compute the distance between two objects
  template<class T>
  double MagDr(T &r0, T &r1) { return Mag(Dr(r0,r1)); };

  /// Enforce the periodic boundary condition on a vector
  void PutInBox(vec<double> &r)
  {
    for (uint32_t d_i=0; d_i<n_d; ++d_i)
      r(d_i) -= nearbyint(r(d_i)*iL)*L;
  }

  /// Get dr, dr_p, and drr_p
  void DrDrpDrrp(const uint32_t b0, const uint32_t b1, const std::shared_ptr<Species>& s0, const std::shared_ptr<Species>& s1, const uint32_t p0, const uint32_t p1, double &r_mag, double &r_p_mag, double &r_r_p_mag)
  {
    vec<double> r = s1->GetBead(p1,b0)->GetR() - s0->GetBead(p0,b0)->GetR();
    vec<double> r_p = s1->GetBead(p1,b1)->GetR() - s0->GetBead(p0,b1)->GetR();
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

  }

  /// Get dr, dr_p, and drr_p and keep vectors
  void DrDrpDrrp(const uint32_t b0, const uint32_t b1, const std::shared_ptr<Species>& s0, const std::shared_ptr<Species>& s1, const uint32_t p0, const uint32_t p1, double &r_mag, double &r_p_mag, double &r_r_p_mag, vec<double> &r, vec<double> &r_p, vec<double> &r_r_p)
  {
    r = s1->GetBead(p1,b0)->GetR() - s0->GetBead(p0,b0)->GetR();
    r_p = s1->GetBead(p1,b1)->GetR() - s0->GetBead(p0,b1)->GetR();
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
  }

  /// Calculate the sign of the current configuration
  int GetSign()
  {
    int sign = 1;
    for (auto& species : species_list)
      sign *= species->GetSign();
    return sign;
  }

  /// Return the importance weight of the current configuration
  const double GetImportanceWeight() { return importance_weight; }

  /// Return the importance weight of the current configuration
  void SetImportanceWeight(const double t_importance_weight) { importance_weight = t_importance_weight; }

  /// Return the time step
  const double GetTau() { return tau; }

  /// Return the dimension
  const uint32_t GetND() { return n_d; }

  /// Return the number of beads
  const uint32_t GetNBead() { return n_bead; }

  /// Return the the number of different species
  const uint32_t GetNSpecies() { return n_species; }

  /// Return if we are using periodic boundary conditions or not
  const bool GetPBC() { return pbc; }
};

#endif // SIMPIMC_PATH_CLASS_H_
