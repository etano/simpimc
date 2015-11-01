#ifndef SIMPIMC_PATH_CLASS_H_
#define SIMPIMC_PATH_CLASS_H_

#include "species_class.h"

/// Main data holding object. Contains all information about the path integral. Holds a vector of Species objects.
class Path
{
private:
  ModeType mode_; ///< Holds the current mode.
  double L_; ///< Cubic box side length
  double iL_; ///< Inverse of cubic box side length
  double vol_; ///< Volume of the simulation cell
  double surface_; ///< Surface area of the simulation cell
  bool pbc_; ///< Whether or not using periodic boundary conditions
  double beta_; ///< Inverse temperature: 1/(k_B T)
  double importance_weight_; ///< Importance weight of current configuration
  double tau_; ///< Time step
  uint32_t n_bead_; ///< Total number of time slices
  uint32_t n_d_; ///< Number of spatial dimensions
  std::vector<std::shared_ptr<Species>> species_; ///< Vector holding pointers to all species objects
public:
  KSpace ks; ///< K space data container

  /// Constructor sets system specific variables
  Path(uint32_t proc_i, Input &in, IO &out, RNG &rng)
   : n_d_(in.GetChild("System").GetAttribute<uint32_t>("n_d")),
     n_bead_(in.GetChild("System").GetAttribute<uint32_t>("n_bead")),
     beta_(in.GetChild("System").GetAttribute<double>("beta")),
     pbc_(in.GetChild("System").GetAttribute<bool>("pbc", true)),
     L_(in.GetChild("System").GetAttribute<double>("L", 0.)),
     importance_weight_(1.)
  {
    // Computed constants
    if (pbc_) {
      iL_ = 1./L_;
      vol_ = pow(L_,n_d_);
      surface_ = 2*n_d_*pow(L_,n_d_-1);
    } else {
      L_ = 0.;
      iL_ = 0.;
      vol_ = 1.;
      surface_ = 1.;
    }
    tau_ = beta_/(1.*n_bead_);

    // Write out constants
    out.CreateGroup("System");
    out.Write("System/n_d",n_d_);
    out.Write("System/n_bead",n_bead_);
    out.Write("System/beta",beta_);
    out.Write("System/pbc",pbc_);
    out.Write("System/L",L_);
    out.Write("System/tau",tau_);

    // Initiate mode
    SetMode(NEW_MODE);

    // Setup k space
    ks.n_d = n_d_;
    ks.L = L_;
    ks.cutoff = 0.;
    double k_cut = in.GetChild("System").GetAttribute<double>("k_cut",2.*M_PI/L_);
    ks.Setup(k_cut);

    // Initialize species
    out.CreateGroup("System/Particles");
    std::vector<Input> species_input = in.GetChild("Particles").GetChildList("Species");
    for (uint32_t s_i=0; s_i<species_input.size(); s_i++)
      species_.push_back(std::make_shared<Species>(species_input[s_i],out,s_i,n_d_,n_bead_,mode_,ks,rng,proc_i));

  }

  /// Return L
  const double GetL() { return L_; }

  /// Return 1/L
  const double GetInverseL() { return iL_; }

  /// Return volume
  const double GetVol() { return vol_; }

  /// Return surface
  const double GetSurface() { return surface_; }

  /// Sets species index by matching the species name string
  std::shared_ptr<Species> GetSpecies(const std::string &species_name)
  {
    for (const auto& species : species_)
      if (species->GetName() == species_name)
        return species;
    std::cerr << "ERROR: No species of name " << species_name << " !" << std::endl;
  }

  /// Sets species index by matching the species name string
  std::vector<std::shared_ptr<Species>>& GetSpecies() { return species_; }

  /// Set the mode to the passed mode
  void SetMode(ModeType mode) { mode_ = mode; }

  /// Returns the current mode
  ModeType GetMode() { return mode_; };

  /// Print the current configuration to screen
  void Print()
  {
    for (const auto& species : species_)
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
    for (uint32_t d_i=0; d_i<n_d_; ++d_i)
      r(d_i) -= nearbyint(r(d_i)*iL_)*L_;
  }

  /// Get dr, dr_p, and drr_p
  void DrDrpDrrp(const uint32_t b0, const uint32_t b1, const std::shared_ptr<Species>& s0, const std::shared_ptr<Species>& s1, const uint32_t p0, const uint32_t p1, double &r_mag, double &r_p_mag, double &r_r_p_mag)
  {
    vec<double> r = s1->GetBead(p1,b0)->GetR() - s0->GetBead(p0,b0)->GetR();
    vec<double> r_p = s1->GetBead(p1,b1)->GetR() - s0->GetBead(p0,b1)->GetR();
    for (uint32_t d_i=0; d_i<n_d_; ++d_i) {
      r(d_i) -= nearbyint(r(d_i)*iL_)*L_;
      r_p(d_i) += nearbyint((r(d_i)-r_p(d_i))*iL_)*L_;
    }
    vec<double> r_r_p = r - r_p;
    for (uint32_t d_i=0; d_i<n_d_; ++d_i)
      r_r_p(d_i) -= nearbyint(r_r_p(d_i)*iL_)*L_;
    r_mag = mag(r);
    r_p_mag = mag(r_p);
    r_r_p_mag = mag(r_r_p);
  }

  /// Get dr, dr_p, and drr_p and keep vectors
  void DrDrpDrrp(const uint32_t b0, const uint32_t b1, const std::shared_ptr<Species>& s0, const std::shared_ptr<Species>& s1, const uint32_t p0, const uint32_t p1, double &r_mag, double &r_p_mag, double &r_r_p_mag, vec<double> &r, vec<double> &r_p, vec<double> &r_r_p)
  {
    r = s1->GetBead(p1,b0)->GetR() - s0->GetBead(p0,b0)->GetR();
    r_p = s1->GetBead(p1,b1)->GetR() - s0->GetBead(p0,b1)->GetR();
    for (uint32_t d_i=0; d_i<n_d_; ++d_i) {
      r(d_i) -= nearbyint(r(d_i)*iL_)*L_;
      r_p(d_i) += nearbyint((r(d_i)-r_p(d_i))*iL_)*L_;
    }
    r_r_p = r - r_p;
    for (uint32_t d_i=0; d_i<n_d_; ++d_i)
      r_r_p(d_i) -= nearbyint(r_r_p(d_i)*iL_)*L_;
    r_mag = mag(r);
    r_p_mag = mag(r_p);
    r_r_p_mag = mag(r_r_p);
  }

  /// Calculate the sign of the current configuration
  int GetSign()
  {
    int sign = 1;
    for (auto& species : species_)
      sign *= species->GetSign();
    return sign;
  }

  /// Return the importance weight of the current configuration
  const double GetImportanceWeight() { return importance_weight_; }

  /// Return the importance weight of the current configuration
  void SetImportanceWeight(const double importance_weight) { importance_weight_ = importance_weight; }

  /// Return the time step
  const double GetTau() { return tau_; }

  /// Return the time step
  const uint32_t GetND() { return n_d_; }

  /// Return the time step
  const uint32_t GetNBead() { return n_bead_; }

  /// Return the time step
  const uint32_t GetNSpecies() { return species_.size(); }

  /// Return if we are using periodic boundary conditions or not
  const bool GetPBC() { return pbc_; }
};

#endif // SIMPIMC_PATH_CLASS_H_
