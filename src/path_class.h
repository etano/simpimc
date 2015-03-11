#ifndef SIMPIMC_PATH_CLASS_H_
#define SIMPIMC_PATH_CLASS_H_

#include "species_class.h"

class Path
{
private:

protected:

public:
  // Constructor
  Path(Communicator &tmp_world_comm, Communicator &tmp_inter_comm, Communicator &tmp_intra_comm)
   : world_comm(tmp_world_comm), inter_comm(tmp_inter_comm), intra_comm(tmp_intra_comm)
  {}

  void Init(Input &in, IO &out, RNG &rng);

  // Parallel communicators
  Communicator& world_comm; // This is the global MPIWORLD communicator.
  Communicator& inter_comm; // This is for communication between the rank 0 procs of each walker group.
  Communicator& intra_comm; // This is for commmunication between procs within a walker group.

  // Given global constants
  uint n_part, n_d, n_bead;
  double beta, L, iL, vol;

  // Calculated global constants
  double tau;
  uint max_level;

  // Species
  uint n_species;
  std::vector<std::shared_ptr<Species>> species_list;
  void GetSpeciesInfo(const std::string &species, uint &species_i);

  // Fast math
  bool approximate;
  inline double FastExp(const double x) { return exp(x); };

  // Mode (use copy or true)
  bool mode;
  void SetMode(uint m) { mode = m; };
  bool GetMode() { return mode; };

  // Beads
  void PrintPath();
  vec<uint> bead_loop;
  std::shared_ptr<Bead> operator() (const uint s_i, const uint p_i, const uint b_i) { return species_list[s_i]->bead(p_i,bead_loop(b_i)); };
  void StoreR(std::vector<std::shared_ptr<Bead>> &affected_beads);
  void RestoreR(std::vector<std::shared_ptr<Bead>> &affected_beads);
  inline vec<double>& GetR(const std::shared_ptr<Bead> b) { return mode ? b->r : b->r_c; };
  inline std::shared_ptr<Bead> GetNextBead(const std::shared_ptr<Bead> b, const uint n) { return mode ? b->NextB(n) : b->NextBC(n); };
  inline std::shared_ptr<Bead> GetPrevBead(const std::shared_ptr<Bead> b, const uint n) { return mode ? b->PrevB(n) : b->PrevBC(n); };
  inline vec<double> Dr(const vec<double> &r0, const vec<double> &r1) { vec<double> dr = r0-r1; PutInBox(dr); return dr; };
  inline vec<double> Dr(const std::shared_ptr<Bead> b0, const vec<double> &r1) { return Dr(GetR(b0), r1); };
  inline vec<double> Dr(const std::shared_ptr<Bead> b0, const std::shared_ptr<Bead> b1) { return Dr(GetR(b0), GetR(b1)); };
  inline vec<double> RBar(const std::shared_ptr<Bead> b0, const std::shared_ptr<Bead> b1) { return GetR(b1) + 0.5*Dr(b0, b1); };
  template<class T>
  inline double MagDr(T &r0, T &r1) { return Mag(Dr(r0,r1)); };

  // Periodic Boundary Condition
  bool pbc;
  void PutInBox(vec<double> &r);

  // k vectors and rho_k
  std::vector<vec<double>> ks;
  std::vector<vec<int>> k_indices;
  std::vector<double> mag_ks;
  field<vec<std::complex<double>>> c_k;
  field<vec<std::complex<double>>> rho_k, rho_k_c;
  vec<double> k_box;
  double k_c;
  vec<int> max_k_index;
  bool Include(const vec<double> &k, const double kCut);
  void SetupKs(const double kCut);
  void InitRhoK();
  void UpdateRhoK();
  void UpdateRhoK(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level);
  void UpdateRhoKP(const uint b0, const uint b1, const std::vector<std::pair<uint,uint>> &particles, const uint level);
  void UpdateRhoKP(const uint b0, const uint b1, const uint s_i, const std::vector<uint> &particles, const uint level);
  void CalcC(const vec<double>& r);
  void AddRhoKP(field<vec<std::complex<double>>> &tmp_rho_k, const uint p_i, const uint b_i, const uint s_i, const int pm);
  inline void CalcRhoKP(const std::shared_ptr<Bead> b);
  inline field<vec<std::complex<double>>>& GetRhoK() { return mode ? (rho_k) : (rho_k_c); };
  inline vec<std::complex<double>>& GetRhoK(const std::shared_ptr<Bead> b) { return mode ? (b->rho_k) : (b->rho_k_c); };
  inline void StoreRhoK(const uint b_i, const uint s_i) { rho_k_c(bead_loop(b_i),s_i) = rho_k(bead_loop(b_i),s_i); };
  inline void RestoreRhoK(const uint b_i, const uint s_i) { rho_k(bead_loop(b_i),s_i) = rho_k_c(bead_loop(b_i),s_i); };
  void StoreRhoKP(std::vector<std::shared_ptr<Bead>> &affected_beads);
  void RestoreRhoKP(std::vector<std::shared_ptr<Bead>> &affected_beads);

  // Importance weight
  double importance_weight;
  int sign;

  // Nodes
  uint ref_bead;

  // Permutations
  int CalcSign();
  struct CompareVecInt
  {
    bool operator() (const std::vector<uint> &a, const std::vector<uint> &b) {
      for (uint i=0; i<a.size(); i++)
        if (a[i] != b[i])
          return (a[i] > b[i]);
      return (a[0] > b[0]);
    }
  };

  std::map<std::vector<uint>,uint,CompareVecInt> poss_perms;
  std::map<std::vector<uint>,uint,CompareVecInt>::const_iterator poss_perms_iterator;
  void SetCycleCount(const uint s_i, std::vector<uint> &cycles);
  uint GetPermSector(const uint s_i);
  uint GetPermSector(const uint s_i, std::vector<uint> &cycles);
  void SetupPermSectors(const uint n, const uint sectors_max);
  bool perm_sectors_setup;


  // Path initialization
  void InitPaths(Input &in, IO &out, RNG &rng);

  // Get dr, dr_p, and drr_p
  inline void DrDrpDrrp(const uint b0, const uint b1, const uint s0, const uint s1, const uint p0, const uint p1, double &r_mag, double &r_p_mag, double &r_r_p_mag)
  {
    //Dr((*this)(s1,p1,b0),(*this)(s0,p0,b0),r);
    //Dr((*this)(s1,p1,b1),(*this)(s0,p0,b1),r_p);

    vec<double> r = GetR((*this)(s1,p1,b0)) - GetR((*this)(s0,p0,b0));
    vec<double> r_p = GetR((*this)(s1,p1,b1)) - GetR((*this)(s0,p0,b1));
    for (uint d_i=0; d_i<n_d; ++d_i) {
      r(d_i) -= nearbyint(r(d_i)*iL)*L;
      r_p(d_i) += nearbyint((r(d_i)-r_p(d_i))*iL)*L;
    }
    vec<double> r_r_p = r - r_p;
    for (uint d_i=0; d_i<n_d; ++d_i)
      r_r_p(d_i) -= nearbyint(r_r_p(d_i)*iL)*L;
    r_mag = mag(r);
    r_p_mag = mag(r_p);
    r_r_p_mag = mag(r_r_p);

  };
};


#endif // SIMPIMC_PATH_CLASS_H_
