#ifndef SIMPIMC_SPECIES_CLASS_H_
#define SIMPIMC_SPECIES_CLASS_H_

#include "bead_class.h"

/// Container for a single species. Holds all the bead objects of that species.
class Species
{
private:
  ModeType& mode_; ///< Holds the current mode.
  KSpace& ks_; ///< K space data container
  bool fermi_; ///< Whether or not species has Fermi statistics
  bool fixed_node_; ///< Whether or not species uses the fixed node approximation
  double lambda_; ///< h_bar^2/2m
  uint32_t n_bead_; ///< Number of time slices
  uint32_t n_d_; ///< Number of physical dimensions
  uint32_t n_part_; ///< Number of particles in species
  uint32_t s_i_; ///< Species index
  uint32_t ref_bead_; ///< Index of the current reference point time slice
  std::string name_; ///< Name of species
  std::vector<uint32_t> cycles_; ///< Vector of current cycle lengths
  std::map<std::vector<uint32_t>,uint32_t,CompareVec<uint32_t>> poss_perms_; ///< map defining all possible permutations for each species
  field<std::shared_ptr<Bead>> bead_; ///< Container of beads
  field<vec<std::complex<double>>> rho_k_; ///< All charge densities
  field<vec<std::complex<double>>> rho_k_c_; ///< All charge densities copy
  bool need_update_rho_k_; ///< Whether or not rho_k needs to be updated

  /// Set the current cycle counts
  void SetCycleCount()
  {
    vec<uint32_t> already_counted(zeros<vec<uint32_t>>(n_part_));
    cycles_.clear();
    for (uint32_t p_i=0; p_i<n_part_; p_i++) {
      if (!already_counted(p_i)) {
        uint32_t cycle_length = 1;
        std::shared_ptr<Bead> b(bead_(p_i,n_bead_-1));
        already_counted(p_i) = 1;
        while (b->GetNextBead(1) != bead_(p_i,0)) {
          cycle_length++;
          b = b->GetNextBead(n_bead_);
          already_counted(b->GetP()) = 1;
        }
        cycles_.push_back(cycle_length);
      }
    }
  }

public:

  /// Constructor
  Species(Input &in, IO &out, const int s_i, const int n_d, const int n_bead, ModeType& mode, KSpace& ks, RNG &rng, const uint32_t proc_i)
    : s_i_(s_i),
      n_d_(n_d),
      n_bead_(n_bead),
      mode_(mode),
      ks_(ks),
      name_(in.GetAttribute<std::string>("name")),
      n_part_(in.GetAttribute<uint32_t>("n_part")),
      lambda_(in.GetAttribute<double>("lambda")),
      fermi_(in.GetAttribute<bool>("fermi",false)),
      fixed_node_(in.GetAttribute<bool>("fixed_node",false)),
      need_update_rho_k_(true),
      ref_bead_(0)
  {
    // Write to file
    out.CreateGroup("System/Particles/"+name_);
    out.Write("System/Particles/"+name_+"/n_part",n_part_);
    out.Write("System/Particles/"+name_+"/lambda",lambda_);
    out.Write("System/Particles/"+name_+"/fermi",fermi_);
    out.Write("System/Particles/"+name_+"/fixed_node",fixed_node_);

    // Initiate beads
    bead_.set_size(n_part_,n_bead_);
    for (uint32_t p_i=0; p_i<n_part_; p_i++)
      for (uint32_t b_i=0; b_i<n_bead_; b_i++)
        bead_(p_i,b_i) = std::make_shared<Bead>(n_d_,s_i_,p_i,b_i,mode_,ks_);

    // Initiate bead connections
    assert(mode_==NEW_MODE);
    for (uint32_t p_i=0; p_i<n_part_; p_i++) {
      bead_(p_i,0)->SetNextBead(bead_(p_i,1));
      bead_(p_i,0)->SetPrevBead(bead_(p_i,n_bead_-1));
      bead_(p_i,0)->StoreNext();
      bead_(p_i,0)->StorePrev();
      for (uint32_t b_i=1; b_i<n_bead_-1; b_i++) {
        bead_(p_i,b_i)->SetNextBead(bead_(p_i,b_i+1));
        bead_(p_i,b_i)->SetPrevBead(bead_(p_i,b_i-1));
        bead_(p_i,b_i)->StoreNext();
        bead_(p_i,b_i)->StorePrev();
      }
      bead_(p_i,n_bead_-1)->SetNextBead(bead_(p_i,0));
      bead_(p_i,n_bead_-1)->SetPrevBead(bead_(p_i,n_bead_-2));
      bead_(p_i,n_bead_-1)->StoreNext();
      bead_(p_i,n_bead_-1)->StorePrev();
    }

    // Initiate the paths
    InitPaths(in,out,rng,proc_i,ks_.L);

    // Initiate rho_k
    InitRhoK();

  };

  /// Cheap way to do beta periodicity
  const uint32_t bead_loop(const uint32_t b_i) { return b_i>=n_bead_ ? b_i-n_bead_ : b_i; }

  /// Returns a specific bead given its species, particle number, and time slice
  const std::shared_ptr<Bead>& GetBead(const uint32_t p_i, const uint32_t b_i) { return bead_(p_i,bead_loop(b_i)); };

  /// Sets a specific bead given its species, particle number, and time slice
  void SetBead(const uint32_t p_i, const uint32_t b_i, const std::shared_ptr<Bead>& b) { bead_(p_i,bead_loop(b_i)) = b; };

  /// Return the number of beads
  const uint32_t GetNBead() const { return n_bead_; }

  /// Return the number of particles
  const uint32_t GetNPart() const { return n_part_; }

  /// Return the species index
  const uint32_t GetId() const { return s_i_; }

  /// Return the name
  const std::string& GetName() const { return name_; }

  /// Return the reference bead
  const uint32_t GetRefBead() const { return ref_bead_; }

  /// Return the reference bead
  void SetRefBead(const uint32_t ref_bead) { ref_bead_ = ref_bead; }

  /// Return whether or not is a fixed node species
  const bool IsFixedNode() const { return fixed_node_; }

  /// Return whether nor is a fermion species
  const bool IsFermi() const { return fermi_; }

  /// Return whether or not rho_k needs updating
  const bool GetNeedUpdateRhoK() const { return need_update_rho_k_; }

  /// Return whether or not rho_k needs updating
  void SetNeedUpdateRhoK(const bool need_update_rho_k) { need_update_rho_k_=need_update_rho_k; }

  /// Return whether nor is a fermion species
  const double GetLambda() const { return lambda_; }

  /// Return the map of possible permutations
  const std::map<std::vector<uint32_t>,uint32_t,CompareVec<uint32_t>>& GetPossPerms() { return poss_perms_; }

  /// Return the size of poss_perms_
  const uint32_t GetNPerm() { return poss_perms_.size(); }

  /// Calculate the sign of the current configuration
  int GetSign()
  {
    int sign = 1;
    if(fermi_){
      SetCycleCount();
      for (const auto& cycle: cycles_)
        sign *= pow(-1,cycle-1);
    }
    return sign;
  }

  /// Get the current cycle counts
  std::vector<uint32_t>& GetCycleCount()
  {
    SetCycleCount();
    return cycles_;
  }

  /// Get the current permutation sector of a given species
  uint32_t GetPermSector()
  {
    SetCycleCount();
    return GetPermSector(cycles_);
  }

  /// Get the current permutation sector of a given species
  uint32_t GetPermSector(std::vector<uint32_t>& cycles)
  {
    sort(cycles.begin(),cycles.end());
    return poss_perms_.find(cycles)->second;
  }

  /// Initialize permutation sectors for a given species
  void SetupPermSectors(const uint32_t sectors_max)
  {
    if (GetNPerm()>0) {
      std::cout << "Setting up permutation sectors..." << std::endl;
      std::vector<int> a;
      a.resize(n_part_);
      for (int i=0; i<n_part_; i++)
        a[i] = 0;
      int k = 1;
      int y = n_part_-1;
      std::vector<std::vector<uint32_t>> poss_perms;
      while (k != 0 && (sectors_max > poss_perms_.size() || !sectors_max)) {
        int x = a[k-1] + 1;
        k -= 1;
        while (2*x <= y) {
          a[k] = x;
          y -= x;
          k += 1;
        }
        int l = k+1;
        while (x <= y && (sectors_max > poss_perms_.size() || !sectors_max)) {
          a[k] = x;
          a[l] = y;
          std::vector<uint32_t> b;
          for (std::vector<int>::size_type j=0; j!=k+2; j++)
            b.push_back(a[j]);
          poss_perms.push_back(b);
          x += 1;
          y -= 1;
        }
        a[k] = x+y;
        y = x+y-1;
        std::vector<uint32_t> c;
        for (std::vector<int>::size_type j=0; j!=k+1; j++)
          c.push_back(a[j]);
        poss_perms.push_back(c);
      }

      int n_sectors = poss_perms.size();
      for (std::vector<int>::size_type j=0; j != n_sectors; j++)
        poss_perms_[poss_perms[j]] = j;
    }
  }

  /// Initiate paths
  void InitPaths(Input &in, IO &out, RNG &rng, const uint32_t proc_i, const double L)
  {
    // Get type
    std::string init_type = in.GetAttribute<std::string>("init_type");
    out.Write("System/Particles/"+name_+"/init_type",init_type);

    // Read configuration from xyz file
    if (init_type == "File") {
      std::string file_name = in.GetAttribute<std::string>("file_name");
      bool all_beads = in.GetAttribute<bool>("all_beads",false);
      out.Write("System/Particles/"+name_+"/file_name",file_name);
      out.Write("System/Particles/"+name_+"/all_beads",all_beads);
      std::fstream init_file;
      init_file.open(file_name.c_str(), std::ios_base::in);
      if (!init_file.fail()) {
        for (uint32_t p_i=0; p_i<n_part_; ++p_i) {
          if (all_beads) {
            for (uint32_t b_i=0; b_i<n_bead_; ++b_i) {
              vec<double> r(n_d_);
              for (uint32_t d_i=0; d_i<n_d_; ++d_i)
                init_file >> r(d_i);
              bead_(p_i,b_i)->SetR(r);
              bead_(p_i,b_i)->StoreR();
            }
          } else {
            vec<double> r(n_d_);
            for (uint32_t d_i=0; d_i<n_d_; ++d_i)
              init_file >> r(d_i);
            for (uint32_t b_i=0; b_i<n_bead_; ++b_i) {
              bead_(p_i,b_i)->SetR(r);
              bead_(p_i,b_i)->StoreR();
            }
          }
        }
        init_file.close();
      } else {
        std::cout << "ERROR: Init file '" << file_name << "' does not exist." << std::endl;
        exit(1);
      }

    // Random configuration
    } else if (init_type == "Random") {
      double cofactor = in.GetAttribute<double>("cofactor",1.);
      for (uint32_t p_i=0; p_i<n_part_; ++p_i) {
        vec<double> r(n_d_);
        for (uint32_t d_i=0; d_i<n_d_; ++d_i)
          r(d_i) = rng.UnifRand(-cofactor*L/2.,cofactor*L/2.);
        for (uint32_t b_i=0; b_i<n_bead_; ++b_i) {
          bead_(p_i,b_i)->SetR(r);
          bead_(p_i,b_i)->StoreR();
        }
      }

    // BCC Lattice
    } else if (init_type == "BCC") {
      int n_part__per_n_d = (int) ceil (pow(0.5*n_part_, 1.0/n_d_));
      double delta = L/n_part__per_n_d;
      for (uint32_t p_i=0; p_i<n_part_; p_i++) {
        int p = p_i/2;
        vec<int> tmp(n_d_);
        if (n_d_ == 2) {
          tmp(0) = p/n_part__per_n_d;
          tmp(1) = p - tmp(0)*n_part__per_n_d;
        } if (n_d_ == 3) {
          tmp(0) = p/(n_part__per_n_d*n_part__per_n_d);
          tmp(1) = (p-(tmp(0)*n_part__per_n_d*n_part__per_n_d))/n_part__per_n_d;
          tmp(2) = p - tmp(0)*n_part__per_n_d*n_part__per_n_d - tmp(1)*n_part__per_n_d;
        }
        vec<double> r(n_d_);
        for (int d_i=0; d_i<n_d_; ++d_i) {
          r(d_i) = delta*tmp(d_i) - 0.5*L;
          if (p_i % 2)
            r(d_i) += 0.5*delta;
        }
        for (int b_i=0; b_i<n_bead_; ++b_i) {
          bead_(p_i,b_i)->SetR(r);
          bead_(p_i,b_i)->StoreR();
        }
      }

    // Cubic Lattice
    } else if (init_type == "Cubic") {
      double cofactor = in.GetAttribute<double>("cofactor",1.);
      int n_part__per_n_d = (int) ceil (pow(0.5*n_part_, 1.0/n_d_));
      double delta = cofactor*L/(1.*n_part__per_n_d);
      for (uint32_t p_i=0; p_i<n_part_; p_i++) {
        vec<int> tmp(n_d_);
        tmp(0) = p_i/(n_part__per_n_d*n_part__per_n_d);
        if (n_d_ > 1)
          tmp(1) = (p_i-(tmp(0)*n_part__per_n_d*n_part__per_n_d))/n_part__per_n_d;
        if (n_d_ > 2)
          tmp(2) = p_i - tmp(0)*n_part__per_n_d*n_part__per_n_d - tmp(1)*n_part__per_n_d;
        vec<double> r(n_d_);
        for (int d_i=0; d_i<n_d_; ++d_i)
          r(d_i) = delta*tmp(d_i) - 0.5*L;
        for (int b_i=0; b_i<n_bead_; ++b_i) {
          bead_(p_i,b_i)->SetR(r);
          bead_(p_i,b_i)->StoreR();
        }
      }

    // Line
    } else if (init_type == "Line") {
      double dr = in.GetAttribute<double>("spacing",L/n_part_);
      vec<double> r(zeros<vec<double>>(n_d_));
      for (uint32_t p_i=0; p_i<n_part_; p_i++) {
        for (int b_i=0; b_i<n_bead_; ++b_i) {
          bead_(p_i,b_i)->SetR(r);
          bead_(p_i,b_i)->StoreR();
        }
        r(0) += dr;
      }

    // Restart from previous run (fixme: This is broken now)
    } else if (init_type == "Restart") {
      std::string prefix = in.GetAttribute<std::string>("prefix");
      std::string path_dump_name = in.GetAttribute<std::string>("path_dump_name","path_dump");
      bool parallel = in.GetAttribute<bool>("parallel",false);
      std::stringstream tmp_ss;
      if (parallel)
        tmp_ss << prefix << "." << proc_i << ".h5";
      else
        tmp_ss << prefix << "." << 0 << ".h5";
      std::string file_name = tmp_ss.str();
      std::cout << "Restarting " << name_ << " from " << file_name << "..." << std::endl;
      IO restart_file;
      restart_file.Load(file_name);
      out.Write("System/Particles/"+name_+"/file_name",file_name);

      // Get number of dumps
      uint32_t n_dump;
      restart_file.Read("Observables/"+path_dump_name+"/"+name_+"/n_dump",n_dump);

      // Get positions
      cube<double> positions(n_d_,n_part_*n_bead_,n_dump);
      restart_file.Read("Observables/"+path_dump_name+"/"+name_+"/positions",positions);
      for (uint32_t p_i=0; p_i<n_part_; ++p_i) {
        for (uint32_t b_i=0; b_i<n_bead_; ++b_i) {
          vec<double> r(n_d_);
          for (uint32_t d_i=0; d_i<n_d_; ++d_i)
            r(d_i) = positions(d_i,p_i*n_bead_ + b_i,n_dump-1);
          bead_(p_i,b_i)->SetR(r);
          bead_(p_i,b_i)->StoreR();
        }
      }

      // Get permutation
      cube<double> permutation(2,n_part_,n_dump);
      restart_file.Read("Observables/"+path_dump_name+"/"+name_+"/permutation",permutation);
      for (uint32_t p_i=0; p_i<n_part_; ++p_i) {
        bead_(p_i,0)->SetPrevBead(bead_(permutation(0,p_i,n_dump-1),n_bead_-1));
        bead_(p_i,0)->StorePrev();
        bead_(p_i,n_bead_-1)->SetNextBead(bead_(permutation(1,p_i,n_dump-1),0));
        bead_(p_i,n_bead_-1)->StoreNext();
      }

    }

  }

  void InitRhoK()
  {
    // Resize rho_k
    rho_k_.set_size(n_bead_);
    rho_k_c_.set_size(n_bead_);

    // Initialize rho_k
    for (uint32_t b_i=0; b_i<n_bead_; b_i++)
      for (uint32_t p_i=0; p_i<n_part_; ++p_i)
        bead_(p_i,b_i)->InitRhoK();

    // Sum up rho_k for each time slice
    for (uint32_t b_i=0; b_i<n_bead_; b_i++) {
      rho_k_(b_i).zeros(ks_.vecs.size());
      for (uint32_t p_i=0; p_i<n_part_; p_i++)
        rho_k_(b_i) += bead_(p_i,b_i)->GetRhoK();
    }

    // Copy over
    for (uint32_t b_i=0; b_i<n_bead_; b_i++) {
      StoreRhoK(b_i);
      for (uint32_t p_i=0; p_i<n_part_; p_i++)
        bead_(p_i,b_i)->StoreRhoK();
    }
  }

  /// Update rho_k_p for specific particles belonging to species s_i between time slices b0 and b1
  void UpdateRhoK(const uint32_t b0, const uint32_t b1, const std::vector<uint32_t> &particles, const uint32_t level)
  {
    // Make sure we're in the NEW_MODE
    assert(mode_==NEW_MODE);

    // Reset to old copy
    uint32_t skip = 1<<level;
    for (uint32_t b_i=b0; b_i<b1; b_i+=skip)
      rho_k_(bead_loop(b_i)) = rho_k_c_(bead_loop(b_i));

    // Update C values of changed particles
    for (const auto& p_i: particles) {
      for (uint32_t b_i=b0; b_i<b1; b_i+=skip) {
        // Calculate new values
        GetBead(p_i,b_i)->CalcRhoK();
        // Add in new values and subtract out old values
        rho_k_(bead_loop(b_i)) += GetBead(p_i,b_i)->GetRhoK(NEW_MODE) - GetBead(p_i,b_i)->GetRhoK(OLD_MODE);
      }
    }
    need_update_rho_k_ = false;
  }

  /// Get rho_k or rho_k_c depending on current mode
  field<vec<std::complex<double>>>& GetRhoK() { return mode_ ? rho_k_ : rho_k_c_; };

  /// Set rho_k_c to rho_k for a specific time slice and species
  void StoreRhoK(const uint32_t b_i) { rho_k_c_(bead_loop(b_i)) = rho_k_(bead_loop(b_i)); };

  /// Set rho_k to rho_k_c for a specific time slice and species
  void RestoreRhoK(const uint32_t b_i) { rho_k_(bead_loop(b_i)) = rho_k_c_(bead_loop(b_i)); };

  /// Print out all beads in species
  void Print()
  {
    std::cout << s_i_ << std::endl;
    for (uint32_t p_i=0; p_i<n_part_; ++p_i) {
      for (uint32_t b_i=0; b_i<n_bead_; ++b_i)
        bead_(p_i,b_i)->Print();
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

};

#endif // SIMPIMC_SPECIES_CLASS_H_
