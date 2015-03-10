#ifndef PathClass_H
#define PathClass_H

#include "SpeciesClass.h"
#include "BeadClass.h"
#include "config.h"

class Path
{
private:

protected:

public:
  // Constructor
  Path(CommunicatorClass& tmpWorldComm, CommunicatorClass& tmpInterComm, CommunicatorClass& tmpIntraComm)
   : WorldComm(tmpWorldComm), InterComm(tmpInterComm), IntraComm(tmpIntraComm)
  {}

  void Init(Input &in, IOClass &out, RNG &rng);

  // Parallel communicators
  CommunicatorClass& WorldComm; // This is the global MPIWORLD communicator.
  CommunicatorClass& InterComm; // This is for communication between the rank 0 procs of each walker group.
  CommunicatorClass& IntraComm; // This is for commmunication between procs within a walker group.

  // Given Global Constants
  uint nPart, nD, nBead;
  double beta, L, iL, vol;

  // Calculated Global Constants
  double tau;
  uint maxLevel;

  // Species
  uint nSpecies;
  vector< std::shared_ptr<Species> > speciesList;
  void GetSpeciesInfo(const string& species, uint &iSpecies);

  // Fast math
  bool approximate;
  inline double fexp(const double x) { return exp(x); };

  // Mode (use copy or true)
  bool mode;
  void SetMode(uint m) { mode = m; };
  bool GetMode() { return mode; };

  // Beads
  void PrintPath();
  vec<uint> beadLoop;
  std::shared_ptr<Bead> operator() (const uint iS, const uint iP, const uint iB) { return speciesList[iS]->bead(iP,beadLoop(iB)); };
  void storeR(std::vector< std::shared_ptr<Bead> > &affBeads);
  void restoreR(std::vector< std::shared_ptr<Bead> > &affBeads);
  inline vec<double>& GetR(const std::shared_ptr<Bead> b) { return mode ? b->r : b->rC; };
  inline std::shared_ptr<Bead> GetNextBead(const std::shared_ptr<Bead> b, const uint n) { return mode ? b->nextB(n) : b->nextBC(n); };
  inline std::shared_ptr<Bead> GetPrevBead(const std::shared_ptr<Bead> b, const uint n) { return mode ? b->prevB(n) : b->prevBC(n); };
  inline vec<double> Dr(const vec<double>& r0, const vec<double>& r1) { vec<double> dr = r0-r1; PutInBox(dr); return dr; };
  inline vec<double> Dr(const std::shared_ptr<Bead> b0, const vec<double>& r1) { return Dr(GetR(b0), r1); };
  inline vec<double> Dr(const std::shared_ptr<Bead> b0, const std::shared_ptr<Bead> b1) { return Dr(GetR(b0), GetR(b1)); };
  inline vec<double> RBar(const std::shared_ptr<Bead> b0, const std::shared_ptr<Bead> b1) { return GetR(b1) + 0.5*Dr(b0, b1); };
  template<class T>
  inline double MagDr(T &r0, T &r1) { return mag(Dr(r0,r1)); };

  // Periodic Boundary Condition
  bool PBC;
  void PutInBox(vec<double> &r);

  // k vectors and rho_k
  vector<vec<double> > ks;
  vector<vec<int> > kIndices;
  vector<double> magKs;
  field<vec<complex<double> > > C;
  field<vec<complex<double> > > rhoK, rhoKC;
  vec<double> kBox;
  double kC;
  vec<int> maxKIndex;
  bool Include(const vec<double>& k, const double kCut);
  void SetupKs(const double kCut);
  void InitRhoK();
  void UpdateRhoK();
  void UpdateRhoK(const uint b0, const uint b1, const vector< pair<uint,uint> >& particles, const uint level);
  void UpdateRhoKP(const uint b0, const uint b1, const vector< pair<uint,uint> >& particles, const uint level);
  void UpdateRhoKP(const uint b0, const uint b1, const uint iS, const vector<uint>& particles, const uint level);
  void CalcC(const vec<double>& r);
  void AddRhoKP(field<vec<complex<double> > >& tmpRhoK, const uint iP, const uint iB, const uint iS, const int pm);
  inline void CalcRhoKP(const std::shared_ptr<Bead> b);
  inline field<vec<complex<double> > >& GetRhoK() { return mode ? (rhoK) : (rhoKC); };
  inline vec<complex<double> >& GetRhoK(const std::shared_ptr<Bead> b) { return mode ? (b->rhoK) : (b->rhoKC); };
  inline void storeRhoK(const uint iB, const uint iS) { rhoKC(beadLoop(iB),iS) = rhoK(beadLoop(iB),iS); };
  inline void restoreRhoK(const uint iB, const uint iS) { rhoK(beadLoop(iB),iS) = rhoKC(beadLoop(iB),iS); };
  void storeRhoKP(vector<std::shared_ptr<Bead>>& affBeads);
  void restoreRhoKP(vector<std::shared_ptr<Bead>>& affBeads);

  // Importance weight
  double importance_weight;
  int sign;

  // Nodes
  uint refBead;

  // Permutations
  int CalcSign();
  struct CompareVecInt
  {
    bool operator() (const vector<uint>& a, const vector<uint>& b) {
      for (uint i = 0; i<a.size(); i++)
        if (a[i] != b[i])
          return (a[i] > b[i]);
      return (a[0]>b[0]);
    }
  };

  map<vector<uint>,uint,CompareVecInt> possPerms;
  map<vector<uint>,uint,CompareVecInt>::const_iterator possPermsIterator;
  void SetCycleCount(const uint iS, vector<uint>& cycles);
  uint GetPermSector(const uint iS);
  uint GetPermSector(const uint iS, vector<uint>& cycles);
  void SetupPermSectors(const uint n, const uint sectorsMax);
  bool permSectorsSetup;


  // Path initialization
  void InitPaths(Input &in, IOClass &out, RNG &rng);

  // Get dr, drP, and drrP
  inline void DrDrPDrrP(const uint b0, const uint b1, const uint s0, const uint s1, const uint p0, const uint p1, double &rMag, double &rPMag, double &rrPMag)
  {
    //Dr((*this)(s1,p1,b0),(*this)(s0,p0,b0),r);
    //Dr((*this)(s1,p1,b1),(*this)(s0,p0,b1),rP);

    vec<double> r = GetR((*this)(s1,p1,b0)) - GetR((*this)(s0,p0,b0));
    vec<double> rP = GetR((*this)(s1,p1,b1)) - GetR((*this)(s0,p0,b1));
    for (uint iD=0; iD<nD; ++iD) {
      r(iD) -= nearbyint(r(iD)*iL)*L;
      rP(iD) += nearbyint((r(iD)-rP(iD))*iL)*L;
    }
    vec<double> rrP = r - rP;
    for (uint iD=0; iD<nD; ++iD)
      rrP(iD) -= nearbyint(rrP(iD)*iL)*L;
    rMag = mag(r);
    rPMag = mag(rP);
    rrPMag = mag(rrP);

  };
};


#endif
