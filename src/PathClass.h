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
  unsigned int nPart, nD, nBead;
  double beta, L, iL, vol;

  // Calculated Global Constants
  double tau;
  unsigned int maxLevel;

  // Species
  unsigned int nSpecies;
  vector<Species*> speciesList;
  void GetSpeciesInfo(string species, int &iSpecies);

  // Fast math
  bool approximate;
  inline double fexp(double x) { return exp(x); };

  // Mode (use copy or true)
  bool mode;
  void SetMode(int m) { mode = m; };
  bool GetMode() { return mode; };

  // Beads
  vector<Bead*>::const_iterator beadIter;
  vec<int> beadLoop;
  Bead* operator() (int iS, int iP, int iB) { return speciesList[iS]->bead(iP,beadLoop(iB)); };
  void storeR(vector<Bead*> &affBeads);
  void restoreR(vector<Bead*> &affBeads);
  inline vec<double>& GetR(Bead* b) { return mode ? b->r : b->rC; };
  inline Bead* GetNextBead(Bead* b, int n) { return mode ? b->nextB(n) : b->nextBC(n); };
  inline Bead* GetPrevBead(Bead* b, int n) { return mode ? b->prevB(n) : b->prevBC(n); };
  inline void Dr(vec<double> &r0, vec<double> &r1, vec<double> &dr) { dr = r0 - r1; PutInBox(dr); };
  inline void Dr(Bead* b0, vec<double> &r1, vec<double> &dr) { Dr(GetR(b0), r1, dr); };
  inline void Dr(Bead* b0, Bead* b1, vec<double> &dr) { Dr(GetR(b0), GetR(b1), dr); };
  inline void RBar(Bead* b0, Bead* b1, vec<double> &rBar) { Dr(b0, b1, rBar); rBar = GetR(b1) + 0.5*rBar; };
  void PrintPath();

  // Periodic Boundary Condition
  bool PBC;
  void PutInBox(vec<double>& r);

  // k vectors and rho_k
  vector< vec<double> > ks;
  vector< vec<int> > kIndices;
  vector<double> magKs;
  field< vec< complex<double> > > C;
  field< vec< complex<double> > > rhoK, rhoKC;
  vec<double> kBox;
  double kC;
  vec<int> maxKIndex;
  bool Include(vec<double> &k, double kCut);
  void SetupKs(double kCut);
  void InitRhoK();
  void UpdateRhoK();
  void UpdateRhoK(int b0, int b1, vector< pair<int,int> > &particles, int level);
  void UpdateRhoKP(int b0, int b1, vector< pair<int,int> > &particles, int level);
  void UpdateRhoKP(int b0, int b1, int iS, vector<int> &particles, int level);
  void CalcC(vec<double> &r);
  void AddRhoKP(field< vec< complex<double> > >& tmpRhoK, int iP, int iB, int iS, int pm);
  inline void CalcRhoKP(Bead* b);
  inline field< vec< complex<double> > >& GetRhoK() { return mode ? (rhoK) : (rhoKC); };
  inline vec< complex<double> >& GetRhoK(Bead* b) { return mode ? (b->rhoK) : (b->rhoKC); };
  inline void storeRhoK(int iB, int iS) { rhoKC(beadLoop(iB),iS) = rhoK(beadLoop(iB),iS); };
  inline void restoreRhoK(int iB, int iS) { rhoK(beadLoop(iB),iS) = rhoKC(beadLoop(iB),iS); };
  void storeRhoKP(vector<Bead*>& affBeads);
  void restoreRhoKP(vector<Bead*>& affBeads);

  // Nodes
  int sign;
  int refBead;
  int CalcSign();

  // Permutations
  struct CompareVecInt
  {
    bool operator() (const vector<int> &a, const vector<int> &b) {
      for (int i = 0; i<a.size(); i++)
        if (a[i] != b[i])
          return (a[i] > b[i]);
      return (a[0]>b[0]);
    }
  };

  map<vector<int>,int,CompareVecInt> possPerms;
  map<vector<int>,int,CompareVecInt>::const_iterator possPermsIterator;
  void SetCycleCount(int iS, vector<int>& cycles);
  int GetPermSector(int iS);
  int GetPermSector(int iS, vector<int>& cycles);
  void SetupPermSectors(int n, int sectorsMax);


  // Path initialization
  void InitPaths(Input &in, IOClass &out, RNG &rng);

  // Get dr, drP, and drrP
  inline void DrDrPDrrP(const int b0, const int b1, const int s0, const int s1, const int p0, const int p1, double &rMag, double &rPMag, double &rrPMag, vec<double> &r, vec<double> &rP, vec<double> &rrP)
  {
    //Dr((*this)(s1,p1,b0),(*this)(s0,p0,b0),r);
    //Dr((*this)(s1,p1,b1),(*this)(s0,p0,b1),rP);
  
    r = GetR((*this)(s1,p1,b0)) - GetR((*this)(s0,p0,b0));
    rP = GetR((*this)(s1,p1,b1)) - GetR((*this)(s0,p0,b1));
    for (int iD=0; iD<nD; ++iD) {
      r(iD) -= nearbyint(r(iD)*iL)*L;
      rP(iD) += nearbyint((r(iD)-rP(iD))*iL)*L;
      rrP(iD) = r(iD) - rP(iD);
      rrP(iD) -= nearbyint(rrP(iD)*iL)*L;
    }
    rMag = mag(r);
    rPMag = mag(rP);
    rrPMag = mag(rrP);
  
  };
};


#endif
