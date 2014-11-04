#ifndef PathClass_H
#define PathClass_H

#include "SpeciesClass.h"
#include "BeadClass.h"
#include "Utils/config.h"
#include "Utils/Communication/Communication.h"
#include "Utils/IO/InputClass.h"
#include "Utils/IO/IOClass.h"
#include "Utils/RNG/RNGClass.h"
#include "Utils/Algorithm/Algorithm.h"

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
  RealType beta, L, iL, vol;

  // Calculated Global Constants
  RealType tau;
  unsigned int maxLevel;

  // Species
  unsigned int nSpecies;
  vector<Species*> speciesList;
  void GetSpeciesInfo(string species, int &iSpecies, int &offset);

  // Fast math
  bool approximate;
  inline RealType fexp(RealType x) { return exp(x); };

  // Mode (use copy or true)
  bool mode;
  void SetMode(int m) { mode = m; };
  bool GetMode() { return mode; };

  // Beads
  field<Bead*> bead;
  vector<Bead*>::const_iterator beadIter;
  Ivector beadLoop;
  Bead* operator() (int iP, int iB) { return bead(iP,beadLoop(iB)); };
  void storeR(vector<Bead*> &affBeads);
  void restoreR(vector<Bead*> &affBeads);
  inline Tvector& GetR(Bead* b) { return mode ? b->r : b->rC; };
  inline Bead* GetNextBead(Bead* b, int n) { return mode ? b->nextB(n) : b->nextBC(n); };
  inline Bead* GetPrevBead(Bead* b, int n) { return mode ? b->prevB(n) : b->prevBC(n); };
  inline void Dr(Tvector &r0, Tvector &r1, Tvector &dr) { dr = r0 - r1; PutInBox(dr); };
  inline void Dr(Bead* b0, Tvector &r1, Tvector &dr) { Dr(GetR(b0), r1, dr); };
  inline void Dr(Bead* b0, Bead* b1, Tvector &dr) { Dr(GetR(b0), GetR(b1), dr); };
  inline void RBar(Bead* b0, Bead* b1, Tvector &rBar) { Dr(b0, b1, rBar); rBar = GetR(b1) + 0.5*rBar; };
  inline void DrDrPDrrP(const int b0, const int b1, const int p0, const int p1, RealType &rMag, RealType &rPMag, RealType &rrPMag, Tvector &r0, Tvector &r1, Tvector &dr);
  void PrintPath();

  // Periodic Boundary Condition
  bool PBC;
  void PutInBox(Tvector& r);

  // k vectors and rho_k
  vector<Tvector> ks;
  vector<Ivector> kIndices;
  vector<RealType> magKs;
  field<Cvector> C;
  field<Cvector> rhoK, rhoKC;
  Tvector kBox;
  RealType kC;
  Ivector maxKIndex;
  bool Include(Tvector &k, RealType kCut);
  void SetupKs(RealType kCut);
  void InitRhoK();
  void UpdateRhoK();
  void UpdateRhoK(int b0, int b1, vector<int> &particles, int level);
  void UpdateRhoKP(int b0, int b1, vector<int> &particles, int level);
  void CalcC(Tvector &r);
  void AddRhoKP(field<Cvector>& tmpRhoK, int iP, int iB, int iS, int pm);
  inline void CalcRhoKP(Bead* b);
  inline field<Cvector>& GetRhoK() { return mode ? (rhoK) : (rhoKC); };
  inline Cvector& GetRhoK(Bead* b) { return mode ? (b->rhoK) : (b->rhoKC); };
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
};

// Get dr, drP, and drrP
inline void Path::DrDrPDrrP(const int b0, const int b1, const int p0, const int p1, RealType &rMag, RealType &rPMag, RealType &rrPMag, Tvector &r, Tvector &rP, Tvector &rrP)
{
  Dr(bead(p1,beadLoop(b0)),bead(p0,beadLoop(b0)),r);
  Dr(bead(p1,beadLoop(b1)),bead(p0,beadLoop(b1)),rP);
  rMag = mag(r);
  rPMag = mag(rP);
  Dr(r,rP,rrP);
  rrPMag = mag(rrP);
}


#endif
