#include "NodalClass.h"

double Nodal::DActionDBeta()
{
  return 0.;
}

double Nodal::GetAction(const int b0, const int b1, const vector< pair<int,int> > &particles, const int level)
{
  // Currently old node should be fine
  if (path.mode == 0)
    return 0.;

  // Decide whether or not to check the node
  if (level > maxLevel || !path.speciesList[iSpecies]->fermi)
    return 0.;
  bool checkNode = false;
  for (auto& p: particles) {
    if (p.first == iSpecies) {
      checkNode = true;
      break;
    }
  }
  if (!checkNode)
    return 0.;

  // Constants
  int skip = 1<<level;
  double levelTau = skip*path.tau;

  // See if ref slice included
  bool checkAll = false;
  if (b1 < path.nBead)
    checkAll = ((b0 <= path.refBead) && (b1 >= path.refBead));
  else
    checkAll = (path.beadLoop(b1) >= path.refBead);

  // Set start and end
  if (checkAll) {
    startB = 0;
    endB = path.nBead-1;
  } else {
    startB = b0;
    endB = b1;
  }

  // Set ref beads and initial beads
  vector< std::shared_ptr<Bead> > refBeads, otherBeads;
  int sliceDiff0 = path.beadLoop(startB) - path.refBead;
  int absSliceDiff0 = abs(sliceDiff0);
  for (int iP=0; iP<nPart; ++iP) {
    refBeads.push_back(path(iSpecies,iP,path.refBead));
    if (absSliceDiff0 >= 0) {
      otherBeads.push_back(path.GetNextBead(refBeads[iP],absSliceDiff0)); // fixme: This may be the only correct form
      //otherBeads.push_back(path(iSpecies,iP,startB));
    } else {
      otherBeads.push_back(path.GetPrevBead(refBeads[iP],absSliceDiff0));
      //otherBeads.push_back(path(iSpecies,iP,startB));
    }
  }

  // Compute action
  vec<double> dr(path.nD);
  mat<double> g(nPart,nPart);
  double tot = 0.;
  for (int iB=startB; iB<=endB; iB+=skip) {
    if (iB != path.refBead) {
      // Form rho_F
      int sliceDiff = path.beadLoop(iB) - path.refBead;
      int absSliceDiff = abs(sliceDiff);
      int invSliceDiff = path.nBead - absSliceDiff;
      int minSliceDiff = min(absSliceDiff, invSliceDiff);

      for (int iP=0; iP<nPart; ++iP) {
        for (int jP=0; jP<nPart; ++jP) {
          path.Dr(refBeads[iP], otherBeads[jP], dr);
          g(iP,jP) = GetGij(dr, minSliceDiff);
        }
      }
      rho_F(iB) = det(g);

      // Check sign
      if (rho_F(iB) < 0.) {
        return 1.e100;
      } else {
        tot += 0.;
      }
    }

    // Move down the line
    for (int iP=0; iP<nPart; ++iP)
      otherBeads[iP] = path.GetNextBead(otherBeads[iP],skip);
  }

  return tot;
}

void Nodal::Accept() {
  //int nCheck = nPart;
  //for (int iB=startB; iB<=endB; ++iB)
  //  rho_F_c(iB) = rho_F(iB);
}
