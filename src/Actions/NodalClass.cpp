#include "NodalClass.h"

void Nodal::Init(Input &in)
{
  // Read in things
  nImages = in.getAttribute<int>("nImages");
  species = in.getAttribute<string>("species");
  cout << "Setting up nodal action for " << species << "..." << endl;
  maxLevel = in.getAttribute<int>("maxLevel",0);
  GetOffset(species,iSpecies,offset);

  // Write things to file
  out.Write("Actions/"+name+"/nImages", nImages);
  out.Write("Actions/"+name+"/species", species);
  out.Write("Actions/"+name+"/maxLevel", maxLevel);
}

RealType Nodal::DActionDBeta()
{
  return 0.;
}

RealType Nodal::GetAction(int b0, int b1, vector<int> &particles, int level)
{
  // Currently old node should be fine
  if (path.mode == 0)
    return 0.;

  // Decide whether or not to check the node
  if (level > maxLevel || !path.speciesList[iSpecies]->fermi)
    return 0.;
  bool checkNode = false;
  for (int p=0; p<particles.size(); ++p) {
    if (path(particles[p],b0)->s == iSpecies) {
      checkNode = true;
      break;
    }
  }
  if (!checkNode)
    return 0.;

  // Constants
  int skip = 1<<level;
  RealType levelTau = skip*path.tau;
  RealType i4LambdaTau = 1./(4.*path.speciesList[iSpecies]->lambda*path.tau);

  // See if ref slice included
  bool checkAll = false;
  if (b1 < path.nBead)
    checkAll = ((b0 <= path.refBead) && (b1 >= path.refBead));
  else
    checkAll = (path.beadLoop(b1) >= path.refBead);
  if (checkAll) {
    b0 = 0;
    b1 = path.nBead;
  }

  // Compute action
  Tmatrix rho_F(path.speciesList[iSpecies]->nPart,path.speciesList[iSpecies]->nPart);
  Tvector dr(path.nD);
  RealType tot = 0.;
  RealType gaussProd, gaussSum, dist;
  for (int iB=b0; iB<b1; iB+=skip) {
    if (iB != path.refBead) {
      // Form rho_F
      int sliceDiff = min(abs(path.beadLoop(iB) - path.refBead), abs(path.refBead - path.beadLoop(iB)));
      for (int iP=offset; iP<offset+path.speciesList[iSpecies]->nPart; ++iP) {
        for (int jP=offset; jP<offset+path.speciesList[iSpecies]->nPart; ++jP) {
          //path.Dr(path.bead(iP,path.refBead), path.bead(jP,path.refBead)->nextB(sliceDiff), dr);
          path.Dr(path.bead(iP,path.refBead), path.bead(jP,iB), dr);
          gaussProd = 1.;
          for (int iD=0; iD<path.nD; iD++) {
            gaussSum = 0.;
            for (int image=-nImages; image<=nImages; image++) {
              dist = dr(iD) + image*path.L;
              gaussSum += path.fexp(-dist*dist*i4LambdaTau/sliceDiff);
            }
            gaussProd *= gaussSum;
          }
          rho_F(iP-offset,jP-offset) = gaussProd;
        }
      }
  
      // Take determinant and check sign
      RealType det_rho_F = det(rho_F);
      if (det_rho_F < 0.) {
        return 1.e100;
      } else {
        tot += 0.;
      }
    }
  }

  return tot;
}

void Nodal::Write() {}
