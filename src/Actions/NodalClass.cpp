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

  // Setup splines
  SetupSpline();
}

// Create spline
void Nodal::SetupSpline()
{
  // Setup grid
  Ugrid r_grid;
  r_grid.start = -path.L/2.;
  r_grid.end = path.L/2.;
  r_grid.num = 1000;
  RealType dr = (r_grid.end - r_grid.start)/(r_grid.num - 1);

  // Resize spline field
  int nSpline = path.nBead/2 + (path.nBead%2) + 1;
  rho_free_r_splines.set_size(nSpline);

  // Create splines
  for (int iSpline=0; iSpline<nSpline; ++iSpline) {
    Tvector rho_free_r(r_grid.num);
    RealType t_i4LambdaTau = 1./(4.*path.speciesList[iSpecies]->lambda*path.tau*(iSpline+1));
    for (int i=0; i<r_grid.num; ++i) {
      RealType r = r_grid.start + i*dr;
      rho_free_r(i) = 0.;
      for (int image=-nImages; image<=nImages; image++) {
        RealType t_r = r + image*path.L;
        rho_free_r(i) += path.fexp(-t_r*t_r*t_i4LambdaTau);
      }
    }
    BCtype_d xBC = {NATURAL, FLAT}; // fixme: Is this correct?
    UBspline_1d_d* rho_free_r_spline = create_UBspline_1d_d(r_grid, xBC, rho_free_r.memptr());
    rho_free_r_splines[iSpline] = rho_free_r_spline;
  }
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

  // See if ref slice included
  bool checkAll = false;
  if (b1 < path.nBead)
    checkAll = ((b0 <= path.refBead) && (b1 >= path.refBead));
  else
    checkAll = (path.beadLoop(b1) >= path.refBead);

  // Set start and end
  int startB, endB;
  if (checkAll) {
    startB = 0;
    endB = path.nBead-1;
  } else {
    startB = b0;
    endB = b1;
  }

  // Set ref beads and initial beads
  vector<Bead*> refBeads, otherBeads;
  int sliceDiff0 = path.beadLoop(startB) - path.refBead;
  int absSliceDiff0 = abs(sliceDiff0);
  for (int iP=0; iP<path.speciesList[iSpecies]->nPart; ++iP) {
    refBeads.push_back(path.bead(iP+offset,path.refBead));
    if (absSliceDiff0 >= 0) {
      //otherBeads.push_back(path.GetNextBead(refBeads[iP],absSliceDiff0)); // fixme: This may be the only correct form
      otherBeads.push_back(path.bead(iP+offset,startB));
    } else {
      //otherBeads.push_back(path.GetPrevBead(refBeads[iP],absSliceDiff0));
      otherBeads.push_back(path.bead(iP+offset,startB));
    }
  }

  // Compute action
  Tmatrix rho_F(path.speciesList[iSpecies]->nPart,path.speciesList[iSpecies]->nPart);
  Tvector dr(path.nD);
  RealType tot = 0.;
  RealType gaussProd, gaussSum, dist;
  for (int iB=startB; iB<=endB; iB+=skip) {
    if (iB != path.refBead) {
      // Form rho_F
      int sliceDiff = path.beadLoop(iB) - path.refBead;
      int absSliceDiff = abs(sliceDiff);
      int invSliceDiff = path.nBead-absSliceDiff;
      int minSliceDiff = min(absSliceDiff, invSliceDiff);
      for (int iP=0; iP<path.speciesList[iSpecies]->nPart; ++iP) {
        for (int jP=0; jP<path.speciesList[iSpecies]->nPart; ++jP) {
          path.Dr(refBeads[iP], otherBeads[jP], dr);
          rho_F(iP,jP) = GetRhoij(dr, minSliceDiff);
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

    // Move down the line
    for (int iP=0; iP<path.speciesList[iSpecies]->nPart; ++iP)
      otherBeads[iP] = path.GetNextBead(otherBeads[iP],skip);
  }

  return tot;
}

RealType Nodal::GetRhoij(Tvector& r, int sliceDiff)
{
  RealType gaussProd = 1.;
  for (int iD=0; iD<path.nD; iD++) {
    RealType gaussSum;
    eval_UBspline_1d_d(rho_free_r_splines[sliceDiff-1],r(iD),&gaussSum);
    gaussProd *= gaussSum;
  }
  return gaussProd;
}

void Nodal::Write() {}
