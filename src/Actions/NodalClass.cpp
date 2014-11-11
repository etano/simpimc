#include "NodalClass.h"

void Nodal::Init(Input &in)
{
  // Read in things
  nImages = in.getAttribute<int>("nImages");
  species = in.getAttribute<string>("species");
  cout << "Setting up nodal action for " << species << "..." << endl;
  maxLevel = in.getAttribute<int>("maxLevel",0);
  GetOffset(species,iSpecies,offset);
  nPart = path.speciesList[iSpecies]->nPart;

  // Write things to file
  out.Write("Actions/"+name+"/nImages", nImages);
  out.Write("Actions/"+name+"/species", species);
  out.Write("Actions/"+name+"/maxLevel", maxLevel);

  // Setup splines
  SetupSpline();

  // Set elements of g matrix
  g.set_size(path.nBead);
  g(0).zeros(nPart,nPart);
  for (int iB=1; iB<path.nBead; iB++) {
    g(iB).set_size(nPart,nPart);
    int sliceDiff = path.beadLoop(iB) - path.refBead;
    int absSliceDiff = abs(sliceDiff);
    int invSliceDiff = path.nBead-absSliceDiff;
    int minSliceDiff = min(absSliceDiff, invSliceDiff);
    for (int iP=0; iP<nPart; ++iP) {
      for (int jP=0; jP<nPart; ++jP) {
        Tvector dr(path.nD);
        path.Dr(path(offset+iP,path.refBead), path(offset+iP,iB), dr);
        g(iB)(iP,jP) = GetGij(dr, minSliceDiff);
      }
    }
  }
  g_c = g;

  // Set elements of c matrix (g inverse transpose)
  c.set_size(path.nBead);
  c_c.set_size(path.nBead);
  c_good.set_size(path.nBead);
  for (int iB=0; iB<path.nBead; ++iB) {
    c(iB).set_size(nPart,nPart);
    c_c(iB).set_size(nPart,nPart);
    c_good(iB) = 0;
    //SetCij(iB);
  }

  // Set up determinants
  rho_F.set_size(path.nBead);
  rho_F_c.set_size(path.nBead);
}

// Set elements of c (g inverse transpose
void Nodal::SetCij(int iB)
{
  if (!c_good(iB) && iB!=path.refBead) {
    c_good(iB) = inv(c_c(iB),g_c(iB));
    c_c(iB) = c_c(iB).t();
  } else
    c_c(iB) = c(iB);
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

    // Make rho_0
    RealType rho0 = 0.;
    for (int image=-nImages; image<=nImages; image++)
      rho0 += path.fexp(-path.L*path.L*t_i4LambdaTau);
    RealType logRho0 = -log(rho0);

    // Make rho_free
    for (int i=0; i<r_grid.num; ++i) {
      RealType r = r_grid.start + i*dr;
      rho_free_r(i) = 0.;
      for (int image=-nImages; image<=nImages; image++) {
        RealType t_r = r + image*path.L;
        rho_free_r(i) += path.fexp(-t_r*t_r*t_i4LambdaTau);
      }
      rho_free_r(i) = -log(rho_free_r(i)) - logRho0;
      if (iSpline == nSpline-1)
        cout << iSpline << " " << r << " " << rho_free_r(i) << endl;
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
  for (int iP=0; iP<nPart; ++iP) {
    refBeads.push_back(path.bead(iP+offset,path.refBead));
    if (absSliceDiff0 >= 0) {
      //otherBeads.push_back(path.GetNextBead(refBeads[iP],absSliceDiff0)); // fixme: This may be the only correct form
      otherBeads.push_back(path(iP+offset,startB));
    } else {
      //otherBeads.push_back(path.GetPrevBead(refBeads[iP],absSliceDiff0));
      otherBeads.push_back(path(iP+offset,startB));
    }
  }

  // Compute action
  Tvector dr(path.nD);
  RealType tot = 0.;
  for (int iB=startB; iB<=endB; iB+=skip) {
    if (iB != path.refBead) {
      // Form rho_F
      int sliceDiff = path.beadLoop(iB) - path.refBead;
      int absSliceDiff = abs(sliceDiff);
      int invSliceDiff = path.nBead-absSliceDiff;
      int minSliceDiff = min(absSliceDiff, invSliceDiff);

      if (0 && c_good(iB) && particles.size()==1) {
        RealType gc = 1.;
        c(iB) = c_c(iB);
        g(iB) = g_c(iB);
        for (int p=0; p<particles.size(); ++p) {
          int jP = particles[p] - offset;

          // Update g
          for (int iP=0; iP<nPart; ++iP) {
            path.Dr(refBeads[iP], otherBeads[jP], dr);
            g(iB)(iP,jP) = GetGij(dr, minSliceDiff);
          }

          // Set factor multiplying new determinant
          gc *= dot(g(iB).col(jP),c(iB).col(jP));

          // Set b vector
          Tvector b = (c(iB).t())*g(iB).col(jP);

          // Update c matrix
          Tmatrix t_c = c(iB) - (c(iB).col(jP)*b.t())/b(jP);
          t_c.col(jP) += c(iB).col(jP)/b(jP);
          c(iB) = t_c;

          if (startB == 0) { // fixme: ref bead may not be 0
            // Update g
            for (int iP=0; iP<nPart; ++iP) {
              path.Dr(refBeads[jP], otherBeads[iP], dr);
              g(iB)(jP,iP) = GetGij(dr, minSliceDiff);
            }

            // Set factor multiplying new determinant
            gc *= dot(g(iB).row(jP),c(iB).row(jP));

            // Set b vector
            b = g(iB).row(jP)*c(iB).t();
            b = b.t();

            // Update c matrix
            t_c = c(iB) - (b*c(iB).row(jP))/b(jP);
            t_c.row(jP) += c(iB).row(jP)/b(jP);
            c(iB) = t_c;
          }
        }
        rho_F(iB) = rho_F_c(iB)*gc;
      } else {
        c_good(iB) = 0;
        for (int iP=0; iP<nPart; ++iP) {
          for (int lP=0; lP<nPart; ++lP) {
            path.Dr(refBeads[iP], otherBeads[lP], dr);
            g(iB)(iP,lP) = GetGij(dr, minSliceDiff);
          }
        }
        rho_F(iB) = det(g(iB));
      }

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

RealType Nodal::GetGij(Tvector& r, int sliceDiff)
{
  RealType gaussProd = 1.;
  for (int iD=0; iD<path.nD; iD++) {
    RealType gaussSum;
    eval_UBspline_1d_d(rho_free_r_splines[sliceDiff-1],r(iD),&gaussSum);
    gaussProd *= exp(-gaussSum);
  }
  return gaussProd;
}

void Nodal::Accept() {
  //for (int iB=startB; iB<=endB; ++iB) {
  //  rho_F_c(iB) = rho_F(iB);
  //  g_c(iB) = g(iB);
  //  SetCij(iB);
  //}
}

void Nodal::Write() {}
