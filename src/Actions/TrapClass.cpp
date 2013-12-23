#include "TrapClass.h"

void Trap::Init(Input &in)
{
  omega = in.getAttribute<RealType>("omega");
}

RealType Trap::DActionDBeta()
{
  RealType tot = 0.;

  for (unsigned int iP=0; iP<path.nPart; iP+=1) {
    for (unsigned int iB=0; iB<path.nBead; iB+=1) {
      tot += dot(path(iP,iB)->r, path(iP,iB)->r);
    }
  }

  return 0.5*omega*omega*(1. + 3.*path.tau*path.tau*omega*omega/12.)*tot;
}

RealType Trap::GetAction(int b0, int b1, vector<int> &particles, int level)
{
  int skip = 1<<level;
  RealType levelTau = skip*path.tau;
  RealType tot = 0.;
  for (int iP=0; iP<particles.size(); ++iP) {
    for (int iB=b0; iB<b1; iB+=skip) {
      Tvector dr(path.nD);
      if(path.mode)
        dr = path(iP,iB)->r;
      else
        dr = path(iP,iB)->rC;
      tot += dot(dr, dr);
    }
  }

  return 0.5*levelTau*omega*omega*(1. + path.tau*path.tau*omega*omega/12.)*tot;
}

void Trap::Write()
{

}
