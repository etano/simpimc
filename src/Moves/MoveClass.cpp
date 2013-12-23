#include "MoveClass.h"

void Move::Init(Input &in)
{
  name = in.getAttribute<string>("name");
}

void Move::Equilibrate()
{
  for (unsigned int iEqSweep = 0; iEqSweep < nEqSweep; iEqSweep += 1) {
    resetCounters();
    for (unsigned int iEqStep = 0; iEqStep < nEqStep; iEqStep += 1)
      MakeMove();
    perAccept = getPerAccept();  // Percentage accepted
    stepSize *= 1.0 - perAcceptDesired + perAccept; // Recalculate step size
    if (stepSize > path.maxLevel) stepSize = path.maxLevel;
  }
  std::cout << name << ": " << stepSize << ", Percent Accepted: " << perAccept << "\n";
  resetCounters();
}

void Move::resetCounters()
{
  nAccept = 0;
  nAttempt = 0;
}

RealType Move::getPerAccept()
{
  perAccept = (nAccept*1.) / (nAttempt*1.);  // Percentage accepted
  return perAccept;
}

// Reassign particle labels
void Move::assignParticleLabels()
{
  Bead *b;

  for (unsigned int iP = 0; iP < path.nPart; iP += 1) {
    b = path.bead(iP,0);
    for (unsigned int iB = 0; iB < path.nBead; iB += 1) {
      path.bead(iP,iB) = b;
      path.bead(iP,iB) -> p = iP;
      b = b -> next;
    }
  }
}

void Move::setPerm( const int permType , int* perm , int* iPerm , const int i , const int j , const int k )
{
  // All 3 particle exchanges (first 1-2 for fermions and bosons, last 3-5 for bosons only, default is identity)
  switch (permType){
    case 1:
      perm[i] = k;
      perm[j] = i;
      perm[k] = j;
      iPerm[i] = j;
      iPerm[j] = k;
      iPerm[k] = i;
      break;
    case 2:
      perm[i] = j;
      perm[j] = k;
      perm[k] = i;
      iPerm[i] = k;
      iPerm[j] = i;
      iPerm[k] = j;
      break;
    case 3:
      perm[i] = j;
      perm[j] = i;
      perm[k] = k;
      iPerm[i] = j;
      iPerm[j] = i;
      iPerm[k] = k;
      break;
    case 4:
      perm[i] = i;
      perm[j] = k;
      perm[k] = j;
      iPerm[i] = i;
      iPerm[j] = k;
      iPerm[k] = j;
      break;
    case 5:
      perm[i] = k;
      perm[j] = j;
      perm[k] = i;
      iPerm[i] = k;
      iPerm[j] = j;
      iPerm[k] = i;
      break;
    default:
      perm[i] = i;
      perm[j] = j;
      perm[k] = k;
      iPerm[i] = i;
      iPerm[j] = j;
      iPerm[k] = k;
      break;
  }
}
