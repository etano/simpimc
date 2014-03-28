#include "ActionClass.h"

void Action::GetOffset(string species, int &iSpecies, int &offset)
{
  int tmpOffset = 0;
  for (unsigned int iS=0; iS<path.nSpecies; iS+=1) {
    if (path.speciesList[iS]->name == species) {
      iSpecies = iS;
      offset = tmpOffset;
    }
    tmpOffset += path.speciesList[iS]->nPart;
  }
}

