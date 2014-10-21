#include "ActionClass.h"

void Action::GetOffset(string species, int &iSpecies, int &offset)
{
  int tmpOffset = 0;
  for (unsigned int iS=0; iS<path.nSpecies; iS+=1) {
    if (path.speciesList[iS]->name == species) {
      iSpecies = iS;
      offset = tmpOffset;
      break;
    }
    tmpOffset += path.speciesList[iS]->nPart;
  }

  if (tmpOffset == path.nPart) {
    // Species doesn't exist
    cout << "Warning: species " << species << " does not exist elsewhere..." << endl;
    iSpecies = -1;
    offset = path.nPart;
  }
}

