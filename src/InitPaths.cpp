#include "PathClass.h"

void Path::InitPaths(Input &in, IOClass &out, RNG &rng)
{
  out.CreateGroup("Init");
  string initType = in.getChild("Init").getAttribute<string>("type");
  out.Write("Init/type",initType);

  // Read configuration from xyz file
  if (initType == "File") {
    string fileName = in.getChild("Init").getAttribute<string>("name");
    int allBeads = in.getChild("Init").getAttribute<bool>("allBeads",false);
    out.Write("Init/fileName",fileName);
    out.Write("Init/allBeads",allBeads);
    fstream initFile;
    initFile.open(fileName.c_str(), std::ios_base::in);
    if (!initFile.fail()) {
      for (int iP=0; iP<nPart; ++iP) {
        if (allBeads) {
          for (int iB=0; iB<nBead; ++iB) {
            initFile >> iP >> iB;
            for (int iD=0; iD<nD; ++iD)
              initFile >> bead(iP,iB)->r(iD);
            bead(iP,iB)->storeR();
          }
        } else {
          initFile >> iP;
          Tvector r(nD);
          for (int iD=0; iD<nD; ++iD)
            initFile >> r(iD);
          for (int iB=0; iB<nBead; ++iB) {
            bead(iP,iB)->r = r;
            bead(iP,iB)->storeR();
          }
        }
      }
      initFile.close();
    } else {
      cout << "ERROR: Init file '" << fileName << "' does not exist." << endl;
      exit(1);
    }

  // Random configuration
  } else if (initType == "Random") {
    for (int iP=0; iP<nPart; ++iP) {
      RealType tmpRand = rng.unifRand();
      for (int iB=0; iB<nBead; ++iB) {
        for (int iD=0; iD<nD; ++iD)
          bead(iP,iB)->r(iD) = tmpRand;
        bead(iP,iB)->storeR();
      }
    }

  // BCC Lattice
  } else if (initType == "BCC") {
    for (int iS=0; iS<nSpecies; iS+=1) {
      int offset;
      GetSpeciesInfo(speciesList[iS]->name,iS,offset);
      int nPartPerND = (int) ceil (pow(0.5*speciesList[iS]->nPart, 1.0/nD)-1.0e-6);
      RealType delta = L/nPartPerND;
      for (unsigned int iP=offset; iP<offset+speciesList[iS]->nPart; iP++) {
        int p = (iP-offset)/2;
        Tvector r(nD);
        r(0) = p/(nPartPerND*nPartPerND);
        if (nD > 1)
          r(1) = (p-(r(0)*nPartPerND*nPartPerND))/nPartPerND;
        if (nD > 2)
          r(2) = p - r(0)*nPartPerND*nPartPerND - r(1)*nPartPerND;
        r *= delta;
        r -= 0.5*L;
        if (iP % 2)
          r += 0.5*delta;
        for (int iB=0; iB<nBead; ++iB) {
          bead(iP,iB)->r = r;
          bead(iP,iB)->storeR();
        }
      }
    }

  // Restart from previous run
  } else if (initType == "Restart") {
    string prefix = in.getChild("Init").getAttribute<string>("prefix");
    int parallel = in.getChild("Init").getAttribute<bool>("parallel",0);
    stringstream tmpSS;
    if (parallel)
      tmpSS << prefix << "." << InterComm.MyProc() << ".h5";
    else
      tmpSS << prefix << "." << 0 << ".h5";
    string fileName = tmpSS.str();
    cout << "Restarting from " << fileName << "..." << endl;
    IOClass restartFile;
    restartFile.load(fileName);
    out.Write("Init/fileName",fileName);

    // Get number of dumps
    int nDump;
    restartFile.Read("Observables/PathDump/nDump",nDump);

    // Get positions
    Tcube pathPositions(nPart*nBead,nD,nDump);
    restartFile.Read("Observables/PathDump/positions",pathPositions);
    for (int iP=0; iP<nPart; ++iP) {
      for (int iB=0; iB<nBead; ++iB) {
        for (int iD=0; iD<nD; ++iD)
          bead(iP,iB)->r(iD) = pathPositions(iP*nBead + iB,iD,nDump-1);
        bead(iP,iB)->storeR();
      }
    }

    // Get permutation
    Tcube pathPermutation(nPart,2,nDump);
    restartFile.Read("Observables/PathDump/permutation",pathPermutation);
    for (int iP=0; iP<nPart; ++iP) {
      bead(iP,0)->prev = bead(pathPermutation(iP,0,nDump-1),nBead-1);
      bead(iP,0)->storePrev();
      bead(iP,nBead-1)->next = bead(pathPermutation(iP,1,nDump-1),0);
      bead(iP,nBead-1)->storeNext();
    }
  }

}
