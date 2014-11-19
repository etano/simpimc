#include "SpeciesClass.h"

void Species::Init(Input &in, IOClass &out)
{
  // Read inputs
  name = in.getAttribute<string>("name");
  nPart = in.getAttribute<int>("nPart");
  lambda = in.getAttribute<RealType>("lambda");
  fermi = in.getAttribute<int>("fermi", 0);
  fixedNode = in.getAttribute<int>("fixedNode", 0);
  initType = in.getAttribute<string>("initType");

  // Write to file
  out.CreateGroup("System/Particles/"+name);
  out.Write("System/Particles/"+name+"/nPart",nPart);
  out.Write("System/Particles/"+name+"/lambda",lambda);
  out.Write("System/Particles/"+name+"/fermi",fermi);
  out.Write("System/Particles/"+name+"/fixedNode",fixedNode);
  out.Write("System/Particles/"+name+"/initType",initType);

  // Set defaults
  needUpdateRhoK = true;
}

void Species::InitPaths(Input &in, IOClass &out, RNG &rng, field<Bead*>& bead, CommunicatorClass& InterComm, int L)
{
  // Read configuration from xyz file
  if (initType == "File") {
    string fileName = in.getAttribute<string>("fileName");
    int allBeads = in.getAttribute<bool>("allBeads",false);
    out.Write("System/Particles/"+name+"/fileName",fileName);
    out.Write("System/Particles/"+name+"/allBeads",allBeads);
    fstream initFile;
    initFile.open(fileName.c_str(), std::ios_base::in);
    if (!initFile.fail()) {
      for (int iP=0; iP<nPart; ++iP) {
        if (allBeads) {
          for (int iB=0; iB<nBead; ++iB) {
            for (int iD=0; iD<nD; ++iD)
              initFile >> bead(iP+offset,iB)->r(iD);
            bead(iP+offset,iB)->storeR();
          }
        } else {
          Tvector r(nD);
          for (int iD=0; iD<nD; ++iD)
            initFile >> r(iD);
          for (int iB=0; iB<nBead; ++iB) {
            bead(iP+offset,iB)->r = r;
            bead(iP+offset,iB)->storeR();
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
    RealType cofactor = in.getAttribute<RealType>("cofactor",1.);
    for (int iP=0; iP<nPart; ++iP) {
      RealType tmpRand = rng.unifRand(-cofactor*L/2.,cofactor*L/2.);
      for (int iB=0; iB<nBead; ++iB) {
        for (int iD=0; iD<nD; ++iD)
          bead(iP+offset,iB)->r(iD) = tmpRand;
        bead(iP+offset,iB)->storeR();
      }
    }

  // BCC Lattice
  } else if (initType == "BCC") {
    int nPartPerND = (int) ceil (pow(0.5*nPart, 1.0/nD)-1.0e-6);
    RealType delta = L/nPartPerND;
    for (unsigned int iP=0; iP<nPart; iP++) {
      int p = iP/2;
      Ivector tmp(nD);
      tmp(0) = p/(nPartPerND*nPartPerND);
      if (nD > 1)
        tmp(1) = (p-(tmp(0)*nPartPerND*nPartPerND))/nPartPerND;
      if (nD > 2)
        tmp(2) = p - tmp(0)*nPartPerND*nPartPerND - tmp(1)*nPartPerND;
      Tvector r(nD);
      for (int iD=0; iD<nD; ++iD)
        r(iD) = 1.*tmp(iD);
      r *= delta;
      r -= 0.5*L;
      if (iP % 2)
        r += 0.5*delta;
      for (int iB=0; iB<nBead; ++iB) {
        bead(iP+offset,iB)->r = r;
        bead(iP+offset,iB)->storeR();
      }
    }

  // Restart from previous run
  } else if (initType == "Restart") {
    string prefix = in.getAttribute<string>("prefix");
    int parallel = in.getAttribute<bool>("parallel",0);
    stringstream tmpSS;
    if (parallel)
      tmpSS << prefix << "." << InterComm.MyProc() << ".h5";
    else
      tmpSS << prefix << "." << 0 << ".h5";
    string fileName = tmpSS.str();
    cout << "Restarting from " << fileName << "..." << endl;
    IOClass restartFile;
    restartFile.load(fileName);
    out.Write("System/Particles/"+name+"/fileName",fileName);

    // Get number of dumps
    int nDump;
    restartFile.Read("Observables/PathDump/"+name+"/nDump",nDump);

    // Get positions
    Tcube pathPositions(nPart*nBead,nD,nDump);
    restartFile.Read("Observables/PathDump/"+name+"/positions",pathPositions);
    for (int iP=0; iP<nPart; ++iP) {
      for (int iB=0; iB<nBead; ++iB) {
        for (int iD=0; iD<nD; ++iD)
          bead(iP+offset,iB)->r(iD) = pathPositions(iP*nBead + iB,iD,nDump-1);
        bead(iP+offset,iB)->storeR();
      }
    }

    // Get permutation
    Tcube pathPermutation(nPart,2,nDump);
    restartFile.Read("Observables/PathDump/"+name+"/permutation",pathPermutation);
    for (int iP=0; iP<nPart; ++iP) {
      bead(iP+offset,0)->prev = bead(pathPermutation(iP,0,nDump-1),nBead-1);
      bead(iP+offset,0)->storePrev();
      bead(iP+offset,nBead-1)->next = bead(pathPermutation(iP,1,nDump-1),0);
      bead(iP+offset,nBead-1)->storeNext();
    }

  }

}
