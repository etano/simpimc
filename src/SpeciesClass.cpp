#include "SpeciesClass.h"

void Species::Init(Input &in, IOClass &out)
{
  // Read inputs
  name = in.getAttribute<string>("name");
  nPart = in.getAttribute<int>("nPart");
  lambda = in.getAttribute<double>("lambda");
  fermi = in.getAttribute<int>("fermi", 0);
  fixedNode = in.getAttribute<int>("fixedNode", 0);
  initType = in.getAttribute<string>("initType","Random");

  // Write to file
  out.CreateGroup("System/Particles/"+name);
  out.Write("System/Particles/"+name+"/nPart",nPart);
  out.Write("System/Particles/"+name+"/lambda",lambda);
  out.Write("System/Particles/"+name+"/fermi",fermi);
  out.Write("System/Particles/"+name+"/fixedNode",fixedNode);
  out.Write("System/Particles/"+name+"/initType",initType);

  // Set defaults
  needUpdateRhoK = true;

  // Initiate beads
  bead.set_size(nPart,nBead);
  unsigned int unique_id = 0;
  for (unsigned int iP=0; iP<nPart; iP++) {
    for (unsigned int iB=0; iB<nBead; iB++) {
      bead(iP,iB) = std::make_shared<Bead>(nD,iS,iP,iB,unique_id);
      unique_id++;
    }
  }

  // Initiate bead connections
  for (unsigned int iP = 0; iP < nPart; iP++) {
    bead(iP,0) -> next = bead(iP,1);
    bead(iP,0) -> nextC = bead(iP,1);
    bead(iP,0) -> prev = bead(iP,nBead-1);
    bead(iP,0) -> prevC = bead(iP,nBead-1);
    for (unsigned int iB = 1; iB < nBead-1; iB++) {
      bead(iP,iB) -> next = bead(iP,iB+1);
      bead(iP,iB) -> nextC = bead(iP,iB+1);
      bead(iP,iB) -> prev = bead(iP,iB-1);
      bead(iP,iB) -> prevC = bead(iP,iB-1);
    }
    bead(iP,nBead-1) -> next = bead(iP,0);
    bead(iP,nBead-1) -> nextC = bead(iP,0);
    bead(iP,nBead-1) -> prev = bead(iP,nBead-2);
    bead(iP,nBead-1) -> prevC = bead(iP,nBead-2);
  }
}

void Species::InitPaths(Input &in, IOClass &out, RNG &rng, CommunicatorClass& InterComm, int L)
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
              initFile >> bead(iP,iB)->r(iD);
            bead(iP,iB)->storeR();
          }
        } else {
          vec<double> r(nD);
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
    double cofactor = in.getAttribute<double>("cofactor",1.);
    for (int iP=0; iP<nPart; ++iP) {
      for (int iD=0; iD<nD; ++iD) {
        double tmpRand = rng.unifRand(-cofactor*L/2.,cofactor*L/2.);
        for (int iB=0; iB<nBead; ++iB)
          bead(iP,iB)->r(iD) = tmpRand;
      }
      for (int iB=0; iB<nBead; ++iB)
        bead(iP,iB)->storeR();
    }

  // BCC Lattice
  } else if (initType == "BCC") {
    int nPartPerND = (int) ceil (pow(0.5*nPart, 1.0/nD));
    double delta = L/nPartPerND;
    for (unsigned int iP=0; iP<nPart; iP++) {
      int p = iP/2;
      vec<int> tmp(nD);
      tmp(0) = p/(nPartPerND*nPartPerND);
      if (nD > 1)
        tmp(1) = (p-(tmp(0)*nPartPerND*nPartPerND))/nPartPerND;
      if (nD > 2)
        tmp(2) = p - tmp(0)*nPartPerND*nPartPerND - tmp(1)*nPartPerND;
      vec<double> r(nD);
      for (int iD=0; iD<nD; ++iD) {
        r(iD) = delta*tmp(iD) - 0.5*L;
        if (iP % 2)
          r(iD) += 0.5*delta;
      }
      for (int iB=0; iB<nBead; ++iB) {
        bead(iP,iB)->r = r;
        bead(iP,iB)->storeR();
      }
    }

  // Cubic Lattice
  } else if (initType == "Cubic") {
    double cofactor = in.getAttribute<double>("cofactor",1.);
    int nPartPerND = (int) ceil (pow(0.5*nPart, 1.0/nD));
    double delta = cofactor*L/(1.*nPartPerND);
    for (unsigned int iP=0; iP<nPart; iP++) {
      vec<int> tmp(nD);
      tmp(0) = iP/(nPartPerND*nPartPerND);
      if (nD > 1)
        tmp(1) = (iP-(tmp(0)*nPartPerND*nPartPerND))/nPartPerND;
      if (nD > 2)
        tmp(2) = iP - tmp(0)*nPartPerND*nPartPerND - tmp(1)*nPartPerND;
      vec<double> r(nD);
      for (int iD=0; iD<nD; ++iD)
        r(iD) = delta*tmp(iD) - 0.5*L;
      for (int iB=0; iB<nBead; ++iB) {
        bead(iP,iB)->r = r;
        bead(iP,iB)->storeR();
      }
    }

  // Restart from previous run (fixme: This is broken now)
  } else if (initType == "Restart") {
    string prefix = in.getAttribute<string>("prefix");
    int parallel = in.getAttribute<bool>("parallel",0);
    stringstream tmpSS;
    if (parallel)
      tmpSS << prefix << "." << InterComm.MyProc() << ".h5";
    else
      tmpSS << prefix << "." << 0 << ".h5";
    string fileName = tmpSS.str();
    cout << "Restarting " << name << " from " << fileName << "..." << endl;
    IOClass restartFile;
    restartFile.load(fileName);
    out.Write("System/Particles/"+name+"/fileName",fileName);

    // Get number of dumps
    int nDump;
    restartFile.Read("Observables/PathDump/"+name+"/nDump",nDump);

    // Get positions
    cube<double> pathPositions(nD,nPart*nBead,nDump);
    restartFile.Read("Observables/PathDump/"+name+"/positions",pathPositions);
    for (int iP=0; iP<nPart; ++iP) {
      for (int iB=0; iB<nBead; ++iB) {
        for (int iD=0; iD<nD; ++iD) {
          bead(iP,iB)->r(iD) = pathPositions(iD,iP*nBead + iB,nDump-1);
        }
        bead(iP,iB)->storeR();
      }
    }

    // Get permutation
    cube<double> pathPermutation(2,nPart,nDump);
    restartFile.Read("Observables/PathDump/"+name+"/permutation",pathPermutation);
    for (int iP=0; iP<nPart; ++iP) {
      bead(iP,0)->prev = bead(pathPermutation(0,iP,nDump-1),nBead-1);
      bead(iP,0)->storePrev();
      bead(iP,nBead-1)->next = bead(pathPermutation(1,iP,nDump-1),0);
      bead(iP,nBead-1)->storeNext();
    }

  }

}
