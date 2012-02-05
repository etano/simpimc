#include "simpimc.h" // Import headers, set global constants, and build system

using namespace std;

int main (int argc, char* argv[])
{
  ////////////////////////////////////
  /* Initialize Simulation Settings */
  ////////////////////////////////////

  // Input File
  char* inputFile = argv[1];
  ifstream inputStream;
  inputStream.open(inputFile);
  if (!inputStream) cout << "There was a problem opening the input file " << inputFile << " for reading.\n";
  string inputFileLabel;
  getline(inputStream, inputFileLabel); // Skip first line, get label
  unsigned int nLineSkip;
  if(argv[2] == NULL) nLineSkip = 1;
  else nLineSkip = atoi(argv[2]);
  cout << "\nRunning simulation " << nLineSkip << " from " << inputFile << " ( " << inputFileLabel << " ) ...\n\n";

  // Inputs
  unsigned int nPart, nD, nBead; // number of particles, dimensions, time slices
  double beta, lambda, L; // Beta, lambda, Box Size
  double duration; // duration of simulation
  int nStep, block, blockOut, nEqSweep, nEqStep; // Frequency of Block Output
  bool fermi; // 0 - Boson, 1 - Fermion
  int halfspace; // -1 - Negative Halfspace, 1 - Positive Halfspace
  int nodeType; // 0 - Exact Nodes, 1 - High T Nodes, 2 - Low T Nodes
  bool useNodeDist; // 0 - No, 1 - Yes
  for (unsigned int iLine = 0; iLine < nLineSkip; iLine += 1) {
    inputStream >> nPart >> nD >> nBead >> beta >> lambda >> L;
    inputStream >> duration >> nStep >> block >> blockOut >> nEqSweep >> nEqStep;
    inputStream >> fermi >> halfspace >> nodeType >> useNodeDist;
  }
  inputStream.close();

  double tau = beta/(1.0*nBead); // tau

  // Output Settings to Screen
  cout << scientific << setprecision(4);
  cout << "\nSIMULATION SETTINGS::\n" 
       << "\nN: " << nPart << "\nD: " << nD << "\nM: " << nBead 
       << "\nBeta: " << beta << "\nLambda: " << lambda << "\nL: " << L 
       << "\nDuration (s): " << duration << "\nnStep: " << nStep << "\nBlock Size: " << block << "\nBlock Output: " << blockOut
       << "\nnEqSweep: " << nEqSweep << "\nnEqStep: " << nEqStep
       << "\nFermions?(1/0): " << fermi << "\nHalfspace(1/0): " << halfspace << "\nNode Type(1/0): " << nodeType << "\nUse Node Distance(1/0): " << useNodeDist
       << endl;

  ///////////////////////////
  /* Initialize Simulation */
  ///////////////////////////

  if(!fermi) useNodeDist = 0;
  Simulation sim(nPart,nD,nBead,beta,lambda,fermi,halfspace,nodeType,useNodeDist,L);

  //////////////////////////
  /* Initialization Moves */
  //////////////////////////

  // ( path , rng ,  perAcceptDesired , nEqSweeps , nEqSteps , moveSkip )
  sim.moves.push_back(new Bisect(sim.path,sim.rng,0.5,nEqSweep,nEqStep,1));
  //sim.moves.push_back(new PermBisect(sim.path,sim.rng,0.5,nEqSweeps,nEqSteps,1));
  //sim.moves.push_back(new DisplaceBead(sim.path,sim.rng,0.5,10,1000,1));
  //sim.moves.push_back(new DisplaceParticle(sim.path,sim.rng,0.5,10,1000,1));
  //sim.moves.push_back(new DisplaceAll(sim.path,sim.rng,0.5,10,1000,1));
  //sim.moves.push_back(new Relabel(sim.path,sim.rng,0.5,10,1000,1));
  //sim.moves.push_back(new SimplePerm(sim.path,sim.rng,0.5,10,1000,1));

  // Equilibrate Moves
  for (vector<Move*>::const_iterator move = sim.moves.begin(); move != sim.moves.end(); ++move) {
    (*move) -> Equilibrate();
  }

  ////////////////////////////
  /* Initialize Observables */
  ////////////////////////////

  // Form Output String
  stringstream outputSuffix;
  outputSuffix << "-" << nPart << "-" << nD << "-" << nBead << "-" << beta << "-" << lambda << "-" << L
               << "-" << duration << "-" << nStep << "-" << block << "-" << blockOut << "-" << nEqSweep << "-" << nEqStep
               << "-" << fermi << "-" << halfspace << "-" << nodeType << "-" << useNodeDist;

  // Permutation type
  int pType;

  // Observables
  // ( path , outputSuffix , outputLabel , skip , block )
  sim.observables.push_back( new Energy(sim.path,outputSuffix.str(),"Energy",1,block) );
  sim.observables.push_back( new R(sim.path,outputSuffix.str(),"R",1,block) );
  sim.observables.push_back( new R2(sim.path,outputSuffix.str(),"R2",1,block) );

  ////////////////////////
  /* Main Loop Settings */
  ////////////////////////

  // Set up time loop
  time_t start, end;
  time (&start);
  time (&end);
  double timeDif = difftime (end,start);

  ///////////////////////////
  /* Main Monte Carlo Loop */
  ///////////////////////////

  cout << "\nRUN SIMULATION::\n\n";

  // Monte Carlo Steps
  int iStep = 1;
  while (timeDif < duration || iStep < nStep) {

    // Make Move
    for (vector<Move*>::const_iterator move = sim.moves.begin(); move != sim.moves.end(); ++move) {
      if (!(iStep % (*move) -> moveSkip)) (*move) -> MakeMove();
    }

    // Make Measurements
    for (vector<Observable*>::const_iterator observable = sim.observables.begin(); observable != sim.observables.end(); ++observable) {
      if (!(iStep % (*observable) -> skip)) {
        pType = sim.path.getPType(); // Get the current permutation configuration
        (*observable) -> Accumulate(pType); // Accumulate Measurements
        if (!(iStep % (*observable) -> block)) {
          (*observable) -> Output(); // Output Block
          (*observable) -> nBlock++; // Increment Block
        }
      }
    }

    // Print block information to screen
    if (!(iStep % blockOut)) {
      time (&end);
      timeDif = difftime (end,start);
      cout << endl << "T- " << duration - timeDif << " s" << endl;
      cout << "#: " << iStep << endl;
      for (vector<Observable*>::const_iterator observable = sim.observables.begin(); observable != sim.observables.end(); ++observable) {
        cout << endl << (*observable) -> observableLabel << " : Block # " << (*observable) -> nBlock << endl;
        (*observable) -> Print();
      }
    }

    // Increment
    time (&end);
    timeDif = difftime (end,start);
    iStep += 1;
  }

  ////////////////////
  /* Output results */
  ////////////////////

  // General Results
  cout << "\n\nFINAL RESULTS::\n\n";
  cout << "Inf: " << sim.path.infCount << " , NaN: " << sim.path.nanCount << " , Err: " << sim.path.errCount << "\n\n";
  cout << "# Monte Carlo sweeps: " << iStep << endl;

  // Move Type Results
  cout << "\nPer Move Type: \n";  
  for (vector<Move*>::const_iterator move = sim.moves.begin(); move != sim.moves.end(); ++move) {
    cout << (*move) -> moveLabel << ": " << (*move) -> getPerAccept() << endl;
  }

  // Per Permutation Configuration
  cout << "\nPer Perm Type: \n";  
  unsigned int nPerm;
  if (fermi) nPerm = 3; // Only 3 particle exchanges
  else nPerm = 6; // Only 2 & 3 particle exchanges

  for (unsigned int iPerm = 0; iPerm < nPerm; iPerm += 1)  {
    cout << iPerm << " : " << sim.path.permCount(iPerm,1) << " " << sim.path.permCount(iPerm,0) << " " << (sim.path.permCount(iPerm,1)*1.0)/(sim.path.permCount(iPerm,0)*1.0) << endl;
  }

  // Compute statistics of measure quantities
  for (vector<Observable*>::const_iterator observable = sim.observables.begin(); observable != sim.observables.end(); ++observable) {
    (*observable) -> Stats();
  }

  cout << endl;
  sim.path.printPerm();
  //sim.path.printBeads();

  return 0;
}
