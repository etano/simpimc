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
  unsigned int nPart; // number of particles
  unsigned int nD; // dimension
  unsigned int nBead; // number of time slices
  double beta; // inverse temperature
  double duration; // duration of simulation
  bool fermi; // 0 - Boson, 1 - Fermion
  int halfspace; // -1 - Negative Halfspace, 1 - Positive Halfspace
  int nodeType; // 0 - Exact Nodes, 1 - High T Nodes, 2 - Low T Nodes  

  for (unsigned int iLine = 0; iLine < nLineSkip; iLine += 1) {
    inputStream >> nPart;
    inputStream >> nD;
    inputStream >> nBead;
    inputStream >> beta;
    inputStream >> duration;
    inputStream >> fermi;
    inputStream >> halfspace;
    inputStream >> nodeType;   
  }
  inputStream.close();

  double tau = beta/(1.0*nBead); // tau

  // Output Settings to Screen
  cout << scientific << setprecision(4);   
  cout << "\nSIMULATION SETTINGS::\n";
  cout << "\nDuration (s): " << duration << "\nBeta: " << beta << "\nN: " << nPart << "\nD: " << nD << "\nM: " << nBead << "\nDTau: " << tau << "\nOmega: " << omega << "\nFermions?(1/0): " << fermi << "\n";   
  

  ///////////////////////////
  /* Initialize Simulation */
  ///////////////////////////
  
  int useNodeDist = 0;
  if(!fermi) useNodeDist = 0;
  Simulation sim(nPart,nD,nBead,beta,fermi,halfspace,nodeType,useNodeDist);
  
  //////////////////////////
  /* Initialization Moves */
  //////////////////////////

  // ( path , perAcceptDesired , nEqSweeps , nEqSteps , moveSkip )
  sim.moves.push_back(new Bisect(sim.path,sim.rng,0.5,10,1000,1));
  sim.moves.push_back(new PermBisect(sim.path,sim.rng,0.5,10,1000,1));
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
  char outputFormat[] = "-%d-%d-%d-%g-%g-%d-%d-%d.dat";
  char outputFile[sizeof outputFormat];
  sprintf(outputFile,outputFormat,nPart,nD,nBead,beta,duration,fermi,halfspace,nodeType); 
  string outputSuffix(outputFile);
  
  // Blocking
  int block = 1000;  
  int blockOut = 100*block;  
  
  // Permutation type
  int pType;

  // Observables
  // ( path , outputSuffix , outputLabel , skip , block )
  sim.observables.push_back( new Energy(sim.path,outputSuffix,"Energy",1,block) );
  sim.observables.push_back( new R(sim.path,outputSuffix,"R",1,block) );
  sim.observables.push_back( new R2(sim.path,outputSuffix,"R2",1,block) );

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
  int iStep = 0;
  while (timeDif < duration) {  

    // Make Move
    for (vector<Move*>::const_iterator move = sim.moves.begin(); move != sim.moves.end(); ++move) {
      if (!(iStep % (*move) -> moveSkip)) (*move) -> MakeMove();
    }

    // Make Measurements
        
    pType = sim.path.getPType(); // Get the current permutation configuration

    for (vector<Observable*>::const_iterator observable = sim.observables.begin(); observable != sim.observables.end(); ++observable) {
      if (!(iStep % (*observable) -> skip)) {
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
      cout << "\nT- " << duration - timeDif << " s" << "\n";
      for (vector<Observable*>::const_iterator observable = sim.observables.begin(); observable != sim.observables.end(); ++observable) {
        cout << "\n" << (*observable) -> observableLabel << " : Block # " << (*observable) -> nBlock << "\n";
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
  cout << "# Monte Carlo sweeps: " << iStep << "\n";
  
  // Move Type Results
  cout << "\nPer Move Type: \n";  
  for (vector<Move*>::const_iterator move = sim.moves.begin(); move != sim.moves.end(); ++move) {
    cout << (*move) -> moveLabel << ": " << (*move) -> getPerAccept() << "\n";
  }
  
  // Per Permutation Configuration
  cout << "\nPer Perm Type: \n";  
  unsigned int nPerm;
  if (fermi) nPerm = 3; // Only 3 particle exchanges
  else nPerm = 6; // Only 2 & 3 particle exchanges
  
  for (unsigned int iPerm = 0; iPerm < nPerm; iPerm += 1)  {
    cout << iPerm << " : " << sim.path.permCount(iPerm,1) << " " << sim.path.permCount(iPerm,0) << " " << (sim.path.permCount(iPerm,1)*1.0)/(sim.path.permCount(iPerm,0)*1.0) << "\n";
  }     
  
  // Compute statistics of measure quantities
  for (vector<Observable*>::const_iterator observable = sim.observables.begin(); observable != sim.observables.end(); ++observable) {
    (*observable) -> Stats();
  }

  cout << "\n"; 
  sim.path.printPerm();
  //sim.path.printBeads();
  
  return 0;
}
