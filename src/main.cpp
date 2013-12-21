#include "main.h"

using namespace std;

int main(int argc, char** argv)
{
  COMM::Init(argc, argv);

  // Get input file
  string inFile = "";
  if( argc == 2 )
    inFile = argv[1];
  else {
    cout << "Usage: ./test InputFile\n";
    return 1;
  }

  Simulation sim;
  sim.SetupIO(inFile);
  sim.BuildMPIModel();
  sim.Run();

  COMM::Finalize();

  return 0;
}

