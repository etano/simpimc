#include "main.h"

int main(int argc, char** argv)
{
  COMM::Init(argc, argv);

  // Get input file
  std::string in_file = "";
  if( argc == 2 )
    in_file = argv[1];
  else {
    std::cout << "Usage: simpimc InputFile.xml\n";
    return 1;
  }

  Simulation sim;
  sim.SetupSimulation(in_file);
  sim.Run();

  COMM::Finalize();

  return 0;
}

