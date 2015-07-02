#include "framework/framework_class.h"

int main(int argc, char** argv)
{
  COMM::Init(argc, argv);

  // Get input file
  std::string in_file = "";
  if( argc == 2 )
    in_file = argv[1];
  else {
    std::cout << "Usage: simpimc input_file.xml\n";
    return 1;
  }

  Framework fw(in_file);
  fw.Run();

  COMM::Finalize();

  return 0;
}

