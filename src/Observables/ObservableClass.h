#ifndef ObservableClass_H
#define ObservableClass_H

#include "../StandardLibraries.h"       // Standard libraries
#include "../GlobalConstants.h"
#include "../PathClass.h"
#include "../BeadClass.h"
#include "../Stats.h" // statistics

class Observable
{
private: 

protected:
  // Path
  Path& path;

  // Label
  std::string outputSuffix;

  // Trace
  std::fstream trace;

  // Constants
  double oneOverNbeadBlock;

public:
  // Constructor
  Observable( Path& pathIn , std::string outputSuffixIn , std::string observableLabelIn , unsigned int skipIn , unsigned int blockIn )
    : path(pathIn) , outputSuffix(outputSuffixIn) , observableLabel(observableLabelIn) , skip(skipIn) , block(blockIn) , nBlock(0)
  {
    outputFile = "data/traces/" + observableLabel + outputSuffix;
    std::cout << "\nOutputting " << observableLabel << " data to " << outputFile << "." << endl;
    trace.open (outputFile.c_str(), ios::out | ios::trunc);

    oneOverNbeadBlock = skip*1.0/(path.nBead * block * 1.0);
    Output();
    Print();
  }

  std::string observableLabel;
  std::string outputFile;
  unsigned int skip;
  unsigned int block;
  unsigned int nBlock;

  // Functions
  virtual void Accumulate( const int iType ) {};
  virtual void Output() {};
  virtual void Print() {};
  virtual void Stats() {};
};

#endif
