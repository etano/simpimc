#include "EventClass.h"

void Event::WriteTime()
{
  double norm = 1./nIterations;
  timeSpent *= norm;

  if (FirstTimeTime) {
    FirstTimeTime = 0;
    out.CreateExtendableDataSet("/TimeAnalysis/", name, timeSpent);
  } else {
    out.AppendDataSet("/TimeAnalysis/", name, timeSpent);
  }

  timeSpent = 0.;
}
