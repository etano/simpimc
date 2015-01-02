#include "MoveClass.h"

void Move::DoEvent()
{
  struct timeval time;
  gettimeofday(&time, NULL); // Start Time
  double start = time.tv_sec + (time.tv_usec / 1000000.);

  // Attempt move
  nAttempt++;
  if(Attempt()) {
    nAccept++;
    Accept();
  } else
    Reject();

  gettimeofday(&time, NULL); //END-TIME
  double end = time.tv_sec + (time.tv_usec / 1000000.);
  timeSpent += end - start;
}

void Move::GenerateActionList(const std::vector<std::string> &species)
{
  for (auto& action: fullActionList) {
    for (auto& sA: species) {
      if (std::find(action->speciesList.begin(), action->speciesList.end(), sA) != action->speciesList.end()) {
        actionList.push_back(action);
        break;
      }
    }
  }
}

void Move::Reset()
{
  nAccept = 0;
  nAttempt = 0;
}

void Move::Write()
{
  double acceptRatio = (1.*nAccept)/(1.*nAttempt);
  if (firstTime) {
    firstTime = 0;
    out.CreateExtendableDataSet("/Moves/"+name+"/", "nAttempt", nAttempt);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "nAccept", nAccept);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "x", acceptRatio);
  } else {
    out.AppendDataSet("/Moves/"+name+"/", "nAttempt", nAttempt);
    out.AppendDataSet("/Moves/"+name+"/", "nAccept", nAccept);
    out.AppendDataSet("/Moves/"+name+"/", "x", acceptRatio);
  }

  Reset();
}
