#include "move_class.h"

void Move::DoEvent()
{
  struct timeval time;
  gettimeofday(&time, NULL); // Start Time
  double start = time.tv_sec + (time.tv_usec / 1000000.);

  // Attempt move
  n_attempt++;
  if(Attempt()) {
    n_accept++;
    Accept();
  } else
    Reject();

  gettimeofday(&time, NULL); //END-TIME
  double end = time.tv_sec + (time.tv_usec / 1000000.);
  time_spent += end - start;
}

void Move::GenerateActionList(const std::vector<std::string> &species)
{
  for (auto& action: full_action_list) {
    for (auto& sA: species) {
      if (std::find(action->species_list.begin(), action->species_list.end(), sA) != action->species_list.end()) {
        action_list.push_back(action);
        break;
      }
    }
  }
}

void Move::Reset()
{
  n_accept = 0;
  n_attempt = 0;
}

void Move::Write()
{
  double accept_ratio = (1.*n_accept)/(1.*n_attempt);
  if (first_time) {
    first_time = 0;
    out.CreateExtendableDataSet("/Moves/"+name+"/", "n_attempt", n_attempt);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "n_accept", n_accept);
    out.CreateExtendableDataSet("/Moves/"+name+"/", "x", accept_ratio);
  } else {
    out.AppendDataSet("/Moves/"+name+"/", "n_attempt", n_attempt);
    out.AppendDataSet("/Moves/"+name+"/", "n_accept", n_accept);
    out.AppendDataSet("/Moves/"+name+"/", "x", accept_ratio);
  }

  Reset();
}
