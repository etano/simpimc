#ifndef SIMPIMC_MOVES_MOVE_CLASS_H_
#define SIMPIMC_MOVES_MOVE_CLASS_H_

#include "../event_class.h"
#include "../../actions/action_class.h"

/// Parent class for all moves
class Move : public Event
{
protected:
  bool first_time; ///< Whether or not the move is called for the first time
  uint32_t n_accept; ///< Number of accepted moves
  uint32_t n_attempt; ///< Number of attempted moves
  Path &path; ///< Reference to path
  RNG &rng; ///< Reference to the random number generator
  IO &out; ///< Reference to the output file
  std::vector<std::shared_ptr<Action>> action_list; ///< Vector of pointers to actions that affect the relevant species
  std::string prefix; ///< Prefix used for all output
  std::string type; ///< Name of the type of move

  /// Accept the move
  virtual void Accept() = 0;

  /// Attempt the move
  virtual bool Attempt() = 0;

  /// Generate a list of actions that affect a given vector of species names
  void GenerateActionList(std::vector<std::shared_ptr<Action>> &t_action_list, const std::vector<std::string> &species)
  {
    for (auto& action: t_action_list)
      for (auto& sA: species)
        if (std::find(action->species_list.begin(), action->species_list.end(), sA) != action->species_list.end()) {
          action_list.push_back(action);
          break; // Only add the action once even if it's found by multiple species
        }
  }

  /// Rejects the move
  virtual void Reject() = 0;

  /// Resets the relevant counters
  virtual void Reset()
  {
    n_accept = 0;
    n_attempt = 0;
  }

public:

  /// Constructor instiates the references to path, rng, and out, as well as create initial output group for the move.
  Move(Path &tmp_path, RNG &tmp_rng, std::vector<std::shared_ptr<Action>> &t_action_list, Input &in, IO &tmp_out)
    : Event(), path(tmp_path), rng(tmp_rng), out(tmp_out)
  {
    // Set move name and type and write them to file
    name = in.GetAttribute<std::string>("name");
    type = in.GetAttribute<std::string>("type");
    prefix = "Moves/"+name;
    out.CreateGroup(prefix);
    out.Write(prefix+"/type",type);
    std::string data_type = "scalar";
    out.Write(prefix+"/data_type",data_type);

    // Reset counters
    first_time = 1;
    Reset();
  }

  /// Attempts the move and determines if it is accepted or rejected
  virtual void DoEvent()
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

  /// Writes relevant information about the move to the output file
  virtual void Write()
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

};

#endif // SIMPIMC_MOVES_MOVE_CLASS_H_
