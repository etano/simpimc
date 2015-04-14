#ifndef SIMPIMC_MOVES_MOVE_CLASS_H_
#define SIMPIMC_MOVES_MOVE_CLASS_H_

#include "../event_class.h"
#include "../actions/action_class.h"

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
  void GenerateActionList(std::vector<std::shared_ptr<Action>> &t_action_list, const std::vector<std::string> &species);

  /// Initializes the move
  virtual void Init(Input &in) = 0;

  /// Rejects the move
  virtual void Reject() = 0;

  /// Resets the relevant counters
  virtual void Reset();
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
  virtual void DoEvent();

  /// Writes relevant information about the move to the output file
  virtual void Write();
};

#endif // SIMPIMC_MOVES_MOVE_CLASS_H_
