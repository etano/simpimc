#ifndef SIMPIMC_MOVES_MOVE_CLASS_H_
#define SIMPIMC_MOVES_MOVE_CLASS_H_

#include "../event_class.h"
#include "../actions/action_class.h"

class Move : public Event
{
private:

protected:
  Path &path;
  RNG &rng;
  IO &out;

  std::vector<std::shared_ptr<Action>> &full_action_list;
  std::vector<std::shared_ptr<Action>> action_list;
  void GenerateActionList(const std::vector<std::string> &species);
  std::string prefix;
public:
  // Constructor
  Move(Path &tmp_path, RNG &tmp_rng, std::vector<std::shared_ptr<Action>> &tmp_action_list, Input &in, IO &tmp_out)
    : Event(), path(tmp_path), rng(tmp_rng), full_action_list(tmp_action_list), out(tmp_out)
  {
    name = in.GetAttribute<std::string>("name");
    type = in.GetAttribute<std::string>("type");
    prefix = "Moves/"+name;
    out.CreateGroup(prefix);
    out.Write(prefix+"/type",type);
    std::string data_type = "scalar";
    out.Write(prefix+"/data_type",data_type);
    first_time = 1;
    Reset();
  }

  std::string type;
  virtual void DoEvent();

  // Moves
  virtual void Init(Input &in) {};
  virtual bool Attempt() {};
  virtual void Accept() {};
  virtual void Reject() {};

  // Acceptance
  bool first_time;
  uint n_attempt, n_accept;
  virtual void Reset();

  // Write
  virtual void Write();

};

#endif // SIMPIMC_MOVES_MOVE_CLASS_H_
