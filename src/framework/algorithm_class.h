#ifndef SIMPIMC_ALGORITHM_CLASS_H_
#define SIMPIMC_ALGORITHM_CLASS_H_

#include "../actions/actions.h"
#include "../events/events.h"

/// Class the actually holds all the events and data objects. This includes the path object which is the main container of all relevant information to the path.
class Algorithm
{
private:
  Loop main_loop; ///< Main loop
  std::vector<std::shared_ptr<Event>> events; /// < Vector of events
  std::vector<std::shared_ptr<Action>> actions; ///< Vector of actions
  Path path; ///< Main datastructure

  /// Initializes all events and objects
  void Init(Input &in, IO &out, RNG &rng, const uint32_t proc_i);
public:
  /// Constructor calls Init.
  Algorithm(Input &in, IO &out, RNG &rng, const uint32_t proc_i)
    : path(proc_i, in, out, rng)
  {
    // Initialize Actions
    out.CreateGroup("Actions");
    for (auto& action_input: in.GetChild("Actions").GetChildList("Action")) {
      std::string type = action_input.GetAttribute<std::string>("type");
      if (type == "Kinetic") {
        actions.push_back(std::make_shared<Kinetic>(path,action_input,out));
      } else if (type == "HarmonicTrap") {
        actions.push_back(std::make_shared<Trap>(path,action_input,out));
      } else if (type == "FreeNodal") {
        actions.push_back(std::make_shared<FreeNodal>(path,action_input,out));
      } else if (type == "OptimizedFreeNodal") {
        actions.push_back(std::make_shared<OptimizedFreeNodal>(path,action_input,out));
      } else if (type == "OptimizedSHONodal") {
        actions.push_back(std::make_shared<OptimizedSHONodal>(path,action_input,out));
      } else if (type == "BarePairAction") {
        actions.push_back(std::make_shared<BarePairAction>(path,action_input,out));
      } else if (type == "DavidPairAction") {
        actions.push_back(std::make_shared<DavidPairAction>(path,action_input,out));
      } else if (type == "IlkkaPairAction") {
        actions.push_back(std::make_shared<IlkkaPairAction>(path,action_input,out));
      } else
        std::cerr << "Warning: Unrecognized Action, " << type << std::endl;
    }

    // Initialize Moves
    out.CreateGroup("Moves");
    for (auto& move_input: in.GetChild("Moves").GetChildList("Move")) {
      std::string type = move_input.GetAttribute<std::string>("type");
      if (type == "Bisect")
        events.push_back(std::make_shared<Bisect>(path,rng,actions,move_input,out));
      else if (type == "DisplaceParticle")
        events.push_back(std::make_shared<DisplaceParticle>(path,rng,actions,move_input,out));
      else if (type == "PermBisectTable")
        events.push_back(std::make_shared<PermBisectTable>(path,rng,actions,move_input,out));
      else if (type == "PermBisectIterative")
        events.push_back(std::make_shared<PermBisectIterative>(path,rng,actions,move_input,out));
      else if (type == "ShiftRefSlice")
        events.push_back(std::make_shared<ShiftRefSlice>(path,rng,actions,move_input,out));
      else if (type == "VaryAction")
        events.push_back(std::make_shared<VaryAction>(path,rng,actions,move_input,out));
      else
        std::cerr << "Warning: Unrecognized Move, " << type << std::endl;
    }

    // Initialize Observables
    out.CreateGroup("Observables");
    std::vector<Input> observable_inputs = in.GetChild("Observables").GetChildList("Observable");
    for (auto& observable_input: observable_inputs) {
      std::string type = observable_input.GetAttribute<std::string>("type");
      if (type == "ContactDensity")
        events.push_back(std::make_shared<ContactDensity>(path,actions,observable_input,out));
      else if (type == "Energy")
        events.push_back(std::make_shared<Energy>(path,actions,observable_input,out));
      else if (type == "ImportanceWeight")
        events.push_back(std::make_shared<ImportanceWeight>(path,actions,observable_input,out));
      else if (type == "PairCorrelation")
        events.push_back(std::make_shared<PairCorrelation>(path,observable_input,out));
      else if (type == "PathDump")
        events.push_back(std::make_shared<PathDump>(path,observable_input,out));
      else if (type == "Permutation")
        events.push_back(std::make_shared<Permutation>(path,observable_input,out));
      else if (type == "RecordOptimizedAction")
        events.push_back(std::make_shared<RecordOptimizedAction>(path,actions,observable_input,out));
      else if (type == "Sign")
        events.push_back(std::make_shared<Sign>(path,observable_input,out));
      else if (type == "StructureFactor")
        events.push_back(std::make_shared<StructureFactor>(path,observable_input,out));
      else if (type == "Time")
        events.push_back(std::make_shared<Time>(path,events,observable_input,out));
      else
        std::cerr << "WARNING: Unrecognized observable, " << type << std::endl;
    }

    // Initialize Write
    events.push_back(std::make_shared<Writes>(out,events,actions,proc_i));

    // Initialize Algorithm
    Input algorithm_input = in.GetChild("Algorithm");
    main_loop.ReadLoop(algorithm_input,events);
  }

  /// Runs the algorithm by calling the DoEvent of the main_loop
  void Run() { main_loop.DoEvent(); };
};

#endif // SIMPIMC_ALGORITHM_CLASS_H_
