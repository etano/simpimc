#ifndef SIMPIMC_SIMULATION_CLASS_H_
#define SIMPIMC_SIMULATION_CLASS_H_

#include "../actions/actions.h"
#include "../events/events.h"

/// Class that actually holds all the events and data objects. This includes the path object which is the main container of all relevant information to the path.
class Simulation
{
private:
  Loop main_loop; ///< Main loop, contains algorithm
  std::vector<std::shared_ptr<Event>> events; /// < Vector of events
  std::vector<std::shared_ptr<Action>> actions; ///< Vector of actions
  Path path; ///< Main data structure
public:
  /// Constructor initializes everything
  Simulation(Input &in, IO &out, RNG &rng, const uint32_t proc_i)
    : path(proc_i, in, out, rng)
  {
    // Initialize Actions
    out.CreateGroup("Actions");
    for (auto& input: in.GetChild("Actions").GetChildList("Action"))
      actions.push_back(ActionFactory(input,out,path));

    // Initialize Moves
    out.CreateGroup("Moves");
    for (auto& input: in.GetChild("Moves").GetChildList("Move"))
      events.push_back(MoveFactory(input,out,path,rng,actions));

    // Initialize Observables
    out.CreateGroup("Observables");
    for (auto& input: in.GetChild("Observables").GetChildList("Observable"))
      events.push_back(ObservableFactory(input,out,path,actions,events));

    // Initialize Write
    events.push_back(std::make_shared<Writes>(out,events,actions,proc_i));

    // Initialize Algorithm
    Input algorithm_input = in.GetChild("Algorithm");
    main_loop.ReadLoop(algorithm_input,events);
  }

  /// Runs the algorithm by calling the DoEvent of the main_loop
  void Run()
  {
    std::cout << "Running main loop..." << std::endl;
    main_loop.DoEvent();
  };
};

#endif // SIMPIMC_SIMULATION_CLASS_H_
