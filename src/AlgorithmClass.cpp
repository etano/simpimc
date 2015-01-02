#include "AlgorithmClass.h"

void Algorithm::Init(Input &in, IOClass &out, RNG &rng)
{

  // Initialize Path
  path.Init(in, out, rng);

  // Initialize Actions
  out.CreateGroup("Actions");
  vector<Input> actionInputs = in.getChild("Actions").getChildList("Action");
  for (auto& actionInput: actionInputs) {
    string type = actionInput.getAttribute<string>("type");
    std::shared_ptr<Action> action;
    if (type == "Kinetic") {
      action.reset(new Kinetic(path,actionInput,out));
    } else if (type == "HarmonicTrap") {
      action.reset(new Trap(path,actionInput,out));
    } else if (type == "Nodal") {
      action.reset(new Nodal(path,rng,actionInput,out));
    } else if (type == "BarePairAction") {
      action.reset(new BarePairAction(path,actionInput,out));
    } else if (type == "DavidPairAction") {
      action.reset(new DavidPairAction(path,actionInput,out));
    } else if (type == "IlkkaPairAction") {
      action.reset(new IlkkaPairAction(path,actionInput,out));
    } else
      std::cerr << "Warning: Unrecognized Action, " << type << endl;
    actions.push_back(action);
  }

  // Initialize Moves
  out.CreateGroup("Moves");
  vector<Input> moveInputs = in.getChild("Moves").getChildList("Move");
  for (auto& moveInput: moveInputs) {
    string type = moveInput.getAttribute<string>("type");
    std::shared_ptr<Event> event;
    if (type == "Bisect")
      event.reset(new Bisect(path,rng,actions,moveInput,out));
    else if (type == "DisplaceParticle")
      event.reset(new DisplaceParticle(path,rng,actions,moveInput,out));
    else if (type == "PermBisect")
      event.reset(new PermBisect(path,rng,actions,moveInput,out));
    else if (type == "PermBisectIterative")
      event.reset(new PermBisectIterative(path,rng,actions,moveInput,out));
    else if (type == "ShiftRefSlice")
      event.reset(new ShiftRefSlice(path,rng,actions,moveInput,out));
    else
      std::cerr << "Warning: Unrecognized Move, " << type << endl;
    events.push_back(event);
  }

  // Initialize Observables
  out.CreateGroup("Observables");
  vector<Input> observableInputs = in.getChild("Observables").getChildList("Observable");
  for (auto& observableInput: observableInputs) {
    string type = observableInput.getAttribute<string>("type");
    std::shared_ptr<Event> event;
    if (type == "Energy")
      event.reset(new Energy(path,actions,observableInput,out));
    else if (type == "PairCorrelation")
      event.reset(new PairCorrelation(path,observableInput,out));
    else if (type == "PathDump")
      event.reset(new PathDump(path,observableInput,out));
    else if (type == "Permutation")
      event.reset(new Permutation(path,observableInput,out));
    else if (type == "Sign")
      event.reset(new Sign(path,observableInput,out));
    else if (type == "StructureFactor")
      event.reset(new StructureFactor(path,observableInput,out));
    else if (type == "Time")
      event.reset(new Time(path,events,observableInput,out));
    else
      std::cerr << "WARNING: Unrecognized observable, " << type << endl;
    events.push_back(event);
  }

  // Initialize Write
  std::shared_ptr<Event> write(new Writes(out,events,InterComm));
  events.push_back(write);

  // Initialize Algorithm
  vector<Input> loopList = in.getChild("Algorithm").getChildList("Loop");
  for (auto& loop: loopList)
    mainLoop.Init(loop,events);
}


void Algorithm::Run()
{
  mainLoop.DoEvent();
}
