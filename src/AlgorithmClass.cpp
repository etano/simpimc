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
    if (type == "Kinetic")
      actions.push_back(new Kinetic(path,actionInput,out));
    else if (type == "HarmonicTrap")
      actions.push_back(new Trap(path,actionInput,out));
    else if (type == "Nodal")
      actions.push_back(new Nodal(path,rng,actionInput,out));
    else if (type == "BarePairAction")
      actions.push_back(new BarePairAction(path,actionInput,out));
    else if (type == "DavidPairAction")
      actions.push_back(new DavidPairAction(path,actionInput,out));
    else if (type == "IlkkaPairAction")
      actions.push_back(new IlkkaPairAction(path,actionInput,out));
    else
      std::cerr << "Warning: Unrecognized Action, " << type << endl;
  }

  // Initialize Moves
  out.CreateGroup("Moves");
  vector<Input> moveInputs = in.getChild("Moves").getChildList("Move");
  for (auto& moveInput: moveInputs) {
    string type = moveInput.getAttribute<string>("type");
    if (type == "Bisect")
      events.push_back(new Bisect(path,rng,actions,moveInput,out));
    else if (type == "DisplaceParticle")
      events.push_back(new DisplaceParticle(path,rng,actions,moveInput,out));
    else if (type == "PermBisect")
      events.push_back(new PermBisect(path,rng,actions,moveInput,out));
    else if (type == "PermBisectIterative")
      events.push_back(new PermBisectIterative(path,rng,actions,moveInput,out));
    else if (type == "ShiftRefSlice")
      events.push_back(new ShiftRefSlice(path,rng,actions,moveInput,out));
    else
      std::cerr << "Warning: Unrecognized Move, " << type << endl;
  }

  // Initialize Observables
  out.CreateGroup("Observables");
  vector<Input> observableInputs = in.getChild("Observables").getChildList("Observable");
  for (auto& observableInput: observableInputs) {
    string type = observableInput.getAttribute<string>("type");
    if (type == "Energy")
      events.push_back(new Energy(path,actions,observableInput,out));
    else if (type == "PairCorrelation")
      events.push_back(new PairCorrelation(path,observableInput,out));
    else if (type == "PathDump")
      events.push_back(new PathDump(path,observableInput,out));
    else if (type == "Permutation")
      events.push_back(new Permutation(path,observableInput,out));
    else if (type == "Sign")
      events.push_back(new Sign(path,observableInput,out));
    else if (type == "StructureFactor")
      events.push_back(new StructureFactor(path,observableInput,out));
    else if (type == "Time")
      events.push_back(new Time(path,events,observableInput,out));
    else
      std::cerr << "WARNING: Unrecognized observable, " << type << endl;
  }

  // Initialize Write
  events.push_back(new Writes(out,events,InterComm));

  // Initialize Algorithm
  vector<Input> loopList = in.getChild("Algorithm").getChildList("Loop");
  for (auto& loop: loopList)
    mainLoop.Init(loop,events);
}


void Algorithm::Run()
{
  mainLoop.DoEvent();
}
