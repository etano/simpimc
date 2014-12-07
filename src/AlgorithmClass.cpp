#include "AlgorithmClass.h"

void Algorithm::Init(Input &in, IOClass &out, RNG &rng)
{

  // Initialize Path
  path.Init(in, out, rng);

  // Initialize Actions
  out.CreateGroup("Actions");
  vector<Input> actionInputs = in.getChild("Actions").getChildList("Action");
  for (int i=0; i<actionInputs.size(); ++i) {
    string type = actionInputs[i].getAttribute<string>("type");
    if (type == "Kinetic")
      actions.emplace_back(new Kinetic(path,actionInputs[i],out));
    else if (type == "HarmonicTrap")
      actions.emplace_back(new Trap(path,actionInputs[i],out));
    else if (type == "Nodal")
      actions.emplace_back(new Nodal(path,rng,actionInputs[i],out));
    else if (type == "BarePairAction")
      actions.emplace_back(new BarePairAction(path,actionInputs[i],out));
    else if (type == "DavidPairAction")
      actions.emplace_back(new DavidPairAction(path,actionInputs[i],out));
    else if (type == "IlkkaPairAction")
      actions.emplace_back(new IlkkaPairAction(path,actionInputs[i],out));
    else
      std::cerr << "Warning: Unrecognized Action, " << type << endl;
  }

  // Initialize Moves
  out.CreateGroup("Moves");
  vector<Input> moveInputs = in.getChild("Moves").getChildList("Move");
  for (int i=0; i<moveInputs.size(); ++i) {
    string type = moveInputs[i].getAttribute<string>("type");
    if (type == "Bisect")
      events.emplace_back(new Bisect(path,rng,actions,moveInputs[i],out));
    else if (type == "DisplaceParticle")
      events.emplace_back(new DisplaceParticle(path,rng,actions,moveInputs[i],out));
    else if (type == "PermBisect")
      events.emplace_back(new PermBisect(path,rng,actions,moveInputs[i],out));
    else if (type == "PermBisectIterative")
      events.emplace_back(new PermBisectIterative(path,rng,actions,moveInputs[i],out));
    else if (type == "ShiftRefSlice")
      events.emplace_back(new ShiftRefSlice(path,rng,actions,moveInputs[i],out));
    else
      std::cerr << "Warning: Unrecognized Move, " << type << endl;
  }

  // Initialize Observables
  out.CreateGroup("Observables");
  vector<Input> observableInputs = in.getChild("Observables").getChildList("Observable");
  for (int i=0; i<observableInputs.size(); ++i) {
    string type = observableInputs[i].getAttribute<string>("type");
    if (type == "Energy")
      events.emplace_back(new Energy(path,actions,observableInputs[i],out));
    else if (type == "PairCorrelation")
      events.emplace_back(new PairCorrelation(path,observableInputs[i],out));
    else if (type == "PathDump")
      events.emplace_back(new PathDump(path,observableInputs[i],out));
    else if (type == "Permutation")
      events.emplace_back(new Permutation(path,observableInputs[i],out));
    else if (type == "Sign")
      events.emplace_back(new Sign(path,observableInputs[i],out));
    else if (type == "StructureFactor")
      events.emplace_back(new StructureFactor(path,observableInputs[i],out));
    else if (type == "Time")
      events.emplace_back(new Time(path,events,observableInputs[i],out));
    else
      std::cerr << "WARNING: Unrecognized observable, " << type << endl;
  }

  // Initialize Write
  events.emplace_back(new Writes(out,events,InterComm));

  // Initialize Algorithm
  vector<Input> loopList = in.getChild("Algorithm").getChildList("Loop");
  for (int i=0; i<loopList.size(); ++i)
    mainLoop.Init(loopList[i],events);
}


void Algorithm::Run()
{
  mainLoop.DoEvent();
}
