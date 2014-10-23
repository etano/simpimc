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
    if (type == "Coulomb")
      actions.push_back(new Coulomb(path,actionInputs[i],out));
    else if (type == "Kinetic")
      actions.push_back(new Kinetic(path,actionInputs[i],out));
    else if (type == "HarmonicTrap")
      actions.push_back(new Trap(path,actionInputs[i],out));
    else if (type == "Nodal")
      actions.push_back(new Nodal(path,actionInputs[i],out));
    else if (type == "DavidPairAction")
      actions.push_back(new DavidPairAction(path,actionInputs[i],out));
    else if (type == "IlkkaPairAction")
      actions.push_back(new IlkkaPairAction(path,actionInputs[i],out));
    else
      std::cerr << "Warning: Unrecognized Action, " << type << endl;
  }

  // Initialize Moves
  out.CreateGroup("Moves");
  vector<Input> moveInputs = in.getChild("Moves").getChildList("Move");
  for (int i=0; i<moveInputs.size(); ++i) {
    string type = moveInputs[i].getAttribute<string>("type");
    if (type == "Bisect")
      events.push_back(new Bisect(path,rng,actions,moveInputs[i],out));
    else if (type == "PermBisect")
      events.push_back(new PermBisect(path,rng,actions,moveInputs[i],out));
    else if (type == "PermBisectIterative")
      events.push_back(new PermBisectIterative(path,rng,actions,moveInputs[i],out));
    else
      std::cerr << "Warning: Unrecognized Move, " << type << endl;
  }

  // Initialize Observables
  out.CreateGroup("Observables");
  vector<Input> observableInputs = in.getChild("Observables").getChildList("Observable");
  for (int i=0; i<observableInputs.size(); ++i) {
    string type = observableInputs[i].getAttribute<string>("type");
    if (type == "Energy")
      events.push_back(new Energy(path,actions,observableInputs[i],out));
    else if (type == "PairCorrelation")
      events.push_back(new PairCorrelation(path,observableInputs[i],out));
    else if (type == "PathDump")
      events.push_back(new PathDump(path,observableInputs[i],out));
    else if (type == "StructureFactor")
      events.push_back(new StructureFactor(path,observableInputs[i],out));
    else if (type == "Time")
      events.push_back(new Time(path,events,observableInputs[i],out));
    else
      std::cerr << "WARNING: Unrecognized observable, " << type << endl;
  }

  // Initialize Write
  events.push_back(new Writes(out,events));

  // Initialize Algorithm
  vector<Input> loopList = in.getChild("Algorithm").getChildList("Loop");
  for (int i=0; i<loopList.size(); ++i)
    mainLoop.Init(loopList[i],events);
}


void Algorithm::Run()
{
  mainLoop.DoEvent();
}
