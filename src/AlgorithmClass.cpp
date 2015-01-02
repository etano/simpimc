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
    if (type == "Kinetic") {
      actions.push_back(std::make_shared<Kinetic>(path,actionInput,out));
    } else if (type == "HarmonicTrap") {
      actions.push_back(std::make_shared<Trap>(path,actionInput,out));
    } else if (type == "FreeNodal") {
      actions.push_back(std::make_shared<FreeNodal>(path,rng,actionInput,out));
    } else if (type == "OptimizedNodal") {
      actions.push_back(std::make_shared<OptimizedNodal>(path,rng,actionInput,out));
    } else if (type == "BarePairAction") {
      actions.push_back(std::make_shared<BarePairAction>(path,actionInput,out));
    } else if (type == "DavidPairAction") {
      actions.push_back(std::make_shared<DavidPairAction>(path,actionInput,out));
    } else if (type == "IlkkaPairAction") {
      actions.push_back(std::make_shared<IlkkaPairAction>(path,actionInput,out));
    } else if (type == "ImportancePairAction") {
      actions.push_back(std::make_shared<ImportancePairAction>(path,actionInput,out));
    } else
      std::cerr << "Warning: Unrecognized Action, " << type << endl;
  }

  // Initialize Moves
  out.CreateGroup("Moves");
  vector<Input> moveInputs = in.getChild("Moves").getChildList("Move");
  for (auto& moveInput: moveInputs) {
    string type = moveInput.getAttribute<string>("type");
    if (type == "Bisect")
      events.push_back(std::make_shared<Bisect>(path,rng,actions,moveInput,out));
    else if (type == "DisplaceParticle")
      events.push_back(std::make_shared<DisplaceParticle>(path,rng,actions,moveInput,out));
    else if (type == "PermBisect")
      events.push_back(std::make_shared<PermBisect>(path,rng,actions,moveInput,out));
    else if (type == "PermBisectIterative")
      events.push_back(std::make_shared<PermBisectIterative>(path,rng,actions,moveInput,out));
    else if (type == "ShiftRefSlice")
      events.push_back(std::make_shared<ShiftRefSlice>(path,rng,actions,moveInput,out));
    else if (type == "VaryOptimizedNodal")
      events.push_back(std::make_shared<VaryOptimizedNodal>(path,rng,actions,moveInput,out));
    else
      std::cerr << "Warning: Unrecognized Move, " << type << endl;
  }

  // Initialize Observables
  out.CreateGroup("Observables");
  vector<Input> observableInputs = in.getChild("Observables").getChildList("Observable");
  for (auto& observableInput: observableInputs) {
    string type = observableInput.getAttribute<string>("type");
    if (type == "Energy")
      events.push_back(std::make_shared<Energy>(path,actions,observableInput,out));
    else if (type == "ImportanceWeight")
      events.push_back(std::make_shared<ImportanceWeight>(path,actions,observableInput,out));
    else if (type == "PairCorrelation")
      events.push_back(std::make_shared<PairCorrelation>(path,observableInput,out));
    else if (type == "PathDump")
      events.push_back(std::make_shared<PathDump>(path,observableInput,out));
    else if (type == "Permutation")
      events.push_back(std::make_shared<Permutation>(path,observableInput,out));
    else if (type == "RecordOptimizedNodal")
      events.push_back(std::make_shared<RecordOptimizedNodal>(path,actions,observableInput,out));
    else if (type == "Sign")
      events.push_back(std::make_shared<Sign>(path,observableInput,out));
    else if (type == "StructureFactor")
      events.push_back(std::make_shared<StructureFactor>(path,observableInput,out));
    else if (type == "Time")
      events.push_back(std::make_shared<Time>(path,events,observableInput,out));
    else
      std::cerr << "WARNING: Unrecognized observable, " << type << endl;
  }

  // Initialize Write
  events.push_back(std::make_shared<Writes>(out,events,InterComm));

  // Initialize Algorithm
  vector<Input> loopList = in.getChild("Algorithm").getChildList("Loop");
  for (auto& loop: loopList)
    mainLoop.Init(loop,events);
}


void Algorithm::Run()
{
  mainLoop.DoEvent();
}
