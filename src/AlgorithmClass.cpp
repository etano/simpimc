#include "AlgorithmClass.h"

void Algorithm::Init(Input &in, IOClass &out, RNG &rng)
{

  string outputPrefix = in.getChild("IO").getAttribute<string>("outputPrefix");

  // Initialize Path
  path.Init(in, out);

  // Initialize Actions
  out.CreateGroup("Actions");
  vector<Input> actionInputs = in.getChild("Actions").getChildList("Action");
  for (int i=0; i<actionInputs.size(); ++i) {
    string type = actionInputs[i].getAttribute<string>("type");
    if (type == "Kinetic")
      actions.push_back(new Kinetic(path,actionInputs[i],out));
    else if (type == "Trap")
      actions.push_back(new Trap(path,actionInputs[i],out));
    else
      std::cerr << "Warning: Unrecognized Action, " << type << endl;
  }

  // Initialize Moves
  vector<Input> moveInputs = in.getChild("Moves").getChildList("Move");
  for (int i=0; i<moveInputs.size(); ++i) {
    string type = moveInputs[i].getAttribute<string>("type");
    if (type == "Bisect")
      events.push_back(new Bisect(path,rng,actions,moveInputs[i],out));
    else
      std::cerr << "Warning: Unrecognized Move, " << type << endl;
  }

  // Initialize Observables
  vector<Input> observableInputs = in.getChild("Observables").getChildList("Observable");
  for (int i=0; i<observableInputs.size(); ++i) {
    string type = observableInputs[i].getAttribute<string>("type");
    if (type == "Energy")
      events.push_back(new Energy(path,actions,observableInputs[i],out));
    else
      std::cerr << "WARNING: Unrecognized observable, " << type << endl;
  }

  // Initialize Write
  events.push_back(new Write(out,events));

  // Initialize Algorithm
  vector<Input> loopList = in.getChild("Algorithm").getChildList("Loop");
  for (int i=0; i<loopList.size(); ++i)
    mainLoop.Init(loopList[i],events);
}


void Algorithm::Run()
{
  ////////////////////////
  /* Main Loop Settings */
  ////////////////////////

  // Set up time loop
  time_t start, end;
  time (&start);
  time (&end);
  RealType timeDif = difftime (end,start);

  mainLoop.DoEvent();
}
