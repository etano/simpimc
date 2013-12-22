#include "AlgorithmClass.h"

void Algorithm::Init(Input &in, IOClass &out, RNG &rng)
{

  string outputPrefix = in.get<string>("Input.IO.OutputPrefix");

  // Initialize Path
  path.Init(in, out);

  // Initialize Actions
  vector<Input> actionInputs = in.getObjectList("Input.Actions");
  for (int i=0; i<actionInputs.size(); ++i) {
    string type = actionInputs[i].get<string>("Type");
    if (type == "Kinetic")
      actions.push_back(new Kinetic(path,actionInputs[i],out));
    else
      std::cerr << "Warning: Unrecognized Action, " << type << endl;
  }

  // Initialize Moves
  vector<Input> moveInputs = in.getObjectList("Input.Moves");
  for (int i=0; i<moveInputs.size(); ++i) {
    string type = moveInputs[i].get<string>("Type");
    if (type == "Bisect")
      events.push_back(new Bisect(path,rng,moveInputs[i],out));
    else
      std::cerr << "Warning: Unrecognized Move, " << type << endl;
  }

  // Initialize Observables
  vector<Input> observableInputs = in.getObjectList("Input.Observables");
  for (int i=0; i<observableInputs.size(); ++i) {
    string type = observableInputs[i].get<string>("Type");
    if (type == "Energy")
      events.push_back(new Energy(path,actions,observableInputs[i],out));
    else
      std::cerr << "WARNING: Unrecognized observable, " << type << endl;
  }

  // Initialize Write
  events.push_back(new Write(out,events));

  // Initialize Algorithm
  vector<Input> loopList = in.getObjectList("Input.Algorithm");
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
  double timeDif = difftime (end,start);

  mainLoop.DoEvent();
}
