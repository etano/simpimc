#include "AlgorithmClass.h"

void Algorithm::Init(Input &in, IOClass &out, RNG &rng)
{

  string outputPrefix = in.get<string>("Input.IO.OutputPrefix");

  // Initialize Path
  path.Init(in, out);

  // Initialization Moves
  vector<Input> moveList = in.getObjectList("Input.Moves");
  for (int i=0; i<moveList.size(); ++i) {
    string type = moveList[i].get<string>("Type");
    if (type == "Bisect")
      events.push_back(new Bisect(path,rng,moveList[i],out));
    else
      std::cerr << "Warning: Unrecognized Move, " << type << endl;
  }

  // Initialize Observables
  vector<Input> observableList = in.getObjectList("Input.Observables");
  for (int i=0; i<observableList.size(); ++i) {
    string type = observableList[i].get<string>("Type");
    if (type == "Energy")
      events.push_back(new Energy(path,observableList[i],out));
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
