#include "InputClass.h"

// Loads settings structure from the specified XML file
void Input::load(const string &filename)
{
  xNode = XMLNode::openFileHelper(filename.c_str(), "Input");
}

