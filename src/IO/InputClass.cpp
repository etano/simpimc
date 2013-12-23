#include "InputClass.h"

// Loads settings structure from the specified XML file
void Input::load(const string &filename)
{
  xNode = XMLNode::openFileHelper(filename.c_str(),"Input");

  // Load the XML file into the property tree. If reading fails
  // (cannot open file, parse error), an exception is thrown.
  //read_xml(filename, pt);
}

