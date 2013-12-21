#include "InputClass.h"

// Loads settings structure from the specified XML file
void Input::load(const std::string &filename)
{
  // Load the XML file into the property tree. If reading fails
  // (cannot open file, parse error), an exception is thrown.
  read_xml(filename, pt);
}

