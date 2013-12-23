#ifndef INPUTFILE_H
#define INPUTFILE_H

#include <string>
#include "../config.h"
#include "xmlParser.h"

template <class T>
inline T convertConstChar(const char * val) { return val; }
template <>
inline string convertConstChar(const char * val) { return string(val); }
template <>
inline RealType convertConstChar(const char * val) { return atof(val); }
template <>
inline int convertConstChar(const char * val) { return atoi(val); }

class Input
{
public:
  void load(const string &filename);
  XMLNode xNode;

  inline string getName() { return string(xNode.getName()); }

  inline Input getNode(string name)
  {
    Input in;
    in.xNode = xNode.getChildNode(name.c_str());
    return in;
  }

  template <class T>
  inline T getAttribute(string name)
  {
    const char* i = xNode.getAttribute(name.c_str());
    if (i == NULL) {
      cerr << "ERROR: " << name << " not found in input!" << endl;
      abort();
    } else
      return convertConstChar<T>(i);
  }

  template <class T>
  inline T getAttribute(string name, T deflt)
  {
    const char* i = xNode.getAttribute(name.c_str());
    if (i == NULL)
      return deflt;
    else
      return convertConstChar<T>(i);
  }

  inline vector<Input> getNodeList(string name) {
    vector<Input> ins;
    int n = xNode.nChildNode(name.c_str());
    for (int i=0; i<n; i++) {
      Input in;
      in.xNode = xNode.getChildNode(name.c_str(), i);
      ins.push_back(in);
    }
    return ins;
  }

  inline vector<Input> getChildList() {
    vector<Input> ins;
    int n = xNode.nChildNode();
    for (int i=0; i<n; i++) {
      Input in;
      in.xNode = xNode.getChildNode(i);
      ins.push_back(in);
    }
    return ins;
  }


};

#endif
