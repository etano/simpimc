#ifndef INPUTFILE_H
#define INPUTFILE_H

#include <string>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include "../config.h"

class Input
{
public:
  void load(const std::string &filename);
  string val;
  boost::property_tree::ptree pt;

  template <class T>
  inline T get(string name) { return pt.get<T>(name); }

  template <class T>
  inline T get(string name, T deflt) { return pt.get<T>(name,deflt); }

  template <class T>
  inline T getValue() { return pt.get_value<T>(); }

  inline string getKey() { return val; }

  template <class T>
  inline std::vector<T> getList(string name) {
    vector<T> vals;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(name))
      vals.push_back(v.second.data());
    return vals;
  }

  inline vector<Input> getObjectList(string name) {
    vector<Input> vals;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(name)) {
      Input tmpIn;
      tmpIn.pt = v.second;
      vals.push_back(tmpIn);
    }
    return vals;
  }

  inline void convertToVector(vector<pair<string,string> > &v) {
    boost::property_tree::ptree::const_iterator end = pt.end();
    for (boost::property_tree::ptree::const_iterator it = pt.begin(); it != end; ++it) {
      Input tmpIn;
      v.push_back(make_pair(it->first,it->second.get_value<std::string>()));
      tmpIn.pt = it->second;
      tmpIn.convertToVector(v);
    }
  }

  inline void getInList(vector<Input*> &inList) {
    boost::property_tree::ptree::const_iterator end = pt.end();
    for (boost::property_tree::ptree::const_iterator it = pt.begin(); it != end; ++it) {
      Input *tmpIn = new Input();
      tmpIn->val = it->first;
      tmpIn->pt = it->second;
      inList.push_back(tmpIn);
    }
  }

  inline Input getObject(string parent, string child) {
    Input tmpIn;
    BOOST_FOREACH(boost::property_tree::ptree::value_type &v, pt.get_child(parent)) {
      if (v.first == child)
        tmpIn.pt = v.second;
    }
    return tmpIn;
  }

};

#endif
