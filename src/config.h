#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <memory>
#include "Utils/scaffold.hpp"

using namespace std;
using namespace scaffold::matrix;
using namespace scaffold::algorithm;
using namespace scaffold::parallel;
using namespace scaffold::IO;
using namespace scaffold::rand;

inline double cmag2 (const complex<double> &z1, const complex<double> &z2) { return (z1.real()*z2.real() + z1.imag()*z2.imag()); }
inline double cmag (const complex<double> &z1, const complex<double> &z2) { return (sqrt(cmag2(z1,z2))); }

#endif
