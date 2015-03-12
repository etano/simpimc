#ifndef SIMPIMC_CONFIG_H_
#define SIMPIMC_CONFIG_H_

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <memory>
#include <sys/time.h>
#include "utils/scaffold.h"

using namespace scaffold::matrix;
using namespace scaffold::algorithm;
using namespace scaffold::parallel;
using namespace scaffold::io;
using namespace scaffold::rand;

inline double CMag2 (const std::complex<double> &z1, const std::complex<double> &z2) { return (z1.real()*z2.real() + z1.imag()*z2.imag()); }
inline double CMag (const std::complex<double> &z1, const std::complex<double> &z2) { return (sqrt(CMag2(z1,z2))); }

#endif // SIMPIMC_CONFIG_H_
