#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>
#include "Utils/types.h"

using namespace std;

inline RealType cmag2 (const ComplexType &z1, const ComplexType &z2) { return (z1.real()*z2.real() + z1.imag()*z2.imag()); }
inline RealType cmag (const ComplexType &z1, const ComplexType &z2) { return (sqrt(cmag2(z1,z2))); }

#endif
