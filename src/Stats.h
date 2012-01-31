#ifndef STATS_H
#define STATS_H

#include "StandardLibraries.h"

double getMean ( const std::vector<double>& data );

double getVar ( const std::vector<double>& data );

double getC( const std::vector<double>& data , unsigned int k , unsigned int N , double mean , double var );

double getKappa( const std::vector<double>& data );

double getError( const std::vector<double>& data );

double jackKnife( const std::vector<double>& data );

void statsEnergy ( const char* energyfile , unsigned int nType , unsigned int nBlock );

void statsR ( const char* rFile , unsigned int nType , unsigned int nBlock , unsigned int nPart , unsigned int nD );
 
void statsRR ( const char* rrFile , unsigned int nType , unsigned int nBlock , unsigned int nPart );

#endif
