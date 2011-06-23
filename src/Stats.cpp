#include "Stats.h"

double getMean ( const std::vector<double>& data ) {
  unsigned int N = data.size();
  double mean = 0.0;
  for (unsigned int i = 0; i < N; i += 1) mean += data[i]/(1.0*N);
  return mean;
}

double getVar ( const std::vector<double>& data ) {
  unsigned int N = data.size();
  std::vector<double> data2;
  for (unsigned int i = 0; i < N; i += 1) data2.push_back(data[i]*data[i]);
  double mean = getMean(data);
  double mean2 = getMean(data2);
  double var = (1.0*N/(N-1.0))*(mean2 - mean*mean);
  return var;
}

double getC( const std::vector<double>& data , unsigned int k , unsigned int N , double mean , double var ) {
  double delta[N];
  double tot = 0.0;
  for (unsigned int i = 0; i < N; i += 1) delta[i] = data[i] - mean;
  for (unsigned int i = 0; i < N-k; i += 1) tot += delta[i]*delta[i+k];
  return tot/((N-k)*var);
}

double getKappa( const std::vector<double>& data ) {
  unsigned int N = data.size();
  bool done = 0;
  unsigned int k = 1;
  double mean = getMean(data);
  double var = getVar(data);
  double c, csum = 1.0;
  while ((!done) && (k < N)) {
    c = getC(data,k,N,mean,var);
    if (c < 0.0) done = 1;
    else csum += 2.0*c;
    k += 1;
  }
  return csum;
}

double getError( const std::vector<double>& data ) {
  int N = data.size();
  double var = getVar(data);
  double kappa = getKappa(data);
  return sqrt(kappa*var/N);
}

void statsEnergy ( const char* energyFile , unsigned int nType , unsigned int nBlock ) {

  std::ifstream energyTrace;

  unsigned int count; 
  double E, KE, VE, NE;
  std::string sE, sKE, sVE, sNE;
  std::vector<double> Evec[nType], KEvec[nType], VEvec[nType], NEvec[nType];
  std::vector<double> Etotvec, KEtotvec, VEtotvec, NEtotvec;
    
  energyTrace.open(energyFile);

  count = 0;
  while ( count < nBlock ) { //!energyTrace.eof() ) { // keep reading until enD-of-file !inData.eof()
    if (count == 0) {
      for (unsigned int iType = 0; iType < nType; iType += 1)
        energyTrace >> sE >> sKE >> sVE >> sNE;
      energyTrace >> sE >> sKE >> sVE >> sNE;
    } else {
      for (unsigned int iType = 0; iType < nType; iType += 1) {
        energyTrace >> E >> KE >> VE >> NE;
        Evec[iType].push_back(E);
        KEvec[iType].push_back(KE);
        VEvec[iType].push_back(VE);
        NEvec[iType].push_back(NE);
      }
      energyTrace >> E >> KE >> VE >> NE;
      Etotvec.push_back(E); 
      KEtotvec.push_back(KE); 
      VEtotvec.push_back(VE); 
      NEtotvec.push_back(NE); 
    }
    
    count += 1;
  }  
  
  std::cout << "\nEnergy estimates:\n";
  std::cout << "E: " << getMean(Etotvec) << " (" << getError(Etotvec) << ")\n";
  std::cout << "T: " << getMean(KEtotvec) << " (" << getError(KEtotvec) << ")\n";
  std::cout << "V: " << getMean(VEtotvec) << " (" << getError(VEtotvec) << ")\n";
  std::cout << "N: " << getMean(NEtotvec) << " (" << getError(NEtotvec) << ")\n";
  
  std::cout << "\nEnergy Per Permutation Configuration, (E, T, V, N):\n";
  for (unsigned int iType = 0; iType < nType; iType += 1) {
    std::cout << getMean(Evec[iType]) << " " << getMean(KEvec[iType]) << " " << getMean(VEvec[iType]) << " " << getMean(NEvec[iType]) << "\n";
  }  
  
  // Close files
  energyTrace.close();
 
}

void statsR ( const char* rFile , unsigned int nType , unsigned int nBlock , unsigned int nPart , unsigned int nD ) {

  std::ifstream rTrace;

  unsigned int count; 
  double r;
  double totr[nType][nPart];
  double avgr[nType][nPart];
  double totavgr[nPart];
  std::string sR;
  
  rTrace.open(rFile);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iType = 0; iType < nType; iType++ ) {
      totr[iType][iPart] = 0;
    }
    totavgr[iPart] = 0;
  }
  
  count = 0;
  while ( count < nBlock ) { // keep reading until enD-of-file !inData.eof()
    if (count == 0) {
      for (unsigned int iType = 0; iType < nType; iType += 1) {
        for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
          rTrace >> sR;
        }
      }
    } else {
      for (unsigned int iType = 0; iType < nType; iType += 1) {
        for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
          rTrace >> r;
          totr[iType][iPart] += r;
        }       
      }
    }
      
    count += 1;  
  }

  for (unsigned int iType = 0; iType < nType; iType++ ) {
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      avgr[iType][iPart] = totr[iType][iPart]/(nBlock*1.0);
      totavgr[iPart] += avgr[iType][iPart];
    } 
  }
  
  std::cout << "\nR estimate:\n";
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    std::cout << iPart << " : " << totavgr[iPart] << "\n";
  }    
  
  // Close files
  rTrace.close();
}
 
void statsRR ( const char* rrFile , unsigned int nType , unsigned int nBlock , unsigned int nPart ) {

  std::ifstream rrTrace;

  unsigned int count; 
  double rr;  
  double totrr[nType][nPart];
  double avgrr[nType][nPart];
  double totavgrr[nPart];
  std::string sRR;
    
  rrTrace.open(rrFile);
  
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    for (unsigned int iType = 0; iType < nType; iType++ ) {
      totrr[iType][iPart] = 0;
    }
    totavgrr[iPart] = 0;
  }
  
  count = 0;
  while ( count < nBlock ) { // keep reading until enD-of-file !inData.eof()
    if (count == 0) {
      for (unsigned int iType = 0; iType < nType; iType += 1) {
        for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
          rrTrace >> sRR;
        }
      }
    } else {
      for (unsigned int iType = 0; iType < nType; iType += 1) {
        for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
          rrTrace >> rr;
          totrr[iType][iPart] += rr;
        }   
           
      }
    }
    count += 1;   
  }
  
  for (unsigned int iType = 0; iType < nType; iType++ ) {
    for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
      avgrr[iType][iPart] = totrr[iType][iPart]/(nBlock*1.0);
      totavgrr[iPart] += avgrr[iType][iPart];
    } 
  }
  
  std::cout << "\nRR estimate:\n";
  for (unsigned int iPart = 0; iPart < nPart; iPart += 1) {
    std::cout << iPart << " : " << totavgrr[iPart] << "\n";
  }    
  
  // Close files
  rrTrace.close();
}
