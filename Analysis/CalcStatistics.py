from numpy import *
import string
import re
from math import *

def C(g,t,mean,sigma2,N):
  if (sigma2 == 0):
    return 0
  total = 0.0
  diff = g - mean
  total = sum(diff[0:N-t]*diff[t:N])
  return (1.0/(sigma2*(N-t)))*total
  
def Kappa(g,mean,sigma2,N):
  total = 0.0
  for i in range(1,len(g)):
	  c = C(g,i,mean,sigma2,N)
	  if (c <= 0):
		  break
	  else:
	    total += c
  return 1.0 + 2.0*total

def Mean(g):
  return sum(g)/len(g)

def Var(g):
  mean = Mean(g)
  mean2 = Mean(g**2)
  return mean2 - mean**2
  
def Sigma(g):
  return sqrt(Var(g))
  
def Error(g):
  kappa = Kappa(g,Mean(g),Var(g),len(g))
  return Sigma(g)/sqrt(len(g)/kappa)
  
def StdError(g):
  return Sigma(g)/sqrt(len(g))

def stats(g):
  N = len(g)
  mean = sum(g)/N
  mean2 = sum(g**2)/N
  var = mean2 - mean**2
  sigma = sqrt(var)
  kappa = Kappa(g,mean,var,N)
  Neff = N/kappa
  return (mean,sigma,kappa,sigma/sqrt(Neff))
  
def getAndOutputStats(myArray, myArrayHeadings):
  col = []
  for i in range(0, len(myArray[0])):
    col.append(zeros(len(myArray))+0.0)
    for j in range(0, len(myArray)):
      col[i][j] = myArray[j,i]
  
  returnVal = []
  for i in range(1, len(col)):
    returnVal.append([])
    print "\n", myArrayHeadings[i], ":"
    (myMean,myStdDev,myKappa,myStdErr) = stats(col[i])  
    print 'Mean     = %8.5f' % (myMean) 
    print 'Standard Deviation     = %8.5f' % (myStdDev) 
    print 'Kappa     = %8.5f' % (myKappa) 
    print 'Standard Error    = %8.5f' % (myStdErr)
    returnVal[i-1].append(myMean)
    returnVal[i-1].append(myStdDev)
    returnVal[i-1].append(myKappa)
    returnVal[i-1].append(myStdErr)
    
  return returnVal
    
#take a set of data and if there are multiple
#x values that are the same, average them together
def HistogramAverageData(theData):
    newData=[]
    theData=theData.tolist()
    theData.sort()
    theData=array(theData)
    counter=0
    (lenX,lenY)=shape(theData)
    while (counter<lenX):
      tempData=[]
      tempData.append(theData[counter,1])
      counter=counter+1
      while (counter<lenX and theData[counter,0]-theData[counter-1,0]<1e-5):
           tempData.append(theData[counter,1])
           counter=counter+1
      newData.append([theData[counter-1,0],Mean(asarray(tempData)),Error(asarray(tempData))])
    return array(newData)
