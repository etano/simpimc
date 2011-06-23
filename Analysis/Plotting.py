import numpy
from numpy import *
from math import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Makes a plot of the pairs of random numbers in 2D  
def Plot2d(xList,yList,plotTitle,plotFileName,xLabel,yLabel,xMin,xMax,yMin,yMax):
    
  plt.xlabel(xLabel)
  plt.ylabel(yLabel)
  plt.plot(xList, yList)
  plt.suptitle(plotTitle, fontsize=12)
  plt.savefig("data/figures/" + plotFileName + ".png")
  plt.clf()
  
  return None

# Makes a plot of the pairs of random numbers in 2D  
def Plot2dError(xList,yList,errList,plotTitle,plotFileName,xLabel,yLabel,xMin,xMax,yMin,yMax):
    
  plt.xlabel(xLabel)
  plt.ylabel(yLabel)
  plt.xlim(xMin,xMax) 
  plt.ylim(yMin,yMax) 
  plt.errorbar(xList, yList, errList)
  plt.suptitle(plotTitle, fontsize=12)
  plt.savefig("data/figures/" + plotFileName + ".png")
  plt.clf()
  
  return None
  
# Makes a plot of the pairs of random numbers in 2D  
def Plot3d(xList,yList,zList,plotTitle,plotFileName,xLabel,yLabel,zLabel):

  fig = plt.figure()
  ax = Axes3D(fig)
      
  ax.scatter(xList, yList, zList)
  ax.set_xlabel(xLabel)
  ax.set_ylabel(yLabel)
  ax.set_zlabel(zLabel)
  fig.suptitle(plotTitle, fontsize=12)
  fig.savefig("data/figures/" + plotFileName + ".png")
  plt.clf()
  
  return None
  
def makePlots(myArray, myArrayHeadings, filename):
  col = []
  for i in range(0, len(myArray[0])):
    col.append(zeros(len(myArray))+0.0)
    for j in range(0, len(myArray)):
      col[i][j] = myArray[j,i]
  
  for i in range(1, len(col)):
    Plot2d(col[0],col[i],myArrayHeadings[i]+" vs "+myArrayHeadings[0],myArrayHeadings[i]+filename,myArrayHeadings[0],myArrayHeadings[i],col[0][0],col[0][len(col)],col[i][0],col[i][len(col)])
    print "\n", myArrayHeadings[i], ": Plot Created!"
    
def make2dErrorPlot(myArray, myArrayHeadings, filename):
  col = []
  for i in range(0, len(myArray[0])):
    col.append(zeros(len(myArray))+0.0)
    for j in range(0, len(myArray)):
      col[i][j] = myArray[j,i]
      
  Plot2dError(col[0],col[1],col[2],myArrayHeadings[1]+" vs "+myArrayHeadings[0],myArrayHeadings[1]+filename,myArrayHeadings[0],myArrayHeadings[1],col[0][0],col[0][len(col)-1],col[1][0],col[1][len(col)-1])
  print "\n", myArrayHeadings[1], ": Plot Created!"
