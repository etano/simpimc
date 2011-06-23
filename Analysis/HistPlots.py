import ReadData 
import Plotting
import CalcStatistics
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import * 
import sys

# Histogram Variable Label
varLabel = str(sys.argv[1])

# Output Figure Label
outputLabel = str(sys.argv[2])

# Input File
inputFile = []
for i in range(0, len(sys.argv)-3):  
  inputFileName = "inputs/" + str(sys.argv[i+3])
  inputFile.append([])
  inputFile[i] = open(inputFileName,'r')

# Input Data
inputFileLabel = []
yValData = []
xValData = []
xVals = []
bins = []
fullHist = []
errHist = []
for i in range(0, len(inputFile)):
  # For every input file
  inputFileLabel.append([])
  yValData.append([])
  xValData.append([])
  xVals.append([])
  bins.append([])
  fullHist.append([])
  errHist.append([])

  firstLine = True
  j = 0
  for line in inputFile[i]:
    if firstLine:
    # First line label
      inputFileLabel[i] = str(line)
      firstLine = False
    else:
    # For every beta point
      inputLine = str(line)
      fileExtension = inputLine.replace(' ','-')
      fileExtension = fileExtension.replace('\n','')
      filePath = "data/traces/" + varLabel + "Trace-" + fileExtension + ".dat"    
      print "\nReading data from " + filePath + "."
      (myArray, myArrayHeadings) = ReadData.loadAscii(filePath)
      yValData[i].append(myArray)
      xValData[i].append(myArrayHeadings)

      # Make histograms and generate error bars
      print "\nSorting Data..."
  
      xVals[i].append([])
      for k in range(0, len(xValData[i][j])):
        # For every x value
        xVals[i][j].append([])          
        xVals[i][j][k] = float(xValData[i][j][k])
    
      bins[i].append([])
      for l in range(0, len(xVals[i][j])):
        # For every bin in the histogram
        bins[i][j].append([])
        for k in range(0, len(yValData[i][j])):
          # For every value in the input file line
          bins[i][j][l].append(yValData[i][j][k][l] + 0.0)
  
      fullHist[i].append([])
      for l in range(0, len(bins[i][j])):
        # For every bin in the histogram
        fullHist[i][j].append([])
        fullHist[i][j][l] = CalcStatistics.Mean(asarray(bins[i][j][l]))
  
      errHist[i].append([])
      for l in range(0, len(bins[i][j])):
        # For every bin in the histogram
        errHist[i][j].append([])
        errHist[i][j][l] = CalcStatistics.StdError(asarray(bins[i][j][l]))

      j += 1

# Generate Plots
print "\nGenerating Plots..."
plt.clf()
plt.xlabel("r")
plt.ylabel(varLabel)
for i in range(0, len(inputFile)):
  # For every input file 
  for j in range(0, len(yValData[i])):
    # For every line in the input file
    plt.errorbar(xVals[i][j], fullHist[i][j], errHist[i][j], label=inputFileLabel[i])
# Legend
leg = plt.legend(loc='best')
for t in leg.get_texts():
    t.set_fontsize('xx-small') 
plt.suptitle(varLabel, fontsize=12)
plt.savefig("data/figures/" + varLabel + "-" + outputLabel + ".png")
plt.clf()
print "\nPlot data/figures/" + varLabel + "-" + outputLabel + ".png Generated!"

print "\nDone.\n"
