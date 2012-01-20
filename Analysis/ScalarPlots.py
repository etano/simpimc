import ReadData 
import Plotting
import CalcStatistics
import math
from numpy import *
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import sys

# Output Figure Label
outputLabel = str(sys.argv[1])

# Input File
inputFile = []
for i in range(0, len(sys.argv)-2):  
  inputFileName = "inputs/" + str(sys.argv[i+2])
  print inputFileName
  inputFile.append([])
  inputFile[i] = open(inputFileName,'r')

# Input Data
inputFileLabel = []
scalarData = []
beta = []
for i in range(0, len(inputFile)):
  # For every input file
  firstLine = True
  inputFileLabel.append([])
  scalarData.append([])
  beta.append([])
  for line in inputFile[i]:
    if firstLine:
    # First line label
      inputFileLabel[i] = str(line)
      firstLine = False
    else:
    # For every beta point
      inputLine = str(line)

      entry = ''
      params = []
      for k in range(len(line)):
        if line[k]!=' ':
          entry += line[k]
        else:
          params += [entry]
          entry=''
      beta[i].append(float(params[3]))
      nPart = float(params[0])
      nD = float(params[1])
      #interaction = float(params[11])

      fileExtension = inputLine.replace(' ','-')
      fileExtension = fileExtension.replace('\n','')
      scalarFilePath = "data/traces/scalarTrace-" + fileExtension + ".dat"    
      scalarFilePath = "data/traces/Energy-7-3-32-2.7366-7200-1-1-1-0.dat"
      print "\nReading data from " + scalarFilePath + "."
      (myArray, myArrayHeadings) = ReadData.loadAscii(scalarFilePath)
      scalarData[i].append(CalcStatistics.getAndOutputStats(myArray, myArrayHeadings))
      #Plotting.makePlots(myArray, myArrayHeadings, fileExtension)

print "\nSorting Data..."
# Rotate Data into Columns
col = []
for i in range(0, len(scalarData)):
  # For every input file
  col.append([])
  for j in range(0, len(scalarData[i][0])):
    # For every observable
    col[i].append([])
    for k in range(0, len(scalarData[i][0][0])):
      # For every statistic
      col[i][j].append([])
      for l in range(0, len(scalarData[i])):
        # For every beta point
        col[i][j][k].append([])
        col[i][j][k][l] = scalarData[i][l][j][k]

# Generate Plots
print "\nGenerating Plots..."

# Normal Plots
for i in range(0, len(col[0])):
  # For every observable
  plt.xlabel("Beta")
  plt.ylabel(myArrayHeadings[i+1])
  for j in range(0, len(col)):
    # For every input file 
    plt.errorbar(beta[j], col[j][i][0], col[j][i][3], label=inputFileLabel[j])

  # Exact Values
  #x = arange(1.0,10.0,0.01)
  #if (interaction==0):  
  #  if ((myArrayHeadings[i+1]=="PE") or (myArrayHeadings[i+1]=="VE")):  
  #    y =  nPart*nD*0.5/tanh(x/2.0)
  #    plt.plot(x, y, label="Exact")
  #  if (myArrayHeadings[i+1]=="R2"):  
  #    y = nD*0.5/tanh(x/2.0)
  #    plt.plot(x, y, label="Exact")

  # Legend
  leg = plt.legend(loc='best')
  # matplotlib.text.Text instances
  for t in leg.get_texts():
      t.set_fontsize('xx-small')    # the legend text fontsize
  plt.suptitle(myArrayHeadings[i+1] + " vs Beta", fontsize=12)
  plt.savefig("data/figures/" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png")
  plt.clf()
  print "\nPlot data/figures/" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png Generated!"

print "\nGenerating Error Plots..."
# Error Plots
for i in range(0, len(col[0])):
  # For every observable
  plt.xlabel("Beta")
  plt.ylabel(myArrayHeadings[i+1])
  for j in range(0, len(col)):
    # For every input file 
    plt.plot(beta[j], col[j][i][3], label=inputFileLabel[j])
  leg = plt.legend(loc='best')
  for t in leg.get_texts():
      t.set_fontsize('xx-small')    # the legend text fontsize
  plt.suptitle("Error in " + myArrayHeadings[i+1] + " vs Beta", fontsize=12)
  plt.savefig("data/figures/Err" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png")
  plt.clf()
  print "\nPlot data/figures/Err" + myArrayHeadings[i+1] + "vBeta-" + outputLabel + ".png Generated!"

print "\nDone.\n"
