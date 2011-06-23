from numpy import *
import string
import re

def removeBlanksFromList(theList):
    newList=[]
    for counter in range(0,len(theList)):
        if theList[counter]!='':
            newList.append(theList[counter])
    return newList

def load(fileName):
    loadAscii(fileName)

def saveAscii(theData,fileName):
    fileOut=open(fileName,'w')
    for counter in range(0,len(theData)):
        for counter2 in range(0,len(theData[counter])):
            print >>fileOut, theData[counter,counter2],
        print >>fileOut

        

#Reads in an ascii file of data much like octave
def loadAscii(fileName):
    infile=open(fileName)
    theLines=infile.readlines()
    splitter=re.compile('\s|,')
    lineSize=len(removeBlanksFromList(splitter.split(theLines[0])))
    numLines=len(theLines)
    if lineSize>1:
        newData=zeros((numLines-1,lineSize))+0.0
    else:
        newData=zeros(numLines-1)+0.0
    
    theLine=removeBlanksFromList(splitter.split(theLines[0]))
    headings = []
    for counter2 in range(0,lineSize):
        headings.append(theLine[counter2])
    for counter in range(1,numLines):
        theLine=removeBlanksFromList(splitter.split(theLines[counter]))
        if len(theLine)!=lineSize:
            print "All lines are not the same size in ",fileName
            #abort()
        else:
            if lineSize>1:
                for counter2 in range(0,lineSize):
                    newData[counter-1,counter2]=string.atof(theLine[counter2])
            else:
                for counter2 in range(0,lineSize):
                    newData[counter-1]=string.atof(theLine[counter2])
    return newData, headings

#Prints a 2d real array of data
def Print2d(theData):
    for counter in range(0,len(theData)):
        for counter2 in range(0,len(theData[counter])):
            print theData[counter,counter2],
        print
