#!/bin/sh
firstVar="1"
for var in "$@"
do
  if [ $firstVar -eq "1" ]
  then
    echo "Using" $1 "Processes"
    firstVar=0
  else
    nLines=$(wc -l < "inputs/$var" )
    nInputLines=$(($nLines-1))
    echo "$var" "has" $nInputLines "input lines"

    line="2"
    echo "$line " > BatchRun.txt
    while [ $line -lt $nLines ]
    do
      line=$[$line+1]
      echo "$line " >> BatchRun.txt
    done

    cat BatchRun.txt | xargs -n1 -P$1 sh Run.sh "$var"
  fi
done
