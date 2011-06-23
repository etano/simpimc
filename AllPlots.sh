#!/bin/sh
varCount="1"
vars=""
for var in "$@"
do
  if [ $varCount -eq "1" ]
  then
    label=$var
    varCount=2
  else
    vars+="$var "
  fi
done
python ScalarPlots.py $label $vars
python HistPlots.py "RDen" $label $vars
python HistPlots.py "Grr" $label $vars
