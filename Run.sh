#!/bin/sh
echo "Running simulation" $(($2-1)) "from" $1 "..."
./shab-pimd "inputs/$1" "$(($2-1))" > "data/output/$1-$(($2-1)).out"

