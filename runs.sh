#!/bin/sh
cat runs.txt | xargs -n8 -P6 sh doit.sh
