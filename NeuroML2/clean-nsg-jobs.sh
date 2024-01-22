#!/bin/bash

# Copyright 2023 Ankur Sinha
# Author: Ankur Sinha <sanjay DOT ankur AT gmail DOT com> 
# File : 
#


nsgr_job -l | cut -f2 -d ' ' > joblist.txt

while IFS="" read -r line ; do
    echo "Deleting job: ${line}"
    nsgr_job -j ${line} -r
done < joblist.txt

echo "All jobs deleted"
echo "Please delete joblist.txt"
