#!/bin/bash

for delta in 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.005
do 
  while [ $(jobs | wc -l) -ge 5 ] 
    do 
      sleep 1
    done 
    ./lz_mass $delta &
done
