#!/bin/bash

#for delta in 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.005
for npart in 100000 1000000 10000000 100000000
do 
  for seed in 1 2 3 4 5 6 7 8 9 10
  do
    while [ $(jobs | wc -l) -ge 30 ] 
      do 
        sleep 1
      done 
      ./lz_mass $npart 0.5 $seed &
  done
done
