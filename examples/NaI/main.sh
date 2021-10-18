#!/bin/bash

for npart in 10000000 100000000 1000000000 10000000000  
do
  for s in 1 2 3 4 5 6 7 8 9 10
  do
    for rate in goddard goddard1 goddard2 goddard3 lasser belayev 
      do 
        while [ $(jobs | wc -l) -ge 20 ] 
        do 
          sleep 1 
        done
        ./many $npart $s $rate & 
      done
  done
done
