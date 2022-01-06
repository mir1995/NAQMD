#!/bin/bash
for s in 1 2 3 4 5 6 7 8 9 10
do
  for delta in 1 0.8 0.6 0.4 0.2 0.1 0.08 0.06 0.04 0.02 0.01  
  do
    for eps in 0.1 0.08 0.06 0.04 0.02
    do
      for p in 2 4 6 
      do
        while [ $(jobs | wc -l) -ge 30 ] 
          do 
            sleep 1
          done 
          ./simple1_mass $eps $delta $p $s &  
      done
    done
  done
done
