#!/bin/bash

for delta in 1 0.5 0.05 
do
  for s in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18
  do 
    ./test_evolution $delta 0.1 4 100000 $s  
  done
done
