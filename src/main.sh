#!/bin/bash

for eps in 0.1 0.05 0.03 0.25 0.02 0.01
do
  for p in 1 2 3 4 5 6 
  do
    while [ $(jobs | wc -l) -ge 28 ]
    do 
      sleep 1
    done
    ./test_variance $eps $p 10000 &  
  done
done
