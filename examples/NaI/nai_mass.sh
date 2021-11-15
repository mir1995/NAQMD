#!/bin/bash

for s in 2 3 4 5 6 7 8 9 10
do
  while [ $(jobs | wc -l) -ge 15 ] 
  do 
    sleep 1 
  done
  ./nai_mass $s & 
done
