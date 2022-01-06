#!/bin/bash

for npart in 10000 100000 1000000 10000000
do
  for rate in lzadia lzdia sa sa1 sa2 sa3 sa12 sa13 sa23 sa123 
  do 
    while [ $(jobs | wc -l) -ge 20 ] 
    do 
      sleep 1 
    done
    ./nai_dynamics2 $npart $rate & 
  done
done
