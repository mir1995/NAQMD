#!/bin/bash

for s in 1 2 3 4 5 6 7 8 9 10
do
  for rate in lzadia lzdia sa sa1 sa2 sa3 sa12 sa13 sa23 sa123 
  do 
    while [ $(jobs | wc -l) -ge 20 ] 
    do 
      sleep 1 
    done
    ./nai_dynamics $rate $s & 
  done
done
