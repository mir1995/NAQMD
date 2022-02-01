#! /bin/bash
#
gcc -c -Wall ../src/Quadrature/hermite_rule.c test_hermite.c
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gcc test_hermite.o hermite_rule.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm test_hermite.o
#
chmod ugo+x a.out
#
echo "Normal end of execution."
