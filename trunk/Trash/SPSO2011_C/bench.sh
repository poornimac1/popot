#!/bin/bash

g++ -o SPSO main.c -O3 `pkg-config --libs --cflags gsl`
rm res.data -f
touch res.data
./SPSO 0 >> res.data
./SPSO 1 >> res.data
./SPSO 2 >> res.data
./SPSO 3 >> res.data
./SPSO 4 >> res.data
./SPSO 5 >> res.data
./SPSO 7 >> res.data
./SPSO 8 >> res.data
./SPSO 9 >> res.data
