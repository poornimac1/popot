#!/bin/bash

g++ -o SPSO Standard_PSO_2006.c -O3
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
