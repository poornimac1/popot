#!/usr/bin/python

import os

#g++ -o bench bench.cc -DGCC_ALGO=spso2011 -DGCC_RNG=CRNG -DGCC_PB=F0 -DGCC_DIM=10 `pkg-config --libs --cflags popot` -O3 -std=c++0x

Algo= ['spso2006','spso2011']
Pbs = ['F0', 'F1','F2', 'F3', 'F4','F5','F7','F8','F9']
Dims = ['30', '30', '30', '30', '2', '30', '40', '30', '40']

for alg in Algo:
    os.system('rm -f res_'+alg+'.data; touch res_'+alg+'.data');
    for pb, d in zip(Pbs, Dims):
        mystr = 'g++ -o bench bench.cc -DGCC_ALGO='+alg+' -DGCC_RNG=CRNG -DGCC_PB='+pb+' -DGCC_DIM='+d+' `pkg-config --libs --cflags popot` -O3 -std=c++0x'
        print mystr
        os.system(mystr)
        os.system('./bench >> res_'+alg+'.data')

