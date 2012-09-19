#!/usr/bin/python

import os



Algo= ['SPSO2006','SPSO2011']
Pbs = ['F0<30>', 'F1<30>','F2<30>', 'F3<30>', 'F4','F5<30>','F7<40>','F8<30>','F9<40>'];
rng = ['CRNG','MersenneRNG'] # KissRNG (only on 32 bits)
for alg in Algo:
    os.system('rm -f res_'+alg+'.data; touch res_'+alg+'.data');
    for rn in rng:
        for pb in Pbs:
            mystr = 'g++ -o bench bench.cc -D\'GCC_RNG=popot::rng::'+rn+'\' -D\'GCC_PB=BenchmarkProblems::'+pb+'\' -D\'GCC_ALGO=popot::PSO::'+alg+'::PSO<Problem>::Type\' `pkg-config --libs --cflags popot` -O3'
            print mystr
            os.system(mystr)
            os.system('./bench >> res_'+alg+'.data')

