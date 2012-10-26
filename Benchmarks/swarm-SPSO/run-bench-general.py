#!/usr/bin/python

import os
from math import * 
metaparameters = [[1.0/(2.0*log(2.0)),0.5 + log(2.0)], [0.729844,1.496180]]

Pbs = ['F0<30>', 'F1<30>','F2<30>', 'F3<30>', 'F4']#,'F5<30>','F7<40>','F8<30>','F9<40>'];

TopologyNames=['Full25','Ring25','VonNeuman25']
Topology = ['Full<25, Particle>', 'Ring<25, Particle>','VonNeuman<5,5,Particle>']
EvaluationMode = ['SYNCHRONOUS_EVALUATION','ASYNCHRONOUS_EVALUATION']

nb_variants = len(metaparameters)*len(Topology)*len(EvaluationMode)
print "I have to test %i variants on %i problems ... go and take some coffees" % (nb_variants, len(Pbs))

index = 0
for pb in Pbs:
    print "Problem : ",pb
    filename = 'res_'+pb.replace('<','_').replace('>','_')+'.data'
    os.system('rm -f '+filename+'; touch '+filename);
    for topo in Topology:
        for evalmode in EvaluationMode:
            for m in metaparameters:
                print(str(index))
                index+=1
                mystr = 'g++ -o bench_general bench_general.cc -D\'GCC_PB='+pb+'\' -DGCC_W=\''+str(m[0])+'\' -DGCC_C=\''+str(m[1])+'\' -DGCC_TOPOLOGY=\''+topo+'\' -DGCC_TOPOLOGY_NAME=\''+TopologyNames[Topology.index(topo)]+'\' -DGCC_EVALUATION_MODE=\''+evalmode+'\' `pkg-config --libs --cflags popot swarm` -O3'
                os.system(mystr)
                os.system('./bench_general >> ' + filename)
