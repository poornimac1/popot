#!/usr/bin/python

import os
from math import * 

Pbs = ['F0<30>', 'F1<30>','F2<30>', 'F3<30>', 'F4']#,'F5<30>','F7<40>','F8<30>','F9<40>'];
rng = 'CRNG'

MetaParameters = [[1.0/(2.0*log(2.0)),0.5 + log(2.0)], [0.729844,1.496180]]
PositionInit = ['PositionUniformRandom']
VelocityInit = ['VelocityZero']#,'VelocityHalfDiff','VelocitySPSO2011']
ParticleType = ['SPSO2006Particle']#,'SPSO2011Particle']
TopologyNames=['Full25_True','Full25_False','Ring25_True','Ring25_False','VonNeuman25_True','VonNeuman25_False'] #,'RandomInformants25','FixedRandomInformants25','AdaptiveRandom25']
Topology = ['Full<25, Particle, true>', 'Full<25, Particle, false>', 'Ring<25, Particle, true>', 'Ring<25, Particle, false>','VonNeuman<25,Particle,true>', 'VonNeuman<25,Particle,false>']#,'RandomInformants<25 , 3, Particle>','FixedRandomInformants<25 , 3, Particle>','AdaptiveRandom<25, 3, Particle>']
Random_Shuffle = ['true','false']
EvaluationMode = ['SYNCHRONOUS_EVALUATION','ASYNCHRONOUS_EVALUATION']

nb_variants = len(MetaParameters)*len(PositionInit)*len(VelocityInit)*len(ParticleType)*len(Topology)*len(Random_Shuffle)*len(EvaluationMode)
print "I have to test %i variants on %i problems ... go and take some coffees" % (nb_variants, len(Pbs))

index = 0
for pb in Pbs:
    print "Problem : ",pb
    filename = 'res_'+pb.replace('<','_').replace('>','_')+'.data'        
    os.system('rm -f '+filename+'; touch '+filename);
    for topo in Topology:
        for particule in ParticleType:
            for pos in PositionInit:
                for vel in VelocityInit:
                    for evalmode in EvaluationMode:
                        for randShuffle in Random_Shuffle:
                            for m in MetaParameters:
                                w_value = m[0]
                                c_value = m[1]
                                print(str(index))
                                index+=1
                                mystr = 'g++ -o bench_general bench_general.cc -D\'GCC_RNG='+rng+'\' -D\'GCC_PB='+pb+'\' -DGCC_POSITION_INIT=\''+pos+'\' -DGCC_VELOCITY_INIT=\''+vel+'\' -DGCC_PARTICLE=\''+particule+'\' -DGCC_W=\''+str(w_value)+'\' -DGCC_C=\''+str(c_value)+'\' -DGCC_TOPOLOGY=\''+topo+'\' -DGCC_TOPOLOGY_NAME=\''+TopologyNames[Topology.index(topo)]+'\' -DGCC_EVALUATION_MODE=\''+evalmode+'\' -DGCC_RANDOM_SHUFFLE=\''+randShuffle+'\' `pkg-config --libs --cflags popot` -O3'
                                os.system(mystr)
                                os.system('./bench_general >> ' + filename)
