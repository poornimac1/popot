Standard PSO 2011, from the Particle Swarm Central http://particleswarm.info
(please, see the ReadMe.txt of the original version for more details).

Maurice.Clerc@WriteMe.com

I have added some options
+ some options. In particular
List Based Optimiser 
Mersenne RNG
Quasi-random numbers


-------
Updates:
2012-09-12	Fixed an inconsistency pointed out by Bernard Leroy.
						The suggested distribution in the hypersphere (param.distrib=0)
						is NOT uniform, as the comments said wrongly.
						The really uniform distribution has been added as an option 
						(param.distrib=-1). See alea_sphere(). Note that it works quite badly.
2012-02-01  Added a Mersenne 64-bits RNG (suggestion of Pascal Mercier)
            because KISS doesn't work on some machines (generates a negative number)
2012-02-01  Added more RNGs (BW[2] option). In particular, it is possible to read the "random"
            numbers as sequences of bits on a file, or directly as real numbers.
            The byte size is also a parameter. Meta-optimisation is possible
            to build a good list of such "random" numbers, for a given set
            of problems (VERY time consuming, though). 
2012-01-25  Frequency modulation sound parameter identification (code 28)



