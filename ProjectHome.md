### Downloads ###

The latest version is [popot-2.12](http://popot.googlecode.com/svn/trunk/Downloads/popot-2.12.tar.gz)



### News ###

popot 2.12:
adding the possibility to load/save particles , check out example-001.cc

popot 2.11:
minor bug correction in pso algorithm .init() (no init of epoch and nb\_new\_neigh)

popot 2.10:
- a python wrapper is provided thanks to H. Glaude. An example of how to use it is provided in examples/example-000.py

popot 2.00:

- major modifications of the interface to avoid having to define all the elements of the problem at compile time; Now, the problem's elements (e.g. its dimension) can be set at execution time. You may also use lambda functions to quickly define your problem, check out the examples.

popot 1.21:

- Some minor bug corrections (unsigned int/int ; return for BaseParticle::operator=)

- a Benchmark class to automatically evaluate the performances of an algorithm

popot 1.10 :

- Added StochasticSPSO2011, StochasticSPSO2006

- The local best within the neighborhoods are recomputed at each step

- The best of the best is now defined as the best of the personal best

- There are some slight differences between the Standards coded in popot and the SPSO provided on swarm central that are being investigated; these are apparently due to numerical roundings (say 1e-13 considered as 0 in one code, not in the other, especially when used to decide if we generate new neighborhoods)


#### Description ####

popot is a C++ Template based library for population based optimization. It aims at providing within the same library various black box optimization algorithms relying on interacting populations such particle swarm, ant colonies, evolutionary strategies, etc...

### Algorithms available ###

So far, the library contains two classes of algorithms :

- Particle Swarm Optimization

- Artificial Bee Colony

(I hope to include also CMAES, Cross Entropy, and other population based metaheuristics. Actually there are tons of swarm based optimization algorithms, just check out for example "A Brief Review of Nature-Inspired Algorithms for Optimization", Fister 2013)


### Documentation ###

The Doxygen documentation is available there : [Doxygen documentation](http://popot.googlecode.com/svn/doc/html/index.html)
(included and readable in the SVN by [fixing the mime-types](http://manjeetdahiya.com/2010/09/29/serving-html-documentation-from-google-code-svn/))


### Particle Swarm Optimization ###


The template based approach actually allows to define various algorithms. There are the standard PSO 2006, 2007 and 2011 implemented but you can also play by combining :

- Different topologies : Full, Ring, VonNeumann, AdaptiveRandom

- Different initialization schemes for the position and velocity : null, half diff, uniform

- Different position or velocity update rules : SPSO2006 , SPSO2011, Barebone and Modified Barebone (more to come : HPSO, FIPS, ...)

- Different confinment methods : position clamping and null velocity, position clamping and inverse velocity, etc..

- Different random number generators (C random number generator, JKiss, Kiss (the implementation only works for 32 bits architecture))

The implementation of the standard PSO is based on the description of M. Clerc available at  : http://clerc.maurice.free.fr/pso/SPSO_descriptions.pdf

The implementation within POPOT of SPSO 2006 and SPSO 2011 is compared to the codes available online at http://www.particleswarm.info/Programs.html on the following functions :

F0 : Sphere in dimension 30, [-100, 100], 75000 Function evaluations, min error of 1e-2

F1 : Griewank in dimension 30, [-600, 600], 75000 Function evaluations, min error of 0.05

F2 : Rosenbrock in dimension 30, [-30, 30] , 75000 function evaluations, min error of 100

F3 : Rastrigin in dimension 30, [-5.12, 5.12], 75000 Function evaluations, min error of 50

F4 : Tripod function (in 2D), [-100, 100] , 100000 Function evaluations, min error of 1e-4

F5 : Ackley in dimension 30 , [-30, 30] , 80000 function evaluations, min error of 0 (therefore run over the 80000 FE)

F7 : Schwefel 1.2 in dimension 40, [-100, 100], 40000 function evaluations, min error of 0

F8 : Schwefel 2.22 in dimension 30, [-10, 10] , 100000 function evaluations, min error of 1e-4

F9 : Neumaier 3 in dimension 40, [-40\*40 ; 40\*40] , 40000 function evaluations, min error of 0


The SPSO 2011 that is benchmarked uses the C random number generator and a uniform initialization of the positions (not the KissRNG nor the Halton sequence). The KissRNG and Halton sequence will be benchmarked soon. These benchmarks were done with the version 2.01 ;

![http://popot.googlecode.com/svn/wiki/performances.png](http://popot.googlecode.com/svn/wiki/performances.png)

### Illustrations in 2D of PSO ###

Below you see an example of SPSO-2011 behavior while optimizing a 2D rastrigin, with the associated best fitness (in log scale for the y-axis) and a messy display of the connections (updated each time the global best is not improved)

<img src='http://popot.googlecode.com/svn/trunk/2DIllustrations/example-2D.gif'>
<img src='http://popot.googlecode.com/svn/trunk/2DIllustrations/example-2D-fitness.png'>
<img src='http://popot.googlecode.com/svn/trunk/2DIllustrations/example-2D-connections.gif'>

<h3>Stochastic PSO</h3>

Below you see an example of StochasticSPSO-2011 behavior while optimizing a 2D gaussian fitness. The center of the gaussian is moving continuously and sometomes abruptly.<br>
<br>
<img src='http://popot.googlecode.com/svn/trunk/2DIllustrations/example-moving.gif'>

<h3>Artificial Bee Colony</h3>

See example-004 for its use;<br>
<br>
<br>
<h2>TODO</h2>

- benchmark the different PSO update rules and the vast amount of variants (RNG, confinment, update rules, meta parameters, population size, topology, ...)<br>
<br>
- bring a class for easily benchmarking an algorithm on a problem (with mean error, std, log progress, etc..)<br>
<br>
- add more PSO versions : FIPS, HPSO, ..<br>
<br>
- extend to other population based heuristics (CMAES, Cross Entropy, Natural Evolutionary Strategy, ...)<br>
<br>
<h2>External resources</h2>

- PyGMO (Python Parallel Global Multiobjective Optimizer) <a href='http://pagmo.sourceforge.net/pygmo/index.html'>here</a>