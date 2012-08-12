# A tutorial on the Cross Entropy Method
# Algorithm 2.2
# An example of optimization problem
# Find the binary vector that maximizes performance_func

import sampling
import CE
import random
import math

Np = 10
N = 50
rho = 0.1

y = [(random.random() < 0.5) for i in range(Np)]
p0 = [1.0/2.0 for i in range(Np)]

def sampling_func(p):
    return [(random.random()<pi) for pi in p]

def performance_func(x):
    global y
    return Np - sum([math.fabs(xi-yi) for xi,yi in zip(x,y)])

pt= CE.algorithm_2_2(p0, performance_func, sampling_func, [N, rho], 0)
print "Optimal value : [%s] \n" % (' '.join('{0:.1f};'.format(k) for k in y))
print "Found         : [%s] \n" % (' '.join('{0:.1f};'.format(k) for k in pt))
print "Performance : %.2f [optimal : %.2f]" % (performance_func(pt), Np)
