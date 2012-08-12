# A tutorial on the Cross Entropy Method
# Algorithm 2.1
# An example of rare event probability estimation

import sampling
import CE
import math

gamma = 2

# Draw X1, ... X5 according to exp_prob
u = [0.250, 0.400, 0.100, 0.300, 0.200]

# Question : what is the probability that the total length of the shortest path
# from A to B is greater than gamma=2 ?


def performance(X):
    # Returns the length of the shortest path
    # of the graph fig. 1 "Tutorial on CE method"
    c1 = X[0] + X[3]
    c2 = X[0] + X[2] + X[4]
    c3 = X[1] + X[4]
    c4 = X[1] + X[2] + X[3]
    return min(c1,c2,c3,c4)

def my_sampling_func(v):
    # The weight of the connections are independently drawn from
    # an exponential distribution with means in v
    # i.e. the parameter of the exponential distributions are 1/v_i
    return [sampling.generate_sample_IT(lambda x: sampling.exp_cdf_1(x,[1.0/vi])) for vi in v]

def my_pdf(x, u):
    # The random variables being independent, the pdf of X
    # is the product of the individual pdfs
    return reduce(lambda x,y: x*y, [sampling.exp_pdf(xi, [1.0/ui]) for xi,ui in zip(x,u)])

##############################################
# MonteCarlo estimation
Nmc = 10**4

# With Monte-Carlo, we compute the mean of the indicator function
[mc_estimation_mu, mc_estimation_var] = sampling.monte_carlo(Nmc, lambda X: performance(X) >= gamma, lambda:my_sampling_func(u))

print "Monte Carlo with %g samples : P(X >= %e) = %.2e (std %.3e)" %(Nmc, gamma, mc_estimation_mu, math.sqrt(mc_estimation_var))



###############################################
# CE Method
N = 1000
Nmc_ce = 10**5
rho = 0.1

[v_ce,gamma_ce] = CE.algorithm_2_1(u, performance, my_sampling_func, my_pdf, [N, rho, gamma],0)

print "Cross Entropy : v = [%s] ; gamma = %.3f" % (' '.join('{0:.3f};'.format(k) for k in v_ce), gamma_ce)

# We now compute the mean of the indicator times the ratio likelihood
[mc_ce_estimation_mu, mc_ce_estimation_var] = sampling.monte_carlo(Nmc_ce, lambda X: (performance(X) >= gamma)*my_pdf(X,u)/my_pdf(X,v_ce), lambda:my_sampling_func(v_ce))
print "Monte Carlo with %g samples : P(X >= %e) = %.2e (std %.3e)" %(Nmc_ce, gamma, mc_ce_estimation_mu, math.sqrt(mc_ce_estimation_var))


