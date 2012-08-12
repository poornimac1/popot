#!/usr/bin/python
from numpy import array
from numpy import ones, zeros
from numpy.random import rand
import numpy as np
from math import *

############### Parameters of the problem 

Nd = 10 # number of dimensions
CS = 10
No = 10 # Number of onlooker beers
Ns = 3 # Number of maintained solutions (employed bees)
Niter = 1000 # Number of iterations
L = (CS*Nd)/2  # Limit before switching to a new random position

def my_func(params):
    # Return the fitness to minimize
    return rastrigin(params)

def my_bounds():
    # Sets up lb and ub 
    rastrigin_bounds()
########################################


# Some global arrays
solutions = zeros((Nd, Ns))
fitnesses = zeros((Ns,1))
probabilities = zeros((Ns,1))
limits = zeros((Ns,1))

lb = zeros((Nd,))
ub = zeros((Nd,))


## Predefinition of some problems
def sphere_bounds():
    global lb, ub
    lb = -5.0 * ones((Nd,))
    ub = 5.0 * ones((Nd,))

def sphere(params):
    val = 0.0
    for i in range(Nd):
        val = val + params[i]**2 
    return val

def rastrigin(params):
    val = 0.0
    for i in range(Nd):
        val = val + params[i]*params[i] + 10.0 * (1.0 - np.cos(2.0 * np.pi * params[i]))
    return val

def rastrigin_bounds():
    global lb, ub
    lb = -5.12 * ones((Nd,))
    ub = 5.12 * ones((Nd,))



## The ABC algorithm

def initialize():
    global solutions, fitnesses
    my_bounds()
    for i in range(Ns):
        for j in range(Nd):
            solutions[j,i] = lb[j] + (ub[j] - lb[j]) * rand()
        fitnesses[i] = fitness(my_func(solutions[:,i]))
        
def fitness(fval):
    if(fval <= 0.0):
        return 1.0 + fabs(fval)
    else:
        return 1.0 / (1.0 + fval)

def fitnesses_to_probabilities():
    # We may avoid scaling and may return the sum value for probabilistic toss
    global probabilities, fitnesses
    s = fitnesses.sum()
    for i in range(Ns):
        probabilities[i] = fitnesses[i] / s

def scout_phase():
    global limits, solutions, fitnesses
    # Check if some limits are above L
    for i in range(Ns):
        if(limits[i] >= L):
            limits[i] = 0
            for j in range(Nd):
                solutions[j,i] = lb[j] + rand()*(ub[j] - lb[j])
            fitnesses[i] = fitness(my_func(solutions[:,i]))

def onlooker_phase():
    global solutions, fitnesses, probabilities
    # For each onlook bee
    # we select one source with a probabistic toss (depending on the fitnesses)
    # and look around it along a direction toward one dimension of another source
    
    for i in range(No):
        rv = rand()
        # Find the selected source
        j = 0
        while((probabilities[0:j].sum() < rv) and (j < Ns - 1)):
            j+=1
        # Select randomly a second source and dimension
        rand_index = int(Ns * rand())
        rand_dim = int(Nd * rand())
        rand_phi = -1.0 + 2.0 * rand()
        # Compute the new position
        v = np.copy(solutions[:,j])
        v[rand_dim] = v[rand_dim] + rand_phi * (v[rand_dim] - solutions[rand_dim, rand_index])
        # evaluate its fitness
        fit = fitness(my_func(v))
        # Perform a greedy selection
        if(fit > fitnesses[j]):
            fitnesses[j] = fit
            for k in range(Nd):
                solutions[k,j] = v[k]
            limits[j] = 0
        else:
            limits[j] = limits[j] + 1

def employed_phase():
    global solutions, fitnesses
    for i in range(Ns):
        rand_index = int(Ns * rand())
        rand_dim = int(Nd * rand())
        rand_phi = -1.0 + 2.0 * rand()
        # Compute the new potential solution
        v = np.copy(solutions[:,i])
        v[rand_dim] = v[rand_dim] + rand_phi * (v[rand_dim] - solutions[rand_dim, rand_index]) 
        fit = fitness(my_func(v))
        # Perform a greedy selection
        if(fit > fitnesses[i]):
            fitnesses[i] = fit
            for j in range(Nd):
                solutions[j,i] = v[j]
            limits[i] = 0
        else:
            limits[i] = limits[i] + 1
 

def main():
    # Initialization of the problem
    initialize()
    best_solution = np.copy(solutions[0])
    best_fit = fitnesses[0,0]
    for j in range(1,Ns):
        if(fitnesses[j,0] > best_fit):
            best_fit = fitnesses[j,0]
            best_solution = np.copy(solutions[:,j])
    # Main loop
    for i in range(Niter):
        employed_phase()
        fitnesses_to_probabilities()
        onlooker_phase()
        # Keep track of the best solution
        for j in range(Ns):
            if(fitnesses[j,0] > best_fit):
                best_fit = fitnesses[j,0]
                best_solution = np.copy(solutions[:,j])
        print "Trial ", i , " : " , best_fit , " ; "
        #print best_solution
        scout_phase()


main()
