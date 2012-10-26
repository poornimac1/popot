
from ctypes import *
import numpy
pb_library = CDLL("Pb_popot.so") 

# Defines all the required bindings to the C++ Problem
get_dimension_func=pb_library.get_dimension
get_dimension_func.restype =c_int

init_problem_func = pb_library.init_problem
free_problem_func = pb_library.free_problem

get_lbound_func = pb_library.get_lbound
get_lbound_func.restype =c_double

get_ubound_func = pb_library.get_ubound
get_ubound_func.restype =c_double

stop_func = pb_library.stop
stop_func.restype =c_bool

evaluate_func = pb_library.evaluate
evaluate_func.restype =c_double

get_nb_f_evaluations_func = pb_library.get_nb_f_evaluations
get_nb_f_evaluations_func.restype = c_int

# # Create the array of parameters
# size = get_size_func()
# params=numpy.zeros((size,),dtype=float)

# def myfunc(x):
#     params=numpy.zeros((size,),dtype=float)
#     for i in range(size):
#         params[i] = x[i]
#     return pb_library.evaluate(params.ctypes.data_as(c_void_p))


# ###### First possibility : probably better with noisy fitness
# import cma, numpy as np

# print cma.Options

# func = myfunc
# es = cma.CMAEvolutionStrategy(np.ones(size), 1)
# logger = cma.CMADataLogger().register(es)
# nh = cma.NoiseHandler(es.N, maxevals=[1, 30])
# while not es.stop():
#     X, fit = es.ask_and_eval(func, evaluations=nh.evaluations) 
#     es.tell(X, fit)  # prepare for next iteration
#     es.sigma *= nh(X, fit, func, es.ask)  # see method __call__
#     es.countevals += nh.evaluations_just_done  # this is a hack, not important though
#     logger.add(more_data = [nh.evaluations, nh.noiseS])  # add a data point 
#     es.disp()
#     # nh.maxevals = ...  it might be useful to start with smaller values and then increase
# print(es.stop())
# print(es.result()[-2])  # take mean value, the best solution is totally off
# #assert sum(es.result()[-2]**2) < 1e-9
# print(X[np.argmin(fit)])  # not bad, but probably worse than the mean
# logger.plot()


# # Second possibility
# #res = cma.fmin(myfunc, size * [1], 1)
# #res[-1].plot()
