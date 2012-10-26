
from ctypes import *
import numpy
pb_library = CDLL("Pb.so") 

get_size_func=pb_library.get_size
get_size_func.restype =c_int

evaluate_func = pb_library.evaluate
evaluate_func.restype =c_double

# Create the array of parameters
size = get_size_func()
params=numpy.zeros((size,),dtype=float)
for i in range(size):
    params[i] = 1;


print("%f"%(pb_library.evaluate(params.ctypes.data_as(c_void_p))))
