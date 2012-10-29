// g++ -shared -O3 -o Problem.so Problem.cc -lboost_python -L/usr/lib `pkg-config --libs --cflags popot python` -I/usr/include


#include <math.h>
#include <cstdlib>

#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include <popot.h>
typedef popot::problems::Ackley<10> Problem;

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(Problem)
{
  class_<Problem>("Problem")
    .def("get_lbound",&Problem::get_lbound).staticmethod("get_lbound")
    .def("get_ubound",&Problem::get_ubound).staticmethod("get_ubound")
    .def("evaluate",&Problem::evaluate).staticmethod("evaluate")
    ;
}
