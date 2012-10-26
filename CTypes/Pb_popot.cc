// g++ -shared -O3 -o Pb_popot.so Pb_popot.cc `pkg-config --libs --cflags popot`

#include <math.h>
#include <cstdlib>

#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include <popot.h>
typedef popot::problems::Ackley<10> Problem;

// Then define the bindings for python
extern "C"
{
  int get_dimension(void) { return Problem::nb_parameters;}
  void init_problem(void) { Problem::init();}
  void free_problem(void) { Problem::free();}
  double get_lbound(int index) { return Problem::get_lbound(index);}
  double get_ubound(int index) { return Problem::get_lbound(index);}
  bool stop(double fitness, int epoch) { return Problem::stop(fitness, epoch);}
  double evaluate(double * params) {return Problem::evaluate(params); }

  int get_nb_f_evaluations(void) { return Problem::count;}

}
