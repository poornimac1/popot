#include <sys/time.h>
#include <typeinfo>
#include "bench.h"

// Example of compilation line :
// Using SPSO2011
//       C random number
//       F0 : Sphere in dimension 30
//g++ -o bench bench.cc -D'GCC_RNG=popot::rng::CRNG' -D'GCC_PB=BenchmarkProblems::F0<30>' -D'GCC_ALGO=popot::PSO::SPSO2011::PSO<Problem>::Type' `pkg-config --libs --cflags popot` -O3


// We first define the generator of random numbers
#include "rng_generators.h"
typedef GCC_RNG RNG_GENERATOR;

// And then include our headers (these make use of the RNG_GENERATOR, so we must include 
// them after the definition of RNG_GENERATOR
#include "popot.h"

// Define our problem
typedef GCC_PB Problem;

// Define our algorithm,
typedef GCC_ALGO PSO;

#define N_RUNS 1000

typedef popot::benchmark::Benchmark<PSO, Problem, N_RUNS> Benchmark;

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {
  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  Benchmark bm;
  bm.run(1);
}
