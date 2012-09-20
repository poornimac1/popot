// Example of compilation line :
// Using SPSO2011
//       C random number
//       F0 : Sphere in dimension 30
//g++ -o bench bench.cc -DGCC_RNG=CRNG -DGCC_PB='F0<30>' -DGCC_ALGO='SPSO2011::PSO<Problem>::Type' `pkg-config --libs --cflags popot` -O3

#include <sys/time.h>
#include <typeinfo>
#include "bench.h"

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::GCC_RNG RNG_GENERATOR;

// And then include our headers (these make use of the RNG_GENERATOR, so we must include 
// them after the definition of RNG_GENERATOR
#include "popot.h"

// Define our problem
typedef BenchmarkProblems::GCC_PB Problem;

// Define our algorithm,
typedef popot::PSO::GCC_ALGO PSO;

#define N_RUNS 500

typedef popot::benchmark::Benchmark<PSO, Problem, N_RUNS> Benchmark;

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {
  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  Benchmark bm;
  bm.run(0);

  //
  
  std::cout << "RNG=" << STRINGIZE_VALUE_OF(GCC_RNG) << ";" 
	    << "PB=" << STRINGIZE_VALUE_OF(GCC_PB) << ";"
	    << "ALGO=" << GCC_ALGO_NAME << ";"
    	    << bm ;
  

}
