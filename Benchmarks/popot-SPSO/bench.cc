// Example of compilation line :
// Using SPSO2011
//       C random number
//       F0 : Sphere in dimension 30
//g++ -o bench bench.cc -DGCC_RNG=CRNG -DGCC_PB=F0 -DGCC_DIM=10 -DGCC_ALGO=spso2011 `pkg-config --libs --cflags popot` -O3 -std=c++0x

#include <sys/time.h>
#include <typeinfo>


// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::GCC_RNG RNG_GENERATOR;

// And then include our headers (these make use of the RNG_GENERATOR, so we must include 
// them after the definition of RNG_GENERATOR
#include "popot.h"
#include "bench.h"
#include "benchmark.h"

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

// Define our problem
typedef BenchmarkProblems::GCC_PB Problem;
// F0, F1, F2, ...

// Define our algorithm,
// GCC_ALGO
//popot::algorithm::spso2011


#define N_RUNS 500

// **************************************** //
// ************** Main ******************** //
// **************************************** //

// difficulté :
// j'aimerais fournir uen référence vers la fonction créant l'algo
// et une référence vers un pb
// mais le compilo n'arrive pas à inférer le type de la fonction de création de l'algo
// parce que son type dépend du pb, et donc sera déterminé quand on l'aura appelé.

int main(int argc, char* argv[]) {
  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  Problem pb(GCC_DIM);
  auto algo = popot::algorithm::GCC_ALGO(pb._dimension,
  					 [&pb] (size_t index) -> double { return pb.get_lbound(index); },
  					 [&pb] (size_t index) -> double { return pb.get_ubound(index); },
  					 [&pb] (double fitness, int epoch) -> bool { return pb.stop(fitness, epoch);},
  					 [&pb] (TVector &pos) -> double { return pb.evaluate(pos);}
					);

  auto bm = popot::benchmark::make_benchmark(algo, pb, N_RUNS);
  bm.run();

  //
  
  std::cout << "RNG=" << STRINGIZE_VALUE_OF(GCC_RNG) << ";" 
	    << "PB=" << STRINGIZE_VALUE_OF(GCC_PB) << ";"
	    << "DIM=" << STRINGIZE_VALUE_OF(GCC_DIM) << ";"
	    << "ALGO=" << STRINGIZE_VALUE_OF(GCC_ALGO) << ";"
	    << bm ;
  

}
