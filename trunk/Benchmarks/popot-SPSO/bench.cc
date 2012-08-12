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

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {
  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  popot::tools::FIFO<N_RUNS> fitnesses;
  popot::tools::FIFO<N_RUNS> epochs;

  double logProgressMean = 0.0;

  // To keep track of the best particle ever found
  PSO::BestType best_particle;

  // To keep track of number of failures (its definition is problem dependent)
  int nb_fails = 0;

  for(int i = 0 ; i < N_RUNS ; ++i)
    {
      // Some initialization of static fields
      Problem::init();
      popot::rng::Halton<Problem::nb_parameters>::init();

      // Let's create our swarm
      PSO* pso = new PSO();
      
      // Run the PSO
      pso->run(0);

      // Some display
      //std::cout << "Run " << i << " f = " << std::scientific <<  pso->getBest()->getFitness() << " with " << Problem::count << " function evaluations" << std::endl;

      // Count a fail if the fitness is larger tha
      if(Problem::has_failed(pso->getBest()->getFitness()))
	nb_fails++;

      // Let's keep track of the best fitness and number of function evaluations
      logProgressMean += -log(pso->getBest()->getFitness());
      fitnesses.insert(pso->getBest()->getFitness());
      epochs.insert(Problem::count);

      //std::cout << i << " " << Problem::count << " " << pso->getBest()->getFitness() << std::endl;

      // Keep track of the best particle ever found
      if(i == 0 || (pso->getBest()->compare(&best_particle) < 0))
	best_particle = *(pso->getBest());

      // Clean up the memory before starting a new trial
      delete pso;
      Problem::free();
      popot::rng::Halton<Problem::nb_parameters>::free();
    }

  logProgressMean /= double(N_RUNS);
  
  // Save the statistics in a file
  std::cout << GCC_PB::name() << "\t" 
	    << N_RUNS << "\t" 
	    << epochs.mean() << "\t" 
	    << fitnesses.mean() << "\t" 
	    << fitnesses.variance() << "\t" 
	    << logProgressMean << "\t" 
	    << nb_fails << "\t" 
	    << (1.0 - nb_fails/double(N_RUNS))*100.0 << "\t"
	    << fitnesses.min() << "\t"
	    << std::endl;

}
