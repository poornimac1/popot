#include <sys/time.h>

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

// And then include our headers (these make use of the RNG_GENERATOR, so we must include 
// them after the definition of RNG_GENERATOR
#include "popot.h"

// Define our problem
typedef popot::problems::Rosenbrock<30> Problem;

// Define our algorithm, here SPSO 2011
//typedef popot::PSO::SPSO2006::PSO<Problem>::Type PSO;
//typedef popot::PSO::SPSO2007::PSO<Problem>::Type PSO;
typedef popot::PSO::SPSO2011::PSO<Problem>::Type PSO;

#define N_RUNS 1000

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {

  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  // Some initialization of static fields
  Problem::init();
  
  // Let's create our swarm
  PSO* pso = new PSO();
  
  // We now run our algorithm
  pso->run(1);

  std::cout << "Nb steps : " << pso->epoch << " ; nb neigh : " << pso->nb_new_neigh << std::endl;

  delete pso;
  Problem::free();

}
