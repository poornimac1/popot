// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"




int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();

  int dim = 10;
  			       
  auto b = popot::PSO::particle::base(dim,   
				      &popot::initializer::position::uniform_random,
				      [] (size_t index) -> double { return -10; },
				      [] (size_t index) -> double { return  10; });

  std::cout << b << std::endl;
  b.init();
  std::cout << b << std::endl;
}
