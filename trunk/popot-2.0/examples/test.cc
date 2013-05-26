// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"


double evaluate(double * x)
{
  return -1.0;
}

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();

  int dim = 10;
  
    /*

  auto b = popot::PSO::particle::make_base(dim,   
					  &popot::initializer::position::uniform_random,
					  [] (size_t index) -> double { return -10; },
					  [] (size_t index) -> double { return  10; },
					  &evaluate);
  */
  /*			       
  auto b = popot::PSO::particle::make_particle(dim,   
					  &popot::initializer::position::uniform_random,
					  &popot::initializer::velocity::half_diff,
					  [] (size_t index) -> double { return -10; },
					  [] (size_t index) -> double { return  10; },
					  &evaluate);
  */
  popot::PSO::particle::make_spso2006_particle(dim,
					  [] (size_t index) -> double { return -10; },
					  [] (size_t index) -> double { return  10; },
							&evaluate, 0.7, 1.2);
  /*
  auto c = b;
  c.setPosition(2,2);

  std::cout << b << std::endl;
  b.init();
  std::cout << b << std::endl;
  std::cout << c << std::endl;
  */
}
