// In this example, we check the validity of the random number generators

#include "rng_generators.h"
typedef popot::rng::JKissRNG RNG_GENERATOR;

#include "popot.h"

int main(int argc, char * argv[])
{
  // // Random number generation
  // Let's initialize the seed
  RNG_GENERATOR::rng_srand();
   for(int i = 0 ; i < 100 ; ++i)
     std::cout << popot::math::uniform_random(0,1) << std::endl;  
}
