// In this example, we check the validity of the random number generators

#include "rng_generators.h"


typedef popot::rng::MersenneRNG RNG_GENERATOR;
#define FILENAME "mrng.data"

/*
typedef popot::rng::CRNG RNG_GENERATOR;
#define FILENAME "crng.data"
*/
#include "popot.h"

typedef popot::problems::Rosenbrock<10> Problem;
typedef popot::PSO::particle::Particle<Problem, popot::PSO::initializer::PositionUniformRandom, popot::PSO::initializer::VelocitySPSO2011> Particle;

int main(int argc, char * argv[])
{
  // // Random number generation
  // Let's initialize the seed
  std::ofstream outfile(FILENAME);

  RNG_GENERATOR::rng_srand();
   for(int i = 0 ; i < 100 ; ++i)
     {
       outfile << i << " " << popot::math::uniform_random(0,1) << std::endl;
       //std::cout << popot::math::uniform_random(0,5) << std::endl;  
     }
   outfile.close();


   Particle p;
   p.init();

   Particle p1 = p;

   std::cout << p << std::endl;

   std::cout << p1 << std::endl;
}
