// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::problems::Rosenbrock Problem;
typedef popot::PSO::initializer::PositionUniformRandom PositionInit;


typedef popot::PSO::particle::BaseParticle< Problem, PositionInit> Particle;

int main(int argc, char * argv[])
{
  Problem prob(10);
  Particle p(prob.dimension);
  p.init(prob);
}
