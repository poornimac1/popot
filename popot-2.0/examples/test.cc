// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::problems::Rosenbrock Problem;

class Params
{
public:
  static double w()  { return 1.0/(2.0*log(2.0));}   // Inertia parameter
  static double c() { return 0.5 + log(2.0);}   // Best particle position weight
};

typedef popot::PSO::particle::SPSO2007Particle< Problem, Params> Particle;

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();

  Problem prob(10);
  Particle p(prob.dimension);
  p.init(prob);
  p.confine(prob);

  std::cout << p << std::endl;
}
