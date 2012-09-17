// Example for comparing to the C code of Clerc et al.

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

// Define our problem
typedef popot::problems::Rosenbrock<5> Problem;

// Define our initializers for the position and velocity
typedef popot::PSO::initializer::PositionUniformRandom PositionInit;
typedef popot::PSO::initializer::VelocityZero VelocityInit;

class ParticleParams
{
public:
  static double w() { return 1. / (2 * log ((double) 2));  }
  static double c() { return 0.5 + log ((double) 2); }
};

// Let's define our particle
typedef popot::PSO::particle::BenchSPSO2006Particle<Problem, ParticleParams> Particle;


// The topology
typedef popot::PSO::topology::AdaptiveRandom<5, 3, Particle> Topology;

// For the algorithm type, we need to mention
// if we use synchronous or asynchronous evaluation
// as well as the criteria to stop iterating
class PSO_Params
{
public:
  static int evaluation_mode() {return popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION;}
};

// For the stopping criteria, we use the one provided by the problem
// Be carefull, the Problem keeps track internally of the number of Function evaluations and this criteria is included in the condition to stop
class StopCriteria
{
public:
  static bool stop(double fitness, int epoch)
  {
    return Problem::stop(fitness, epoch);
  }
};


// We can now define our algorithm
typedef popot::PSO::algorithm::Base<PSO_Params, Particle, Topology, StopCriteria> PSO;

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {

  RNG_GENERATOR::rng_srand(0);
  RNG_GENERATOR::rng_warm_up();
  
  // To keep track of the best particle ever found
  PSO::BestType best_particle;

  Problem::init();
  
  // Let's create our swarm



  std::cout << "Random before creating : " << popot::math::uniform_random(0,1) << std::endl;
  RNG_GENERATOR::print();

  PSO pso;
  pso.print(0);


  // We now run our algorithm
  for(int i = 0 ; i< 2; ++i)
    {
      std::cout << "########## Step " << i << "##############" << std::endl;
      pso.step();
      pso.print(0);
    }
  
  // Some display
  std::cout << "Best particle : " << *(pso.getBest()) << std::endl;

  Problem::free();
}
