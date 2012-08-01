// In this example, we completely define a swarm
// made of Modified Barebone particles

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::JKissRNG RNG_GENERATOR;

// And then include our headers (these make use of the RNG_GENERATOR, so we must include 
// them after the definition of RNG_GENERATOR
#include "popot.h"

// Define our problem
// here Rosenbrock in 30 dimensions (see problems.h)
typedef popot::problems::SPSO2011::F2<30> Problem;

// Define our initializers for the position and velocity
typedef popot::PSO::initializer::PositionUniformRandom PositionInit;
typedef popot::PSO::initializer::VelocityZero VelocityInit;

// Let's define our particle
typedef popot::PSO::particle::ModifiedBareboneParticle<Problem, PositionInit, VelocityInit> Particle;


// The topology
// Let's use 25 particles 
typedef popot::PSO::topology::AdaptiveRandom<25, 3, Particle> Topology;
//typedef popot::PSO::topology::Full<25, Particle> Topology;
//typedef popot::PSO::topology::Ring<25, Particle> Topology;
//typedef popot::PSO::topology::VonNeuman<5, 5, Particle> Topology;

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

  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  // To keep track of the best particle ever found
  PSO::BestType best_particle;

  Problem::init();
  
  // Let's create our swarm
  PSO pso;
  
  // We now run our algorithm
  pso.run(1);
  
  // Some display
  std::cout << "Best particle : " << *(pso.getBest()) << std::endl;

  Problem::free();
}
