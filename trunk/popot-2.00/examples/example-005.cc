// Example for comparing to the C code of Clerc et al.

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;
//typedef popot::rng::MersenneRNG RNG_GENERATOR;

#include "popot.h"

// Define our problem
typedef popot::problems::SPSO2011Bench::Rosenbrock<5> Problem;

class ParticleParams
{
public:
  static double w() { return 1. / (2 * log ((double) 2));  }
  static double c() { return 0.5 + log ((double) 2); }
};

// Let's define our particle
//typedef popot::PSO::particle::SPSO2011Particle<Problem, ParticleParams> Particle;
typedef popot::PSO::particle::ModifiedBareboneParticle<Problem> Particle;
 

// The topology
typedef popot::PSO::topology::VonNeuman<25, Particle> Topology;
//typedef popot::PSO::topology::AdaptiveRandom<5, 3, Particle> Topology;

// For the algorithm type, we need to mention
// if we use synchronous or asynchronous evaluation
// as well as the criteria to stop iterating
class PSO_Params
{
public:
  static bool random_shuffle() { return false;}
  static int evaluation_mode() { return popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION;}
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

//typedef popot::PSO::SPSO2006::PSO<Problem>::Type PSO;
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
  //pso.print(0);

  std::ofstream outfile("f_trace.data");

  // We now run our algorithm
  double delta;  
  for(int i = 0 ; i<= 100; ++i)
    {
      delta = pso.step();
      outfile << i << " " << pso.nb_new_neigh << " " << pso.getBest().getFitness() << " " << std::scientific << delta << std::endl;
    }
  std::cout << RNG_GENERATOR::nb_calls << " random numbers used " << std::endl;
  outfile.close();
  
  
  // Some display
  std::cout << "Best particle : " << pso.getBest() << std::endl;


  std::cout << "Nb of function evaluations : " << Problem::count << std::endl;
  std::cout << "Nb of new neighbors : " << pso.nb_new_neigh << std::endl;

  Problem::free();
}
