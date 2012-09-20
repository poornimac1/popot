// Example of compilation line :
// Using SPSO2011
//       C random number
//       F0 : Sphere in dimension 30
//g++ -o bench bench.cc -DGCC_RNG=CRNG -DGCC_PB='F0<30>' -DGCC_ALGO='SPSO2011::PSO<Problem>::Type' `pkg-config --libs --cflags popot` -O3

#include <sys/time.h>
#include <typeinfo>
#include "bench.h"

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::GCC_RNG RNG_GENERATOR;

// And then include our headers (these make use of the RNG_GENERATOR, so we must include 
// them after the definition of RNG_GENERATOR
#include "popot.h"

// Define our problem
typedef BenchmarkProblems::GCC_PB Problem;

//*********************************************************

// Define our initialization methods
typedef popot::PSO::initializer::GCC_POSITION_INIT PositionInitializer;
typedef popot::PSO::initializer::GCC_VELOCITY_INIT VelocityInitializer;


// Define our particle

class Particle_Params
{
public:
  static double w()  { return GCC_W;}   // Inertia parameter
  static double c() { return GCC_C;}   // Best particle position weight
};

typedef popot::PSO::particle::GCC_PARTICLE<Problem, Particle_Params, PositionInitializer, VelocityInitializer> Particle;

// Define our topology
typedef popot::PSO::topology::GCC_TOPOLOGY Topology;

class PSO_params
{
public:
  static bool random_shuffle() { return GCC_RANDOM_SHUFFLE;}
  static int evaluation_mode() { return popot::PSO::algorithm::GCC_EVALUATION_MODE;}
};
  
class Stop_Criteria
{
public:
  static bool stop(double fitness, int epoch)
  {
    return  Problem::stop(fitness,epoch);
  }
};

// Define our algorithm,
typedef popot::PSO::algorithm::Base<PSO_params, Particle, Topology, Stop_Criteria > PSO;



//*********************************************************

#define N_RUNS 100

typedef popot::benchmark::Benchmark<PSO, Problem, N_RUNS> Benchmark;

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {
  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  Benchmark bm;
  bm.run(0);

  //
  
  std::cout << "RNG=" << STRINGIZE_VALUE_OF(GCC_RNG) << ";" 
	    << "PB=" << STRINGIZE_VALUE_OF(GCC_PB) << ";"
	    << "PARTICULE=" << STRINGIZE_VALUE_OF(GCC_PARTICLE) << ";"
	    << "POS=" << STRINGIZE_VALUE_OF(GCC_POSITION_INIT) << ";"
	    << "VEL=" << STRINGIZE_VALUE_OF(GCC_VELOCITY_INIT) << ";"
	    << "TOPO=" << STRINGIZE_VALUE_OF(GCC_TOPOLOGY_NAME) << ";"
	    << "evalMode=" << STRINGIZE_VALUE_OF(GCC_EVALUATION_MODE) << ";"
	    << "randShuff=" << STRINGIZE_VALUE_OF(GCC_RANDOM_SHUFFLE) << ";"
    	    << bm ;
  
  
}
