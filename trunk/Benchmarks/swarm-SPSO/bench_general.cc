// Example of compilation line :
// Using SPSO2011
//       C random number
//       F0 : Sphere in dimension 30
//g++ -o bench bench.cc -DGCC_RNG=CRNG -DGCC_PB='F0<30>' -DGCC_ALGO='SPSO2011::PSO<Problem>::Type' `pkg-config --libs --cflags popot` -O3

#include <sys/time.h>
#include <typeinfo>
#include "bench_general.h"
#include "pso.h"

#define STRINGIZE(x) #x
#define STRINGIZE_VALUE_OF(x) STRINGIZE(x)

// Define our problem
typedef BenchmarkProblems::GCC_PB Problem;

//*********************************************************

// Define our particle
class Particle_PSO_Params
{
public:
    static double w()  { return GCC_W;}   // Inertia parameter
    static double c1() { return GCC_C;}   // Best particle position weight
    static double c2() { return GCC_C;}   // Best swarm position weight
};

typedef swarm::particle::NonStochasticTraditionalParticle< Problem, Particle_PSO_Params > Particle;


typedef swarm::topology::GCC_TOPOLOGY Topology;

class Swarm_PSO_params
{
public:
    static int evaluation_mode() { return swarm::algorithm::GCC_EVALUATION_MODE;}
};


class Swarm_Stop_Criteria
{
public:
    static bool stop(double fitness, int epoch)
    {
      return Problem::stop(fitness,epoch);
    }
};


// Define our algorithm,
typedef swarm::algorithm::PSO<Swarm_PSO_params, Particle, Topology, Swarm_Stop_Criteria> SwarmPSO;

class PSO : public SwarmPSO
{
public:
  PSO() : SwarmPSO()
  {}

  double getBestFitness(void)
  {
    return this->getBest()->getFitness();
  }

  void getBestPosition(double * pos)
  {
  }
  
  void run(void)
  {
    while(!Swarm_Stop_Criteria::stop(getBestFitness(), getEpoch()))
      {
	this->step();
      }
  }

  int getEpoch(void)
  {
    return epoch;
  }
};


//*********************************************************

#define N_RUNS 500

#include "benchmark.h"
typedef popot::benchmark::Benchmark<PSO, Problem, N_RUNS> Benchmark;

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {
  srand(time(NULL));
 
  Benchmark bm;
  bm.run(0);

  
  std::cout << "PB=" << STRINGIZE_VALUE_OF(GCC_PB) << ";"
	    << "TOPO=" << STRINGIZE_VALUE_OF(GCC_TOPOLOGY_NAME) << ";"
	    << "evalMode=" << STRINGIZE_VALUE_OF(GCC_EVALUATION_MODE) << ";"
	    << "w=" << GCC_W << ";"
	    << "c=" << GCC_C << ";"
	    << bm; 
  
}
