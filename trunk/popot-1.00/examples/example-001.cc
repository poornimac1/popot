#include <sys/time.h>

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

// And then include our headers (these make use of the RNG_GENERATOR, so we must include 
// them after the definition of RNG_GENERATOR
#include "popot.h"

// Define our problem
typedef popot::problems::SPSO2011::F1<30> Problem;

// Define our algorithm, here SPSO 2011
//typedef popot::PSO::SPSO2006::PSO<Problem>::Type PSO;
//typedef popot::PSO::SPSO2007::PSO<Problem>::Type PSO;
typedef popot::PSO::SPSO2011::PSO<Problem>::Type PSO;

#define N_RUNS 100

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {

  RNG_GENERATOR::rng_srand();
  RNG_GENERATOR::rng_warm_up();
  
  popot::tools::FIFO<N_RUNS> fitnesses;
  popot::tools::FIFO<N_RUNS> epochs;

  double logProgressMean = 0.0;

  // To keep track of the best particle ever found
  PSO::BestType best_particle;

  // To keep track of number of failures (its definition is problem dependent)
  int nb_fails = 0;

  // To measure the execution time
  struct timeval before, after;
  double total_time;
  total_time = 0;
  
  for(int i = 0 ; i < N_RUNS ; ++i)
    {
      // Some initialization of static fields
      Problem::init();

      // Let's create our swarm
      PSO* pso = new PSO();
      
      gettimeofday(&before, NULL);

      // We now run our algorithm
      pso->run(0);

      // Get time after execution
      gettimeofday(&after, NULL);

      // Evaluation the exectuion time
      total_time += (after.tv_sec + after.tv_usec * 1E-6) - (before.tv_sec + before.tv_usec * 1E-6);

      // Some display
      std::cout << "Run " << i << " f = " << std::scientific <<  pso->getBest()->getFitness() << " with " << Problem::count << " function evaluations" << std::endl;

      // Count a fail if the fitness is larger tha
      if(Problem::has_failed(pso->getBest()->getFitness()))
	nb_fails++;

      // Let's keep track of the best fitness and number of function evaluations
      logProgressMean += -log(pso->getBest()->getFitness());
      fitnesses.insert(pso->getBest()->getFitness());
      epochs.insert(Problem::count);

      // Keep track of the best particle ever found
      if(i == 0 || (pso->getBest()->compare(&best_particle) < 0))
	best_particle = *(pso->getBest());

      // Clean up the memory before starting a new trial
      delete pso;
      Problem::free();
    }

  logProgressMean /= double(N_RUNS);
  
  std::cout << "Statistics over " << N_RUNS << " runs : " << std::endl;
  std::cout << "Eval. (mean) : " << epochs.mean() << std::endl;
  std::cout << "Error (mean) : " << fitnesses.mean() << std::endl;
  std::cout << "Std. dev. : " << fitnesses.variance() << std::endl;
  std::cout << "Log_progress (mean) : " << logProgressMean << std::endl;
  std::cout << "Failure(s) " << nb_fails << " Success rate = " << std::fixed << double(N_RUNS - nb_fails)*100.0 / N_RUNS << " %" << std::endl;
  std::cout << "Best min value : " << std::scientific << fitnesses.min() << std::endl;
  std::cout << "Position of the extremum : " << best_particle << std::endl;
  std::cout << "errMax : " << fitnesses.max() << std::endl;
  std::cout << "Mean time per run : " << total_time / N_RUNS << " s." << std::endl;

}
