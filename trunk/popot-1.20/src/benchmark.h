#ifndef POPOT_BENCHMARK
#define POPOT_BENCHMARK

namespace popot
{
  namespace benchmark
  {
    /**
     * Benchmarking class; Given an algorithm and a number of trials for this algorithm
     * it computes the mean error, std, etc... and various statistics on the performances
     * We suppose ALGO provides :
     * - a constructor without arguments, 
     * - a step method, 
     * - a getBestFitness method,
     * - a getBestPosition method
     * - a getEpoch() 
     * The Problem is also the one used as template for ALGO; It provides the stop method and the function 
     * evaluations counter
     */
    template<typename ALGO, typename PROBLEM, int NB_TRIALS>
    class Benchmark
    {
    private:
      double _mean_error;
      double _std_error;
      double _mean_nb_f_evaluations;
      double _log_progress;
      int    _nb_failures;
      double * _best_position;
      double _best_fitness;

    public:
      Benchmark() : _mean_error(0.0),_std_error(0.0), _mean_nb_f_evaluations(0.0),
		    _log_progress(0.0), _nb_failures(0), _best_position(0), _best_fitness(0)
      {
	_best_position = new double[PROBLEM::nb_parameters];
      }

      ~Benchmark()
      {
	delete[] _best_position;
      }

      /**
       * Returns the mean error over NB_TRIALS
       */
      double getMeanError(void) const
      {
	return _mean;
      }

      /**
       * Returns the std error over NB_TRIALS
       */
      double getStdError(void) const
      {
	return _std_error;
      }

      /**
       * Returns the mean number of function evaluations to meet the stop criteria over NB_TRIALS
       */
      double getMeanFEvaluations(void) const
      {
	return _mean_nb_f_evaluations;
      }

      /**
       * Returns the log progress over NB_TRIALS
       * NOT FUNCTIONAL YET
       */
      double getLogProgress(void) const
      {
	return _log_progress;
      } 

      /**
       * Returns the success rate over NB_TRIALS
       */
      double getSuccessRate(void) const
      {
	return double(NB_TRIALS-_nb_failures)/double(NB_TRIALS);
      } 

      /**
       * Returns the best position over NB_TRIALS
       */
      void getBestPosition(double * pos) const
      {
	memcpy(pos, _best_position, PROBLEM::nb_parameters*sizeof(double));
      } 

      /**
       * Returns the best fitness over NB_TRIALS
       */
      double getBestFitness(void) const
      {
	return _best_fitness;
      }
 
      /**
       * Runs the benchmarks and computes the statistics
       */
      void run(int mode=0)
      {
	ALGO* algo;
	for(int i = 0 ; i < NB_TRIALS ; ++i)
	  {
	    PROBLEM::init();
	    
	    algo = new ALGO();
	    algo.run();

	    // Get the best 


	    PROBLEM::free();
	    delete algo;

	    // Eventually print the statistics
	    if(mode)
	      print(i);
	  }
      }

      /**
       * Print the statistics
       */
      void print(int trial_nb)
      {
	std::cout << "Trial " << trial_nb;

	std::cout << std::endl;
      }
    };
  } // benchmark
} // popot

#endif
