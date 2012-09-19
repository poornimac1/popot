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
      double _success_rate;

    public:
      Benchmark() : _mean_error(0.0),_std_error(0.0), _mean_nb_f_evaluations(0.0),
		    _log_progress(0.0), _nb_failures(0), _best_position(0), _best_fitness(0), _success_rate(0.0)
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
	return _mean_error;
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
	return _success_rate;
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
	// Some temporary variables used to compute online the mean, std, .. 
	double sum_error = 0;
	double sum_square_errors=0;
	double sum_fe = 0;
	 
	for(int i = 1 ; i <= NB_TRIALS ; ++i)
	  {
	    PROBLEM::init();
	    
	    algo = new ALGO();
	    algo->run();

	    // Update the statistics
	    sum_error += algo->getBestFitness();
	    sum_square_errors += algo->getBestFitness()*algo->getBestFitness();
	    sum_fe += PROBLEM::count;
	    if(PROBLEM::has_failed(algo->getBestFitness()))
	      _nb_failures++;
	    if(i == 1 || algo->getBestFitness() < _best_fitness)
	      {
		_best_fitness = algo->getBestFitness();
		algo->getBestPosition(_best_position);
	      }
	      
	    // Eventually print the statistics
	    if(mode)
	      {
		_mean_error = sum_error/double(i);
		_std_error = sqrt((sum_square_errors - 2.0 * _mean_error * sum_error + double(i) * _mean_error*_mean_error)/double(i));
		_mean_nb_f_evaluations = sum_fe/double(i);
		_log_progress = 0.0;
		_success_rate = 100.*double(i-_nb_failures)/double(i);

		print(i);
	      }
	    
	    // Free the memory
	    PROBLEM::free();
	    delete algo;
	  }
		     
	// Before leaving, ensure that the statistics are up-to-date;
	_mean_error = sum_error/double(NB_TRIALS);
	_std_error = sqrt((sum_square_errors - 2.0 * _mean_error * sum_error + double(NB_TRIALS) * _mean_error*_mean_error)/double(NB_TRIALS));
	_mean_nb_f_evaluations = sum_fe/double(NB_TRIALS);
	_log_progress = 0.0;
	print(NB_TRIALS);
      }

      /**
       * Print the statistics
       */
      void print(int trial_nb)
      {
	std::cout << "Trial " << trial_nb << ";";
	std::cout << "Error (mean) = " << getMeanError() << " ";
	std::cout << "Error (std) = " << _std_error << " ";
	std::cout << "FE (mean) = " << _mean_nb_f_evaluations << " ";
	std::cout << "Log progress (mean) = " << getLogProgress() << " ";
	std::cout << "Best min = " << getBestFitness() << " ";
	std::cout << "Success rate = " << getSuccessRate() << "% ";
	std::cout << std::endl;
      }
    };
  } // benchmark
} // popot

#endif