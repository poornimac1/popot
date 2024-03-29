#ifndef POPOT_MATH_H
#define POPOT_MATH_H

#include <math.h>
#include <stdlib.h>

//Warning:
// Before inclusion of this file, it is mandatory to use define the random number generator
// e.g.  #define popot::rng::KissRNG RNG_GENERATOR

namespace popot
{
  namespace math
  {

    /**
     * Returns a random value in [min; max]
     * uniformely sampled
     */
    double uniform_random(double min, double max)
    {
      return min + (max - min) * double(RNG_GENERATOR::rng_rand()) / RNG_GENERATOR::RNG_RAND_MAX;
    }

    /**
     * Returns a random value in [min; max[
     * uniformely sampled
     */
    double uniform_random_exclusive(double min, double max)
    {
      return min + (max - min) * double(RNG_GENERATOR::rng_rand()) / (RNG_GENERATOR::RNG_RAND_MAX+1);
    }

    /**
     * Returns a value from a normal distribution
     */
    double normal(double mean, double std)
    {
      // We generate a sample following a normal distribution
      // with the Box-Muller method
      // u and v are randomly generated in ] 0 ; 1]
      
      double u = 1.0 - uniform_random_exclusive(0,1);
      double v = 1.0 - uniform_random_exclusive(0,1);
      double x = sqrt(-2.0 * log(u)) * cos(2.0 * M_PI * v);
      
      return mean + std * x;
    }

    /**
     * Returns a random value in [|min; max|]
     * uniformely sampled
     */
    int uniform_integer(int imin, int imax)
    {
      return floor(uniform_random(double(imin), double(imax+1)));
    }

    /**
     * Returns a random index from a discrete set of probabilities
     * probabilities are supposed to sum to 1.0 !
     */
    int random_from_array(double* proba)
    {
      double val = uniform_random(0.0, 1.0);
      int i = 0;
      double sum = proba[0];
      while(sum < val)
	{
	  i++;
	  sum += proba[i];
	}
      return i;
    }

    /**
     * Returns a permutation of indexes [|0 .. max|]
     * It supposes indexes to be allocated to size max + 1
     **/
    void random_shuffle_indexes(int * indexes, int S)
    {
      int * index_temp = new int[S];
      for(int i = 0 ; i < S ; i++)
	index_temp[i] = i;

      int rank;
      int length = S;
      int i ,t;
      for (i=0;i < S;i++)
	{
	  rank=uniform_integer(0,length-1);
	  indexes[i]=index_temp[rank];
	  
	  if (rank<length-1)	// Compact
	    {
	      for (t=rank;t<length-1;t++)
		index_temp[t]=index_temp[t+1];
	    }					
	  length=length-1;
	}
      delete[] index_temp;
    }

  }// namespace math
} // namespace popot

#endif // POPOT_MATH_H
