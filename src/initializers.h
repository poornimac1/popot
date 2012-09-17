#ifndef POPOT_INITIALIZER_H
#define POPOT_INITIALIZER_H

// For various initialization procedures and bounds handling, see
// Helwig, Wanka (2008) "Theoretical Analysis of Initial Particle Swarm Behavior"
// In Proceedings of the 10th International Conference on Parallel Problem Solving from Nature, G. Rudolph et al. (Eds.),LNCS 5199, pages 889–898

namespace popot
{
  namespace PSO
  {
    namespace initializer
    {

      class PositionZero
      {
      public:
	static void init(double (*lowerBound)(int), double (*upperBound)(int), int size, double * positions)
	{
	  // params = [min, max]
	  // returns 0
	  for(int i = 0 ; i < size ; ++i)
	    positions[i] = 0.0;
	}
      };

      class PositionUniformRandom
      {
      public:
	static void init(double (*lowerBound)(int), double (*upperBound)(int), int size, double * positions)
	{
	  for(int i = 0 ; i < size ; ++i)
	    positions[i] = popot::math::uniform_random(lowerBound(i),upperBound(i));
	}
	static double init(double lbound, double ubound)
	{
	  return popot::math::uniform_random(lbound,ubound);
	}
      };

      template< int D>
      class PositionHalton
      {
      public:
	static void init(double (*lowerBound)(int), double (*upperBound)(int), int size, double * positions)
	{
	  if(D != size)
	    std::cerr << "WARNING : Using Halton<" << D << "> for generating " << size << " dimension quasi random numbers " << std::endl;

	  // Get the halton sequence, in [0,1]
	  popot::rng::Halton<D>::rng_rand(positions);
	  
	  // Scale it into min, max
	  for(int i = 0 ; i < size ; ++i)
	    positions[i] = lowerBound(i)+ (upperBound(i)-lowerBound(i))* positions[i];
	}
      };

      class VelocityZero
      {
      public:
	static double init(void * params)
	{
	  // params = [min, max, xi]
	  // returns 0
	  return 0.0;
	}
      };

      class VelocityHalfDiff
      {
      public:
	static double init(void * params)
	{
	  // params = [min, max, xi]
	  // returns (U(min,max) - xi)/2.0
	  double * dparams = (double*) params;
	  return 0.5*(popot::math::uniform_random(dparams[0], dparams[1])-dparams[2]);
	}
      };

      class VelocitySPSO2011
      {
      public:
	static double init(void * params)
	{
	  // dparams = [min, max, xi]
	  // returns U(min,max) - xi
	  double * dparams = (double*) params;
	  return popot::math::uniform_random(dparams[0], dparams[1])-dparams[2];
	}
      };

    }
  }
}

#endif
