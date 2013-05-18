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
      template<typename PROBLEM>
      class PositionZero
      {
      public:
	static void init(PROBLEM& p, double (PROBLEM::*lowerBound)(size_t), double (PROBLEM::*upperBound)(size_t), double * positions)
	{
	  // params = [min, max]
	  // returns 0
	  for(size_t i = 0 ; i < p.dimension ; ++i)
	    positions[i] = 0.0;
	}
	/*
	static double init(double lbound, double ubound)
	{
	  return 0.0;
	}
	*/
      };

      template<typename PROBLEM>
      class PositionUniformRandom
      {
      public:
	static void init(PROBLEM& p, double (PROBLEM::*lowerBound)(size_t), double (PROBLEM::*upperBound)(size_t), double * positions)
	{
	  for(size_t i = 0 ; i < p.dimension ; ++i)
	    positions[i] = popot::math::uniform_random((p.*(lowerBound))(i),(p.*(upperBound))(i));
	}
	/*
	static double init(double lbound, double ubound)
	{
	  return popot::math::uniform_random(lbound,ubound);
	}
	*/
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
