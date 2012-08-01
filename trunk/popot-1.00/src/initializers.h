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
	static double init(void * params)
	{
	  // params = [min, max]
	  // returns 0
	  return 0.0;
	}
      };

      class PositionUniformRandom
      {
      public:
	static double init(void * params)
	{
	  // params = [min, max]
	  // returns U(min, max)
	  double * bounds = (double*) params;
	  return popot::math::uniform_random(bounds[0], bounds[1]);
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
