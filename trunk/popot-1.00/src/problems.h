#ifndef POPOT_PROBLEMS_H
#define POPOT_PROBLEMS_H

namespace popot
{
  /**
   * Namespace where all the concepts are defined
   */
  namespace concept
  {
    /**
     * A problem defines the cost function to minimize
     * This class gives the required methods to implement 
     * in your problem
     * @brief
     */
    class Problem
    {
    public:     
      /**
       * Dimension of the parameter vector
       */
      static const int nb_parameters;

      /**
       * This gives the lower bound for a given parameter component
       * @param index The index of the component you want to get the lower bound.
       */
      static double get_lbound(int index);

      /**
       * This gives the upper bound for a given parameter component
       * @param index The index of the component you want to get the upper bound.
       */
      static double get_ubound(int index);

      /**
       * Determines if we should stop the optimization algorithm
       */
      static bool stop(double fitness, int epoch);
	  
      /**
       * Returns the value of the cost function for the given parameters
       * @param params A vector of parameters.
       */
      static double evaluate(void * params);

    };

  } // namespace concept

  /**
   * Definition of some benchmarking problems
   */
  namespace problems
  {
    /**
     * N-dimensional Ackley function
     * @brief \f$ 20 (1 - \exp(-0.2 * \sqrt{\frac{1}{N} \sum_{i=1}^{N} x_i^2})) + \exp(0) - \exp(\frac{1}{N} \sum_{i=1}^{N} \cos(2\pi x_i) )\f$
     *  Bounds [-30; 30]
     */
    template< int dimension>
    class Ackley
    {
    public:
      static const int nb_parameters = dimension;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	return 32.0;
      }

      static double get_ubound(int index)
      {
	return -32.0;
      }

      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 0.0;
	double cos_x = 0.0;
	double sq_x = 0.0;
	for(int i = 0 ; i < nb_parameters ; ++i)
	  {
	    sq_x += pow(params[i],2.0);
	    cos_x += cos(2.0 * M_PI * params[i]);
	  }
	fit = 20.0 * (1.0 - exp(-0.2 * sqrt(1.0/double(nb_parameters) * sq_x)))
	  + exp(1) - exp(1.0 / double(nb_parameters) * cos_x);
	return fit;
      }
    };
    template<int dimension> int Ackley<dimension>::count;

    /**
     * N-dimensional Quadric function
     * @brief \f$ \sum_{i=1}^{N} (\sum_{j=1}^{i} x_i)^2 \f$
     *  Bounds [-100; 100]
     */
    template< int dimension>
    class Quadric
    {
    public:
      static const int nb_parameters = dimension;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	return -100;
      }

      static double get_ubound(int index)
      {
	return 100;
      }

      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 0.0;
	double tmp_fit = 0.0;
	for(int i = 0 ; i < nb_parameters ; ++i)
	  {
	    tmp_fit += params[i];
	    fit += tmp_fit*tmp_fit;
	  }
	return fit;
      }
    };
    template<int dimension> int Quadric<dimension>::count;


    /**
     * N-dimensional Griewank function
     * @brief \f$ 1 + \frac{1}{4000} \sum_{i=1}^{N} x_i^2 - \prod_i cos(\frac{x_i}{\sqrt{i}})\f$
     *        Bounds [-600;600]
     */
    template< int dimension>
    class Griewank
    {
    public:
      static const int nb_parameters = dimension;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	return -600;
      }

      static double get_ubound(int index)
      {
	return 600;
      }

      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 0.0;
	double cos_x = 1.0;
	double sq_x = 0.0;
	for(int i = 0 ; i < nb_parameters ; ++i)
	  {
	    sq_x += pow(params[i],2.0);
	    cos_x *= cos(params[i]/sqrt(double(i+1)));
	  }
	fit = 1.0 + 1.0/4000.0 * sq_x - cos_x;
	return fit;
      }
    };
    template<int dimension> int     Griewank<dimension>::count;

    /**
     * N-dimensional Sphere function
     * @brief \f$ \sum_{i=1}^{N} x_i^2\f$
     * Bounds [-100;100]
     */
    template< int dimension>
    class Sphere
    {
    public:
      static const int nb_parameters = dimension;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	return -100;
      }

      static double get_ubound(int index)
      {
	return 100;
      }


      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 0.0;
	for(int i = 0 ; i < nb_parameters ; ++i)
	  fit += pow(params[i],2.0);
	return fit;
      }
    };
    template<int dimension> int     Sphere<dimension>::count;

    /**
     * N-dimensional Quartic function with uniform noise
     * @brief \f$ \sum_{i=1}^{N} (i x_i^4 + n_i)\f$
     * \f$ n_i \sim U(0,1)\f$ , Bounds [-1.28,1.28]
     */
    template< int dimension>
    class QuarticNoise
    { 
    public:
      static const int nb_parameters = dimension;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	return 1.28;
      }

      static double get_ubound(int index)
      {
	return -1.28;
      }

      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 0.0;
	for(int i = 0 ; i < nb_parameters ; ++i)
	  fit += double(i+1)*pow(params[i],4.0) + popot::math::uniform_random(0.0,1.0);
	return fit;
      }
    };
    template<int dimension> int     QuarticNoise<dimension>::count;

    /**
     * N-dimensional Rastrigin function
     * @brief \f$ \sum_{i=1}^{N} (x_i^2 + 10 (1 - cos(2\pi x_i)))\f$
     * Bounds [-5.12,5.12]
     */
    template< int dimension>
    class Rastrigin
    {
    public:
      static const int nb_parameters = dimension;
      static const int suggested_swarm_size;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	if(index < 0 || index >= nb_parameters)
	  throw popot::Exception::IndexOutOfRange(index, nb_parameters);
	      
	return -5.12;
      }

      static double get_ubound(int index)
      {
	if(index < 0 || index >= nb_parameters)
	  throw popot::Exception::IndexOutOfRange(index, nb_parameters);
	      
	return 5.12;
      }

      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 0.0;
	for(int i = 0 ; i < nb_parameters ; ++i)
	  fit += pow(params[i],2.0) + 10.0*(1.0 - cos(2.0*M_PI*params[i]));
	return fit;
      }
    };
    template<int dimension> int     Rastrigin<dimension>::count;
    template<int dimension> const int Rastrigin<dimension>::suggested_swarm_size =  int(10 + 2.0 * sqrt(dimension));

    /**
     * N-dimensional Rosenbrock banana function
     * @brief \f$ \sum_{i=1}^{N-1} (100 (y_{i+1} - y_i^2)^2 + (y_i - 1)^2)\f$
     *        with \f$ y_i = x_i + o_i - 1\f$
     * Bounds [-30,30]
     */
    template< int dimension>
    class Rosenbrock
    {
    public:
      static const int nb_parameters = dimension;
      static const int suggested_swarm_size ;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	return -30;
      }

      static double get_ubound(int index)
      {
	return 30;
      }
	
      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 0.0;
	double y_i, y_i_1;
	for(int i = 0 ; i < nb_parameters-1 ; ++i)
	  {
	    y_i = params[i];
	    y_i_1 = params[i+1];
	    fit += 100 * pow(y_i_1 - pow(y_i,2.0),2.0)+pow(y_i - 1.0,2.0);
	  }
	return fit;
      }
    };
    template<int dimension> int     Rosenbrock<dimension>::count;

    /**
     * N-dimensional Schwefel1_2 function
     * @brief \f$ \sum_{i=1}^{N} (\sum_{j=1}^{i} x_i)^2\f$
     * Bounds [-500,500]
     */
    template< int dimension>
    class Schwefel1_2
    {
    public:
      static const int nb_parameters = dimension;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	return -500;
      }

      static double get_ubound(int index)
      {
	return 500;
      }

      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 0.0;
	double sum = 0.0;
	for(int i = 0 ; i < nb_parameters ; ++i)
	  {
	    sum += params[i];
	    fit += pow(sum, 2.0);
	  }
	return fit;
      }
    };
    template<int dimension> int     Schwefel1_2<dimension>::count;

    /**
     * N-dimensional Schwefel function
     * @brief \f$ 418.9829 N - \sum_{i=1}^{N} x_i \sin(\sqrt{|x_i|}) \f$
     * Bounds [-500,500]
     */
    template< int dimension>
    class Schwefel
    {
    public:
      static const int nb_parameters = dimension;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	return -500;
      }

      static double get_ubound(int index)
      {
	return 500;
      }

      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 418.9829 * nb_parameters;
	for(int i = 0 ; i < nb_parameters ; ++i)
	  {
	    fit += -params[i]*sin(sqrt(fabs(params[i])));
	  }
	return fit;
      }
    };
    template<int dimension> int     Schwefel<dimension>::count;

    /**
     * N-dimensional Salomon function
     * @brief \f$ 1 + 0.1 \sqrt{\sum_{i=1}^{N} x_i^2} - cos(2\pi \sum_{i=1}^{N} x_i^2)\f$
     * Bounds [-600,600]
     */
    template< int dimension>
    class Salomon
    {
    public:
      static const int nb_parameters = dimension;
      static int count;

      static void init(void)
      {
	count = 0;
      }

      static void free(void)
      {
      }

      static double get_lbound(int index)
      {
	return -600;
      }

      static double get_ubound(int index)
      {
	return 600;
      }

      static bool stop(double fitness, int epoch)
      {
	return (fitness <= 1e-4) || (count >= 10000*dimension);
      }

      static double evaluate(void * x)
      {
	double * params = (double*) x;
	count++;
	double fit = 0.0;
	double sum_sq = 0.0;
	for(int i = 0 ; i < nb_parameters ; ++i)
	  {
	    sum_sq += pow(params[i],2.0);
	  }
	fit = -cos(2.0 * M_PI * sum_sq) + 0.1 * sqrt(sum_sq) + 1.0;
	return fit;
      }
    };
    template<int dimension> int     Salomon<dimension>::count;

    /**
     * 2-dimensional Dropwave function
     * @brief \f$ -\frac{1 + \cos(12 \sqrt{x_0^2 + x_1^2} + x_1^2)}{0.5 (x_0^2 + x_1^2) + 2} \f$
     * Bounds [-100,100]
     */
    class Dropwave
    {
    public:
      static const int nb_parameters = 2;
      static int count;
      static void init(void);
      static void free(void);
      static double get_lbound(int index);
      static double get_ubound(int index);
      static bool stop(double fitness, int epoch);
      static double evaluate(void * x);
    };


    // The test functions of the CEC2005 benchmark
    // "Problem Definitions and Evaluation Criteria for the CEC 2005 Special Session on Real-Parameter Optimization", 2005
    // P. N. Suganthan, N. Hansen, J. J. Liang, K. Deb, Y. -P. Chen, A. Auger, S. Tiwari
    namespace CEC2005
    {
    }

    /**
     * Functions from BBOB : http://coco.lri.fr/downloads/download11.06/bbobc.tar.gz
     */
    namespace BBOB
    {
    }

    /**
     * Functions from the C code of SPSO 2011
     */
    namespace SPSO2011
    {

      /**
       * F0 : Sphere
       */
      template <int dimension>
      class F0
      {
      public:
	static const int nb_parameters = dimension;
	static int count;

	static void init(void)
	{
	  count = 0;
	}

	static void free(void)
	{
	}

	static double get_lbound(int index)
	{
	  return -100;
	}

	static double get_ubound(int index)
	{
	  return 100;
	}
	
	static std::string name(void)
	{
	  return "SPSO2011-F0";
	}

	static bool has_failed(double fitness)
	{
	  return fitness > 1e-2;
	}

	static bool stop(double fitness, int epoch)
	{
	  return (fitness <= 1e-2) || (count >= 75000);
	}

	static double evaluate(void * params)
	{
	  double * dparams = (double*)params;
	  count ++;
	  double f = 0.0;
	  for(int i = 0 ; i < dimension ; ++i)
	    f += dparams[i]*dparams[i];
	  return f;
	}
      };
      template<int dimension> int      F0<dimension>::count;

      /**
       * F1 : N-dimensional Griewank function
       * @brief \f$ 1 + \frac{1}{4000} \sum_{i=1}^{N} x_i^2 - \prod_i cos(\frac{x_i}{\sqrt{i}})\f$
       *        Bounds [-600;600]
       */
      template< int dimension>
      class F1
      {
      public:
	static const int nb_parameters = dimension;
	static int count;

	static void init(void)
	{
	  count = 0;
	}

	static void free(void)
	{
	}

	static double get_lbound(int index)
	{
	  return -600;
	}

	static double get_ubound(int index)
	{
	  return 600;
	}
	
	static std::string name(void)
	{
	  return "SPSO2011-F1";
	}

	static bool has_failed(double fitness)
	{
	  return fitness > 0.05;
	}

	static bool stop(double fitness, int epoch)
	{
	  return (fitness <= 0.05) || (count >= 75000);
	}

	static double evaluate(void * x)
	{
	  double * params = (double*) x;
	  count++;
	  double fit = 0.0;
	  double cos_x = 1.0;
	  double sq_x = 0.0;
	  for(int i = 0 ; i < nb_parameters ; ++i)
	    {
	      sq_x += pow(params[i],2.0);
	      cos_x *= cos(params[i]/sqrt(double(i+1)));
	    }
	  fit = 1.0 + 1.0/4000.0 * sq_x - cos_x;
	  return fit;
	}
      };
      template<int dimension> int F1<dimension>::count;

      /**
       * F2 : N-dimensional Rosenbrock banana function
       * @brief \f$ \sum_{i=1}^{N-1} (100 (y_{i+1} - y_i^2)^2 + (y_i - 1)^2)\f$
       *        with \f$ y_i = x_i + o_i - 1\f$
       * Bounds [-30,30]
       */

      template< int dimension>
      class F2
      {
      public:
	static const int nb_parameters = dimension;
	static int count;

	static void init(void)
	{
	  count = 0;
	}

	static void free(void)
	{
	}

	static double get_lbound(int index)
	{
	  return -30;
	}

	static double get_ubound(int index)
	{
	  return 30;
	}
	
	static std::string name(void)
	{
	  return "SPSO2011-F2";
	}

	static bool has_failed(double fitness)
	{
	  return fitness > 100;
	}

	static bool stop(double fitness, int epoch)
	{
	  return (fitness <= 100) || (count >= 75000);
	}

	static double evaluate(void * x)
	{
	  double * params = (double*) x;
	  count++;
	  double fit = 0.0;
	  double y_i, y_i_1;
	  for(int i = 0 ; i < nb_parameters-1 ; ++i)
	    {
	      y_i = params[i];
	      y_i_1 = params[i+1];
	      fit += 100 * pow(y_i_1 - pow(y_i,2.0),2.0)+pow(y_i - 1.0,2.0);
	    }
	  return fit;
	}
      };
      template<int dimension> int F2<dimension>::count;

      /**
       * F3 : N-dimensional Rastrigin function
       * @brief \f$ \sum_{i=1}^{N} (x_i^2 + 10 (1 - cos(2\pi x_i)))\f$
       * Bounds [-5.12,5.12]
       */

      template< int dimension>
      class F3
      {
      public:
	static const int nb_parameters = dimension;
	static int count;

	static void init(void)
	{
	  count = 0;
	}

	static void free(void)
	{
	}

	static double get_lbound(int index)
	{
	  return -5.12;
	}

	static double get_ubound(int index)
	{
	  return 5.12;
	}
	
	static std::string name(void)
	{
	  return "SPSO2011-F3";
	}

	static bool has_failed(double fitness)
	{
	  return fitness > 50;
	}

	static bool stop(double fitness, int epoch)
	{
	  return (fitness <= 50) || (count >= 75000);
	}

	static double evaluate(void * x)
	{
	  double * params = (double*) x;
	  count++;
	  double fit = 0.0;
	  for(int i = 0 ; i < nb_parameters ; ++i)
	    fit += pow(params[i],2.0) + 10.0*(1.0 - cos(2.0*M_PI*params[i]));
	  return fit;
	}
      };
      template<int dimension> int F3<dimension>::count;


      /**
       * F4
       * 2-dimensional Tripod function
       * @brief \f$  \f$
       * Bounds [-100,100]
       */
      class F4
      {
      public:
	static const int nb_parameters = 2;
	static int count;
	static void init(void);
	static void free(void);
	static double get_lbound(int index);
	static double get_ubound(int index);
	
	static std::string name(void);
	static bool has_failed(double fitness);
	static bool stop(double fitness, int epoch);
	static int sign(double x);
	static double evaluate(void * x);
      };

      /**
       * F5 : N-dimensional Ackley function
       * @brief \f$ 20 (1 - \exp(-0.2 * \sqrt{\frac{1}{N} \sum_{i=1}^{N} x_i^2})) + \exp(0) - \exp(\frac{1}{N} \sum_{i=1}^{N} \cos(2\pi x_i) )\f$
       *  Bounds [-30; 30]
       */
      
      template< int dimension>
      class F5
      {
      public:
	static const int nb_parameters = dimension;
	static int count;
	  
	static void init(void)
	{
	  count = 0;
	}
	  
	static void free(void)
	{
	}
	  
	static double get_lbound(int index)
	{
	  return -30;
	}

	static double get_ubound(int index)
	{
	  return 30;
	}
	
	static std::string name(void)
	{
	  return "SPSO2011-F5";
	}

	static bool has_failed(double fitness)
	{
	  return fitness > 0;
	}

	static bool stop(double fitness, int epoch)
	{
	  return (fitness <= 0) || (count >= 80000);
	}

	static double evaluate(void * x)
	{
	  double * params = (double*) x;
	  count++;
	  double fit = 0.0;
	  double cos_x = 0.0;
	  double sq_x = 0.0;
	  for(int i = 0 ; i < nb_parameters ; ++i)
	    {
	      sq_x += pow(params[i],2.0);
	      cos_x += cos(2.0 * M_PI * params[i]);
	    }
	  fit = 20.0 * (1.0 - exp(-0.2 * sqrt(1.0/double(nb_parameters) * sq_x)))
	    + exp(1) - exp(1.0 / double(nb_parameters) * cos_x);
	  return fit;
	}
      };
      template<int dimension> int F5<dimension>::count;

      /**
       * F6 : N-dimensional Schwefel function
       *  Bounds [-500; 500]
       */
      template< int dimension>
      class F6
      {
      public:
	static const int nb_parameters = dimension;
	static const double fobj = -12569.5;
	static const double epsilon = 2569.5;
	static int count;
	  
	static void init(void)
	{
	  count = 0;
	}
	  
	static void free(void)
	{
	}
	  
	static double get_lbound(int index)
	{
	  return -500;
	}

	static double get_ubound(int index)
	{
	  return 500;
	}
	
	static std::string name(void)
	{
	  return "SPSO2011-F6";
	}

	static bool has_failed(double fitness)
	{
	  return fitness > epsilon;
	}

	static bool stop(double fitness, int epoch)
	{
	  return (fitness <= epsilon) || (count >= 300000);
	}

	static double evaluate(void * x)
	{
	  double * params = (double*) x;
	  count++;
	  double fit = 0.0;
	  for(int i = 0 ; i < dimension ; ++i)
	    fit -= params[i] * sin(sqrt(fabs(params[i])));
	  return fabs(fit-fobj);
	}
      };
      template<int dimension> int F6<dimension>::count;

      /**
       * F7 : N-dimensional Schwefel 1.2 function
       *  Bounds [-100; 100]
       */
      template< int dimension>
      class F7
      {
      public:
	static const int nb_parameters = dimension;
	static int count;
	  
	static void init(void)
	{
	  count = 0;
	}
	  
	static void free(void)
	{
	}
	  
	static double get_lbound(int index)
	{
	  return -100;
	}

	static double get_ubound(int index)
	{
	  return 100;
	}
	
	static std::string name(void)
	{
	  return "SPSO2011-F7";
	}

	static bool has_failed(double fitness)
	{
	  return fitness > 0;
	}

	static bool stop(double fitness, int epoch)
	{
	  return (fitness <= 0) || (count >= 40000);
	}

	static double evaluate(void * x)
	{
	  double * params = (double*) x;
	  count++;
	  double fit = 0.0;
	  double sum = 0.0;
	  for(int i = 0 ; i < dimension ; ++i)
	    {
	      sum = 0.0;
	      for(int j = 0 ; j <= i ; ++j)
		sum += params[j];
	      fit += sum * sum;
	    }
	  return fit;
	}
      };
      template<int dimension> int F7<dimension>::count;


      /**
       * F8 : N-dimensional Schwefel 2.22 function
       *  Bounds [-10; 10]
       */
      template< int dimension>
      class F8
      {
      public:
	static const int nb_parameters = dimension;
	static int count;
	  
	static void init(void)
	{
	  count = 0;
	}
	  
	static void free(void)
	{
	}
	  
	static double get_lbound(int index)
	{
	  return -10;
	}

	static double get_ubound(int index)
	{
	  return 10;
	}
	
	static std::string name(void)
	{
	  return "SPSO2011-F8";
	}

	static bool has_failed(double fitness)
	{
	  return fitness > 1e-4;
	}

	static bool stop(double fitness, int epoch)
	{
	  return (fitness <= 1e-4) || (count >= 100000);
	}

	static double evaluate(void * x)
	{
	  double * params = (double*) x;
	  count++;
	  double sum = 0.0;
	  double prod = 1.0;
	  for(int i = 0 ; i < dimension ; ++i)
	    {
	      sum += fabs(params[i]);
	      prod *= fabs(params[i]);
	    }
	  return sum+prod;
	}
      };
      template<int dimension> int F8<dimension>::count;

      /**
       * F9 : N-dimensional Neumaier 3 function
       *  Bounds [-d^2; d^2]
       */
      template< int dimension>
      class F9
      {
      public:
	static const int nb_parameters = dimension;
	static int count;
	  
	static void init(void)
	{
	  count = 0;
	}
	  
	static void free(void)
	{
	}
	  
	static double get_lbound(int index)
	{
	  return -dimension*dimension;
	}

	static double get_ubound(int index)
	{
	  return dimension*dimension;
	}
	
	static std::string name(void)
	{
	  return "SPSO2011-F9";
	}

	static bool has_failed(double fitness)
	{
	  return fitness > 0;
	}

	static bool stop(double fitness, int epoch)
	{
	  return (fitness <= 0) || (count >= 40000);
	}

	static double evaluate(void * x)
	{
	  double * params = (double*) x;
	  count++;
	  double fit = 0.0;
	  for(int i = 0 ; i < dimension ; ++i)
	    fit += (params[i]-1)*(params[i]-1);
	  for(int i = 1 ; i < dimension ; ++i)
	    fit += (params[i])*(params[i-1]);  
	  return fit;
	}
      };
      template<int dimension> int F9<dimension>::count;

    }

  } // namespace problems
} // namespace popot

#endif // POPOT_PROBLEMS_H
