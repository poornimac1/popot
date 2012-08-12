#include <iostream>
#include <cmath>

// Definition of the problems for the benchmark
namespace BenchmarkProblems
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
      fit = sq_x/4000.0 - cos_x + 1.0;
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
	fit += params[i]*params[i] - 10.0*cos(2.0*M_PI*params[i]);
      fit += 10.0 * nb_parameters;
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
    static void init(void)
    {
      count = 0;
    }
    static void free(void){}
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
      return "Tripod";
    }
    static bool has_failed(double fitness)
    {
      return fitness > 1e-4;
    }
    static bool stop(double fitness, int epoch)
    {
      return (fitness <= 1e-4) || (count >= 1e4);
    }

    static int sign(double x)
    {
      if(x > 0)
	return 1;
      else
	return -1;
    }

    static double evaluate(void * x)
    {
      double * params = (double*) x;
      count++;
      double s11 = (1.0 - sign(params[0]))/2.0;
      double s12 = (1.0 + sign(params[0]))/2.0;
      double s21 = (1.0 - sign(params[1]))/2.0;
      double s22 = (1.0 + sign(params[1]))/2.0;

      double f;
      //f = s21 * (fabs (params[0]) - params[1]); // Solution on (0,0)
      f = s21 * (fabs (params[0]) +fabs(params[1]+50)); // Solution on (0,-50)  
      f = f + s22 * (s11 * (1 + fabs (params[0] + 50) +
			    fabs (params[1] - 50)) + s12 * (2 +
							    fabs (params[0] - 50) +
							    fabs (params[1] - 50)));
      //std::cout << "Evaluatation at " << params[0] << ";" << params[1] << " -> " << f << std::endl;

      return fabs(f);
    }
  };
  int F4::count;

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
