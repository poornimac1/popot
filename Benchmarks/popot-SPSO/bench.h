#include <iostream>
#include <cmath>

typedef popot::algorithm::ParticleSPSO::VECTOR_TYPE TVector;

// Definition of the problems for the benchmark
namespace BenchmarkProblems
{

  /**
   * F0 : N-dimensional Sphere
   */
  class F0
  {
  public:
    int _count;
    size_t _dimension;

    F0(size_t dimension) : _dimension(dimension), _count(0)  {}
    F0(const F0& p) : _dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -100;            }
    double get_ubound(int index)          { return 100;             }	
    std::string name(void)                { return "SPSO2011-F0";   }
    bool has_failed(double fitness)       { return fitness > 1e-2;  }
    bool stop(double fitness, int epoch)  { return (fitness <= 1e-2) || (_count >= 75000);   }

    double evaluate(TVector &pos)
    {
      _count ++;
      double f = 0.0;
      for(int i = 0 ; i < _dimension ; ++i)
	f += pos[i] * pos[i];
      return f;
    }
  };

  /**
   * F1 : N-dimensional Griewank function
   * @brief \f$ 1 + \frac{1}{4000} \sum_{i=1}^{N} x_i^2 - \prod_i cos(\frac{x_i}{\sqrt{i}})\f$
   *        Bounds [-600;600]
   */

  class F1
  {
  public:
    int _count;
    size_t _dimension;

    F1(size_t dimension) : _dimension(dimension), _count(0)  {}
    F1(const F1& p) : _dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -600;            }
    double get_ubound(int index)          { return 600;             }	
    std::string name(void)                { return "SPSO2011-F1";   }
    bool has_failed(double fitness)       { return fitness > 0.05;  }
    bool stop(double fitness, int epoch)  { return (fitness <= 0.05) || (_count >= 75000);   }

    double evaluate(TVector &pos)
    {
      _count++;
      double fit = 0.0;
      double cos_x = 1.0;
      double sq_x = 0.0;
      for(int i = 0 ; i < _dimension ; ++i)
	{
	  sq_x += pow(pos[i],2.0);
	  cos_x *= cos(pos[i]/sqrt(double(i+1)));
	}
      fit = sq_x/4000.0 - cos_x + 1.0;
      return fit;
    }
  };

  /**
   * F2 : N-dimensional Rosenbrock banana function
   * @brief \f$ \sum_{i=1}^{N-1} (100 (y_{i+1} - y_i^2)^2 + (y_i - 1)^2)\f$
   *        with \f$ y_i = x_i + o_i - 1\f$
   * Bounds [-30,30]
   */

  class F2
  {
  public:
    int _count;
    size_t _dimension;

    F2(size_t dimension) : _dimension(dimension), _count(0)  {}
    F2(const F2& p) : _dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -30;            }
    double get_ubound(int index)          { return 30;             }	
    std::string name(void)                { return "SPSO2011-F2";  }
    bool has_failed(double fitness)       { return fitness > 100;  }
    bool stop(double fitness, int epoch)  { return (fitness <= 100) || (_count >= 75000);  }

    double evaluate(TVector &pos)
    {
      _count++;
      double fit = 0.0;
      double y_i, y_i_1;
      for(int i = 0 ; i < _dimension-1 ; ++i)
	{
	  y_i = pos[i];
	  y_i_1 = pos[i+1];
	  fit += 100 * pow(y_i_1 - pow(y_i,2.0),2.0)+pow(y_i - 1.0,2.0);
	}
      return fit;
    }
  };

  /**
   * F3 : N-dimensional Rastrigin function
   * @brief \f$ \sum_{i=1}^{N} (x_i^2 + 10 (1 - cos(2\pi x_i)))\f$
   * Bounds [-5.12,5.12]
   */
  class F3
  {
  public:
    int _count;
    size_t _dimension;

    F3(size_t dimension) : _dimension(dimension), _count(0)  {}
    F3(const F3& p) : _dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -5.12;         }
    double get_ubound(int index)          { return 5.12;          }	
    std::string name(void)                { return "SPSO2011-F3"; }
    bool has_failed(double fitness)       { return fitness > 50;  }
    bool stop(double fitness, int epoch)  { return (fitness <= 50) || (_count >= 75000);  }

    double evaluate(TVector &pos)
    {
      _count++;
      double fit = 0.0;
      for(int i = 0 ; i < _dimension ; ++i)
	fit += pos[i]*pos[i] - 10.0*cos(2.0*M_PI*pos[i]);
      fit += 10.0 * _dimension;
      return fit;
    }
  };

  /**
   * F4
   * 2-dimensional Tripod function
   * @brief \f$  \f$
   * Bounds [-100,100]
   */
  class F4
  {
  public:
    int _count;
    size_t _dimension;

    F4(size_t dimension) : _dimension(2), _count(0)  {
      if(dimension != 2)
	std::cerr << "Warning, Tripod is working in dimension 2 not " << dimension << " !! " << std::endl;
    }
    F4(const F4& p) : _dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -100;            }
    double get_ubound(int index)          { return 100;             }	
    std::string name(void)                { return "SPSO2011-F4";   }
    bool has_failed(double fitness)       { return fitness > 1e-4;  }
    bool stop(double fitness, int epoch)  { return  (fitness <= 1e-4) || (_count >= 1e4);  }

    int sign(double x)
    {
      if(x > 0)
	return 1;
      else
	return -1;
    }

    double evaluate(TVector &pos)
    {      
      _count++;
      double s11 = (1.0 - sign(pos[0]))/2.0;
      double s12 = (1.0 + sign(pos[0]))/2.0;
      double s21 = (1.0 - sign(pos[1]))/2.0;
      double s22 = (1.0 + sign(pos[1]))/2.0;

      double f;
      //f = s21 * (fabs (pos[0]) - pos[1]); // Solution on (0,0)
      f = s21 * (fabs (pos[0]) +fabs(pos[1]+50)); // Solution on (0,-50)  
      f = f + s22 * (s11 * (1 + fabs (pos[0] + 50) +
			    fabs (pos[1] - 50)) + s12 * (2 +
							 fabs (pos[0] - 50) +
							 fabs (pos[1] - 50)));
      //std::cout << "Evaluatation at " << pos[0] << ";" << pos[1] << " -> " << f << std::endl;

      return fabs(f);

    }
  };


  /**
   * F5 : N-dimensional Ackley function
   * @brief \f$ 20 (1 - \exp(-0.2 * \sqrt{\frac{1}{N} \sum_{i=1}^{N} x_i^2})) + \exp(0) - \exp(\frac{1}{N} \sum_{i=1}^{N} \cos(2\pi x_i) )\f$
   *  Bounds [-30; 30]
   */
  class F5
  {
  public:
    int _count;
    size_t _dimension;

    F5(size_t dimension) : _dimension(dimension), _count(0)  {}
    F5(const F5& p) : _dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -30;            }
    double get_ubound(int index)          { return 30;             }	
    std::string name(void)                { return "SPSO2011-F5";  }
    bool has_failed(double fitness)       { return fitness > 0;    }
    bool stop(double fitness, int epoch)  { return (fitness <= 0) || (_count >= 80000);  }

    double evaluate(TVector &pos)
    {
      _count++;
      double fit = 0.0;
      double cos_x = 0.0;
      double sq_x = 0.0;
      for(int i = 0 ; i < _dimension ; ++i)
	{
	  sq_x += pow(pos[i],2.0);
	  cos_x += cos(2.0 * M_PI * pos[i]);
	}
      fit = 20.0 * (1.0 - exp(-0.2 * sqrt(1.0/double(_dimension) * sq_x)))
	+ exp(1) - exp(1.0 / double(_dimension) * cos_x);
      return fit;
    }
  };   

  /**
   * F6 : N-dimensional Schwefel function
   *  Bounds [-500; 500]
   */

  class F6
  {
  public:
    const double fobj;
    const double epsilon;
    int _count;
    size_t _dimension;

    F6(size_t dimension) : _dimension(dimension), _count(0), fobj(-12569.5), epsilon(2569.5)  {}
    F6(const F6& p) : fobj(p.fobj), epsilon(p.epsilon),_dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -500;            }
    double get_ubound(int index)          { return 500;             }	
    std::string name(void)                { return "SPSO2011-F6";   }
    bool has_failed(double fitness)       { return fitness > epsilon;  }
    bool stop(double fitness, int epoch)  { return (fitness <= epsilon) || (_count >= 300000);   }

    double evaluate(TVector &pos)
    {
      _count++;
      double fit = 0.0;
      for(int i = 0 ; i < _dimension ; ++i)
	fit -= pos[i] * sin(sqrt(fabs(pos[i])));
      return fabs(fit-fobj);
    }
  };

  /**
   * F7 : N-dimensional Schwefel 1.2 function
   *  Bounds [-100; 100]
   */
  class F7
  {
  public:
    int _count;
    size_t _dimension;

    F7(size_t dimension) : _dimension(dimension), _count(0)  {}
    F7(const F7& p) : _dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -100;            }
    double get_ubound(int index)          { return 100;             }	
    std::string name(void)                { return "SPSO2011-F7";   }
    bool has_failed(double fitness)       { return fitness > 0;     }
    bool stop(double fitness, int epoch)  { return (fitness <= 0) || (_count >= 40000);   }

    double evaluate(TVector &pos)
    {
      _count++;
      double fit = 0.0;
      double sum = 0.0;
      for(int i = 0 ; i < _dimension ; ++i)
	{
	  sum = 0.0;
	  for(int j = 0 ; j <= i ; ++j)
	    sum += pos[j];
	  fit += sum * sum;
	}
      return fit;
    }
  };


  /**
   * F8 : N-dimensional Schwefel 2.22 function
   *  Bounds [-10; 10]
   */
  class F8
  {
  public:
    int _count;
    size_t _dimension;

    F8(size_t dimension) : _dimension(dimension), _count(0)  {}
    F8(const F8& p) : _dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -10;            }
    double get_ubound(int index)          { return 10;             }	
    std::string name(void)                { return "SPSO2011-F8";  }
    bool has_failed(double fitness)       { return fitness > 1e-4; }
    bool stop(double fitness, int epoch)  { return (fitness <= 1e-4) || (_count >= 100000);   }

    double evaluate(TVector &pos)
    {
      _count++;
      double sum = 0.0;
      double prod = 1.0;
      for(int i = 0 ; i < _dimension ; ++i)
	{
	  sum += fabs(pos[i]);
	  prod *= fabs(pos[i]);
	}
      return sum+prod;
    }
  };

  /**
   * F9 : N-dimensional Neumaier 3 function
   *  Bounds [-d^2; d^2]
   */
  class F9
  {
  public:
    int _count;
    size_t _dimension;

    F9(size_t dimension) : _dimension(dimension), _count(0)  {}
    F9(const F9& p) : _dimension(p._dimension), _count(p._count) {}
    void init(void) { _count = 0;}
    double get_lbound(int index)          { return -_dimension*_dimension;            }
    double get_ubound(int index)          { return _dimension*_dimension;             }	
    std::string name(void)                { return "SPSO2011-F9";   }
    bool has_failed(double fitness)       { return fitness > 0;     }
    bool stop(double fitness, int epoch)  { return (fitness <= 0) || (_count >= 40000);   }

    double evaluate(TVector &pos)
    {
      _count++;
      double fit = 0.0;
      for(size_t i = 0 ; i < _dimension ; ++i)
	fit += (pos[i]-1)*(pos[i]-1);
      for(size_t i = 1 ; i < _dimension ; ++i)
	fit -= (pos[i])*(pos[i-1]);  
      return fit;
    }
  };
}
