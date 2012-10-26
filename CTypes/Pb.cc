// gcc -shared -O3 -o Pb.so Pb.cc
#include <math.h>
#include <cstdlib>

class MyPb
{
public:
  static double testfunc(double x)
  {
    return (x-1)*(x-1);
  }
};


extern "C"
{
  int nb_params = 10;

  int get_size() { return nb_params;}

  double evaluate(double * params)
  {
    double val = 0.0;
    for(int i = 0 ; i < nb_params ; ++i)
      {
	val += MyPb::testfunc(params[i]);
      }
    return val+(2.0*rand()/double(RAND_MAX)-1.0);
  }
}
