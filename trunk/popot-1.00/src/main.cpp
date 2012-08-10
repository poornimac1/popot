#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"

// Dropwave 
void popot::problems::Dropwave::init(void)
{
  count = 0;
}

void popot::problems::Dropwave::free(void)
{
}

double popot::problems::Dropwave::get_lbound(int index)
{
  return -5.12;
}

double popot::problems::Dropwave::get_ubound(int index)
{
  return 5.12;
}

bool popot::problems::Dropwave::stop(double fitness, int epoch)
{
  return (fitness <= -1.0+1e-4) || (count >= 20000);
}

double popot::problems::Dropwave::evaluate(void * x)
{
  double * params = (double*) x;
  count++;
  return -(1.0 + cos(12.0*sqrt(pow(params[0],2.0) + pow(params[1],2.0))))
    / (0.5 * (pow(params[0],2.0) + pow(params[1],2.0))+2.0);
}
int     popot::problems::Dropwave::count;

// F4 Tripod
void popot::problems::SPSO2011::F4::init(void)
{
  count = 0;
}

void popot::problems::SPSO2011::F4::free(void)
{
}

double popot::problems::SPSO2011::F4::get_lbound(int index)
{
  return -100;
}

double popot::problems::SPSO2011::F4::get_ubound(int index)
{
  return 100;
}
std::string popot::problems::SPSO2011::F4::name(void)
{
  return "SPSO2011-F4";
}

bool popot::problems::SPSO2011::F4::has_failed(double fitness)
{
  return fitness > 1e-4;
}

bool popot::problems::SPSO2011::F4::stop(double fitness, int epoch)
{
  return (fitness <= 1e-4) || (count >= 1e4);
}

int popot::problems::SPSO2011::F4::sign(double x)
{
  if(x > 0)
    return 1;
  else
    return -1;
}

double popot::problems::SPSO2011::F4::evaluate(void * x)
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

int popot::problems::SPSO2011::F4::count;


int main(int argc, char * argv[])
{

}
