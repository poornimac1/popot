// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::JKissRNG RNG_GENERATOR;

#include "popot.h"

typedef popot::Vector<double> TVector;

double evaluate(TVector& x)
{
  double sum = 0.0;
  for(size_t i = 0 ; i < x.size() ; ++i)
    sum += x[i] * x[i];
  return sum;
}

bool stop(double f, int epoch)
{
  return (f <= 1e-10) || (epoch >= 2000);
}

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();

  auto algo = popot::algorithm::spso2006(50,
					 [] (size_t index) -> double { return -10; },
					 [] (size_t index) -> double { return  10; },
					 stop, evaluate);

  algo.generateGraph("graph.dot");
  algo.run(1);
}
