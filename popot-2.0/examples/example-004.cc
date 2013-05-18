

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"
// Define our problem
// Here we use a predefined parametrized problem
typedef popot::problems::Dropwave Problem;

class ABC_Params
{
public:
  static int ColonySize() { return 50;};
};

class Stop_Criteria
{
public:
  static bool stop(double f, int epoch)
  {
    return (f <= -1 + 1e-4) || (epoch >= 100000);
  }
};

typedef popot::ABC::algorithm::Base<ABC_Params, Problem, Stop_Criteria> ABC;

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();

  /*
  Problem p(50);
  ABC_Params params;
  Stop_Criteria stop;
  auto algo = popot::ABC::algorithm::base(params, p, stop);
  */

  Problem p;
  ABC algo(p);
  
  algo.init();
  algo.run();

  std::cout << "Best minimum found :" << algo.getBest().getFValue() << " in " << algo.getEpoch() << " steps " << std::endl;
  std::cout << "Position of the optimum : " << std::endl;
  for(size_t i = 0 ; i < p.dimension ; ++i)
    std::cout << algo.getBest().getValueAt(i) << " ";
  std::cout << std::endl;
  std::cout << (int)p.count << " Function Evaluations" << std::endl;
  
}
