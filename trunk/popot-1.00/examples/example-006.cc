#include "popot.h"


// Define our problem
// Here we use a predefined parametrized problem
typedef popot::problems::Rosenbrock<50> Problem;

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
    return (epoch >= 1e4) || (f <= 1e-15);
  }
};

typedef popot::ABC::algorithm::Base<ABC_Params, Problem, Stop_Criteria> ABC;

typedef popot::Vector<double, Problem::nb_parameters> TVector;

int main(int argc, char * argv[])
{
  srand(time(NULL));
  Problem::init();

  ABC algo;
  algo.init();
  algo.run();
  std::cout << "Best minimum found :" << algo.getBest().getFValue() << " in " << algo.getEpoch() << " steps " << std::endl;


  Problem::free();
  
}
