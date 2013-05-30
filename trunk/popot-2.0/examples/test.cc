// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::JKissRNG RNG_GENERATOR;

#include "popot.h"


typedef popot::algorithm::ParticleSPSO::VECTOR_TYPE TVector;
typedef popot::problems::Ackley Problem;

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();

  size_t dimension = 10,
  Problem prob(dimension);
  
  auto algo = popot::algorithm::spso2011(dimension,
					 [&prob] (size_t index) -> double { return prob.get_lbound(index); },
					 [&prob] (size_t index) -> double { return  prob.get_ubound(index); },
					 [&prob] (double fitness, int epoch) -> bool { return prob.stop(fitness, epoch);},
					 [&prob] (TVector &pos) -> double { return prob.evaluate(pos.getValuesPtr());});

  algo.run(1);
  std::cout << "best particle " << algo.getBest() << std::endl;
}
