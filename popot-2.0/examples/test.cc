// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

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
  return (f <= 1e-2) || (epoch >= 10000);
}

int main(int argc, char * argv[])
{
  RNG_GENERATOR::rng_srand();

  typedef popot::PSO::particle::Particle<> PARTICLE_TYPE ;
  PARTICLE_TYPE b(10);

  std::cout << b << std::endl;
  
  auto lbound = [] (size_t index) -> double { return -10; };
  auto ubound = [] (size_t index) -> double { return  10; };
  auto init_func = [lbound, ubound] (PARTICLE_TYPE& p) -> void { popot::initializer::position::uniform_random<PARTICLE_TYPE::VECTOR_TYPE>(p.getPosition(), lbound, ubound);};

  //popot::initializer::position::zero(b.getPosition());
  popot::initializer::position::zero(b.getPosition());
  popot::initializer::position::zero(b.getBestPosition().getPosition());
  //popot::PSO::particle::init_position(b, popot::initializer::position::zero<PARTICLE_TYPE::VECTOR_TYPE>);
  init_func(b);
  b.evaluateFitness(evaluate);

  std::cout << b << std::endl;
  
  auto algo = popot::algorithm::spso2006(5,
					 [] (size_t index) -> double { return -10; },
					 [] (size_t index) -> double { return  10; },
					 stop, evaluate);
  algo.print();
  //
  auto position_update = popot::PSO::particle::updatePosition<popot::algorithm::ParticleSPSO2006>;
  position_update(algo.getParticles()[0]);
  //decltype(popot::PSO::topology::ring_fillNeighborhoods<popot::PSO::particle::Particle<> >) titi = popot::PSO::topology::ring_fillNeighborhoods<PARTICLE_TYPE>;

  algo.step();

}
