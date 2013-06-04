#include <iostream>
#include <array>

typedef std::array<double, 3>  Param;

typedef pso::by_default::Op<Param> Op;

double ubound(size_t index)
{
  
}
double lbound(size_t index)
{
  
}

double fitness(const Param& p)
{
}

int main(int argc, char * argv[])
{
  Op op;
  Param p1, p2;
  op.add(p1, p2);

  /*
  auto initp = pso::init::random(Param(), lbound, ubound);

  auto part = pso::particle(initp, fitness); // decltype(initp()) toto;
  auto spedpart = pso::speed(part);
  auto swarm = pso::grid(5,5, Particle(), initp);
  */

  pso::synchronous_step(swarm, op, fitness);
  pso::asynchronous_step();




  auto swarm = pso::make_std2011(lbound, ubound, fitness);

}

