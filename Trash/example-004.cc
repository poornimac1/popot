// In this example, we use a standard PSO on a predefined problem

#include <stdio.h>
#include <cmath>
#include "popot.h"

// Define our problem
// Here we use a predefined parametrized problem
typedef popot::problems::Rastrigin<2> Problem;

// Define the parameters of our swarm
// A rule of thumb is to use
// w = 0.729844
// c1 = 1.496180
// c2 = 1.496180
// a result coming from the work of M. Clerc.
class Particle_PSO_Params
{
public:
  static double w()  { return 0.729844;}   // Inertia parameter
  static double c1() { return 1.496180;}   // Best particle position weight
  static double c2() { return 1.496180;}   // Best swarm position weight
};

// From the previous definition of the problem
// and parameters, we can now define our particles
// we here use non stochastic standard particles
typedef popot::PSO::particle::TraditionalParticle< Problem, Particle_PSO_Params > Particle;

// The parameters of the swarm algorithm
// for the algorithm we use below, we just need to provide
// the evaluation mode (synchronous or asynchronous)
class Swarm_PSO_params
{
public:
  static int evaluation_mode() { return popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION;}
};

// The topology of the particles' network
typedef popot::PSO::topology::Full<16, Particle> Topology;
//typedef popot::PSO::topology::Ring<6, Particle> Topology;
//typedef popot::PSO::topology::VonNeuman<4 ,4, Particle> Topology;

// Define the stopping criteria
class Swarm_Stop_Criteria
{
public:
    static bool stop(double fitness, int epoch)
    {
        return (fitness <= 1e-5) || (epoch >= 10000);
    }
};

// With the above definitions of the parameters, particles, topology and stopping criteria
// we can now set up our PSO algorithm
typedef popot::PSO::algorithm::Base<Swarm_PSO_params, Particle, Topology, Swarm_Stop_Criteria> PSO;


void presave_function_values()
{
  std::ofstream outfile("PlotExample004/function_values.data");
  
  int N = 200;
  double params[2];

  double x_min, x_max, y_min, y_max;
  x_min = Problem::get_lbound(0);
  x_max = Problem::get_ubound(0);
  y_min = Problem::get_lbound(1);
  y_max = Problem::get_ubound(1);

  /*x_min = -2.0;
  x_max = 2.0;
  y_min = -0.5;
  y_max = 3.0;*/

  double z_min, z_max, z_tmp;

  params[0] = x_min;
  params[1] = y_min;
  z_min = Problem::evaluate(params);
  z_max = Problem::evaluate(params);
  for(int i = 0 ; i < N ; ++i)
    {
    for(int j = 0 ; j < N ; ++j)
      {
	params[0] = x_min + (x_max - x_min) * i / double(N-1);
	params[1] = y_min + (y_max - y_min) * j / double(N-1);
	z_tmp = Problem::evaluate(params);
	z_min = z_min < z_tmp ? z_min : z_tmp;
	z_max = z_max > z_tmp ? z_max : z_tmp;	
      }
    }
  std::cout << "Function min = " << z_min << " ;  max = " << z_max << std::endl;
  
  for(int i = 0 ; i < N ; ++i)
    {
    for(int j = 0 ; j < N ; ++j)
      {
	params[0] = x_min + (x_max - x_min) * i / double(N-1);
	params[1] = y_min + (y_max - y_min) * j / double(N-1);
	outfile << params[0] << " " 
		<< params[1] << " "
		<< (Problem::evaluate(params)-z_min)/(z_max - z_min) << std::endl;
      }
    outfile << std::endl;
    }
  outfile.close();

}

void save_particle_positions(int epoch, PSO & p)
{
  std::ostringstream filename_plot;
  std::ofstream outfile_plot;
  
  std::ostringstream root_filename;
  root_filename << std::setw(5) << std::setfill('0') 
		<< epoch 
		<< std::setfill(' ');
  
  filename_plot << "PlotExample004/" << root_filename.str()
		<< ".plot";

  outfile_plot.open(filename_plot.str().c_str());

  outfile_plot << "set term gif" << std::endl
	       << "set output \"" << root_filename.str() << ".gif\"" << std::endl
	       << "set size ratio 1" << std::endl
	       << "set xrange [" << Problem::get_lbound(0) << ":" 
	       << Problem::get_ubound(0) << "];"<< std::endl
	       << "set yrange [" << Problem::get_lbound(1) << ":" 
	       << Problem::get_ubound(1) << "];"<< std::endl
	       << "set zrange [0:1];"<< std::endl
	       << "set cbrange [0:1];"<< std::endl
	       << "set title \"PSO on Rastrigin " << epoch << "\";" << std::endl
	       << "set palette defined ( 0 \"white\", 1 \"black\");"<< std::endl
	       << "set xlabel \"x\";"<< std::endl
	       << "set ylabel \"y\";" << std::endl
	       << "set cblabel \"z\";"<< std::endl
	       << "set view map;"<< std::endl
	       << "set pm3d at s;"<< std::endl
	       << "splot 'function_values.data' with pm3d notitle, \
  '-' with points notitle pt 7 ps 1.5 lc rgb \"red\""<< std::endl;


    // Put the position of the particles with z-value = ...
  std::vector<Particle*>::iterator it, it_end;

  for(it = p.getParticles()->begin() , it_end = p.getParticles()->end(); it != it_end ; ++it)
    {
      outfile_plot << (*it)->getPosition(0) << " " 
		   << (*it)->getPosition(1) << " " 
		   << "1" << std::endl;
    }
    
    outfile_plot.close();


}

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {

  // Initialize the random number generator 
  // to a specific seed so that we keep the same initial conditions
  //srand(0);

    // Initialize our problem
    // this actually allocates memory and initializes the boundaries
    Problem::init();

    // Let's create our swarm
    PSO pso;

    // Print the fitness of our particles
    printf("Before learning : \n");
    pso.print(1);
    std::cout << std::endl;

    presave_function_values();

    // We now iterate the algorithm
    // We can iterate step by step
    for(int i = 0 ; i < 100 ; ++i)
    {
        pso.step();
	save_particle_positions(i, pso);

        std::cout << '\r' << std::setw(6) << std::setfill('0') << i << " " << pso.getBest()->getFitness() << std::setw(5) << std::setfill(' ') << ' ' << std::flush;
    }
    // Or run the algorithm until the stopping criteria is met
    //pso.run();
    std::cout << "epoch : " << pso.epoch << std::endl;
    std::cout << "\n" << std::endl;

    // And display the best fitness we got
    printf("After learning : \n");
    pso.print(1);

    std::cout << " To get an animated gif of the behavior of the PSO, and if you used the source package, go into PlotExample004 and type make all . It requires convert and gnuplot" << std::endl;

    // Free the memory used by the problem (e.g. the bounds)
    Problem::free();
}
