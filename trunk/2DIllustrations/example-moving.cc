// In this example, we use a standard PSO on a predefined problem

// Compile with :
// g++ -o example-moving example-moving.cc `pkg-config --libs --cflags popot` -O3

#include "rng_generators.h"
typedef popot::rng::JKissRNG RNG_GENERATOR;

#include <stdio.h>
#include <cmath>
#include "popot.h"

// Define our problem
// Here we use a new problem
class Problem
{
public:
  static const int nb_parameters = 2;
  static int epoch;
public:
  
  static double get_lbound(int index) { return -10;}
  static double get_ubound(int index) { return 10;}

  static void step()
  {
    epoch ++;
  }

  static double evaluate(void * params)
  {
    double * dparams = (double*) params;
    // A single moving bump making a turn in 100 steps
    double x0 = 5.0 * cos(double(epoch) * 2.0 * M_PI / 100.);
    double x1 = 5.0 * sin(double(epoch) * 2.0 * M_PI / 100.);
    double d2 = (dparams[0] - x0)*(dparams[0] -x0)
      + (dparams[1] - x1) * (dparams[1] - x1);
    return 1.0 - exp(-d2/(2.0 * 5.0 * 5.0)); 
  }
};
int Problem::epoch = 0;

// We use SPSO 2011
typedef popot::PSO::SPSO2011::PSO<Problem>::Type PSO;


void save_function_values(int epoch)
{
  std::ostringstream str;
  str.str("");
  str << "PlotExampleMoving/function_values" << std::setw(5) << std::setfill('0') << epoch << ".data";
  std::ofstream outfile(str.str().c_str());
  
  int N = 200;
  double params[2];

  double x_min, x_max, y_min, y_max;
  x_min = Problem::get_lbound(0);
  x_max = Problem::get_ubound(0);
  y_min = Problem::get_lbound(1);
  y_max = Problem::get_ubound(1);

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
  //std::cout << "Function min = " << z_min << " ;  max = " << z_max << std::endl;
  
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
  
  filename_plot << "PlotExampleMoving/" << root_filename.str()
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
	       << "splot 'function_values" << std::setw(5) << std::setfill('0') << epoch << ".data' with pm3d notitle, \
  '-' with points notitle pt 7 ps 1.5 lc rgb \"red\""<< std::endl;


  // Put the position of the particles with z-value = ...
  popot::PSO::SPSO2011::Particle<Problem>::Type* it, *it_end;

  for(it = &((p.getParticles())[0]) , it_end = it + p.getSize(); it != it_end ; ++it)
    {
      outfile_plot << it->getPosition(0) << " " 
		   << it->getPosition(1) << " " 
		   << "1" << std::endl;
    }
    
  outfile_plot.close();

}

void save_connectivity(int epoch, PSO&p)
{
  popot::PSO::SPSO2011::Particle<Problem>::Type* it, *it_end;
  
  // Connectivity matrix : m[i*width+j] : particule i informs particule j
  int* connectivity_matrix = new int[p.getSize()*p.getSize()];
  for(int i = 0 ; i < p.getSize()*p.getSize() ; ++i)
    connectivity_matrix[i] = 0;

  std::vector<int>* neigh;
  for(int i = 0 ; i < p.getSize() ; ++i)
    {
      neigh = &(p.getNeighborhoodMembership()->at(i));
      for(int j = 0 ; j < neigh->size() ; ++ j)
	{
	  connectivity_matrix[neigh->at(j) * p.getSize() + i] = 1;
	}
    }



  // Generate the Graphviz file
  std::ostringstream filename;
  std::ostringstream root_filename;
  root_filename << std::setw(5) << std::setfill('0') 
		<< epoch 
		<< std::setfill(' ');
  
  filename << "PlotExampleMoving/" << root_filename.str()
		<< ".dot";
  std::ofstream outfile(filename.str().c_str());

  // Define the nodes :
  outfile << "digraph \"connectivity\" { " << std::endl;
  double x, y;
  for(int i = 0 ; i < p.getSize() ; ++i)
    {
      x = 5.0 * cos(2.0 * M_PI * i / p.getSize());
      y = 5.0 * sin(2.0 * M_PI * i / p.getSize());
      outfile << "PART_" << i << "[pos=\""<< x << "," << y << "!\", label=\""<<i << "\"];" << std::endl;
    }
  // Define the connections
  for(int i = 0 ; i < p.getSize() ; ++i)
    {
      if(connectivity_matrix[i*p.getSize() + i])
	{
	  outfile << "PART_" << i << " -> PART_" << i << ";" << std::endl;
	}

      for(int j = 0 ; j < i ; ++j)
	{
	  // If i informs j
	  if(connectivity_matrix[i*p.getSize() + j] && connectivity_matrix[j*p.getSize() + i])
	    {
	      outfile << "PART_" << i << " -> PART_" << j << " [dir=both];" << std::endl;
	    }
	  else if(connectivity_matrix[i*p.getSize() + j])
	    {
	      outfile << "PART_" << i << " -> PART_" << j <<";" << std::endl;
	    }
	  else if(connectivity_matrix[j*p.getSize() + i])
	    {
	      outfile << "PART_" << j << " -> PART_" << i <<";" << std::endl;
	    }
	}
      outfile << std::endl;
    }

  outfile << "}" << std::endl;


  outfile.close();


  outfile.open("file.connec");
  for(int i = 0 ; i < p.getSize() ; ++i)
    {
      for(int j = 0 ; j < p.getSize() ; ++j)
	{
	  outfile << connectivity_matrix[i*p.getSize() + j] << " ";
	}
      outfile << std::endl;
    }
  outfile.close();

  delete[] connectivity_matrix;
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
  popot::rng::Halton<Problem::nb_parameters>::init();

  std::ofstream outfile("PlotExampleMoving/fitness.data");
  
  // Let's create our swarm
  PSO pso;

  // Print the fitness of our particles
  printf("Before learning : \n");
  pso.print(1);
  std::cout << std::endl;

  double * bestpos = new double[2];

  // We now iterate the algorithm
  // We can iterate step by step
  for(int i = 0 ; i < 200 ; ++i)
    {
      save_function_values(i);
      save_particle_positions(i, pso);
      //save_connectivity(i,pso);
      
      bestpos[0] = pso.getBest()->getPosition(0);
      bestpos[1] = pso.getBest()->getPosition(1);
      outfile << pso.getBest()->getFitness() << " " << Problem::evaluate(bestpos) << std::endl;

      
      pso.step();

      std::cout << '\r' << std::setw(6) << std::setfill('0') << i << " " << pso.getBest()->getFitness() << std::setw(5) << std::setfill(' ') << ' ' << Problem::evaluate(bestpos) << std::flush;

      Problem::step();
      
    }
  outfile.close();

  // Or run the algorithm until the stopping criteria is met
  std::cout << "epoch : " << pso.epoch << std::endl;
  std::cout << "\n" << std::endl;

  // And display the best fitness we got
  printf("After learning : \n");
  pso.print(1);

  std::cout << " To get an animated gif of the behavior of the PSO, and if you used the source package, go into PlotExampleMoving and type make all . It requires convert and gnuplot" << std::endl;
  std::cout << " The fitnesses are saved in PlotExampleMoving/fitness.data" << std::endl;

  // Free the memory used by the problem (e.g. the bounds)
  popot::rng::Halton<Problem::nb_parameters>::free();
}
