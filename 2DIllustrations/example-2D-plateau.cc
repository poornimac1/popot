// In this example, we use a standard PSO on a predefined problem

// Compile with :
// g++ -o example-2D-plateau example-2D-plateau.cc `pkg-config --libs --cflags popot` -O3

#include "rng_generators.h"
typedef popot::rng::JKissRNG RNG_GENERATOR;

#include <stdio.h>
#include <cmath>
#include "popot.h"

// Define our problem
class Plateau
{
public:
  static const int nb_parameters = 2;
  static int count;
  
  static void init(void)
  {
    count = 0;
  }
  
  static void free(void)
  {
  }
  
  static double get_lbound(int index)
  {
    return -10.0;
  }

  static double get_ubound(int index)
  {
    return 10.0;
  }
  
  static bool stop(double fitness, int epoch)
  {
    return (epoch >= 100);
  }

  static double evaluate(void * x)
  {
    return .5;
  }
};
int Plateau::count;

typedef Plateau Problem;

// We use SPSO 2006 or 2011
typedef popot::PSO::SPSO2006::PSO<Problem>::Type PSO;
//typedef popot::PSO::SPSO2011::PSO<Problem>::Type PSO;


void presave_function_values()
{
  std::ofstream outfile("PlotExample2DPlateau/function_values.data");
  
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
  
  filename_plot << "PlotExample2DPlateau/" << root_filename.str()
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
	       << "set title \"PSO on plateau " << epoch << "\";" << std::endl
	       << "set palette defined ( 0 \"white\", 1 \"black\");"<< std::endl
	       << "set xlabel \"x\";"<< std::endl
	       << "set ylabel \"y\";" << std::endl
	       << "set cblabel \"z\";"<< std::endl
	       << "set view map;"<< std::endl
	       << "set pm3d at s;"<< std::endl
	       << "splot 'function_values.data' with pm3d notitle, \
  '-' with points notitle pt 7 ps 1.5 lc rgb \"red\""<< std::endl;


  // Put the position of the particles with z-value = ...
  PSO::ParticleType*it, *it_end;
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
  PSO::ParticleType* it, *it_end;
  
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
  
  filename << "PlotExample2DPlateau/" << root_filename.str()
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
  Problem::init();
  popot::rng::Halton<Problem::nb_parameters>::init();

  std::ofstream outfile("PlotExample2DPlateau/fitness.data");
  
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
      save_particle_positions(i, pso);
      save_connectivity(i,pso);
      outfile << pso.getBest()->getFitness() << std::endl;
      pso.step();

      std::cout << '\r' << std::setw(6) << std::setfill('0') << i << " " << pso.getBest()->getFitness() << std::setw(5) << std::setfill(' ') << ' ' << std::flush;
    }
  outfile.close();

  // Or run the algorithm until the stopping criteria is met
  std::cout << "epoch : " << pso.epoch << std::endl;
  std::cout << "\n" << std::endl;

  // And display the best fitness we got
  printf("After learning : \n");
  pso.print(1);

  std::cout << " To get an animated gif of the behavior of the PSO, and if you used the source package, go into PlotExample2DPlateau and type make all . It requires convert and gnuplot" << std::endl;
  std::cout << " The fitnesses are saved in PlotExample2DPlateau/fitness.data" << std::endl;

  // Free the memory used by the problem (e.g. the bounds)
  Problem::free();
  popot::rng::Halton<Problem::nb_parameters>::free();
}
