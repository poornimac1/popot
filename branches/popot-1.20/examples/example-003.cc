// In this example, we show how ones can define its own problem with a PSO
// The problem we try to solve is optimizing a Multilayer Perceptron to detect if a number is odd or even
// The input is represented by 7 segments indexed the following way
//         ----- 0 -----
//         |           |
//         5           1
//         |           |
//         ----- 6 -----
//         |           |
//         4           2
//         |           |
//         ----- 3 -----
// Therefore, we code :
//   = [ 6 5 4 3 2 1 0]  <- Index of the segments
// 0 = [ 0 1 1 1 1 1 1] ; 0x3F
// 1 = [ 0 0 0 0 1 1 0] ; 0x06
// 2 = [ 1 0 1 1 0 1 1] ; 0x5B
// 3 = [ 1 0 0 1 1 1 1] ; 0x4F
// 4 = [ 1 1 0 0 1 1 0] ; 0x66
// 5 = [ 1 1 0 1 1 0 1] ; 0x6D
// 6 = [ 1 1 1 1 1 0 1] ; 0x7D
// 7 = [ 0 0 0 0 1 1 1] ; 0x07
// 8 = [ 1 1 1 1 1 1 1] ; 0x7F
// 9 = [ 1 1 0 1 1 1 1] ; 0x6F
// In the code below, these are stored in inputs_mlp in the order : 0 1 2 3 4 5 6 ; 0 1 2 3 4 5 6 ....

#include <stdio.h>
#include <cmath>
#include <list>

// We first define the generator of random numbers
#include "rng_generators.h"
typedef popot::rng::CRNG RNG_GENERATOR;

#include "popot.h"


// Define our problem
// We define it as a MLP tested on our dataset
class MLPClassifier
{
public:
  static const int nb_inputs = 7;
  static const int nb_outputs = 1;
  static const int nb_hidden = 4;
  static const int nb_digits = 10;

  static const int nb_parameters = (nb_inputs+1)*nb_hidden + (nb_hidden+1)*nb_outputs;
  static int *hexa_codes;
  static int *inputs_mlp;

  static void init(void)
  {
    hexa_codes = new int[nb_digits];
    hexa_codes[0] = 0x3F;
    hexa_codes[1] = 0x06;
    hexa_codes[2] = 0x5B;
    hexa_codes[3] = 0x4F;
    hexa_codes[4] = 0x66;
    hexa_codes[5] = 0x6D;
    hexa_codes[6] = 0x7D;
    hexa_codes[7] = 0x07;
    hexa_codes[8] = 0x7F;
    hexa_codes[9] = 0x6F;

    inputs_mlp = new int[nb_inputs*nb_digits];
    int y;
    for(int x = 0 ; x < nb_digits ; ++x)
      {
	y = hexa_codes[x];
	for(int i = 0 ; i < nb_inputs ; ++i)
	  {
	    inputs_mlp[nb_inputs*x + i] = (y & 0x01);
	    y = y >> 1;
	  }
      }
  }

  static double get_lbound(int index)
  {  
    return -10;
  }

  static double get_ubound(int index)
  {    
    return 10;
  }

  static void display_digits(std::vector<int> &list_digits)
  {
    if(list_digits.size() == 0)
      {
	printf("Nothing to display .. \n");
	return;
      }
    // We want to display several digits on a row under their 7 segments format
    // We therefore display one after the other the rows for all the digits

    // Horizontal segments span 13 characters
    // Vertical segments span 3 characters
    int span_horiz = 13;
    int span_vert = 3;
    int space_digit = 3;

    // Segment 0
    for(int i = 0 ; i < list_digits.size() ; ++i)
      {
	for(int j = 0 ; j < span_horiz ; ++j)
	  {
	    if(inputs_mlp[list_digits[i]*nb_inputs+0])
	      printf("-");
	    else
	      printf(" ");
	  }
	for(int j = 0 ; j < space_digit ; ++j)
	  printf(" ");
      }
    printf("\n");
    // Segments 5 and 1
    for(int i = 0 ; i < span_vert ; ++i)
      {
	for(int j = 0 ; j < list_digits.size() ; ++j)
	  {
	    if(inputs_mlp[list_digits[j]*nb_inputs+5])
	      printf("|");
	    else
	      printf(" ");
	    for(int k = 0 ; k < span_horiz - 2 ; ++k)
	      printf(" ");
	    if(inputs_mlp[list_digits[j]*nb_inputs+1])
	      printf("|");
	    else
	      printf(" ");
	    for(int k = 0 ; k < space_digit ; ++k)
	      printf(" ");
	  }
	printf("\n");
      }

    // Segment 6
    for(int i = 0 ; i < list_digits.size() ; ++i)
      {
	for(int j = 0 ; j < span_horiz ; ++j)
	  {
	    if(inputs_mlp[list_digits[i]*nb_inputs+6])
	      printf("-");
	    else
	      printf(" ");
	  }
	for(int j = 0 ; j < space_digit ; ++j)
	  printf(" ");
      }
    printf("\n");
    // Segments 4 2
    for(int i = 0 ; i < span_vert ; ++i)
      {
	for(int j = 0 ; j < list_digits.size() ; ++j)
	  {
	    if(inputs_mlp[list_digits[j]*nb_inputs+4])
	      printf("|");
	    else
	      printf(" ");
	    for(int k = 0 ; k < span_horiz - 2 ; ++k)
	      printf(" ");
	    if(inputs_mlp[list_digits[j]*nb_inputs+2])
	      printf("|");
	    else
	      printf(" ");
	    for(int k = 0 ; k < space_digit ; ++k)
	      printf(" ");
	  }
	printf("\n");
      }
    // Segment 3
    for(int i = 0 ; i < list_digits.size() ; ++i)
      {
	for(int j = 0 ; j < span_horiz ; ++j)
	  {
	    if(inputs_mlp[list_digits[i]*nb_inputs+3])
	      printf("-");
	    else
	      printf(" ");
	  }
	for(int j = 0 ; j < space_digit ; ++j)
	  printf(" ");
      }

    printf("\n");
  }

  static void free(void)
  {
  }

  static double transfer_function(double x)
  {
    // We use a sigmoidal transfer function
    return 1.0 / (1.0 + exp(-x));
  }

  static double compute_output(int input_index, double * params)
  {
    double act_input[nb_inputs];
    double act_hidden[nb_hidden];
    double act_output;

    // Set up the input
    for(int i = 0 ; i < nb_inputs ; ++i)
      act_input[i] = inputs_mlp[nb_inputs * input_index + i];

    // Compute the activity of the hidden nodes
    int params_index = 0;
    for(int i = 0 ; i < nb_hidden ; ++i)
      {
	act_hidden[i] = params[params_index++];
	for(int j = 0 ; j < nb_inputs ; ++j)
	  {
	    act_hidden[i] += params[params_index++] * act_input[j];
	  }
	// Apply the transfer function
	act_hidden[i] = transfer_function(act_hidden[i]);
      }
    // Compute the activity of the output node
    act_output = params[params_index++];
    for(int j = 0 ; j < nb_hidden ; ++j)
      act_output += params[params_index++] * act_hidden[j];
    act_output = transfer_function(act_output);

    return act_output;
  }

  static bool stop(double fitness, int epoch)
  {
    return (fitness <= 1e-10) || (epoch >= 1000);
  }

  static double evaluate(double * params)
  {
    double fitness = 0.0;
    double act_output;
    // Test over all the inputs
    for(int k = 0 ; k < nb_digits ; ++k)
      {
	act_output = compute_output(k, params);

	fitness += pow(act_output - (k%2),2.0);
      }
    return fitness;
  }
};
int *MLPClassifier::hexa_codes;
int *MLPClassifier::inputs_mlp;

// Let's typedef the problem so that the following is identical to example-001
typedef MLPClassifier Problem;

// With the above definitions of the parameters, particles, topology and stopping criteria
// we can now set up our PSO algorithm
typedef popot::PSO::SPSO2011::PSO<Problem>::Type PSO;

// **************************************** //
// ************** Main ******************** //
// **************************************** //

int main(int argc, char* argv[]) {

  RNG_GENERATOR::rng_srand();

  // For testing
  std::vector<int> odd_numbers;
  std::vector<int> even_numbers;
  std::vector<int> unclassified_numbers;
  double * params = new double[Problem::nb_parameters];

  // Initialize our problem
  // this actually allocates memory and initializes the boundaries
  Problem::init();

  popot::rng::Halton<Problem::nb_parameters>::init();

  // Let's create our swarm
  PSO pso;

  // Let's generate the graph of the connections within the swarm
  pso.generateGraph("connections.dot");


  ////////////////////////////////////////:
  // Test before learning :
  printf("--------------------------------------------- \n Before learning : \n");
  printf("Best fitness : %f \n", pso.getBest()->getFitness());
  
  for(int i = 0 ; i < Problem::nb_parameters ; ++i)
    params[i] = pso.getBest()->getPosition(i);

  for(int i = 0 ; i < Problem::nb_digits; ++i)
    {
      if(Problem::compute_output(i, params) >= 0.85)
	even_numbers.push_back(i);
      else if(Problem::compute_output(i, params) <= 0.15)
	odd_numbers.push_back(i);
      else
	unclassified_numbers.push_back(i);
    }

  printf("I classified %i numbers as EVEN : \n", even_numbers.size());
  Problem::display_digits(even_numbers);
  printf("I classified %i numbers as ODD : \n", odd_numbers.size());
  Problem::display_digits(odd_numbers);
  printf("I was not able to classify %i numbers : \n", unclassified_numbers.size());
  Problem::display_digits(unclassified_numbers);  

  ////////////////////////////////////////:
  // We now iterate the algorithm
  pso.run();
  std::cout << "epoch : " << pso.epoch << std::endl;
  std::cout << "\n" << std::endl;

  ////////////////////////////////////////:
  // Test after learning :
  printf("--------------------------------------------- \n After learning : \n");
  printf("Best fitness : %f \n", pso.getBest()->getFitness());

  // Get the best parameters
  for(int i = 0 ; i < Problem::nb_parameters ; ++i)
    params[i] = pso.getBest()->getPosition(i);

  // And test them on our inputs
  odd_numbers.clear();
  even_numbers.clear();
  unclassified_numbers.clear();
  //for(int i = 0 ; i < Problem::nb_digits ; ++i)
  //    printf("%i : %f \n", i, Problem::compute_output(i, params));
  for(int i = 0 ; i < Problem::nb_digits; ++i)
    {
      if(Problem::compute_output(i, params) >= 0.85)
	even_numbers.push_back(i);
      else if(Problem::compute_output(i, params) <= 0.15)
	odd_numbers.push_back(i);
      else
	unclassified_numbers.push_back(i);
    }

  printf("I classified %i numbers as EVEN : \n", even_numbers.size());
  Problem::display_digits(even_numbers);
  printf("I classified %i numbers as ODD : \n", odd_numbers.size());
  Problem::display_digits(odd_numbers);
  printf("I was not able to classify %i numbers : \n", unclassified_numbers.size());
  Problem::display_digits(unclassified_numbers);

  // Free the memory used by the problem (e.g. the bounds)
  Problem::free();

  popot::rng::Halton<Problem::nb_parameters>::free();
}

