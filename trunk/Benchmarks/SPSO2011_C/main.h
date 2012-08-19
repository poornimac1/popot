#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include <time.h>

#define	DMax 100 			// 114		// Max number of dimensions of the search space
#define fMax 6				// Max number of constraints +1
#define funcMax 13		// Max number of functions sequentially solved
#define R_max 100000		// Max number of runs
#define	S_max 100 		// Max swarm size
#define zero  1.0e-40	// To avoid numerical instabilities on some machines
#define infinity 1.e+100

# define ulong unsigned long // To generate pseudo-random numbers with KISS
# define RAND_MAX_KISS ((unsigned long) 4294967295) //  Needs ISO C90

// For quasi-random numbers generation (option)
#include <gsl/gsl_qrng.h> // Do not forget to link to gsl and gslcblas
#define halton gsl_qrng_halton
#define sobol gsl_qrng_sobol

# define randNbMax 200 //8000 // When using a list of previously generated random numbers
											// DMax * S_max

#define ERROR(a) {printf("\nERROR: %s\n", a);exit(-1);}

// Structures
struct quantum 
{
	double q[DMax];
	int size; 
};

struct SS
{ 
	int D;							// Dimension of the search space (D-rectangle)
	double max[DMax];		// Maximum value on each dimension 
	double min[DMax]; 	// Minimum value on each dimension
	struct quantum q;		// Quantisation step size. 0 => continuous problem
	double normalise; 	//  >0  if normalisation needed on [0,normalise]^D
	int quantisation; 	// Flag. Set to 1 if quantisation needed
};

struct param 
{
	int BW[3];			// For "bells and whistles" that are not really
									//			part of the standard
	double c;				// Acceleration coefficient
	int confin;			// Confinement option
	int distrib;		// Distribution option
	int K;					// Max number of particles informed by a given one
	double mean;		// Mean in case of Gaussian distribution
	double p;				// Probability threshold for random topology	
									// (is actually computed as p(S,K) )
	int S;					// Swarm size. Optionnaly (BW[0]=1) only the _mean_ swarm size
									// is given by the user
	double sigma;		// variance in case of Gaussian distribution
	int topology;		// option to modify the topology of the information links
	double w;				// Inertia weight
};

struct fitness 
{
	int size; 
	double f[fMax];
};

struct position 
{ 

struct vector 
{  
	int size;  
	double v[DMax]; 
};

struct vectorList
{
	int size;
	struct vector V[randNbMax];
};

struct problem 
{ 
	//int constraintNb;						// Number of constraints
	double epsilon; 					// Admissible error
	int evalMax; 							// Maximum number of fitness evaluations
	int function; 						// Function code
	double objective; 				// Objective value
	struct SS SS;							// Search space
};

struct swarm 
{ 
	int best; 					// rank of the best particle
	struct position P[S_max];	// Previous best positions found by each particle
	int S; 						// Swarm size 
	struct vector V[S_max];	// Velocities
	struct position X[S_max];	// Positions 
};

struct result 
{  
	double nEval; 		// Number of evaluations  
	struct swarm SW;	// Final swarm
	double error;		// Numerical result of the run
};


// --------------- Sub-programs (lexicographic order)

// In alea.c
double alea (double a, double b);
void aleaIndex(int index[], int S);
int alea_integer (int a, int b);
double alea_normal (double mean, double stdev);
struct vector alea_sphere(int D, double radius, int option, double mean,double sigma);
double alea_stable (double alpha, double beta, double nu, double delta);
struct vectorList quasiRand(int D, int nRand, int option);

// For the pseudo-random number generator KISS
ulong	rand_kiss(); 
void	seed_rand_kiss(ulong seed); 

// In lennard_jones.c (included in perf.c)
double lennard_jones (struct position x); 

// In perf.c
struct fitness constraint(struct position x, int functCode);
double perf_bench (struct position x, int function, double objective);
double perf (struct position x, int function,struct SS SS, double objective);	

// In problemDef.c
struct problem problemDefBench(int functionCode);
struct problem problemDef(int functionCode);
	// See also cec2005.c, cec2005data.c, cec2005pb.c
	
// In PSO.c
struct result PSO ( struct param param, struct problem problem);
struct position quantis (struct position x, struct SS SS);

// In tools.c
static int compareDoubles (void const *a, void const *b);
double distanceL(struct position x1, struct position x2, double L);
double Gamma(double u);
double max(double a, double b);
double min(double a, double b);
int sign (double x);


// --------------- Global variables	
long double	E ;  // exp(1). Useful for some test functions
double errMax,errMin;
long double	pi ; // Useful for some test functions

double randNumber[randNbMax];
int randRank;
double randChaos;

// For Network problem
int bcsNb;
int btsNb;

// --------------- Files
//FILE * f_run;
//FILE * f_synth;
//FILE * f_rand;