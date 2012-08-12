/*
  Standard PSO 2011

  Contact for remarks, suggestions etc.:
  pso@writeme.com

  For more details, see ReadMe.txt		
*/

// g++ -o main main.c -O3 `pkg-config --libs --cflags gsl`

#include "main.h"

// =================================================
int main (int argc, char * argv[]) 
{ 

  if(argc != 2)
    {
      printf("Usage : %s <func> \n", argv[0]);
      return -1;
    }
  int funcNumber = atoi(argv[1]);

  struct position bestBest; // Best position over all runs
  int d;			// Current dimension
  double D;
  double error;			// Current error
  double errorMean;		// Average error
  double errorMin;		// Best result over all runs
  double errorMeanBest[R_max]; 
  double evalMean;		// Mean number of evaluations
  int func[funcMax]; // List of functions (codes) to optimise
  int indFunc;

  int nFailure;		// Number of unsuccessful runs
  double logProgressMean;
  struct param param;
  struct problem pb; 

  int run, runMax; 
  struct result result; 

  time_t seconds; 
	
  double Smean;
  double success[funcMax];
  double successRate;
  int t;
  double variance;
  float z;
  double zz;

  E = exp ((long double) 1);  
  pi = acos ((long double) -1); 
  errMax=0;

  //------------------------------------------------ PARAMETERS
  // Bells and Whistles
  // Not really part of the standard
  // May improve the performance (not always)
  // Sometimes not well mathematically founded (rules of thumbs)
  // * => suggested value
				
  param.BW[0]=0; 	//	0 => same swarm size for each run 	
  //	1 => random swarm size around the given mean
									
  param.BW[1]=0; 	//*	0 => when P=G, use the "standard" method	
  // 	1 => when P=G, a specific probabilistic method 
  //	2 => when P=G, a more conservative method
  //  3 => when P=G, just look around
  // 4 =>  different weights for X, P, G  (TEST)
									
  param.BW[2]=0;	// Randomness options
  // 0 => pseudo-random number generator KISS
  // 1 => quasi-random numbers for initialisation, Sobol sequences
  //* 2 => quasi-random numbers for initialisation, Halton sequences
  // <integer> nnn =>  **** TO DO. 
  //--------
	
  param.confin=0; 	// 0 => keep inside the search space (supposed to be a D-rectangle)
  // 1 => no confinement 
  //   WARNING: may be very slow (and bad) for discrete problems
										
  param.distrib=0; 	// 0 => uniform
  // 1 => Gaussian (Box-Muller method)
  // 2 =>	Gaussian (CMS method)
  // 3 => Other stable (CMS, experimental parameters)
  // 4 => Slash distribution (Gaussian BM/Gaussian BM)
  // Useful only if param.distrib>0;
  //param.mean=0.5; //Default: 0.5. For some functions 0 is better, though 
  //	Example: shifted Rosenbrock (code 102)
  //param.sigma=1./12; // Default: 1./12 (standard deviation of U(0,1))			

  Smean=40; //Swarm size or Mean swarm size (if BW[0]=1). Default: 40 

  param.K=3; 	// Parameter to compute the probability p for a particle to be an
  // external informant. You may also directly define p (see below),
  // but K is about the mean number of the informants of a particle. 
  // Default: 3

  // Confidence coefficients. Default:
  param.w = 1. / (2 * log ((double) 2)); // 0.721
  param.c = 0.5 + log ((double) 2); // 1.193
  //param.w=0.7215; // A bit beyond of the theoretical convergence area
  //param.w=-0.55; param.c=1.8;
  param.topology = 0; // 0 => information links as in SPSO 2007 (quasi-random)
  // 1 => variable random ring (EXPERIMENTAL)

  // ----------------------------------------------- PROBLEM
  // Functions to optimise

  runMax = 1000; // Numbers of runs

  // Some information
  //printf ("\n Function %i ", funcNumber);

  // Define the problem
  pb=problemDefBench(funcNumber);			

  // ----------------------------------------------- RUNS	
  errorMean = 0;	    
  evalMean = 0;	    
  nFailure = 0;	
  D=pb.SS.D;

  int new_neigh = 0;

  srand(time(NULL));
  //seed_rand_kiss(1234567890);

  /*
    seconds=time(NULL); // Initialise the RNG KISS more randomly
    printf("\n time %ld",seconds);
    seed_rand_kiss(time(NULL)); 
  */
  for (t=0;t<10000;t++) 
    {
      zz=alea(0,1); // "Warm up" the RNG
      //printf("\n%f",zz);
    }
		
  for (run = 0; run < runMax; run++)  
    {
      //printf("%i \n", run);

      param.S=Smean; // Constant swarm size
      param.p=1-pow(1-1./param.S,param.K); 

      result = PSO (param, pb);
      error = result.error;

      if (error > pb.epsilon) // Failure
	nFailure = nFailure + 1;	

      // Memorize the best (useful if more than one run)
      if(run==0) bestBest=result.SW.P[result.SW.best];
      else
	if(error<bestBest.f) bestBest=result.SW.P[result.SW.best];

      // Result display
      errorMean=errorMean+error;

      // Compute/save some statistical information
      if (run == 0)
	errorMin = error;
      else if (error < errorMin)
	errorMin = error;
	
      new_neigh += result.new_neigh;
			
      evalMean = evalMean + result.nEval;	
      errorMeanBest[run] = error;
      logProgressMean  = logProgressMean - log(error);		

      //printf("%i %f %i  \n", run, result.nEval, result.new_neigh);
    }		// End loop on "run"
      

  // ---------------------END 
		
  // Display some statistical information
  evalMean = evalMean / (double) runMax;   
  errorMean = errorMean / (double) runMax;
  logProgressMean = logProgressMean/(double) runMax;

  printf ("\n Eval. (mean)= %f", evalMean);	
  printf ("\n Error (mean) = %e", errorMean);
  printf ("\n New neigh : %i , mean = %f \n", new_neigh, new_neigh/double(runMax));
  // Variance
  variance = 0;

  for (run = 0; run < runMax; run++)
    variance = variance + pow (errorMeanBest[run] - errorMean, 2);

  variance = sqrt (variance / runMax);	    
  printf ("\n Std. dev. %e", variance); 
  printf("\n Log_progress (mean) = %f", logProgressMean);	

  // Success rate and minimum value
  printf("\n Failure(s) %i  ",nFailure);

  successRate = 100 * (1 - nFailure / (double) runMax);			
  printf ("Success rate = %.2f%%", successRate);
  success[0]=successRate;

  printf ("\n Best min value = %1.20e", errorMin);

  printf("\n errMax : %f",errMax);
  // Repeat informations
  printf("\n---------");
  printf("\n Function(s): %i ", funcNumber);
  printf("\n Confinement: "); 
  if(param.confin==0) printf("YES"); else printf("NO");
  printf("\n Distribution: ");
  switch(param.distrib)
    {
    case 0: printf(" uniform"); break;
    case 1: printf(" Gaussian (%f,%f), Box-Muller",param.mean,param.sigma); break;
    case 2: printf(" Gaussian (%f,%f), CMS",param.mean,param.sigma); break;
    case 3: printf(" Stable (%f,%f)",param.mean,param.sigma); break;
    case 4: printf(" Slash (%f,%f)",param.mean,param.sigma); break;
    }

  printf("\n BW = (%i, %i, %i )",param.BW[0],param.BW[1],param.BW[2]);
  printf("\n Swarm size: ");
  if(param.BW[0]==0) printf("%i",(int)Smean); else printf(" mean %i",(int)Smean);
  printf("\n K = %i",param.K);
  printf("\n w = %f",param.w);
  printf("\n c = %f",param.c);
  printf("\n");
  return 0; // End of main program
}
// ===============================================================
#include "alea.c"
#include "perf.c"
#include "problemDef.c"
#include "PSO.c"
#include "tools.c"


