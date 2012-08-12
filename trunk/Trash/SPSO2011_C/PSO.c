struct result PSO (struct param param, struct problem pb) 
{  
  int d; 
  double error;   
  double errorPrev;
  int g;  
  struct position Gr;
  struct vector GX={0};   
  int index[param.S];     
  int initLinks;	// Flag to (re)init or not the information links
  int iter; 		// Iteration number (time step)
  int iterBegin;
  int LINKS[param.S][param.S];	// Information links
  int m; 	
  int noStop;
  double out;
  double p;
  struct vector PX;	
  struct vectorList qRand;
  struct result R;
  double rad;
  int s0, s,s1; 
  struct vector V1,V2; 
  double w1,w2,w3;
  double xMax,xMin;
  //struct position XPrev;
  double zz;

  PX.size=pb.SS.D; 
  GX.size=pb.SS.D; 
  V1.size=pb.SS.D;
  V2.size=pb.SS.D;
  Gr.size=pb.SS.D;
  // -----------------------------------------------------
  // INITIALISATION
  p=param.p; // Probability threshold for random topology

  R.SW.S = param.S; // Size of the current swarm

  // Position and velocity
  for (s = 0; s < R.SW.S; s++)   
    {
      R.SW.X[s].size = pb.SS.D;
      R.SW.V[s].size = pb.SS.D;
    }

  switch(param.BW[2])
    {
    default:
      for (s = 0; s < R.SW.S; s++)   
	{
	  for (d = 0; d < pb.SS.D; d++)  
	    {  
	      if(pb.SS.normalise==0) { xMin=pb.SS.min[d]; xMax=pb.SS.max[d];}
	      R.SW.X[s].x[d] = alea (xMin,xMax);		
	      R.SW.V[s].v[d] = alea( xMin-R.SW.X[s].x[d], xMax-R.SW.X[s].x[d] ); 
	      // So that  xMin<= x+v <= xMax
	    }
	}
      break;
	
    case 1:
    case 2:
      qRand=quasiRand(pb.SS.D,R.SW.S,param.BW[2]);
	
      for (s = 0; s < R.SW.S; s++)   
	{
	  for (d = 0; d < pb.SS.D; d++)  
	    {  
	      if(pb.SS.normalise==0) { xMin=pb.SS.min[d]; xMax=pb.SS.max[d];}
	      R.SW.X[s].x[d] = xMin+(xMax-xMin)*qRand.V[s].v[d];		
	      R.SW.V[s].v[d] = alea( xMin-R.SW.X[s].x[d], xMax-R.SW.X[s].x[d] ); 
	    }
	}
	
      break;
	
    }

  // First evaluations
  errMax=0;
  errMin=infinity;
  for (s = 0; s < R.SW.S; s++) 
    {	
      R.SW.X[s].f = perf_bench (R.SW.X[s], pb.function,pb.objective);

      R.SW.P[s] = R.SW.X[s];	// Best position = current one
    }



  //
  // If the number max of evaluations is smaller than 
  // the swarm size, just keep evalMax particles, and finish
  if (R.SW.S>pb.evalMax) R.SW.S=pb.evalMax;	
  R.nEval = R.SW.S;

  // Find the best
  R.SW.best = 0;

  errorPrev =R.SW.P[R.SW.best].f; // "distance" to the wanted f value (objective)

  for (s = 1; s < R.SW.S; s++)     
    {
      zz=R.SW.P[s].f;
      if (zz < errorPrev)
	{
	  R.SW.best = s;
	  errorPrev=zz;
	} 
    }

  // ---------------------------------------------- ITERATIONS	
  noStop = 0;	
  error=errorPrev;	
  iter=0; iterBegin=0;
  initLinks = 1;		// So that information links will be initialized

  int new_neigh = 0;

  // Each particle informs itself
  for (m = 0; m < R.SW.S; m++) LINKS[m][m] = 1;	

  while (noStop == 0) 
    {
      iter=iter+1;
		
      aleaIndex(index, R.SW.S); // Random numbering of the particles
				
      if (initLinks==1)	// Modify topology
	{
	  // Who informs who, at random
	  for (s = 0; s < R.SW.S; s++)
	    {				
	      for (m = 0; m < R.SW.S; m++)
		{		    
		  if(m==s) continue;		

		  if (alea (0, 1)<p)
		    { 
		      LINKS[m][s] = 1;	// Probabilistic method
		    }
		  else LINKS[m][s] = 0;
		}	
	    }					 
	}

      // Loop on particles
      for (s0 = 0; s0 < R.SW.S; s0++)	// For each particle ...
	{	
	  s=index[s0];

	  // Find the best informant			
	  g = 0;	
	  for (m = 0; m < R.SW.S; m++) 
	    {	    
	      if (LINKS[m][s] == 1 && R.SW.P[m].f < R.SW.P[g].f)
		g = m;
	    }	

	  //.. compute the new velocity, and move

	  // Exploration tendency
	  for (d = 0; d < pb.SS.D; d++)
	    {
	      R.SW.V[s].v[d]=param.w *R.SW.V[s].v[d];
	    }

	  // Prepare Exploitation tendency  p-x
	  for (d = 0; d < pb.SS.D; d++)
	    {
	      PX.v[d]= R.SW.P[s].x[d] - R.SW.X[s].x[d];
	    }					

	  if(g!=s) // If the particle is not its own local best, prepare g-x
	    {
	      for (d = 0; d < pb.SS.D; d++) 
		GX.v[d]= R.SW.P[g].x[d] - R.SW.X[s].x[d];
	    }


	  // Gravity centre Gr
	  w1=1; w2=1; w3=1; 

	  if(g==s) w3=0; 

	  zz=w1+w2+w3;
	  w1=w1/zz; w2=w2/zz; w3=w3/zz;

	  for (d = 0; d < pb.SS.D; d++) 
	    {
	      Gr.x[d]=w1*R.SW.X[s].x[d] + w2*(R.SW.X[s].x[d] + param.c*PX.v[d]);
	
	      if(g!=s) 
		Gr.x[d]=Gr.x[d]+ w3*(R.SW.X[s].x[d] +param.c*GX.v[d]) ;	

	      V1.v[d]= Gr.x[d]-R.SW.X[s].x[d]; // Vector X-Gr
	    }

	  // Random point around
	  rad=distanceL(R.SW.X[s],Gr,2);

	  V2=alea_sphere(pb.SS.D,rad,param.distrib, param.mean,param.sigma); 

	  // New "velocity"

	  for (d = 0; d < pb.SS.D; d++) 
	    R.SW.V[s].v[d]=R.SW.V[s].v[d]+V1.v[d] + V2.v[d]; // New "velocity"

	  // New position

	  for (d = 0; d < pb.SS.D; d++) 
	    R.SW.X[s].x[d] = R.SW.X[s].x[d] + R.SW.V[s].v[d];			

	  if (R.nEval >= pb.evalMax)  goto end;

	  // Confinement			
	  //	  out=0;

	  for (d = 0; d < pb.SS.D; d++)
	    {	
	      xMin=pb.SS.min[d]; 
	      xMax=pb.SS.max[d];

	      if (R.SW.X[s].x[d] < xMin)
		{	
		  R.SW.X[s].x[d] = xMin;
		  R.SW.V[s].v[d] = -0.5*R.SW.V[s].v[d];
		  //out=1;
		}
	      else
		{
		  if (R.SW.X[s].x[d] > xMax)
		    {			
		      R.SW.X[s].x[d] = xMax;
		      R.SW.V[s].v[d] = -0.5*R.SW.V[s].v[d];
		      //out=1;
		    }
		}				
	    }	
	  
	  // Evaluation
	  //R.SW.X[s].f =perf(R.SW.X[s],pb.function, pb.SS,pb.objective);
	  R.SW.X[s].f =perf_bench(R.SW.X[s],pb.function, pb.objective);
	  R.nEval = R.nEval + 1;				
	  
	  // ... update the best previous position		
	  if (R.SW.X[s].f < R.SW.P[s].f)	// Improvement
	    {		 		
	      R.SW.P[s] = R.SW.X[s];
	      
	      // ... update the best of the bests
	      if (R.SW.P[s].f < R.SW.P[R.SW.best].f)
		{		
		  R.SW.best = s;			
		}			
	    }
	}			// End of "for (s0=0 ...  "	

      // Check if finished
      error = R.SW.P[R.SW.best].f;

      if (error < errorPrev)	// Improvement of the global best
	initLinks = 0;							
      else			// No global improvement
	{
	  initLinks = 1;	// Information links will be	reinitialized	
	  new_neigh++;
	}
      
      errorPrev = error;
    end:

      if (error > pb.epsilon && R.nEval < pb.evalMax)
	noStop = 0;	// Won't stop
      else
	noStop = 1;	// Will stop

    } // End of "while nostop ...

  // printf( "\n and the winner is ... %i", R.SW.best );			
  R.error = error;
  R.new_neigh = new_neigh;
  return R;  
}


