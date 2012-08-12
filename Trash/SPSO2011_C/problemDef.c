
struct problem problemDefBench(int functionCode)
{

  struct problem pb;
  int d;
  pb.function=functionCode;	
  // Default values
  // Can be modified below for each function
  pb.epsilon = 0.00000;	// Acceptable error (default). May be modified below
  pb.objective = 0;       // Objective value (default). May be modified below
  //	pb.constraintNb=0;
  pb.SS.normalise=0; // Set to a value x.
  
  switch(functionCode)
    {
    case 0:
      // Sphere
      pb.SS.D = 30;
      pb.epsilon = 1e-2;
      pb.objective = 0;
      pb.evalMax = 75000;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -100;
	  pb.SS.max[d] = 100;
	}
      break;
    case 1:
      // Griewank
      pb.SS.D = 30;
      pb.epsilon = 0.05;
      pb.objective = 0.0;
      pb.evalMax = 75000;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -600;
	  pb.SS.max[d] = 600;
	}
      break;
    case 2:
      // Rosenbrock
      pb.SS.D = 30;
      pb.epsilon =   100;
      pb.objective = 0;
      pb.evalMax = 75000;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -30;
	  pb.SS.max[d] = 30;
	}
      break;
    case 3:
      // Rastrigin
      pb.SS.D = 30;
      pb.epsilon =   50;
      pb.objective = 0;
      pb.evalMax = 75000;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -5.12;
	  pb.SS.max[d] = 5.12;
	}
      break;
    case 4:
      // Tripod
      pb.SS.D = 2;
      pb.epsilon =   1e-4;
      pb.objective = 0;
      pb.evalMax = 1e4;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -100;
	  pb.SS.max[d] = 100;
	}
      break;
    case 5:
      // Ackley
      pb.SS.D = 30;
      pb.epsilon =   0;
      pb.objective = 0;
      pb.evalMax = 80000;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -30;
	  pb.SS.max[d] = 30;
	}
      break;
    case 6:
      // Schwefel
      pb.SS.D = 30;
      pb.epsilon =   2569.5;
      pb.objective = -12569.5;
      pb.evalMax = 300000;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -500;
	  pb.SS.max[d] = 500;
	}
      break;
    case 7:
      // Schwefel 1.2
      pb.SS.D = 40;
      pb.epsilon =   0;
      pb.objective = 0;
      pb.evalMax = 40000;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -100;
	  pb.SS.max[d] = 100;
	}
      break;
    case 8:
      // Schwefel 2.22
      pb.SS.D = 30;
      pb.epsilon =   1e-4;
      pb.objective = 0;
      pb.evalMax = 100000;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -10;
	  pb.SS.max[d] = 10;
	}
      break;
    case 9:
      // Neumaier 3
      pb.SS.D = 40;
      pb.epsilon =   0;
      pb.objective = 0;
      pb.evalMax = 40000;
      // Boundaries
      for (d = 0; d < pb.SS.D; d++)
	{
	  pb.SS.min[d] = -pb.SS.D*pb.SS.D;
	  pb.SS.max[d] = pb.SS.D*pb.SS.D;
	}
      break;
    default:
      printf("Unknonw function !!! \n");
    }


  return pb;



}
