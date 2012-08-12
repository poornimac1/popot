
double perf_bench(struct position x, int function, double objective)
{
  struct position xs = x;

  double f;
  int d,d1;
  double sum1, sum2;
  double prod;
  double p;
  double xd;
  double t0, t1, tt, x1, x2;
  int k;
  double y;

  switch ( function )
    {
    case 0: // Parabola (Sphere)
      f = 0;
      p = 0; // Shift
      for ( d = 0; d < xs.size; d++ )
        {
          xd = xs.x[d] - p;
          f = f + xd * xd;
        }
      break;
    case 1: // Griewank
      f = 0;
      p = 1;
      for ( d = 0; d < xs.size; d++ )
        {
          xd = xs.x[d];
          f = f + xd * xd;
          p = p * cos( xd / sqrt( d + 1 ) );
        }
      f = f / 4000 - p + 1;
      break;
    case 2: // Rosenbrock
      f = 0;
      t0 = xs.x[0];
      for ( d = 1; d < xs.size; d++ )
        {
          t1 = xs.x[d];
          tt = 1 - t0;
          f += tt * tt;
          tt = t1 - t0 * t0;
          f += 100 * tt * tt;
          t0 = t1;
        }
      break;
    case 3: // Rastrigin. Minimum value 0. Solution (0,0 ...0)
      k = 10;
      f = 0;
      for ( d = 0; d < xs.size; d++ )
        {
          xd = xs.x[d];
          f += xd * xd - k * cos( 2 * pi * xd );
        }
      f += xs.size * k;
      break;
    case 4: // Tripod
      // Search [-100, 100]^2. min 0 on (0  -50)
      x1 = xs.x[0];
      x2 = xs.x[1];
      if ( x2 < 0 )
        {
          f = fabs( x1 ) + fabs( x2 + 50 );
        }
      else
        {
          if ( x1 < 0 )
            f = 1 + fabs( x1 + 50 ) + fabs( x2 - 50 );
          else
            f = 2 + fabs( x1 - 50 ) + fabs( x2 - 50 );
        }
      break;
    case 5: // Ackley
      sum1 = 0;
      sum2 = 0;
      for ( d = 0; d < xs.size; d++ )
        {
          xd = xs.x[d];
          sum1 += xd * xd;
          sum2 += cos( 2 * pi * xd );
        }
      y = xs.size;
      f = ( -20 * exp( -0.2 * sqrt( sum1 / y ) ) - exp( sum2 / y ) + 20 + E );
      break;
    case 6: // Schwefel
      f = 0;
      for ( d = 0; d < xs.size; d++ )
        {
          xd = xs.x[d];
          f -= xd * sin(sqrt(fabs(xd)));
        }
      break;
    case 7: // Schwefel 1.2
      f = 0;
      sum1 = 0;
      for ( d = 0; d < xs.size; d++ )
        {
	  sum1 = 0;
	  for( d1 = 0 ; d1 <= d ; d1++)
	    {
	      xd = xs.x[d1];
	      sum1 += xd;
	    }
	  f += sum1 * sum1;
        }
      break;
    case 8: // Schwefel 2.22
      f = 0;
      sum1 = 0;
      prod = 1;
      for ( d = 0; d < xs.size; d++ )
        {
	  xd = xs.x[d];
	  sum1 += fabs(xd);
	  prod *= fabs(xd);
	}
      f = sum1 + prod;
      break;
    case 9: // Neumaier 3
      f = 0;
      for ( d = 0; d < xs.size; d++ )
        {
	  xd = xs.x[d];
	  f += (xd - 1)*(xd-1);
	}
      for(d = 1 ; d < xs.size ; d++)
	f += xs.x[d] * xs.x[d-1];
      break;
    }


  
  return fabs(f-objective);
}

