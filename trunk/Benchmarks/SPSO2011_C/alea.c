
double alea (double a, double b) 
{				// random number (uniform distribution) in  [a b]
	double r;
	int ir;
// Read from afile
		//r=randNumber[randRank]; r=a+r*(b-a);
		
	// Chaos
/*
	randRank=randRank+1; if (randRank>=randNbMax) {randRank=0;randChaos=0.02;};

	ir=floor(randChaos*10000); r=(double)ir/10000;
	
	randChaos=4*r*(1-r); 
	r=a+randChaos*(b-a);
*/
// KISS	
//r=a+(double)rand_kiss()*(b-a)/RAND_MAX_KISS;	

// Standard C
r=a+(double)rand()*(b-a)/RAND_MAX; 
	return r; 
}

//==================================================
void aleaIndex(int index[], int S)
{
	int indexTemp[S_max];
	int length;
	int rank;
	int s;
	int t;
	
	length=S;
			for (s=0;s<S;s++) indexTemp[s]=s; //=index[s];

			for (s=0;s<S;s++)
			{
				rank=alea_integer(0,length-1);
				index[s]=indexTemp[rank];
				
				if (rank<length-1)	// Compact
				{
					for (t=rank;t<length-1;t++)
						indexTemp[t]=indexTemp[t+1];
				}					
				length=length-1;
			}
}
// ===========================================================
int alea_integer (int a, int b) 
{				// Integer random number in [a b]
	int ir;
	double r;

	//r = alea (0, 1);
	//ir = (int) (a + r * (b + 1 - a));
	//if (ir > b)	ir = b;
	
	r=alea((double)a,(double)b+1); ir=floor(r);
	return ir;  
}
// ===========================================================
double alea_stable (double alpha, double beta, double nu, double delta)
{	// CMS algorithm (J.M. Chambers, C.L. Mallows and B.W. Stuck)
	// alpha 	:= 	stability parameter. Tail thickness (the higher the thinner)
	//						 Must be in ]0,2]. 
	//						For normal distribution, alpha=2
	// beta		:=	skweness.  0 => symmetric distribution
	// nu			:=	dispersion. Must be >0. For normal distribution, nu=standard dev.
	// delta	:=	mean (mesure of centrality)
	double betaPrime; 
	double d;
	double eps;
	double kAlpha;
	double phi0, phi;
	double r;
	double s;
	double t1,t2,t3;
	double tau;
	double temp;
	double w;
	double z;

	if(alpha<0 || alpha>2)
	{ printf("\n alpha %f ",alpha);
		ERROR("alea_levy. alpha must be in ]0,2]");
	}
	if(nu<0)
	{
		printf("\n nu %f ",nu);
		ERROR("alea_levy. nu must be positive");
	}
	//--------------------------------------------
	// Define k(α) = 1 − |1−α|. Thus if α≤1 then k(α)=α and if α≥1 then k(α)=2−α.   
	if(alpha<1) kAlpha=alpha; else kAlpha=2-alpha;

	// Compute φ0 = −½β(k(α)/α).  
	phi0=0.5*beta*(kAlpha/alpha);

	temp=tan(alpha*phi0);

	// Transform β to β' by β' = β if α=1 and β'= −tan(½π(1−α))tan(αφ0) otherwise.

	if(fabs(alpha-1)<zero) betaPrime=beta;
	else
		betaPrime=-tan(0.5*pi*(1-alpha))*temp;

	// Generate a random variable u uniformly distributed on the interval [0, 1] and compute φ = π(u−½). 		
	phi=pi*(alea(0,1)-0.5);

	//Compute ε=1−α and then τ = −εtan(αφ0) 		
	eps=1-alpha; 
	tau=-eps*temp;

	//Compute tan(½φ), tan(½εφ) and tan(½εφ)/(½εφ). 
	t1=tan(0.5*phi);
	t2=tan(0.5*eps*phi);
	t3=2*t2/(eps*phi);

	// Generate a random variable v which has a uniform distribution 
	//	on the interval [0,1] and then compute w=−ln(v)
	w=-log(alea(0,1));

	// Compute z = (cos(εφ)−tan(αφ0)sin(εφ)/(wcos(φ))
	z=cos(eps*phi)-tan(alpha*phi0)*sin(eps*phi)/(w*cos(phi));

	// Compute d = z^ε/α /ε 
	temp=pow(z,eps/alpha);
	d=temp/eps;

	// Compute s = tan(αφ0) + z^ε/α (sin(αφ)−tan(αφ0)cos(αφ))/cos(φ) 
	s=tan(alpha*phi0) + temp*(sin(alpha*phi) - tan(alpha*phi0)*cos(alpha*phi))/cos(phi);	

	// Multiply by the dispersion, and add the mean
	r=s*nu + delta;

	return r;
}
// ===========================================================
double alea_normal (double mean, double std_dev) 
{ 
	// Use the polar form of the Box-Muller transformation to obtain a pseudo
	// random number from a Gaussian distribution 

	double x1, x2, w, y1;  

	do  
	{
		x1 = 2.0 * alea (0, 1) - 1.0;
		x2 = 2.0 * alea (0, 1) - 1.0;
		w = x1 * x1 + x2 * x2;     
	}
	while (w >= 1.0);

	w = sqrt (-2.0 * log (w) / w);
	y1 = x1 * w;

	if(alea(0,1)<0.5) y1=-y1; 
	y1 = y1 * std_dev + mean;
	return y1;  
}

// ===========================================================
struct vector alea_sphere(int D, double radius, int distrib, double mean, double sigma)
{
	/*  ******* Random point in a hypersphere ********
	 Maurice Clerc 2003-07-11
	 Last update: 2011-01-01

	 Put  a random point inside the hypersphere S(center 0, radius 1), 
	 or on its surface
		 */

	int 	j;
	double   length;
	double      pw;
	double      r;
	struct	vector	x;

	x.size=D;
	pw=1/(double)D;

	// ----------------------------------- Step 1.  Direction
	length=0;
	for (j=0;j<D;j++)
	{
		if(distrib<2) x.v[j]=alea_normal(0,1); // Gaussian (Box-Muller method)
		else x.v[j]=alea_stable(2,0,1,0); // Gaussian (CMS method)
		
		length=length+  x.v[j]*x.v[j];
	}

	length=sqrt(length);
//printf("\nalea 189 length %f",length);
	//----------------------------------- Step 2. 

	switch(distrib) 
	{
		default: // Random radius, uniform distribution
			r=alea(0,1); 
		break;

		case 1: //Gaussian (Box-Muller)
			r=fabs(alea_normal(mean,sigma)); 
		break; 

		case 2: //Gaussian (CMS)
			r=fabs(alea_stable(2,0,sigma,mean)); 
		break; 
		
		case 3: //TEST
			//r=fabs(alea_stable(0.2,0,sigma,alea(0,1))); 
			r=fabs(alea_stable(0.2,0,sigma,0.5*(alea(0, 1)+alea(0,1)))); 
		break; 
		
		case 4: // Slash distribution
		r=fabs(alea_normal(mean,sigma))/fabs(alea_normal(0,1));
		break;

		case 99: // Constant radius, for some specific uses (random ON the sphere)
			r=1; 
		break;
	} 
//printf("\nalea 219 r %f",r);
	for (j=0;j<D;j++)
	{
		x.v[j]=radius*r*x.v[j]/length;
	}
	return x;
}
//================================================== KISS
/*
 A good pseudo-random numbers generator

 The idea is to use simple, fast, individually promising
 generators to get a composite that will be fast, easy to code
 have a very long period and pass all the tests put to it.
 The three components of KISS are
 x(n)=a*x(n-1)+1 mod 2^32
 y(n)=y(n-1)(I+L^13)(I+R^17)(I+L^5),
 z(n)=2*z(n-1)+z(n-2) +carry mod 2^32
 The y's are a shift register sequence on 32bit binary vectors
 period 2^32-1;
 The z's are a simple multiply-with-carry sequence with period
 2^63+2^32-1.  The period of KISS is thus
 2^32*(2^32-1)*(2^63+2^32-1) > 2^127
 */

static ulong kiss_x; //= 1;
static ulong kiss_y; //= 2;
static ulong kiss_z; //= 4;
static ulong kiss_w; //= 8;
static ulong kiss_carry = 0;
static ulong kiss_k;
static ulong kiss_m;



void seed_rand_kiss(ulong seed) 
{
	kiss_x = seed | 1;
	kiss_y = seed | 2;
	kiss_z = seed | 4;
	kiss_w = seed | 8;
	kiss_carry = 0;
}

ulong rand_kiss() 
{
	kiss_x = kiss_x * 69069 + 1;
	kiss_y ^= kiss_y << 13;
	kiss_y ^= kiss_y >> 17;
	kiss_y ^= kiss_y << 5;
	kiss_k = (kiss_z >> 2) + (kiss_w >> 3) + (kiss_carry >> 2);
	kiss_m = kiss_w + kiss_w + kiss_z + kiss_carry;
	kiss_z = kiss_w;
	kiss_w = kiss_m;
	kiss_carry = kiss_k >> 30;
	//printf("\n%f ",(double) (kiss_x + kiss_y + kiss_w));
	return kiss_x + kiss_y + kiss_w;
}
//==========================================================================    
struct vectorList quasiRand(int D, int nRand, int option)
{   
	/*
	 Generate nRand vectors of size D that are the coordinates of
	 nRand quasi-random points in [0,1]^D

		 For C under Linux, you need to use the GSL library:
#include <gsl/gsl_qrng.h> // Do not forget to link to gsl and gslcblast

		 */

	int i;	
	gsl_qrng *qrng_q;    
	struct vectorList qRand;
	
	switch(option)
	{
		case 1: // Sobol
			if(D>40) 
		{
			printf("\nThe embedded Sobol sequences generator can not be used for dimensions greater than 40");
			printf("\nYou should use Halton sequences (for dimensions up to 1229)");
			ERROR("\n I stop here");
			
		}
		qrng_q=gsl_qrng_alloc (sobol, D);
		break;
		
		default: // Halton
			if(D>1229) 
		{printf("\nThe embedded Halton sequences generator can not be used for dimensions greater than 1229");
			printf("\nSorry");
			ERROR("\n I stop here");
		}
		qrng_q=gsl_qrng_alloc (halton, D);
		break;
	}

	gsl_qrng * q = qrng_q;
	
	for (i = 0; i < nRand; i++)
	{
		gsl_qrng_get (q, qRand.V[i].v);       
	}

	gsl_qrng_free (q); 

	return qRand;
}
