#ifndef POPOT_INDIVIDUALS_H
#define POPOT_INDIVIDUALS_H

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "maths.h"
#include "neighborhood.h"
#include "exceptions.h"
#include "tools.h"
#include "initializers.h"

// For particle param values see
// For standard PSO 
// Clerc, M., et al.: Standard PSO 2007
// Bratton, D., Kennedy, J.: Defining a Standard for Particle Swarm Optimization.

namespace popot
{

  /**
   * Basic type representing a position in a search space
   * @brief The template parameters are the type of values (e.g. double, or bool)
   *        and the dimension of the vector
   */
  template< typename VALUE_TYPE, int SIZE>
    class Vector
  {
  protected:
    VALUE_TYPE * values;
    int dimension;
      
  public:

    /**
     * Constructor
     */
    Vector(void)
      {
	values = new VALUE_TYPE[SIZE];
	dimension = SIZE;
      }

    /**
     * Copy constructor
     */
    Vector(const Vector& other)
      {
	values = new VALUE_TYPE[SIZE];
	dimension = SIZE;
	for(int i = 0 ; i < dimension ; ++i)
	  values[i] = other.values[i];
      }
      

    /**
     * Copy Constructor
     */
    Vector & operator=(const Vector& other)
      {
	for(int i = 0 ; i < dimension ; ++i)
	  values[i] = other.values[i];

	return *this;
      }


    /**
     * Destructor
     */
    virtual ~Vector(void)
      {
	delete[] values;
      }


    /**
     * Setter for a value
     * @param index position
     * @param value value to set
     */
    void setValueAt(int index, VALUE_TYPE value)
    {
      values[index] = value;
    }

    /**
     * Returns a component of the vector
     * @param index position
     */
    VALUE_TYPE getValueAt(int index) const
    {
      if(index >= 0 && index < SIZE)
	return values[index];
      else
	throw popot::Exception::IndexOutOfRange(index, SIZE);
    }

    /**
     * Returns the size of the vector
     */
    int size(void) const
    {
      return SIZE;
    }

    /**
     * Returns a pointer to the raw data
     */
    VALUE_TYPE * getValuesPtr(void)
    {
      return values;
    }

    /**
     * Sum of two vectors
     */
    Vector operator+ (const Vector& b) const
    {
      // No need to check the dimensions
      // they would be of different types
      Vector res;
      for(int i = 0 ; i < SIZE ; ++i)
	res.setValueAt(i, getValueAt(i) + b.getValueAt(i));
      return res;
    }

    /**
     * In-place Sum of two vectors
     */  
    Vector& operator+=(const Vector& b)
      {
	for(int i = 0 ; i < SIZE ; ++i)
	  setValueAt(i, getValueAt(i) + b.getValueAt(i));
	return *this;
      }
      
    /**
     * Right Multiplication by a constant
     */ 
    Vector operator* (double b) const
      {
	Vector res;
	for(int i = 0 ; i < SIZE ; ++i)
	  res.setValueAt(i, b * getValueAt(i));
	return res;
      }
      
    /**
     * Left-Multiplication by a constant
     */ 
    friend Vector operator*(double a, const Vector&b) 
      {
	return b * a;
      }

    /**
     * Method for displaying a vector, fills in the stream
     */ 
    virtual void print(std::ostream & os) const
    {
      os << "[";
      for(int i = 0 ; i < size() ; ++i)
	{
	  os << getValueAt(i);
	  if(i != size() - 1)
	    os << ";";
	}
      os << "]";
    }

    /**
     * Serialization operator, for display, based on the print method
     * to modify the display of derived class, just overload the print method
     */
    friend std::ostream & operator <<(std::ostream & os, const Vector &v)
    {
      v.print(os);
      return os;
    }
  };


  namespace PSO
  {
    namespace particle
    {

      /**
       * BaseParticle is a vector with a fitness and an associated Problem
       */
      template< typename PROBLEM , typename POSITION_INITIALIZER>
	class BaseParticle : public Vector<double, PROBLEM::nb_parameters>
	{
	private:
	  typedef Vector<double, PROBLEM::nb_parameters> TSuper;

	protected:
	  double _fitness;
	  static const int _dimension = PROBLEM::nb_parameters;

	public:	 

	  typedef PROBLEM Problem;

	BaseParticle(void) : TSuper()
	    {
	    }
	  
	  virtual ~BaseParticle(void) 
	    {}


	  /**
	   * Random intialization of the position followed by a fitness evaluation
	   */
	  virtual void init(void)
	  {
	    // Initialize the position
	    initPosition();

	    // Evaluate the fitness
	    evaluateFitness();
	  }

	  /**
	   * Initialization of the position
	   */
	  virtual void initPosition()
	  {
	    double bounds[2];
	    for(int i = 0 ; i < _dimension ; ++i)
	      {
		bounds[0] = PROBLEM::get_lbound(i);
		bounds[1] = PROBLEM::get_ubound(i);
		// The position is initialized within the bounds
		setPosition(i, POSITION_INITIALIZER::init(bounds));
	      }
	  }

	  /**
	   * Getter on the position
	   */
	  double getPosition(int i)
	  {
	    if(i >= 0 && i < _dimension)
	      return this->getValueAt(i);
	    else
	      throw popot::Exception::IndexOutOfRange(i, _dimension);
	  }
	  /**
	   * Setter on the position
	   */
	  void setPosition(int i, double val)
	  {
	    if(i >= 0 && i < _dimension)
	      this->setValueAt(i, val);
	    else
	      throw popot::Exception::IndexOutOfRange(i, _dimension);
	  }

	  /**
	   * Ensures the position is within the boundaries
	   */
	  virtual void confine(void)
	  {
	    for(int i = 0 ; i < _dimension ; ++i)
	      {
		if(getPosition(i) < PROBLEM::get_lbound(i))
		  setPosition(i, PROBLEM::get_lbound(i));
		else if(getPosition(i) > PROBLEM::get_ubound(i))
		  setPosition(i, PROBLEM::get_ubound(i));
	      }
	  }

	  /**
	   * Returns the currently known fitness
	   */
	  double getFitness(void) const
	  {
	    return _fitness;
	  }

	  /**
	   * Set the fitness
	   */
	  void setFitness(double f)
	  {
	    _fitness = f;
	  }
	
	  /**
	   * Recompute the fitness
	   */
	  virtual double evaluateFitness(void)
	  {
	    _fitness = PROBLEM::evaluate(this->getValuesPtr());
	    return _fitness;
	  }
	  /**
	   * Comparison of the fitness of two particles
	   */
	  bool operator<(const BaseParticle &p) const
	  {
	    return (compare(p) < 0);
	  }
	  /**
	   * Comparison of the fitness of two particles p1.compare(p2)
	   * @return -1 if p1.f < p2.f
	   * @return 1 if p1.f > p2.f
	   * @return 0 otherwise
	   */
	  virtual int compare(BaseParticle * p)
	  {
	    double myfitness = getFitness();
	    double otherfitness = p->getFitness();
	    if(myfitness < otherfitness)
	      return -1;
	    else if(myfitness > otherfitness)
	      return 1;
	    else
	      return 0;
	  }  

	  virtual void print(std::ostream & os) const
	  {
	    TSuper::print(os);
	    os << " ; Fitness : " << this->getFitness();
	  }
	
	};


      // Particle introduces the notion
      // of neighborhood and velocity
      // and best position
      template< typename PROBLEM, typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER>
	class Particle : public BaseParticle<PROBLEM, POSITION_INITIALIZER>
	{
	private:
	  typedef BaseParticle<PROBLEM, POSITION_INITIALIZER> TSuper;
	  typedef Particle<PROBLEM, POSITION_INITIALIZER,VELOCITY_INITIALIZER > ThisParticleType;

	public:
	  typedef BaseParticle<PROBLEM, POSITION_INITIALIZER> BestType;
	  typedef popot::PSO::neighborhood::Neighborhood< ThisParticleType > NeighborhoodType;
	
	
	protected:
	  double *_velocity;
	  BestType _best_position;
	  NeighborhoodType *_neighborhood;
	
	public:

	Particle(void) : TSuper()
	    {
	      _velocity = new double[TSuper::_dimension];
	      memset(_velocity, 0, TSuper::_dimension * sizeof(double));

	      // Create an empty neighborhood
	      _neighborhood = new NeighborhoodType();
	    }

	  virtual ~Particle(void)
	    {
	      delete[] _velocity;
	      delete _neighborhood;
	    }

	  /**
	   * Initialization of the position, velocity and fitness of the particle and set the best_position as the current position
	   */
	  virtual void init(void)
	  {
	    // Initialize the position
	    // and computes the fitness
	    TSuper::init();

	    // Initialize the velocity
	    initVelocity();

	    // Set the best particle to the current position
	    _best_position = *this;
	  }


	  /**
	   * Initialization of the velocity
	   */
	  virtual void initVelocity(void)
	  {
	    double params[3];
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      {
		params[0] = PROBLEM::get_lbound(i);
		params[1] = PROBLEM::get_ubound(i);
		params[2] = this->getPosition(i);
		setVelocity(i, VELOCITY_INITIALIZER::init(params));
	      }
	  }

	  /**
	   * Getter on the velocity
	   */
	  double getVelocity(int i)
	  {
	    if(i >= 0 && i < TSuper::_dimension)
	      return _velocity[i];
	    else
	      throw popot::Exception::IndexOutOfRange(i, TSuper::_dimension);
	  }

	  /**
	   * Setter on the velocity
	   */
	  void setVelocity(int i, double value)
	  {
	    if(i >= 0 && i < TSuper::_dimension)
	      _velocity[i] = value;
	    else
	      throw popot::Exception::IndexOutOfRange(i, TSuper::_dimension);
	  }

	  /**
	   * Bounds the velocity and position of the particle
	   */
	  virtual void confine(void)
	  {
	    // In case the position is out of the bounds
	    // we reset the velocities
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      {
		if((this->getPosition(i) < PROBLEM::get_lbound(i)))
		  {
		    this->setPosition(i, PROBLEM::get_lbound(i));
		    this->setVelocity(i,0);
		  }
		else if(this->getPosition(i) > PROBLEM::get_ubound(i))
		  {		  
		    this->setPosition(i, PROBLEM::get_ubound(i));
		    this->setVelocity(i,0);
		  }
	      }
	  }

	  /**
	   * Bounds the velocity of the particle
	   */
	  virtual void boundVelocity(void)
	  {
	    return;
	  }

	  /**
	   * Updates the best position
	   */
	  virtual void updateBestPosition(void)
	  {
	    // Update the best position the particle ever had
	    if(this->compare(&_best_position) < 0)
	      _best_position = *this;
	  }

	  /**
	   * Updates the velocity of the particle
	   */
	  virtual void updatePosition(void)
	  {
	    std::cout << "No updatePosition rule is provided for a particle of type Particle " << std::endl;
	  }

	  /**
	   * Updates the velocity of the particle
	   */
	  virtual void updateVelocity(void)
	  {
	    std::cout << "No updateVelocity rule is provided for a particle of type Particle " << std::endl;
	  }

	  virtual BestType * getBestPosition(void)
	  {
	    return &_best_position;
	  }

	  /**
	   * Returns a pointer to the neighborhood of the particle
	   */
	  NeighborhoodType* getNeighborhood(void)
	  {
	    return _neighborhood;
	  }

	  /**
	   * Set the pointer to the neighborhood
	   */
	  void setNeighborhood(NeighborhoodType * neighborhood)
	  {
	    _neighborhood = neighborhood;
	  }

	  virtual void print(std::ostream & os) const
	  {
	    TSuper::print(os);
	    os << "; Best position : " << _best_position;
	  }

	};

      /**
       * Definition of a Standard PSO 2006 particle
       * @brief The template parameters is the problem you want to solve and a class providing the parameters w and c
       */
      template< typename PROBLEM, typename PARTICLE_PARAMS>
	class SPSO2006Particle : public Particle<PROBLEM, popot::PSO::initializer::PositionUniformRandom, popot::PSO::initializer::VelocityHalfDiff >
	{
	  typedef Particle<PROBLEM, popot::PSO::initializer::PositionUniformRandom, popot::PSO::initializer::VelocityHalfDiff > TSuper;

	public:

	SPSO2006Particle(void) : TSuper()
	    {
	    }

	  /**
	   * Destructor
	   */
	  virtual ~SPSO2006Particle(void)
	    {
	    }

	  /**
	   * Update the position of the particle
	   */
	  virtual void updatePosition(void)
	  {
	    // Here it is simply : p_{k+1} = p_k + v_k
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i));
	  }

	  /**
	   * Updates the velocity of the particle
	   */
	  virtual void updateVelocity(void)
	  {
	    // The update of the velocity is done according to the equation :
	    // v = w * v + c r1 (best_p - p) + c r2 (best_g - p)
	    // with :
	    // r_p, r_g two random real numbers in [0.0,1.0]
	    // w, c1, c2 : user defined parameters
	    // best_p : the best position the particle ever had
	    // best_g : the best position the neighborhood ever had
	    double r1,r2;
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      {
		r1 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c());
		r2 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c());
		this->setVelocity(i, PARTICLE_PARAMS::w() * this->getVelocity(i)
				  + r1 * (this->getBestPosition()->getPosition(i) - this->getPosition(i))
				  + r2 * (this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i)));
	      }
	  }
	};

      /**
       * Definition of a Standard PSO 2007 particle
       * @brief The template parameters is the problem you want to solve and a class providing the parameters w and c
       */
      template< typename PROBLEM, typename PARTICLE_PARAMS>
	class SPSO2007Particle : public Particle<PROBLEM, popot::PSO::initializer::PositionUniformRandom, popot::PSO::initializer::VelocityHalfDiff >
	{
	  typedef Particle<PROBLEM, popot::PSO::initializer::PositionUniformRandom, popot::PSO::initializer::VelocityHalfDiff > TSuper;

	public:

	SPSO2007Particle(void) : TSuper()
	    {
	    }

	  /**
	   * Destructor
	   */
	  virtual ~SPSO2007Particle(void)
	    {
	    }

	  /**
	   * Update the position of the particle
	   */
	  virtual void updatePosition(void)
	  {
	    // Here it is simply : p_{k+1} = p_k + v_k
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i));
	  }

	  /**
	   * Updates the velocity of the particle
	   */
	  virtual void updateVelocity(void)
	  {
	    // The update of the velocity is done according to the equation :
	    // v = w * v + r1 (best_p - p) + r2 (best_g - p)
	    // with :
	    // r_1, r_2 two random real numbers in [0.0,c]
	    // w, c : user defined parameters
	    // best_p : the best position the particle ever had
	    // best_g : the best position the neighborhood ever had

	    bool li_equals_pi = (this->getBestPosition() == this->getNeighborhood()->getBest());
	    double r1,r2;
	    if(li_equals_pi)
	      {
		for(int i = 0 ; i < TSuper::_dimension ; ++i)
		  {
		    r1 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c());
		    this->setVelocity(i, PARTICLE_PARAMS::w() * this->getVelocity(i)
				      + r1 * (this->getBestPosition()->getPosition(i) - this->getPosition(i)));
		  }
	      }
	    else
	      {
		for(int i = 0 ; i < TSuper::_dimension ; ++i)
		  {
		    r1 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c());
		    r2 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c());
		    this->setVelocity(i, PARTICLE_PARAMS::w() * this->getVelocity(i)
				      + r1 * (this->getBestPosition()->getPosition(i) - this->getPosition(i))
				      + r2 * (this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i)));
		  }
	      }
	  }
	};


      /**
       * Definition of a Standard PSO 2011 particle
       * @brief The template parameters is the problem you want to solve,
       *  and the parameters (inertia, accelaration) of the velocity update rule
       * Compared to the base type BaseParticle, this type adds the some specific position and velocity update rules
       * These come from the Standard PSO 2011 
       */
      template< typename PROBLEM, typename PARTICLE_PARAMS>
	class SPSO2011Particle : public Particle<PROBLEM, popot::PSO::initializer::PositionUniformRandom, popot::PSO::initializer::VelocitySPSO2011 >
	{
	  typedef Particle<PROBLEM, popot::PSO::initializer::PositionUniformRandom, popot::PSO::initializer::VelocitySPSO2011 > TSuper;

	  double * xpi;
	  double * p1;
	  double * p2;
	  double * gr;
	public:

	SPSO2011Particle(void) : TSuper()
	    {
	      xpi = new double[TSuper::_dimension];
	      p1 = new double[TSuper::_dimension];
	      p2 = new double[TSuper::_dimension];
	      gr = new double[TSuper::_dimension];
	    }

	  /**
	   * Destructor
	   */
	  virtual ~SPSO2011Particle(void)
	    {
	      delete[] xpi;
	      delete[] p1;
	      delete[] p2;
	      delete[] gr;
	    }

	  /**
	   * Update the position of the particle
	   */
	  virtual void updatePosition(void)
	  {
	    // Here it is simply : p_{k+1} = p_k + v_k
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i));
	  }

	  /**
	   * Updates the velocity of the particle
	   */
	  virtual void updateVelocity(void)
	  {
	    // The update of the velocity is done according to the equation :
	    // v_i(t+1) = w v_i(t) + x'_i(t) - x_i(t)
	    // where x'_i is normally sampled from the hypersphere (Gi, Ri)
	    // with Gi the center and Ri the radius and
	    // Gi = 1/3 ( x_i + (x_i + c(p_i - x_i)) + (x_i + c(l_i - x_i)))
	    //    = 1/3 ( x_i +          p1          +            p2     )))
	    // or 
	    // Gi = 1/2 (x_i + (x_i + c(p_i - x_i))) if p_i == l_i
	    //    = 1/2 (x_i +          p1         )
	    // Ri = ||G_i - x_i||

	    // We first check if the local best and personal best are identical by computing the norm of the difference
	    // between both positions
	    // This can be done with pointer comparison
	    // as the neighborhood holds a pointer to the best personal best
	    bool li_equals_pi = (this->getBestPosition() == this->getNeighborhood()->getBest());

	    // First position
	    // p1 = xi + c * (pi - xi)
	    // where pi is the personal best position
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      p1[i] = this->getPosition(i) + PARTICLE_PARAMS::c() *(this->getBestPosition()->getPosition(i) - this->getPosition(i));
	  
	    // Second position
	    // p2 = xi + c * (li - xi)
	    // where li is the local best position (within the neighborhood)
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      p2[i] = this->getPosition(i) + PARTICLE_PARAMS::c() *(this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i));
	  
	    // Compute the gravity center of p1, p2 and xi
	    if(!li_equals_pi)
	      {
		// We here consider p1, p2 and xi
		for(int i = 0 ; i < TSuper::_dimension ; ++i)
		  gr[i] = 1.0/3.0 * (this->getPosition(i) + p1[i] + p2[i]);
	      }
	    else
	      {
		// We here consider only p1 (or p2, since they are equal) and xi
		for(int i = 0 ; i < TSuper::_dimension ; ++i)
		  gr[i] = 1.0/2.0 * (this->getPosition(i) + p1[i]);
	      }

	    // And the radius of the hypersphere
	    double ri = 0.0;
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      ri += (gr[i] - this->getPosition(i))*(gr[i] - this->getPosition(i));
	    ri = sqrt(ri);

	    // Compute the auxiliary position x'_i randomly within the hypersphere (gr, ri);
	    // To uniformely sample from the hypersphere we uniformely sample a direction
	    double norm = 0.0;
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      {
		xpi[i] = popot::math::normal(0.0,1.0);
		norm += xpi[i] * xpi[i];
	      }
	    norm = sqrt(norm);
	  
	    // And then scale by a random radius
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      xpi[i] =  gr[i] + popot::math::uniform_random(0.0,1.0) * ri * xpi[i] / norm;

	    // And then update the velocity
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      this->setVelocity(i, PARTICLE_PARAMS::w() * this->getVelocity(i)
				+ xpi[i] - this->getPosition(i));
	  }


	  /**
	   * Bounds the velocity and position of the particle
	   */
	  /*
	  virtual void confine(void)
	  {
	    // In case the position is out of the bounds
	    // we reset the velocities
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      {
		if((this->getPosition(i) < PROBLEM::get_lbound(i)))
		  {
		    this->setPosition(i, PROBLEM::get_lbound(i));
		    this->setVelocity(i,-0.5*this->getVelocity(i));
		  }
		else if(this->getPosition(i) > PROBLEM::get_ubound(i))
		  {		  
		    this->setPosition(i, PROBLEM::get_ubound(i));
		    this->setVelocity(i,-0.5*this->getVelocity(i));
		  }
	      }
	  }
	  */

	};

      /* /\** */
      /*  * Definition of a particle with stochastic fitness */
      /*  * @brief The template parameters is the problem you want to solve, */
      /*  *  and the parameters (inertia, accelaration) of the velocity update rule */
      /*  * It is a standard PSO, the only difference is that the fitness of the best position is reevaluated before being updated */
      /*  *\/ */
      /* template< typename PROBLEM, typename PARTICLE_PARAMS> */
      /* class StochasticTraditionalParticle : public TraditionalParticle<PROBLEM, PARTICLE_PARAMS> */
      /* { */
      /* 	typedef TraditionalParticle<PROBLEM, PARTICLE_PARAMS> TSuper; */

      /* public: */

      /* 	StochasticTraditionalParticle(void) : TSuper() */
      /* 	{ */
      /* 	} */

      /* 	/\** */
      /* 	 * Destructor */
      /* 	 *\/ */
      /* 	virtual ~StochasticTraditionalParticle(void) */
      /* 	{ */
      /* 	} */

      /* 	/\** */
      /* 	 * Updates the best position */
      /* 	 *\/ */
      /* 	virtual void updateBestPosition(void) */
      /* 	{ */
      /* 	  // Before comparing the current and best positions */
      /* 	  // we reevaluate the fitness of the best */
      /* 	  this->getBestPosition()->evaluateFitness(); */

      /* 	  // Update the best position the particle ever had */
      /* 	  TSuper::updateBestPosition(); */
      /* 	} */

      /* }; */

      /* /\** */
      /*  * Definition of a particle with Dynamic Constriction */
      /*  * @brief The template parameters is the problem you want to solve, */
      /*  *  and the parameters (inertia, accelaration) of the velocity update rule */
      /*  * When updating the best position of the particle, we reevaluate its fitness to handle some stochasticity in the fitness */
      /*  *\/ */
      /* template< typename PROBLEM, typename PARTICLE_PARAMS> */
      /* class DynamicConstrictionParticle : public Particle<PROBLEM> */
      /* { */
      /* 	typedef Particle<PROBLEM> TSuper; */

      /* public: */

      /* 	DynamicConstrictionParticle(void) : TSuper() */
      /* 	{ */
      /* 	} */

      /* 	/\** */
      /* 	 * Destructor */
      /* 	 *\/ */
      /* 	virtual ~DynamicConstrictionParticle(void) */
      /* 	{ */
      /* 	} */

      /* 	/\** */
      /* 	 * Update the position of the particle */
      /* 	 *\/ */
      /* 	virtual void updatePosition(void) */
      /* 	{ */
      /* 	  // Here it is simply : p_{k+1} = p_k + v_k */
      /* 	  for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 	    this->setPosition(i, this->getPosition(i) + this->getVelocity(i)); */
      /* 	} */

      /* 	/\** */
      /* 	 * Updates the velocity of the particle */
      /* 	 *\/ */
      /* 	virtual void updateVelocity(void) */
      /* 	{ */
      /* 	  double r1,r2; */
      /* 	  double xi,phi1,phi2,phi,kpa; */
      /* 	  kpa = 0.8; */
      /* 	  for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 	    { */
      /* 	      r1 = popot::math::uniform_random(0.0,1.0); */
      /* 	      r2 = popot::math::uniform_random(0.0,1.0); */
      /* 	      phi1 = r1 * PARTICLE_PARAMS::c1(); */
      /* 	      phi2 = r2 * PARTICLE_PARAMS::c2(); */
      /* 	      phi = phi1+phi2; */
      /* 	      if(phi <= 4) */
      /* 		xi = sqrt(kpa); */
      /* 	      else */
      /* 		xi = 2.0 * kpa / (phi-2.0 +sqrt(phi*(phi-4.0))); */
      /* 	      this->setVelocity(i, xi*( this->getVelocity(i) */
      /* 					+ phi1 * (this->getBestPosition()->getPosition(i) - this->getPosition(i)) */
      /* 					+ phi2 * (this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i)))); */
      /* 	    } */
      /* 	} */
      /* }; */


      /**
       * Definition of the Barebone particles
       * @brief The template parameters is the problem you want to solve,
       *  and the parameters (inertia, accelaration) of the velocity update rule
       */
      template< typename PROBLEM, typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER>
	class BareboneParticle : public Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER>
	{
	  typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> TSuper;

	public:

      	BareboneParticle(void) : TSuper()
	    {
	    }

	  /**
	   * Destructor
	   */
	  virtual ~BareboneParticle(void)
	    {
	    }

	  /**
	   * Update the position of the particle
	   */
	  virtual void updatePosition(void)
	  {
	    // Here it is simply : p_{k+1} = v_{k+1}
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      this->setPosition(i , this->getVelocity(i));
	  }

	  /**
	   * Updates the velocity of the particle
	   */
	  virtual void updateVelocity(void)
	  {
	    // The update of the velocity is done according to the equation :
	    // v_{k+1} = N( (best_p+best_g)/2), sigma)
	    // with :
	    // sigma : the variance defined as |best_p_ij - best_g_ij| with
	    // best_p : the best position the particle ever had
	    // best_g : the best position the neighborhood ever had

	    double mean,var;
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      {
		mean = 0.5*(this->getBestPosition()->getPosition(i) + this->getNeighborhood()->getBest()->getPosition(i));
		var = fabs(this->getBestPosition()->getPosition(i) - this->getNeighborhood()->getBest()->getPosition(i));
		this->setVelocity(i, popot::math::normal(mean,var));
	      }
	  }
	};

      /**
       * Definition of a modified barebone particle
       * @brief The template parameters is the problem you want to solve,
       *  and the parameters (inertia, accelaration) of the velocity update rule
       * In this implementation, the center of the normal distribution is the mean between the personal best
       * and local best positions (other variations exist, using a random personal best in the neighbor instead
       * of the local best for example)
       */
      template< typename PROBLEM , typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER>
	class ModifiedBareboneParticle : public Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER>
	{
	  typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> TSuper;
	
	public:

      	ModifiedBareboneParticle(void) : TSuper()
	    {
	    }

	  /**
	   * Destructor
	   */
	  virtual ~ModifiedBareboneParticle(void)
	    {
	    }

	  /**
	   * Update the position of the particle
	   */
	  virtual void updatePosition(void)
	  {
	    // Here it is simply : p_{k+1} = v_{k+1}
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      this->setPosition(i, this->getVelocity(i));
	  }

	  /**
	   * Updates the velocity of the particle
	   */
	  virtual void updateVelocity(void)
	  {
	    // The update of the velocity is done according to the equation :
	    // v_{k+1} = N( (best_p+best_g)/2), sigma)
	    // with :
	    // sigma : the variance defined as |best_p_ij - best_g_ij| with
	    // best_p : the best position the particle ever had
	    // best_g : the best position the neighborhood ever had
	    double mean,var;
	        
	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
	      {
		if(popot::math::uniform_random(0.0,1.0) < 0.5)
		  {
		    this->setVelocity(i, this->getBestPosition()->getPosition(i));
		  }
		else
		  {
		    mean = 0.5*(this->getBestPosition()->getPosition(i) + this->getNeighborhood()->getBest()->getPosition(i));
		    var = fabs(this->getBestPosition()->getPosition(i) - this->getNeighborhood()->getBest()->getPosition(i));
		    this->setVelocity(i, popot::math::normal(mean,var));
		  }
	      }
	  }
	}; // Modified Barebone

      /**
       * Definition of a fully informed particle
       * @brief The template parameters is the problem you want to solve,
       *  and the parameters (k, phi) of the velocity update rule
       */
      /* template< typename PARAMS, typename PROBLEM, typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER> */
      /* 	class FullyInformedParticle : public Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> */
      /* { */
      /* 	typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> TSuper; */

      /* public: */
      /* 	FullyInformedParticle(void) : TSuper() */
      /* 	{} */

      /* 	virtual ~FullyInformedParticle(void) */
      /* 	{ */
      /* 	} */

      /* 	/\** */
      /* 	 * Overloaded update velocity */
      /* 	 *\/  */
      /* 	virtual void updatePosition(void) */
      /* 	{ */
      /* 	  // Here it is simply : p_{k+1} = p_k + v_k */
      /* 	  for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 	    this->setPosition(i, this->getPosition(i) + this->getVelocity(i));	   */
      /* 	} */

      /* 	/\** */
      /* 	 * Overloaded update velocity */
      /* 	 *\/  */
      /* 	virtual void updateVelocity(void) */
      /* 	{ */
      /* 	  // Two possibilities are considered, */
      /* 	  // See "Fundamentals of computational swarm intelligence", p. 185 */
	  
      /* 	  // 1) The contribution of the neighborhood */
      /* 	  // is a mean velocity computed with the best personal */
      /* 	  // of all my neighboors. */
      /* 	  typename TSuper::NeighborhoodType * myneigh = this->getNeighborhood(); */
      /* 	  double r = 0.0; */
      /* 	  double r_sum = 0.0; */
      /* 	  double dv = 0.0; */
      /* 	  double kappa = 2 * PARAMS::k() / fabs(2.0 - PARAMS::phi() - sqrt(PARAMS::phi() * PARAMS::phi() - 4.0 * PARAMS::phi())); */

      /* 	  for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 	    { */
      /* 	      // Compute the mean velocity contribution */
      /* 	      dv = 0.0; */
      /* 	      for(int j = 0 ; j < myneigh->size() ; ++j) */
      /* 		{ */
      /* 		  r = popot::math::uniform_random(0.0, PARAMS::phi()); */
      /* 		  dv += r * (myneigh->get(j)->getPosition(i) - this->getPosition(i)); */
      /* 		} */
      /* 	      dv /= double(myneigh->size()); */
      /* 	      this->setVelocity(i, kappa * (this->getVelocity(i) + dv)); */
      /* 	    } */
      /* 	} */
      /* }; // Fully informed particle */

      /* /\** */
      /*  * Definition of a heterogeneous Particle */
      /*  * @brief The template parameters is the problem you want to solve, */
      /*  *  and the parameters (inertia, accelaration, ..) of the different velocity update rule */
      /*  *\/ */
      /* template< typename PROBLEM, typename PARTICLE_PARAMS, int FITNESS_HISTORY_SIZE> */
      /* class HeterogeneousParticle : public Particle<PROBLEM> */
      /* { */
      /* private: */
      /* 	typedef Particle<PROBLEM> TSuper; */
      /* 	popot::tools::FIFO<FITNESS_HISTORY_SIZE> _fitness_history; */

      /* 	// In the paper of Engelbrecht, 5 different update rules are considered */
      /* 	// - Traditional PSO */
      /* 	// - Cognitive only PSO */
      /* 	// - Social only PSO */
      /* 	// - Barebone PSO */
      /* 	// - Modified barebone */
      /* 	typedef enum { */
      /* 	  TRADITIONAL = 0, */
      /* 	  COGNITIVE = 1, */
      /* 	  SOCIAL = 2, */
      /* 	  BAREBONE = 3, */
      /* 	  MODIFIED_BAREBONE = 4, */
      /* 	  NB_RULES = 5 } RuleMode; */

      /* 	RuleMode _mode; */

      /* public: */
      /* 	HeterogeneousParticle(void) : TSuper() */
      /* 	{ */
      /* 	} */

      /* 	virtual ~HeterogeneousParticle(void) */
      /* 	{ */
      /* 	} */

      /* 	virtual void init(void) */
      /* 	{ */
      /* 	  TSuper::init(); */
      /* 	  // Random initialization of the update mode */
      /* 	  _mode = RuleMode((int)(popot::math::uniform_random(0.0,NB_RULES))); */
      /* 	} */

      /* 	/\** */
      /* 	 * Recompute the fitness and stores it in the FIFO */
      /* 	 *\/ */
      /* 	virtual double evaluateFitness(void) */
      /* 	{ */
      /* 	  _fitness_history.insert(TSuper::evaluateFitness()); */
      /* 	  if(_fitness_history.size() == FITNESS_HISTORY_SIZE) */
      /* 	    if(_fitness_history.variance() <= 1e-2) */
      /* 	      switchBehavior(); */
      /* 	} */

      /* 	void switchBehavior(void) */
      /* 	{ */
      /* 	  //std::cout << "behavior change" << std::endl; */
      /* 	  // Random initialization of the update mode */
      /* 	  _mode = RuleMode((int)(popot::math::uniform_random(0.0,NB_RULES))); */
      /* 	  _fitness_history.clear(); */
      /* 	} */

      /* 	virtual void updatePosition(void) */
      /* 	{ */
      /* 	  switch(_mode) */
      /* 	    { */
      /* 	    case TRADITIONAL: */
      /* 	    case COGNITIVE: */
      /* 	    case SOCIAL: */
      /* 	      { */
      /* 		// Here it is simply : p_{k+1} = p_k + v_k */
      /* 		for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 		  this->setPosition(i, this->getPosition(i) + this->getVelocity(i)); */
      /* 		break; */
      /* 	      } */
      /* 	    case BAREBONE: */
      /* 	    case MODIFIED_BAREBONE: */
      /* 	      { */
      /* 		// Here it is simply : p_{k+1} = p_k + v_k */
      /* 		for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 		  this->setPosition(i, this->getVelocity(i)); */
      /* 		break; */
      /* 	      } */
      /* 	    } */
      /* 	} */

      /* 	virtual void updateVelocity(void) */
      /* 	{ */
      /* 	  double r1,r2; */
      /* 	  double mean,var; */
      /* 	  switch(_mode) */
      /* 	    { */
      /* 	    case TRADITIONAL: */
      /* 	      for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 		{ */
      /* 		  r1 = popot::math::uniform_random(0.0,1.0); */
      /* 		  r2 = popot::math::uniform_random(0.0,1.0); */
      /* 		  this->setVelocity(i, PARTICLE_PARAMS::w_trad() * this->getVelocity(i) */
      /* 				    + r1 * PARTICLE_PARAMS::c1_trad() * (this->getBestPosition()->getPosition(i) - this->getPosition(i)) */
      /* 				    + r2 * PARTICLE_PARAMS::c2_trad() * (this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i))); */
      /* 		} */
      /* 	      break; */
      /* 	    case COGNITIVE: */
      /* 	      for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 		{ */
      /* 		  r1 = popot::math::uniform_random(0.0,1.0); */
      /* 		  this->setVelocity(i, PARTICLE_PARAMS::w_cognitive() * this->getVelocity(i)+ r1 * PARTICLE_PARAMS::c1_cognitive() * (this->getBestPosition()->getPosition(i) - this->getPosition(i))); */
      /* 		} */
      /* 	      break; */
      /* 	    case SOCIAL: */
      /* 	      for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 		{ */
      /* 		  r2 = popot::math::uniform_random(0.0,1.0); */
      /* 		  this->setVelocity(i, PARTICLE_PARAMS::w_social() * this->getVelocity(i)+ r2 * PARTICLE_PARAMS::c2_social() * (this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i))); */
      /* 		} */
      /* 	      break; */
      /* 	    case BAREBONE: */
      /* 	      for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 		{ */
      /* 		  mean = 0.5*(this->getBestPosition()->getPosition(i) + this->getNeighborhood()->getBest()->getPosition(i)); */
      /* 		  var = fabs(this->getBestPosition()->getPosition(i) - this->getNeighborhood()->getBest()->getPosition(i)); */
      /* 		  this->setVelocity(i, popot::math::normal(mean,var)); */
      /* 		} */
      /* 	      break; */
      /* 	    case MODIFIED_BAREBONE: */
      /* 	      if(popot::math::uniform_random(0.0,1.0) < 0.5) */
      /* 		{ */
      /* 		  for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 		    this->setVelocity(i, this->getBestPosition()->getPosition(i)); */
      /* 		} */
      /* 	      else */
      /* 		{ */
      /* 		  for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 		    { */
      /* 		      mean = 0.5*(this->getBestPosition()->getPosition(i) + this->getNeighborhood()->getBest()->getPosition(i)); */
      /* 		      var = fabs(this->getBestPosition()->getPosition(i) - this->getNeighborhood()->getBest()->getPosition(i)); */
      /* 		      this->setVelocity(i, popot::math::normal(mean,var)); */
      /* 		    } */
      /* 		} */
      /* 	      break; */
      /* 	    } */
      /* 	} */
      /* }; // Heterogenous PSO */

    } // namespace particle
  } // namespace PSO

  // Artificial Bee Colony
  namespace ABC
  {
    namespace individuals
    {
      

      template< typename PROBLEM>
	class FoodSource : public Vector<double, PROBLEM::nb_parameters>
	{
	  typedef Vector<double, PROBLEM::nb_parameters> TVector;

	private:
	  double fitness;
	  double fvalue;
	  int counter;

	public:
	FoodSource(void): TVector()
	    {}

	  /**
	   * Copy constructor
	   */
	FoodSource(const FoodSource & other) : TVector(other), fitness(other.fitness), fvalue(other.fvalue), counter(other.counter)
	    {
	    }

	  virtual ~FoodSource(void)
	    {}

	  void init(void)
	  {
	    for(int i = 0 ; i < PROBLEM::nb_parameters ; ++i)
	      setValueAt(i, popot::math::uniform_random(PROBLEM::get_lbound(i), PROBLEM::get_ubound(i)));

	    computeFitness();

	    // Initialize the counter
	    counter = 0;
	  }

	  void combine(const FoodSource & other_source)
	  {
	    // Combination operator
	    int change_dim;
	    double phi;
	    FoodSource new_source;
	    double new_param_value;

	    // Randomly select a dimension to change
	    change_dim = (int) popot::math::uniform_random(0, PROBLEM::nb_parameters);
	  
	    // Random combination coefficient
	    phi = popot::math::uniform_random(-1.0, 1.0);

	    // Compute the new source
	    new_source = *this;
	    new_param_value = this->getValueAt(change_dim) + phi * (this->getValueAt(change_dim) - other_source.getValueAt(change_dim));

	    // Bound the parameter value
	    if(new_param_value < PROBLEM::get_lbound(change_dim))
	      new_param_value = PROBLEM::get_lbound(change_dim);
	    else if(new_param_value > PROBLEM::get_ubound(change_dim))
	      new_param_value = PROBLEM::get_ubound(change_dim);
	  
	    new_source.setValueAt(change_dim, new_param_value);

	    // Evaluate the fitness of the new solution
	    new_source.computeFitness();

	    // Perform a greedy selection
	    if(new_source.getFitness() > getFitness())
	      {
		// If we get a better solution, it takes the place of the previous one
		*this = new_source;
		resetCounter();
	      }
	    else
	      {
		// otherwise we increment its counter
		(*this)++;
	      }
	  }

	  void resetCounter(void)
	  {
	    counter = 0;
	  }


	  int getCounter(void)
	  {
	    return counter;
	  }

	  void operator++(int) { 
	    counter ++;
	  }

	  void computeFitness(void)
	  {
	    // Compute the fitness
	    fvalue = PROBLEM::evaluate(this->getValuesPtr());
	    fitness = fitnessFunction(fvalue);
	  }

	  double fitnessFunction(double x)
	  {
	    if(x >= 0.0)
	      return 1.0 / (1.0 + x);
	    else
	      return 1.0 + fabs(x);
	  }

	  double getFitness(void)
	  {
	    return fitness;
	  }

	  double getFValue(void)
	  {
	    return fvalue;
	  }

	  virtual void print(std::ostream & os) const
	  {
	    TVector::print(os);
	    os << " Fitness : " << fitness << " Fvalue : " << fvalue << " Count : " << counter;
	  }

	  /**
	   * Serialization operator, for display
	   */
	  friend std::ostream & operator <<(std::ostream & os, const FoodSource &b)
	  {
	    b.print(os);
	    return os;
	  }

	}; // FoodSource     
    } // namespace individuals

  } // namespace ABC

}//namespace popot

#endif // POPOT_INDIVIDUALS_H
