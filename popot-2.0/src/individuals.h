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
   * @brief The template parameters is the type of values (e.g. double, or bool)
   */
  template< typename VALUE_TYPE>
  class Vector
  {
  protected:
    size_t _dimension;
    VALUE_TYPE * _values;

  public:

    /**
     * Default constructor, dimension=0 is assumed
     */
    Vector(void) : _dimension(0), _values(0)
    {}
    
    /**
     * Constructor
     */
    Vector(size_t dimension) 
      : _dimension(dimension)
    {
      _values = new VALUE_TYPE[_dimension];
      for(size_t i = 0 ; i < dimension; ++i)
	_values[i] = 0;
    }

    /**
     * Copy constructor
     */
    Vector(const Vector& other)
    {
      _dimension = other._dimension;
      _values = new VALUE_TYPE[_dimension];
      for(size_t i = 0 ; i < _dimension ; ++i)
	_values[i] = other._values[i];
    }

    /**
     * Copy Constructor
     */
    Vector & operator=(const Vector& other)
    {
      if(_dimension != other._dimension){
	delete[] _values;
	_dimension = other._dimension;
	_values = new VALUE_TYPE[_dimension];
      }
      for(size_t i = 0 ; i < _dimension ; ++i)
	_values[i] = other._values[i];

      return *this;
    }

    /**
     * Destructor
     */
    virtual ~Vector(void)
    {
      delete[] _values;
    }

    /**
     * Setter for a value
     * @param index position
     * @param value value to set
     */
    /*
    void setValueAt(size_t index, VALUE_TYPE value)
    {
      if(index >= 0 && index < _dimension)
	_values[index] = value;
      else
	throw popot::Exception::IndexOutOfRange(index, _dimension);
    }
    */
    /**
     * Returns a component of the vector
     * @param index position
     */
    /*
    VALUE_TYPE getValueAt(size_t index) const
    {
      if(index >= 0 && index < _dimension)
	return _values[index];
      else
	throw popot::Exception::IndexOutOfRange(index, _dimension);
    }
    */
    VALUE_TYPE& operator[](size_t index) const 
    {
      if(index >= 0 && index < _dimension)
	return _values[index];
      else
	throw popot::Exception::IndexOutOfRange(index, _dimension);
    }

    /**
     * Returns the size of the vector
     */
    size_t size(void) const
    {
      return _dimension;
    }

    /**
     * Returns a pointer to the raw data
     */
    VALUE_TYPE * getValuesPtr(void)
    {
      return _values;
    }

    /**
     * Prints the current position + fitness
     */
    virtual void print(std::ostream & os) const
    {
      os << "[";
      for(size_t i = 0 ; i < size() ; ++i)
	{
	  os << (*this)[i];
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

      template<typename TVECTOR_TYPE=Vector<double> >
      class Base
      {
      protected:
	TVECTOR_TYPE _position;
	double _fitness;

      public:
	typedef TVECTOR_TYPE VECTOR_TYPE;

	/**
	 * Default constructor, dimension=0 is assumed
	 */
	Base()
	  : _position(),
	    _fitness(0)
	{}
	
	/**
	 * Default constructor
	 */
	Base(size_t dimension)
	  : _position(dimension),
	    _fitness(0)
	{}
	
	/**
	 * Copy constructor
	 */
	Base(const Base & other) 
	  : _position(other._position),
	    _fitness(other._fitness)
	{}

	/**
	 * Destructor
	 */
	virtual ~Base(void)
	{}

	/**
	 * Assignement operator
	 */
	Base & operator=(const Base &other)
	{
	  if (this == &other) 
	    return *this;

	  _position = other._position;
	  _fitness = other._fitness;
	  return *this;
	}

	/**
	 * Random intialization of the position followed by a fitness evaluation
	 */
      /*
      template<typename INIT_FUNCTION, typename COST_FUNCTION>
      void init(const INIT_FUNCTION& init_function, const COST_FUNCTION& cost_function)
	{
	  // Initialize the position
	  init_function(*this);

	  // Evaluate the fitness
	  evaluateFitness(cost_function);
	}
	*/
	/**
	 * Getter on the position
	 */
	TVECTOR_TYPE& getPosition()
	{
	  return _position;
	}

	/**
	 * Ensures the position is within the boundaries
	 */
	/*
	template<typename CONFINE_FUNCTION>
	void confine(const CONFINE_FUNCTION& confine_function)
	{
	  confine_function(*this);
	}
	*/
	

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
	
	template< typename COST_FUNCTION>
	double evaluateFitness(const COST_FUNCTION& cost_function)
	{
	  _fitness = cost_function(this->getPosition());
	  return _fitness;
	}
	

	/**
	 * Comparison of two particles through the fitness
	 */
	bool operator<(const Base &p) const
	{
	  return (compare(p) < 0);
	}

	/**
	 * Comparison of the fitness of two particles p1.compare(p2)
	 * @return -1 if p1.f < p2.f
	 * @return 1 if p1.f > p2.f
	 * @return 0 otherwise
	 */
	virtual int compare(const Base& p) const
	{
	  double myfitness = getFitness();
	  double otherfitness = p.getFitness();
	  if(myfitness < otherfitness)
	    return -1;
	  else if(myfitness > otherfitness)
	    return 1;
	  else
	    return 0;
	}

	/**
	 * Prints the current position + fitness
	 */
	virtual void print(std::ostream & os) const
	{
	  os << "Position : " << this->_position;
	  os << " ; Fitness : " << this->getFitness();
	}

	/**
	 * Serialization operator 
	 */
	friend std::ostream & operator <<(std::ostream & os, const Base &v)
	{
	  v.print(os);
	  return os;
	}
      };


      // Particle introduces the notion
      // of neighborhood, velocity and best position
      template<typename TVECTOR_TYPE=Vector<double> >
      class Particle : public Base<TVECTOR_TYPE>
      {
      private:
	typedef Base<TVECTOR_TYPE> TSuper;
	typedef Particle<TVECTOR_TYPE> ThisParticleType;

      public:
	typedef TSuper BestType;
	typedef popot::PSO::neighborhood::Neighborhood< ThisParticleType > NeighborhoodType;

	TVECTOR_TYPE _velocity;
	BestType _best_position;
	NeighborhoodType _neighborhood;
	  
      public:

	/**
	 * Default constructor, dimension=0 is assumed
	 */
	Particle() 
	  : TSuper(),
	    _best_position(),
	    _neighborhood()
	{}

	/**
	 * Default constructor
	 */
	Particle(size_t dimension) 
	  : TSuper(dimension), 
	    _velocity(dimension),
	    _best_position(dimension),
	    _neighborhood()
	{}

	/**
	 * Copy constructor
	 */
	Particle(const Particle& other) 
	  : TSuper(other), 
	    _velocity(other._velocity),
	    _best_position(other._best_position),
	    _neighborhood(other._neighborhood)
	{}

	/**
	 * Destructor
	 */
	virtual ~Particle(void)
	{}

	/**
	 * Getter on the velocity
	 */
	TVECTOR_TYPE& getVelocity(void)
	{
	  return _velocity;
	}

	/**
	 * Bounds the velocity and position of the particle
	 * If the position is out of the boundaries, set the position on the
	 * boundaries and the velocity to zero
	 */
	
	template<typename LBOUND_FUNC, typename UBOUND_FUNC>
	void confine(const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
	{
	  // In case the position is out of the bounds
	  // we reset the velocities
	  TVECTOR_TYPE& pos = this->getPosition();
	  for(size_t i = 0 ; i < pos.size() ; ++i)
	    {
	      if(pos[i] < lbound(i))
		{
		  pos[i] = lbound(i);
		  this->getVelocity()[i] = 0;
		}
	      else if(pos[i] > ubound(i))
		{
		  pos[i] = ubound(i);
		  this->getVelocity()[i] = 0;
		}
	    }
	}


      /**
       * Initialization of the best position to the current position
       */
      void initBestPosition()
      {
	_best_position = *this;
      }

	/**
	 * Returns a reference to the current best position
	 */
	BestType& getBestPosition(void)
	{
	  return this->_best_position;
	}

	/**
	 * Returns a pointer to the neighborhood of the particle
	 */
	NeighborhoodType& getNeighborhood(void)
	{
	  return _neighborhood;
	}
	  
	/**
	 * Display the current position+fitness
	 * and the current best position + fitness
	 */
	virtual void print(std::ostream & os) const
	{
	  TSuper::print(os);
	  os << "; Velocity : " << _velocity;
	  os << "; Best position : " << _best_position;
	}

      };

      /*
      template<typename PARTICLE, typename POSITION_INITIALIZER>
      void init_position(PARTICLE& p, const POSITION_INITIALIZER& pinit)
      {
	pinit(p.getPosition());
      }

      template<typename PARTICLE, typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER>
      void init_position_velocity(PARTICLE& p, const POSITION_INITIALIZER& pinit, const VELOCITY_INITIALIZER& vinit)
      {
	pinit(p.getPosition());
	vinit(p.getPosition(), p.getVelocity());
      }
      */


      template<typename PARTICLE>
      void updateBestPosition(PARTICLE& p)
      {
	// Update the best position the particle ever had
	// with a copy of the current position
	if(p.compare(p.getBestPosition()) < 0)
	  p._best_position = p;
      }

      template<typename PARTICLE>
      void updateBestPosition_reevaluateBest(PARTICLE& p)
      {
	// Reevalute the fitness of the personal best before
	// possibly changing the personal best
	p._best_position.evaluateFitness();

	updateBestPosition(p);
      }


      /**
       * Basic update of the position of a particle
       */
      template<typename PARTICLE>
      void updatePosition(PARTICLE& p)
      {
	// Here it is simply : p_{k+1} = p_k + v_k      
	for(size_t i = 0 ; i < p.getPosition().size() ; ++i)
	  p.getPosition()[i] = p.getPosition()[i] + p.getVelocity()[i];
	
      }

      /**
       * Updates the velocity of the particle
       */
      template<typename PARTICLE, typename PARAMS>
      void updateVelocity_spso2006(PARTICLE& p)
      {
	// The update of the velocity is done according to the equation :
	// v = w * v + c r1 (best_p - p) + c r2 (best_g - p)
	// with :
	// r_p, r_g two random real numbers in [0.0,1.0]
	// w, c1, c2 : user defined parameters
	// best_p : the best position the particle ever had
	// best_g : the best position the neighborhood ever had
	
	double r1,r2;
        typename PARTICLE::VECTOR_TYPE& pos = p.getPosition();
	typename PARTICLE::VECTOR_TYPE& vel = p.getVelocity();

	for(size_t i = 0 ; i < pos.size(); ++i)
	  {
	    r1 = popot::math::uniform_random(0.0, PARAMS::c());
	    r2 = popot::math::uniform_random(0.0, PARAMS::c());
	    vel[i] = PARAMS::w() * vel[i]
	      + r1 * (p.getBestPosition().getPosition()[i] - pos[i])
	      + r2 * (p.getNeighborhood().getBest()->getPosition()[i] - pos[i]);
	  }
      }


 
      /**
       * Definition of a Stochastic Standard PSO 2006 particle
       * @brief The template parameters is the problem you want to solve,
       *  and the parameters (inertia, accelaration) of the velocity update rule
       * Compared to the base type BaseParticle, this type adds the some specific position and velocity update rules
       * These come from the Standard PSO 2006 . The stochastic code simply reevaluates the fitness of the personal best
       * position before possibly changing it
       */
      
      // template<typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER,
      // 	       typename LBOUND_FUNC, typename UBOUND_FUNC, typename COST_FUNC>
      // class StochasticSPSO2006Particle : public SPSO2006Particle< POSITION_INITIALIZER, VELOCITY_INITIALIZER,
      // 								  LBOUND_FUNC, UBOUND_FUNC, COST_FUNC>
 
      // 	{
      // 	  typedef SPSO2006Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC> TSuper;

      // 	public:
      // 	  StochasticSPSO2006Particle(size_t dimension,
      // 				     const LBOUND_FUNC& lbound, 
      // 				     const UBOUND_FUNC& ubound,
      // 				     const COST_FUNC& cost_func,
      // 				     double w, double c,
      // 				     const POSITION_INITIALIZER& pinit=popot::initializer::position::uniform_random,
      // 				     const VELOCITY_INITIALIZER& vinit=popot::initializer::velocity::half_diff) 
      // 	    : TSuper(dimension, lbound, ubound, cost_func, w, c, pinit, vinit)
      // 	    {
      // 	    }

      // 	  /**
      // 	   * Updates the best position
      // 	   */
      // 	  virtual void updateBestPosition(void)
      // 	  {
      // 	    // Reevalute the fitness of the personal best before
      // 	    // possibly changing the personal best
      // 	    this->_best_position.evaluateFitness();

      // 	    // Update the best position the particle ever had
      // 	    // with a copy of the current position
      // 	    if(this->compare(this->getBestPosition()) < 0)
      // 	      this->_best_position = *this;
      // 	  }
      // 	};

      /**
       * Builds up a particle, its types is deduced by the compiler
       */
    /*
      template<typename LBOUND_FUNC, typename UBOUND_FUNC, 
	       typename COST_FUNC>
      StochasticSPSO2006Particle<void(*)(size_t, double(*)(size_t), double(*)(size_t), double *),
				 void(*)(size_t, double(*)(size_t), double(*)(size_t), double *, double *),
				 LBOUND_FUNC, UBOUND_FUNC, COST_FUNC> 
      make_stochastic_spso2006_particle(size_t dimension,
					const LBOUND_FUNC& lbound, 
					const UBOUND_FUNC& ubound,
					const COST_FUNC& cost_func, double w, double c) 
				 {
				   return StochasticSPSO2006Particle<void(*)(size_t, double(*)(size_t), double(*)(size_t), double *),
								     void(*)(size_t, double(*)(size_t), double(*)(size_t), double *, double *), 
								     LBOUND_FUNC, UBOUND_FUNC, COST_FUNC>(dimension, lbound, ubound, cost_func, w, c);
				 };

    */


//   template<typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER,
// 	   typename LBOUND_FUNC, typename UBOUND_FUNC, typename COST_FUNC>
//   class SPSO2007Particle : public Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC >
//   {

//   private:
//     typedef Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC > TSuper;
//     double _w, _c;

//   public:

//     SPSO2007Particle(size_t dimension,
// 		     const LBOUND_FUNC& lbound, 
// 		     const UBOUND_FUNC& ubound,
// 		     const COST_FUNC& cost_func,
// 		     double w, double c,
// 		     const POSITION_INITIALIZER& pinit=popot::initializer::position::uniform_random,
// 		     const VELOCITY_INITIALIZER& vinit=popot::initializer::velocity::half_diff) 
//       : TSuper(dimension, pinit, vinit, lbound, ubound, cost_func),
// 	_w(w), _c(c)
//     {
//     }

//     /**
//      * Destructor
//      */
//     virtual ~SPSO2007Particle(void)
//     {
//     }

//     /**
//      * Update the position of the particle
//      */
//     virtual void updatePosition(void)
//     {
//       // Here it is simply : p_{k+1} = p_k + v_k
//       for(size_t i = 0 ; i < this->_dimension ; ++i)
// 	this->setPosition(i, this->getPosition(i) + this->getVelocity(i));
//     }

//     /**
//      * Updates the velocity of the particle
//      */
//     virtual void updateVelocity(void)
//     {
//       // The update of the velocity is done according to the equation :
//       // v = w * v + r1 (best_p - p) + r2 (best_g - p)
//       // with :
//       // r_1, r_2 two random real numbers in [0.0,c]
//       // w, c : user defined parameters
//       // best_p : the best position the particle ever had
//       // best_g : the best position the neighborhood ever had

//       bool li_equals_pi = (&(this->getBestPosition()) == this->getNeighborhood()->getBest());
//       double r1,r2;
//       if(li_equals_pi)
// 	{
// 	  for(size_t i = 0 ; i < this->_dimension ; ++i)
// 	    {
// 	      r1 = popot::math::uniform_random(0.0, _c);
// 	      this->setVelocity(i, _w * this->getVelocity(i)
// 				+ r1 * (this->getBestPosition().getPosition(i) - this->getPosition(i)));
// 	    }
// 	}
//       else
// 	{
// 	  for(size_t i = 0 ; i < this->_dimension ; ++i)
// 	    {
// 	      r1 = popot::math::uniform_random(0.0, _c);
// 	      r2 = popot::math::uniform_random(0.0, _c);
// 	      this->setVelocity(i, _w * this->getVelocity(i)
// 				+ r1 * (this->getBestPosition().getPosition(i) - this->getPosition(i))
// 				+ r2 * (this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i)));
// 	    }
// 	}
//     }
//   };

//   /**
//    * Builds up a particle, its types is deduced by the compiler
//    */
//   template<typename LBOUND_FUNC, typename UBOUND_FUNC, 
// 	   typename COST_FUNC>
//   SPSO2007Particle<void(*)(size_t, double(*)(size_t), double(*)(size_t), double *),
// 		   void(*)(size_t, double(*)(size_t), double(*)(size_t), double *, double *),
// 		   LBOUND_FUNC, UBOUND_FUNC, COST_FUNC> 
//   make_spso2007_particle(size_t dimension,
// 			 const LBOUND_FUNC& lbound, 
// 			 const UBOUND_FUNC& ubound,
// 			 const COST_FUNC& cost_func, double w, double c) 
// 		   {
// 		     return SPSO2007Particle<void(*)(size_t, double(*)(size_t), double(*)(size_t), double *),
// 					     void(*)(size_t, double(*)(size_t), double(*)(size_t), double *, double *), 
// 					     LBOUND_FUNC, UBOUND_FUNC, COST_FUNC>(dimension, lbound, ubound, cost_func, w, c);
// };


//       /**
//        * Definition of a Standard PSO 2011 particle
//        * @brief The template parameters is the problem you want to solve,
//        *  and the parameters (inertia, accelaration) of the velocity update rule
//        * Compared to the base type BaseParticle, this type adds the some specific position and velocity update rules
//        * These come from the Standard PSO 2011
//        */
//       template<typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER,
// 	       typename LBOUND_FUNC, typename UBOUND_FUNC, 
// 	       typename COST_FUNC>
//       class SPSO2011Particle : public Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC >
 
//       {
//       private:

// 	typedef Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC > TSuper;
// 	double _w, _c;
// 	double * xpi;
// 	double * p1;
// 	double * p2;
// 	double * gr;
//       public:

//       	SPSO2011Particle(size_t dimension,
// 			 const LBOUND_FUNC& lbound, 
// 			 const UBOUND_FUNC& ubound,
// 			 const COST_FUNC& cost_func,
// 			 double w, double c,
// 			 const POSITION_INITIALIZER& pinit=popot::initializer::position::uniform_random,
// 			 const VELOCITY_INITIALIZER& vinit=popot::initializer::velocity::spso2011) 
// 	  : TSuper(dimension, pinit, vinit, lbound, ubound, cost_func),
// 	    _w(w), _c(c)
// 	{
// 	  xpi = new double[this->_dimension];
// 	  p1 = new double[this->_dimension];
// 	  p2 = new double[this->_dimension];
// 	  gr = new double[this->_dimension];
// 	}

// 	/**
// 	 * Destructor
// 	 */
// 	virtual ~SPSO2011Particle(void)
// 	{
// 	  delete[] xpi;
// 	  delete[] p1;
// 	  delete[] p2;
// 	  delete[] gr;
// 	}

// 	/**
// 	 * Update the position of the particle
// 	 */
// 	virtual void updatePosition(void)
// 	{
// 	  // Here it is simply : p_{k+1} = p_k + v_k
// 	  for(size_t i = 0 ; i < this->_dimension ; ++i)
// 	    this->setPosition(i, this->getPosition(i) + this->getVelocity(i));
// 	}

// 	/**
// 	 * Updates the velocity of the particle
// 	 */
// 	virtual void updateVelocity(void)
// 	{
// 	  // The update of the velocity is done according to the equation :
// 	  // v_i(t+1) = w v_i(t) + x'_i(t) - x_i(t)
// 	  // where x'_i is normally sampled from the hypersphere (Gi, Ri)
// 	  // with Gi the center and Ri the radius and
// 	  // Gi = 1/3 ( x_i + (x_i + c(p_i - x_i)) + (x_i + c(l_i - x_i)))
// 	  //    = 1/3 ( x_i +          p1          +            p2     )))
// 	  // or
// 	  // Gi = 1/2 (x_i + (x_i + c(p_i - x_i))) if p_i == l_i
// 	  //    = 1/2 (x_i +          p1         )
// 	  // Ri = ||G_i - x_i||

// 	  size_t i;

// 	  // We first check if the local best and personal best are identical by computing the norm of the difference
// 	  // between both positions
// 	  // This can be done with pointer comparison
// 	  // as the neighborhood holds a pointer to the best personal best
// 	  bool li_equals_pi = (&(this->getBestPosition()) == this->getNeighborhood()->getBest());

// 	  // First position
// 	  // p1 = xi + c * (pi - xi)
// 	  // where pi is the personal best position
// 	  for(i = 0 ; i < this->_dimension ; ++i)
// 	    p1[i] = this->getPosition(i) + _c *(this->getBestPosition().getPosition(i) - this->getPosition(i));
	  
// 	  // Second position
// 	  // p2 = xi + c * (li - xi)
// 	  // where li is the local best position (within the neighborhood)
// 	  for(i = 0 ; i < this->_dimension ; ++i)
// 	    p2[i] = this->getPosition(i) + _c *(this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i));
	  
// 	  // Compute the gravity center of p1, p2 and xi
// 	  if(!li_equals_pi)
// 	    {
// 	      // We here consider p1, p2 and xi
// 	      for(i = 0 ; i < this->_dimension ; ++i)
// 		gr[i] = 1.0/3.0 * (this->getPosition(i) + p1[i] + p2[i]);
// 	    }
// 	  else
// 	    {
// 	      // We here consider only p1 (or p2, since they are equal) and xi
// 	      for(i = 0 ; i < this->_dimension ; ++i)
// 		gr[i] = 1.0/2.0 * (this->getPosition(i) + p1[i]);
// 	    }

// 	  // And the radius of the hypersphere
// 	  double ri = 0.0;
// 	  for(i = 0 ; i < this->_dimension ; ++i)
// 	    ri += (gr[i] - this->getPosition(i))*(gr[i] - this->getPosition(i));
// 	  ri = sqrt(ri);

// 	  // Compute the auxiliary position x'_i randomly within the hypersphere (gr, ri);
// 	  // To uniformely sample from the hypersphere we uniformely sample a direction
// 	  double norm = 0.0;
// 	  for(i = 0 ; i < this->_dimension ; ++i)
// 	    {
// 	      xpi[i] = popot::math::normal(0.0,1.0);
// 	      norm += xpi[i] * xpi[i];
// 	    }
// 	  norm = sqrt(norm);
	  
// 	  // And then scale by a random radius
// 	  double r = popot::math::uniform_random(0.0,1.0);
// 	  for(i = 0 ; i < this->_dimension ; ++i)
// 	    xpi[i] =  gr[i] + r * ri * xpi[i] / norm;

// 	  // And then update the velocity
// 	  for(i = 0 ; i < this->_dimension ; ++i)
// 	    this->setVelocity(i, _w * this->getVelocity(i)
// 			      + xpi[i] - this->getPosition(i));
// 	}


// 	/**
// 	 * Bounds the velocity and position of the particle
// 	 */
	  
// 	virtual void confine(void)
// 	{
// 	  // In case the position is out of the bounds
// 	  // we reset the velocities
// 	  for(size_t i = 0 ; i < this->_dimension ; ++i)
// 	    {
// 	      if((this->getPosition(i) < this->_lbound(i)))
// 		{
// 		  this->setPosition(i, this->_lbound(i));
// 		  this->setVelocity(i,-0.5*this->getVelocity(i));
// 		}
// 	      else if(this->getPosition(i) > this->_ubound(i))
// 		{
// 		  this->setPosition(i, this->_ubound(i));
// 		  this->setVelocity(i,-0.5*this->getVelocity(i));
// 		}
// 	    }
// 	}
//       };
      
//       /**
//        * Builds up a particle, its types is deduced by the compiler
//        */
//       template<typename LBOUND_FUNC, typename UBOUND_FUNC, 
// 	       typename COST_FUNC>
//       SPSO2011Particle<void(*)(size_t, double(*)(size_t), double(*)(size_t), double *),
// 		       void(*)(size_t, double(*)(size_t), double(*)(size_t), double *, double *),
// 		       LBOUND_FUNC, UBOUND_FUNC, COST_FUNC> 
//       make_spso2011_particle(size_t dimension,
// 			     const LBOUND_FUNC& lbound, 
// 			     const UBOUND_FUNC& ubound,
// 			     const COST_FUNC& cost_func, double w, double c) 
// 		       {
// 			 return SPSO2011Particle<void(*)(size_t, double(*)(size_t), double(*)(size_t), double *),
// 						 void(*)(size_t, double(*)(size_t), double(*)(size_t), double *, double *), 
// 						 LBOUND_FUNC, UBOUND_FUNC, COST_FUNC>(dimension, lbound, ubound, cost_func, w, c);
// 		       };
      
//       /**
//        * Definition of a Stochastic Standard PSO 2011 particle
//        * @brief The template parameters is the problem you want to solve,
//        *  and the parameters (inertia, accelaration) of the velocity update rule
//        * Compared to the base type BaseParticle, this type adds the some specific position and velocity update rules
//        * These come from the Standard PSO 2011 . The stochastic code simply reevaluates the fitness of the personal best
//        * position before possibly changing it
//        */
//       template<typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER,
// 	       typename LBOUND_FUNC, typename UBOUND_FUNC, 
// 	       typename COST_FUNC>
//       class StochasticSPSO2011Particle : public SPSO2011Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC>
 
//       	{
//       	  typedef SPSO2011Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC> TSuper;

//       	public:
//      	StochasticSPSO2011Particle(size_t dimension,
// 				   const LBOUND_FUNC& lbound, 
// 				   const UBOUND_FUNC& ubound,
// 				   const COST_FUNC& cost_func,
// 				   double w, double c,
// 				   const POSITION_INITIALIZER& pinit=popot::initializer::position::uniform_random,
// 				   const VELOCITY_INITIALIZER& vinit=popot::initializer::velocity::spso2011)
// 	  : TSuper(dimension, lbound, ubound, cost_func, w, c, pinit, vinit)
//       	    {}

//       	  /**
//       	   * Updates the best position
//       	   */
//       	  virtual void updateBestPosition(void)
//       	  {
//       	    // Reevalute the fitness of the personal best before
//       	    // possibly changing the personal best
//       	    this->_best_position.evaluateFitness();

//       	    // Update the best position the particle ever had
//       	    // with a copy of the current position
//       	    if(this->compare(this->getBestPosition()) < 0)
//       	      this->_best_position = *this;
//       	  }
//       	};

//       /**
//        * Builds up a particle, its types is deduced by the compiler
//        */
//       template<typename LBOUND_FUNC, typename UBOUND_FUNC, 
// 	       typename COST_FUNC>
//       StochasticSPSO2011Particle<void(*)(size_t, double(*)(size_t), double(*)(size_t), double *),
// 		       void(*)(size_t, double(*)(size_t), double(*)(size_t), double *, double *),
// 		       LBOUND_FUNC, UBOUND_FUNC, COST_FUNC> 
//       make_stochastic_spso2011_particle(size_t dimension,
// 					const LBOUND_FUNC& lbound, 
// 					const UBOUND_FUNC& ubound,
// 					const COST_FUNC& cost_func, double w, double c) 
// 		       {
// 			 return StochasticSPSO2011Particle<void(*)(size_t, double(*)(size_t), double(*)(size_t), double *),
// 							   void(*)(size_t, double(*)(size_t), double(*)(size_t), double *, double *), 
// 							   LBOUND_FUNC, UBOUND_FUNC, COST_FUNC>(dimension, lbound, ubound, cost_func, w, c);
// 		       };

      
    } // namespace particle
  } // namespace PSO
  
  // Artificial Bee Colony
  namespace ABC
  {
    namespace individuals
    {
      class FoodSource : public Vector<double>
      {
	typedef Vector<double> TVector;

      private:
	double _fitness;
	double _fvalue;
	size_t _counter;

      public:

	FoodSource(size_t dimension)
	  : TVector(dimension),
	    _fitness(0), 
	    _fvalue(0), 
	    _counter(0)
	{}

	FoodSource(): TVector(), _fitness(0), _fvalue(0), _counter(0)
	{}

	/**
	 * Copy constructor
	 */
	FoodSource(const FoodSource & other) 
	  : TVector(other), 
	    _fitness(other._fitness), 
	    _fvalue(other._fvalue), 
	    _counter(other._counter)
	{
	}

	virtual ~FoodSource(void)
	{}

	template<typename LBOUND_FUNC, typename UBOUND_FUNC, typename COST_FUNCTION>
	void init(const LBOUND_FUNC& lbound, const UBOUND_FUNC &ubound, const COST_FUNCTION& cost_function)
	{
	  for(size_t i = 0 ; i < this->_dimension ; ++i)
	    (*this)[i] = popot::math::uniform_random(lbound(i), ubound(i));

	  computeFitness(cost_function);

	  // Initialize the counter
	  _counter = 0;
	}

	template<typename LBOUND_FUNC, typename UBOUND_FUNC, typename COST_FUNCTION>
	void combine(const FoodSource & other_source, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound, const COST_FUNCTION& cost_function)
	{
	  // Combination operator
	  int change_dim;
	  double phi;
	  FoodSource new_source;
	  double new_param_value;

	  // Randomly select a dimension to change
	  change_dim = (int) popot::math::uniform_random(0, this->_dimension);
	  
	  // Random combination coefficient
	  phi = popot::math::uniform_random(-1.0, 1.0);

	  // Compute the new source
	  new_source = *this;
	  new_param_value = (*this)[change_dim] + phi * ((*this)[change_dim] - other_source[change_dim]);

	  // Bound the parameter value
	  if(new_param_value < lbound(change_dim))
	    new_param_value = lbound(change_dim);
	  else if(new_param_value > ubound(change_dim))
	    new_param_value = ubound(change_dim);
	  
	  new_source[change_dim] = new_param_value;

	  // Evaluate the fitness of the new solution
	  new_source.computeFitness(cost_function);

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

	void resetCounter(void)  {
	  _counter = 0;
	}


	size_t getCounter(void) {
	  return _counter;
	}

	void operator++(int) { 
	  _counter ++;
	}

	template<typename COST_FUNCTION>
	void computeFitness(const COST_FUNCTION& cost_function)  {
	  // Compute the fitness
	  _fvalue = cost_function(this->getValuesPtr());
	  _fitness = fitnessFunction(_fvalue);
	}

	double fitnessFunction(double x)  {
	  if(x >= 0.0)
	    return 1.0 / (1.0 + x);
	  else
	    return 1.0 + fabs(x);
	}

	double getFitness(void)
	{
	  return _fitness;
	}

	double getFValue(void)
	{
	  return _fvalue;
	}

	virtual void print(std::ostream & os) const
	{
	  TVector::print(os);
	  os << " Fitness : " << _fitness << " Fvalue : " << _fvalue << " Count : " << _counter;
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
