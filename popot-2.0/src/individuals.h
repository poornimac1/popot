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
       * Constructor
       */
    Vector(size_t dimension) 
      : _dimension(dimension)
      {
	_values = new VALUE_TYPE[_dimension];
      }

    Vector(void) : _dimension(0), _values(0)
	{}

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
      void setValueAt(size_t index, VALUE_TYPE value)
      {
	if(index >= 0 && index < _dimension)
	  _values[index] = value;
	else
	  throw popot::Exception::IndexOutOfRange(index, _dimension);
      }

      /**
       * Returns a component of the vector
       * @param index position
       */
      VALUE_TYPE getValueAt(size_t index) const
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
       * Method for displaying a vector, fills in the stream
       */
      virtual void print(std::ostream & os) const
      {
	os << "[";
	for(size_t i = 0 ; i < size() ; ++i)
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

      template<typename POSITION_INITIALIZER, typename LBOUND_FUNC, typename UBOUND_FUNC, typename COST_FUNC>
      class Base: public Vector<double>
      {

      private:
	typedef Vector<double> TSuper;

      protected:
	const POSITION_INITIALIZER& _pinit;
	const LBOUND_FUNC& _lbound;
	const UBOUND_FUNC& _ubound;
	const COST_FUNC& _cost_func;

	double _fitness;

      public:

	Base() 
	  : TSuper(),
	    _pinit(0),
	    _lbound(0),
	    _ubound(0),
	    _cost_func(0)
	{}

	Base(const POSITION_INITIALIZER & pinit, 
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound, 
	     const COST_FUNC& cost_func)
	  : TSuper() ,
	    _pinit(pinit),
	    _lbound(lbound),
	    _ubound(ubound),
	    _cost_func(cost_func),
	    _fitness(0)
	{}
	
	Base(size_t dimension, const POSITION_INITIALIZER & pinit, 
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const COST_FUNC& cost_func)
	  : TSuper(dimension),
	    _pinit(pinit),
	    _lbound(lbound),
	    _ubound(ubound),
	    _cost_func(cost_func),
	    _fitness(0)
	{}
	
	virtual ~Base(void)
	{}

      	Base(const Base & other) 
	  : TSuper(other),
	    _pinit(other._pinit),
	    _lbound(other._lbound),
	    _ubound(other._ubound),
	    _cost_func(other._cost_func),
	    _fitness(other._fitness)
	{}

	Base & operator=(const Base &other)
	{
	  this->TSuper::operator=(other);
	  _pinit = other._pinit;
	  _lbound = other._lbound;
	  _ubound = other._bound;
	  _cost_func = other._cost_func;
	  _fitness = other._fitness;
	  return *this;
	}

	/**
	 * Random intialization of the position followed by a fitness evaluation
	 */
	virtual void init()
	{
	  // Initialize the position
	  initPosition();

	  // Evaluate the fitness
	  evaluateFitness();
	}

	/**
	 * Initialization of the position
	 */
	void initPosition()
	{
	  _pinit(this->_dimension, this->_lbound, this->_ubound, this->getValuesPtr());
	}

	/**
	 * Getter on the position
	 */
	double getPosition(int i) const
	{
	  return this->getValueAt(i);
	}
	/**
	 * Setter on the position
	 */
	void setPosition(int i, double val)
	{
	  this->setValueAt(i, val);
	}

	/**
	* Ensures the position is within the boundaries
	*/
	virtual void confine()
	{
	  for(size_t i = 0 ; i < this->_dimension ; ++i)
	    {
	      if(getPosition(i) < this->_lbound(i))
		setPosition(i, this->_lbound(i));
	      else if(getPosition(i) > this->_ubound(i))
		setPosition(i, this->_ubound(i));
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
	      virtual double evaluateFitness()
	{
	  _fitness = _cost_func(this->getValuesPtr());
	  return _fitness;
	}

	/**
	 * Comparison of the fitness of two particles
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

	virtual void print(std::ostream & os) const
	{
	  TSuper::print(os);
	  os << " ; Fitness : " << this->getFitness();
	}

      };

      template<typename POSITION_INITIALIZER, typename LBOUND_FUNC, typename UBOUND_FUNC, typename COST_FUNC>
      Base<POSITION_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC> base(size_t dimension,
								const POSITION_INITIALIZER& pinit,
								const LBOUND_FUNC& lbound, 
								const UBOUND_FUNC& ubound,
								const COST_FUNC& cost_func) 
      {
	return Base<POSITION_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC>(dimension, pinit, lbound, ubound, cost_func);
      };

      // Particle introduces the notion
      // of neighborhood and velocity
      // and best position
      template< typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER, 
		typename LBOUND_FUNC, typename UBOUND_FUNC,
		typename COST_FUNC>
      class Particle : public Base<POSITION_INITIALIZER,
				   LBOUND_FUNC, UBOUND_FUNC,
				   COST_FUNC>
      	{
      	private:
      	  typedef Base<POSITION_INITIALIZER,
		       LBOUND_FUNC, UBOUND_FUNC,
		       COST_FUNC> TSuper;
      	  typedef Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER,
			   LBOUND_FUNC, UBOUND_FUNC,
			   COST_FUNC> ThisParticleType;

	public:
	  typedef TSuper BestType;
      	  typedef popot::PSO::neighborhood::Neighborhood< ThisParticleType > NeighborhoodType;

      	protected:
      	  double *_velocity;
      	  BestType _best_position;
      	  NeighborhoodType *_neighborhood;
	  const VELOCITY_INITIALIZER& _vinit;
	  
      	public:

	  Particle(const POSITION_INITIALIZER& pinit, const VELOCITY_INITIALIZER &vinit,
		   const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound, 
	     const COST_FUNC& cost_func) 
	    : TSuper(pinit, lbound, ubound, cost_func),
	      _vinit(vinit)
      	    {}

	  Particle(size_t dimension,
		   const POSITION_INITIALIZER& pinit, const VELOCITY_INITIALIZER &vinit,
		   const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound, 
		   const COST_FUNC& cost_func) 
	    : TSuper(dimension, pinit, lbound, ubound, cost_func), _vinit(vinit)
      	    {
      	      _velocity = new double[TSuper::dimension];
      	      memset(_velocity, 0, TSuper::dimension * sizeof(double));

      	      // Create an empty neighborhood
      	      _neighborhood = new NeighborhoodType();
      	    }

      	Particle(const Particle& other) 
	  : TSuper(other), _vinit(other._vinit)
      	    {
      	      // Copy the velocity
      	      _velocity = new double[this->_dimension];
      	      for(size_t i = 0 ; i < this->_dimension ; ++i)
      		_velocity[i] = other.getVelocity(i);

      	      _neighborhood = new NeighborhoodType(*(other.getNeighborhood()));
      	    }

      	  virtual ~Particle(void)
      	    {
      	      delete[] _velocity;
      	      delete _neighborhood;
      	    }

      	  /**
      	   * Initialization of the position, velocity and fitness of the particle and set the best_position as the current position
      	   */
      	  virtual void init()
      	  {
      	    // Initialize the position
      	    // and computes the fitness
      	    TSuper::init();

      	    // Initialize the velocity
      	    initVelocity();

      	    // Set the best particle to the current position
      	    // Just copies the position and fitness
      	    _best_position = *this;
      	  }


      	  /**
      	   * Initialization of the velocity
      	   */
      	  virtual void initVelocity()
      	  {
	    _vinit(this->_dimension, this->_lbound, this->_ubount, this->getValuesPtr(), _velocity);
	    /*
      	    double params[3];
      	    for(size_t i = 0 ; i < TSuper::dimension ; ++i)
      	      {
      		params[0] = p.get_lbound(i);
      		params[1] = p.get_ubound(i);
      		params[2] = this->getPosition(i);
      		setVelocity(i, VELOCITY_INITIALIZER::init(params));
      	      }
	    */
      	  }

      	  /**
      	   * Getter on the velocity
      	   */
      	  double getVelocity(size_t i) const
      	  {
      	    if(i >= 0 && i < TSuper::dimension)
      	      return _velocity[i];
      	    else
      	      throw popot::Exception::IndexOutOfRange(i, TSuper::dimension);
      	  }

      	  /**
      	   * Setter on the velocity
      	   */
      	  void setVelocity(size_t i, double value)
      	  {
      	    if(i >= 0 && i < this->_dimension)
      	      _velocity[i] = value;
      	    else
      	      throw popot::Exception::IndexOutOfRange(i, this->_dimension);
      	  }

      	  /**
      	   * Bounds the velocity and position of the particle
      	   * If the position is out of the boundaries, set the position on the
      	   * boundaries and the velocity to zero
      	   */
      	  virtual void confine()
      	  {
      	    // In case the position is out of the bounds
      	    // we reset the velocities
      	    for(size_t i = 0 ; i < this->_dimension ; ++i)
      	      {
      		if((this->getPosition(i) < this->_lbound(i)))
      		  {
      		    this->setPosition(i, this->_lbound(i));
      		    this->setVelocity(i,0);
      		  }
      		else if(this->getPosition(i) > this->_ubound(i))
      		  {
      		    this->setPosition(i, this->_ubound(i));
      		    this->setVelocity(i,0);
      		  }
      	      }
      	  }

      	  /**
      	   * Updates the best position
      	   */
      	  virtual void updateBestPosition(void)
      	  {
      	    if(VERBOSE_BENCH)
      	      {
      		std::cout << "Previous best : " << _best_position.getFitness() << " : ";
      		for(size_t i = 0 ; i < this->_dimension ; ++i)
      		  std::cout << _best_position.getPosition(i) << " ";
      		std::cout << std::endl;
      	      }
      	    // Update the best position the particle ever had
      	    if(this->compare(_best_position) < 0)
      	      _best_position = *this;
	    
      	    if(VERBOSE_BENCH)
      	      {
      		std::cout << "New best : " << _best_position.getFitness() << " : ";
      		for(size_t i = 0 ; i < this->_dimension ; ++i)
      		  std::cout << _best_position.getPosition(i) << " ";
      		std::cout << std::endl;
      	      }
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

      	  virtual BestType & getBestPosition(void) const
      	  {
      	    return _best_position;
      	  }

      	  /**
      	   * Returns a pointer to the neighborhood of the particle
      	   */
      	  NeighborhoodType* getNeighborhood(void) const
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
      	    /*
      	    printf("%f - ", this->getFitness());
      	    for(int i = 0 ; i < TSuper::_dimension ; ++i)
      	      {
      		printf("(%f,%f,%f) ", this->getPosition(i),this->getVelocity(i),_best_position.getPosition(i));
      	      }
      	    printf("\n");
      	    */
      	    TSuper::print(os);
      	    os << "; Best position : " << _best_position;
      	  }

      	};


      template<typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER,
	       typename LBOUND_FUNC, typename UBOUND_FUNC, 
	       typename COST_FUNC>
      Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC> 
      particle(size_t dimension,
	       const POSITION_INITIALIZER& pinit,
	       const VELOCITY_INITIALIZER& vinit,
	       const LBOUND_FUNC& lbound, 
	       const UBOUND_FUNC& ubound,
	       const COST_FUNC& cost_func) 
      {
	return Particle<POSITION_INITIALIZER, VELOCITY_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC, COST_FUNC>(dimension, pinit, vinit, lbound, ubound, cost_func);
      };


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
	    this->setValueAt(i, popot::math::uniform_random(lbound(i), ubound(i)));

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
	  new_param_value = this->getValueAt(change_dim) + phi * (this->getValueAt(change_dim) - other_source.getValueAt(change_dim));

	  // Bound the parameter value
	  if(new_param_value < lbound(change_dim))
	    new_param_value = lbound(change_dim);
	  else if(new_param_value > ubound(change_dim))
	    new_param_value = ubound(change_dim);
	  
	  new_source.setValueAt(change_dim, new_param_value);

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
