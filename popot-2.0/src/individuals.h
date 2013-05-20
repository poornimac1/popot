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

      template<typename POSITION_INITIALIZER, typename LBOUND_FUNC, typename UBOUND_FUNC>
      class Base: public Vector<double>
      {
      private:
	typedef Vector<double> TSuper;

	const POSITION_INITIALIZER& _pinit;
	const LBOUND_FUNC& _lbound;
	const UBOUND_FUNC& _ubound;

      protected:
	double _fitness;

      public:

	Base(const POSITION_INITIALIZER & pinit, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
	  : TSuper() ,
	    _pinit(pinit),
	    _lbound(lbound),
	    _ubound(ubound),
	    _fitness(0)
	{}
	
	Base(size_t dimension, const POSITION_INITIALIZER & pinit, const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound)
	  : TSuper(dimension),
	    _pinit(pinit),
	    _lbound(lbound),
	    _ubound(ubound),
	    _fitness(0)
	{}
	
	virtual ~Base(void)
	{}

      	Base(const Base & other) 
	  : TSuper(other),
	    _pinit(other._pinit),
	    _lbound(other._lbound),
	    _ubound(other._ubound),
	    _fitness(other._fitness)
	{}

	Base & operator=(const Base &other)
	{
	  this->TSuper::operator=(other);
	  _pinit = other._pinit;
	  _lbound = other._lbound;
	  _ubound = other._bound;
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

	  /* // Evaluate the fitness */
	  /* evaluateFitness(p); */
	}

	/**
	 * Initialization of the position
	 */
	void initPosition()
	{
	  _pinit(this->_dimension, _lbound, _ubound, this->getValuesPtr());
	}
      };

      template<typename POSITION_INITIALIZER, typename LBOUND_FUNC, typename UBOUND_FUNC>
      Base<POSITION_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC> base(size_t dimension,
								const POSITION_INITIALIZER& pinit,
								const LBOUND_FUNC& lbound, 
								const UBOUND_FUNC& ubound) 
      {
	return Base<POSITION_INITIALIZER, LBOUND_FUNC, UBOUND_FUNC>(dimension, pinit, lbound, ubound);
      };

      /**
       * BaseParticle is a vector with a fitness and an associated Problem
       */
      /* template<typename PROBLEM, typename POSITION_INITIALIZER> */
      /* 	class BaseParticle : public Vector<double> */
      /* 	{ */
      /* 	private: */
      /* 	  typedef Vector<double> TSuper; */

      /* 	protected: */
      /* 	  double _fitness; */

      /* 	public: */

      /* 	BaseParticle(): TSuper() */
      /* 	    {} */

      /* 	BaseParticle(size_t dimension) : TSuper(dimension) */
      /* 	    {} */
	  
      /* 	  virtual ~BaseParticle(void) */
      /* 	    {} */

      /* 	BaseParticle(const BaseParticle & other) : TSuper(other) */
      /* 	    { */
      /* 	      _fitness = other._fitness; */
      /* 	    } */

      /* 	  BaseParticle & operator=(const BaseParticle &other) */
      /* 	    { */
      /* 	      this->TSuper::operator=(other); */
      /* 	      _fitness = other._fitness; */
      /* 	      return *this; */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Random intialization of the position followed by a fitness evaluation */
      /* 	   *\/	   */
      /* 	  virtual void init(PROBLEM& p) */
      /* 	  { */
      /* 	    // Initialize the position */
      /* 	    initPosition(p); */

      /* 	    // Evaluate the fitness */
      /* 	    evaluateFitness(p); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Initialization of the position */
      /* 	   *\/ */
      /* 	  void initPosition(PROBLEM& p) */
      /* 	  { */
      /* 	    POSITION_INITIALIZER::init(p,&PROBLEM::get_lbound, &PROBLEM::get_ubound, this->getValuesPtr()); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Getter on the position */
      /* 	   *\/ */
      /* 	  double getPosition(int i) const */
      /* 	  { */
      /* 	    return this->getValueAt(i); */
      /* 	  } */
      /* 	  /\** */
      /* 	   * Setter on the position */
      /* 	   *\/ */
      /* 	  void setPosition(int i, double val) */
      /* 	  { */
      /* 	    this->setValueAt(i, val); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Ensures the position is within the boundaries */
      /* 	   *\/ */
      /* 	  virtual void confine(PROBLEM &p) */
      /* 	  { */
      /* 	    for(size_t i = 0 ; i < p.dimension ; ++i) */
      /* 	      { */
      /* 		if(getPosition(i) < p.get_lbound(i)) */
      /* 		  setPosition(i, p.get_lbound(i)); */
      /* 		else if(getPosition(i) > p.get_ubound(i)) */
      /* 		  setPosition(i, p.get_ubound(i)); */
      /* 	      } */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Returns the currently known fitness */
      /* 	   *\/ */
      /* 	  double getFitness(void) const */
      /* 	  { */
      /* 	    return _fitness; */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Set the fitness */
      /* 	   *\/ */
      /* 	  void setFitness(double f) */
      /* 	  { */
      /* 	    _fitness = f; */
      /* 	  } */
	
      /* 	  /\** */
      /* 	   * Recompute the fitness */
      /* 	   *\/ */
      /* 	  virtual double evaluateFitness(PROBLEM& p) */
      /* 	  { */
      /* 	    _fitness = p.evaluate(this->getValuesPtr()); */
      /* 	    return _fitness; */
      /* 	  } */
      /* 	  /\** */
      /* 	   * Comparison of the fitness of two particles */
      /* 	   *\/ */
      /* 	  bool operator<(const BaseParticle &p) const */
      /* 	  { */
      /* 	    return (compare(p) < 0); */
      /* 	  } */
      /* 	  /\** */
      /* 	   * Comparison of the fitness of two particles p1.compare(p2) */
      /* 	   * @return -1 if p1.f < p2.f */
      /* 	   * @return 1 if p1.f > p2.f */
      /* 	   * @return 0 otherwise */
      /* 	   *\/ */
      /* 	  virtual int compare(const BaseParticle& p) const */
      /* 	  { */
      /* 	    double myfitness = getFitness(); */
      /* 	    double otherfitness = p.getFitness(); */
      /* 	    if(myfitness < otherfitness) */
      /* 	      return -1; */
      /* 	    else if(myfitness > otherfitness) */
      /* 	      return 1; */
      /* 	    else */
      /* 	      return 0; */
      /* 	  } */

      /* 	  virtual void print(std::ostream & os) const */
      /* 	  { */
      /* 	    TSuper::print(os); */
      /* 	    os << " ; Fitness : " << this->getFitness(); */
      /* 	  } */
	
      /* 	}; */


      // Particle introduces the notion
      // of neighborhood and velocity
      // and best position
      /* template< typename PROBLEM, typename POSITION_INITIALIZER, typename VELOCITY_INITIALIZER> */
      /* 	class Particle : public BaseParticle<PROBLEM, POSITION_INITIALIZER> */
      /* 	{ */
      /* 	private: */
      /* 	  typedef BaseParticle<PROBLEM, POSITION_INITIALIZER> TSuper; */
      /* 	  typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER > ThisParticleType; */

      /* 	public: */
      /* 	  typedef BaseParticle<PROBLEM, POSITION_INITIALIZER> BestType; */
      /* 	  typedef popot::PSO::neighborhood::Neighborhood< ThisParticleType > NeighborhoodType; */
	
	
      /* 	protected: */
      /* 	  double *_velocity; */
      /* 	  BestType _best_position; */
      /* 	  NeighborhoodType *_neighborhood; */
	
      /* 	public: */

      /* 	Particle() : TSuper() */
      /* 	    {} */

      /* 	Particle(size_t dimension) : TSuper(dimension) */
      /* 	    { */
      /* 	      _velocity = new double[TSuper::dimension]; */
      /* 	      memset(_velocity, 0, TSuper::dimension * sizeof(double)); */

      /* 	      // Create an empty neighborhood */
      /* 	      _neighborhood = new NeighborhoodType(); */
      /* 	    } */

      /* 	Particle(const Particle& other) : TSuper(other), _best_position(other.getBestPosition()) */
      /* 	    { */
      /* 	      // Copy the velocity */
      /* 	      _velocity = new double[TSuper::dimension]; */
      /* 	      for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		_velocity[i] = other.getVelocity(i); */

      /* 	      _neighborhood = new NeighborhoodType(*(other.getNeighborhood())); */
      /* 	      //_neighborhood = new NeighborhoodType(); */
      /* 	      //for(int i = 0 ; i < other.getNeighborhood()->size() ; ++i) */
      /* 	      // _neighborhood->add(other.getNeighborhood()->get(i)); */
      /* 	    } */

      /* 	  virtual ~Particle(void) */
      /* 	    { */
      /* 	      delete[] _velocity; */
      /* 	      delete _neighborhood; */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Initialization of the position, velocity and fitness of the particle and set the best_position as the current position */
      /* 	   *\/ */
      /* 	  virtual void init(PROBLEM &p) */
      /* 	  { */
      /* 	    // Initialize the position */
      /* 	    // and computes the fitness */
      /* 	    TSuper::init(p); */

      /* 	    // Initialize the velocity */
      /* 	    initVelocity(p); */

      /* 	    // Set the best particle to the current position */
      /* 	    // Just copies the position and fitness */
      /* 	    _best_position = *this; */
      /* 	  } */


      /* 	  /\** */
      /* 	   * Initialization of the velocity */
      /* 	   *\/ */
      /* 	  virtual void initVelocity(PROBLEM& p) */
      /* 	  { */
      /* 	    double params[3]; */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		params[0] = p.get_lbound(i); */
      /* 		params[1] = p.get_ubound(i); */
      /* 		params[2] = this->getPosition(i); */
      /* 		setVelocity(i, VELOCITY_INITIALIZER::init(params)); */
      /* 	      } */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Getter on the velocity */
      /* 	   *\/ */
      /* 	  double getVelocity(size_t i) const */
      /* 	  { */
      /* 	    if(i >= 0 && i < TSuper::dimension) */
      /* 	      return _velocity[i]; */
      /* 	    else */
      /* 	      throw popot::Exception::IndexOutOfRange(i, TSuper::dimension); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Setter on the velocity */
      /* 	   *\/ */
      /* 	  void setVelocity(size_t i, double value) */
      /* 	  { */
      /* 	    if(i >= 0 && i < TSuper::dimension) */
      /* 	      _velocity[i] = value; */
      /* 	    else */
      /* 	      throw popot::Exception::IndexOutOfRange(i, TSuper::dimension); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Bounds the velocity and position of the particle */
      /* 	   * If the position is out of the boundaries, set the position on the */
      /* 	   * boundaries and the velocity to zero */
      /* 	   *\/ */
      /* 	  virtual void confine(PROBLEM& p) */
      /* 	  { */
      /* 	    // In case the position is out of the bounds */
      /* 	    // we reset the velocities */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		if((this->getPosition(i) < p.get_lbound(i))) */
      /* 		  { */
      /* 		    this->setPosition(i, p.get_lbound(i)); */
      /* 		    this->setVelocity(i,0); */
      /* 		  } */
      /* 		else if(this->getPosition(i) > p.get_ubound(i)) */
      /* 		  { */
      /* 		    this->setPosition(i, p.get_ubound(i)); */
      /* 		    this->setVelocity(i,0); */
      /* 		  } */
      /* 	      } */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Updates the best position */
      /* 	   *\/ */
      /* 	  virtual void updateBestPosition(void) */
      /* 	  { */
      /* 	    if(VERBOSE_BENCH) */
      /* 	      { */
      /* 		std::cout << "Previous best : " << _best_position.getFitness() << " : "; */
      /* 		for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  std::cout << _best_position.getPosition(i) << " "; */
      /* 		std::cout << std::endl; */
      /* 	      } */
      /* 	    // Update the best position the particle ever had */
      /* 	    if(this->compare(_best_position) < 0) */
      /* 	      _best_position = *this; */
	    
      /* 	    if(VERBOSE_BENCH) */
      /* 	      { */
      /* 		std::cout << "New best : " << _best_position.getFitness() << " : "; */
      /* 		for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  std::cout << _best_position.getPosition(i) << " "; */
      /* 		std::cout << std::endl; */
      /* 	      } */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Updates the velocity of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    std::cout << "No updatePosition rule is provided for a particle of type Particle " << std::endl; */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Updates the velocity of the particle */
      /* 	   *\/ */
      /* 	  virtual void updateVelocity(void) */
      /* 	  { */
      /* 	    std::cout << "No updateVelocity rule is provided for a particle of type Particle " << std::endl; */
      /* 	  } */

      /* 	  virtual const BestType & getBestPosition(void) const */
      /* 	  { */
      /* 	    return _best_position; */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Returns a pointer to the neighborhood of the particle */
      /* 	   *\/ */
      /* 	  NeighborhoodType* getNeighborhood(void) const */
      /* 	  { */
      /* 	    return _neighborhood; */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Set the pointer to the neighborhood */
      /* 	   *\/ */
      /* 	  void setNeighborhood(NeighborhoodType * neighborhood) */
      /* 	  { */
      /* 	    _neighborhood = neighborhood; */
      /* 	  } */

      /* 	  virtual void print(std::ostream & os) const */
      /* 	  { */
      /* 	    /\* */
      /* 	    printf("%f - ", this->getFitness()); */
      /* 	    for(int i = 0 ; i < TSuper::_dimension ; ++i) */
      /* 	      { */
      /* 		printf("(%f,%f,%f) ", this->getPosition(i),this->getVelocity(i),_best_position.getPosition(i)); */
      /* 	      } */
      /* 	    printf("\n"); */
      /* 	    *\/ */
      /* 	    TSuper::print(os); */
      /* 	    os << "; Best position : " << _best_position; */
      /* 	  } */

      /* 	}; */


      /**
       * Definition of a Standard PSO 2006 particle
       * @brief The template parameters is the problem you want to solve and a class providing the parameters w and c
       * popot::PSO::initializer::PositionUniformRandom, popot::PSO::initializer::VelocityHalfDiff
       */
      /* template< typename PROBLEM, typename PARTICLE_PARAMS, typename POSITION_INITIALIZER=popot::PSO::initializer::PositionUniformRandom<PROBLEM>, typename VELOCITY_INITIALIZER=popot::PSO::initializer::VelocityHalfDiff> */
      /* 	class SPSO2006Particle : public Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER > */
      /* 	{ */
      /* 	  typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER > TSuper; */

      /* 	public: */

      /* 	SPSO2006Particle(size_t dimension) : TSuper(dimension) */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Destructor */
      /* 	   *\/ */
      /* 	  virtual ~SPSO2006Particle(void) */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Update the position of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = p_k + v_k */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i)); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Updates the velocity of the particle */
      /* 	   *\/ */
      /* 	  virtual void updateVelocity(void) */
      /* 	  { */
      /* 	    // The update of the velocity is done according to the equation : */
      /* 	    // v = w * v + c r1 (best_p - p) + c r2 (best_g - p) */
      /* 	    // with : */
      /* 	    // r_p, r_g two random real numbers in [0.0,1.0] */
      /* 	    // w, c1, c2 : user defined parameters */
      /* 	    // best_p : the best position the particle ever had */
      /* 	    // best_g : the best position the neighborhood ever had */
      /* 	    double r1,r2; */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		r1 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c()); */
      /* 		r2 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c()); */
      /* 		this->setVelocity(i, PARTICLE_PARAMS::w() * this->getVelocity(i) */
      /* 				  + r1 * (this->getBestPosition().getPosition(i) - this->getPosition(i)) */
      /* 				  + r2 * (this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i))); */
      /* 	      } */
      /* 	  } */
      /* 	}; */

      /**
       * Definition of a Standard PSO 2007 particle
       * @brief The template parameters is the problem you want to solve and a class providing the parameters w and c
       */
      /* template< typename PROBLEM, typename PARTICLE_PARAMS, typename POSITION_INITIALIZER=popot::PSO::initializer::PositionUniformRandom<PROBLEM>, typename VELOCITY_INITIALIZER=popot::PSO::initializer::VelocityHalfDiff > */
      /* 	class SPSO2007Particle : public Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER > */
      /* 	{ */
      /* 	  typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> TSuper; */

      /* 	public: */

      /* 	SPSO2007Particle(size_t dimension) : TSuper(dimension) */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Destructor */
      /* 	   *\/ */
      /* 	  virtual ~SPSO2007Particle(void) */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Update the position of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = p_k + v_k */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i)); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Updates the velocity of the particle */
      /* 	   *\/ */
      /* 	  virtual void updateVelocity(void) */
      /* 	  { */
      /* 	    // The update of the velocity is done according to the equation : */
      /* 	    // v = w * v + r1 (best_p - p) + r2 (best_g - p) */
      /* 	    // with : */
      /* 	    // r_1, r_2 two random real numbers in [0.0,c] */
      /* 	    // w, c : user defined parameters */
      /* 	    // best_p : the best position the particle ever had */
      /* 	    // best_g : the best position the neighborhood ever had */

      /* 	    bool li_equals_pi = (&(this->getBestPosition()) == this->getNeighborhood()->getBest()); */
      /* 	    double r1,r2; */
      /* 	    if(li_equals_pi) */
      /* 	      { */
      /* 		for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  { */
      /* 		    r1 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c()); */
      /* 		    this->setVelocity(i, PARTICLE_PARAMS::w() * this->getVelocity(i) */
      /* 				      + r1 * (this->getBestPosition().getPosition(i) - this->getPosition(i))); */
      /* 		  } */
      /* 	      } */
      /* 	    else */
      /* 	      { */
      /* 		for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  { */
      /* 		    r1 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c()); */
      /* 		    r2 = popot::math::uniform_random(0.0,PARTICLE_PARAMS::c()); */
      /* 		    this->setVelocity(i, PARTICLE_PARAMS::w() * this->getVelocity(i) */
      /* 				      + r1 * (this->getBestPosition().getPosition(i) - this->getPosition(i)) */
      /* 				      + r2 * (this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i))); */
      /* 		  } */
      /* 	      } */
      /* 	  } */
      /* 	}; */


      /**
       * Definition of a Standard PSO 2011 particle
       * @brief The template parameters is the problem you want to solve,
       *  and the parameters (inertia, accelaration) of the velocity update rule
       * Compared to the base type BaseParticle, this type adds the some specific position and velocity update rules
       * These come from the Standard PSO 2011
       */
      /* template< typename PROBLEM, typename PARTICLE_PARAMS, typename POSITION_INITIALIZER=popot::PSO::initializer::PositionUniformRandom<PROBLEM>, typename VELOCITY_INITIALIZER=popot::PSO::initializer::VelocitySPSO2011> */
      /* 	class SPSO2011Particle : public Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER > */
 
      /* 	{ */
      /* 	  typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER > TSuper; */

      /* 	  double * xpi; */
      /* 	  double * p1; */
      /* 	  double * p2; */
      /* 	  double * gr; */
      /* 	public: */

      /* 	SPSO2011Particle(size_t dimension) : TSuper(dimension) */
      /* 	    { */
      /* 	      xpi = new double[TSuper::dimension]; */
      /* 	      p1 = new double[TSuper::dimension]; */
      /* 	      p2 = new double[TSuper::dimension]; */
      /* 	      gr = new double[TSuper::dimension]; */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Destructor */
      /* 	   *\/ */
      /* 	  virtual ~SPSO2011Particle(void) */
      /* 	    { */
      /* 	      delete[] xpi; */
      /* 	      delete[] p1; */
      /* 	      delete[] p2; */
      /* 	      delete[] gr; */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Update the position of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = p_k + v_k */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i)); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Updates the velocity of the particle */
      /* 	   *\/ */
      /* 	  virtual void updateVelocity(void) */
      /* 	  { */
      /* 	    // The update of the velocity is done according to the equation : */
      /* 	    // v_i(t+1) = w v_i(t) + x'_i(t) - x_i(t) */
      /* 	    // where x'_i is normally sampled from the hypersphere (Gi, Ri) */
      /* 	    // with Gi the center and Ri the radius and */
      /* 	    // Gi = 1/3 ( x_i + (x_i + c(p_i - x_i)) + (x_i + c(l_i - x_i))) */
      /* 	    //    = 1/3 ( x_i +          p1          +            p2     ))) */
      /* 	    // or */
      /* 	    // Gi = 1/2 (x_i + (x_i + c(p_i - x_i))) if p_i == l_i */
      /* 	    //    = 1/2 (x_i +          p1         ) */
      /* 	    // Ri = ||G_i - x_i|| */

      /* 	    size_t i; */

      /* 	    // We first check if the local best and personal best are identical by computing the norm of the difference */
      /* 	    // between both positions */
      /* 	    // This can be done with pointer comparison */
      /* 	    // as the neighborhood holds a pointer to the best personal best */
      /* 	    bool li_equals_pi = (&(this->getBestPosition()) == this->getNeighborhood()->getBest()); */

      /* 	    // First position */
      /* 	    // p1 = xi + c * (pi - xi) */
      /* 	    // where pi is the personal best position */
      /* 	    for(i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      p1[i] = this->getPosition(i) + PARTICLE_PARAMS::c() *(this->getBestPosition().getPosition(i) - this->getPosition(i)); */
	  
      /* 	    // Second position */
      /* 	    // p2 = xi + c * (li - xi) */
      /* 	    // where li is the local best position (within the neighborhood) */
      /* 	    for(i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      p2[i] = this->getPosition(i) + PARTICLE_PARAMS::c() *(this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i)); */
	  
      /* 	    // Compute the gravity center of p1, p2 and xi */
      /* 	    if(!li_equals_pi) */
      /* 	      { */
      /* 		// We here consider p1, p2 and xi */
      /* 		for(i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  gr[i] = 1.0/3.0 * (this->getPosition(i) + p1[i] + p2[i]); */
      /* 	      } */
      /* 	    else */
      /* 	      { */
      /* 		// We here consider only p1 (or p2, since they are equal) and xi */
      /* 		for(i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  gr[i] = 1.0/2.0 * (this->getPosition(i) + p1[i]); */
      /* 	      } */

      /* 	    // And the radius of the hypersphere */
      /* 	    double ri = 0.0; */
      /* 	    for(i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      ri += (gr[i] - this->getPosition(i))*(gr[i] - this->getPosition(i)); */
      /* 	    ri = sqrt(ri); */

      /* 	    // Compute the auxiliary position x'_i randomly within the hypersphere (gr, ri); */
      /* 	    // To uniformely sample from the hypersphere we uniformely sample a direction */
      /* 	    double norm = 0.0; */
      /* 	    for(i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		xpi[i] = popot::math::normal(0.0,1.0); */
      /* 		norm += xpi[i] * xpi[i]; */
      /* 	      } */
      /* 	    norm = sqrt(norm); */
	  
      /* 	    // And then scale by a random radius */
      /* 	    double r = popot::math::uniform_random(0.0,1.0); */
      /* 	    for(i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      xpi[i] =  gr[i] + r * ri * xpi[i] / norm; */

      /* 	    // And then update the velocity */
      /* 	    for(i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setVelocity(i, PARTICLE_PARAMS::w() * this->getVelocity(i) */
      /* 				+ xpi[i] - this->getPosition(i)); */
      /* 	  } */


      /* 	  /\** */
      /* 	   * Bounds the velocity and position of the particle */
      /* 	   *\/ */
	  
      /* 	  virtual void confine(PROBLEM& p) */
      /* 	  { */
      /* 	    // In case the position is out of the bounds */
      /* 	    // we reset the velocities */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		if((this->getPosition(i) < p.get_lbound(i))) */
      /* 		  { */
      /* 		    this->setPosition(i, p.get_lbound(i)); */
      /* 		    this->setVelocity(i,-0.5*this->getVelocity(i)); */
      /* 		  } */
      /* 		else if(this->getPosition(i) > p.get_ubound(i)) */
      /* 		  { */
      /* 		    this->setPosition(i, p.get_ubound(i)); */
      /* 		    this->setVelocity(i,-0.5*this->getVelocity(i)); */
      /* 		  } */
      /* 	      } */
      /* 	  } */
	  

      /* 	}; */

      /* /\** */
      /*  * Definition of a Standard PSO 2011 particle, respecting the order of the operations of SPSO on swarm central */
      /*  * basically, in the initialization */
      /*  *\/ */
      /* template< typename PROBLEM, typename PARTICLE_PARAMS, typename POSITION_INITIALIZER=popot::PSO::initializer::PositionUniformRandom, typename VELOCITY_INITIALIZER=popot::PSO::initializer::VelocitySPSO2011> */
      /* 	class BenchSPSO2011Particle : public Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER > */
 
      /* 	{ */
      /* 	  typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> TSuper; */

      /* 	  double * xpi; */
      /* 	  double * p1; */
      /* 	  double * p2; */
      /* 	  double * gr; */
      /* 	public: */

      /* 	BenchSPSO2011Particle(void) : TSuper() */
      /* 	    { */
      /* 	      xpi = new double[TSuper::dimension]; */
      /* 	      p1 = new double[TSuper::dimension]; */
      /* 	      p2 = new double[TSuper::dimension]; */
      /* 	      gr = new double[TSuper::dimension]; */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Destructor */
      /* 	   *\/ */
      /* 	  virtual ~BenchSPSO2011Particle(void) */
      /* 	    { */
      /* 	      delete[] xpi; */
      /* 	      delete[] p1; */
      /* 	      delete[] p2; */
      /* 	      delete[] gr; */
      /* 	    } */

      /* 	  virtual void init(void) */
      /* 	  { */
      /* 	    double params[3]; */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		this->setPosition(i,POSITION_INITIALIZER::init(PROBLEM::get_lbound(i),PROBLEM::get_ubound(i))); */

      /* 		params[0] = PROBLEM::get_lbound(i); */
      /* 		params[1] = PROBLEM::get_ubound(i); */
      /* 		params[2] = this->getPosition(i); */
      /* 		this->setVelocity(i, VELOCITY_INITIALIZER::init(params)); */
      /* 	      } */

      /* 	    this->evaluateFitness(); */

      /* 	    // Set the best particle to the current position */
      /* 	    this->_best_position = *this; */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Update the position of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = p_k + v_k */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i)); */
      /* 	    if(VERBOSE_BENCH) */
      /* 	      { */
      /* 		std::cout << "New position : " ; */
      /* 		for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  std::cout << this->getPosition(i) << " "; */
      /* 		std::cout << std::endl; */
      /* 	      } */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Updates the velocity of the particle */
      /* 	   *\/ */
      /* 	  virtual void updateVelocity(void) */
      /* 	  { */

      /* 	    if(VERBOSE_BENCH) */
      /* 	      { */
      /* 		std::cout << "Position of the best informant (over " << this->getNeighborhood()->size() << ") : " << this->getNeighborhood()->getBest()->getFitness() << " " ; */
      /* 		for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  std::cout << this->getNeighborhood()->getBest()->getPosition(i) << " " ; */
      /* 		std::cout << std::endl; */

      /* 		std::cout << "Old velocity : " ; */
      /* 		for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  std::cout << this->getVelocity(i) << " "; */
      /* 		std::cout << std::endl; */
      /* 	      } */

      /* 	    // The update of the velocity is done according to the equation : */
      /* 	    // v_i(t+1) = w v_i(t) + x'_i(t) - x_i(t) */
      /* 	    // where x'_i is normally sampled from the hypersphere (Gi, Ri) */
      /* 	    // with Gi the center and Ri the radius and */
      /* 	    // Gi = 1/3 ( x_i + (x_i + c(p_i - x_i)) + (x_i + c(l_i - x_i))) */
      /* 	    //    = 1/3 ( x_i +          p1          +            p2     ))) */
      /* 	    // or */
      /* 	    // Gi = 1/2 (x_i + (x_i + c(p_i - x_i))) if p_i == l_i */
      /* 	    //    = 1/2 (x_i +          p1         ) */
      /* 	    // Ri = ||G_i - x_i|| */

      /* 	    // We first check if the local best and personal best are identical by computing the norm of the difference */
      /* 	    // between both positions */
      /* 	    // This can be done with pointer comparison */
      /* 	    // as the neighborhood holds a pointer to the best personal best */
      /* 	    bool li_equals_pi = (&(this->getBestPosition()) == this->getNeighborhood()->getBest()); */

      /* 	    if(VERBOSE_BENCH) */
      /* 	      { */
      /* 		if(li_equals_pi) */
      /* 		  std::cout << "The particle IIISSSS its own local best" << std::endl; */
      /* 	      } */

      /* 	    if(VERBOSE_BENCH) */
      /* 	      { */
      /* 		std::cout << "BW = .... nb_calls = " << RNG_GENERATOR::nb_calls << std::endl; */
      /* 	      } */

      /* 	    // First position */
      /* 	    // p1 = xi + c * (pi - xi) */
      /* 	    // where pi is the personal best position */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      p1[i] = this->getPosition(i) + PARTICLE_PARAMS::c() *(this->getBestPosition().getPosition(i) - this->getPosition(i)); */
	  
      /* 	    // Second position */
      /* 	    // p2 = xi + c * (li - xi) */
      /* 	    // where li is the local best position (within the neighborhood) */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      p2[i] = this->getPosition(i) + PARTICLE_PARAMS::c() *(this->getNeighborhood()->getBest()->getPosition(i) - this->getPosition(i)); */
	  
      /* 	    // Compute the gravity center of p1, p2 and xi */
      /* 	    if(!li_equals_pi) */
      /* 	      { */
      /* 		// We here consider p1, p2 and xi */
      /* 		for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  gr[i] = 1.0/3.0 * (this->getPosition(i) + p1[i] + p2[i]); */
      /* 	      } */
      /* 	    else */
      /* 	      { */
      /* 		// We here consider only p1 (or p2, since they are equal) and xi */
      /* 		for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  gr[i] = 1.0/2.0 * (this->getPosition(i) + p1[i]); */
      /* 	      } */

      /* 	    // And the radius of the hypersphere */
      /* 	    double ri = 0.0; */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      ri += (gr[i] - this->getPosition(i))*(gr[i] - this->getPosition(i)); */
      /* 	    ri = sqrt(ri); */

      /* 	    // Compute the auxiliary position x'_i randomly within the hypersphere (gr, ri); */
      /* 	    // To uniformely sample from the hypersphere we uniformely sample a direction */
      /* 	    double norm = 0.0; */
      /* 	    if(VERBOSE_BENCH) */
      /* 	      { */
      /* 		std::cout << "Nb RNG calls before alea normal: " << RNG_GENERATOR::nb_calls << std::endl; */
      /* 	      } */

      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		xpi[i] = popot::math::normal_spso2011(0.0,1.0); */
      /* 		norm += xpi[i] * xpi[i]; */
      /* 	      } */
      /* 	    norm = sqrt(norm); */

      /* 	    if(VERBOSE_BENCH) */
      /* 	      { */
      /* 		std::cout << "Alea normal : " ; */
      /* 		for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  std::cout << xpi[i] << " "; */
      /* 		std::cout << std::endl; */
	      
      /* 		std::cout << "Rad = " << ri << std::endl; */
      /* 	      } */

      /* 	    // And then scale by a random radius */
      /* 	    double r = popot::math::uniform_random(0.0,1.0); */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      xpi[i] =  gr[i] +  r * ri * xpi[i] / norm; */

      /* 	    // And then update the velocity */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setVelocity(i, PARTICLE_PARAMS::w() * this->getVelocity(i) */
      /* 				+ xpi[i] - this->getPosition(i)); */
      /* 	    if(VERBOSE_BENCH) */
      /* 	      { */
      /* 		std::cout << "New velocity : " ; */
      /* 		for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 		  std::cout << this->getVelocity(i) << " "; */
      /* 		std::cout << std::endl; */
      /* 		std::cout << "Nb RNG calls after new velocity: " << RNG_GENERATOR::nb_calls << std::endl; */
      /* 	      } */

      /* 	  } */


      /* 	  /\** */
      /* 	   * Bounds the velocity and position of the particle */
      /* 	   *\/ */
	  
      /* 	  virtual void confine(void) */
      /* 	  { */
      /* 	    // In case the position is out of the bounds */
      /* 	    // we reset the velocities */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		if((this->getPosition(i) < PROBLEM::get_lbound(i))) */
      /* 		  { */
      /* 		    this->setPosition(i, PROBLEM::get_lbound(i)); */
      /* 		    this->setVelocity(i,-0.5*this->getVelocity(i)); */
      /* 		  } */
      /* 		else if(this->getPosition(i) > PROBLEM::get_ubound(i)) */
      /* 		  { */
      /* 		    this->setPosition(i, PROBLEM::get_ubound(i)); */
      /* 		    this->setVelocity(i,-0.5*this->getVelocity(i)); */
      /* 		  } */
      /* 	      } */
      /* 	  } */
      /* 	}; */


      /* /\** */
      /*  * Definition of a Stochastic Standard PSO 2006 particle */
      /*  * @brief The template parameters is the problem you want to solve, */
      /*  *  and the parameters (inertia, accelaration) of the velocity update rule */
      /*  * Compared to the base type BaseParticle, this type adds the some specific position and velocity update rules */
      /*  * These come from the Standard PSO 2006 . The stochastic code simply reevaluates the fitness of the personal best */
      /*  * position before possibly changing it */
      /*  *\/ */
      /* template< typename PROBLEM, typename PARTICLE_PARAMS, typename POSITION_INITIALIZER=popot::PSO::initializer::PositionUniformRandom, typename VELOCITY_INITIALIZER=popot::PSO::initializer::VelocityHalfDiff> */
      /* 	class StochasticSPSO2006Particle : public SPSO2006Particle<PROBLEM, PARTICLE_PARAMS, POSITION_INITIALIZER, VELOCITY_INITIALIZER> */
 
      /* 	{ */
      /* 	  typedef SPSO2006Particle<PROBLEM, PARTICLE_PARAMS, POSITION_INITIALIZER, VELOCITY_INITIALIZER> TSuper; */

      /* 	public: */
      /* 	StochasticSPSO2006Particle() : TSuper() */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Updates the best position */
      /* 	   *\/ */
      /* 	  virtual void updateBestPosition(void) */
      /* 	  { */
      /* 	    // Reevalute the fitness of the personal best before */
      /* 	    // possibly changing the personal best */
      /* 	    this->_best_position.evaluateFitness(); */

      /* 	    // Update the best position the particle ever had */
      /* 	    // with a copy of the current position */
      /* 	    if(this->compare(this->getBestPosition()) < 0) */
      /* 	      this->_best_position = *this; */
      /* 	  } */
      /* 	}; */

      /* /\** */
      /*  * Definition of a Stochastic Standard PSO 2011 particle */
      /*  * @brief The template parameters is the problem you want to solve, */
      /*  *  and the parameters (inertia, accelaration) of the velocity update rule */
      /*  * Compared to the base type BaseParticle, this type adds the some specific position and velocity update rules */
      /*  * These come from the Standard PSO 2011 . The stochastic code simply reevaluates the fitness of the personal best */
      /*  * position before possibly changing it */
      /*  *\/ */
      /* template< typename PROBLEM, typename PARTICLE_PARAMS, typename POSITION_INITIALIZER=popot::PSO::initializer::PositionUniformRandom, typename VELOCITY_INITIALIZER=popot::PSO::initializer::VelocitySPSO2011> */
      /* 	class StochasticSPSO2011Particle : public SPSO2011Particle<PROBLEM, PARTICLE_PARAMS, POSITION_INITIALIZER, VELOCITY_INITIALIZER> */
 
      /* 	{ */
      /* 	  typedef SPSO2011Particle<PROBLEM, PARTICLE_PARAMS, POSITION_INITIALIZER, VELOCITY_INITIALIZER> TSuper; */

      /* 	public: */
      /* 	StochasticSPSO2011Particle() : TSuper() */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Updates the best position */
      /* 	   *\/ */
      /* 	  virtual void updateBestPosition(void) */
      /* 	  { */
      /* 	    // Reevalute the fitness of the personal best before */
      /* 	    // possibly changing the personal best */
      /* 	    this->_best_position.evaluateFitness(); */

      /* 	    // Update the best position the particle ever had */
      /* 	    // with a copy of the current position */
      /* 	    if(this->compare(this->getBestPosition()) < 0) */
      /* 	      this->_best_position = *this; */
      /* 	  } */
      /* 	}; */


      /* /\** */
      /*  * Definition of the Barebone particles */
      /*  * @brief The template parameters is the problem you want to solve, */
      /*  *  and the parameters (inertia, accelaration) of the velocity update rule */
      /*  *\/ */
      /* template< typename PROBLEM, typename POSITION_INITIALIZER=popot::PSO::initializer::PositionUniformRandom, typename VELOCITY_INITIALIZER=popot::PSO::initializer::VelocityHalfDiff> */
      /* 	class BareboneParticle : public Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> */
      /* 	{ */
      /* 	  typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> TSuper; */

      /* 	public: */

      /* 	BareboneParticle(void) : TSuper() */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Destructor */
      /* 	   *\/ */
      /* 	  virtual ~BareboneParticle(void) */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Update the position of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = v_{k+1} */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i , this->getVelocity(i)); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Updates the velocity of the particle */
      /* 	   *\/ */
      /* 	  virtual void updateVelocity(void) */
      /* 	  { */
      /* 	    // The update of the velocity is done according to the equation : */
      /* 	    // v_{k+1} = N( (best_p+best_g)/2), sigma) */
      /* 	    // with : */
      /* 	    // sigma : the variance defined as |best_p_ij - best_g_ij| with */
      /* 	    // best_p : the best position the particle ever had */
      /* 	    // best_g : the best position the neighborhood ever had */

      /* 	    double mean,var; */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		mean = 0.5*(this->getBestPosition().getPosition(i) + this->getNeighborhood()->getBest()->getPosition(i)); */
      /* 		var = fabs(this->getBestPosition().getPosition(i) - this->getNeighborhood()->getBest()->getPosition(i)); */
      /* 		this->setVelocity(i, popot::math::normal(mean,var)); */
      /* 	      } */
      /* 	  } */
      /* 	}; */

      /* /\** */
      /*  * Definition of a modified barebone particle */
      /*  * @brief The template parameters is the problem you want to solve, */
      /*  *  and the parameters (inertia, accelaration) of the velocity update rule */
      /*  * In this implementation, the center of the normal distribution is the mean between the personal best */
      /*  * and local best positions (other variations exist, using a random personal best in the neighbor instead */
      /*  * of the local best for example) */
      /*  *\/ */
      /* template< typename PROBLEM , typename POSITION_INITIALIZER=popot::PSO::initializer::PositionUniformRandom, typename VELOCITY_INITIALIZER=popot::PSO::initializer::VelocityHalfDiff> */
      /* 	class ModifiedBareboneParticle : public Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> */
      /* 	{ */
      /* 	  typedef Particle<PROBLEM, POSITION_INITIALIZER, VELOCITY_INITIALIZER> TSuper; */
	
      /* 	public: */

      /* 	ModifiedBareboneParticle(void) : TSuper() */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Destructor */
      /* 	   *\/ */
      /* 	  virtual ~ModifiedBareboneParticle(void) */
      /* 	    { */
      /* 	    } */

      /* 	  /\** */
      /* 	   * Update the position of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = v_{k+1} */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i, this->getVelocity(i)); */
      /* 	  } */

      /* 	  /\** */
      /* 	   * Updates the velocity of the particle */
      /* 	   *\/ */
      /* 	  virtual void updateVelocity(void) */
      /* 	  { */
      /* 	    // The update of the velocity is done according to the equation : */
      /* 	    // v_{k+1} = N( (best_p+best_g)/2), sigma) */
      /* 	    // with : */
      /* 	    // sigma : the variance defined as |best_p_ij - best_g_ij| with */
      /* 	    // best_p : the best position the particle ever had */
      /* 	    // best_g : the best position the neighborhood ever had */
      /* 	    double mean,var; */
	        
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      { */
      /* 		if(popot::math::uniform_random(0.0,1.0) < 0.5) */
      /* 		  { */
      /* 		    this->setVelocity(i, this->getBestPosition().getPosition(i)); */
      /* 		  } */
      /* 		else */
      /* 		  { */
      /* 		    mean = 0.5*(this->getBestPosition().getPosition(i) + this->getNeighborhood()->getBest()->getPosition(i)); */
      /* 		    var = fabs(this->getBestPosition().getPosition(i) - this->getNeighborhood()->getBest()->getPosition(i)); */
      /* 		    this->setVelocity(i, popot::math::normal(mean,var)); */
      /* 		  } */
      /* 	      } */
      /* 	  } */
      /* 	}; // Modified Barebone */
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
