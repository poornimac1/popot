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

      /* 	  /** */
      /* 	   * Destructor */
      /* 	   */
      /* 	  virtual ~SPSO2006Particle(void) */
      /* 	    { */
      /* 	    } */

      /* 	  /** */
      /* 	   * Update the position of the particle */
      /* 	   */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = p_k + v_k */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i)); */
      /* 	  } */

      /* 	  /** */
      /* 	   * Updates the velocity of the particle */
      /* 	   */
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

      /* 	  /** */
      /* 	   * Destructor */
      /* 	   * */
      /* 	  virtual ~SPSO2007Particle(void) */
      /* 	    { */
      /* 	    } */

      /* 	  /**
      /* 	   * Update the position of the particle */
      /* 	   */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = p_k + v_k */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i)); */
      /* 	  } */

      /* 	  /** */
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

      /* 	  /** */
      /* 	   * Destructor */
      /* 	   *\/ */
      /* 	  virtual ~SPSO2011Particle(void) */
      /* 	    { */
      /* 	      delete[] xpi; */
      /* 	      delete[] p1; */
      /* 	      delete[] p2; */
      /* 	      delete[] gr; */
      /* 	    } */

      /* 	  /** */
      /* 	   * Update the position of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = p_k + v_k */
      /* 	    for(size_t i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i, this->getPosition(i) + this->getVelocity(i)); */
      /* 	  } */

      /* 	  /** */
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


      /* 	  /** */
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

      /* /** */
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

      /* 	  /** */
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

      /* 	  /** */
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

      /* 	  /** */
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


      /* 	  /** */
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


      /* /** */
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

      /* 	  /** */
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

      /* /** */
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

      /* 	  /** */
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


      /* /** */
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

      /* 	  /** */
      /* 	   * Destructor */
      /* 	   *\/ */
      /* 	  virtual ~BareboneParticle(void) */
      /* 	    { */
      /* 	    } */

      /* 	  /** */
      /* 	   * Update the position of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = v_{k+1} */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i , this->getVelocity(i)); */
      /* 	  } */

      /* 	  /** */
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

      /* /** */
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

      /* 	  /** */
      /* 	   * Destructor */
      /* 	   *\/ */
      /* 	  virtual ~ModifiedBareboneParticle(void) */
      /* 	    { */
      /* 	    } */

      /* 	  /** */
      /* 	   * Update the position of the particle */
      /* 	   *\/ */
      /* 	  virtual void updatePosition(void) */
      /* 	  { */
      /* 	    // Here it is simply : p_{k+1} = v_{k+1} */
      /* 	    for(int i = 0 ; i < TSuper::dimension ; ++i) */
      /* 	      this->setPosition(i, this->getVelocity(i)); */
      /* 	  } */

      /* 	  /** */
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
