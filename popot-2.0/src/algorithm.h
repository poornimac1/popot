#ifndef POPOT_ALGORITHM_H
#define POPOT_ALGORITHM_H

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

#include "neighborhood.h"
#include "individuals.h"
#include "topology.h"
#include <fstream>


// TODO
// --- rendre init générique <-- mais comment ? les lambda fonction sont des types particuliers, je n'arrive pas à le
//   mettre dans le decltype ou type de retour de mon spso2006(...), et pour le coup les adapter en fonction de l'algo
// --- paramétrer les fonctions de confinenemt : les confinements sont différents pour la 2006, 2007 et 2011

namespace popot
{
  namespace PSO
  {
    namespace algorithm
    {
      typedef enum
	{
	  SYNCHRONOUS_EVALUATION,
	  ASYNCHRONOUS_EVALUATION,
	  ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION
	} EvaluationMode;

      /**
       * Standard particle swarm optimization
       * @param PARTICLE The particle type you want to use, e.g. swarm::particle::TraditionalParticle
       * @param TOPOLOGY Defines the topology, e.g. swarm::topology::VonNeuman
       * @param STOP_CRITERIA Defines the condition to stop the evolution of the swarm
       */

      template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION,
		typename TOPOLOGY, typename UPDATE_POSITION_RULE, typename UPDATE_VELOCITY_RULE, typename UPDATE_BEST_POSITION_RULE, 
		typename PARTICLE>
      class Base
      {
      public:
	typedef typename PARTICLE::BestType BestType;
	typedef PARTICLE  ParticleType;
	typedef typename PARTICLE::NeighborhoodType NeighborhoodType;

      private:
	std::vector<PARTICLE> particles;
	std::vector<NeighborhoodType * > neighborhoods;
	std::map< size_t , std::vector<size_t> > neighborhood_membership;

	size_t *particles_indexes;

	BestType _best_particle;

	size_t _dimension;
	size_t _swarm_size;

	const LBOUND_FUNC& _lbound;
	const UBOUND_FUNC& _ubound;
	const STOP_CRITERIA& _stop_criteria;
	const COST_FUNCTION& _cost_function;
	const TOPOLOGY _topology;
	const UPDATE_POSITION_RULE _update_position_rule;
	//const INIT_FUNCTION& _init_function;
	const UPDATE_VELOCITY_RULE _update_velocity_rule;
	const UPDATE_BEST_POSITION_RULE _update_best_position_rule;

	const bool _reevaluate_best_before_updating;

	EvaluationMode _evaluation_mode;
      public:
	int epoch;
	int nb_new_neigh;

      public:
	Base(size_t swarm_size,
	     size_t dimension,
	     const LBOUND_FUNC& lbound,
	     const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop,
	     const COST_FUNCTION& cost_function,    
	     const TOPOLOGY topology,
	     const UPDATE_POSITION_RULE update_position_rule,
	     const UPDATE_VELOCITY_RULE update_velocity_rule,
	     const UPDATE_BEST_POSITION_RULE update_best_position_rule,
	     const PARTICLE& p,
	     EvaluationMode evaluation_mode,
	     bool reevaluate_best_before_updating) 
	  : _best_particle(dimension),
	    _dimension(dimension),
	    _swarm_size(swarm_size),
	    _lbound(lbound),
	    _ubound(ubound),
	    _stop_criteria(stop),
	    _cost_function(cost_function),
	    _topology(topology),
	    _update_position_rule(update_position_rule),
	    _update_velocity_rule(update_velocity_rule),
	    _update_best_position_rule(update_best_position_rule),
	    _reevaluate_best_before_updating(reevaluate_best_before_updating),
	    _evaluation_mode(evaluation_mode),
	    epoch(0),
	    nb_new_neigh(0)
	{
	  //
	  particles_indexes = new size_t[swarm_size];

	  // Declare our particles
	  particles.clear();
	  for(size_t i = 0 ; i < swarm_size ; ++i)
	    {
	      // Create the particle
	      particles.push_back(PARTICLE(_dimension));

	      // Initialize the position, best position and velocity
	      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      // TODO :  This part must be added in the template !!!!
	      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	      //_position_initializer(particles[i]);

	      popot::initializer::position::uniform_random<typename PARTICLE::VECTOR_TYPE>(particles[i].getPosition(), lbound, ubound);
	      popot::initializer::velocity::half_diff<typename PARTICLE::VECTOR_TYPE>(particles[i].getPosition(), particles[i].getVelocity(), lbound, ubound);
				       
	      // Evaluate the initial fitnesses
	      particles[i].evaluateFitness(_cost_function);
		
	      // And initialize the best position
	      particles[i].initBestPosition();

	      // Set up the array of indices, used to random shuffle 
	      particles_indexes[i] = i;
	    }
	    
	  // We now form groups of particles depending on the topology
	  _topology(particles, neighborhoods, neighborhood_membership);
	  nb_new_neigh ++;

	  // We now browse all the neighborhoods and find the best particles
	  // within each of them
	  // and initialize the best best particle of the whole swarm
	  _best_particle = *(neighborhoods[0]->findBest());
	  for(size_t i = 1 ; i < neighborhoods.size() ; ++i)
	    if(neighborhoods[i]->findBest()->compare(_best_particle) < 0)
	      _best_particle = *(neighborhoods[i]->getBest());
	}

	virtual ~Base(void)
	{
	  delete[] particles_indexes;
	}


	/**
	 * Returns the size of the swarm
	 */
	size_t getSize(void)
	{
	  return _swarm_size;
	}

	/**
	 * One step of the swarm
	 */
	double step(void)
	{
	  // For each particle,
	  // - Update the particle's velocity
	  // - Enforce velocity boundaries
	  // - Move the particle to its new position
	  // - Enforce position boundaries 
	  // - Update the particle's best position
	  // - Update the swarm's best position

	  // if(VERBOSE_BENCH)
	  //   std::cout << "Random before looping : " << RNG_GENERATOR::nb_calls << " calls " << std::endl;

	  // If we use an asynchronous update, we first shuffle the particles
	  if(_evaluation_mode == ASYNCHRONOUS_EVALUATION || _evaluation_mode == ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION)
	    {
	      // If comparing to SPSO, the following should be commented
	      if(_evaluation_mode == ASYNCHRONOUS_EVALUATION)
	      	popot::math::random_shuffle_indexes(particles_indexes, _swarm_size);

	      size_t particle_index;
	      for(size_t i = 0 ; i < _swarm_size ; ++i)
	   	{
	   	  particle_index = particles_indexes[i];

	   	  // First find the best informant within the neighborhood
	   	  particles[particle_index].getNeighborhood().findBest();

	   	  // Update the velocities and position
		  _update_velocity_rule(particles[particle_index]);
	   	  _update_position_rule(particles[particle_index]);

	   	  // Confine the positions and velocities if required
	   	  particles[particle_index].confine(_lbound, _ubound);

	   	  // Compute the fitness of the new position
	   	  particles[particle_index].evaluateFitness(_cost_function);

	   	  // And see if we update the personal best
		  if(_reevaluate_best_before_updating)
		    particles[i].getBestPosition().evaluateFitness(_cost_function);

		  _update_best_position_rule(particles[particle_index]);


	   	  if(VERBOSE_BENCH)
	   	    std::cout << std::endl;

	   	}
	    }
	  else // Synchronous evaluation
	    {
	      // In synchronous mode
	      // we first update all the current positions and evaluate their fitness
	      for(size_t i = 0 ; i < _swarm_size ; ++i)
	  	{
	   	  // Update the velocities and position
		  _update_velocity_rule(particles[i]);
	   	  _update_position_rule(particles[i]);

	   	  // Confine the positions and velocities if required
	   	  particles[i].confine(_lbound, _ubound);

	   	  // Compute the fitness of the new position
	   	  particles[i].evaluateFitness(_cost_function);
	  	}
	      // before changing the personal bests


	      // The personal best position is updated after all the positions of the particles
	      // have been updated because the neighborhoods hold a pointer to a personal best
	      // As this best in the neighborhood is used by the particles when they
	      // update their velocity, we must ensure that, for a synchronous evaluation
	      // updating the best is done after all the positions updates
	      for(size_t i = 0 ; i < _swarm_size ; ++i)
		{
		  if(_reevaluate_best_before_updating)
		    particles[i].getBestPosition().evaluateFitness(_cost_function);
		  _update_best_position_rule(particles[i]);
		}

	      // We now update the best particle of all the neighborhoods
	      for(size_t i = 0 ; i < neighborhoods.size() ; ++i)
	  	neighborhoods[i]->findBest();
	    }

	  // Update the best particle the whole swarm ever had
	  double old_fitness = _best_particle.getFitness();

	  _best_particle = *(neighborhoods[0]->findBest());
	  for(size_t i = 1 ; i < neighborhoods.size() ; ++i)
	    if(neighborhoods[i]->findBest()->compare(_best_particle) < 0)
	      _best_particle = *(neighborhoods[i]->getBest());

	  if(_best_particle.getFitness() >= old_fitness)
	    {
	      // We consider that there is no improvement in the best particle
	      // and ask the topology if it wants to regenerate its topology
	      _topology(particles, neighborhoods, neighborhood_membership);
	      nb_new_neigh ++;
	    }
	  
	  epoch++;
	  return _best_particle.getFitness() - old_fitness;
	}

	/**
	 * Current epoch
	 */
	size_t getEpoch(void) const
	{
	  return epoch;
	}

	/**
	 * Run until stop criteria is met. 
	 * @param verbose set it to 1 to print the best fitness at each iteration
	 */
	void run(int verbose=0)
	{
	  while(!_stop_criteria(getBest().getFitness(), epoch))
	    {
	      step();
	      if(verbose) std::cout << '\r' << std::setw(6) << std::setfill('0') << epoch << " " << getBest().getFitness() << std::setw(5) << std::setfill(' ') << ' ' << std::flush;
	    }
	  if(verbose) std::cout << std::endl;
	}

	/**
	 * Returns the current best particle of the whole swarm (best of the personal best)
	 */
	BestType& getBest(void)
	{
	  return _best_particle;
	}

	/**
	 * Returns the fitness of the best individual (Required for benchmarking)
	 */ 
	double getBestFitness() const
	{
	  return _best_particle.getFitness();
	}

	/**
	 * Returns a reference to the vector of particles
	 */
	std::vector<PARTICLE>& getParticles(void)
	{
	  return particles;
	}

	/**
	 * Print in the console the status of the swarm
	 * @param mode 0 displays all the particles, 1 display only the best
	 */
	void print(int mode=0)
	{
	  if(_evaluation_mode == SYNCHRONOUS_EVALUATION)
	    printf("Using a synchronous evaluation mode \n");
	  else if(_evaluation_mode == ASYNCHRONOUS_EVALUATION)
	    printf("Using an asynchronous evaluation mode \n");
	  else if(_evaluation_mode == ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION)
	    printf("Using an asynchronous evaluation mode with random shuffling the particles\n");
	  else
	    printf("WARNING : Unrecognized evaluation mode \n");

	  switch(mode)
	    {
	    case 0:
	      {
		for(size_t i = 0 ; i < _swarm_size ; ++i)
		  {
		    std::cout << "Particle " << std::setw(3) << std::setfill(' ') << i
			      << " : " << "\n" << particles[i] << std::endl;
		  }
		// Display the best particle
		std::cout << "Best particle : \n" << _best_particle << std::endl;
	      }
	      break;
	    case 1:
	      {
		// Display the best particle
		std::cout << "Best particle : \n" << _best_particle << std::endl;
	      }
	      break;
	    case 2:
	      {
		for(size_t i = 0 ; i < _swarm_size ; ++i)
		  {
		    std::cout << "Particle " << std::setw(3) << std::setfill(' ') << i
			      << " : " << particles[i].getFitness() << " ; Best : " << particles[i].getBestPosition().getFitness() << std::endl;
		  }
		// Display the best particle
		std::cout << "Best particle : \n" << _best_particle.getFitness() << std::endl;
	      }
	      break;
	    default:
	      break;
	    }
	}

	
	std::map< size_t , std::vector<size_t> >& getNeighborhoodMembership(void)
	{
	  return neighborhood_membership;
	}

	/**
	 * Generate a graph of the connections between the particles using DOT
	 */
	void generateGraph(std::string filename)
	{
	  printf("----------- \n Generating a graph of the connectivity \n");
	  size_t nb_particles = _swarm_size;

	  // From this connectivity matrix, we generate the file
	  // to be used with the DOT command
	  std::ofstream outfile(filename.c_str());
	  outfile << "digraph Connections {" << std::endl;
	  outfile << "node [shape=circle,fixedsize=true,width=0.6];  ";
	  for(size_t i = 0 ; i < nb_particles ; ++i)
	    {
	      outfile << i << ";";
	    }
	  outfile << std::endl;

	  // Generate a connectivity matrix
	  int *connectivity_matrix = new int[nb_particles*nb_particles];
	  for(size_t i = 0 ; i < nb_particles*nb_particles ; ++i)
	    connectivity_matrix[i] = 0;

	  for(size_t i = 0 ; i < nb_particles ; ++i)
	    for(size_t j = 0 ; j < neighborhood_membership[i].size(); ++j)
	      connectivity_matrix[i*nb_particles+neighborhood_membership[i][j]] += 1;

	  // Put in the connections
	  for(size_t i = 0 ; i < nb_particles ; ++i)
	    {
	      for(size_t j = 0 ; j < i ; ++j)
		{
		  if(connectivity_matrix[i*nb_particles+j] && !connectivity_matrix[j*nb_particles+i])
		    outfile << j << "->" << i << ";" << std::endl;
		  else if(connectivity_matrix[i*nb_particles+j] && connectivity_matrix[j*nb_particles+i])
		    outfile << j << "->" << i << " [dir=both];" << std::endl;
		}
	      if(connectivity_matrix[i*nb_particles+i])
		outfile << i << "->" << i << ";" << std::endl;
	    }

	  // The title and the end of the graphic
	  outfile << "overlap=false" << std::endl;
	  outfile << "label=\"Swarm's connections\"" << std::endl;
	  outfile << "fontsize=12;" << std::endl;
	  outfile << "}" << std::endl;
	  outfile.close();

	  printf("You can generate the graph by calling e.g.: \n");
	  printf(" neato -Tpng %s > graph.png\n", filename.c_str());
	  printf("Other export formats can be used with dot/neato : svg, eps, dia, xfig, ...\n");

	  // Compute the distribution of the informants
	  // We first compute the number of informants for each particle
	  // The informants of particle i are at line i
	  size_t * nb_informants = new size_t[nb_particles];
	  size_t max_nb_informants = 0;
	  for(size_t i = 0 ; i < nb_particles ; ++i)
	    {
	      nb_informants[i] = 0;
	      for(size_t j = 0 ; j < nb_particles ; ++j)
		nb_informants[i] += connectivity_matrix[i*nb_particles+j];
	      if(nb_informants[i] > max_nb_informants)
		max_nb_informants = nb_informants[i];
	    }
	  printf("\n");
	  printf("Max number of informants : %i \n", max_nb_informants);

	  outfile.open((filename + ".conn").c_str());
	  size_t nb_part = 0;
	  for(size_t i = 0 ; i <= max_nb_informants ; ++i)
	    {
	      // How many particles have nb_informants <= i ?
	      nb_part = 0;
	      for(size_t j = 0 ; j < nb_particles ; ++j)
		if(nb_informants[j] == i)
		  nb_part ++;
	      outfile << i << "\t" << nb_part/double(nb_particles) << std::endl;
	    }
	  outfile.close();
	  printf("-------------\n");

	  delete[] connectivity_matrix;
	  delete[] nb_informants;
	}
      };

      template< typename LBOUND_FUNC, typename UBOUND_FUNC,
		typename STOP_CRITERIA, typename COST_FUNCTION,
		typename TOPOLOGY,
		typename UPDATE_POSITION_RULE, typename UPDATE_VELOCITY_RULE,typename UPDATE_BEST_POSITION_RULE, 
		typename PARTICLE>
      Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION, 
	   TOPOLOGY, UPDATE_POSITION_RULE, UPDATE_VELOCITY_RULE, UPDATE_BEST_POSITION_RULE, PARTICLE>
      base(size_t swarm_size,
	   size_t dimension,
	   const LBOUND_FUNC& lbound,
	   const UBOUND_FUNC& ubound,
	   const STOP_CRITERIA& stop,
	   const COST_FUNCTION& cost_function,
	   const TOPOLOGY& topology,
	   const UPDATE_POSITION_RULE& update_position_rule,
	   const UPDATE_VELOCITY_RULE& update_velocity_rule,
	   const UPDATE_BEST_POSITION_RULE& update_best_position_rule,
	   const PARTICLE& p,
	   EvaluationMode evaluation_mode,
	   bool reevaluate_best_before_updating)
      {
	return Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION, 
		    TOPOLOGY, UPDATE_POSITION_RULE, UPDATE_VELOCITY_RULE, UPDATE_BEST_POSITION_RULE, PARTICLE>
	  (swarm_size, dimension, lbound, ubound, stop, cost_function, topology, 
	   update_position_rule, update_velocity_rule, update_best_position_rule, p, evaluation_mode, reevaluate_best_before_updating);
      }


    } // namespace algorithm
  } // namespace PSO

  namespace ABC
  {
    namespace algorithm
    {
      template<typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
      class Base
      {
	typedef popot::ABC::individuals::FoodSource FoodSourceType;

      private:
	size_t _epoch;
	size_t _CS;
	size_t _dimension;
	size_t _limitForScout;
	size_t _nb_employed;
	size_t _nb_onlookers;
	double * _probabilities;

	FoodSourceType * _foodSources;
	FoodSourceType _bestSource;

	const LBOUND_FUNC& _lbound;
	const UBOUND_FUNC& _ubound;
	const STOP_CRITERIA& _stop_criteria;
	const COST_FUNCTION& _cost_function;

	
      private:
	void findBestSource(void)
	{
	  double bestf = _bestSource.getFitness();
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    if(bestf < _foodSources[i].getFitness())
	      {
		bestf = _foodSources[i].getFitness();
		_bestSource = _foodSources[i];
	      }
	}

	void employedPhase(void)
	{
	  // In the employed bees phase, we generate new solutions
	  // around each nectar source
	  //int change_dim;
	  size_t other_source;
	  //double phi;
	  FoodSourceType new_source;
	  //double new_param_value;
	  double sum_fitnesses = 0;
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    {
	      // Randomly select another source, different from the current source
	      other_source = (size_t) popot::math::uniform_random(0, _nb_employed);
	      while(other_source == i)
		other_source = (size_t) popot::math::uniform_random(0, _nb_employed);

	      _foodSources[i].combine(_foodSources[other_source], _lbound, _ubound, _cost_function);

	      _probabilities[i] = (_foodSources[i]).getFitness();
	      sum_fitnesses += _probabilities[i];
	    }

	  // At the end, we normalize the probabilities for each source
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    _probabilities[i] /= sum_fitnesses;
	}

	void onlookerPhase(void)
	{
	  // Onlooker phase
	  // Each onlooker bee selects a food source
	  // based on its probability (reflecting its relative fitness)
	  size_t selected_source=0;
	  size_t other_source=0;
	  for(size_t i = 0 ; i < _nb_onlookers ; ++i)
	    {
	      // Select a source based on its fitness
	      selected_source = popot::math::random_from_array(_probabilities);

	      // Randomly select another source, different from the current source
	      other_source = (size_t) popot::math::uniform_random(0, _nb_employed);
	      while(other_source == i)
		other_source = (size_t) popot::math::uniform_random(0, _nb_employed);

	      _foodSources[selected_source].combine(_foodSources[other_source], _lbound, _ubound, _cost_function);
	      
	    }
	}

	void scoutPhase(void)
	{
	  // Scout phase
	  // We browse all the food sources
	  // If a source has a counter higher than the limit
	  // we reset it
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    if(_foodSources[i].getCounter() >= _limitForScout)
	      _foodSources[i].init(_lbound, _ubound, _cost_function);
	}

      public:

	Base(int colony_size, int dimension, const LBOUND_FUNC &lbound, const UBOUND_FUNC &ubound, const STOP_CRITERIA &stop_criteria, const COST_FUNCTION &cost_function)
	  : _epoch(0), 
	    _CS(colony_size), 
	    _dimension(dimension), 
	    _limitForScout(colony_size * dimension / 2), 
	    _nb_employed(colony_size/2), 
	    _nb_onlookers(colony_size/2), 
	    _probabilities(0), 
	    _foodSources(0), 
	    _bestSource(), 
	    _lbound(lbound), 
	    _ubound(ubound), 
	    _stop_criteria(stop_criteria), 
	    _cost_function(cost_function)
	{
	  // Initialize our populations
	  _foodSources = new FoodSourceType[_nb_employed];
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    _foodSources[i] = FoodSourceType(_dimension);

	  // And the probabilities of their solutions
	  _probabilities = new double[_nb_employed];

	}

	virtual ~Base(void)
	{
	  delete[] _foodSources;
	  delete[] _probabilities;
	}

	void init(void)
	{
	  // Initialize the positions and fitnesses
	  // of the employed bees
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    _foodSources[i].init(_lbound, _ubound, _cost_function);

	  _bestSource = _foodSources[0];

	  _epoch = 0;

	  // Keep track of the best solution
	  findBestSource();
	}

	int getEpoch(void)
	{
	  return _epoch;
	}

	void step(void)
	{
	  // Employed bees phase
	  employedPhase();
	  
	  // Onlooker bees phase
	  onlookerPhase();

	  // Scout bees phase
	  scoutPhase();

	  // Memorize the best solution
	  findBestSource();

	  _epoch ++;
	}

	FoodSourceType& getBest()
	{
	  return _bestSource;
	}

	void run(void)
	{
	  while(!_stop_criteria(_bestSource.getFValue(),_epoch))
	    {
	      step();
	    }
	}

	void print(void)
	{
	  std::cout << "Artificial Bee Colony " << std::endl;
	  std::cout << "Food Sources : " << std::endl;
	  for(size_t i = 0 ; i < _nb_employed ; ++i)
	    std::cout << _foodSources[i] << std::endl;
	  std::cout << "Best source : " << _bestSource << std::endl;

	}
      }; // class Base
    } // namespace algorithm
  } // namespace ABC


  namespace algorithm
  {

    /**
     * ABC algorithm
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
    popot::ABC::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION> 
    abc(size_t colony_size, size_t dimension,
	const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	const STOP_CRITERIA& stop, const COST_FUNCTION& func)  {
      return popot::ABC::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION>(colony_size, dimension, lbound, ubound, stop, func);
    };

    // The parameters for updating the velocity of the particles
    class SPSO2006_Params
    {
    public:
      static double w(void) { return 1.0/(2.0*log(2.0));};
      static double c(void) { return 0.5 + log(2.0);};
    };
    

    // The particle type, shortcut simplifying the expression of the types in spso2006
    typedef popot::PSO::particle::Particle<> ParticleSPSO;

    /**
     * Builds the SPSO2006 algorithm
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
    popot::PSO::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION, 
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				ParticleSPSO>
    spso2006(size_t dimension,
	     const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
	     const STOP_CRITERIA& stop, const COST_FUNCTION& cost_function) {
      size_t swarm_size = 10 + int(2.0 * sqrt(dimension));

      // Particle type
      ParticleSPSO p;

      // Position and velocity updates
      auto position_update = popot::PSO::particle::updatePosition<ParticleSPSO>;
      auto velocity_update = popot::PSO::particle::updateVelocity_spso2006<ParticleSPSO, SPSO2006_Params>;
     
      // Initialization functions
      
      auto init_function = [lbound, ubound] (ParticleSPSO& p) -> void { 
	popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), lbound, ubound);
	popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), p.getVelocity(), lbound, ubound);
      };
      

      // The rule to update the best position
      auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>;


      // Topology
      //auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::ring_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::vonNeuman_fillNeighborhoods<ParticleSPSO>;
      auto topology = popot::PSO::topology::randomInformants_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::adaptiveRandom_fillNeighborhoods<ParticleSPSO>;


      auto algo = popot::PSO::algorithm::base(swarm_size, dimension, 
					      lbound, ubound, stop, cost_function, 
					      topology, position_update, velocity_update, best_position_update, 
					      p, popot::PSO::algorithm::ASYNCHRONOUS_WITHOUT_SHUFFLE_EVALUATION, false); // This is the only difference from spso2006
      return algo;
    }


    /**
     * Builds the "stochastic" SPSO2006 algorithm : we just reevaluate the best position before updating it
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
    popot::PSO::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION, 
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				ParticleSPSO>
    stochastic_spso2006(size_t dimension,
			const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
			const STOP_CRITERIA& stop, const COST_FUNCTION& cost_function) {
      size_t swarm_size = 10 + int(2.0 * sqrt(dimension));

      // Particle type
      ParticleSPSO p;

      // Position and velocity updates
      auto position_update = popot::PSO::particle::updatePosition<ParticleSPSO>;
      auto velocity_update = popot::PSO::particle::updateVelocity_spso2006<ParticleSPSO, SPSO2006_Params>;
     
      // Initialization functions
      /*
	auto init_function = [lbound, ubound] (ParticleSPSO& p) -> void { 
	popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), lbound, ubound);
	popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), p.getVelocity(), lbound, ubound);
	};
      */

      // The rule to update the best position
      auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>;


      // Topology
      //auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::ring_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::vonNeuman_fillNeighborhoods<ParticleSPSO>;
      auto topology = popot::PSO::topology::randomInformants_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::adaptiveRandom_fillNeighborhoods<ParticleSPSO>;


      auto algo = popot::PSO::algorithm::base(swarm_size, dimension, 
					      lbound, ubound, stop, cost_function, 
					      topology, position_update, velocity_update, best_position_update, 
					      p, popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION, true);
      return algo;
    }



    /**
     * Builds the SPSO2007 algorithm
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
    popot::PSO::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION, 
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				ParticleSPSO>
    spso2007(size_t dimension,
			const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
			const STOP_CRITERIA& stop, const COST_FUNCTION& cost_function) {
      size_t swarm_size = 10 + int(2.0 * sqrt(dimension));

      // Particle type
      ParticleSPSO p;

      // Position and velocity updates
      auto position_update = popot::PSO::particle::updatePosition<ParticleSPSO>;
      auto velocity_update = popot::PSO::particle::updateVelocity_spso2007<ParticleSPSO, SPSO2006_Params>;
     
      // Initialization functions
      /*
	auto init_function = [lbound, ubound] (ParticleSPSO& p) -> void { 
	popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), lbound, ubound);
	popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), p.getVelocity(), lbound, ubound);
	};
      */

      // The rule to update the best position
      auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>;


      // Topology
      //auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::ring_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::vonNeuman_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::randomInformants_fillNeighborhoods<ParticleSPSO>;
      auto topology = popot::PSO::topology::adaptiveRandom_fillNeighborhoods<ParticleSPSO>;


      auto algo = popot::PSO::algorithm::base(swarm_size, dimension, 
					      lbound, ubound, stop, cost_function, 
					      topology, position_update, velocity_update, best_position_update, 
					      p, popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION, false);
      return algo;
    }
    
    /**
     * Builds the SPSO2011 algorithm
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
    popot::PSO::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION, 
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				ParticleSPSO>
    spso2011(size_t dimension,
			const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
			const STOP_CRITERIA& stop, const COST_FUNCTION& cost_function) {
      size_t swarm_size = 40;

      // Particle type
      ParticleSPSO p;

      // Position and velocity updates
      auto position_update = popot::PSO::particle::updatePosition<ParticleSPSO>;
      auto velocity_update = popot::PSO::particle::updateVelocity_spso2011<ParticleSPSO, SPSO2006_Params>;
     
      // Initialization functions
      /*
	auto init_function = [lbound, ubound] (ParticleSPSO& p) -> void { 
	popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), lbound, ubound);
	popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), p.getVelocity(), lbound, ubound);
	};
      */

      // The rule to update the best position
      auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>;


      // Topology
      //auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::ring_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::vonNeuman_fillNeighborhoods<ParticleSPSO>;
      //auto topology = popot::PSO::topology::randomInformants_fillNeighborhoods<ParticleSPSO>;
      auto topology = popot::PSO::topology::adaptiveRandom_fillNeighborhoods<ParticleSPSO, 3>;


      auto algo = popot::PSO::algorithm::base(swarm_size, dimension, 
					      lbound, ubound, stop, cost_function, 
					      topology, position_update, velocity_update, best_position_update, 
					      p, popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION, false);
      return algo;
    }

    /**
     * Builds the barebone algorithm
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
    popot::PSO::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION, 
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				ParticleSPSO>
    barebone(size_t dimension,
			const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
			const STOP_CRITERIA& stop, const COST_FUNCTION& cost_function) {
      size_t swarm_size = 40;

      // Particle type
      ParticleSPSO p;

      // Position and velocity updates
      auto position_update = popot::PSO::particle::updatePosition_barebone<ParticleSPSO>;
      auto velocity_update = popot::PSO::particle::updateVelocity_barebone<ParticleSPSO>;
     
      // Initialization functions
      /*
	auto init_function = [lbound, ubound] (ParticleSPSO& p) -> void { 
	popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), lbound, ubound);
	popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), p.getVelocity(), lbound, ubound);
	};
      */

      // The rule to update the best position
      auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>;

      // Topology
      auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>;

      auto algo = popot::PSO::algorithm::base(swarm_size, dimension, 
					      lbound, ubound, stop, cost_function, 
					      topology, position_update, velocity_update, best_position_update, 
					      p, popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION, false);
      return algo;
    }


    /**
     * Builds the barebone algorithm
     */
    template< typename LBOUND_FUNC, typename UBOUND_FUNC, typename STOP_CRITERIA, typename COST_FUNCTION>
    popot::PSO::algorithm::Base<LBOUND_FUNC, UBOUND_FUNC, STOP_CRITERIA, COST_FUNCTION, 
				void(*)(std::vector<ParticleSPSO >&, 
					std::vector< typename ParticleSPSO::NeighborhoodType *> &, 
					std::map< size_t, std::vector<size_t> > &),
				void(*)(ParticleSPSO&),
				void(*)(ParticleSPSO&), 
				void(*)(ParticleSPSO&),
				ParticleSPSO>
    modified_barebone(size_t dimension,
			const LBOUND_FUNC& lbound, const UBOUND_FUNC& ubound,
			const STOP_CRITERIA& stop, const COST_FUNCTION& cost_function) {
      size_t swarm_size = 40;

      // Particle type
      ParticleSPSO p;

      // Position and velocity updates
      auto position_update = popot::PSO::particle::updatePosition_barebone<ParticleSPSO>;
      auto velocity_update = popot::PSO::particle::updateVelocity_modifiedBarebone<ParticleSPSO>;
     
      // Initialization functions
      /*
	auto init_function = [lbound, ubound] (ParticleSPSO& p) -> void { 
	popot::initializer::position::uniform_random<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), lbound, ubound);
	popot::initializer::velocity::half_diff<ParticleSPSO::VECTOR_TYPE>(p.getPosition(), p.getVelocity(), lbound, ubound);
	};
      */

      // The rule to update the best position
      auto best_position_update = popot::PSO::particle::updateBestPosition<ParticleSPSO>;

      // Topology
      auto topology = popot::PSO::topology::full_fillNeighborhoods<ParticleSPSO>;

      auto algo = popot::PSO::algorithm::base(swarm_size, dimension, 
					      lbound, ubound, stop, cost_function, 
					      topology, position_update, velocity_update, best_position_update, 
					      p, popot::PSO::algorithm::ASYNCHRONOUS_EVALUATION, false);
      return algo;
    }

} // namespace algorithm

  } // namespace popot

#endif // POPOT_ALGORITHM_H
