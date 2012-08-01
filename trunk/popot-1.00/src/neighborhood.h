#ifndef POPOT_NEIGHBORHOOD_H
#define POPOT_NEIGHBORHOOD_H

#include <vector>
#include "exceptions.h"

namespace popot
{
  namespace PSO
  {

    namespace particle
    {
      template<typename PROBLEM,typename POSITION_INITIALIZER> class BaseParticle;
    }

    namespace neighborhood
    {

      /**
       * Neighborhood
       * @param PARTICLE The particle type you want to use, e.g. swarm::particle::TraditionalParticle
       * \brief A neighborhood contains a vector of particles belonging to the neighborhood as well as a pointer (not a copy!)
       *        to the best particle
       */
      template<typename PARTICLE>
        class Neighborhood
        {
	public:
	  typedef typename PARTICLE::BestType BestType;
	  typedef PARTICLE InNeighborhoodType;

        protected:
	  std::vector<InNeighborhoodType *> _particles;
	  BestType *_best_particle;

        public:
	  Neighborhood(void){};
	  virtual ~Neighborhood(void){
	    clear();
	  };

	  void add(InNeighborhoodType * p)
	  {
	    _particles.push_back(p);
	  }

	  int size() const
	  {
	    return _particles.size();
	  }

	  void clear()
	  {
	    _particles.clear();
	  }

	  InNeighborhoodType * get(int i)
	  {
	    if(i >= 0 && i < size())
	      return _particles[i];
	    else
	      throw popot::Exception::IndexOutOfRange(i, size());
	  }

	  BestType * findBest(void)
	  {
	    if(_particles.size() == 0)
	      throw popot::Exception::FindBestFromEmptyNeighborhood();

	    _best_particle = _particles[0]->getBestPosition();
	    for(unsigned int i = 1 ; i < _particles.size() ; ++i)
	      {
		if(_particles[i]->getBestPosition()->compare(_best_particle) < 0)
		  _best_particle = _particles[i]->getBestPosition();
	      }
	    return _best_particle;
	  }

	  void updateBest(InNeighborhoodType * p)
	  {
	    if(p->getBestPosition()->compare(_best_particle) < 0)
	      _best_particle = p->getBestPosition();
	  }

	  BestType * getBest(void)
	  {
	    return _best_particle;
	  }
        };
    } // namespace neighborhood
  } // namespace PSO
} // namespace popot

#endif // POPOT_NEIGHBORHOOD_H
