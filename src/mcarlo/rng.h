/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
/** \file rng_gsl.h
    \brief File for definition of \ref o2scl::rng_gsl
*/
#ifndef O2SCL_RNG_H
#define O2SCL_RNG_H

#include <random>
#include <ctime>
#include <o2scl/err_hnd.h>

namespace o2scl {

  /** \brief Simple C++11 random number generator class
   */
  template<class fp_t=double> class rng {
    
  protected:

    /// Random number engine
    std::default_random_engine def_engine;

    /// Distribution for \f$ [0,1) \f$
    std::uniform_real_distribution<fp_t> dist;

    /// Seed for the random number generator
    unsigned int seed;

  public:
    
    rng() {
      seed=time(0);
      def_engine.seed(seed);
    }

    /// Set the seed using <tt>time(0)</tt>
    void clock_seed() {
      seed=time(0);
      def_engine.seed(seed);
      return;
    }

    /// Get the seed
    unsigned int get_seed() {
      return seed;
    }
    
    /// Return a random number in \f$(0,1]\f$
    fp_t random() {
      return dist(def_engine);
    }
  
    /// Set the seed
    void set_seed(unsigned int s) { 
      seed=s;
      def_engine.seed(seed);
      return;
    }

    /// Return random integer in \f$[0,\mathrm{max}-1]\f$.
    unsigned long int random_int(unsigned long int max=1) {
      if (max<1) {
        O2SCL_ERR("Cannot specify max less than 1.",o2scl::exc_einval);
      }
      std::uniform_int_distribution<unsigned long int> uid(0,max-1);
      return uid(def_engine);
    }
    
  };
  
};

#endif
