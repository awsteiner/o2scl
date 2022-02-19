/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_RNG_GSL_H
#define O2SCL_RNG_GSL_H

#include <iostream>

// For memcpy()
#include <cstring>

#ifdef O2SCL_NEEDS_TIME_H
#include "time.h"
#endif

#include <gsl/gsl_rng.h>

namespace o2scl {

  /** \brief Random number generator (GSL)

      If \c seed is zero, or is not given, then the default seed
      specific to the particular random number generator is used. 
      
      \future This is being replaced by \ref o2scl::rng now that
      C++11 support is more widespread.
  */
  class rng_gsl : public gsl_rng {

  public:

    /// Desc
    typedef unsigned long int result_type;

    /** \brief Initialize the random number generator with type \c gtype 
	and the default seed 
    */
    rng_gsl(const gsl_rng_type *gtype=gsl_rng_mt19937);
    
    /// Initialize the random number generator with \c seed 
    rng_gsl(unsigned long int seed, 
	     const gsl_rng_type *gtype=gsl_rng_mt19937);
    
    ~rng_gsl();
    
    /// Return generator type
    const gsl_rng_type *get_type() { return rng; }

    /** \brief Return a random number in \f$(0,1]\f$
     */
    result_type operator()() {
      return gsl_rng_get(this);
    }
    
    /// Return a random number in \f$(0,1]\f$
    double random() {
      return (this->type->get_double)(this->state);
    }

    /** \brief Return the entropy (0.0 since not applicable for
	pseudo-random engines */
    double entropy() {
      return 0.0;
    }

    /// Return the maximum integer for random_int()
    unsigned long int max() {
      return gsl_rng_max(this);
    }

    /// Return the minimum integer for random_int()
    unsigned long int min() {
      return gsl_rng_min(this);
    }

    /// Return random integer in \f$[0,\mathrm{max}-1]\f$.
    unsigned long int random_int(unsigned long int max=0);
    
    /// Set the seed
    void set_seed(unsigned long int s) { 
      seed=s;
      gsl_rng_set(this,seed);
    }

    /// Get the seed
    unsigned long int get_seed() {
      return seed;
    }
    
    /// Set the seed
    void clock_seed() {
      seed=time(0);
      gsl_rng_set(this,seed);
    }
    
    /// Copy constructor with equals operator
    rng_gsl& operator=(const rng_gsl &rg) {
      if (this!=&rg) {
	seed=rg.seed;
	rng=rg.rng;
	this->state=malloc(rg.type->size);
	this->type=rg.type;
	memcpy(this->state,rg.state,this->type->size);
      }
      return *this;
    }

    /// Copy constructor
    rng_gsl(const rng_gsl &rg) {
      seed=rg.seed;
      rng=rg.rng;
      this->state=malloc(rg.type->size);
      this->type=rg.type;
      memcpy(this->state,rg.state,this->type->size);
    }

#ifndef DOXYGEN_INTERNAL
    
  protected:

    /// The seed
    unsigned long int seed;

    /// The GSL random number generator type
    const gsl_rng_type *rng;

#endif
  
  };

  /** \brief An alternative RNG type used to design templates 
      for either GSL or STL random number generators
   */
  class rng_gsl_uniform_real {
  public:
    double operator()(rng_gsl &r) {
      return (r.type->get_double)(r.state);
    }
  };

}

#endif
