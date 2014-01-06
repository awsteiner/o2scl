/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

#ifdef O2SCL_NEEDS_TIME_H
#include "time.h"
#endif

#include <gsl/gsl_rng.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Random number generator (GSL)

      If \c seed is zero, or is not given, then the default seed
      specific to the particular random number generator is used. 

      An interesting application of this class to generate an
      arbitrary distribution through a Markov chain Monte Carlo
      method is in \ref ex_markov_sect.
  */

  class rng_gsl {

  public:

    /** \brief Initialize the random number generator with type \c gtype 
	and the default seed 
    */
    rng_gsl(const gsl_rng_type *gtype=gsl_rng_mt19937);

    rng_gsl(double ig1, double ig2) {
      rng=gsl_rng_mt19937;
      gr=gsl_rng_alloc(gsl_rng_mt19937);
      gsl_rng_set(gr,0);
      seed=time(0);
    }
    
    /// Initialize the random number generator with \c seed 
    rng_gsl(unsigned long int seed, 
	     const gsl_rng_type *gtype=gsl_rng_mt19937);
    
    ~rng_gsl();
    
    /// Return generator type
    const gsl_rng_type *get_type() { return rng; }

    /** \brief Return a random number in \f$(0,1]\f$
     */
    double operator()(int ignored) {
      return random();
    }
    
    /// Return a random number in \f$(0,1]\f$
    double random() {
      return (gr->type->get_double)(gr->state);
    }

    /// Return the maximum integer for random_int()
    unsigned long int get_max() {
      return gsl_rng_max(gr)-gsl_rng_min(gr);
    }

    /// Return random integer in \f$[0,\mathrm{max}-1]\f$.
    unsigned long int random_int(unsigned long int n=0);
    
    /// Set the seed
    void set_seed(unsigned long int s) { 
      seed=s;
      gsl_rng_set(gr,seed);
    }

    /// Set the seed
    void clock_seed() {
      seed=time(0);
      gsl_rng_set(gr,seed);
    }

  protected:

#ifndef DOXYGEN_NO_O2NS

    rng_gsl(const rng_gsl &);
    rng_gsl& operator=(const rng_gsl&);

#endif
#ifndef DOXYGEN_INTERNAL

    /// The seed
    unsigned long int seed;

    /// The GSL random number generator
    gsl_rng *gr;

    /// The GSL random number generator type
    const gsl_rng_type *rng;

#endif
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
