/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Random number generator (GSL)

      This object is built on the <tt>gsl_rng</tt> struct and modeled
      to look like a <tt>std::random_device</tt> object.

      If \c seed is zero, or is not given, then the default seed
      specific to the particular random number generator is used. 
      
      \future This will likely be completely replaced by the random
      number generators in the standard library.
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

  /** \brief Swap function for vector_shuffle()
      
      \note This function is based on the static GSL swap function in
      <tt>randist/shuffle.c</tt>.
   */
  template<class data_t>
  void shuffle_swap(data_t *base, size_t size, size_t i, size_t j) {
    register char * a = size * i + (char *) base ;
    register char * b = size * j + (char *) base ;
    register size_t s = size ;
    
    if (i == j) {
      return;
    }
    
    do {                                       
      char tmp = *a;                            
      *a++ = *b;                                
      *b++ = tmp;                               
    } while (--s > 0);
    
    return;
  }

  /** \brief Copy function for vector_choose()

      \note This function is based on the static GSL copy function in
      <tt>randist/shuffle.c</tt>.
   */
  template<class data_t>
  void choose_copy(data_t *dest, size_t i, data_t *src, size_t j,
		    size_t size) {
    register char * a = size * i + (char *) dest;
    register char * b = size * j + (char *) src;
    register size_t s = size ;
    
    do {                                          
      *a++ = *b++;                              
    } 
    while (--s > 0);

    return;
  }

  /** \brief Shuffle the first \c n elements of vector \c data

      \note This function is based on the GSL function
      <tt>gsl_ran_shuffle()</tt> in <tt>randist/shuffle.c</tt>.

      \note This function works only on vector types which
      guarantee adjacent storage, such as <tt>vector<double></tt>.
      
      \note If \c n is 0, this function silently does nothing.
   */
  template<class vec_t, class data_t>
  void vector_shuffle(rng_gsl &r, size_t n, vec_t &data) {
    if (n==0) return;
    
    for (size_t i = n - 1; i > 0; i--) {
      size_t j = r.random_int(i+1);
      shuffle_swap(&data[0],sizeof(data_t),i,j);
    }
    
    return;
  }

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
