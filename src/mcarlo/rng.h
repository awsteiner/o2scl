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
/** \file rng.h
    \brief File for definition of \ref o2scl::rng
*/
#ifndef O2SCL_RNG_H
#define O2SCL_RNG_H

#include <boost/random.hpp>
#include <random>
#include <ctime>
#include <o2scl/err_hnd.h>

namespace o2scl {

  /** \brief Simple C++11 random number generator class
   */
  template<class fp_t=double> class rng {
    
  protected:

    /// Distribution for \f$ [0,1) \f$
    mutable std::uniform_real_distribution<fp_t> dist;

    /// Seed for the random number generator
    unsigned int seed;

  public:
    
    /// Random number engine
    mutable std::mt19937 def_engine;

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
    unsigned int get_seed() const {
      return seed;
    }
    
    /// Return a random number in \f$(0,1]\f$
    fp_t random() const {
      return dist(def_engine);
    }
  
    /// Set the seed
    void set_seed(unsigned int s) { 
      seed=s;
      def_engine.seed(seed);
      return;
    }

    /// Return random integer in \f$[0,\mathrm{max}-1]\f$.
    unsigned long int random_int(unsigned long int max=1) const {
      if (max<1) {
        O2SCL_ERR("Cannot specify max less than 1.",o2scl::exc_einval);
      }
      std::uniform_int_distribution<unsigned long int> uid(0,max-1);
      return uid(def_engine);
    }

    /// Copy constructor with equals operator
    rng& operator=(const rng &rg) {
      if (this!=&rg) {
	seed=rg.seed;
	dist=rg.dist;
	def_engine=rg.def_engine;
      }
      return *this;
    }

    /// Copy constructor
    rng(const rng &rg) {
      seed=rg.seed;
      dist=rg.dist;
      def_engine=rg.def_engine;
    }

    
  };

  /** \brief Shuffle the first \c n elements of vector \c data

      \note If \c n is 0, this function silently does nothing.
  */
  template<class vec_t, class data_t>
  void vector_shuffle(rng<> &r, size_t n, vec_t &data) {
    if (n==0) return;
    
    for (int i=0;i<((int)n)-1;i++) {
      int j=r.random_int(n-i)+i;
      if (i!=j) {
        data_t tmp=data[i];
        data[i]=data[j];
        data[j]=tmp;
      }
    }
    
    return;
  }

};

#endif
