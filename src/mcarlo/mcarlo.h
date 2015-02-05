/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Andrew W. Steiner
  
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
#ifndef O2SCL_MCARLO_H
#define O2SCL_MCARLO_H

/** \file mcarlo.h
    \brief File defining \ref o2scl::mcarlo
*/

#include <iostream>
#include <random>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/inte_multi.h>
#include <o2scl/rng_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Monte-Carlo integration [abstract base]
      
      This class provides the generic Monte Carlo parameters and the
      random number generator. The default type for the random number
      generator is a \ref rng_gsl object. 
  */
  template<class func_t=multi_funct11, 
    class vec_t=boost::numeric::ublas::vector<double>,
    class rng_t=int, class rng_dist_t=rng_gsl>
    //class rng_t=std::mt19937, 
    //    class rng_dist_t=std::uniform_real_distribution<double> > 
    class mcarlo : public inte_multi<func_t,vec_t> {

  public:
  
  mcarlo() : rng_dist(0.0,1.0) {
      n_points=1000;
    }

    virtual ~mcarlo() {}
  
    /** \brief Number of integration points (default 1000)
     */
    unsigned long n_points;

    /// The random number distribution
    rng_dist_t rng_dist;
      
    /// The random number generator
    rng_t rng;
  
    /// Return string denoting type ("mcarlo")
    virtual const char *type() { return "mcarlo"; }
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

