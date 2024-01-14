/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
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
#include <o2scl/rng.h>

namespace o2scl {

  /** \brief Monte-Carlo integration [abstract base]
      
      This class provides the generic Monte Carlo parameters and the
      random number generator. The default type for the random number
      generator is a \ref rng object. 
  */
  template<class func_t=multi_funct, 
    class vec_t=boost::numeric::ublas::vector<double>,
           class rng_t=rng<> >
    class mcarlo : public inte_multi<func_t,vec_t> {

  public:
  
  mcarlo() {
      n_points=1000;
    }

    virtual ~mcarlo() {}
  
    /** \brief Number of integration points (default 1000)
     */
    unsigned long n_points;
  
    /// The random number generator
    rng_t rng;
  
    /// Return string denoting type ("mcarlo")
    virtual const char *type() { return "mcarlo"; }
  
  };

}

#endif

