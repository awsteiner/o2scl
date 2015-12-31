/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
/* monte/plain.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
#ifndef O2SCL_MCARLO_PLAIN_H
#define O2SCL_MCARLO_PLAIN_H

/** \file mcarlo_plain.h
    \brief File defining \ref o2scl::mcarlo_plain
*/
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/multi_funct.h>
#include <o2scl/mcarlo.h>
#include <o2scl/rng_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Multidimensional integration using plain Monte Carlo (GSL)
   */
  template<class func_t=multi_funct11, 
    class vec_t=boost::numeric::ublas::vector<double>,
    class rng_t=int, class rng_dist_t=rng_gsl>
    //class rng_t=std::mt19937, 
    //    class rng_dist_t=std::uniform_real_distribution<double> > 
    class mcarlo_plain : public mcarlo<func_t,vec_t,rng_t,rng_dist_t> {
    
  public:
  
  virtual ~mcarlo_plain() {}
  
  /// Integrate function \c func from x=a to x=b.
  virtual int minteg_err(func_t &func, size_t ndim, const vec_t &a, 
			 const vec_t &b, double &res, double &err) {

    double r;
      
    vec_t cc(ndim);
      
    double vol=1.0, m=0.0, q=0.0;
    for(size_t i=0;i<ndim;i++) vol*=b[i]-a[i];
      
    for(size_t n=0;n<this->n_points;n++) {
      for(size_t i=0;i<ndim;i++) {
	do {
	  r=this->rng_dist(this->rng);
	} while (r==0.0);
	cc[i]=a[i]+r*(b[i]-a[i]);
      }
      double fval;
      fval=func(ndim,cc);
      double d=fval-m;
      m+=d/(n+1.0);
      q+=d*d*(n/(n+1.0));
    }

    res=vol*m;
      
    if (this->n_points<2) {
      err=0.0;
    } else {
      err=vol*sqrt(q/(this->n_points*(this->n_points-1.0)));
    }
	  
    return 0;

  }
      
  /** \brief Integrate function \c func over the hypercube from
      \f$ x_i=a_i \f$ to \f$ x_i=b_i \f$ for \f$ 0<i< \f$ ndim-1
  */
  virtual double minteg(func_t &func, size_t ndim, const vec_t &a, 
			const vec_t &b) {
    double res;
    minteg_err(func,ndim,a,b,res,this->interror);
    return res;
  }

  /// Return string denoting type ("mcarlo_plain")
  virtual const char *type() { return "mcarlo_plain"; }
      
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

