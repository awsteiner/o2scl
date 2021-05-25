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
#ifndef O2SCL_ANNEAL_H
#define O2SCL_ANNEAL_H

/** \file anneal.h
    \brief File defining \ref o2scl::anneal_base
*/

#include <iostream>

#include <boost/config.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/multi_funct.h>
#include <o2scl/mmin.h>
#include <o2scl/rng.h>
#include <o2scl/prob_dens_func.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Simulated annealing base

      The seed of the generator is not fixed initially by calls to
      mmin(), so if successive calls should reproduce the same
      results, then the random seed should be set by the user before
      each call.
      
      For the algorithms here, it is important that all of the inputs
      <tt>x[i]</tt> to the function are scaled similarly relative to
      the temperature. For example, if the inputs <tt>x[i]</tt> are
      all of order 1, one might consider a temperature schedule which
      begins with \f$ T=1 \f$ .

      The number of iterations at each temperature is controlled by
      \ref o2scl::mmin_base::ntrial which defaults to 100.

      \verbatim embed:rst
      .. todo:: 
         
         In class anneal_base: I'm having trouble with
         std::uniform_real_distribution on clang at the moment, so
         this class uses \ref o2scl::prob_dens_uniform for the moment.

      \endverbatim
  */
  template<class func_t=multi_funct,
    class vec_t=boost::numeric::ublas::vector<double>,
           class rng_t=rng<> > class anneal_base :
    public mmin_base<func_t,func_t,vec_t> {
      
#ifdef O2SCL_NEVER_DEFINED
    }
  {
#endif
    
  public:
    
  anneal_base() : dist(0.0,1.0) {
      this->ntrial=100;
    }
      
    virtual ~anneal_base() {}
      
    /** \brief Calculate the minimum \c fmin of \c func w.r.t the 
        array \c x of size \c nvar.
    */
    virtual int mmin(size_t nvar, vec_t &x, double &fmin, 
                     func_t &func)=0;
    
    /** \brief Print out iteration information.

        Depending on the value of the variable verbose, this prints out
        the iteration information. If verbose=0, then no information is
        printed, while if verbose>1, then after each iteration, the
        present values of x and y are output to std::cout along with the
        iteration number. If verbose>=2 then each iteration waits for a
        character.  
    */
    virtual int print_iter(size_t nv, vec_t &x, double y, int iter,
                           double tptr, std::string comment) 
    {
      if (this->verbose<=0) return 0;
        
      size_t i;
      char ch;
      
      (*this->outs) << comment << " Iteration: " << iter << std::endl;
      std::cout << "x: ";
      for(i=0;i<nv;i++) std::cout << x[i] << " ";
      std::cout << std::endl;
      (*this->outs) << "y: " << y << " Tptr: " << tptr << std::endl;
      if (this->verbose>1) {
        (*this->outs) << "Press a key and type enter to continue. ";
        (*this->ins) >> ch;
      }
        
      return 0;
    }

    /// The default random number generator
    rng_t rng;

    /// The random distribution object
    o2scl::prob_dens_uniform dist;
    //std::uniform_real_distribution<> dist;

    /// Return string denoting type, \c "anneal_base".
    virtual const char *type() { return "anneal_base"; }

    /** \brief Copy constructor
     */
    anneal_base<func_t,vec_t,rng_t>
      (const anneal_base<func_t,vec_t,rng_t> &ab) : 
    mmin_base<func_t,func_t,vec_t>() {
      
      this->rng=ab.rng;
      
    }
    
    /** \brief Copy constructor from operator=
     */
    anneal_base<func_t,vec_t,rng_t>& operator=
      (const anneal_base<func_t,vec_t,rng_t> &ab) {

      if (this != &ab) {
        mmin_base<func_t,func_t,vec_t>::operator=(ab);
        this->rng=ab.rng;
      }
      return *this;
    }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

