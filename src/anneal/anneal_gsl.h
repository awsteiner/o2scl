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
/* siman/siman.c
 * 
 * Copyright (C) 2007 Brian Gough
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Mark Galassi
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
#ifndef O2SCL_ANNEAL_GSL_H
#define O2SCL_ANNEAL_GSL_H

/** \file anneal_gsl.h
    \brief File defining \ref o2scl::anneal_gsl
*/
#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/anneal.h>
#include <o2scl/multi_funct.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Multidimensional minimization by simulated annealing (GSL)
      
      This class is a modification of simulated annealing as
      implemented in GSL in the function \c gsl_siman_solve(). It acts
      as a generic multidimensional minimizer for any function given a
      generic temperature schedule specified by the user.

      There are a large variety of strategies for choosing the
      temperature evolution. To offer the user the largest possible
      flexibility, the temperature evolution is controlled by a 
      the virtual functions start() and next() which can be freely
      changed by creating a child class which overwrites these 
      functions. 

      The simulated annealing algorithm proposes a displacement of one 
      coordinate of the previous point by
      \f[
      x_{i,\mathrm{new}} = \mathrm{step\_size}_i 
      (2 u_i - 1) + x_{i,\mathrm{old}}
      \f]
      where the \f$u_i\f$ are random numbers between 0 and 1. The
      displacement is accepted or rejected based on the Metropolis
      method. The random number generator is set in the parent,
      anneal.

      The default behavior is as follows: Initially, the step sizes
      are chosen to be 1.0 (or whatever was recently specified in \ref
      set_step() ) and the temperature to be \ref T_start (default
      1.0). Each iteration decreases the temperature by a factor of
      \ref T_dec (default 1.5) for each step, and the minimizer is
      finished when the next decrease would bring the temperature
      below \ref o2scl::mmin_base::tol_abs. If none of the
      mmin_base::ntrial steps in a particular iteration changes the
      value of the minimum, and the step sizes are greater than \ref
      min_step_ratio (default 100) times \ref
      o2scl::mmin_base::tol_abs, then the step sizes are decreased by
      a factor of \ref step_dec (default 1.5) for the next iteration.

      If \ref o2scl::mmin_base::verbose is greater than zero, then
      \ref mmin() will print out information and/or request a keypress
      after the function iterations for each temperature.
      
      \verbatim embed:rst
      An example demonstrating the usage of this class is given in
      ``examples/ex_anneal.cpp`` and in the 
      :ref:`Simulated annealing example`.
      \endverbatim

      \comment
      The form of the user-specified function is as in \ref
      multi_funct has a "function value" which is the value of the
      function (given in the third argument as a number of type \c
      double), and a "return value" (the integer return value). The
      initial function evaluation which is performed at the
      user-specified initial guess must give 0 as the return value. If
      subsequent function evaluations have a non-zero return value,
      then the resulting point is ignored and a new point is selected.

      This class thus can sometimes handle constrained minimization
      problems. If the user ensures that the function's return value
      is non-zero when the function is evaluated outside the allowed
      region, the minimizer will not accept any steps which take the
      minimizer outside the allowed region. Note that this should be
      done with care, however, as this approach may cause convergence
      problems with sufficiently difficult functions or constraints.
      \endcomment

      \comment
      \future There's x0, old_x, new_x, and x? There's probably
      some duplication here which could be avoided.
      9/19/14: Some better documentation now, and it looks like
      all four have some utility. 
      \endcomment

      \future Implement a more general simulated annealing routine
      which would allow the solution of discrete problems like the
      Traveling Salesman problem.
      \comment
      This could probably be done just by creating a parent abstract 
      base class which doesn't have any reference to step() and
      step_vec.
      \endcomment
      
  */
  template<class func_t=multi_funct,
    class vec_t=boost::numeric::ublas::vector<double>,
           class rng_t=o2scl::rng<> > class anneal_gsl :
    public anneal_base<func_t,vec_t,rng_t> {
    
  public:
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  
  anneal_gsl() {
    boltz=1.0;
    step_vec.resize(1);
    step_vec[0]=1.0;
    T_start=1.0;
    T_dec=1.5;
    step_dec=1.5;
    min_step_ratio=100.0;
  }

  virtual ~anneal_gsl() {
  }
  
  /** \brief Calculate the minimum \c fmin of \c func w.r.t the 
      array \c x0 of size \c nvar.
  */
  virtual int mmin(size_t nvar, vec_t &x0, double &fmin, 
		   func_t &func) {
    
    if (nvar==0) {
      O2SCL_ERR2("Tried to minimize over zero variables ",
		 " in anneal_gsl::mmin().",exc_einval);
    }
    
    fmin=0.0;

    step_norm=1.0;
    
    double E, new_E, T, old_E;
    int i, iter=0;
    size_t j;
    
    // Current point
    vec_t x(nvar);
    // Proposed next point
    vec_t new_x(nvar);
    // Last point from previous iteration
    vec_t old_x(nvar);

    for(j=0;j<nvar;j++) {
      x[j]=x0[j];
    }
    
    E=func(nvar,x);

    // We use x0 and fmin to store the best point found 
    fmin=E;

    // Setup initial temperature and step sizes
    start(nvar,T);

    bool done=false;

    while (!done) {

      // Copy old value of x for next() function
      for(j=0;j<nvar;j++) old_x[j]=x[j];
      old_E=E;
	  
      size_t nmoves=0;

      for (i=0;i<this->ntrial;++i) {
	for (j=0;j<nvar;j++) new_x[j]=x[j];
      
	step(new_x,nvar);
	new_E=func(nvar,new_x);
	
	// Store best value obtained so far
	if(new_E<=fmin) {
	  for(j=0;j<nvar;j++) x0[j]=new_x[j];
	  fmin=new_E;
	}
	
	// Take the crucial step: see if the new point is accepted
	// or not, as determined by the Boltzmann probability
	if (new_E<E) {
	  for(j=0;j<nvar;j++) x[j]=new_x[j];
	  E=new_E;
	  nmoves++;
	} else {
	  double r=this->local_rng.random();
          if (this->verbose>=3) {
            std::cout << x[0] << " " << new_x[0] << " "
                      << r << " " << exp(-(new_E-E)/(boltz*T)) << " "
                      << (r < exp(-(new_E-E)/(boltz*T))) << " "
                      << nmoves << std::endl;
          }
          if (r < exp(-(new_E-E)/(boltz*T))) {
	    for(j=0;j<nvar;j++) x[j]=new_x[j];
	    E=new_E;
	    nmoves++;
	  }
	}

      }
	  
      if (this->verbose>0) {
	this->print_iter(nvar,x0,fmin,iter,T,"anneal_gsl");
	iter++;
      }
      
      // See if we're finished and proceed to the next step
      next(nvar,old_x,old_E,x,E,T,nmoves,x0,fmin,done);
      
    }
  
    return 0;
  }
      
  /// Return string denoting type ("anneal_gsl")
  virtual const char *type() { return "anneal_gsl"; }

  /** \brief Boltzmann factor (default 1.0).
   */
  double boltz;
      
  /// Set the step sizes 
  template<class vec2_t> int set_step(size_t nv, vec2_t &stepv) {
    if (nv>0) {
      step_vec.resize(nv);
      for(size_t i=0;i<nv;i++) step_vec[i]=stepv[i];
    }
    return 0;
  }

  /// Initial temperature (default 1.0)
  double T_start;

  /// Factor to decrease temperature by (default 1.5)
  double T_dec;

  /// Factor to decrease step size by (default 1.5)
  double step_dec;

  /// Ratio between minimum step size and \ref tol_abs (default 100.0)
  double min_step_ratio;

  /** \brief Copy constructor
   */
  anneal_gsl<func_t,vec_t,rng_t>
  (const anneal_gsl<func_t,vec_t,rng_t> &ag) :
  anneal_base<func_t,vec_t,rng_t>() {
    
    boltz=ag.boltz;
    T_start=ag.T_start;
    T_dec=ag.T_dec;
    step_dec=ag.step_dec;
    min_step_ratio=ag.min_step_ratio;
    step_vec=ag.step_vec;
  }
  
  /** \brief Copy constructor from operator=
   */
  anneal_gsl<func_t,vec_t,rng_t>& operator=
  (const anneal_gsl<func_t,vec_t,rng_t> &ag) {
    if (this != &ag) {
      anneal_base<func_t,vec_t,rng_t>::operator=(ag);
      boltz=ag.boltz;
      T_start=ag.T_start;
      T_dec=ag.T_dec;
      step_dec=ag.step_dec;
      min_step_ratio=ag.min_step_ratio;
      step_vec=ag.step_vec;
    }
    return *this;
  }
      
#ifndef DOXYGEN_INTERNAL
      
  protected:
      
  /// \name Storage for points in parameter space
  //@{
  //@}

  /// Normalization for step
  double step_norm;
  
  /// Vector of step sizes
  ubvector step_vec;
      
  /// Determine how to change the minimization for the next iteration
  virtual int next(size_t nvar, vec_t &x_old, double min_old, 
		   vec_t &x_new, double min_new, double &T, 
		   size_t n_moves, vec_t &best_x, double best_E, 
		   bool &finished) {
    
    if (T/T_dec<this->tol_abs) {
      finished=true;
      return success;
    }

    if (n_moves==0) {
      // If we haven't made any progress, shrink the step by
      // decreasing step_norm
      step_norm/=step_dec;
      if (this->verbose>=3) {
        std::cout << "Shrinking step " << step_norm << " " 
                  << step_dec << std::endl;
      }
      // Also reset x to best value so far
      for(size_t i=0;i<nvar;i++) {
	x_new[i]=best_x[i];
      }
      min_new=best_E;
    }
    T/=T_dec;
    return success;
  }

  /// Setup initial temperature and stepsize
  virtual int start(size_t nvar, double &T) {
    T=T_start;
    return success;
  }

  /** \brief Make a step to a new attempted minimum
   */
  virtual int step(vec_t &sx, int nvar) {
    size_t nstep=step_vec.size();
    for(int i=0;i<nvar;i++) {
      double u=this->local_rng.random();

      // Construct the step in the ith direction
      double step_i=step_norm*step_vec[i%nstep];
      // Fix if the step is too small
      if (step_i<this->tol_abs*min_step_ratio) {
	step_i=this->tol_abs*min_step_ratio;
      }

      sx[i]=(2.0*u-1.0)*step_i+sx[i];
    }
    return 0;
  }
  
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
