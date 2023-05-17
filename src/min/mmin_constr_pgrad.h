/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
/*───────────────────────────────────────────────────────────────────-------*
 * Open Optimization Library - Projected Gradient Method
 * 
 * pgrad/pgrad.c
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * Iara da Cunha
 * since: June, 29, 2004 
 *
 * $Id: pgrad.c,v 1.14 2005/05/19 19:37:07 biloti Exp $
 *───────────────────────────────────────────────────────────────────-------*/
#ifndef O2SCL_OOL_MMIN_PGRAD_H
#define O2SCL_OOL_MMIN_PGRAD_H

/** \file mmin_constr_pgrad.h
    \brief File defining \ref o2scl::mmin_constr_pgrad 
*/

#include <gsl/gsl_math.h>

#include <o2scl/multi_funct.h>
#include <o2scl/mmin_constr.h>

namespace o2scl {

  /** \brief Constrained minimization by the projected gradient method (OOL)

      This is the simple extension of the steepest descent algorithm
      to constrained minimization. Each step is a line search is
      performed along the projected gradient direction subject to 
      the specified constraints. 

      \verbatim embed:rst
      This algorithm is likely not ideal for most problems and is
      provided mostly for demonstration and educational purposes.
      Based on implementation of [Kelley99]_ in OOL. 
      \endverbatim

      Default template arguments
      - \c func_t - \ref multi_funct
      - \c dfunc_t - \ref grad_funct
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>

      \future Replace the explicit norm computation below with 
      the more accurate dnrm2 from linalg
      \future Replace the generic variable 'tol' with 'tolf'
      or 'tolx' from \ref o2scl::mmin_base.
  */
  template<class func_t=multi_funct, class dfunc_t=grad_funct,
    class vec_t=boost::numeric::ublas::vector<double> > 
    class mmin_constr_pgrad : 
    public mmin_constr<func_t,dfunc_t,ool_hfunct,vec_t> {
    
#ifndef DOXYGEN_INTERNAL

    protected:

    /// A convenient typedef for the unused Hessian product type
    typedef ool_hfunct hfunc_t;

    /// Temporary vector
    vec_t xx;
    
    /// Project into feasible region
    int proj(vec_t &xt) {
      size_t ii, n=this->dim;
      
      for(ii=0;ii<n;ii++) {
	double dtmp;
	if (xt[ii]<this->U[ii]) dtmp=xt[ii];
	else dtmp=this->U[ii];
	if (this->L[ii]>dtmp) xt[ii]=this->L[ii];
	else xt[ii]=dtmp;
	//xt[ii]=GSL_MAX(this->L[ii],GSL_MIN(xt[ii],this->U[ii]));
      }
      return 0;
    }

    /// Line search
    int line_search() {

      double t, tnew, fx, fxx, dif2, gtd;

      fx=this->f;

      tnew=1;
      
      do {
	
	t = tnew;
	/* xx = x - t df */
	o2scl::vector_copy(this->dim,this->x,xx);
	for(size_t i=0;i<this->dim;i++) xx[i]-=t*this->gradient[i];
	proj(xx);
	
	fxx=(*this->func)(this->dim,xx);
	
	o2scl::vector_copy(this->dim,xx,this->dx);
	for(size_t i=0;i<this->dim;i++) this->dx[i]-=this->x[i];

	dif2=0.0;
	for(size_t i=0;i<this->dim;i++) dif2+=this->dx[i]*this->dx[i];
	
	gtd=0.0;
	for(size_t i=0;i<this->dim;i++) {
	  gtd+=this->gradient[i]*this->dx[i];
	}
	
	tnew = - t * t * gtd / (2*(fxx -  this->f - gtd));

	double arg1=sigma1*t;
	double arg2;
	if (sigma2*t<tnew) arg2=sigma2*t;
	else arg2=tnew;
	if (arg1>arg2) tnew=arg1;
	else tnew=arg2;
	
	// sufficient decrease condition (Armijo)
      } while(fxx > fx - (alpha/t) * dif2 );
      
      o2scl::vector_copy(this->dim,xx,this->x);
      this->f = fxx;
      
      (*this->dfunc)(this->dim,this->x,this->gradient);
      
      return 0;
    }
    
#endif

    public:

    mmin_constr_pgrad() {
      fmin=-1.0e99;
      tol=1.0e-4;
      alpha=1.0e-4;
      sigma1=0.1;
      sigma2=0.9;
      this->ntrial=1000;
    }

    /** \brief Minimum function value (default \f$ 10^{-99} \f$)
	
	If the function value is below this value, then the algorithm
	assumes that the function is not bounded and exits.
    */
    double fmin;

    /// Tolerance on infinite norm
    double tol;

    /** \brief Constant for the sufficient decrease condition 
	(default \f$ 10^{-4} \f$)
    */
    double alpha;

    /// Lower bound to the step length reduction
    double sigma1;

    /// Upper bound to the step length reduction
    double sigma2;
    
    /// Allocate memory
    virtual int allocate(const size_t n) {
      if (this->dim>0) free();
      xx.resize(n);
      return mmin_constr<func_t,dfunc_t,hfunc_t,vec_t>::allocate(n);
    }
    
    /// Free previously allocated memory
    virtual int free() {
      xx.clear();
      return 0;
    }
    
    /// Set the function, the initial guess, and the parameters
    virtual int set(func_t &fn, dfunc_t &dfn, vec_t &init) {
      
      int ret=mmin_constr<func_t,dfunc_t,hfunc_t,vec_t>::set(fn,dfn,init);
      
      // Turn initial guess into a feasible point
      proj(this->x);

      // Evaluate function and gradient
      this->f=(*this->func)(this->dim,this->x);
      (*this->dfunc)(this->dim,this->x,this->gradient);
      
      return 0;
    }

    /// Restart the minimizer
    virtual int restart() {
      // Turn x into a feasible point
      proj(this->x);

      // Evaluate function and gradient
      this->f=(*this->func)(this->dim,this->x);
      (*this->dfunc)(this->dim,this->x,this->gradient);

      return 0;
    }

    /// Perform an iteration
    virtual int iterate() {
      double t;

      line_search();
      
      /* Infinite norm of g1 = d <- [P(x - g) - x] */
      o2scl::vector_copy(this->dim,this->x,xx);
      for(size_t i=0;i<this->dim;i++) {
	xx[i]-=this->gradient[i];
      }
      proj(xx);
      for(size_t i=0;i<this->dim;i++) xx[i]-=this->x[i];
      
      double maxval=fabs(xx[0]);
      for(size_t i=1;i<this->dim;i++) {
	if (fabs(xx[i])>maxval) {
	  maxval=fabs(xx[i]);
	}
      }
      this->size=fabs(maxval);

      return success;
    }

    /// See if we're finished
    virtual int is_optimal() {
      if (this->size>tol && this->f>fmin) {
	return gsl_continue;
      }
      return success;
    }

    /// Return string denoting type ("mmin_constr_pgrad")
    const char *type() { return "mmin_constr_pgrad"; }

#ifndef DOXYGEN_INTERNAL

  private:
  
  mmin_constr_pgrad<func_t,dfunc_t,vec_t>
  (const mmin_constr_pgrad<func_t,dfunc_t,vec_t> &);
  mmin_constr_pgrad<func_t,dfunc_t,vec_t>& operator=
  (const mmin_constr_pgrad<func_t,dfunc_t,vec_t>&);

#endif
      
  };
  
}

#endif

