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
/*───────────────────────────────────────────────────────────────────-----*
 * Open Optimization Library - Spectral Projected Gradient Method
 * 
 * spg/spg.c
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
 * Ricardo Biloti
 * since: May 3rd, 2005
 *
 * $Id: spg.c,v 1.4 2005/05/10 20:24:27 biloti Exp $
 *───────────────────────────────────────────────────────────────────-----*/
#ifndef O2SCL_OOL_MMIN_SPG_H
#define O2SCL_OOL_MMIN_SPG_H

/** \file mmin_constr_spg.h
    \brief File defining \ref o2scl::mmin_constr_spg 
*/

#include <o2scl/multi_funct.h>
#include <o2scl/mmin_constr.h>
#include <gsl/gsl_math.h>

namespace o2scl {

  /** \brief Constrained minimization by the spectral projected
      gradient method (OOL)

      This class applies a non-monotone line search strategy 
      to the classical projected gradient method. 
      
      \verbatim embed:rst
      As in [Birgin00]_, this class applies a nonmonotone Armijo
      sufficient decrease condition for accepting trial points as an
      improvement over the classical spectral projected gradient
      method. This method may be competitive with large problems
      because it has low memory requirements.
      \endverbatim

      Default template arguments
      - \c func_t - \ref multi_funct
      - \c dfunc_t - \ref grad_funct
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>

      \future There is some memory allocation which isn't
      deallocated until the destructor, and should be handled
      a bit more elegantly.
  */
  template<class func_t=multi_funct, class dfunc_t=grad_funct, 
    class vec_t=boost::numeric::ublas::vector<double> > class mmin_constr_spg : 
    public mmin_constr<func_t,dfunc_t,ool_hfunct,vec_t> {

#ifndef DOXYGEN_INTERNAL

    protected:

    /// A convenient typedef for the unused Hessian product type
    typedef ool_hfunct hfunc_t;

    /// Armijo parameter
    double alpha;

    /// Temporary vector
    vec_t xx;

    /// Temporary vector
    vec_t d;

    /// Temporary vector
    vec_t s;

    /// Temporary vector
    vec_t y;

    /// Temporary vector
    vec_t fvec;

    /// Non-monotone parameter
    size_t m;

    /// Desc
    int tail;
    
    /// Line search
    int line_search() {

      double lambda, lambdanew;
      double fmax, faux, fxx, dTg;
      short armijo_failure;
      size_t ii;

      // Saving the previous gradient and compute 
      // d = P(x - alpha * g) - x
      for(size_t i=0;i<this->dim;i++) {
	y[i]=this->gradient[i];
	d[i]=this->x[i];
      }
      for(size_t i=0;i<this->dim;i++) {
	d[i]-=alpha*this->gradient[i];
      }
      proj(d);
      for(size_t i=0;i<this->dim;i++) {
	d[i]-=this->x[i];
      }

      dTg=0.0;
      for(size_t i=0;i<this->dim;i++) dTg+=d[i]*this->gradient[i];

      lambda=1;
      armijo_failure=1;
      while (armijo_failure) {

	/* x trial */
	for(size_t i=0;i<this->dim;i++) xx[i]=this->x[i];
	for(size_t i=0;i<this->dim;i++) xx[i]+=lambda*d[i];
	fxx=(*this->func)(this->dim,xx);

	fmax=0;
	for(ii=0;ii<m;ii++) {
	  faux=fvec[ii]+gamma*lambda*dTg;
	  fmax=GSL_MAX(fmax,faux);
	}
	
	if (fxx>fmax) {
	  armijo_failure=1;
	  lambdanew=-lambda*lambda*dTg/(2*(fxx-this->f-lambda*dTg));
	  double dtmp;
	  if (sigma2*lambda<lambdanew) dtmp=sigma2*lambda;
	  else dtmp=lambdanew;
	  if (sigma1*lambda>dtmp) lambda=sigma1*lambda;
	  else lambda=dtmp;
	} else {
	  armijo_failure=0;
	}
      }
      
      // st->s = x_{k+1} - x_k 
      for(size_t i=0;i<this->dim;i++) s[i]=xx[i];
      for(size_t i=0;i<this->dim;i++) s[i]-=this->x[i];

      // M->x = x_{k+1} 
      for(size_t i=0;i<this->dim;i++) this->x[i]=xx[i];
      (*this->dfunc)(this->dim,this->x,this->gradient);
      this->f=fxx;
      
      // Infinite norm of g1 = d <- [P(x - g) - x] 
      for(size_t i=0;i<this->dim;i++) d[i]=this->x[i];
      for(size_t i=0;i<this->dim;i++) d[i]-=this->gradient[i];
      proj(d);
      for(size_t i=0;i<this->dim;i++) d[i]-=this->x[i];
      
      double maxval=fabs(d[0]);
      for(size_t i=1;i<this->dim;i++) {
	if (fabs(d[i])>maxval) {
	  maxval=fabs(d[i]);
	}
      }
      this->size=fabs(maxval);
      
      // st->y = - (g(x_{k+1}) - g(x_k)) 
      for(size_t i=0;i<this->dim;i++) y[i]-=this->gradient[i];

      m++;
      if (M<m) m=M;
      tail++;
      tail=tail % M;
      fvec[tail]=this->f;

      return 0;
    }
    
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

#endif

    public:

    mmin_constr_spg() {
      fmin=-1.0e99;
      tol=1.0e-4;
      alphamin=1.0e-30;
      alphamax=1.0e30;
      gamma=1.0e-4;
      sigma1=0.1;
      sigma2=0.9;
      M=10;
    }

    ~mmin_constr_spg() {
      fvec.clear();
    }

    /** \brief Minimum function value (default \f$ 10^{-99} \f$)
	
	If the function value is below this value, then the algorithm
	assumes that the function is not bounded and exits.
    */
    double fmin;

    /// Tolerance on infinite norm (default \f$ 10^{-4} \f$)
    double tol;

    /// Lower bound to spectral step size (default \f$ 10^{-30} \f$)
    double alphamin;

    /// Upper bound to spectral step size (default \f$ 10^{30} \f$)
    double alphamax;

    /// Sufficient decrease parameter (default \f$ 10^{-4} \f$)
    double gamma;

    /// Lower bound to the step length reduction (default 0.1)
    double sigma1;

    /// Upper bound to the step length reduction (default 0.9)
    double sigma2;

    /// Monotonicity parameter (M=1 forces monotonicity) (default 10)
    size_t M;
    
    /// Allocate memory
    virtual int allocate(const size_t n) {
      if (this->dim>0) free();
      xx.resize(n);
      d.resize(n);
      s.resize(n);
      y.resize(n);
      return mmin_constr<func_t,dfunc_t,hfunc_t,vec_t>::allocate(n);
    }
    
    /// Free previously allocated memory
    virtual int free() {
      xx.clear();
      d.clear();
      s.clear();
      y.clear();
      return 0;
    }
    
    /// Set the function, the initial guess, and the parameters
    virtual int set(func_t &fn, dfunc_t &dfn, vec_t &init) {
      
      int ret=mmin_constr<func_t,dfunc_t,hfunc_t,
      vec_t>::set(fn,dfn,init);
      
      // Turn initial guess into a feasible point
      proj(this->x);
      
      // Evaluate function and gradient
      this->f=(*this->func)(this->dim,this->x);
      (*this->dfunc)(this->dim,this->x,this->gradient);
      
      /* Infinite norm of d <- g1 = [P(x - g) - x] */
      o2scl::vector_copy(this->dim,this->x,d);
      for(size_t i=0;i<this->dim;i++) {
	d[i]-=this->gradient[i];
      }
      proj(d);
      for(size_t i=0;i<this->dim;i++) d[i]-=this->x[i];

      double maxval=fabs(d[0]);
      for(size_t i=1;i<this->dim;i++) {
	if (fabs(d[i])>maxval) {
	  maxval=fabs(d[i]);
	}
      }
      this->size=fabs(maxval);
      
      /* alpha_0 */
      maxval=fabs(this->gradient[0]);
      for(size_t i=1;i<this->dim;i++) {
	if (fabs(this->gradient[i])>maxval) {
	  maxval=fabs(this->gradient[i]);
	}
      }
      alpha=1.0/fabs(maxval);

      // FIXME: at the moment, this memory allocation is 
      // left hanging. It's taken care of the destructor,
      // but this isn't so elegant.
      fvec.resize(M);
      m=1;
      tail=0;
      fvec[tail]=this->f;

      return 0;
    }

    /// Restart the minimizer
    virtual int restart() {

      proj(this->x);

      // Evaluate function and gradient
      this->f=(*this->func)(this->dim,this->x);
      (*this->dfunc)(this->dim,this->x,this->gradient);

      // alpha_0 
      double maxval=fabs(this->gradient[0]);
      for(size_t i=1;i<this->dim;i++) {
	if (fabs(this->gradient[i])>maxval) {
	  maxval=fabs(this->gradient[i]);
	}
      }
      alpha=1.0/fabs(maxval);

      // Infinite norm of d <- g1 = [P(x - g) - x]
      o2scl::vector_copy(this->dim,this->x,d);
      for(size_t i=0;i<this->dim;i++) {
	d[i]-=this->gradient[i];
      }
      proj(d);
      for(size_t i=0;i<this->dim;i++) d[i]-=this->x[i];

      maxval=fabs(d[0]);
      for(size_t i=1;i<this->dim;i++) {
	if (fabs(d[i])>maxval) {
	  maxval=fabs(d[i]);
	}
      }
      this->size=fabs(maxval);

      m=1;
      tail=0;
      fvec[tail]=this->f;

      return 0;
    }

    /// Perform an iteration
    virtual int iterate() {
      double t;
      
      line_search();

      double bk=0.0;
      for(size_t i=0;i<this->dim;i++) bk+=s[i]*y[i];

      if (bk>=0.0) {
	alpha=alphamax;
      } else {
	double ak=0.0;
	for(size_t i=0;i<this->dim;i++) ak+=s[i]*s[i];
	ak=-ak/bk;
	double dtmp;
	if (alphamin>ak) dtmp=alphamin;
	else dtmp=ak;
	if (alphamax<dtmp) alpha=alphamax;
	else alpha=dtmp;
	//alpha=GSL_MIN(alphamax,GSL_MAX(alphamin,ak));
      }

      return 0;
    }

    /// See if we're finished
    virtual int is_optimal() {
      if (this->size>tol && this->f>fmin) {
	return GSL_CONTINUE;
      }
      return GSL_SUCCESS;
    }

    /// Return string denoting type ("mmin_constr_spg")
    const char *type() { return "mmin_constr_spg"; }

#ifndef DOXYGEN_INTERNAL

  private:
  
  mmin_constr_spg<func_t,dfunc_t,vec_t>
  (const mmin_constr_spg<func_t,dfunc_t,vec_t> &);
  mmin_constr_spg<func_t,dfunc_t,vec_t>& operator=
  (const mmin_constr_spg<func_t,dfunc_t,vec_t>&);

#endif
      
  };
  
}

#endif

