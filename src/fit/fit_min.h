/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
*/
#ifndef O2SCL_FIT_MIN_H
#define O2SCL_FIT_MIN_H

/** \file fit_min.h
    \brief File defining \ref o2scl::fit_min
*/

#include <o2scl/mmin.h>
#include <o2scl/multi_funct.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/fit_base.h>
#include <o2scl/fit_nonlin.h>

namespace o2scl {

  /** \brief Non-linear least-squares fitting class with generic minimizer
    
      This minimizes a generic fitting function using any \ref
      o2scl::mmin_base object, and then uses the GSL routines to
      calculate the uncertainties in the parameters and the covariance
      matrix.

      This can be useful for fitting problems which might be better
      handled by more complex minimizers than those that are used in
      \ref o2scl::fit_nonlin. For problems with many local minima near
      the global minimum, using a \ref o2scl::anneal_base object with
      this class can sometimes produce better results than \ref
      o2scl::fit_nonlin.

      Default template arguments
      - \c func_t - \ref gen_fit_funct\<\>
      - \c vec_t - \ref boost::numeric::ublas::vector \<double \>
      - \c mat_t - \ref boost::numeric::ublas::matrix \<double \>
  */
  template<class func_t=gen_fit_funct<>, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class fit_min : 
    public fit_base<func_t,vec_t,mat_t>, public fit_nonlin_b<vec_t,mat_t> {    

  public:
  
    fit_min() {
      min_set=false;
      mmp=&def_mmin;
    }

    virtual ~fit_min() {}

    /** \brief Fit the data specified in (xdat,ydat) to
	the function \c fitfun with the parameters in \c par.
	
	The covariance matrix for the parameters is returned in \c
	covar and the value of \f$ \chi^2 \f$ is returned in \c chi2.
    */
    virtual int fit(size_t npar, vec_t &par, mat_t &covar, double &chi2,
		    func_t &fitfun) {

      func=&fitfun;

      size_t ndata=fitfun.get_ndata();
      fval.resize(ndata);

      // ---------------------------------------------------
      // First minimize with the specified mmin object
      // ---------------------------------------------------
      
      multi_funct mmf=std::bind(std::mem_fn<double(size_t, const vec_t &)>
				  (&fit_min::min_func),
				  this,std::placeholders::_1,
				  std::placeholders::_2);
      //multi_funct_mfptr<fit_min> 
      //mmf(this,&fit_min::min_func);
      
      double dtemp;
      mmp->mmin(npar,par,dtemp,mmf);

      // ---------------------------------------------------
      // Now compute the Jacobian and do the QR decomposition 
      // ---------------------------------------------------

      // Allocate space
      int signum;
      permutation perm(npar);
      vec_t diag, tau, norm;
      diag.resize(npar);
      if (ndata<npar) {
	tau.resize(ndata);
      } else {
	tau.resize(npar);
      }
      norm.resize(npar);
      mat_t J, r;
      J.resize(ndata,npar);
      r.resize(ndata,npar);
      
      // Evaluate function and Jacobian at optimal parameter values
      fitfun(npar,par,ndata,fval);
      fitfun.jac(npar,par,ndata,fval,J);

      double fnorm=o2scl_cblas::dnrm2(ndata,fval);

      // Use scaled version
      this->compute_diag(npar,ndata,J,diag);

      double xnorm=this->scaled_enorm(diag,npar,fval);
      double delta=this->compute_delta(diag,npar,fval);
      matrix_copy(ndata,npar,J,r);

      o2scl_linalg::QRPT_decomp(ndata,npar,r,tau,perm,signum,norm);
      
      // ---------------------------------------------------
      // Compute the covariance matrix
      // ---------------------------------------------------

      this->covariance(ndata,npar,J,covar,norm,r,tau,perm,
		       this->tol_rel_covar);

      chi2=fnorm*fnorm;

      // Free previously allocated space

      diag.clear();
      tau.clear();
      norm.clear();
      J.clear();
      r.clear();
      fval.clear();

      return 0;
    }
    
    /** \brief Set the mmin object to use (default is 
	of type \ref o2scl::mmin_simp2)
     */
    int set_mmin(mmin_base<multi_funct> &mm) {
      mmp=&mm;
      min_set=true;
      return 0;
    }
    
    /// The default minimizer
    mmin_simp2<multi_funct> def_mmin;

    /// Return string denoting type ("fit_min")
    virtual const char *type() { return "fit_min"; }

#ifndef DOXYGEN_INTERNAL

    protected:

    /// Storage for the function values
    vec_t fval;

    /// Pointer to the user-specified fitting function
    func_t *func;

    /// The minimizer
    mmin_base<multi_funct> *mmp;
    
    /// True if the minimizer has been set by the user
    bool min_set;
    
    /// The function to minimize, \f$ \chi^2 \f$
    double min_func(size_t np, const vec_t &xp) {
      double ret=0.0;
      (*func)(np,xp,func->get_ndata(),fval);
      for(size_t i=0;i<func->get_ndata();i++) {
	ret+=fval[i]*fval[i];
      }
      return ret;
    }
    
#endif
  
  };

}

#endif
