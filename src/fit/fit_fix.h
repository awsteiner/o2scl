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
#ifndef O2SCL_FIT_FIX_H
#define O2SCL_FIT_FIX_H

/** \file fit_fix.h
    \brief File defining \ref o2scl::fit_fix_pars
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/fit_base.h>
#include <o2scl/fit_nonlin.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multidimensional fitting class fixing some parameters and 
      varying others

      The number of trials used in the fit can be specified in
      the data member of the parent class \ref fit_base::ntrial
      associated with the fit_fix_pars object. Similarly for the
      verbosity parameter in \ref fit_base::verbose, the absolute
      tolerance in \ref fit_base::tol_abs, and the relative tolerance
      in \ref fit_base::tol_abs. These values are copied to the
      minimizer used by <tt>fit_fix_pars::mmin()</tt> during each call.
      After the minimizer is called, the value of \ref fit_base::ntrial
      associated with the \ref fit_fix_pars object is filled
      with the last number of trials required for the last
      minimization.
      \comment
      For some reason the reference to mmin() above doesn't work
      in doxygen.
      \endcomment

      Default template arguments
      - \c func_t - \ref gen_fit_funct\<\>
      - \c vec_t - \ref boost::numeric::ublas::vector \<double \>
      - \c mat_t - \ref boost::numeric::ublas::matrix \<double \>
  */
  template<class bool_vec_t, class func_t=gen_fit_funct<>, 
    class vec_t=boost::numeric::ublas::vector<double>,
    class mat_t=boost::numeric::ublas::matrix<double> >
    class fit_fix_pars : public fit_base<func_t,vec_t,mat_t>,
    public gen_fit_funct<vec_t,mat_t> {
    
  public:

  /// The generic fitter type
  typedef fit_base<fit_fix_pars<bool_vec_t,func_t,vec_t,mat_t>,
  vec_t,mat_t> base_fit_t;

  /// The default fitter type
  typedef fit_nonlin<fit_fix_pars<bool_vec_t,func_t,vec_t,mat_t>,
  vec_t,mat_t> def_fit_t;
  
  /** \brief Specify the member function pointer
   */
  fit_fix_pars() {
    fitp=&def_fit;
    expand_covar=false;
  }
    
  virtual ~fit_fix_pars() {}

  /** \brief If true, expand the covariance matrix to the
      larger space by filling with the identity matrix (default false)

      If this varable is false (the default), then the covariance
      matrix is computed in the smaller space which enumerates only
      the parameters which are not fixed. If this variable is true,
      then the covariance matrix must be a full \c np by \c np
      matrix (where \c np is the argument to \ref fit() or \ref
      fit_fix() ) and rows and columns which correspond with
      parameters which are fixed are replaced by elements from the
      identity matrix.

      The optimal parameters and \f$ \chi^2 \f$ reported by the
      fit are unchanged. 
  */
  bool expand_covar;
    
  /** \brief Fit the data specified in (xdat,ydat) to
      the function \c fitfun with the parameters in \c par.

      The covariance matrix for the parameters is returned in \c covar
      and the value of \f$ \chi^2 \f$ is returned in \c chi2.
  */
  virtual int fit(size_t np, vec_t &par, mat_t &covar, double &chi2,
		  func_t &fitfun) {
      
    x_par=&par;
    funcp=&fitfun;
    u_np=np;
    u_np_new=np;
    size_t nd=fitfun.get_ndata();

    u_par.resize(np);
    J.resize(nd,np);

    fitp->verbose=this->verbose;
    fitp->ntrial=this->ntrial;
    fitp->tol_rel=this->tol_rel;
    fitp->tol_abs=this->tol_abs;
    fitp->fit(u_np_new,par,covar,chi2,*this);

    u_par.clear();
    J.clear();

    return 0;
  }

  /** \brief Fit function \c func while fixing some parameters as
      specified in \c fix.

      \note When some parameters are fixed, each fixed parameter 
      leads to a row and column in the covariance matrix which
      is unused by this function. If there are <tt>npar2<npar</tt>
      unfixed parameters, then only the first <tt>npar2</tt> rows
      and columns of \c covar are filled.
  */
  virtual int fit_fix(size_t np, vec_t &par, mat_t &covar, double &chi2,
		      func_t &fitfun, bool_vec_t &fix) {

    x_par=&par;
    funcp=&fitfun;
    u_np=np;
    fix_par=&fix;
    size_t nd=fitfun.get_ndata();

    // Count number of new parameters
    u_np_new=0;
    for(size_t i=0;i<np;i++) {
      if (fix[i]==false) u_np_new++;
    }
    if (u_np_new==0) return 0;

    // Allocate space
    u_par.resize(np);
    u_par_new.resize(u_np_new);
    J.resize(nd,np);

    // Copy to smaller vector containing only parameters to be varied
    size_t j=0;
    for(size_t i=0;i<np;i++) {
      if (fix[i]==false) {
	u_par_new[j]=par[i];
	j++;
      }
    }

    fitp->verbose=this->verbose;
    fitp->ntrial=this->ntrial;
    fitp->tol_rel=this->tol_rel;
    fitp->tol_abs=this->tol_abs;
    fitp->fit(u_np_new,u_par_new,covar,chi2,*this);

    // Copy optimal parameter set back to initial guess object
    j=0;
    for(size_t i=0;i<np;i++) {
      if (fix[i]==false) {
	par[i]=u_par_new[j];
	j++;
      }
    }

    // If required, expand covariance matrix
    if (expand_covar && u_np_new<u_np) {
      int i_new=((int)u_np_new)-1;
      for(int i=((int)np)-1;i>=0;i--) {
	int k_new=((int)u_np_new)-1;
	for(int k=((int)np)-1;k>=0;k--) {
	  if (fix[i]==false && fix[k]==false) {
	    if (i_new<0 || k_new<0 || 
		i_new>=((int)u_np_new) || k_new>=((int)u_np_new)) {
	      O2SCL_ERR2("Covariance alignment in ",
			 "fit_fix::fit_fix().",exc_esanity);
	    }
	    covar(i,k)=covar(i_new,k_new);
	    k_new--;
	  } else if (i==k) {
	    covar(i,k)=1.0;
	  } else {
	    covar(i,k)=0.0;
	  }
	}
	if (fix[i]==false) i_new--;
      }
    }

    // Free space
    u_par_new.clear();
    u_par.clear();
    J.clear();

    return 0;
  }

  /// Change the base fitter
  int set_fit(base_fit_t &fitter) {
    fitp=&fitter;
    return 0;
  }
      
  /// The default base fitter
  def_fit_t def_fit;

  /// \name Reimplementation of \ref gen_fit_funct
  //@{
  /// The function to return the number of data points
  virtual size_t get_ndata() {
    return (*funcp).get_ndata();
  }
      
  /** \brief The function computing deviations

      This function operates in the smaller space of size np_new.
  */
  virtual void operator()(size_t np_new, const vec_t &par_new, 
			  size_t nd, vec_t &f) {

    // Variable 'i' runs over the user space and 'i_new' runs
    // over the smaller space of size np_new.
    size_t i_new=0;
    for(size_t i=0;i<u_np;i++) {
      if (u_np_new<u_np && (*fix_par)[i]==true) {
	u_par[i]=(*x_par)[i];
      } else {
	u_par[i]=par_new[i_new];
	i_new++;
      }
    }
    if (i_new!=np_new) {
      O2SCL_ERR("Alignment failure 1 in fit_fix::operator().",
		exc_esanity);
    }

    // Call the function in the larger space of size np
    (*funcp)(u_np,u_par,nd,f);

    return; 
  }

  /** \brief The function computing the Jacobian
   */
  virtual void jac(size_t np_new, vec_t &par_new, size_t nd, vec_t &f, 
		   mat_t &J_new) {
    
    size_t i_new=0;
    for(size_t i=0;i<u_np;i++) {
      if (u_np_new<u_np && (*fix_par)[i]==true) {
	u_par[i]=(*x_par)[i];
      } else {
	u_par[i]=par_new[i_new];
	i_new++;
      }
    }
    if (i_new!=np_new) {
      O2SCL_ERR("Alignment failure 2 in fit_fix::operator().",
		exc_esanity);
    }
    
    // Call the Jacobian in the larger space of size (nd,u_np)
    (*funcp).jac(u_np,u_par,nd,f,J);

    // Copy the Jacobian over to the smaller space, ignoring
    // rows corresponding to parameters which are fixed
    for(size_t id=0;id<nd;id++) {
      i_new=0;
      for(size_t i=0;i<u_np;i++) {
	if (u_np_new==u_np || (*fix_par)[i]==false) {
	  J_new(id,i_new)=J(id,i);
	  i_new++;
	}
      }
    }

    return;
  }
  //@}
    
#ifndef DOXYGEN_INTERNAL
      
  protected:

  /// Temporary vector to store full parameter list of size u_np
  vec_t u_par;

  /// Vector for smaller parameter list of size u_np_new
  vec_t u_par_new;
  
  /// Jacobian in the user space of size (nd,u_np)
  mat_t J;

  /// The fitter
  base_fit_t *fitp;

  /// The user-specified function
  func_t *funcp;

  /// The user-specified number of variables
  size_t u_np;

  /// The new number of variables
  size_t u_np_new;
    
  /// Specify which parameters to fix (vector of size u_np)
  bool_vec_t *fix_par;
    
  /// The user-specified initial vector of size u_np
  vec_t *x_par;

  private:
    
  fit_fix_pars(const fit_fix_pars &);
  fit_fix_pars& operator=(const fit_fix_pars&);

#endif

  };


#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
