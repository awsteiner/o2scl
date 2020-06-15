/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2020, Andrew W. Steiner

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
#ifndef O2SCL_MMIN_FIX_H
#define O2SCL_MMIN_FIX_H

/** \file mmin_fix.h
    \brief File defining \ref o2scl::mmin_fix_params
*/

#include <vector>

#include <o2scl/mmin.h>
#include <o2scl/mmin_simp2.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multidimensional minimizer fixing some parameters and 
      varying others

      This class allows one to min a function after having fixed
      some of the parameters. The parameters which should be fixed are
      specified through a <tt>bool</tt> vector. This class performs the
      extra bookkeeping associated with reordering the parameters and
      performs the minimization with a separate minimizer object. This
      class is most useful for minimization problems which do not use
      information about the gradient.
      
      The number of trials used in the minimizer can be specified in
      the data member of the parent class \ref mmin_base::ntrial
      associated with the \ref o2scl::mmin_fix_params object.
      Similarly for the verbosity parameter in \ref
      mmin_base::verbose, the absolute tolerance in \ref
      mmin_base::tol_abs, and the relative tolerance in \ref
      mmin_base::tol_abs. These values are copied to the minimizer
      used by \ref mmin_fix_params::mmin() during each call. After the
      minimizer is called, the value of \ref mmin_base::ntrial
      associated with the \ref mmin_fix_params object is filled with
      the last number of trials required for the last minimization.

      \verbatim embed:rst
      See an example for the usage of this class in 
      :ref:`Minimizer fixing variables example`.
      \endverbatim

      \comment
      We cannot really do a version of mmin_de() for this class
      because there's no good way to rewrite the gradient
      without the contribution from the parameters which are 
      fixed.
      \endcomment
  */
  template<class func_t=multi_funct, 
    class vec_t=boost::numeric::ublas::vector<double> > 
    class mmin_fix_params : public mmin_base<func_t,func_t,vec_t> {
    
  public:
  
  typedef boost::numeric::ublas::vector<double> ubvector;

  /// The generic minimizer type
  typedef mmin_base<mmin_fix_params<func_t,vec_t>,
  mmin_fix_params<func_t,vec_t>,vec_t> base_mmin_t;

  /// The default minimizer type
  typedef mmin_simp2<mmin_fix_params<func_t,vec_t>,vec_t> def_mmin_t;

  /** \brief Specify the member function pointer
   */
  mmin_fix_params() {
    mmp=&def_mmin;
  }
    
  virtual ~mmin_fix_params() {}
    
  /** \brief Calculate the minimum \c min of \c func w.r.t. the
      array \c x of size \c nvar.
  */
  virtual int mmin(size_t nvar, vec_t &x, double &fmin, 
		   func_t &func) {

    if (fixp.size()<nvar) {
      fixp.resize(nvar);
      for(size_t i=0;i<nvar;i++) fixp[i]=false;
    }
    
    // Copy arguments for later use
    xp=&x;
    funcp=&func;
    unv=nvar;

    nv_new=0;
    for(size_t i=0;i<nvar;i++) {
      if (fixp[i]==false) nv_new++;
    }
    if (nv_new==0) return 0;

    // Copy initial guess to new format
    size_t j=0;
    ubvector xnew(nv_new);
    for(size_t i=0;i<nvar;i++) {
      if (fixp[i]==false) {
	xnew[j]=x[i];
	j++;
      }
    }

    // Perform minimization
    mmp->verbose=this->verbose;
    mmp->ntrial=this->ntrial;
    mmp->tol_rel=this->tol_rel;
    mmp->tol_abs=this->tol_abs;

    int ret=mmp->mmin(nv_new,xnew,fmin,*this);
    if (ret!=0) {
      O2SCL_ERR("Minimizer failed in mmin_fix_params::mmin_fix().",ret);
    }

    // Copy final result back to original format
    j=0;
    for(size_t i=0;i<nvar;i++) {
      if (fixp[i]==false) {
	x[i]=xnew[j];
	j++;
      }
    }
      
    this->last_ntrial=mmp->last_ntrial;

    return ret;
  }

  template<class bool_vec_t> void set_fix(size_t n, bool_vec_t &fix) {
    fixp.resize(n);
    for(size_t i=0;i<n;i++) fixp[i]=fix[i];
    return;
  }

  /** \brief Calculate the minimum of \c func while fixing 
      some parameters as specified in \c fix.

      If all of entries <tt>fix[0], fix[1], ... fix[nvar-1]</tt>
      are true, then this function assumes all of the parameters
      are fixed and that there is no minimization to be performed.
      In this case, it will return 0 for success without calling
      the error handler. 
  */
  template<class bool_vec_t>
  int mmin_fix(size_t nvar, ubvector &x, double &fmin, 
	       bool_vec_t &fix, multi_funct &func) {
    
    fixp.resize(nvar);
    for(size_t i=0;i<nvar;i++) fixp[i]=fix[i];

    return mmin(nvar,x,fmin,func);
  }

  /// Change the base minimizer
  int set_mmin(base_mmin_t &min) {
		 
    mmp=&min;
    return 0;
  }

  /// The default base minimizer
  def_mmin_t def_mmin;
    
  /** \brief The new function to send to the minimizer
   */
  virtual double operator()(size_t nv, const vec_t &x) {
    
    ubvector tmp(unv);
    size_t j=0;
    for(size_t i=0;i<unv;i++) {
      if (nv_new<unv && fixp[i]==true) {
	tmp[i]=(*xp)[i];
      } else {
	tmp[i]=x[j];
	j++;
      }
    }
    return (*funcp)(unv,tmp);
  }
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
  /// The minimizer
  base_mmin_t *mmp;

  /// The user-specified function
  func_t *funcp;

  /// The user-specified number of variables
  size_t unv;

  /// The new number of variables
  size_t nv_new;
    
  /// Specify which parameters to fix
  std::vector<bool> fixp;
    
  /// The user-specified initial vector
  vec_t *xp;
  
  private:
  
  mmin_fix_params(const mmin_fix_params &);
  mmin_fix_params& operator=(const mmin_fix_params&);

#endif

  };


#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
