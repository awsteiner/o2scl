/*
  -------------------------------------------------------------------
   
  Copyright (C) 2006-2014, Andrew W. Steiner

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
/*------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * tools/tools_diff.c
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
 * since: May, 23, 2004
 *
 * $Id: tools_diff.c,v 1.2 2004/06/04 21:26:23 akiles Exp $
 *------------------------------------------------------------------------*/
#ifndef O2SCL_OOL_CONSTR_MMIN_H
#define O2SCL_OOL_CONSTR_MMIN_H

/** \file mmin_constr.h
    \brief File defining \ref o2scl::mmin_constr and associated function
    objects
*/
#include <o2scl/multi_funct.h>
#include <o2scl/mmin.h>
#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Hessian product function for mmin_constr [abstract base]
   */
  template<class vec_t=boost::numeric::ublas::vector<double> > 
    class ool_hfunct {
    
    public:  

    ool_hfunct() {}

    virtual ~ool_hfunct() {}

    /** \brief Evaluate \f$ H(x) \cdot v \f$
     */
    virtual int operator()(size_t nv, const vec_t &x, const vec_t &v, 
			   vec_t &hv)=0;

#ifndef DOXYGEN_NO_O2NS

    private:

    ool_hfunct(const ool_hfunct &);
    ool_hfunct& operator=(const ool_hfunct&);

#endif

  };

  /** \brief A hessian product supplied by a function pointer
   */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class ool_hfunct_fptr : public ool_hfunct<vec_t> {
    
    public:

    /** \brief Specify the function pointer
     */
    ool_hfunct_fptr(int (*fp)(size_t nv, const vec_t &x, const vec_t &v, 
			      vec_t &hv)) {
      fptr=fp;
    }
    
    
    virtual ~ool_hfunct_fptr() {};

    /** \brief Evaluate \f$ H(x) \cdot v \f$
     */
    virtual int operator()(size_t nv, const vec_t &x, const vec_t &v, 
			   vec_t &hv) {
      return fptr(nv,x,v,hv);
    }
    
#ifndef DOXYGEN_INTERNAL

    protected:
    
    /// Store the function pointer
    int (*fptr)(size_t nv, const vec_t &x, const vec_t &v, vec_t &hv);
    
    ool_hfunct_fptr() {}

#ifndef DOXYGEN_NO_O2NS
#endif

    private:

    ool_hfunct_fptr(const ool_hfunct_fptr &);
    ool_hfunct_fptr& operator=(const ool_hfunct_fptr&);

#endif

  };

  /** \brief A hessian product supplied by a member function pointer
   */
  template<class tclass, class vec_t=boost::numeric::ublas::vector<double> >
    class ool_hfunct_mfptr : public ool_hfunct<vec_t> {
    public:
  
    /** \brief Specify the class instance and member function pointer
     */
    ool_hfunct_mfptr(tclass *tp, int (tclass::*fp)
		     (size_t nv, const vec_t &x, const vec_t &v, 
		      vec_t &hv)) {
      tptr=tp;
      fptr=fp;
    }
    
    virtual ~ool_hfunct_mfptr() {}
    
    /** \brief Evaluate \f$ H(x) \cdot v \f$
     */
    virtual int operator()(size_t nv, const vec_t &x, const vec_t &v, 
			   vec_t &hv) {
      return (*tptr.*fptr)(nv,x,v,hv);
    }

#ifndef DOXYGEN_INTERNAL
  
    protected:
  
    /// Store the function pointer
    int (tclass::*fptr)(size_t nv, const vec_t &x, const vec_t &v, 
			vec_t &hv);

    /// Store a pointer to the class instance
    tclass *tptr;
  
#ifndef DOXYGEN_NO_O2NS
#endif

    private:

    ool_hfunct_mfptr(const ool_hfunct_mfptr &);
    ool_hfunct_mfptr& operator=(const ool_hfunct_mfptr&);

#endif
    
  };
  
  /** \brief Constrained multidimensional minimization (OOL) [abstract base]

      \future Implement automatic computations of \gradient and Hessian
      \future Construct a more difficult example for the "examples" directory
      \future Finish mmin() interface
      \future Implement a direct computation of the hessian as the 
      jacobian of the gradient
  */
  template<class func_t, class dfunc_t=func_t, 
    class hfunc_t=func_t, class vec_t=boost::numeric::ublas::vector<double> >
    class mmin_constr :  public mmin_base<func_t,dfunc_t,vec_t> {
    
    public:
    
#ifndef DOXYGEN_INTERNAL

    protected:
    
    /// The current function value
    double f;
    /// Desc
    double size;
    
    /// The current minimum vector
    vec_t x;
    /// The current gradient vector
    vec_t gradient;
    /// Desc
    vec_t dx;
    
    /// Number of function evaluations
    size_t fcount;
    /// Number of gradient evaluations
    size_t gcount;
    /// Number of Hessian evaluations
    size_t hcount;

    /// Number of parameters
    size_t dim;
    /// Number of constraints
    size_t nconstr;
    /// User-supplied function
    func_t *func;
    /// Gradient function
    dfunc_t *dfunc;
    /// Hessian function
    hfunc_t *hfunc;

    /// Lower bound constraints
    vec_t L;
    /// Upper bound constraints
    vec_t U;

    /// If true, the algorithm requires the hessian vector product
    bool requires_hess;

    /// Shrink vector \c V from the full to the reduced space
    void shrink(const size_t nind, gsl_vector_uint *Ind, 
		const vec_t &V) {
      size_t ii, indii;
      
      for( ii = 0; ii < nind; ii++ ) {
	indii = gsl_vector_uint_get(Ind,ii);
	
	double tmp=V[ii];
	V[ii]=V[indii];
	V[indii]=tmp;
      }
    }

    /// Expand vector \c V from the reduced to the full space
    void expand(const size_t nind, gsl_vector_uint *Ind, 
		const vec_t &V) {

      size_t jj, ii, indii;
      
      ii = nind;
      for( jj = 0; jj < nind; jj++ ) {
	ii--;
	indii = gsl_vector_uint_get( Ind, ii );
	
	double tmp=V[ii];
	V[ii]=V[indii];
	V[indii]=tmp;
      }
    }
    
    /// Evaluate the objective function from the reduced space
    double calc_f(const size_t nind, gsl_vector_uint *Ind, 
		  vec_t &X, vec_t &Xc) {
      
      const size_t missing=this->dim-nind;
      double ftmp;
      
      if (missing>0) {
	
	// Complete X with values from Xc 
	for(size_t i=nind;i<nind+missing;i++) {
	  X[i]=Xc[i];
	}
	
	// Expand to full space
	expand(nind,Ind,X);
      }

      // Compute f 
      func(dim,X,ftmp);
      
      if(missing>0)	{
	// Shrink to reduced space 
	shrink(nind,Ind,X);
      }

      return f;
    }
    
    /// Compute gradient in the reduced space
    int calc_g(const size_t nind, gsl_vector_uint *Ind, vec_t &X,
	       vec_t &Xc, vec_t &G) {
      
      const size_t missing=this->dim-nind;
      
      if( missing > 0 )	{

	// Complete X with values from Xc 
	for(size_t i=nind;i<nind+missing;i++) {
	  X[i]=Xc[i];
	}
	
	// Expand to full space
	expand(nind,Ind,X);
      }

      // Compute gradient 
      dfunc_t(dim,X,G);
      
      if (missing>0) {
	// Shrink to reduced space 
	shrink(nind,Ind,X);
	shrink(nind,Ind,G);
      }

      return 0;
    }

    /// Evaluate a hessian times a vector from the reduced space
    int calc_Hv(const size_t nind, gsl_vector_uint *Ind,
		vec_t &X, vec_t &Xc, vec_t &V, vec_t &Hv) {
      
      const size_t missing=this->ndim-nind;
      
      if (missing>0) {
	// Complete X with values from Xc and set V to zero
	for(size_t i=nind;i<nind+missing;i++) {
	  X[i]=Xc[i];
	  V[i]=0.0;
	}
	/// Expand to full space
	expand(nind,Ind,X);
	expand(nind,Ind,V);
      }

      // Compute gradient
      hfunc(this->dim,X,V,Hv);

      if (missing>0) {
	// Shrink to reduced space
	shrink(nind,Ind,X);
	shrink(nind,Ind,V);
	shrink(nind,Ind,Hv);
      }

      return 0;
    }
    
#endif

    public:
    
    mmin_constr() {
      dim=0;
      nconstr=0;
      requires_hess=false;
    }
    
    virtual ~mmin_constr() {
    }
    
    /// Allocate memory
    virtual int allocate(const size_t n) {
      int status;
      
      x.resize(n);
      gradient.resize(n);
      dx.resize(n);
      dim=n;

      for(size_t i=0;i<n;i++) {
	x[i]=0.0;
	gradient[i]=0.0;
	dx[i]=0.0;
      }

      return 0;
    }
    
    /// Restart the minimizer
    virtual int restart() {
   
      // Reset dx vector 
      for(size_t i=0;i<dim;i++) dx[i]=0.0;
      
      // Restart evaluation counters 
      fcount = 0;
      gcount = 0;
      hcount = 0;
      
      return 0;
    }
    
    /// Set the function, the gradient, and the initial guess
    virtual int set(func_t &fn, dfunc_t &dfn, vec_t &init) {
      
      o2scl::vector_copy(dim,init,x);
      for(size_t i=0;i<dim;i++) dx[i]=0.0;
      
      func=&fn;
      dfunc=&dfn;
      
      // Start evaluation conunters 
      fcount = 0;
      gcount = 0;
      hcount = 0;

      return 0;
    }
    
    /** \brief Set the function, the gradient, the Hessian product,
	and the initial guess
    */
    virtual int set_hess(func_t &fn, dfunc_t &dfn, hfunc_t &hfn,
			 vec_t &init) {
      
      o2scl::vector_copy(dim,init,x);
      for(size_t i=0;i<dim;i++) dx[i]=0.0;
      
      func=&fn;
      dfunc=&dfn;
      hfunc=&hfn;
      
      // Start evaluation conunters 
      fcount = 0;
      gcount = 0;
      hcount = 0;

      return 0;
    }
    
    /// Set the constraints
    virtual int set_constraints(size_t nc, vec_t &lower, vec_t &upper) {
      L.resize(nc);
      U.resize(nc);
      nconstr=nc;
      o2scl::vector_copy(nc,lower,L);
      o2scl::vector_copy(nc,upper,U);
      return 0;
    }
    
    /// Perform an iteration
    virtual int iterate()=0;
    
    /// See if we're finished
    virtual int is_optimal()=0;
    
    /** \brief Calculate the minimum \c min of \c func w.r.t. the
	array \c x of size \c nvar.

	\note This is unimplemented.
    */
    virtual int mmin(size_t nvar, vec_t &xx, double &fmin, 
		     func_t &ff) 
    {
      O2SCL_ERR_RET("Not yet implemented mmin_constr::mmin().",
		    exc_eunimpl);
    }

    /** \brief Calculate the minimum \c min of \c ff
	w.r.t. the array \c x of size \c nvar with gradient
	\c df and hessian vector product \c hf
    */
    virtual int mmin_hess(size_t nvar, vec_t &xx, double &fmin, 
			  func_t &ff, dfunc_t &df, hfunc_t &hf)
    {
      
      int status;
      allocate(nvar);
      if (requires_hess) set_hess(ff,df,hf,xx);
      else set(ff,df,xx);
      int it=0;
      do {
	status=iterate();
	status=is_optimal();
	it++;
      } while (it<this->ntrial && status==gsl_continue);
      
      for(size_t i=0;i<nvar;i++) xx[i]=this->x[i];
      fmin=this->f;
      this->last_ntrial=it;

      return 0;
    }

    /** \brief Calculate the minimum \c min of \c func
	w.r.t. the array \c x of size \c nvar with gradient
	\c dfunc
    */
    virtual int mmin_de(size_t nvar, vec_t &xx, double &fmin, 
			func_t &ff, dfunc_t &df) {
      
      int status;
      allocate(nvar);
      set(ff,df,xx);
      int it=0;
      do {
	status=iterate();
	status=is_optimal();
	it++;
      } while (it<this->ntrial && status==gsl_continue);

      for(size_t i=0;i<nvar;i++) xx[i]=this->x[i];
      fmin=this->f;
      this->last_ntrial=it;

      return 0;
    }
      
    /// Return string denoting type ("mmin_constr")
    const char *type() { return "mmin_constr"; }

#ifndef DOXYGEN_INTERNAL

  private:
  
  mmin_constr<func_t,dfunc_t,hfunc_t,vec_t>
  (const mmin_constr<func_t,dfunc_t,hfunc_t,vec_t> &);
  mmin_constr<func_t,dfunc_t,hfunc_t,vec_t>& operator=
  (const mmin_constr<func_t,dfunc_t,hfunc_t,vec_t>&);

#endif
      
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

