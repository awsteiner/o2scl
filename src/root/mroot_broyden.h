/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2023, Andrew W. Steiner
  
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

#ifndef MROOT_BROYDEN_H
#define MROOT_BROYDEN_H

/** \file mroot_broyden.h
    \brief File defining \ref o2scl::mroot_broyden 
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/mroot.h>
#include <o2scl/permutation.h>
#include <o2scl/lu.h>

namespace o2scl {
  
  /** \brief Multidimensional root-finding using Broyden's method (GSL)

      Experimental.

      \verbatim embed:rst
      Based on [Broyden65]_.
      \endverbatim
  */
  template<class func_t=mm_funct, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double>, 
    class jfunc_t=jac_funct> 
    class mroot_broyden : public mroot<func_t,vec_t,jfunc_t> {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

    protected:

    /// Desc
    ubmatrix H;
  
    /// LU decomposition
    ubmatrix lu;

    /// Permutation object for the LU decomposition
    permutation perm;

    /// Desc
    ubvector v;

    /// Desc
    ubvector w;

    /// Desc
    ubvector y;

    /// Desc
    ubvector p;

    /// Desc
    ubvector fnew;

    /// Desc
    ubvector x_trial;

    /// Desc
    double phi;

    /// Stepsize vector
    vec_t dx_int;

    /// Function value vector
    vec_t f_int;

    /// A pointer to the user-specified function
    func_t *user_func;

    /// Function values
    vec_t *user_f;

    /// Initial guess and current solution
    vec_t *user_x;

    /// Initial and current step
    vec_t *user_dx;

    /// Number of variables
    size_t user_nvar;

    /// Size of memory allocated
    size_t mem_size;

    /// Jacobian
    jacobian<func_t,vec_t,mat_t> *ajac;

    /** \brief Clear allocated vectors and matrices

        This function is called by set() before each solve.
     */
    void clear() {
      if (mem_size>0) {
	for(size_t i=0;i<lu.size1();i++) {
	  for(size_t j=0;j<lu.size2();j++) {
	    lu(i,j)=0.0;
	  }
	}
	for(size_t i=0;i<mem_size;i++) {
	  perm[i]=0;
	}
	for(size_t i=0;i<H.size1();i++) {
	  for(size_t j=0;j<H.size2();j++) {
	    H(i,j)=0.0;
	  }
	}
	for(size_t i=0;i<v.size();i++) v[i]=0.0;
	for(size_t i=0;i<w.size();i++) w[i]=0.0;
	for(size_t i=0;i<y.size();i++) y[i]=0.0;
	for(size_t i=0;i<fnew.size();i++) fnew[i]=0.0;
	for(size_t i=0;i<x_trial.size();i++) x_trial[i]=0.0;
	for(size_t i=0;i<p.size();i++) p[i]=0.0;
      }
      return;
    }

    public:

    mroot_broyden() {
      mem_size=0;
      ajac=&def_jac;
      def_jac.set_epsrel(sqrt(std::numeric_limits<double>::epsilon()));
    }
    
    /// Default Jacobian object
    jacobian_gsl<func_t,vec_t,mat_t> def_jac;
    
    virtual ~mroot_broyden() {
    }

    /// Allocate memory 
    void allocate(size_t n) {

      if (n!=mem_size) {

	lu.resize(n,n);
	perm.resize(n);
	H.resize(n,n);
	v.resize(n);
	w.resize(n);
	y.resize(n);
	fnew.resize(n);
	x_trial.resize(n);
	p.resize(n);
	dx_int.resize(n);
	f_int.resize(n);

	mem_size=n;

      }

      return;
    }

    /// Euclidean norm
    double enorm(size_t nvar, const vec_t &ff) {
      double e2=0;
      size_t i;
      for (i=0;i<nvar;i++) {
	double fi=ff[i];
	e2+=fi*fi;
      }
      return sqrt(e2);
    }

    /** \brief Set the function, initial guess, and provide vectors to
	store function values and stepsize

	The initial values of \c f and \c dx are ignored.
    */
    void set(func_t &func, size_t nvar, vec_t &x, vec_t &f, vec_t &dx) {

      if (nvar!=mem_size) allocate(nvar);
      clear();

      user_func=&func;
      user_nvar=nvar;
      user_x=&x;
      user_f=&f;
      user_dx=&dx;

      size_t i, j;
      int signum=0;
      
      func(nvar,x,f);
      
      ajac->set_function(func);
      (*ajac)(nvar,x,nvar,f,lu);
      
      o2scl_linalg::LU_decomp<>(nvar,lu,perm,signum);
      o2scl_linalg::LU_invert<ubmatrix,ubmatrix,ubmatrix_column>
      (nvar,lu,perm,H);
      
      for (i=0;i<nvar;i++) {
	for (j=0;j<nvar;j++) {
	  H(i,j)=-H(i,j);
	}
      }
      
      for(size_t iiq=0;iiq<dx.size();iiq++) dx[iiq]=0.0;
      
      phi=enorm(nvar,f);
      
      return;
    }

    /// Perform an iteration
    int iterate() {

      func_t &func=*user_func;
      size_t nvar=user_nvar;
      vec_t &x=*user_x;
      vec_t &f=*user_f;
      vec_t &dx=*user_dx;
    
      double phi0, phi1, t, lambda;
      
      size_t i, j, iter;
      
      /* p=H f */
      
      for (i=0;i<nvar;i++) {
	double sum=0;
	for (j=0;j<nvar;j++) {
	  sum+=H(i,j)*f[j];
	}
	p[i]=sum;
      }
      
      t=1.0;
      iter=0;
      
      phi0=phi;

      bool new_step=true;
      while (new_step) {
	new_step=false;
      
	for (i=0;i<nvar;i++) {
	  x_trial[i]=x[i]+t*p[i];
	}
      
	{ 
	  int status=func(nvar,x_trial,fnew);
	  if (status != success) {
	    return exc_ebadfunc;
	  }
	}
      
	phi1=enorm(nvar,fnew);
      
	iter++;
      
	if (phi1>phi0 && iter<10 && t>0.1) {

	  // [GSL] Full step goes uphill, take a reduced step instead.
	  double theta=phi1/phi0;
	  t*=(sqrt(1.0+6.0*theta)-1.0)/(3.0*theta);
	  new_step=true;
	}
      
	if (new_step==false && phi1>phi0) {

	  // [GSL] Need to recompute Jacobian
	  int signum=0;
	  
	  (*ajac)(nvar,x,nvar,f,lu);
	  
	  for (i=0;i<nvar;i++) {
	    for (j=0;j<nvar;j++) {
	      lu(i,j)=-lu(i,j);
	    }
	  }
	  
	  o2scl_linalg::LU_decomp(nvar,lu,perm,signum);
	  o2scl_linalg::LU_invert<ubmatrix,ubmatrix,ubmatrix_column>
	    (nvar,lu,perm,H);
	  o2scl_linalg::LU_solve(nvar,lu,perm,f,p);
	  
	  t=1.0;
	  
	  for (i=0;i<nvar;i++) {
	    x_trial[i]=x[i]+t*p[i];
	  }
	  
	  {
	    int status=func(nvar,x_trial,fnew);
	    if (status != success) {
	      return exc_ebadfunc;
	    }
	  }
	  
	  phi1=enorm(nvar,fnew);
	}

      }
      
      /* y=f'-f */
      for (i=0;i<nvar;i++) {
	y[i]=fnew[i]-f[i];
      }
      
      /* v=H y */
      for (i=0;i<nvar;i++) {
	double sum=0.0;
	for (j=0;j<nvar;j++) {
	  sum+=H(i,j)*y[j];
	}
	v[i]=sum;
      }
      
      /* lambda=p dot v */
      lambda=0.0;
      for (i=0;i<nvar;i++) {
	lambda+=p[i]*v[i];
      }
      if (lambda==0) {
	O2SCL_ERR2("Approximation to Jacobian has collapsed in ",
		   "mroot_broyden::iterate().",exc_ezerodiv);
      }
      
      /* v'=v+t*p */
      for (i=0;i<nvar;i++) {
	v[i]=v[i]+t*p[i];
      }
      
      /* w^T=p^T H */
      for (i=0;i<nvar;i++) {
	double sum=0;
	for (j=0;j<nvar;j++) {
	  sum+=H(j,i)*p[j];
	}
	w[i]=sum;
      }
      
      /* Hij -> Hij-(vi wj/lambda) */
      for (i=0;i<nvar;i++) {
	double vi=v[i];
	for (j=0;j<nvar;j++) {
	  double wj=w[j];
	  H(i,j)-=vi*wj/lambda;
	}
      }
      
      /* copy fnew into f */
      vector_copy(nvar,fnew,f);
      
      /* copy x_trial into x */
      vector_copy(nvar,x_trial,x);
      for (i=0;i<nvar;i++) {
	dx[i]=t*p[i];
      }
      
      phi=phi1;
      
      return success;
    }

    /// Desc
    virtual int msolve(size_t n, vec_t &x, func_t &func) {

      int status;

      set(func,n,x,f_int,dx_int);
      
      int iter=0;

      do {
	iter++;
	
	status=iterate();
	if (status) break;

	// ------------------------------------------------------
	// The equivalent of the statement:
	// 
	// status=gsl_multiroot_test_residual(f,this->tol_rel);

	double resid=0.0;
	for(size_t i=0;i<n;i++) {
	  resid+=fabs(f_int[i]);
	}
	if (resid<this->tol_rel) status=success;
	else status=gsl_continue;
	
	// ------------------------------------------------------
	
	if (this->verbose>0) {
	  this->print_iter(n,x,f_int,iter,resid,this->tol_rel,
			   "mroot_broyden");
	}
	
      } while (status==gsl_continue && iter<this->ntrial);

      this->last_ntrial=iter;
    
      if (iter>=this->ntrial) {
	O2SCL_CONV2_RET("Function mroot_broyden::msolve() ",
			"exceeded max. number of iterations.",
			exc_emaxiter,this->err_nonconv);
      }
      
      if (status!=0) {
	O2SCL_ERR2("Function iterate() failed in ",
		       "mroot_broyden::solve_set().",exc_efailed);
      }

      return success;
    }

#ifndef DOXYGEN_INTERNAL

 private:
  
  mroot_broyden<func_t,vec_t,mat_t,jfunc_t>
  (const mroot_broyden<func_t,vec_t,mat_t,jfunc_t> &);
  mroot_broyden<func_t,vec_t,mat_t,jfunc_t>& operator=
  (const mroot_broyden<func_t,vec_t,mat_t,jfunc_t>&);

#endif

  };

}

#endif
