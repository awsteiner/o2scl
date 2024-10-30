/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
/* multiroots/hybrid.c
 * multiroots/dogleg.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
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
/** \file mroot_hybrids.h
    \brief File defining \ref o2scl::mroot_hybrids and 
    specializations
*/
#ifndef O2SCL_MROOT_HYBRIDS_H
#define O2SCL_MROOT_HYBRIDS_H

#include <string>

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_linalg.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/vector.h>
#include <o2scl/mroot.h>
#include <o2scl/jacobian.h>
#include <o2scl/qr.h>
#include <o2scl/table.h>
// For matrix_out() below
#include <o2scl/columnify.h>

namespace o2scl {

  /** \brief Base functions for \ref mroot_hybrids
   */
  template<class fp_t=double>
  class mroot_hybrids_base {

  public:
    
    typedef boost::numeric::ublas::vector<fp_t> ubvector;
    typedef boost::numeric::ublas::matrix<fp_t> ubmatrix;

    /// Compute the actual reduction
    fp_t compute_actual_reduction(fp_t fnorm0, fp_t fnorm1) {
      fp_t actred;
      if (fnorm1<fnorm0) {
	fp_t u=fnorm1/fnorm0;
	actred=1-u*u;
      } else {
	actred=-1;
      }
      return actred;
    }

    /// Compute the predicted reduction phi1p=|Q^T f + R dx| 
    fp_t compute_predicted_reduction(fp_t fnorm0, fp_t fnorm1) {
      fp_t prered;
      if (fnorm1<fnorm0) {
	fp_t u=fnorm1/fnorm0;
	prered=1-u*u;
      } else {
	prered=0;
      }
      return prered;
    }

    /// Compute \f$ R \cdot g \f$ where \c g is the \gradient
    template<class vec2_t, class mat_t>
    void compute_Rg(size_t N, const mat_t &r2, 
                    const ubvector &gradient2, vec2_t &Rg) {
    
      for (size_t i=0;i<N;i++) {
	fp_t sum=0;
	for (size_t j=i;j<N;j++) {
	  sum+=r2(i,j)*gradient2[j];
	}
	Rg[i]=sum;
      }

      return;
    }

    /// Compute \c w and \c v
    template<class vec2_t>
    void compute_wv(size_t n, const ubvector &qtdf2, 
                    const ubvector &rdx2, const vec2_t &dx2, 
                    const ubvector &diag2, fp_t pnorm, 
                    ubvector &w2, ubvector &v2) {
    
      for (size_t i=0;i<n;i++) {
	fp_t diagi=diag2[i];
	w2[i]=(qtdf2[i]-rdx2[i])/pnorm;
	v2[i]=diagi*diagi*dx2[i]/pnorm;
      }
    
      return;
    }

    /// Compute \f$ R \cdot \mathrm{dx} \f$
    template<class vec2_t, class mat_t>
    void compute_rdx(size_t N, const mat_t &r2, 
                     const vec2_t &dx2, ubvector &rdx2) {
    
      for (size_t i=0;i<N;i++) {
	fp_t sum=0;
	for (size_t j=i;j<N;j++) {
	  sum+=r2(i,j)*dx2[j];
	}
	rdx2[i]=sum;
      }

      return;
    }

    /** \brief Compute the norm of the vector \f$ \vec{v} \f$ defined
	by \f$ v_i = d_i ff_i \f$
    */
    template<class vec2_t, class vec3_t>
    fp_t scaled_enorm(size_t n, const vec2_t &d, const vec3_t &ff) {
      fp_t e2=0;
      for (size_t i=0;i<n;i++) {
	fp_t u=d[i]*ff[i];
	e2+=u*u;
      }
      return sqrt(e2);
    }
  
    /// Compute delta
    template<class vec2_t>
    fp_t compute_delta(size_t n, ubvector &diag2, vec2_t &x2) {
      fp_t Dx=scaled_enorm(n,diag2,x2);
      fp_t factor=100;
      return (Dx > 0) ? factor*Dx : factor;
    }

    /** \brief Compute the Euclidean norm of \c ff
      
        \verbatim embed:rst      
        
        .. todo::
        
        class mroot_hybrids
        
        Future: Replace this with \c dnrm2 from \ref cblas_base.h
        
        \endverbatim
    */
    template<class vec2_t> fp_t enorm(size_t N, const vec2_t &ff) {
      fp_t e2=0;
      for (size_t i=0;i<N;i++) {
	fp_t fi=ff[i];
	e2+=fi*fi;
      }
      return sqrt(e2);
    }

    /// Compute the Euclidean norm of the sum of \c a and \c b
    fp_t enorm_sum(size_t n, const ubvector &a, 
		     const ubvector &b) {
      fp_t e2=0;
      for (size_t i=0;i<n;i++) {
	fp_t u=a[i]+b[i];
	e2+=u*u;
      }
      return sqrt(e2);
    }

    /** \brief Compute a trial step and put the result in \c xx_trial
	
        \verbatim embed:rst      
        
        .. todo::
        
        class mroot_hybrids
        
        Future: Replace this function with daxpy?
        
        \endverbatim
    */
    template<class vec2_t>
    int compute_trial_step(size_t N, vec2_t &xl, vec2_t &dxl, 
                           vec2_t &xx_trial) {
      for (size_t i=0;i<N;i++) {
	xx_trial[i]=xl[i]+dxl[i];
      }
      return 0;
    }

    /** \brief Compute the change in the function value
	
        \verbatim embed:rst      
        
        .. todo::
        
        class mroot_hybrids
        
        Future: Replace this function with daxpy?
        
        \endverbatim
    */
    template<class vec2_t>
    int compute_df(size_t n, const vec2_t &ff_trial, 
                   const vec2_t &fl, ubvector &dfl) {
      for (size_t i=0;i<n;i++) {
	dfl[i]=ff_trial[i]-fl[i];
      }
    
      return 0;
    }
   
    /// Compute \c diag, the norm of the columns of the Jacobian
    template<class mat2_t>
    void compute_diag(size_t n, const mat2_t &J2, ubvector &diag2) {
      for (size_t j=0;j<n;j++) {
	fp_t sum=0;
	for (size_t i=0;i<n;i++) {
	  const fp_t Jij=J2(i,j);
	  sum+=Jij*Jij;
	}
	if (sum == 0) {
	  sum=1;
	}
	diag2[j]=sqrt(sum);
      }

      return;
    }
  
    /** \brief Compute \f$ Q^{T} f \f$

        \verbatim embed:rst      
        
        .. todo::
        
        class mroot_hybrids
        
        Future: This function is just right-multiplication, so we
        could use the O2scl cblas routines instead.
        
        \endverbatim
    */
    template<class vec2_t, class vec3_t, class vec4_t>
    void compute_qtf(size_t N, const vec2_t &q2, 
                     const vec3_t &f2, vec4_t &qtf2) {
      for (size_t j=0;j<N;j++) {
	fp_t sum=0;
	for (size_t i=0;i<N;i++) {
	  sum+=q2(i,j)*f2[i];
	}
	qtf2[j]=sum;
      }

      return;
    }
  
    /// Update \c diag 
    template<class mat2_t>
    void update_diag(size_t n, const mat2_t &J2, ubvector &diag2) {
      for (size_t j=0;j<n;j++) {
	fp_t cnorm, diagj, sum=0;
	for (size_t i=0;i<n;i++) {
	  fp_t Jij=J2(i,j);
	  sum+=Jij*Jij;
	}
	if (sum == 0) {
	  sum=1;
	}
	cnorm=sqrt(sum);
	diagj=diag2[j];
	if (cnorm > diagj) {
	  diag2[j]=cnorm;
	}
      }

      return;
    }
    
    /** \brief Form appropriate convex combination of the Gauss-Newton
	direction and the scaled gradient direction.

	Using the Gauss-Newton direction given in \c newton (a vector of
	size N), and the gradient direction in \c gradient (also a
	vector of size N), this computes
	\f[
	\mathrm{pp}=\alpha \mathrm{newton}+\beta \mathrm{gradient}
	\f]
    */
    template<class vec2_t>
    void scaled_addition(size_t N, fp_t alpha, ubvector &newton2, 
                         fp_t beta, ubvector &gradient2, vec2_t &pp) {
      for (size_t i=0;i<N;i++) {
	pp[i]=alpha*newton2[i]+beta*gradient2[i];
      }

      return;
    }

    /// Compute the Gauss-Newton direction
    template<class mat_t>
    int newton_direction(const size_t N, const mat_t &r2, 
                         const ubvector &qtf2, ubvector &p) {
      size_t i;
      int status=0;

      vector_copy(N,qtf2,p);
      o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
			 o2scl_cblas::o2cblas_Upper,
			 o2scl_cblas::o2cblas_NoTrans,
			 o2scl_cblas::o2cblas_NonUnit,N,N,r2,p);
      
      for (i=0;i<N;i++) {
	p[i]*=-1;
      }

      return status;
    }

    /// Compute the gradient direction
    template<class mat_t>
    void gradient_direction(const size_t M, const size_t N,
                            const mat_t &r2, const ubvector &qtf2,
                            const ubvector &diag2, ubvector &g) {
      for (size_t j=0;j<M;j++) {
	fp_t sum=0;
	for (size_t i=0;i<N;i++) {
	  sum+=r2(i,j)*qtf2[i];
	}
	g[j]=-sum/diag2[j];
      }

      return;
    }

    /// Compute the point at which the gradient is minimized
    void minimum_step(const size_t N, fp_t gnorm, 
		      const ubvector &diag2, ubvector &g) {
      for (size_t i=0;i<N;i++) {
	g[i]=(g[i]/gnorm)/diag2[i];
      }

      return;
    }
  
    /** \brief Take a dogleg step
      
	Given the QR decomposition of an \c n by \c n matrix "A",
	a diagonal matrix \c diag, a vector \c b, and a positive
	number \c delta, this function determines the convex 
	combination \c x of the Gauss-Newton and scaled gradient
	directions which minimizes \f$ A x-b \f$ in the least
	squares sense, subject to the restriction that the 
	Euclidean norm of \f$ d x \f$ is smaller than the 
	value of \c delta.
    */
    template<class vec2_t, class mat_t>
    int dogleg(size_t n, const mat_t &r2, const ubvector &qtf2,
               const ubvector &diag2, fp_t delta2,
               ubvector &newton2, ubvector &gradient2, 
               vec2_t &p) {

      fp_t qnorm, gnorm, sgnorm, bnorm, temp;
      
      // Compute the Gauss-Newton direction
      newton_direction(n,r2,qtf2,newton2);
      qnorm=scaled_enorm(n,diag2,newton2);

      // Test whether the gauss-newton direction is acceptable
      if (qnorm <= delta2) {
	for(size_t i=0;i<n;i++) p[i]=newton2[i];
	return success;
      }
      
      // Compute the gradient direction
      gradient_direction(n,n,r2,qtf2,diag2,gradient2);
      gnorm=enorm(n,gradient2);
      
      // The special case that the scaled gradient is zero
      if (gnorm == 0) {
	fp_t alpha=delta2/qnorm;
	fp_t beta=0;
	scaled_addition(n,alpha,newton2,beta,gradient2,p);
	return success;
      }

      // -----------------------------------------------------
      // The next set of statements calculates the point along the
      // scaled gradient at which the quadratic is minimized.

      minimum_step(n,gnorm,diag2,gradient2);
      
      // Use p as temporary space to compute Rg
      compute_Rg(n,r2,gradient2,p); 
      
      temp=enorm(n,p);
      sgnorm=(gnorm/temp)/temp;

      // -----------------------------------------------------
      // Test whether the scaled gradient direction is acceptable

      if (sgnorm>delta2)  {

	fp_t alpha=0;
	fp_t beta=delta2;  

	scaled_addition(n,alpha,newton2,beta,gradient2,p);

	return success;
      }
      
      // If not, calculate the point along the dogleg step
      // at which the quadratic is minimized
      bnorm=enorm(n,qtf2);
      
      fp_t bg=bnorm/gnorm;
      fp_t bq=bnorm/qnorm;
      fp_t dq=delta2/qnorm;
      fp_t dq2=dq*dq;
      fp_t sd=sgnorm/delta2;
      fp_t sd2=sd*sd;
      
      fp_t t1=bg*bq*sd;
      fp_t u=t1-dq;
      fp_t t2=t1-dq*sd2+sqrt(u*u+(1-dq2)*(1-sd2));
      
      fp_t alpha=dq*(1-sd2)/t2;
      fp_t beta=(1-alpha)*sgnorm;

      // Form the appropriate convex combination of the gauss-newton
      // direction and the scaled gradient direction.

      scaled_addition(n,alpha,newton2,beta,gradient2,p);
      
      return success;
    }
    
  };

  /** \brief Multidimensional root-finding algorithm using
      Powell's Hybrid method (GSL)
      
      \verbatim embed:rst
      This is a recasted version of the GSL routines which use a
      modified version of Powell's Hybrid method as implemented in the
      HYBRJ algorithm in MINPACK [Garbow80]_.  
      \endverbatim

      Both the scaled and
      unscaled options are available by setting \ref int_scaling (the
      scaled version is the default). If derivatives are not provided,
      they will be computed automatically. This class provides the
      GSL-like interface using allocate(), set() (or set_de() in case
      where derivatives are available), iterate(), and the
      higher-level interfaces, msolve() and msolve_de(), which perform
      the solution and the memory allocation automatically. Some
      additional checking is performed in case the user calls the
      functions out of order (i.e. set() without allocate()).
      
      The functions msolve() and msolve_de() use the condition \f$
      \sum_i |f_i|<\f$ \ref mroot::tol_rel to determine if the solver has
      succeeded.
      
      \verbatim embed:rst
      See the :ref:`Multi-dimensional solvers` section of the User's
      guide for general information about the O2scl solvers. There is
      an example for the usage of the multidimensional solver classes
      given in ``examples/ex_mroot.cpp``, see the
      :ref:`Multi-dimensional solver example`.
      \endverbatim

      \note The \ref set() and \ref set_de() functions store a pointer
      to the function object and the user must ensure that the object
      is still valid for a later call to \ref iterate(). 

      The original GSL algorithm has been modified to shrink the
      stepsize if a proposed step causes the function to return a
      non-zero value. This allows the routine to automatically try to
      avoid regions where the function is not defined. The algorithm
      is also modified to check that it is not sending non-finite
      values to the user-specified function. To return to the default
      GSL behavior, set \ref shrink_step and \ref extra_finite_check
      to false.

      The default method for numerically computing the Jacobian is
      from \ref jacobian_gsl. This default is identical to the GSL
      approach, except that the default value of \ref
      jacobian_gsl::epsmin is non-zero. See \ref jacobian_gsl
      for more details.

      By default convergence failures result in calling the exception
      handler, but this can be turned off by setting \ref
      mroot::err_nonconv to false. If \ref mroot::err_nonconv is
      false, then the functions \ref iterate(), \ref msolve() and \ref
      msolve_de() will return a non-zero value if convergence fails.
      Note that the default Jacobian object, \ref def_jac also has a
      data member \ref jacobian_gsl::err_nonconv which separately
      handles the case where the one row of the Jacobian is all zero.

      \verbatim embed:rst      
        
      .. todo::
       
      class mroot_hybrids
        
      Future: 
        
      - Is all the setting of vectors and matrices to zero really
      necessary? Do they need to be executed even if memory hasn't
      been recently allocated?
      - Convert more ubvectors to vec_t.
      - Some more of the element-wise vector manipulation could be
      converted to BLAS routines.
      - It's kind of strange that set() sets jac_given to false
      and set_de() has to reset it to true. Can this be simplified?
      - Many of these minpack functions could be put in their
      own "minpack_tools" class, or possibly moved to be
      linear algebra routines instead.
      - There are still some numbers in here which the user
      could have control over, for example, the ``nslow2``
      threshold which indicates failure.
         
      \endverbatim
  */
  template<
    class func_t=mm_funct,
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double>,
    class jfunc_t=jac_funct > class mroot_hybrids : 
    public mroot<func_t,vec_t,jfunc_t>, mroot_hybrids_base<double> {
    
  protected:

    /// Number of iterations
    int iter;
    /// Compute the number of failures
    size_t ncfail;
    /// Compute the number of successes
    size_t ncsuc;
    /// The number of times the actual reduction is less than 0.001
    size_t nslow1;
    /// The number of times the actual reduction is less than 0.1
    size_t nslow2;
    /// The norm of the current function value
    double fnorm;
    /// The limit of the Nuclidean norm
    double delta;
    /// Jacobian
    mat_t J;
    /// Q matrix from QR decomposition
    mat_t q;
    /// R matrix from QR decomposition
    mat_t r;
    /// The diagonal elements
    ubvector diag;
    /// The value of \f$ Q^T f \f$
    ubvector qtf;
    /// The Newton direction
    ubvector newton;
    /// The gradient direction
    ubvector gradient;
    /// The change in the function value
    ubvector df;
    /// The value of \f$ Q^T \cdot \mathrm{df} \f$
    ubvector qtdf;
    /// The value of \f$ R \cdot \mathrm{dx} \f$
    ubvector rdx;
    /// The value of \f$ w=(Q^T df-R dx)/|dx| \f$
    ubvector w;
    /// The value of \f$ v=D^2 dx/|dx| \f$
    ubvector v;
  
    /// The user-specified Jacobian
    jfunc_t *jac;
    
    /// The automatic Jacobian
    jacobian<func_t,vec_t,mat_t> *ajac;

    /// The value of the derivative
    vec_t dx;

    /// Trial root
    vec_t x_trial;
    
    /// Trial function value
    vec_t f_trial;

    /// The number of equations and unknowns
    size_t dim;

    /// True if the jacobian has been given
    bool jac_given;

    /// The user-specified function
    func_t *fnewp;

    /// True if "set" has been called
    bool set_called;

    /// Finish the solution after set() or set_de() has been called
    virtual int solve_set(size_t nn, vec_t &xx, func_t &ufunc) {

      int status;
      iter=0;

      do {
        iter++;
	
        if (iterate()!=0) {
          O2SCL_CONV2_RET("Function iterate() failed in mroot_hybrids::",
                          "solve_set().",exc_efailed,this->err_nonconv);
        }

        // ------------------------------------------------------
        // The equivalent of the statement:
        // 
        // status=gsl_multiroot_test_residual(f,this->tol_rel);

        double resid=0.0;
        for(size_t i=0;i<nn;i++) {
          resid+=fabs(f[i]);
        }
        if (resid<this->tol_rel) status=success;
        else status=gsl_continue;
	
        // ------------------------------------------------------
	
        if (this->verbose>0) {
          this->print_iter(nn,x,f,iter,resid,this->tol_rel,
                           "mroot_hybrids");
        }
	
      } while (status==gsl_continue && iter<this->ntrial);

      for(size_t i=0;i<nn;i++) {
        xx[i]=x[i];
      }

      this->last_ntrial=iter;
    
      if (((int)iter)>=this->ntrial) {
        O2SCL_CONV2_RET("Function mroot_hybrids::msolve() ",
                        "exceeded max. number of iterations.",
                        exc_emaxiter,this->err_nonconv);
      }
    
      return success;
    }

  public:
      
    mroot_hybrids() {
      shrink_step=true;
      dim=0;
      ajac=&def_jac;
      //def_jac.set_epsrel(sqrt(std::numeric_limits<double>::epsilon()));
      int_scaling=true;
      jac_given=false;
      set_called=false;
      extra_finite_check=true;
      store_funcs=false;
    }
  
    virtual ~mroot_hybrids() {
    }

    /** \brief If true, iterate() will shrink the step-size automatically if
        the function returns a non-zero value (default true)

        The original GSL behavior can be obtained by setting 
        this to \c false.
    */
    bool shrink_step;

    /** \brief If true, double check that the input function values are
        finite (default true)
    */
    bool extra_finite_check;

    /// If true, use the internal scaling method (default true)
    bool int_scaling;

    /// Default automatic Jacobian object
    jacobian_gsl<func_t,vec_t,mat_t> def_jac;

    /// If true, store function evaluations
    bool store_funcs;
    
    /// Set the automatic Jacobian object
    virtual int set_jacobian(jacobian<func_t,vec_t,mat_t> &j) {
      ajac=&j;
      return 0;
    }
    
    /** \brief The value of the function at the present iteration
	
        \comment
        We need this to be public so that the user can see if 
        iterate() is working
        \endcomment
    */
    vec_t f;

    /// The present solution
    vec_t x;

    /** \brief Perform an iteration
	
        At the end of the iteration, the current value of the solution 
        is stored in \ref x.
    */
    int iterate() {
	
      if (!set_called) {
        O2SCL_ERR2("Function set() not called or most recent call ",
                   "failed in mroot_hybrids::iterate().",
                   exc_efailed);
      }

      double prered, actred;
      double pnorm, fnorm1, fnorm1p;
      double ratio;
      double p1=0.1, p5=0.5, p001=0.001, p0001=0.0001;
    
      /* Compute qtf=Q^T f */
      
      compute_qtf(dim,q,f,qtf);
  
      /* Compute dogleg step */
      
      dogleg(dim,r,qtf,diag,delta,newton,gradient,dx);

      /* Take a trial step */
      
      compute_trial_step(dim,x,dx,x_trial);
      pnorm=scaled_enorm(dim,diag,dx);

      if (iter==1) {
        if (pnorm<delta) {
          delta=pnorm;
        }
      } 

      if (extra_finite_check) {
        // AWS, 11/14/13: This is not strictly necessary, because if 
        // x_trial is not finite the iteration will fail to converge
        // anyway, but it is disconcerting to the user to have 
        // non-finite values sent to the user-specified function for 
        // no apparent reason. In reality what appears to be happening 
        // is that pnorm was previously zero (because of a vanishingly 
        // small step size), and the call to compute_wv() below then 
        // leads to non-finite values. On the other hand, checking all 
        // the vector values is time consuming, so I perform this 
        // check only if extra_finite_check is true.
        for(size_t ik=0;ik<dim;ik++) {
          if (!std::isfinite(x_trial[ik])) {
            O2SCL_CONV2_RET("Iteration lead to non-finite values in ",
                            "mroot_hybrids::iterate().",exc_efailed,
                            this->err_nonconv);
          }
        }
      }

      /* Evaluate function at x+p */
      
      int status;
    
      if (shrink_step==false) {
	  
        // Evaluate the function at the new point, exit if it fails
	
        status=(*fnewp)(dim,x_trial,f_trial);
	  
        if (status != success) {
          std::string str="Function returned non-zero value ("+itos(status)+
            ") in mroot_hybrids::iterate().";
          O2SCL_CONV_RET(str.c_str(),o2scl::exc_ebadfunc,this->err_nonconv);
        }
      
      } else {
	  
        // Evaluate the function at the new point, try to recover
        // if it fails
	  
        status=(*fnewp)(dim,x_trial,f_trial);
	  
        int bit=0;
        while(status!=0 && bit<20) {
          for(size_t ib=0;ib<dim;ib++) {
            x_trial[ib]=(x_trial[ib]+x[ib])/2.0;
          }
          status=(*fnewp)(dim,x_trial,f_trial);
          bit++;
        }
    
        // Exit if we weren't able to find a new good point

        if (status != success) {
          std::string str="No suitable step found, function returned ("+
            itos(status)+") in mroot_hybrids::iterate().";
          O2SCL_CONV_RET(str.c_str(),o2scl::exc_ebadfunc,this->err_nonconv);
        }
    
      }

      /* Set df=f_trial-f */

      compute_df(dim,f_trial,f,df);
      
      /* Compute the scaled actual reduction */
  
      fnorm1=enorm(dim,f_trial);
      actred=compute_actual_reduction(fnorm,fnorm1);
      
      /* Compute rdx=R dx */
  
      compute_rdx(dim,r,dx,rdx);

      /* Compute the scaled predicted reduction phi1p=|Q^T f+R dx| */
  
      fnorm1p=enorm_sum(dim,qtf,rdx);
      prered=compute_predicted_reduction(fnorm,fnorm1p);

      /* Compute the ratio of the actual to predicted reduction */
  
      if (prered > 0) {
        ratio=actred/prered;
      } else {
        ratio=0;
      }
  
      /* Update the step bound */
  
      if (ratio<p1) {
        ncsuc=0;
        ncfail++;
        delta*=p5;
      } else {
        ncfail=0;
        ncsuc++;
      
        if (ratio >= p5 || ncsuc > 1) {
          delta=GSL_MAX(delta,pnorm/p5);
        }
        if (fabs (ratio-1) <= p1) {
          delta=pnorm/p5;
        }
      }
  
      /* Test for successful iteration */

      if (ratio >= p0001) {
        for(size_t i=0;i<dim;i++) {
          x[i]=x_trial[i];
          f[i]=f_trial[i];
        }
        fnorm=fnorm1;
        iter++;
      }
      
      /* Determine the progress of the iteration */
  
      nslow1++;
      if (actred >= p001) {
        nslow1=0;
      }
  
      if (actred >= p1) {
        nslow2=0;
      }
    
      if (ncfail==2) {

        int jac_ret;
      
        if (jac_given) jac_ret=(*jac)(dim,x,dim,f,J);
        else jac_ret=(*ajac)(dim,x,dim,f,J);
      
        if (jac_ret!=0) {
          std::string str="Jacobian failed and returned ("+
            itos(jac_ret)+") in mroot_hybrids::iterate().";
          O2SCL_CONV_RET(str.c_str(),exc_efailed,this->err_nonconv);
        }
      
        nslow2++;
      
        if (iter==1) {
          if (int_scaling) {
            compute_diag(dim,J,diag);
          }
          delta=compute_delta(dim,diag,x);
        } else {
          if (int_scaling) {
            update_diag(dim,J,diag);
          }
        }
      
        o2scl_linalg::QR_decomp_unpack<mat_t,mat_t,mat_t,double>
          (dim,dim,this->J,this->q,this->r);

        return success;
      }
  
      /* Compute qtdf=Q^T df, w=(Q^T df-R dx)/|dx|, v=D^2 dx/|dx| */
  
      compute_qtf(dim,q,df,qtdf);
      compute_wv(dim,qtdf,rdx,dx,diag,pnorm,w,v);

      /* Rank-1 update of the jacobian Q'R'=Q(R+w v^T) */
    
      o2scl_linalg::QR_update<mat_t,mat_t,ubvector,
        ubvector,double>(dim,dim,q,r,w,v);

      /* No progress as measured by jacobian evaluations */

      if (nslow2==5) {
        O2SCL_CONV2_RET("No progress in Jacobian in mroot_hybrids::",
                        "iterate().",exc_enoprogj,this->err_nonconv);
      }
  
      /* No progress as measured by function evaluations */

      if (nslow1==10) {
        O2SCL_CONV2_RET("No progress in function in mroot_hybrids::",
                        "iterate().",exc_enoprog,this->err_nonconv);
      }
      
      return success;
    }
      
    /// Allocate the memory
    void allocate(size_t n) {

      q.resize(n,n);
      r.resize(n,n);
      diag.resize(n);
      qtf.resize(n);
      newton.resize(n);
      gradient.resize(n);
      df.resize(n);
      qtdf.resize(n);
      rdx.resize(n);
      w.resize(n);
      v.resize(n);

      for(size_t ii=0;ii<n;ii++) {
        for(size_t jj=0;jj<n;jj++) {
          q(ii,jj)=0.0;
        }
      }
      for(size_t ii=0;ii<n;ii++) {
        for(size_t jj=0;jj<n;jj++) {
          r(ii,jj)=0.0;
        }
      }

      for(size_t ii=0;ii<diag.size();ii++) diag[ii]=0.0;
      for(size_t ii=0;ii<qtf.size();ii++) qtf[ii]=0.0;
      for(size_t ii=0;ii<newton.size();ii++) newton[ii]=0.0;
      for(size_t ii=0;ii<gradient.size();ii++) gradient[ii]=0.0;
      for(size_t ii=0;ii<df.size();ii++) df[ii]=0.0;
      for(size_t ii=0;ii<qtdf.size();ii++) qtdf[ii]=0.0;
      for(size_t ii=0;ii<rdx.size();ii++) rdx[ii]=0.0;
      for(size_t ii=0;ii<w.size();ii++) w[ii]=0.0;
      for(size_t ii=0;ii<v.size();ii++) v[ii]=0.0;

      J.resize(n,n);

      x.resize(n);
      dx.resize(n);
      f.resize(n);
      x_trial.resize(n);
      f_trial.resize(n);

      for(size_t i=0;i<n;i++) {
        x[i]=0.0;
        dx[i]=0.0;
        f[i]=0.0;
      }
      
      dim=n;

      return;
    }
    
    /// Return the type,\c "mroot_hybrids".
    virtual const char *type() { return "mroot_hybrids"; }
    
    /** \brief Solve \c func with derivatives \c dfunc using \c x as 
        an initial guess, returning \c x.

    */
    virtual int msolve_de(size_t nn, vec_t &xx, func_t &ufunc,
                          jfunc_t &dfunc) {

      int ret=set_de(nn,xx,ufunc,dfunc);
      if (ret!=success) {
        O2SCL_CONV2_RET("Function set_de() failed in mroot_hybrids::",
                        "msolve_de().",exc_efailed,this->err_nonconv);
      }
      
      return solve_set(nn,xx,ufunc);
    }

    /// Solve \c ufunc using \c xx as an initial guess, returning \c xx.
    virtual int msolve(size_t nn, vec_t &xx, func_t &ufunc) {
      
      int ret=set(nn,xx,ufunc);
      if (ret!=success) {
        O2SCL_CONV2_RET("Function set() failed in mroot_hybrids::",
                        "msolve().",exc_efailed,this->err_nonconv);
      }
      int status=solve_set(nn,xx,ufunc);
      for(size_t i=0;i<nn;i++) xx[i]=x[i];

      return status;
    }
      
    /** \brief Set the function, the parameters, and the initial guess 
     */
    int set(size_t nn, vec_t &ax, func_t &ufunc) {

      int status;
  
      if (nn!=dim) { 
        allocate(nn);
      }
      
      fnewp=&ufunc;

      // Specify function for automatic Jacobian
      ajac->set_function(ufunc);
      
      if (dim==0) {
        O2SCL_ERR2("No memory allocated in ",
                   "mroot_hybrids::set().",o2scl::exc_ebadlen);
      }
    
      // Copy the user-specified solution
      for(size_t i=0;i<dim;i++) x[i]=ax[i];
      
      status=ufunc(dim,ax,f);
      if (status!=0) {
        O2SCL_CONV2_RET("Function returned non-zero in ",
                        "mroot_hybrids::set().",exc_ebadfunc,this->err_nonconv);
      }
    
      if (jac_given) status=(*jac)(dim,ax,dim,f,J);
      else status=(*ajac)(dim,ax,dim,f,J);

      if (status!=0) {
        O2SCL_CONV2_RET("Jacobian failed in ",
                        "mroot_hybrids::set().",exc_efailed,this->err_nonconv);
      }

      iter=1;
      fnorm=enorm(dim,f);
      ncfail=0;
      ncsuc=0;
      nslow1=0;
      nslow2=0;
      
      for(size_t i=0;i<dim;i++) dx[i]=0.0;
      
      /* Store column norms in diag */
  
      if (int_scaling) {
        compute_diag(dim,J,diag);
      } else {
        for(size_t ii=0;ii<diag.size();ii++) diag[ii]=0.0;
      }
	
      /* Set delta to factor |D x| or to factor if |D x| is zero */
	
      delta=compute_delta(dim,diag,x);
  
      /* Factorize J into QR decomposition */
      o2scl_linalg::QR_decomp_unpack<mat_t,mat_t,mat_t,double>
        (dim,dim,this->J,this->q,this->r);
      set_called=true;
      jac_given=false;

      return 0;
    }

    /** \brief Set the function, the Jacobian, the parameters,
        and the initial guess 
    */
    int set_de(size_t nn, vec_t &ax, func_t &ufunc, jfunc_t &dfunc) {

      fnewp=&ufunc;
      jac=&dfunc;

      // Make sure set() uses the right Jacobian
      jac_given=true;
    
      int ret=set(nn,ax,ufunc);
    
      // Reset jac_given since set() will set it back to false
      jac_given=true;
      set_called=true;
    
      return ret;
    }

  private:
  
    mroot_hybrids<func_t,vec_t,mat_t,jfunc_t>
    (const mroot_hybrids<func_t,vec_t,mat_t,jfunc_t> &);
    mroot_hybrids<func_t,vec_t,mat_t,jfunc_t>& operator=
    (const mroot_hybrids<func_t,vec_t,mat_t,jfunc_t>&);

  };

}

#if defined (O2SCL_COND_FLAG) || defined (DOXYGEN)
#if defined (O2SCL_SET_ARMA) || defined (DOXYGEN)
#include <armadillo>
namespace o2scl {
  /** \brief A version of \ref mroot_hybrids which uses 
      Armadillo for the QR decomposition

      \note This class template is only defined if Armadillo
      was enabled when \o2 was installed.
  */
  template<class func_t, class vec_t, class mat_t, class jfunc_t>
  class mroot_hybrids_arma_qr_econ :
    public mroot_hybrids<func_t,vec_t,mat_t,jfunc_t> {

    virtual void qr_decomp_unpack() {
      arma::qr_econ(this->q,this->r,this->J);
      return;
    }
  
  };
}
#endif
#if defined (O2SCL_HAVE_EIGEN) || defined (DOXYGEN)
#include <eigen3/Eigen/Dense>
namespace o2scl {
  /** \brief A version of \ref mroot_hybrids
      which uses Eigen for the QR decomposition

      \note This class template is only defined if Eigen 
      was enabled when \o2 was installed.
  */
  template<class func_t, class vec_t, class mat_t, class jfunc_t>
  class mroot_hybrids_eigen :
    public mroot_hybrids<func_t,vec_t,mat_t,jfunc_t> {
  
    virtual void qr_decomp_unpack() {
      Eigen::HouseholderQR<Eigen::MatrixXd> hqr(this->J);
      this->q=hqr.householderQ();
      this->r=hqr.matrixQR().triangularView<Eigen::Upper>();
      return;
    }
  
  };
}
#endif
#else
#include <o2scl/mroot_special.h>
#endif
  
#endif
