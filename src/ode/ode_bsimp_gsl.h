/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
/* ode-initval/bsimp.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA 02110-1301, USA.
 */

#ifndef O2SCL_GSL_BSIMP_H
#define O2SCL_GSL_BSIMP_H

/** \file ode_bsimp_gsl.h
    \brief File defining \ref o2scl::ode_bsimp_gsl 
*/

#include <gsl/gsl_math.h>

#include <o2scl/err_hnd.h>
#include <o2scl/ode_step.h>
#include <o2scl/ode_jac_funct.h>
#include <o2scl/lu.h>
#include <o2scl/permutation.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Bulirsch-Stoer implicit ODE stepper (GSL)
      
      Bader-Deuflhard implicit extrapolative stepper (\ref Bader83). 

      \note The variable <tt>h_next</tt> was defined in the original
      GSL version has been removed here, as it was unused by the
      stepper routine. 

      \note At the moment, this class retains the GSL approach to
      handling non-integer return values in the user-specified
      derivative function. If the user-specified function returns \c
      exc_efailed, then the stepper attempts to decrease the stepsize
      to continue. If the user-specified function returns a non-zero
      value other than \c exc_efailed, or if the Jacobian evaluation
      returns any non-zero value, then the stepper aborts and returns
      the error value without calling the error handler. This behavior
      may change in the future.

      There is an example for the usage of this class in
      <tt>examples/ex_stiff.cpp</tt> documented in the \ref
      ex_stiff_sect section.

      \future More detailed documentation about the algorithm

      \future Figure out if this should be a child of ode_step or
      astep. The function step_local() is actually its own ODE
      stepper and could be reimplemented as an object of type ode_step.

      \future I don't like setting yerr to GSL_POSINF, there should
      be a better way to force an adaptive stepper which is calling
      this stepper to readjust the stepsize.

      \future The functions deuf_kchoice() and compute_weights() can
      be moved out of this header file.

      \future Rework internal arrays

      \future Rework linear solver to be amenable to using
      a sparse matrix solver
  */
  template<class func_t=ode_funct, class jac_func_t=ode_jac_funct, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class ode_bsimp_gsl {

  public:

  typedef boost::numeric::ublas::unbounded_array<double> ubarr;
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

#ifndef DOXYGEN_INTERAL

  protected:

  /// Size of allocated vectors
  size_t dim;

  /// Function specifying derivatives
  func_t *funcp;

  /// Jacobian
  jac_func_t *jfuncp;

  /// Number of entries in the Bader-Deuflhard extrapolation sequence
  static const int sequence_count=8;

  /// Largest index of the Bader-Deuflhard extrapolation sequence
  static const int sequence_max=7;

  /// Workspace for extrapolation 
  ubmatrix d;

  /// Workspace for linear system matrix 
  ubmatrix a_mat;            

  /** \brief Workspace for extrapolation 
	
      (This state variable was named 'x' in GSL.)
  */
  double ex_wk[sequence_max];       

  /// \name State info
  //@{
  size_t k_current;
  size_t k_choice;
  double eps;
  //@}

  /// \name Workspace for extrapolation step 
  //@{
  ubarr yp;
  ubarr y_save;
  ubarr yerr_save;
  ubarr y_extrap_save;
  vec_t y_extrap_sequence;
  ubarr extrap_work;
  vec_t dfdt;
  vec_t y_temp;
  vec_t delta_temp;
  ubarr weight;
  mat_t dfdy;
  //@}
    
  /// \name Workspace for the basic stepper
  //@{
  vec_t rhs_temp;
  ubarr delta;
  //@}

  /// Order of last step 
  size_t order;

  /// Compute weights
  int compute_weights(const ubarr &y, ubarr &w, size_t n) {
    size_t i;
    // w[i] is 1 if y[i] is zero and the absolute value of y[i]
    // otherwise
    for (i=0;i<n;i++) {
      double u=fabs(y[i]);
      w[i]=(u>0.0) ? u : 1.0;
    }
    return 0;
  }
  
  /** \brief Calculate a choice for the "order" of the method, using the
      Deuflhard criteria. 

      Used in the allocate() function.
  */
  size_t deuf_kchoice(double eps2, size_t dimension) {

    const double safety_f=0.25;
    const double small_eps=safety_f*eps2;
      
    double a_work[sequence_count];
    double alpha[sequence_max][sequence_max];
    
    /* Bader-Deuflhard extrapolation sequence */
    static const int bd_sequence[sequence_count]=
    {2,6,10,14,22,34,50,70};
    
    int i, k;
    
    a_work[0]=bd_sequence[0]+1.0;
    
    for (k=0;k<sequence_max;k++) {
      a_work[k+1]=a_work[k]+bd_sequence[k+1];
    }

    for (i=0;i<sequence_max;i++) {
      alpha[i][i]=1.0;
      for (k=0;k<i;k++) {
	const double tmp1=a_work[k+1]-a_work[i+1];
	const double tmp2=(a_work[i+1]-a_work[0]+1.0)*(2*k+1);
	alpha[k][i]=pow(small_eps,tmp1/tmp2);
      }
    }

    a_work[0]+=dimension;

    for (k=0;k<sequence_max;k++) {
      a_work[k+1]=a_work[k]+bd_sequence[k+1];
    }

    for (k=0;k<sequence_max-1;k++) {
      if (a_work[k+2]>a_work[k+1]*alpha[k][k+1]) {
	break;
      }
    }

    return k;
  }

  /** \brief Polynomial extrapolation

      Compute the step of index <tt>i_step</tt> using polynomial
      extrapolation to evaulate functions by fitting a polynomial
      to estimates <tt>(x_i,y_i)</tt> and output the result to
      <tt>y_0</tt> and <tt>y_0_err</tt>. 

      The index <tt>i_step</tt> begins with zero.
	
  */
  int poly_extrap(ubmatrix &dloc, const double x[],
		  const unsigned int i_step, const double x_i,
		  const vec_t &y_i, vec_t &y_0, vec_t &y_0_err, 
		  ubarr &work) {
    size_t j, k;
      
    o2scl::vector_copy(dim,y_i,y_0_err);
    o2scl::vector_copy(dim,y_i,y_0);

    if (i_step == 0) {

      for (j=0;j<dim;j++) {
	dloc(0,j)=y_i[j];
      }

    } else {
	
      o2scl::vector_copy(dim,y_i,work);
	
      for (k=0;k<i_step;k++) {
	double deltaloc=1.0/(x[i_step-k-1]-x_i);
	const double f1=deltaloc*x_i;
	const double f2=deltaloc*x[i_step-k-1];
	  
	for (j=0;j<dim;j++) {
	  const double q_kj=dloc(k,j);
	  dloc(k,j)=y_0_err[j];
	  deltaloc=work[j]-q_kj;
	  y_0_err[j]=f1*deltaloc;
	  work[j]=f2*deltaloc;
	  y_0[j]+=y_0_err[j];
	}
      }
	
      for (j=0;j<dim;j++) {
	dloc(i_step,j)=y_0_err[j];
      }
    }
    return 0;
  }

  /** \brief Basic implicit Bulirsch-Stoer step
	
      Divide the step <tt>h_total</tt> into <tt>n_step</tt> smaller
      steps and do the Bader-Deuflhard semi-implicit iteration. This
      function starts at <tt>t0</tt> with function values
      <tt>y</tt>, derivatives <tt>yp_loc</tt>, and information from
      the Jacobian to compute the final value <tt>y_out</tt>.
  */
  int step_local(const double t0, const double h_total, 
		 const unsigned int n_step, const ubarr &y, 
		 const ubarr &yp_loc, const vec_t &dfdt_loc,
		 const mat_t &dfdy_loc, vec_t &y_out) {
      
    const double h=h_total/n_step;
    double t=t0+h;

    double sum;

    /* This is the factor sigma referred to in equation 3.4 of the
       paper.  A relative change in y exceeding sigma indicates a
       runaway behavior. According to the authors suitable values
       for sigma are >>1.  I have chosen a value of 100*dim. BJG
    */
    const double max_sum=100.0*dim;

    int signum, status;
    size_t i, j;
    size_t n_inter;

    /* Calculate the matrix for the linear system. */
    for (i=0;i<dim;i++) {
      for (j=0;j<dim;j++) {
	a_mat(i,j)=-h*dfdy_loc(i,j);
      }
      a_mat(i,i)=a_mat(i,i)+1.0;
    }

    /* LU decomposition for the linear system. */
      
    o2scl::permutation p_vec(dim);

    o2scl_linalg::LU_decomp(dim,a_mat,p_vec,signum);

    /* Compute weighting factors */
    compute_weights(y,weight,dim);

    /* Initial step. */

    for (i=0;i<dim;i++) {
      y_temp[i]=h*(yp_loc[i]+h*dfdt_loc[i]);
    }

    o2scl_linalg::LU_solve(dim,a_mat,p_vec,y_temp,delta_temp);
      
    sum=0.0;
    for (i=0;i<dim;i++) {
      const double di=delta_temp[i];
      delta[i]=di;
      y_temp[i]=y[i]+di;
      sum+=fabs(di)/weight[i];
    }
    if (sum>max_sum) {
      return exc_efailed;
    }

    /* Intermediate steps. */

    status=(*funcp)(t,dim,y_temp,y_out);
    if (status) {
      return status;
    }
      
    for (n_inter=1;n_inter<n_step;n_inter++) {

      for (i=0;i<dim;i++) {
	rhs_temp[i]=h*y_out[i]-delta[i];
      }
	
      o2scl_linalg::LU_solve(dim,a_mat,p_vec,rhs_temp,delta_temp);

      sum=0.0;
      for (i=0;i<dim;i++) {
	delta[i]+=2.0*delta_temp[i];
	y_temp[i]+=delta[i];
	sum+=fabs(delta[i])/weight[i];
      }
      if (sum>max_sum) {
	return exc_efailed;
      }
	
      t+=h;
	
      status=(*funcp)(t,dim,y_temp,y_out);
      if (status) {
	return status;
      }
    }


    /* Final step. */

    for (i=0;i<dim;i++) {
      rhs_temp[i]=h*y_out[i]-delta[i];
    }
      
    o2scl_linalg::LU_solve(dim,a_mat,p_vec,rhs_temp,delta_temp);
      
    sum=0.0;
    for (i=0;i<dim;i++) {
      y_out[i]=y_temp[i]+delta_temp[i];
      sum+=fabs(delta_temp[i])/weight[i];
    }

    if (sum>max_sum) {
      return exc_efailed;
    }
      
    return success;
  }

  /// Allocate memory for a system of size \c n
  int allocate(size_t n)  {

    if (dim>0) free();

    dim=n;
      
    d.resize(sequence_max,n);
    a_mat.resize(n,n);
      
    yp.resize(n);

    // AWS, 12/2/08 - This was added to ensure memory reallocation 
    // resets the stepper just like reset() does
    for(size_t i=0;i<n;i++) yp[i]=0.0;

    y_save.resize(n);
    yerr_save.resize(n);
    y_extrap_save.resize(n);
    extrap_work.resize(n);
    weight.resize(n);

    dfdy.resize(n,n);
    dfdt.resize(n);
    y_extrap_sequence.resize(n);
    y_temp.resize(n);
    rhs_temp.resize(n);
    delta_temp.resize(n);

    delta.resize(n);

    double sqrt_dbl_eps=sqrt(std::numeric_limits<double>::epsilon());

    // This choice of epsilon is not necessarily optimal, it has
    // a "FIXME" comment in the original GSL code
    size_t k_choice_loc=deuf_kchoice(sqrt_dbl_eps,n);
    k_choice=k_choice_loc;
    k_current=k_choice_loc;
    order=2*k_choice_loc;
      
    return 0;
  }

  /// Free allocated memory
  void free() {
    if (dim>0) {
      delta.resize(0);
      weight.resize(0);
      extrap_work.resize(0);
      y_extrap_save.resize(0);
      y_save.resize(0);
      yerr_save.resize(0);
      yp.resize(0);
      d.resize(0,0);
      dim=0;
    }
  }

#endif
    
  public:
    
  ode_bsimp_gsl() {
    dim=0;
    verbose=0;
  }
      
  virtual ~ode_bsimp_gsl() {
    if (dim>0) free();
  }

  /// Verbose parameter
  int verbose;
  
  /// Reset stepper
  int reset() {
    for(size_t i=0;i<dim;i++) yp[i]=0.0;
    return success;
  }
    
  /** \brief Perform an integration step
	
      Given initial value of the n-dimensional function in \c y and
      the derivative in \c dydx (which must generally be computed
      beforehand) at the point \c x, take a step of size \c h giving
      the result in \c yout, the uncertainty in \c yerr, and the new
      derivative in \c dydx_out (at \c x+h) using function \c derivs
      to calculate derivatives.  Implementations which do not
      calculate \c yerr and/or \c dydx_out do not reference these
      variables so that a blank \c vec_t can be given. All of the
      implementations allow \c yout=y and \c dydx_out=dydx if
      necessary.
  */
  virtual int step(double x, double h, size_t n, vec_t &y, vec_t &dydx, 
		   vec_t &yout, vec_t &yerr, vec_t &dydx_out, 
		   func_t &derivs, jac_func_t &jac) {

    int ret;

    funcp=&derivs;
    jfuncp=&jac;

    if (n!=dim) allocate(n);

    /* Bader-Deuflhard extrapolation sequence */
    static const int bd_sequence[sequence_count]={2,6,10,14,22,34,50,70};
    
    double t_local=x;

    size_t i, k;

    if (h+t_local==t_local) {
      // This section is labeled with a "FIXME" comment in the
      // original GSL code. I'm not sure why, but an error is
      // sensible here.
      O2SCL_ERR("Stepsize underflow in ode_bsimp_gsl::step().",
		    exc_eundrflw);
    }

    /* Save inputs */
    o2scl::vector_copy(dim,y,y_extrap_save);
    o2scl::vector_copy(dim,y,y_save);
    o2scl::vector_copy(dim,yerr,yerr_save);

    // Copy derivative
    o2scl::vector_copy(dim,dydx,yp);
      
    // Evaluate the Jacobian for the system. */
    ret=jac(t_local,dim,y,dfdy,dfdt);
    if (ret!=success) {
      return ret;
    }

    /* Make a series of refined extrapolations, up to the specified
       maximum order, which was calculated based on the Deuflhard
       criterion in the deuf_kchoice() function (which is called by
       allocate() ).
    */
    for (k=0;k<=k_current;k++) {

      const unsigned int N=bd_sequence[k];
      const double r=(h/N);
      const double x_k=r*r;

      // Each step computes a value of y_extrap_sequence,
      // using the number of sub-steps dictated by
      // the BD sequence
      int status=step_local(t_local,h,N,y_extrap_save,yp,
			    dfdt,dfdy,y_extrap_sequence);

      if (verbose>=2) {
	std::cout << "ode_bsimp_gsl: " << k << "/" << k_current << " "
	  << "status=" << status << std::endl;
      }

      if (status==exc_efailed) {
	/* If the local step fails, set the error to infinity in
	   order to force a reduction in the step size 
	*/
	for (i=0;i<dim;i++) {
	  yerr[i]=GSL_POSINF;
	}
	break;
      } else if (status!=success) {
	return status;
      }
	
      // Use the information in y_extrap_sequence to compute 
      // the new value of y and yerr .
      ex_wk[k]=x_k;
      poly_extrap(d,ex_wk,k,x_k,y_extrap_sequence,y,yerr,extrap_work);
    }

    /* Evaluate dydt_out[]. */

    ret=derivs(t_local+h,dim,y,dydx_out);
      
    // If we failed, copy the old values back to y and yerr
    if (ret!=success) {
      o2scl::vector_copy(dim,y_save,y);
      o2scl::vector_copy(dim,yerr_save,yerr);
      return ret;
    }

    return success;
  }
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
