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
/* roots/steffenson.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, 
 * Brian Gough
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

#ifndef O2SCL_ROOT_STEF_H
#define O2SCL_ROOT_STEF_H

/** \file root_stef.h
    \brief File defining \ref o2scl::root_stef 
*/

#include <string>

#include <o2scl/misc.h>
#include <o2scl/root.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Steffenson equation solver (GSL)

      This is Newton's method with an Aitken "delta-squared"
      acceleration of the iterates. This can improve the convergence
      on multiple roots where the ordinary Newton algorithm is slow.

      Defining the next iteration with
      \f[
      x_{i+1} = x_i - f(x_i) / f^{\prime}(x_i)
      \f]
      the accelerated value is
      \f[
      x_{\mathrm{acc},i} = x_i - (x_{i+1}-x_i)^2 / (x_{i+2} - 
      2 x_{i+1} + x_i)
      \f]
      We can only use the accelerated estimate after three iterations,
      and use the unaccelerated value until then.

      This class finds a root of a function a derivative. If the
      derivative is not analytically specified, it is most likely
      preferable to use of the alternatives, \ref o2scl::root_brent_gsl,
      \ref o2scl::root_bkt_cern, or \ref o2scl::root_cern. The function
      solve_de() performs the solution automatically, and a
      lower-level GSL-like interface with set() and iterate() is also
      provided.

      By default, this solver compares the present value of the root
      (\f$ \mathrm{root} \f$) to the previous value (\f$ \mathrm{x}
      \f$), and returns success if \f$ | \mathrm{root} - \mathrm{x} |
      < \mathrm{tol} \f$, where \f$ \mathrm{tol} = \mathrm{tol\_abs} +
      \mathrm{tol\_rel2}~\mathrm{root} \f$ .

      If \ref test_residual is set to true, then the solver
      additionally requires that the absolute value of the function is
      less than \ref root::tol_rel.

      The original variable \c x_2 has been removed as it was unused
      in the original GSL code.

      See the \ref onedsolve_subsect section of the User's guide for
      general information about \o2 solvers. 

      \future There's some extra copying here which can probably
      be removed.
      \future Compare directly to GSL.
      \future This can probably be modified to shorten the step
      if the function goes out of bounds as in exc_mroot_hybrids.
      
  */
  template<class func_t=funct, class dfunc_t=func_t> class root_stef : 
  public root_de<func_t,dfunc_t> {
    
  protected:
      
    /// Function value
    double f;

    /// Derivative value
    double df;

    /// Previous value of root
    double x_1;

    /// Root
    double x;

    /// Number of iterations
    int count;

    /// The function to solve
    func_t *fp;

    /// The derivative
    dfunc_t *dfp;
      
  public:
  
    root_stef() {
      test_residual=false;
      tol_rel2=1.0e-12;
    }
    
    /// Return the type, \c "root_stef".
    virtual const char *type() { return "root_stef"; }

    /// The present solution estimate
    double root;

    /** \brief The relative tolerance for subsequent solutions 
	(default \f$ 10^{-12} \f$)
    */
    double tol_rel2;

    /** \brief Perform an iteration
	
	After a successful iteration, \ref root contains the
	most recent value of the root. 
    */
    int iterate() {
      
      double x_new, f_new, df_new;
	
      double x_1t=x_1;
      double xt=x;
	
      if (df == 0.0) {
	O2SCL_CONV_RET("Derivative is zero in root_stef::iterate().",
		       o2scl::exc_ezerodiv,this->err_nonconv);
      }
	
      x_new=xt-(f/df);
      
      // It is important that the derivative be evaluated first here,
      // because we want the last function evaluation to be an
      // evaluation for the returned root
      df_new=(*dfp)(x_new);
      f_new=(*fp)(x_new);
    
      x_1=xt;
      x=x_new;
	
      f=f_new;
      df=df_new;
      
      if (!o2scl::is_finite(f_new)) {
	std::string str="Function not finite (returned "+dtos(f_new)+
	  ") in root_stef::iterate().";
	O2SCL_ERR(str.c_str(),o2scl::exc_ebadfunc);
      }
      
      if (count<3) {
	
	root=x_new;
	count++;
	
      } else {

	double u=(xt-x_1t);
	double v=(x_new-2*xt+x_1t);
	
	if (v == 0) {
	  // Avoid division by zero
	  root=x_new; 
	} else {
	  // Accelerated value
	  root=x_1t-u*u/v; 
	}
      }
      
      if (!o2scl::is_finite(df_new)) {
	std::string str="Derivative not finite (returned "+dtos(df_new)+
	  ") in root_stef::iterate().";
	O2SCL_ERR(str.c_str(),o2scl::exc_ebadfunc);
      }
      
      return o2scl::success;
    }

    /** \brief Solve \c func using \c x as an initial
	guess using derivatives \c df.
    */
    virtual int solve_de(double &xx, func_t &fun, dfunc_t &dfun) {
      
      int status1, status2=gsl_continue, iter=0;
	
      set(fun,dfun,xx);

      while (status1==success && status2==gsl_continue && 
	     iter<this->ntrial) {
	iter++;

	status1=iterate();

	// Compare present value to previous value
	status2=gsl_root_test_delta(root,xx,this->tol_abs,tol_rel2);

	if (test_residual && status2==success) {
	  double y;
	  y=fun(root);
	  if (fabs(y)>=this->tol_rel) status2=gsl_continue;
	}

	if (this->verbose>0) {
	  double fval;
	  fval=fun(root);
	  this->print_iter(root,fval,iter,fabs(root-xx),this->tol_abs*root,
			   "root_stef");
	}
	xx=root;
      }
      
      this->last_ntrial=iter;

      if (status1!=success || status2!=success) {
	int ret=o2scl::err_hnd->get_errno();
	return ret;
      }
      if (iter>=this->ntrial) {
	std::string str=((std::string)"Function solve_de() exceeded the ")+
	  "maximum number of iterations, "+itos(this->ntrial)+".";
	O2SCL_CONV_RET(str.c_str(),exc_emaxiter,this->err_nonconv);
      }

      return o2scl::success;
    }
    
    /// True if we should test the residual also (default false)
    bool test_residual;

    /** \brief Set the information for the solver

	Set the function, the derivative, the initial guess and
	the parameters.
    */
  void set(func_t &fun, dfunc_t &dfun, double guess) {
    
    fp=&fun;
    dfp=&dfun;
    root=guess;
    
    df=dfun(root);
    f=fun(root);
    
    x=root;
    x_1=0.0;
    count=1;
    
    return;
  }
  
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

