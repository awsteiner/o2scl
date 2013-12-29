/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
/* ode-initval/evolve.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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
#ifndef O2SCL_GSL_ASTEP_OLD_H
#define O2SCL_GSL_ASTEP_OLD_H

#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <o2scl/adapt_step.h>
#include <o2scl/gsl_astep.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Adaptive ODE stepper (GSL)

      [This is the stepper from \o2 version 0.908 which is 
      kept around for diagnostic purposes and some legacy code.]
      
      This class performs an adaptive step of a system of ODEs. To
      modify the ODE stepper which is used, use the function
      \ref adapt_step::set_step().

      There is an example for the usage of this class in
      <tt>examples/ex_ode.cpp</tt> documented in the \ref ex_ode_sect
      section.

      Default template arguments
      - \c func_t - \ref ode_funct
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
  */
  template<class func_t=ode_funct<>, 
    class vec_t=ubvector, class alloc_vec_t=ubvector, 
    class alloc_t=ubvector_alloc> class gsl_astep_old : 
  public adapt_step<func_t,vec_t,alloc_vec_t,alloc_t> {
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// Memory allocator for objects of type \c alloc_vec_t
    alloc_t ao;
    
    /// Temporary storage for yout
    alloc_vec_t yout;

    /// Internal storage for dydx
    alloc_vec_t dydx_int;
    
    /// The size of the last step
    double last_step;

    /// The number of steps
    unsigned long int count;

    /// The number of failed steps
    unsigned long int failed_steps;

    /// The size of the allocated vectors
    size_t msize;

    /// Apply the evolution for the next adaptive step
    int evolve_apply(double t0, double &t, double &h, double t1, size_t nvar,
		     vec_t &y, vec_t &dydx, vec_t &yout2, vec_t &yerr,
		     vec_t &dydx_out2, func_t &derivs) {
      
      double h0 = h;
      int step_status;
      int final_step = 0;
      /* remaining step size, possibly less than h */
      double dt = t1 - t0;  
    
      if ((dt < 0.0 && h0 > 0.0) || (dt > 0.0 && h0 < 0.0)) {
	std::string str="Interval direction (t1-t0="+dtos(dt)+
	  ") does not match step direction (h="+dtos(h)+
	  " in gsl_astep_old::evolve_apply().";
	O2SCL_ERR_RET(str.c_str(),gsl_einval);
      }
      
      bool try_step_again=true;
      while(try_step_again) {

	if ((dt >= 0.0 && h0 > dt) || (dt < 0.0 && h0 < dt)) {
	  h0 = dt;
	  final_step = 1;
	} else{
	  final_step = 0;
	}
	
	step_status=this->stepp->step(t0,h0,nvar,y,dydx,yout2,yerr,
				      dydx_out2,derivs);
	
	/* Check for stepper internal failure */
	if (step_status != gsl_success)	{
	  // Notify user which step-size caused the failure
	  h=h0;
	  return step_status;
	}
  
	count++;
	last_step = h0;
	if (final_step) {
	  t = t1;
	} else{
	  t = t0 + h0;
	}
  
	/* Check error and attempt to adjust the step. */
	const int hadjust_status=con.hadjust(nvar,this->stepp->get_order(),
					     yout2,yerr,dydx_out2,h0);
	
	if (hadjust_status == GSL_ODEIV_HADJ_DEC) {
	  
	  /* Step was decreased. Undo and go back to try again. */
	  failed_steps++;
	  try_step_again=true;
	} else {
	  try_step_again=false;
	}

      }
	
      /* suggest step size for next time-step */
      h = h0;  
  
      // This used to be return step_status, but the code never 
      // reaches this point if step_status is anything other than zero.
      return gsl_success;
    }
      
#endif
    
    public:
      
    gsl_astep_old() {
      this->verbose=0;
      msize=0;
    }
      
    virtual ~gsl_astep_old() {
      if (msize>0) {
	ao.free(yout);
	ao.free(dydx_int);
      }
    }
    
    /// Control specification
    gsl_ode_control<vec_t> con;
    
    /** \brief Make an adaptive integration step of the system 
	\c derivs with derivatives
	  
	This attempts to take a step of size \c h from the point \c x
	of an \c n-dimensional system \c derivs starting with \c y and
	given the initial derivatives \c dydx. On exit, \c x, \c y and
	\c dydx contain the new values at the end of the step, \c h
	contains the size of the step, \c dydx contains the derivative
	at the end of the step, and \c yerr contains the estimated
	error at the end of the step.

	If the base stepper returns a non-zero value, the step is
	aborted and \c y is unmodified. The error handler is never
	called.
    */
    virtual int astep_derivs(double &x, double xmax, double &h, 
			     size_t n, vec_t &y, vec_t &dydx, 
			     vec_t &yerr, func_t &derivs) {
      count=0;
      failed_steps=0;
      last_step=0.0;

      if (n!=msize) {
	ao.allocate(yout,n);
	ao.allocate(dydx_int,n);
	msize=n;
      }

      int ret=evolve_apply(x,x,h,xmax,n,y,dydx,yout,yerr,dydx_int,
			   derivs);

      if (ret==gsl_success) {
	for(size_t i=0;i<n;i++) {
	  y[i]=yout[i];
	  dydx[i]=dydx_int[i];
	}
      }
	
      if (this->verbose>0) {
	std::cout << x << " ";
	for(size_t j=0;j<n;j++) std::cout << y[j] << " ";
	std::cout << std::endl;
      }

      return ret;
    }
  
    /** \brief Make an adaptive integration step of the system 
	\c derivs

	This function performs an adaptive integration step with the
	\c n-dimensional system \c derivs and parameter \c pa. It
	Begins at \c x with initial stepsize \c h, ensuring that the
	step goes no farther than \c xmax. At the end of the step, the
	size of the step taken is \c h and the new value of \c x is in
	\c x_out. Initially, the function values and derivatives
	should be specified in \c y and \c dydx. The function values,
	derivatives, and the error at the end of the step are given in
	\c yout, \c yerr, and \c dydx_out. Unlike in \c ode_step
	objects, the objects \c y, \c yout, \c dydx, and \c dydx_out
	must all be distinct.

	If the base stepper returns a non-zero value, the step is
	aborted. The error handler is never called.

	This adaptive stepper function is faster than astep() or
	astep_derivs() because it does not require any copying of
	vectors.
    */
    virtual int astep_full(double x, double xmax, double &x_out, double &h, 
			   size_t n, vec_t &y, vec_t &dydx, 
			   vec_t &yout2, vec_t &yerr, vec_t &dydx_out, 
			   func_t &derivs) {

      count=0;
      failed_steps=0;
      last_step=0.0;
    
      int ret=evolve_apply(x,x_out,h,xmax,n,y,dydx,yout2,yerr,dydx_out,
			   derivs);
    
      if (this->verbose>0) {
	std::cout << x_out << " ";
	for(size_t j=0;j<n;j++) std::cout << yout2[j] << " ";
	std::cout << std::endl;
      }

      return ret;
    }
      
    /** \brief Make an adaptive integration step of the system 
	\c derivs 
	  
	This attempts to take a step of size \c h from the point \c
	x of an \c n-dimensional system \c derivs starting with \c
	y. On exit, \c x and \c y contain the new values at the end
	of the step, \c h contains the size of the step, \c dydx_out
	contains the derivative at the end of the step, and \c yerr
	contains the estimated error at the end of the step.

	If the system \c derivs or the base stepper return a non-zero 
	value, the adaptive step is aborted and \c y is unmodified. 
	The error handler is never called.
    */
    virtual int astep(double &x, double xmax, double &h, 
		      size_t n, vec_t &y, vec_t &dydx_out,
		      vec_t &yerr, func_t &derivs) {
      count=0;
      failed_steps=0;
      last_step=0.0;

      if (n!=msize) {
	ao.allocate(yout,n);
	ao.allocate(dydx_int,n);
	msize=n;
      }
      
      int ret=derivs(x,n,y,dydx_int);
      if (ret!=0) return ret;
      
      ret=evolve_apply(x,x,h,xmax,n,y,dydx_int,yout,yerr,dydx_out,
		       derivs);

      if (ret==gsl_success) {
	for(size_t i=0;i<n;i++) {
	  y[i]=yout[i];
	}
      }

      if (this->verbose>0) {
	std::cout << x << " ";
	for(size_t j=0;j<n;j++) std::cout << y[j] << " ";
	std::cout << std::endl;
      }

      return ret;
    }
      
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
