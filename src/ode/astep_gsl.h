/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
/* ode-initval2/evolve.c
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
#ifndef O2SCL_ASTEP_GSL_H
#define O2SCL_ASTEP_GSL_H

/** \file astep_gsl.h
    \brief File defining \ref o2scl::astep_gsl
*/

#include <string>

#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv.h>

#include <o2scl/astep.h>
#include <o2scl/vector.h>
#include <o2scl/string_conv.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Control structure for astep_gsl

      This class implements both the "standard" and "scaled" step
      control methods from GSL. The standard control method is the
      default. To use the scaled control, set \ref standard to
      <tt>false</tt> and set the scale for each component using
      set_scale().
      
      The control object is a four parameter heuristic based on
      absolute and relative errors \ref eps_abs and \ref eps_rel, and
      scaling factors \ref a_y and \ref a_dydt for the system state
      \f$ y(t) \f$ and derivatives \f$ y^{\prime}(t) \f$ respectively.

      The step-size adjustment procedure for this method begins by
      computing the desired error level \f$ D_i \f$ for each
      component. In the unscaled version,
      \f[
      D_i = \mathrm{eps\_abs}+\mathrm{eps\_rel} \times
      \left( \mathrm{a\_y} | y_i| + \mathrm{a\_dydt}~h 
      | y_i^{\prime}| \right)
      \f]
      while in the scaled version the user specifies the scale for
      each component, \f$ s_i \f$,
      \f[
      D_i = \mathrm{eps\_abs}~s_i+\mathrm{eps\_rel} \times
      \left( \mathrm{a\_y} | y_i| + \mathrm{a\_dydt}~h 
      | y_i^{\prime}| \right)
      \f]

      The desired error level \f$ D_i \f$ is compared to then observed
      error \f$ E_i = |\mathrm{yerr}_i| \f$. If the observed error \f$
      E \f$ exceeds the desired error level \f$ D \f$ by more than 10
      percent for any component then the method reduces the step-size
      by an appropriate factor,
      \f[
      h_{\mathrm{new}} = S~h_{\mathrm{old}} 
      \left(\frac{E}{D}\right)^{-1/q}
      \f]
      where \f$ q \f$ is the consistency order of the method (e.g. \f$
      q=4 \f$ for 4(5) embedded RK), and \f$ S \f$ is a safety factor
      of 0.9. The ratio \f$ E/D \f$ is taken to be the maximum of the
      ratios \f$ E_i/D_i \f$.
      
      If the observed error E is less than 50 percent of the desired
      error level \f$ D \f$ for the maximum ratio \f$ E_i/D_i \f$ then
      the algorithm takes the opportunity to increase the step-size to
      bring the error in line with the desired level,
      \f[
      h_{\mathrm{new}} = S~h_{\mathrm{old}} 
      \left(\frac{E}{D}\right)^{-1/(q+1)}
      \f]
      This encompasses all the standard error scaling methods. To
      avoid uncontrolled changes in the stepsize, the overall scaling
      factor is limited to the range 1/5 to 5.

      If the user specified fewer scaling parameters than the number
      of ODEs, then the scaling parameters are reused as follows. If
      there are \f$ N \f$ ODEs and \f$ M \f$ scaling parameters, then
      for \f$ i>M \f$, the ith scaling parameter \f$ s_i \f$ is set to
      be \f$ s_{i\%M} \f$ . If the user selects the scaled control by
      setting \ref standard to <tt>false</tt> and no scale parameters
      are specified, this class reverts to the standard control.

      \todo Double check that the improvements in the ode-initval2
      routines are available here
  */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t>
    class ode_control_gsl {

  public:  

  typedef boost::numeric::ublas::vector<double> ubvector;
  
  protected:
  
  /// Scalings
  ubvector scale_abs;
  
  public:
  
  /// \name Adjustment specification
  //@{
  /// No adjustment required
  static const size_t hadj_nil=0;
  /// Recommend step decrease
  static const size_t hadj_dec=1;
  /// Recommend step increase
  static const size_t hadj_inc=2;
  //@}
  
  /// Absolute precision (default \f$ 10^{-6} \f$)
  double eps_abs;
  /// Relative precision (default 0)
  double eps_rel;
  /// Function scaling factor (default 1)
  double a_y;
  /// Derivative scaling factor (default 0)
  double a_dydt;
  /// Use standard or scaled algorithm (default true)
  bool standard;

  ode_control_gsl() {
    eps_abs=1.0e-6;
    eps_rel=0.0;
    a_y=1.0;
    a_dydt=0.0;
    standard=true;
  }
    
  virtual ~ode_control_gsl() {
  }      

  /** \brief Set the scaling for each differential equation
   */
  template<class svec_t> int set_scale(size_t nscal, const svec_t &scale) {
    scale_abs.resize(nscal);
    for(size_t i=0;i<nscal;i++) {
      scale_abs[i]=scale[i];
    }
    return 0;
  }

  virtual int hadjust(size_t dim, unsigned int ord, const vec_y_t &y,
		      vec_yerr_t &yerr, vec_dydx_t &yp, double &h) {
      
    const double S=0.9;
    const double h_old=h;
      
    double rmax=DBL_MIN;
    size_t i;
      
    for(i=0;i<dim;i++) {
      double D0;
      if (standard || scale_abs.size()==0) {
	D0=eps_rel*(a_y*fabs(y[i])+a_dydt*fabs(h_old*yp[i]))+eps_abs;
      } else {
	D0=eps_rel*(a_y*fabs(y[i])+a_dydt*fabs(h_old*yp[i]))+
	  eps_abs*scale_abs[i%scale_abs.size()];
      }
      const double r=fabs(yerr[i]) / fabs(D0);
      rmax=GSL_MAX_DBL(r, rmax);
    }
	
    if (rmax > 1.1) {

      /* decrease step, no more than factor of 5, but a fraction S more
	 than scaling suggests (for better accuracy) */
      double r= S / pow(rmax, 1.0/ord);
	  
      if (r < 0.2) {
	r=0.2;
      }
      h=r*h_old;
      return hadj_dec;

    } else if (rmax < 0.5) {

      /* increase step, no more than factor of 5 */
      double r=S / pow(rmax, 1.0/(ord+1.0));
      if (r > 5.0) {
	r=5.0;
      }
      /* don't allow any decrease caused by S<1 */
      if (r < 1.0) {
	r=1.0;
      }
      h=r*h_old;
      return hadj_inc;

    } else {

      /* no change */
      return hadj_nil;

    }
    
  }
  
  };
  
  /** \brief Adaptive ODE stepper (GSL)
      
      This class performs an adaptive step of a system of ODEs. To
      modify the ODE stepper which is used, use the function \ref
      astep_base::set_step().

      Note, this has been updated to correspond to the
      <tt>ode-initval2</tt> functions in GSL.

      There is an example for the usage of this class in
      <tt>examples/ex_ode.cpp</tt> documented in the \ref ex_ode_sect
      section.

      \todo Document what happens when the stepper function returns
      a non-zero value, as it's different now with the ode-initval2
      function.
      \todo Document count, failed_steps, etc.
      \future Compare more directly to GSL

      Default template arguments
      - \c func_t - \ref ode_funct
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
  */
  template<class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t, class vec_yerr_t=vec_y_t, 
    class func_t=ode_funct > class astep_gsl : 
    public astep_base<vec_y_t,vec_dydx_t,vec_yerr_t,func_t> {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;

#ifndef DOXYGEN_INTERNAL
    
  protected:
    
  /// Temporary storage for yout
  vec_y_t yout_int;

  /// Internal storage for dydx
  vec_dydx_t dydx_int;
    
  /// The size of the last step
  double last_step;

  /// The number of steps
  unsigned long int count;

  /// The number of failed steps
  unsigned long int failed_steps;

  /// The size of the allocated vectors
  size_t msize;

  /** \brief Apply the evolution for the next adaptive step

      This function is based on <tt>gsl_odeiv2_evolve_apply</tt>.
      
      \note This function requres that \c y, \c yout, \c dydx and \c
      dydx_out are all distinct vectors. 
  */
  int evolve_apply(double t0, double t1, double &t, double &h, size_t nvar,
		   vec_y_t &y, vec_dydx_t &dydx, vec_y_t &yout, 
		   vec_yerr_t &yerr, vec_dydx_t &dydx_out, 
		   func_t &derivs) {

    double h0=h;
    int step_status;
    bool final_step=false;
    // Remaining step size, possibly less than h
    double dt=t1-t0; 
    
    if ((dt < 0.0 && h0 > 0.0) || (dt > 0.0 && h0 < 0.0)) {
      std::string str="Interval direction (t1-t0="+o2scl::dtos(dt)+
	") does not match step direction (h="+o2scl::dtos(h)+
	" in astep_gsl::evolve_apply().";
      O2SCL_ERR(str.c_str(),exc_einval);
    }

    // Copy the initial point from y to yout. The stepper won't modify
    // y, so we can undo the step by going back to the values in y.
    vector_copy(nvar,y,yout);
      
    bool try_step_again=true;
    while (try_step_again) {
      try_step_again=false;
      
      if ((dt >= 0.0 && h0 > dt) || (dt < 0.0 && h0 < dt)) {
	h0=dt;
	final_step=true;
      } else{
	final_step=false;
      }
	
      step_status=this->stepp->step(t0,h0,nvar,yout,dydx,yout,yerr,
				    dydx_out,derivs);

      // [GSL] Return if the stepper indicates a pointer or user function
      // failure
      if (step_status == exc_efault || step_status==exc_ebadfunc) {
	return step_status;
      }

      // [GSL] Check for stepper internal failure
      if (step_status!=success) {

	// [GSL] Stepper was unable to calculate step. Try decreasing
	// stepsize
	double h_old=h0;
	h0*=0.5;

	// [GSL] Check that an actual decrease in h0 occured and that
	// the suggested h0 will change the independent variable
	// will change by at least 1ulp
	//double t_curr=GSL_COERCE_DBL(t);
	//double t_next=GSL_COERCE_DBL(t+h0);
	double t_curr=t;
	double t_next=t+h0;
	
	if (fabs(h0) < fabs(h_old) && t_next != t_curr) {
	  // [GSL] Step was decreased. Undo step, and try again with 
	  // new h0. 
	  vector_copy(nvar,y,yout);
	  failed_steps++;
	  try_step_again=true;
	} else {
	  // [GSL] notify user of step-size which caused the failure 
	  h=h0;            
	  // [GSL] restore original t value
	  t=t0;            
	  return step_status;
	}
      }

      if (try_step_again==false) {
      
	count++;
	last_step=h0;
	if (final_step) {
	  t=t1;
	} else{
	  t=t0+h0;
	}
  
	// [GSL] Check error and attempt to adjust the step.
	double h_old=h0;
	const size_t hadjust_status=con.hadjust
	  (nvar,this->stepp->get_order(),yout,yerr,dydx_out,h0);
      
	if (hadjust_status == 
	    ode_control_gsl<vec_y_t,vec_dydx_t,vec_yerr_t>::hadj_dec) {
	
	  // [GSL] Check that the reported status is correct (i.e. an actual
	  // decrease in h0 occured) and the suggested h0 will change
	  // the independent variable by at least 1 ulp

	  //double t_curr=GSL_COERCE_DBL (*t);
	  //double t_next=GSL_COERCE_DBL ((*t) + h0);
	  double t_curr=t;
	  double t_next=t+h0;
	  
	  if (fabs(h0) < fabs(h_old) && t_next != t_curr) {
	    // [GSL] Step was decreased. Undo step, and try again with 
	    // new h0.
	    vector_copy(nvar,y,yout);
	    failed_steps++;
	    try_step_again=true;
	  } else {
	    /* [GSL] Can not obtain required error tolerance, and can
	       not decrease step-size any further, so give up and return
	       gsl_failure.
	    */
	    // [GSL] notify user of step-size which caused the failure
	    h=h0;         
	    return gsl_failure;
	  }
	}
      }
      
    }

    if (final_step==false) {
      // [GSL] Suggest step size for next time-step. Change of step 
      // size is not suggested in the final step, because that step 
      // can be very small compared to previous step, to reach t1. 
      h=h0; 
    }
  
    // This used to be return step_status, but the code never 
    // reaches this point if step_status is anything other than zero.
    return step_status;
  }
      
#endif
    
  public:
      
  astep_gsl() {
    this->verbose=0;
    msize=0;
  }
      
  virtual ~astep_gsl() {
  }
    
  /// Control specification
  ode_control_gsl<vec_y_t,vec_dydx_t,vec_yerr_t> con;
    
  /** \brief Make an adaptive integration step of the system 
      \c derivs 
	  
      This attempts to take a step of size \c h from the point \c
      x of an \c n-dimensional system \c derivs starting with \c
      y. On exit, \c x and \c y contain the new values at the end
      of the step, \c h contains the size of the step, \c dydx_out
      contains the derivative at the end of the step, and \c yerr
      contains the estimated error at the end of the step.
  */
  virtual int astep(double &x, double xmax, double &h, 
		    size_t n, vec_y_t &y, vec_dydx_t &dydx_out,
		    vec_yerr_t &yerr, func_t &derivs) {
    count=0;
    failed_steps=0;
    last_step=0.0;

    if (n!=msize) {
      yout_int.resize(n);
      dydx_int.resize(n);
      msize=n;
    }
      
    int ret=derivs(x,n,y,dydx_int);
    if (ret!=0) return ret;
      
    ret=evolve_apply(x,xmax,x,h,n,y,dydx_int,yout_int,yerr,dydx_out,derivs);

    if (ret==success) {
      for(size_t i=0;i<n;i++) {
	y[i]=yout_int[i];
      }
    }

    if (this->verbose>0) {
      std::cout << "astep_gsl step:";
      std::cout << x << " ";
      for(size_t j=0;j<n;j++) std::cout << y[j] << " ";
      std::cout << std::endl;
    }

    return ret;
  }

  /** \brief Make an adaptive integration step of the system 
      \c derivs with derivatives
	  
      This attempts to take a step of size \c h from the point \c x
      of an \c n-dimensional system \c derivs starting with \c y and
      given the initial derivatives \c dydx. On exit, \c x, \c y and
      \c dydx contain the new values at the end of the step, \c h
      contains the size of the step, \c dydx contains the derivative
      at the end of the step, and \c yerr contains the estimated
      error at the end of the step.
  */
  virtual int astep_derivs(double &x, double xmax, double &h, 
			   size_t n, vec_y_t &y, vec_dydx_t &dydx, 
			   vec_yerr_t &yerr, func_t &derivs) {
    count=0;
    failed_steps=0;
    last_step=0.0;

    if (n!=msize) {
      yout_int.resize(n);
      dydx_int.resize(n);
      msize=n;
    }

    int ret=evolve_apply(x,xmax,x,h,n,y,dydx,yout_int,yerr,dydx_int,derivs);

    if (ret==success) {
      for(size_t i=0;i<n;i++) {
	y[i]=yout_int[i];
	dydx[i]=dydx_int[i];
      }
    }
	
    if (this->verbose>0) {
      std::cout << "astep_gsl step:";
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

      This adaptive stepper function is faster than astep() or
      astep_derivs() because it does not require any copying of
      vectors.
  */
  virtual int astep_full(double x, double xmax, double &x_out, double &h, 
			 size_t n, vec_y_t &y, vec_dydx_t &dydx, 
			 vec_y_t &yout, vec_yerr_t &yerr, 
			 vec_dydx_t &dydx_out, func_t &derivs) {

    count=0;
    failed_steps=0;
    last_step=0.0;
    
    int ret=evolve_apply(x,xmax,x_out,h,n,y,dydx,yout,yerr,dydx_out,
			 derivs);
    
    if (this->verbose>0) {
      std::cout << "astep_gsl step:";
      std::cout << x_out << " ";
      for(size_t j=0;j<n;j++) std::cout << yout[j] << " ";
      std::cout << std::endl;
    }

    return ret;
  }
            
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
