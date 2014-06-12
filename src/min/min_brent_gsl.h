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
/* min/brent.c
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
#ifndef O2SCL_MIN_BRENT_GSL_H
#define O2SCL_MIN_BRENT_GSL_H

/** \file min_brent_gsl.h
    \brief File defining \ref o2scl::min_brent_gsl
*/
#include <limits>

// For gsl_min_test_interval()
#include <gsl/gsl_min.h>

#include <o2scl/min.h>


#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
   
  /** \brief One-dimensional minimization using Brent's method (GSL)

      The minimization in the function min_bkt() is complete when the
      bracketed interval is smaller than 
      \f$ \mathrm{tol} = \mathrm{tol\_abs} + \mathrm{tol\_rel} \cdot
      \mathrm{min} \f$, where \f$ \mathrm{min} = 
      \mathrm{min}(|\mathrm{lower}|,|\mathrm{upper}|) \f$.

      Note that this algorithm requires that the initial guess
      already brackets the minimum, i.e. \f$ x_1 < x_2 < x_3 \f$,
      \f$ f(x_1) > f(x_2) \f$ and \f$ f(x_3) > f(x_2) \f$. This is
      different from \ref min_cern, where the initial value
      of the first parameter to min_cern::min_bkt() is
      ignored.

      The set functions throw an error if the initial bracket is not
      correctly specified. The function \ref iterate() never calls the
      error handler. The function min_bkt() calls the error handler if
      the tolerances are negative or if the number of iterations is
      insufficient to give the specified tolerance.

      Setting \ref min_base::err_nonconv to false will force
      min_bkt() not to call the error handler when the number of
      iterations is insufficient.

      Note that if there are more than 1 local minima in the specified
      interval, there is no guarantee that this method will find the
      global minimum.

      See also \ref root_brent_gsl for a similar algorithm
      applied as a solver rather than a minimizer.

      \note There was a bug in this minimizer which was fixed for
      GSL-1.11 which has also been fixed here. 
  */
  template<class func_t=funct11> class min_brent_gsl : 
  public min_bkt_base<func_t> {
    
#ifndef DOXYGEN_INTERNAL
    
    protected:

    /// The function
    func_t *uf;
      
    /// \name Temporary storage
    //@{
    double d, e, v, w, f_v, f_w;
    //@}
  
    /// Compute the function values at the various points
    int compute_f_values(func_t &func, double xminimum, double *fminimum, 
			 double xlower, double *flower, double xupper,
			 double *fupper) {
    
      *flower=func(xlower);
      *fupper=func(xupper);
      *fminimum=func(xminimum);
    
      return success;
    }
      
#endif
      
    public:
 
    /// Location of minimum
    double x_minimum;
    /// Lower bound
    double x_lower;
    /// Upper bound
    double x_upper;
    /// Minimum value
    double f_minimum;
    /// Value at lower bound
    double f_lower;
    /// Value at upper bound
    double f_upper;

    min_brent_gsl() {
    }
      
    virtual ~min_brent_gsl() {}

    /// Set the function and the initial bracketing interval
    int set(func_t &func, double xmin, double lower, double upper) {
      
      double fmin, fl, fu;
      int status=compute_f_values(func,lower,&fl,xmin,&fmin,upper,&fu);
    
      status=set_with_values(func,xmin,fmin,lower,fl,upper,fu);

      return status;
    }

    /** \brief Set the function, the initial bracketing interval,
	and the corresponding function values.
    */
    int set_with_values(func_t &func, double xmin, double fmin,
			double lower, double fl, double upper,
			double fu) {

      if (lower > upper) {
	std::string tmp=((std::string)"Invalid interval (lower > upper) ")+
	"in min_brent_gsl::set_with_values().";
	O2SCL_ERR(tmp.c_str(),exc_einval);
      }
      if (xmin>=upper || xmin<=lower) {
	std::string tmp=((std::string)"'xmin' was not inside interval ")+
	"in min_brent_gsl::set_with_values().";
	O2SCL_ERR(tmp.c_str(),exc_einval);
      }
      if (fmin>=fl || fmin>=fu) {
	std::string tmp=((std::string)"Endpoints don't enclose minimum ")+
	"in min_brent_gsl::set_with_values().";
	O2SCL_ERR(tmp.c_str(),exc_einval);
      }
    
      x_lower=lower;
      x_minimum=xmin;
      x_upper=upper;
      f_lower=fl;
      f_minimum=fmin;
      f_upper=fu;

      uf=&func;

      /* golden=(3-sqrt(5))/2 */
      const double golden=0.3819660;      
    
      v=x_lower+golden*(x_upper-x_lower);
      w=v;
      d=0;
      e=0;
      f_v=func(v);
      f_w=f_v;
    
      return success;
    }
  
    /** \brief Perform an iteration
	
        \future It looks like x_left and x_right can be removed. Also,
	it would be great to replace the one-letter variable names
	with something more meaningful.
     */
    int iterate() {
    
      const double x_left=x_lower;
      const double x_right=x_upper;
	
      const double z=x_minimum;
      double u, f_u;
      const double f_z=f_minimum;
	
      /* golden=(3-sqrt(5))/2 */
      const double golden=0.3819660;      
    
      const double w_lower=(z-x_left);
      const double w_upper=(x_right-z);
    
#ifndef O2SCL_NO_CPP11
      double sqrt_dbl_eps=sqrt(std::numeric_limits<double>::epsilon());
#else 
      double sqrt_dbl_eps=GSL_SQRT_DBL_EPSILON;
#endif

      const double tolerance=sqrt_dbl_eps*fabs(z);
    
      double p=0, q=0, r=0;
    
      const double midpoint=0.5*(x_left+x_right);
    
      if (fabs (e) > tolerance) {
	/* fit parabola */
	  
	r=(z-w)*(f_z-f_v);
	q=(z-v)*(f_z-f_w);
	p=(z-v)*q-(z-w)*r;
	q=2*(q-r);
	  
	if (q > 0) {
	  p=-p;
	} else {
	  q=-q;
	}
	  
	r=e;
	e=d;
      }

      if (fabs (p) < fabs (0.5*q*r) && p < q*w_lower && 
	  p < q*w_upper) {

	double t2=2*tolerance;
	  
	d=p / q;
	u=z+d;
	  
	if ((u-x_left) < t2 || (x_right-u) < t2) {
	  d=(z < midpoint) ? tolerance : -tolerance ;
	}

      } else {

	e=(z < midpoint) ? x_right-z : -(z-x_left) ;
	d=golden*e;

      }
    
    
      if (fabs (d) >= tolerance) {
	u=z+d;
      } else {
	u=z+((d > 0) ? tolerance : -tolerance);
      }
	
      f_u=(*uf)(u);
      
      if (f_u<=f_z) {

	if (u<z) {

	  x_upper=z;
	  f_upper=f_z;

	} else {

	  x_lower=z;
	  f_lower=f_z;

	}

	v=w;
	f_v=f_w;
	w=z;
	f_w=f_z;
	x_minimum=u;
	f_minimum=f_u;

	return success;

      } else {

	if (u<z) {

	  x_lower=u;
	  f_lower=f_u;

	} else {

	  x_upper=u;
	  f_upper=f_u;

	}

	if (f_u<=f_w || w==z) {

	  v=w;
	  f_v=f_w;
	  w=u;
	  f_w=f_u;
	  return success;

	} else if (f_u<=f_v || v==z || v==w) {

	  v=u;
	  f_v=f_u;
	  return success;
	}

      }

      return success;
    }
  
    /** \brief Calculate the minimum \c fmin of \c func 
	with \c x2 bracketed between \c x1 and \c x3.

	Note that this algorithm requires that the initial guess
	already brackets the minimum, i.e. \f$ x_1 < x_2 < x_3 \f$,
	\f$ f(x_1) > f(x_2) \f$ and \f$ f(x_3) > f(x_2) \f$. This is
	different from \ref min_cern, where the initial value
	of the first parameter to min_cern::min_bkt() is
	ignored.
    */
    virtual int min_bkt(double &x2, double x1, double x3, double &fmin,
			func_t &func) {
    
      int status, iter=0;
    
      int rx=set(func,x2,x1,x3);
      if (rx!=0) {
	O2SCL_ERR2("Function set() failed in ",
		       "min_brent_gsl::min_bkt().",rx);
      }
	
      do {
	iter++;
	status=iterate();
	x2=x_minimum;
	if (status) {
	  // This should never actually happen in the current
	  // version, but we retain it for now
	  O2SCL_CONV2_RET("Function iterate() failed in gsl_min_",
			  "brent::min_bkt().",status,this->err_nonconv);
	  break;
	}
	  
	status=gsl_min_test_interval(x_lower,x_upper,this->tol_abs,
				     this->tol_rel);
	if (status>0) {
	  // The function gsl_min_test_interval() fails if the
	  // tolerances are negative or if the lower bound is larger
	  // than the upper bound
	  std::string s="Function gsl_min_test_interval() failed in ";
	  s+="min_brent_gsl::min_bkt().";
	  O2SCL_ERR(s.c_str(),status);
	}
	  
	if (this->verbose>0) {
	  if (x_lower*x_upper<0.0) {
	    if (x_lower<x_upper) {
	      this->print_iter(x2,f_minimum,iter,fabs(x_lower-x_upper)-
			       this->tol_rel*x_lower,this->tol_abs,
			       "min_brent_gsl");
	    } else {
	      this->print_iter(x2,f_minimum,iter,fabs(x_lower-x_upper)-
			       this->tol_rel*x_upper,this->tol_abs,
			       "min_brent_gsl");
	    }
	  } else {
	    this->print_iter(x2,f_minimum,iter,fabs(x_lower-x_upper),
			     this->tol_abs,"min_brent_gsl");
	  }
	}
      
      } while (status == gsl_continue && iter<this->ntrial);
      
      if (iter>=this->ntrial) {
	O2SCL_CONV2_RET("Exceeded maximum number of ",
			"iterations in min_brent_gsl::min_bkt().",
			o2scl::exc_emaxiter,this->err_nonconv);
      }

      this->last_ntrial=iter;
      fmin=f_minimum;
	
      return status;
    }
  
    /// Return string denoting type ("min_brent_gsl")
    virtual const char *type() { return "min_brent_gsl"; }
      
    };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
