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
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*  quad_golden.c                                                           */
/*                                                                          */
/*  Copyright (C) 2007 James Howse                                          */
/*  Copyright (C) 2009 Brian Gough                                          */
/*                                                                          */
/* This program is free software; you can redistribute it and/or modify     */
/* it under the terms of the GNU General Public License as published by     */
/* the Free Software Foundation; either version 3 of the License, or (at    */
/* your option) any later version.                                          */
/*                                                                          */
/* This program is distributed in the hope that it will be useful, but      */
/* WITHOUT ANY WARRANTY; without even the implied warranty of               */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        */
/* General Public License for more details.                                 */
/*                                                                          */
/* You should have received a copy of the GNU General Public License        */
/*  along with this program; if not, write to the Free Software             */
/*  Foundation, Inc., 51 Franklin Street, Fifth Floor,                      */
/*  Boston, MA 02110-1301, USA.                                             */
/*--------------------------------------------------------------------------*/
#ifndef O2SCL_MIN_QUAD_GOLDEN_H
#define O2SCL_MIN_QUAD_GOLDEN_H

/** \file min_quad_golden.h
    \brief File defining \ref o2scl::min_quad_golden
*/

#include <limits>
#include <gsl/gsl_min.h>
#include <o2scl/min.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
   
  /** \brief Minimization of a function using the safeguarded step-length
      algorithm of Gill and Murray [GSL] 

      This class is unfinished. 
      
      Documentation from GSL:
      \verbatim
                                                                            
      This algorithm performs univariate minimization (i.e., line
      search). It requires only objective function values g(x) to
      compute the minimum. The algorithm maintains an interval of
      uncertainty [a,b] and a point x in the interval [a,b] such that
      a < x < b, and g(a) > g(x) and g(x) < g(b). The algorithm also
      maintains the three points with the smallest objective values x,
      v and w such that g(x) < g(v) < g(w). The algorithm terminates
      when max( x-a, b-x ) < 2(r |x| + t) where r and t are small
      positive reals. At a given iteration, the algorithm first fits a
      quadratic through the three points (x, g(x)), (v, g(v)) and (w,
      g(w)) and computes the location of the minimum u of the
      resulting quadratic. If u is in the interval [a,b] then g(u) is
      computed. If u is not in the interval [a,b], and either v < x
      and w < x, or v > x and w > x (i.e., the quadratic is
      extrapolating), then a point u' is computed using a safeguarding
      procedure and g(u') is computed. If u is not in the interval
      [a,b], and the quadratic is not extrapolating, then a point u''
      is computed using approximate golden section and g(u'') is
      computed. After evaluating g() at the appropriate new point, a,
      b, x, v, and w are updated and the next iteration is performed.
      The algorithm is based on work presented in the following
      references.
                                                                            
      Algorithms for Minimization without derivatives
      Richard Brent
      Prentice-Hall Inc., Englewood Cliffs, NJ, 1973
      
      Safeguarded Steplength Algorithms for Optimization using Descent Methods
      Philip E. Gill and Walter Murray
      Division of Numerical Analysis and Computing
      National Physical Laboratory, Teddington, United Kingdom
      NPL Report NAC 37, August 1974
      \endverbatim

      \future Take common elements of this and min_brent and
      move to a generic GSL minimizer type? 
  */
  template<class func_t=funct> class min_quad_golden : 
  public min_bkt_base<func_t> {
	
#ifndef DOXYGEN_INTERNAL
	
    protected:

    /// The function
    func_t *uf;
      
    double x_prev_small;
    double f_prev_small;
    double x_small;
    double f_small;
    double step_size;
    double stored_step;
    double prev_stored_step;
    size_t num_iter;

    double rel_err_val;
    double abs_err_val;
    double golden_mean;
    double golden_ratio;

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

    min_quad_golden() {
      rel_err_val=1.0e-06;
      abs_err_val=1.0e-10;
      golden_mean=0.3819660112501052;
      golden_ratio=1.6180339887498950;
    }
      
    virtual ~min_quad_golden() {}

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
	"in min_quad_golden::set_with_values().";
	O2SCL_ERR(tmp.c_str(),exc_einval);
      }
      if (xmin>=upper || xmin<=lower) {
	std::string tmp=((std::string)"'xmin' was not inside interval ")+
	"in min_quad_golden::set_with_values().";
	O2SCL_ERR(tmp.c_str(),exc_einval);
      }
      if (fmin>= fl || fmin>=fu) {
	std::string tmp=((std::string)"Endpoints don't enclose minimum ")+
	"in min_quad_golden::set_with_values().";
	O2SCL_ERR(tmp.c_str(),exc_einval);
      }
    
      x_lower=lower;
      x_minimum=xmin;
      x_upper=upper;
      f_lower=fl;
      f_minimum=fmin;
      f_upper=fu;

      uf=&func;

      x_prev_small=xmin;
      x_small=xmin;
      f_prev_small=fmin;
      f_small=fmin;
      
      step_size=0.0;
      stored_step=0.0;
      prev_stored_step=0.0;
      num_iter=0;
    
      return success;
    }
  
    /** \brief Perform an iteration
     */
    int iterate() {

      const double x_m=x_minimum;
      const double f_m=f_minimum;

      const double x_l=x_lower;
      const double x_u=x_upper;
  
      double quad_step_size=prev_stored_step;
  
      double x_trial;
      double x_eval, f_eval;
      
      double x_midpoint=0.5*(x_l+x_u);
      double tol=rel_err_val*fabs(x_m)+abs_err_val; 

      double dbl_eps=std::numeric_limits<double>::epsilon();

      if (fabs(stored_step)-tol>-2.0*dbl_eps) {
	
	// [GSL] Fit quadratic
	double c3=(x_m-x_small)*(f_m-f_prev_small);
	double c2=(x_m-x_prev_small)*(f_m-f_small);
	double c1=(x_m-x_prev_small)*c2-(x_m-x_small)*c3;

	c2=2.0*(c2-c3);

	// [GSL] if (c2!=0)
	if (fabs(c2)>dbl_eps) {
	  if (c2>0.0) {
	    c1=-c1;
	  }
	  c2=fabs(c2);
	  quad_step_size=c1/c2;
	} else {
	  // [GSL] Handle case where c2 is nearly 0. Insure that the
	  // line search will NOT take a quadratic interpolation step
	  // in this iteration
	  quad_step_size=stored_step;
	}

	prev_stored_step=stored_step;
	stored_step=step_size;
      }

      x_trial=x_m+quad_step_size;

      if (fabs(quad_step_size)<fabs(0.5*prev_stored_step) && 
	  x_trial>x_l && x_trial<x_u) {

	// [GSL] Take quadratic interpolation step 
	step_size=quad_step_size;

	// [GSL] Do not evaluate function too close to x_l or x_u 
	if ((x_trial-x_l)<2.0*tol || (x_u-x_trial)<2.0*tol) {
	  step_size=(x_midpoint >= x_m ? +1.0 : -1.0)*fabs(tol);
	}
	
      } else if ((x_small!=x_prev_small && x_small<x_m && 
		  x_prev_small<x_m) ||
		 (x_small!=x_prev_small && x_small>x_m && 
		  x_prev_small>x_m)) {

	// [GSL] Take safeguarded function comparison step
	double outside_interval, inside_interval;

	if (x_small<x_m) {
	  outside_interval=x_l-x_m;
	  inside_interval=x_u-x_m;
	} else {
	  outside_interval=x_u-x_m;
	  inside_interval=x_l-x_m;
	}

	if (fabs(inside_interval) <= tol) {
	  // [GSL] Swap inside and outside intervals
	  double tmp=outside_interval;
	  outside_interval=inside_interval;
	  inside_interval=tmp;
	}

	{
	  double step=inside_interval;
	  double scale_factor;

	  if (fabs(outside_interval)<fabs(inside_interval)) {
	    scale_factor=0.5*sqrt (-outside_interval/inside_interval);
	  } else {
	    scale_factor=(5.0/11.0)*
	    (0.1-inside_interval/outside_interval);
	  }

	  stored_step=step;
	  step_size=scale_factor*step;
	}

      } else {

	// [GSL] Take golden section step
	double step;

	if (x_m<x_midpoint) {
	  step=x_u-x_m;
	} else {
	  step=x_l-x_m;
	}
	  
	stored_step=step;
	step_size=golden_mean*step;

      }

      // [GSL] Do not evaluate function too close to x_minimum 
      if (fabs(step_size)>tol) {
	x_eval=x_m+step_size;
      } else {
	x_eval=x_m+(step_size >= 0 ? +1.0 : -1.0)*fabs(tol);
      }
      
      
      // [GSL] Evaluate function at the new point x_eval
      f_eval=(*uf)(x_eval);

      // [GSL] Update {x,f}_lower, {x,f}_upper, {x,f}_prev_small, 
      // {x,f}_small, and {x,f}_minimum 
      if (f_eval <= f_m) {
	if (x_eval<x_m) {
	  x_upper=x_m;
	  f_upper=f_m;
	} else {
	  x_lower=x_m;
	  f_upper=f_m;
	}

	x_prev_small=x_small;
	f_prev_small=f_small;

	x_small=x_m;
	f_small=f_m;

	x_minimum=x_eval;
	f_minimum=f_eval;
      } else {
	if (x_eval<x_m) {
	  x_lower=x_eval;
	  f_lower=f_eval;
	} else {
	  x_upper=x_eval;
	  f_upper=f_eval;
	}

	if (f_eval <= f_small || fabs(x_small-x_m)<2.0*dbl_eps) {
	  x_prev_small=x_small;
	  f_prev_small=f_small;

	  x_small=x_eval;
	  f_small=f_eval;
	} else if (f_eval <= f_prev_small ||
		   fabs(x_prev_small-x_m)<2.0*dbl_eps ||
		   fabs(x_prev_small-x_small)<2.0*dbl_eps) {
	  x_prev_small=x_eval;
	  f_prev_small=f_eval;
	}
      }

      // Update for next iteration
      num_iter++;

      return success;
    }
  
    /** \brief Calculate the minimum \c fmin of \c func 
	with \c x2 bracketed between \c x1 and \c x3.

	Note that this algorithm requires that the initial guess
	already brackets the minimum, i.e. \f$ x_1<x_2<x_3 \f$,
	\f$ f(x_1)>f(x_2) \f$ and \f$ f(x_3)>f(x_2) \f$. This is
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
		       "min_quad_golden::min_bkt().",rx);
      }
	
      do {
	iter++;
	status=iterate();
	x2=x_minimum;
	if (status) {
	  // This should never actually happen in the current
	  // version, but we retain it for now
	  O2SCL_CONV2_RET("Function iterate() failed in gsl_min_",
			  "quad_golden::min_bkt().",status,this->err_nonconv);
	  break;
	}
	  
	status=gsl_min_test_interval(x_lower,x_upper,this->tol_abs,
				     this->tol_rel);
	if (status>0) {
	  // The function gsl_min_test_interval() fails if the
	  // tolerances are negative or if the lower bound is larger
	  // than the upper bound
	  std::string s="Function gsl_min_test_interval() failed in ";
	  s+="min_quad_golden::min_bkt().";
	  O2SCL_ERR(s.c_str(),status);
	}
	  
	if (this->verbose>0) {
	  if (x_lower*x_upper<0.0) {
	    if (x_lower<x_upper) {
	      this->print_iter(x2,f_minimum,iter,fabs(x_lower-x_upper)-
			       this->tol_rel*x_lower,this->tol_abs,
			       "min_quad_golden");
	    } else {
	      this->print_iter(x2,f_minimum,iter,fabs(x_lower-x_upper)-
			       this->tol_rel*x_upper,this->tol_abs,
			       "min_quad_golden");
	    }
	  } else {
	    this->print_iter(x2,f_minimum,iter,fabs(x_lower-x_upper),
			     this->tol_abs,"min_quad_golden");
	  }
	}
	
      } while (status == gsl_continue && iter<this->ntrial);
      
      if (iter>=this->ntrial) {
	O2SCL_CONV2_RET("Exceeded maximum number of ",
			"iterations in min_quad_golden::min_bkt().",
			o2scl::exc_emaxiter,this->err_nonconv);
      }

      this->last_ntrial=iter;
      fmin=f_minimum;
	
      return status;
    }
  
    /// Return string denoting type ("min_quad_golden")
    virtual const char *type() { return "min_quad_golden"; }
      
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
