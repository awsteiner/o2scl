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
#ifndef O2SCL_MINIMIZE_H
#define O2SCL_MINIMIZE_H

// For fabs()
#include <cmath>

#include <o2scl/err_hnd.h>
#include <o2scl/funct.h>
 
/** \file min.h
    \brief One-dimensional minimization base class and associated functions
*/

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief One-dimensional minimization [abstract base]
   */
  template<class func_t=funct, class dfunc_t=func_t> class min_base {
    
  public:
  
  min_base() {
    verbose=0;
    ntrial=100;
    tol_rel=1.0e-4;
    tol_abs=1.0e-4;
    last_ntrial=0;
    bracket_iter=20;
    err_nonconv=true;
  }
      
  virtual ~min_base() {}
      
  /// Output control
  int verbose;
    
  /// Maximum number of iterations
  int ntrial;
    
  /// The tolerance for the minimum function value
  double tol_rel;

  /// The tolerance for the location of the minimum
  double tol_abs;
    
  /// The number of iterations used in the most recent minimization
  int last_ntrial;
    
  /** \brief The number of iterations for automatically 
      bracketing a minimum (default 20)
  */
  int bracket_iter;

  /// If true, call the error handler if the routine does not "converge"
  bool err_nonconv;
    
  /** \brief Print out iteration information.
	
      Depending on the value of the variable \ref verbose, this
      prints out the iteration information. If verbose=0, then no
      information is printed, while if verbose>1, then after each
      iteration, the present values of x and y are output to
      std::cout along with the iteration number. If verbose>=2 then
      each iteration waits for a character.
  */
  virtual int print_iter(double x, double y, int iter, double value=0.0,
			 double limit=0.0, std::string comment="") {
    
    if (verbose<=0) return success;
    
    char ch;
    
    std::cout << comment << " Iteration: " << iter << std::endl;
    if (x<0) std::cout << x << " ";
    else std::cout << " " << x << " ";
    if (y<0) std::cout << y << " ";
    else std::cout << " " << y << " ";
    if (value<0) std::cout << value << " ";
    else std::cout << " " << value << " ";
    if (limit<0) std::cout << limit << std::endl;
    else std::cout << " " << limit << std::endl;
    if (verbose>1) {
      std::cout << "Press a key and type enter to continue. ";
      std::cin >> ch;
    }
    
    return success;
  }
  
  /** \brief Calculate the minimum \c min of \c func w.r.t 'x'.

      If this is not overloaded, it attempts to bracket the 
      minimum using bracket() and then calls min_bkt() with
      the newly bracketed minimum.
  */
  virtual int min(double &x, double &fmin, func_t &func)=0;
    
  /** \brief Calculate the minimum \c min of \c func with x2
      bracketed between x1 and x3
	  
      If this is not overloaded, it ignores the bracket and calls min().  
  */
  virtual int min_bkt(double &x2, double x1, double x3, double &fmin,
		      func_t &func)=0;
    
  /** \brief Calculate the minimum \c min of \c func with
      derivative \c dfunc w.r.t 'x'.
	
      If this is not overloaded, it attempts to bracket the 
      minimum using bracket() and then calls min_bkt_de() with
      the newly bracketed minimum.
  */
  virtual int min_de(double &x, double &fmin, func_t &func,
		     dfunc_t &df)=0;
  
  /** \brief Given interval <tt>(ax,bx)</tt>, attempt to bracket a
      minimum for function \c func.
	
      Upon success, <tt>fa=func(ax)</tt>, <tt>fb=func(bx)</tt>, and
      <tt>fc=func(cx)</tt> with <tt>fb<fa</tt>, <tt>fb<fc</tt> and
      <tt>ax<bx<cx</tt>. The initial values of \c cx, \c fa, 
      \c fb, and \c fc are all ignored.

      The number of iterations is controlled by \ref bracket_iter.

      \note This algorithm can fail if there's a minimum which has a
      much smaller size than \f$ bx-ax \f$, or if the function has the
      same value at \c ax, \c bx, and the midpoint <tt>(ax+bx)/2</tt>.
      
      \future Improve this algorithm with the golden ratio
      method in gsl/min/bracketing.c?
  */
  virtual int bracket(double &ax, double &bx, double &cx, double &fa, 
		      double &fb, double &fc, func_t &func) {
      
    double x=ax, x2=bx, x3=(ax+bx)/2.0;
    double fx, fx3, fx2;
    int i=0;
      
    bool done=false;
    while(done==false && i<bracket_iter) {
      fx=func(x);
      fx2=func(x2);
      fx3=func(x3);

      if (verbose>0) {
	std::cout << "Function min::bracket(), Iteration: " 
		  << i << std::endl;
	std::cout << " " << x << " " << x3 << " " << x2 << std::endl;
	std::cout << " " << fx << " " << fx3 << " " << fx2 << std::endl;
	if (verbose>1) {
	  char ch;
	  std::cout << "Press a key and type enter to continue. ";
	  std::cin >> ch;
	}
      }

      if (fx3>=fx2 && fx3>=fx) {
	// If the middle value is higher than the endpoints, 
	// try again with one side or the other
	if (fx2>fx) {
	  x2=x3;
	  x3=(x+x2)/2.0;
	} else {
	  x=x3;
	  x3=(x+x2)/2.0;
	}
      } else if (fx<=fx3 && fx3<=fx2) {
	// If we're increasing, extend to the left
	x3=x;
	x=x2-2.0*(x2-x);
      } else if (fx3<fx2 && fx3<fx) {
	// If we've succeeded, copy the results over
	done=true;
	ax=x;
	bx=x3;
	cx=x2;
	fa=fx;
	fb=fx3;
	fc=fx2;
      } else {
	// Otherwise we're decreasing, extend to the right
	x3=x2;
	x2=x+2.0*(x2-x);
      }
      i++;
    }
    
    if (done==false) {
      O2SCL_ERR("Too many iterations in min::bracket().",
		    exc_emaxiter);
    }
    
    return 0;
  }
  
  /// Return string denoting type ("min")
  virtual const char *type() { return "min"; }
  
  };

  /** \brief One-dimensional bracketing minimization [abstract base]
  */
  template<class func_t, class dfunc_t=func_t> 
    class min_bkt_base : public min_base<func_t,dfunc_t> {

  public:
    
  min_bkt_base() {
    bracket_iter=20;
  }
  
  virtual ~min_bkt_base() {}
      
  /** \brief The number of iterations for automatically 
      bracketing a minimum (default 20)
  */
  int bracket_iter;

  /** \brief Calculate the minimum \c min of \c func w.r.t 'x'.

      If this is not overloaded, it attempts to bracket the 
      minimum using bracket() and then calls min_bkt() with
      the newly bracketed minimum.
  */
  virtual int min(double &x, double &fmin, func_t &func) {
    double xl, xr, f, fl, fr;
    xl=x*0.9;
    xr=x*1.1;
    if (this->bracket(xl,x,xr,fl,f,fr,func)!=0) {
      O2SCL_CONV_RET("Failed to bracket in min().",exc_efailed,
			     this->err_nonconv);
    }
    return min_bkt(x,xl,xr,fmin,func);
  }
    
  /** \brief Calculate the minimum \c min of \c func with x2
      bracketed between x1 and x3
	  
      If this is not overloaded, it ignores the bracket and calls min().  
  */
  virtual int min_bkt(double &x2, double x1, double x3, double &fmin,
		      func_t &func)=0;
    
  /** \brief Calculate the minimum \c min of \c func with
      derivative \c dfunc w.r.t 'x'.
	
      If this is not overloaded, it attempts to bracket the 
      minimum using bracket() and then calls min_bkt_de() with
      the newly bracketed minimum.
  */
  virtual int min_de(double &x, double &fmin, func_t &func,
		     dfunc_t &df) {
    double xl, xr, f, fl, fr;
    xl=x*0.9;
    xr=x*1.1;
    if (this->bracket(xl,x,xr,fl,f,fr,func)!=0) {
      O2SCL_CONV_RET("Failed to bracket in min_de().",exc_efailed,
			     this->err_nonconv);
    }
    return min_bkt(x,xl,xr,fmin,func);
  }
  
  /// Return string denoting type ("min_bkt")
  virtual const char *type() { return "min_bkt"; }
  
  };

  /** \brief One-dimensional minimization using derivatives [abstract base]

      At the moment there are no minimizers of this type implemented in
      \o2. 

      \future Create a version of \ref o2scl::mmin_conf which
      implements a minimizer with this interface.
  */
  template<class func_t, class dfunc_t=func_t> 
    class min_de_base : public min_base<func_t,dfunc_t> {

  public:
    
  min_de_base() {
  }
      
  virtual ~min_de_base() {}
      
  /** \brief Calculate the minimum \c min of \c func w.r.t 'x'.

      If this is not overloaded, it attempts to bracket the 
      minimum using bracket() and then calls min_bkt() with
      the newly bracketed minimum.
  */
  virtual int min(double &x, double &fmin, func_t &func)=0;
    
  /** \brief Calculate the minimum \c min of \c func with x2
      bracketed between x1 and x3
	  
      If this is not overloaded, it ignores the bracket and calls min().  
  */
  virtual int min_bkt(double &x2, double x1, double x3, double &fmin,
		      func_t &func)=0;
  
  /** \brief Calculate the minimum \c min of \c func with
      derivative \c dfunc w.r.t 'x'.
      
      If this is not overloaded, it attempts to bracket the 
      minimum using bracket() and then calls min_bkt_de() with
      the newly bracketed minimum.
  */
  virtual int min_de(double &x, double &fmin, func_t &func,
		     dfunc_t &df)=0;
  
  /// Return string denoting type ("min_de")
  virtual const char *type() { return "min_de"; }
  
  };

  /** \brief Constrain \c x to be within \c width
      of the value given by \c center
      
      Defining \f$ c= \f$ \c center, \f$ w= \f$ \c width, \f$ h= \f$
      \c height, this returns the value \f$ h (1+|x-c-w|/w) \f$ if \f$
      x>c+w \f$ and \f$ h (1+|x-c+w|/w) \f$ if \f$ x<c-w \f$ . The
      value near \f$ x=c-w \f$ or \f$ x=c+w \f$ is \f$ h \f$ (the
      value of the function exactly at these points is zero) and the
      value at \f$ x=c-2w \f$ or \f$ x=c+2w \f$ is \f$ 2 h \f$ .

      It is important to note that, for large distances of \c x 
      from \c center, this only scales linearly. If you are trying to
      constrain a function which decreases more than linearly by
      making \c x far from \c center, then a minimizer may
      ignore this constraint.
  */
  inline double constraint(double x, double center, double width, 
			   double height) {
    double ret=0.0;
    if (x>center+width) {
      ret=height*(1.0+fabs(x-center-width)/width);
    } else if (x<center-width) {
      ret=height*(1.0+fabs(x-center+width)/width);
    }
    return ret;
  }

  /** \brief Constrain \c x to be within \c width of the value given
      by \c center

      Defining \f$ c= \f$ \c center, \f$ w= \f$ \c width, \f$ h= \f$
      \c height, \f$ t= \f$ \c tightness, and \f$ \ell= \f$ \c
      exp_arg_limit, this returns the value
      \f[
      h \left(\frac{x-c}{w}\right)^2 \left[ 
      1+ e^{t\left(x-c+w\right)\left(c+w-x\right)/w^2}
      \right]^{-1}
      \f]
      
      This function is continuous and differentiable.  Note that if
      \f$ x=c \f$ , then the function returns zero.

      The exponential is handled gracefully by assuming that anything
      smaller than \f$ \exp(-\ell) \f$ is zero. This creates a small
      discontinuity which can be removed with the sufficiently large
      value of \f$ \ell \f$. 
      
      It is important to note that, for large distances of \c x from
      \c center, this scales quadratically. If you are trying to
      constrain a function which decreases faster than quadratically
      by making \c x far from \c center, then a minimizer may ignore
      this constraint.
      
      In the limit \f$ t \rightarrow \infty \f$, this function
      converges towards the squared value of \ref constraint(), except
      exactly at the points \f$ x=c-w \f$ and \f$ x=c+w \f$.
  */
  inline double cont_constraint(double x, double center, double width, 
				double height, double tightness=40.0, 
				double exp_arg_limit=50.0) {
    double ret, wid2=width*width;
    double arg=tightness/wid2*(x-center+width)*(center+width-x);
    if (arg<-exp_arg_limit) {
      ret=(x-center)*(x-center)/wid2;
    } else {
      ret=(x-center)*(x-center)/wid2/(1.0+exp(arg));
    }
    return ret*height;
  }
  
  /** \brief Constrain \c x to be greater than the value given by \c
      center
      
      Defining \f$ c= \f$ \c center, \f$ w= \f$ \c width, \f$ h= \f$
      \c height, this returns \f$ h(1+|x-c|/w) \f$ if \f$ x<c \f$ and
      zero otherwise.  The value at \f$ x=c \f$ is \f$ h \f$ , while
      the value at \f$ x=c-w \f$ is \f$ 2 h \f$ .

      It is important to note that, for large distances of \c x from
      \c center, this only scales linearly. If you are trying to
      constrain a function which decreases more than linearly by
      making \c x far from \c center, then a minimizer may ignore this
      constraint.
  */
  inline double lower_bound(double x, double center, double width, 
			    double height) {
    double ret=0.0;
    if (x<center) ret=height*(1.0+fabs(x-center)/width);
    return ret;
  }
  
  /** \brief Constrain \c x to be greater than the value given by \c
      center 

      Defining \f$ c= \f$ \c center, \f$ w= \f$ \c width, \f$ h= \f$
      \c height, \f$ t= \f$ \c tightness, and \f$ \ell= \f$ \c
      exp_arg_limit, this returns \f$ h(c-x+w)/(w+w\exp(t(x-c)/w)) \f$
      and has the advantage of being a continuous and differentiable
      function.  The value of the function exactly at \f$ x=c \f$ is
      \f$ h/2 \f$, but for \f$ x \f$ just below \f$ c \f$ the function
      is \f$ h \f$ and just above \f$ c \f$ the function is quite
      small. 

      The exponential is handled gracefully by assuming that anything
      smaller than \f$ \exp(-\ell) \f$ is zero. This creates a small
      discontinuity which can be removed with the sufficiently large
      value of \f$ \ell \f$. 

      It is important to note that, for large distances of \c x 
      from \c center, this only scales linearly. If you are trying to
      constrain a function which decreases more than linearly by
      making \c x far from \c center, then a minimizer may
      ingore this constraint.

      In the limit \f$ t \rightarrow \infty \f$, this function
      converges towards \ref lower_bound(), except exactly at the
      point \f$ x=c \f$.
  */
  inline double cont_lower_bound(double x, double center, double width, 
				 double height, double tightness=40.0, 
				 double exp_arg_limit=50.0) {
    double ret, arg=tightness*(x-center)/width;
    if (arg>exp_arg_limit) {
      ret=0.0;
    } else if (arg<-exp_arg_limit) {
      ret=height*(center-x+width)/width;
    } else {
      ret=height*(center-x+width)/width/(1.0+exp(arg));
    }
    return ret;
  }

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
