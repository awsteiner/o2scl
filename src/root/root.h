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
/** \file root.h
    \brief File for one-dimensional solver base class \ref o2scl::root
*/
#ifndef O2SCL_ROOT_H
#define O2SCL_ROOT_H

#include <iostream>
#include <cmath>
#include <o2scl/err_hnd.h>
#include <o2scl/funct.h>
#include <o2scl/misc.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief One-dimensional solver [abstract base]
      
      See the \ref onedsolve_subsect section of the User's guide for
      general information about \o2 solvers. 

      \future Maybe consider allowing the user to specify
      the stream to which 'verbose' information is sent.
  */
  template<class func_t=funct11, class dfunc_t=func_t> class root {
    
  public:
  
  root() {
    ntrial=100;
    tol_rel=1.0e-8;
    verbose=0;
    tol_abs=1.0e-12;
    last_ntrial=0;
    err_nonconv=true;
  }

  virtual ~root() {}

  /** \brief The maximum value of the functions for success 
      (default \f$ 10^{-8} \f$ )
  */
  double tol_rel;
      
  /** \brief The minimum allowable stepsize 
      (default \f$ 10^{-12} \f$ )
  */
  double tol_abs;
    
  /// Output control (default 0)
  int verbose;
    
  /// Maximum number of iterations (default 100)
  int ntrial;

  /** \brief If true, call the error handler if the solver does 
      not converge (default true)
  */
  bool err_nonconv;
    
  /// The number of iterations used in the most recent solve
  int last_ntrial;

  /// Return the type, \c "root".
  virtual const char *type() { return "root"; }

  /** \brief Print out iteration information.
         
      Depending on the value of the variable verbose, this prints
      out the iteration information. If verbose=0, then no
      information is printed, while if verbose>1, then after each
      iteration, the present values of \c x and \c y are
      output to std::cout along with the iteration number. If
      verbose>=2 then each iteration waits for a character before
      continuing.
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
    
  /** \brief Solve \c func using \c x as an initial guess
  */
  virtual int solve(double &x, func_t &func)=0;

  /** \brief Solve \c func in region \f$ x_1<x<x_2 \f$  
      returning \f$ x_1 \f$ .
  */
  virtual int solve_bkt(double &x1, double x2, func_t &func) {
    return solve(x1,func);
  }

  /** \brief Solve \c func using \c x as an initial
      guess using derivatives \c df.
  */
  virtual int solve_de(double &x, func_t &func, dfunc_t &df) {
    return solve(x,func);
  }

  };

  /** \brief One-dimensional bracketing solver [abstract base]
      
  */
  template<class func_t=funct11, class dfunc_t=func_t> class root_bkt :
  public root<func_t,dfunc_t> {

  public:

  root_bkt() {
    bracket_step=0.0;
    bracket_iters=10;
  }
      
  virtual ~root_bkt() {}

  /** \brief The stepsize for automatic bracketing 
      (default \f$ 10^{-4} \f$)

      If this is exactly zero, it will be reset to
      \f$ 10^{-4} \f$ by solve().
  */
  double bracket_step;

  /// The number of iterations in attempt to bracket root (default 10)
  size_t bracket_iters;
 
  /// Return the type, \c "root_bkt".
  virtual const char *type() { return "root_bkt"; }

  /** \brief Solve \c func in region \f$ x_1<x<x_2 \f$  
      returning \f$ x_1 \f$ .
  */
  virtual int solve_bkt(double &x1, double x2, func_t &func)=0; 
    
  /** \brief Solve \c func using \c x as an initial guess
  */
  virtual int solve(double &x, func_t &func) {
    
    if (bracket_step<=0.0) {
      bracket_step=fabs(x)*1.0e-4;
      if (bracket_step<1.0e-15) bracket_step=1.0e-15;
    }

    double x2=0.0, dx, fx, fx2;
    size_t i=0;
    bool done=false;
    // Use function to try to bracket a root
    while(done==false && i<bracket_iters) {
      
      fx=func(x);
      fx2=func(x*(1.0+bracket_step));

      if (this->verbose>0) {
	std::cout << "root_bkt::solve(): Iteration " << i << " of "
		  << bracket_iters << "." << std::endl;
	std::cout << "\tx1: " << x << " f(x1): " << fx << std::endl;
	std::cout << "\tx2: " << x*(1.0+bracket_step) 
		  << " f(x2): " << fx2 << std::endl;
      }
      
      dx=(fx2-fx)/(bracket_step*x);
      if (dx==0.0) {
	O2SCL_CONV_RET("Failed to bracket (dx==0) in root_bkt::solve().",
		       o2scl::exc_emaxiter,this->err_nonconv);
      }
      x2=x-2.0*fx/dx;
      
      fx2=func(x2);
      if (this->verbose>0) {
	std::cout << "\tx_new: " << x2 << " f(x_new): " << fx2 << std::endl;
	if (this->verbose>1) {
	  std::cout << "Press a key and type enter to continue. ";
	  char ch;
	  std::cin >> ch;
	}
      }
	  
      if (fx*fx2<0.0) {
	done=true;
      } else {
	x=(x2+x)/2.0;
      }

      i++;
    }
    if (done==false) {
      O2SCL_CONV_RET("Failed to bracket (iters>max) in root_bkt::solve().",
		     o2scl::exc_emaxiter,this->err_nonconv);
    }
    if (this->verbose>0) {
      std::cout << "root_bkt::solve(): Going to solve_bkt()." 
		<< std::endl;
    }
    return solve_bkt(x,x2,func);
  }

  /** \brief Solve \c func using \c x as an initial
      guess using derivatives \c df.
  */
  virtual int solve_de(double &x, func_t &func, dfunc_t &df) {

    double x2=0.0, dx, fx, fx2;
    int i=0;
    bool done=false;
	
    // Use derivative information to try to bracket a root
    while(done==false && i<10) {
	  
      fx=func(x);
      dx=df(x);
	  
      x2=x-2.0*fx/dx;
	  
      fx2=func(x2);
	    
      if (fx*fx2<0.0) {
	done=true;
      } else {
	x=(x2+x)/2.0;
      }
      i++;
    }
	
    if (done==false) {
      O2SCL_CONV_RET("Failed to bracket function in root_bkt::solve_de().",
		     o2scl::exc_emaxiter,this->err_nonconv);
    }
    
    return solve_bkt(x,x2,func);
  }
  
  };

  /** \brief One-dimensional with solver with derivatives [abstract base]

      \todo At the moment, the functions solve() and solve_bkt() 
      are not implemented for derivative solvers.
  */
  template<class func_t=funct11, class dfunc_t=func_t> 
    class root_de : public root<func_t> {
    
  public:
  
  root_de() {
  }
  
  virtual ~root_de() {}
  
  /// Return the type, \c "root_de".
  virtual const char *type() { return "root_de"; }

  /** \brief Solve \c func in region \f$ x_1<x<x_2 \f$  
      returning \f$ x_1 \f$ .
  */
  virtual int solve_bkt(double &x1, double x2, func_t &func) {
    O2SCL_ERR("Function solve_bkt() not implemented.",exc_eunimpl);
    return 0;
  }

  /** \brief Solve \c func using \c x as an initial guess
  */
  virtual int solve(double &x, func_t &func) {
    O2SCL_ERR("Function solve() not implemented.",exc_eunimpl);
    return 0;
  }

  /** \brief Solve \c func using \c x as an initial
      guess using derivatives \c df.
  */
  virtual int solve_de(double &x, func_t &func, dfunc_t &df)=0;

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

