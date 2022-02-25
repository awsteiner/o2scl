/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_ODE_IV_SOLVE_H
#define O2SCL_ODE_IV_SOLVE_H

/** \file ode_iv_solve.h
    \brief File defining \ref o2scl::ode_iv_solve 
*/

#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/astep.h>
#include <o2scl/astep_gsl.h>

namespace o2scl {

  /** \brief Solve an initial-value ODE problems given an adaptive ODE
      stepper

      This class gives several functions which solve an initial
      value ODE problem. The functions \ref solve_final_value()
      gives only the final value of the functions at the end
      of the ODE integration and is relatively fast. 

      The function solve_store() stores the solution of the ODE over
      the full range into a set of vectors and matrices which are
      allocated and specified by the user. This function is designed
      to give exactly the same results (though this cannot be
      guaranteed) as solve_final_value() and additionally records some
      or all of the results from the adaptive steps which were taken.

      All of these functions automatically evaluate the derivatives
      from the specified function at the initial point and
      user-specified initial derivatives are ignored. The total
      number of steps taken is limited by \ref ntrial and \ref
      nsteps stores the number of steps taken by the most recent
      solution. The variable \ref nsteps_out is the maximum number
      of points in the interval for which verbose output will be
      given when \ref o2scl::ode_iv_solve::verbose is greater than zero.

      \verbatim embed:rst
      There is an example for the usage of this class in
      ``examples/ex_ode.cpp<`` documented in the
      :ref:`Ordinary differential equations example`.
      \endverbatim

      <b>Convergence error handling</b>
      
      There are two different convergence errors which can 
      be controlled separately in this class. 
      - The adaptive stepper may require too many steps. If this
      happens, then the solver immediately stops. The solver
      calls the error handler if \ref err_nonconv is true, and
      otherwise it returns a non-zero value.
      - The adaptive stepper may fail. If \ref exit_on_fail
      is true, then the error handler is called. Otherwise,
      the solver proceeds to continue computing the whole solution.
      So long as the number of adaptive steps required is less
      than \ref ntrial, then the full solution is computed and
      a non-zero value is returned to indicate the accuracy of 
      the solution may be impacted. If the number of adaptive
      steps required after a failure of the adaptive stepper
      is larger than \ref ntrial, then the behavior of the
      solver is controlled by \ref err_nonconv as described
      above.

      Documentation links for default template arguments
      - \c func_t - \ref ode_funct
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
      
      The default adaptive stepper is an object of type \ref astep_gsl.

      \future The form of solve_final_value() is very similar to that
      of astep_base::astep_full(), but not quite the same. Maybe
      these functions should be consistent with each other?
  */
  template<class func_t=ode_funct, 
	   class vec_t=boost::numeric::ublas::vector<double> > 
  class ode_iv_solve {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
   
#ifndef DOXYGEN_INTERNAL
    
  protected:

    /// \name Vectors for temporary storage 
    //@{
    vec_t vtemp, vtemp2, vtemp3, vtemp4;
    //@}

    /// The size of the temporary vectors
    size_t mem_size;

    /// The adaptive stepper
    astep_base<vec_t,vec_t,vec_t,func_t> *astp;
    
    /// Print out iteration information
    virtual int print_iter(double x, size_t nv, vec_t &y) {
      std::cout << type() << " x: " << x << " y: ";
      for(size_t i=0;i<nv;i++) std::cout << y[i] << " ";
      std::cout << std::endl;
      if (verbose>1) {
	char ch;
	std::cin >> ch;
      }
      return 0;
    }

    /// Free allocated memory
    void free() {
      if (mem_size>0) {
	vtemp.clear();
	vtemp2.clear();
	vtemp3.clear();
	vtemp4.clear();
      }
    }

    /** \brief Allocate space for temporary vectors
     */
    void allocate(size_t n) {
      if (n!=mem_size) {
	free();
	vtemp.resize(n);
	vtemp2.resize(n);
	vtemp3.resize(n);
	vtemp4.resize(n);
	mem_size=n;
      }
    }
  
#endif

  public:
      
    ode_iv_solve() {
      verbose=0;
      ntrial=1000;
      nsteps_out=10;
      astp=&gsl_astp;
      exit_on_fail=true;
      mem_size=0;
      err_nonconv=true;
    }
      
    virtual ~ode_iv_solve() {
      free();
    }

    /** \brief If true, call the error handler if the solution does 
	not converge (default true)
    */
    bool err_nonconv;
  
    /// \name Main solver functions
    //@{
    /** \brief Solve the initial-value problem to get the final value

	Given the \c n initial values of the functions in \c ystart,
	this function integrates the ODEs specified in \c derivs over
	the interval from \c x0 to \c x1 with an initial stepsize of \c
	h. The final values of the function are given in \c yend and the
	initial values of \c yend are ignored.

	If \ref verbose is greater than zero, The solution at less than
	or approximately equal to \ref nsteps_out points will be written
	to \c std::cout. If \ref verbose is greater than one, a
	character will be required after each selected point.
    */
    int solve_final_value(double x0, double x1, double h, size_t n,
			  vec_t &ystart, vec_t &yend, func_t &derivs) {

      allocate(n);
      return solve_final_value(x0,x1,h,n,ystart,yend,vtemp2,
			       vtemp3,derivs);
    }
      
    /** \brief Solve the initial-value problem to get the final value
	with errors

	Given the \c n initial values of the functions in \c ystart,
	this function integrates the ODEs specified in \c derivs over
	the interval from \c x0 to \c x1 with an initial stepsize of \c
	h. The final values of the function are given in \c yend and the
	associated errors are given in \c yerr. The initial values of \c
	yend and \c yerr are ignored.

	If \ref verbose is greater than zero, The solution at less
	than or approximately equal to \ref nsteps_out points will be
	written to \c std::cout. If \ref verbose is greater than one,
	a character will be required after each selected point.
    */
    int solve_final_value(double x0, double x1, double h, size_t n,
			  vec_t &ystart, vec_t &yend, vec_t &yerr,
			  func_t &derivs) {

      allocate(n);
      return solve_final_value(x0,x1,h,n,ystart,yend,yerr,
			       vtemp2,derivs);
    }
      
    /** \brief Solve the initial-value problem to get the final value,
	derivatives, and errors

	Given the \c n initial values of the functions in \c ystart,
	this function integrates the ODEs specified in \c derivs over
	the interval from \c x0 to \c x1 with an initial stepsize of \c
	h. The final values of the function are given in \c yend, the
	derivatives in \c dydx_end, and the associated errors are given
	in \c yerr. The initial values of \c yend and \c yerr are
	ignored.

	This function is designed to be relatively fast,
	avoiding extra copying of vectors back and forth. 

	If \ref verbose is greater than zero, The solution at less
	than or approximately equal to \ref nsteps_out points will be
	written to \c std::cout. If \ref verbose is greater than one,
	a character will be required after each selected point.

	This function computes \c dydx_start automatically and the
	values given by the user are ignored. 

	The solution fails if more than \ref ntrial steps are required.
	This function will also fail if <tt>x1>x0</tt> and <tt>h<0</tt>
	or if <tt>x1-x0</tt> and <tt>h>0</tt> do not have the same sign.
    */
    int solve_final_value(double x0, double x1, double h, size_t n,
			  vec_t &ystart, vec_t &yend, vec_t &yerr, 
			  vec_t &dydx_end, func_t &derivs) {
    
      if ((x1>x0 && h<=0.0) || (x0>x1 && h>=0.0)) {
	std::string str="Interval direction (x1-x0="+o2scl::dtos(x1-x0)+
	  ") does not match step direction (h="+o2scl::dtos(h)+
	  " in ode_iv_solve::solve_final_value().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }
      if (x0==x1) {
	O2SCL_ERR2("Starting and final endpoints identical in ",
		   "ode_iv_solve::solve_final_value().",exc_einval);
      }

      // Allocate for temporary vectors
      allocate(n);
    
      // The variable 'x' is the current independent variable, xnext is 
      // the next value of x, xverb is the next value of x for which 
      // verbose output should be provided, and dxverb is the stepsize 
      // for xverb.
      double x=x0, xverb=0.0, dxverb=0.0, xnext;
      int ret=0, first_ret=0;

      nsteps=0;

      // Create a reference to make the code easier to read
      vec_t &dydx_start=vtemp;

      // Compute stepsize for verbose output
      if (verbose>0) {
	print_iter(x0,n,ystart);
	if (verbose>1) {
	  char ch;
	  std::cin >> ch;
	}
	// Ensure that 'dxverb' is positive
	if (x1>x0) {
	  dxverb=(x1-x0)/((double)nsteps_out);
	  xverb=x0+dxverb;
	} else {
	  dxverb=(x0-x1)/((double)nsteps_out);
	  xverb=x0-dxverb;
	}
      }

      // Use yend as workspace
      for(size_t i=0;i<n;i++) yend[i]=ystart[i];

      // Compute initial derivative
      derivs(x,n,yend,dydx_start);
      
      bool done=false;
      while (done==false) {
	
	// Take a step of the first type
	ret=astp->astep_full(x,x1,xnext,h,n,ystart,dydx_start,
			     yend,yerr,dydx_end,derivs);
      
	if (ret!=0) {
	  if (exit_on_fail) {
	    O2SCL_ERR2("Adaptive stepper failed in ",
		       "ode_iv_solve::solve_final_value()",ret);
	  } else if (first_ret!=0) {
	    first_ret=ret;
	  }
	}

	// Print out verbose info 
	if (verbose>0 && xnext>xverb) {
	  print_iter(xnext,n,yend);
	  if (verbose>1) {
	    char ch;
	    std::cin >> ch;
	  }
	  xverb+=dxverb;
	}

	// Check number of iterations
	nsteps++;
	if (nsteps>ntrial) {
	  std::string str="Too many steps required (nsteps="+
	    szttos(nsteps)+", ntrial="+szttos(ntrial)+
	    ", x="+o2scl::dtos(x)+") in ode_iv_solve::solve_final_value().";
	  O2SCL_CONV_RET(str.c_str(),exc_emaxiter,err_nonconv);
	}

	if (ret!=0) {
	  done=true;
	} else {
	  if (x1>x0) {
	    if (xnext>=x1) done=true;
	  } else {
	    if (xnext<=x1) done=true;
	  }
	}

	if (done==false) {

	  // Take a step of the second type
	  ret=astp->astep_full(xnext,x1,x,h,n,yend,dydx_end,ystart,yerr,
			       dydx_start,derivs);
	
	  if (ret!=0) {
	    if (exit_on_fail) {
	      O2SCL_ERR2("Adaptive stepper failed in ",
			 "ode_iv_solve::solve_final_value()",ret);
	    } else if (first_ret!=0) {
	      first_ret=ret;
	    }
	  }

	  // Print out verbose info 
	  if (verbose>0 && x>xverb) {
	    print_iter(x,n,ystart);
	    if (verbose>1) {
	      char ch;
	      std::cin >> ch;
	    }
	    xverb+=dxverb;
	  }

	  // Check number of iterations
	  nsteps++;
	  if (nsteps>ntrial) {
	    std::string str="Too many steps required (ntrial="+itos(ntrial)+
	      ", x="+o2scl::dtos(x)+") in ode_iv_solve::solve_final_value().";
	    O2SCL_ERR(str.c_str(),exc_emaxiter);
	  }

	  if (ret!=0) {
	    done=true;
	  } else {
	    if (x1>x0) {
	      if (x>=x1) {
		done=true;
		for(size_t j=0;j<n;j++) {
		  yend[j]=ystart[j];
		  dydx_end[j]=dydx_start[j];
		}
	      }
	    } else {
	      if (x<=x1) {
		done=true;
		for(size_t j=0;j<n;j++) {
		  yend[j]=ystart[j];
		  dydx_end[j]=dydx_start[j];
		}
	      }
	    }
	  }

	}

	// End of while loop
      }
	
      // Print out final verbose info
      if (verbose>0) {
	print_iter(x1,n,yend);
	if (verbose>1) {
	  char ch;
	  std::cin >> ch;
	}
      }
	
      return first_ret;
    }
      
    /** \brief Solve the initial-value problem and store the 
	associated output

	Initially, \c x_sol should be a vector of size \c n_sol, and \c
	y_sol, \c dydx_sol, and \c yerr_sol should be matrices with size
	\c [n_sol][n].  On exit, \c n_sol will will be number of points
	store, less than or equal to the original value of \c
	n_sol. This function avoids performing extra calls to the
	adaptive stepper, and the solution will be approximately evenly
	spaced.

	This function is also designed to give the exactly the same
	results as solve_final_value(). This cannot be strictly
	guaranteed, but will likely hold in most applications.

	This template function works with any matrix class \c mat_t
	which can be accessed using successive applications of
	operator[] and which has an associated class \c mat_row_t
	which returns a row of a matrix of type \c mat_t as
	an object with type \c vec_t. 

	If \ref verbose is greater than zero, The solution at each
	internal point will be written to \c std::cout.  If \ref verbose
	is greater than one, a character will be required after each
	point.

	\todo Document \c istart
    */
    template<class mat_t>
    int solve_store(double x0, double x1, double h, size_t n, size_t &n_sol, 
		    vec_t &x_sol, mat_t &y_sol, mat_t &yerr_sol, 
		    mat_t &dydx_sol, func_t &derivs, size_t istart=0) {
    
      int ret=0;
      int first_ret=0;
      size_t nmax=n_sol-1;
      nsteps=0;

      // Stepsize for next verbose output. Use nsteps_out-1 instead of
      // nsteps_out since the first point is always output below.
      double dx_verb=(x1-x0)/((double)(nsteps_out-1));
      // Stepsize for next point for storage
      double dx_tab=(x1-x0)/((double)(n_sol-istart-1));

      double x_verb=x0+dx_verb;
      double x_tab=x0+dx_tab;
    
      // Allocate space for errors and derivatives and extra storage
      allocate(n);

      // Create some references just to make the code easier
      // to read
      vec_t &ystart=vtemp4;
      vec_t &dydx_start=vtemp2;

      // Copy first point to ystart
      for(size_t j=0;j<n;j++) ystart[j]=y_sol(istart,j);

      // Output first point 
      if (verbose>0) {
	print_iter(x0,n,ystart);
	if (verbose>1) {
	  char ch;
	  std::cin >> ch;
	}
      }

      // Initial derivative evaulation
      derivs(x0,n,ystart,dydx_start);

      // Add first derivatives to storage, and set the errors to zero
      x_sol[istart]=x0;
      for(size_t j=0;j<n;j++) {
	dydx_sol(istart,j)=dydx_start[j];
	yerr_sol(istart,j)=0.0;
      }

      // Copy first point to storage again for first step
      size_t icurr=istart+1;
      x_sol[icurr]=x0;
      for(size_t j=0;j<n;j++) y_sol(icurr,j)=ystart[j];

      double xnext;
    
      bool done=false;
      while (done==false) {

	// Create some references just to make the code easier
	// to read
	vec_t &yerr=vtemp;
	vec_t &dydx_out=vtemp3;

	// Use ystart as temporary storage for the end of the current
	// adaptive step
	vec_t yrow(n);
	for(size_t i=0;i<n;i++) yrow[i]=y_sol(icurr,i);
	ret=astp->astep_full(x0,x1,xnext,h,n,yrow,dydx_start,ystart,yerr,
			     dydx_out,derivs);
      
	if (ret!=0) {
	  if (exit_on_fail) {
	    n_sol=icurr+1;
	    // call error handler
	    O2SCL_ERR2("Adaptive stepper returned non-zero in ",
		       "ode_iv_solve::solve_store().",exc_efailed);
	  } else if (first_ret==0) {
	    first_ret=ret;
	  }
	}
      
	// Update step counter and abort if necessary
	nsteps++;
	if (nsteps>ntrial) {
	  std::string str="Too many steps required (ntrial="+itos(ntrial)+
	    ", x="+o2scl::dtos(x0)+") in ode_iv_solve::solve_store().";
	  O2SCL_ERR(str.c_str(),exc_emaxiter);
	}
      
	// If we've made enough progress, do verbose output
	if (xnext>=x_verb && verbose>0) {
	  print_iter(xnext,n,ystart);
	  if (verbose>1) {
	    char ch;
	    std::cin >> ch;
	  }
	  x_verb+=dx_verb;
	}

	// If we're at the end
	if (xnext>=x1) {

	  // Exit the loop
	  done=true;

	  // Store the final entry
	  x_sol[icurr]=xnext;
	  for(size_t j=0;j<n;j++) {
	    y_sol(icurr,j)=ystart[j];
	    dydx_sol(icurr,j)=dydx_out[j];
	    yerr_sol(icurr,j)=yerr[j];
	  }

	  // Update the solution size
	  n_sol=icurr+1;

	} else {

	  if (xnext>=x_tab && icurr<nmax) {
	  
	    // If we've made enough progress, store an entry in the table
	    x_sol[icurr]=xnext;
	    for(size_t j=0;j<n;j++) {
	      y_sol(icurr,j)=ystart[j];
	      dydx_sol(icurr,j)=dydx_out[j];
	      yerr_sol(icurr,j)=yerr[j];
	    
	      // Also prepare for the next step
	      y_sol(icurr+1,j)=ystart[j];
	      dydx_start[j]=dydx_out[j];
	    }
	    x_tab+=dx_tab;
	  
	    // Update x0 and the current and next indicies
	    x0=xnext;
	    icurr++;
	  
	  } else {

	    // Otherwise, just prepare for the next step
	    // without updating icurr
	    x0=xnext;
	    for(size_t j=0;j<n;j++) {
	      dydx_start[j]=dydx_out[j];
	      y_sol(icurr,j)=ystart[j];
	    }
	  
	  }
	
	}
     
	// End of loop 'while (done==false)'
      }

      return first_ret;
    }
    //@}

    /// Set output level
    int verbose;
  
    /** \brief Number of output points for verbose output (default 10)
      
	This is used in functions solve_store() and solve_final_value()
	to control how often steps are output when verbose is greater
	than zero.
    */
    size_t nsteps_out;

    /// Maximum number of applications of the adaptive stepper (default 1000)
    size_t ntrial;

    /// Number of adaptive ste!ps employed
    size_t nsteps;

    /// \name The adaptive stepper
    //@{
    /// Set the adaptive stepper to use
    int set_astep(astep_base<vec_t,vec_t,vec_t,func_t> &as) {
      astp=&as;
      return 0;
    }

    /** \brief If true, stop the solution if the adaptive stepper fails
	(default true)
    */
    bool exit_on_fail;

    /// The default adaptive stepper
    astep_gsl<vec_t,vec_t,vec_t,func_t> gsl_astp;
    //@}

    /// Return the type, \c "ode_iv_solve".
    virtual const char *type() { return "ode_iv_solve"; }
  
  };

  typedef matrix_row_gen_ctor<boost::numeric::ublas::matrix<double> >
    solve_grid_mat_row;

  typedef std::function<int(double,size_t,const solve_grid_mat_row &,
			    solve_grid_mat_row &)>
    ode_funct_solve_grid;
 
  /** \brief Solve an initial-value ODE problems on a grid 
      given an adaptive ODE stepper
      
      This class works as similar to ode_iv_solve::solve_store()
      except that the solution is stored on a grid of points in the
      independent variable specified by the user, at the cost of
      taking extra steps to ensure that function values, derivatives,
      and errors are computed at each grid point.

      \verbatim embed:rst
      There is an example for the usage of this class in
      ``examples/ex_ode.cpp<`` documented in the
      :ref:`Ordinary differential equations example`.
      \endverbatim
   */
  template<class func_t=ode_funct_solve_grid,
    class mat_row_t=solve_grid_mat_row> 
    class ode_iv_solve_grid {
    
#ifndef DOXYGEN_INTERNAL
    
  protected:

    /// The adaptive stepper
    astep_base<mat_row_t,mat_row_t,mat_row_t,func_t> *astp;
    
    /// Print out iteration information
    virtual int print_iter(double x, size_t nv, mat_row_t &y) {
      std::cout << type() << " x: " << x << " y: ";
      for(size_t i=0;i<nv;i++) std::cout << y[i] << " ";
      std::cout << std::endl;
      if (verbose>1) {
	char ch;
	std::cin >> ch;
      }
      return 0;
    }

#endif

  public:
      
    ode_iv_solve_grid() {
      verbose=0;
      ntrial=1000;
      astp=&gsl_astp;
      exit_on_fail=true;
      err_nonconv=true;
    }
      
    virtual ~ode_iv_solve_grid() {
    }

    /** \brief If true, call the error handler if the solution does 
	not converge (default true)
    */
    bool err_nonconv;
  
    /// \name Main solver function
    //@{
    /** \brief Solve the initial-value problem from \c x0 to \c x1 over
	a grid storing derivatives and errors

	Initially, \c xsol should be an array of size \c nsol, and \c
	ysol should be a \c ubmatrix of size \c [nsol][n]. This function
	never takes a step larger than the grid size.

	If \ref verbose is greater than zero, The solution at each grid
	point will be written to \c std::cout. If \ref verbose is
	greater than one, a character will be required after each point.

	\future Consider making a version of grid which gives the 
	same answers as solve_final_value(). After each proposed step,
	it would go back and fill in the grid points if the proposed
	step was past the next grid point.
    */
    template<class vec_t, class mat_t>
    int solve_grid(double h, size_t n, size_t nsol, vec_t &xsol, 
		   mat_t &ysol, mat_t &err_sol, mat_t &dydx_sol, 
		   func_t &derivs) {
    
      double x0=xsol[0];
      double x1=xsol[nsol-1];

      double x=x0, xnext;
      int ret=0, first_ret=0;
      nsteps=0;

      // Copy the initial point to the first row
      xsol[0]=x0;

      mat_row_t y_start(n);
      mat_row_t dydx_start(n);

      // Copy ysol[0] to ystart
      for(size_t j=0;j<n;j++) {
	y_start[j]=ysol(0,j);
      }
    
      // Verbose output for the first row
      if (verbose>0) print_iter(xsol[0],n,y_start);

      // Evaluate the first derivative
      derivs(x0,n,y_start,dydx_start);

      // Set the derivatives and the uncertainties for the first row
      for(size_t j=0;j<n;j++) {
	dydx_sol(0,j)=dydx_start[j];
	err_sol(0,j)=0.0;
      }

      for(size_t i=1;i<nsol && ret==0;i++) {
      
	mat_row_t y_row(ysol,i);
	mat_row_t dydx_row(dydx_sol,i);
	mat_row_t yerr_row(err_sol,i);

	// Step until we reach the next grid point
	bool done=false;
	while(done==false && ret==0) {
	
	  ret=astp->astep_full(x,xsol[i],xnext,h,n,y_start,dydx_start,
			       y_row,yerr_row,dydx_row,derivs);

	  nsteps++;
	  if (ret!=0) {
	    if (exit_on_fail) {
	      O2SCL_ERR2("Adaptive stepper failed in ",
			 "ode_iv_solve_grid::solve_grid()",ret);
	    } else if (first_ret!=0) {
	      first_ret=ret;
	    }
	  }
	
	  if (nsteps>ntrial) {
	    std::string str="Too many steps required (ntrial="+itos(ntrial)+
	      ", x="+o2scl::dtos(x)+") in ode_iv_solve_grid::solve_grid().";
	    O2SCL_ERR(str.c_str(),exc_emaxiter);
	  }

	  // Adjust independent variable for next step
	  x=xnext;
	  for(size_t j=0;j<n;j++) {
	    y_start[j]=y_row[j];
	    dydx_start[j]=dydx_row[j];
	  }
	
	  // Decide if we have reached the grid point
	  if (x1>x0) {
	    if (x>=xsol[i]) done=true;
	  } else {
	    if (x<=xsol[i]) done=true;
	  }
	    
	}
      
	if (verbose>0) print_iter(xsol[i],n,y_row);

      }
    
      return first_ret;
    }
    //@}

    /// Set output level
    int verbose;
  
    /// Maximum number of applications of the adaptive stepper (default 1000)
    size_t ntrial;

    /// Number of adaptive steps employed
    size_t nsteps;

    /// \name The adaptive stepper
    //@{
    /// Set the adaptive stepper to use
    int set_astep(astep_base<mat_row_t,mat_row_t,mat_row_t,func_t> &as) {
      astp=&as;
      return 0;
    }

    /** \brief If true, stop the solution if the adaptive stepper fails
	(default true)
    */
    bool exit_on_fail;

    /// The default adaptive stepper
    astep_gsl<mat_row_t,mat_row_t,mat_row_t,func_t> gsl_astp;
    //@}

    /// Return the type, \c "ode_iv_solve".
    virtual const char *type() { return "ode_iv_solve_grid"; }
  
  };

}

#endif
