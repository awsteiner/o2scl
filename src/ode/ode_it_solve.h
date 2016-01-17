/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#ifndef O2SCL_ODE_IT_SOLVE_H
#define O2SCL_ODE_IT_SOLVE_H

/** \file ode_it_solve.h
    \brief File defining \ref o2scl::ode_it_solve 
*/

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/misc.h>
#include <o2scl/test_mgr.h>
#include <o2scl/linear_solver.h>
#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /// Function for iterative solving of ODEs
  typedef std::function<double
    (size_t,double,boost::numeric::ublas::matrix_row
     <boost::numeric::ublas::matrix<double> > &)> ode_it_funct11;
  
  /// Function derivatives for iterative solving of ODEs
  typedef std::function<double
    (size_t,size_t,double,boost::numeric::ublas::matrix_row
     <boost::numeric::ublas::matrix<double> > &)> ode_it_dfunct11;
  
  /** \brief ODE solver using a generic linear solver to solve 
      finite-difference equations

      \future Set up convergence error if it goes beyond max iterations
      \future Create a GSL-like set() and iterate() interface
      \future Implement as a child of ode_bv_solve ?
      \future Max and average tolerance?
      \future Allow the user to ensure that the solver doesn't
      apply the full correction
  */
  template <class func_t=ode_it_funct11,
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double>, 
    class matrix_row_t=boost::numeric::ublas::matrix_row
    <boost::numeric::ublas::matrix<double> >, 
    class solver_vec_t=boost::numeric::ublas::vector<double> , 
    class solver_mat_t=boost::numeric::ublas::matrix<double> >
    class ode_it_solve {
    
  public:

  bool make_mats;
  
  ode_it_solve() {
    h=1.0e-4;
    niter=30;
    tol_rel=1.0e-8;
    verbose=0;
    solver=&def_solver;
    alpha=1.0;
    make_mats=false;
    fd=0;
    fl=0;
    fr=0;
  }

  virtual ~ode_it_solve() {}

  /// Set level of output (default 0)
  int verbose;

  /** \brief Stepsize for finite differencing (default \f$ 10^{-4} \f$)
   */
  double h;

  /// Tolerance (default \f$ 10^{-8} \f$)
  double tol_rel;
  
  /// Maximum number of iterations (default 30)
  size_t niter;
  
  /// Set the linear solver
  int set_solver(o2scl_linalg::linear_solver<solver_vec_t,solver_mat_t> &ls) {
    solver=&ls;
    return 0;
  }
    
  /// Size of correction to apply (default 1.0)
  double alpha;

  /** \brief Solve \c derivs with boundary conditions \c left and 
      \c right

      Given a grid of size \c n_grid and \c n_eq differential equations,
      solve them by relaxation. The grid is specified in \c x, which
      is a vector of size \c n_grid. The differential equations are
      given in \c derivs, the boundary conditions on the left hand
      side in \c left, and the boundary conditions on the right hand
      side in \c right. The number of boundary conditions on the left
      hand side is \c nb_left, and the number of boundary conditions on
      the right hand side should be <tt>n_eq-nb_left</tt>. The initial
      guess for the solution, a matrix of size <tt>[n_grid][n_eq]</tt>
      should be given in \c y. Upon success, \c y will contain an
      approximate solution of the differential equations. The matrix
      \c mat is workspace of size <tt>[n_grid*n_eq][n_grid*n_eq]</tt>, and
      the vectors \c rhs and \c y are workspace of size
      <tt>[n_grid*n_eq]</tt>.
  */
  int solve(size_t n_grid, size_t n_eq, size_t nb_left, vec_t &x, 
	    mat_t &y, func_t &derivs, func_t &left, func_t &right,
	    solver_mat_t &mat, solver_vec_t &rhs, solver_vec_t &dy) {

    // Store the functions for simple derivatives
    fd=&derivs;
    fl=&left;
    fr=&right;
    
    ode_it_dfunct11 d2_derivs=std::bind
      (std::mem_fn<double(size_t,size_t,double,matrix_row_t &)>
       (&ode_it_solve::fd_derivs),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);
    ode_it_dfunct11 d2_left=std::bind
      (std::mem_fn<double(size_t,size_t,double,matrix_row_t &)>
       (&ode_it_solve::fd_left),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);
    ode_it_dfunct11 d2_right=std::bind
      (std::mem_fn<double(size_t,size_t,double,matrix_row_t &)>
       (&ode_it_solve::fd_right),this,std::placeholders::_1,
       std::placeholders::_2,std::placeholders::_3,std::placeholders::_4);

    return solve_derivs(n_grid,n_eq,nb_left,x,y,derivs,left,right,
			d2_derivs,d2_left,d2_right,mat,rhs,dy);
  }

  /** \brief Solve \c derivs with boundary conditions \c left and 
      \c right

      Given a grid of size \c n_grid and \c n_eq differential equations,
      solve them by relaxation. The grid is specified in \c x, which
      is a vector of size \c n_grid. The differential equations are
      given in \c derivs, the boundary conditions on the left hand
      side in \c left, and the boundary conditions on the right hand
      side in \c right. The number of boundary conditions on the left
      hand side is \c nb_left, and the number of boundary conditions on
      the right hand side should be <tt>n_eq-nb_left</tt>. The initial
      guess for the solution, a matrix of size <tt>[n_grid][n_eq]</tt>
      should be given in \c y. Upon success, \c y will contain an
      approximate solution of the differential equations. The matrix
      \c mat is workspace of size <tt>[n_grid*n_eq][n_grid*n_eq]</tt>, and
      the vectors \c rhs and \c y are workspace of size
      <tt>[n_grid*n_eq]</tt>.
  */
  template<class dfunc_t>
  int solve_derivs(size_t n_grid, size_t n_eq, size_t nb_left, vec_t &x, 
		   mat_t &y, func_t &derivs, func_t &left, func_t &right,
		   dfunc_t &d_derivs, dfunc_t &d_left, dfunc_t &d_right,
		   solver_mat_t &mat, solver_vec_t &rhs, solver_vec_t &dy) {

    // Variable index
    size_t ix;

    // Number of RHS boundary conditions
    size_t nb_right=n_eq-nb_left;

    // Number of variables
    size_t nvars=n_grid*n_eq;
    
    bool done=false;
    for(size_t it=0;done==false && it<niter;it++) {
      
      ix=0;
      
      for(size_t i=0;i<nvars;i++) {
	for(size_t j=0;j<nvars;j++) {
	  mat(i,j)=0.0;
	}
      }

      // Construct the entries corresponding to the LHS boundary. 
      // This makes the first nb_left rows of the matrix.
      for(size_t i=0;i<nb_left;i++) {
	matrix_row_t yk=o2scl::matrix_row<mat_t,matrix_row_t>(y,0);
	rhs[ix]=-left(i,x[0],yk);
	for(size_t j=0;j<n_eq;j++) {
	  mat(ix,j)=d_left(i,j,x[0],yk);
	}
	ix++;
      }

      // Construct the matrix entries for the internal points
      // This loop adds n_grid-1 sets of n_eq rows
      for(size_t k=0;k<n_grid-1;k++) {
	size_t kp1=k+1;
	double tx=(x[kp1]+x[k])/2.0;
	double dx=x[kp1]-x[k];
	matrix_row_t yk=o2scl::matrix_row<mat_t,matrix_row_t>(y,k);
	matrix_row_t ykp1=o2scl::matrix_row<mat_t,matrix_row_t>(y,k+1);
	
	for(size_t i=0;i<n_eq;i++) {
	  
	  rhs[ix]=y(k,i)-y(kp1,i)+(x[kp1]-x[k])*
	    (derivs(i,tx,ykp1)+derivs(i,tx,yk))/2.0;
	  
	  size_t lhs=k*n_eq;
	  for(size_t j=0;j<n_eq;j++) {
	    mat(ix,lhs+j)=-d_derivs(i,j,tx,yk)*dx/2.0;
	    mat(ix,lhs+j+n_eq)=-d_derivs(i,j,tx,ykp1)*dx/2.0;
	    if (i==j) {
	      mat(ix,lhs+j)=mat(ix,lhs+j)-1.0;
	      mat(ix,lhs+j+n_eq)=mat(ix,lhs+j+n_eq)+1.0;
	    }
	  }

	  ix++;

	}
	
      }
      
      // Construct the entries corresponding to the RHS boundary
      // This makes the last nb_right rows of the matrix.
      for(size_t i=0;i<nb_right;i++) {
	matrix_row_t ylast=o2scl::matrix_row<mat_t,matrix_row_t>(y,n_grid-1);
	size_t lhs=n_eq*(n_grid-1);
	
	rhs[ix]=-right(i,x[n_grid-1],ylast);
	
	for(size_t j=0;j<n_eq;j++) {
	  mat(ix,lhs+j)=d_right(i,j,x[n_grid-1],ylast);
	}
	  
	ix++;

      }
      
      // Compute correction by calling the linear solver

      if (verbose>3) {
	std::cout << "Matrix: " << std::endl;
	for(size_t i=0;i<nvars;i++) {
	  for(size_t j=0;j<nvars;j++) {
	    std::cout << mat(i,j) << " ";
	  }
	  std::cout << std::endl;
	}
	std::cout << "Deviations:" << std::endl;
	for(size_t i=0;i<nvars;i++) {
	  std::cout << rhs[i] << std::endl;
	}
      }

      if (make_mats) return 0;

      solver->solve(ix,mat,rhs,dy);

      if (verbose>3) {
	std::cout << "Corrections:" << std::endl;
	for(size_t i=0;i<nvars;i++) {
	  std::cout << dy[i] << std::endl;
	}
	std::cout << "Press a key and press enter to continue: " 
		  << std::endl;
	char ch;
	std::cin >> ch;
      }
      
      // Apply correction and compute its size

      double res=0.0;
      ix=0;

      for(size_t igrid=0;igrid<n_grid;igrid++) {
	for(size_t ieq=0;ieq<n_eq;ieq++) {
	  y(igrid,ieq)+=alpha*dy[ix];
	  res+=dy[ix]*dy[ix];
	  ix++;
	}
      }

      if (verbose>0) {
	// Since we're in the o2scl namespace, we explicitly
	// specify std::sqrt() here
	std::cout << "ode_it_solve: " << it << " " << std::sqrt(res) << " " 
		  << tol_rel << std::endl;
	if (verbose>1) {
	  char ch;
	  std::cout << "Press a key and type enter to continue. ";
	  std::cin >> ch;
	}
      }
      
      // If the correction has become small enough, we're done
      if (std::sqrt(res)<=tol_rel) done=true;
    }

    if (done==false) {
      O2SCL_ERR("Exceeded number of iterations in solve().",
		    o2scl::exc_emaxiter);
    }

    return 0;
  }
  
  /// Default linear solver
  o2scl_linalg::linear_solver_HH<solver_vec_t,solver_mat_t> def_solver;
  
  protected:
  
  /// \name Storage for functions
  //@{
  func_t *fl, *fr, *fd;
  //@}

  /// Solver
  o2scl_linalg::linear_solver<solver_vec_t,solver_mat_t> *solver;

  /** \brief Compute the derivatives of the LHS boundary conditions

      This function computes \f$ \partial f_{left,\mathrm{ieq}} / \partial
      y_{\mathrm{ivar}} \f$
  */
  virtual double fd_left(size_t ieq, size_t ivar, double x, matrix_row_t &y) {

    double ret, dydx;
    
    y[ivar]+=h;
    ret=(*fl)(ieq,x,y);
    
    y[ivar]-=h;
    ret-=(*fl)(ieq,x,y);
    
    ret/=h;
    return ret;
  }
  
  /** \brief Compute the derivatives of the RHS boundary conditions
	
      This function computes \f$ \partial f_{right,\mathrm{ieq}} / \partial
      y_{\mathrm{ivar}} \f$
  */
  virtual double fd_right(size_t ieq, size_t ivar, double x, matrix_row_t &y) {

    double ret, dydx;
    
    y[ivar]+=h;
    ret=(*fr)(ieq,x,y);
    
    y[ivar]-=h;
    ret-=(*fr)(ieq,x,y);
    
    ret/=h;
    return ret;
  }
  
  /** \brief Compute the finite-differenced part of the 
      differential equations

      This function computes \f$ \partial f_{\mathrm{ieq}} / \partial
      y_{\mathrm{ivar}} \f$
  */
  virtual double fd_derivs(size_t ieq, size_t ivar, double x, matrix_row_t &y) {

    double ret, dydx;
    
    y[ivar]+=h;
    ret=(*fd)(ieq,x,y);
    
    y[ivar]-=h;
    ret-=(*fd)(ieq,x,y);
    
    ret/=h;
    
    return ret;
  }
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
