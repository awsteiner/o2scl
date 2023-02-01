/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifndef O2SCL_ODE_BV_SOLVE_H
#define O2SCL_ODE_BV_SOLVE_H

/** \file ode_bv_solve.h
    \brief File defining \ref o2scl::ode_bv_solve 
*/

#include <string>
#include <o2scl/astep.h>
#include <o2scl/astep_gsl.h>
#include <o2scl/ode_iv_solve.h>
#include <o2scl/mroot_hybrids.h>

namespace o2scl {

  /** \brief Base class for boundary-value ODE solvers

      This class is experimental.
  */
  class ode_bv_solve {

    public:

    ode_bv_solve() {
      verbose=0;
    }

    virtual ~ode_bv_solve() {}
    
    /** \name Values for index arrays */
    //@{
    /// Unknown on both the left and right boundaries
    static const int unk=0;
    /// Known on the right boundary
    static const int right=1;
    /// Known on the left boundary
    static const int left=2;
    /// Known on both the left and right boundaries
    static const int both=3;
    //@}

    /// Set output level
    int verbose;
    
  };
  
  /** \brief Solve boundary-value ODE problems by shooting from one 
      boundary to the other

      This class is experimental.

      Documentation links for default template arguments
      - \c func_t - \ref ode_funct 
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
      - \c vec_int_t - \ref boost::numeric::ublas::vector \< int \>
  */
  template<class func_t=ode_funct, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class vec_int_t=boost::numeric::ublas::vector<int> > 
    class ode_bv_shoot : public ode_bv_solve {
    
    public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    
    ode_bv_shoot() {
      oisp=&def_ois;
      mrootp=&def_mroot;
      mem_size=0;
    }
    
    virtual ~ode_bv_shoot() {
    }

    /// Allocate internal storage
    void allocate(size_t n) {
      if (n!=mem_size) {
	sy.resize(n);
	sy2.resize(n);
	syerr.resize(n);
	sdydx.resize(n);
	mem_size=n;
      }
      return;
    }
    
    /** \brief Solve the boundary-value problem and store the 
	solution

	Given the \c n initial values of the functions in \c ystart,
	this function integrates the ODEs specified in \c derivs over
	the interval from \c x0 to \c x1 with an initial stepsize of
	\c h. The final values of the function are given in \c yend,
	the derivatives in \c dydx_end, and the associated errors are
	given in \c yerr. The initial values of \c yend and \c yerr
	are ignored.

    */
    int solve_final_value(double x0, double x1, double h, size_t n,
			  vec_t &ystart, vec_t &yend, vec_int_t &index,
			  vec_t &yerr, vec_t &dydx_end, func_t &derivs) {
      
      // Get pointers to inputs for access by solve function
      this->l_index=&index;
      this->l_ystart=&ystart;
      this->l_yend=&yend;
      this->l_yerr=&yerr;
      this->l_dydx_end=&dydx_end;
      this->l_x0=x0;
      this->l_x1=x1;
      this->l_h=h;
      this->l_derivs=&derivs;
      this->l_n=n;
  
      // lhs_unks is the number of unknowns on the LHS
      // rhs_conts is the number of constraints on the RHS
      int lhs_unks=0, rhs_conts=0;

      // Count the number of variables we'll need to solve for
      for (size_t i=0;i<n;i++) {
	if (index[i]<2) lhs_unks++;
	if (index[i]%2==1) rhs_conts++;
      }

      // Make sure that the boundary conditions make sense
      if (lhs_unks!=rhs_conts) {
	O2SCL_ERR2("Incorrect boundary conditions in ",
		       "ode_bv_shoot::solve()",gsl_einval);
      } 

      // Make space for the solution
      ubvector tx(lhs_unks);

      // Copy initial guesses from ystart
      int j=0;
      for (size_t i=0;i<n;i++) {
	if (index[i]<2) {
	  tx[j]=ystart[i];
	  j++;
	}
      }

      allocate(n);
  
      // The function object
      mm_funct_mfptr<ode_bv_shoot<func_t,vec_t,vec_int_t> > 
      mfm(this,&ode_bv_shoot<func_t,vec_t,vec_int_t>::solve_fun);

      // Solve
      int ret=this->mrootp->msolve(lhs_unks,tx,mfm);
      if (ret!=0) {
	O2SCL_ERR("Solver failed in ode_bv_shoot::solve().",ret);
      }

      // Copy the solution back to ystart
      j=0;
      for (size_t i=0;i<n;i++) {
	if (index[i]<2) {
	  ystart[i]=tx[j];
	  j++;
	}
      }

      return 0;
    }

    /** \brief Solve the boundary-value problem and store the 
	solution
     */
    template<class mat_t, class mat_row_t> 
    int solve_store(double x0, double x1, double h, size_t n,
		    vec_t &ystart, vec_t &yend, 
		    vec_int_t &index, size_t &n_sol, vec_t &x_sol, 
		    mat_t &y_sol, mat_t &yerr_sol, mat_t &dydx_sol, 
		    func_t &derivs) {

      mat_row_t yerr(yerr_sol,0);
      mat_row_t dydx_end(dydx_sol,0);

      // Solve for the final value
      solve_final_value(x0,x1,h,n,ystart,yend,index,yerr,dydx_end,derivs);

      // Copy ystart to y_sol[0]
      for(size_t i=0;i<n;i++) {
	y_sol[0][i]=ystart[i];
      }

      // Evaluate the function one more time to create the table
      oisp->solve_store<mat_t,mat_row_t>(this->l_x0,this->l_x1,this->l_h,
					 this->l_n,n_sol,x_sol,
					 y_sol,yerr_sol,dydx_sol,derivs);

      // Copy the stored solution back to ystart and yend
      for (size_t i=0;i<n;i++) {
	if (index[i]%2==0) {
	  yend[i]=y_sol[n_sol-1][i];
	}
	ystart[i]=y_sol[0][i];
      }
  
      return 0;
    }
    
    /** \brief Set initial value solver
     */
    int set_iv(ode_iv_solve<func_t,vec_t> &ois) {
      oisp=&ois;
      return 0;
    }
    
    /** \brief Set the equation solver */
    int set_mroot(mroot<mm_funct<> > &root) {
      mrootp=&root;
      return 0;
    }

    /// The default initial value solver
    ode_iv_solve<func_t,vec_t> def_ois;

    /// The default equation solver
    gsl_mroot_hybrids<mm_funct<> > def_mroot;

#ifndef DOXYGEN_INTERNAL

    protected:
    
    /// The solver for the initial value problem
    ode_iv_solve<func_t,vec_t> *oisp;

    /// The equation solver
    mroot<mm_funct<> > *mrootp;

    /// The index defining the boundary conditions
    vec_int_t *l_index;

    /// Storage for the starting vector
    vec_t *l_ystart;

    /// Storage for the ending vector
    vec_t *l_yend;

    /// Storage for the starting vector
    vec_t *l_yerr;

    /// Storage for the ending vector
    vec_t *l_dydx_end;

    /// Storage for the starting point
    double l_x0;

    /// Storage for the ending abcissa
    double l_x1;

    /// Storage for the stepsize
    double l_h;

    /// The functions to integrate
    func_t *l_derivs;

    /// The number of functions
    size_t l_n;

    /// \name Temporary storage for \ref solve_fun()
    //@{
    vec_t sy, sy2, syerr, sdydx;
    //@}
    
    /// Size of recent allocation
    size_t mem_size;

    /// The shooting function to be solved by the multidimensional solver
    int solve_fun(size_t nv, const vec_t &tx, vec_t &ty) {

      int j;

      // Create leftmost point from combination of 
      // starting array and proposed solution
      j=0;
      for(size_t i=0;i<this->l_n;i++) {
	if ((*this->l_index)[i]<2) {
	  sy[i]=tx[j];
	  j++;
	} else {
	  sy[i]=(*this->l_ystart)[i];
	}
      }
  
      // Shoot across. We specify memory for the derivatives and
      // errors to prevent ode_iv_solve from having to continuously
      // allocate and reallocate it.
      oisp->solve_final_value(this->l_x0,this->l_x1,this->l_h,
			      this->l_n,sy,sy2,*this->l_yerr,
			      *this->l_dydx_end,*this->l_derivs);
  
      j=0;
      for(size_t i=0;i<this->l_n;i++) {
	if ((*this->l_index)[i]%2==1) {
	  // Construct the equations from the rightmost point
	  if ((*this->l_yend)[i]==0.0) {
	    ty[j]=sy2[i];
	  } else {
	    ty[j]=(sy2[i]-(*this->l_yend)[i])/(*this->l_yend)[i];
	  }
	  j++;
	} else {
	  // Otherwise copy the final values from y2 to *l_yend
	  (*this->l_yend)[i]=sy2[i];
	}
      }
  
      return 0;
    }

#endif

  };

  /** \brief Solve boundary-value ODE problems by shooting from one 
      boundary to the other on a grid

      This class is experimental.

      Default template arguments
      - \c func_t - \ref ode_funct
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
      - \c vec_int_t - \ref boost::numeric::ublas::vector \< int \>
  */
  template<class mat_t=boost::numeric::ublas::matrix<double>, 
    class mat_row_t=boost::numeric::ublas::matrix_row<
    boost::numeric::ublas::matrix<double> >, 
    class func_t=ode_funct<>, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class vec_int_t=boost::numeric::ublas::vector<int> > 
    class ode_bv_shoot_grid : 
    public ode_bv_shoot<func_t,vec_t,vec_int_t> {
    
    public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<int> ubvector_int;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    ode_bv_shoot_grid() {
    }
    
    virtual ~ode_bv_shoot_grid() {
    }
    
    /// Desc
    int solve_grid(double x0, double x1, double h, size_t n, 
		   vec_t &ystart, vec_t &yend, vec_int_t &index,
		   size_t nsol, vec_t &xsol, mat_t &ysol, 
		   mat_t &err_sol, mat_t &dydx_sol, func_t &derivs) {

      // lhs_unks is the number of unknowns on the LHS
      // rhs_conts is the number of constraints on the RHS
      int lhs_unks=0, rhs_conts=0;

      // Count the number of variables we'll need to solve for
      for (size_t i=0;i<n;i++) {
	if (index[i]<2) lhs_unks++;
	if (index[i]%2==1) rhs_conts++;
      }

      // Make sure that the boundary conditions make sense
      if (lhs_unks!=rhs_conts) {
	O2SCL_ERR2("Incorrect boundary conditions in ",
		       "ode_bv_shoot_grid::solve_grid()",gsl_einval);
      } 

      // Get pointers to inputs for access by solve function
      this->l_n=n;
      this->l_index=&index;
      this->l_derivs=&derivs;
      this->l_ystart=&ystart;
      this->l_yend=&yend;
      this->l_h=h;
      l_nsol=nsol;
      l_xsol=&xsol;
      l_ysol=&ysol;
      l_dydxsol=&dydx_sol;
      l_errsol=&err_sol;

      // Make space for the solution
      ubvector tx(lhs_unks);

      // Copy initial guesses from ysol
      int j=0;
      for (size_t i=0;i<n;i++) {
	if (index[i]<2) {
	  tx[j]=ystart[i];
	  j++;
	}
      }

      // Allocate internal storage
      this->allocate(n);

      // The function object
      mm_funct_mfptr<ode_bv_shoot_grid<mat_t,mat_row_t,
      func_t,vec_t,vec_int_t> > 
      mfm(this,&ode_bv_shoot_grid<mat_t,mat_row_t,func_t,vec_t,
	  vec_int_t>::solve_grid_fun);
  
      // Solve
      int ret=this->mrootp->msolve(lhs_unks,tx,mfm);
      if (ret!=0) {
	O2SCL_ERR("Solver failed in ode_bv_shoot_grid::solve_grid().",ret);
      }

      // Copy the solution back to ysol
      j=0;
      for (size_t i=0;i<n;i++) {
	if (index[i]<2) {
	  ysol[0][i]=tx[j];
	  j++;
	}
      }

      return 0;
    }

#ifndef DOXYGEN_INTERNAL

    protected:
    
    /// Desc
    size_t l_nsol;

    /// Desc
    vec_t *l_xsol;

    /// Desc
    mat_t *l_ysol;

    /// Desc
    mat_t *l_dydxsol;

    /// Desc
    mat_t *l_errsol;

    /// The shooting function to be solved by the multidimensional solver    
    int solve_grid_fun(size_t nv, const vec_t &tx, vec_t &ty) {

      int j;

      // Create leftmost point from combination of 
      // starting array and proposed solution
      j=0;
      for(size_t i=0;i<this->l_n;i++) {
	if ((*this->l_index)[i]<2) {
	  (*l_ysol)[0][i]=tx[j];
	  j++;
	}
      }

      // Perform the solution
      this->oisp->template solve_grid<mat_t,mat_row_t>
	(this->l_h,this->l_n,l_nsol,*l_xsol,*l_ysol,*l_errsol,
	 *l_dydxsol,*this->l_derivs);

      // The last vector
      mat_row_t yend2(*l_ysol,l_nsol-1);

      j=0;
      for(size_t i=0;i<this->l_n;i++) {
	if ((*this->l_index)[i]%2==1) {
	  // Construct the equations from the rightmost point
	  if ((*this->l_yend)[i]==0.0) {
	    ty[j]=yend2[i];
	  } else {
	    ty[j]=((*this->l_yend)[i]-yend2[i])/(*this->l_yend)[i];
	  }
	  j++;
	}
      }

      return 0;
    }

#endif

  };

}

#endif
