 /*
  -------------------------------------------------------------------
  
  Copyright (C) 2008-2022, Julien Garaud and Andrew W. Steiner
  
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
#ifndef O2SCL_ODE_BV_MSHOOT_H
#define O2SCL_ODE_BV_MSHOOT_H

/** \file ode_bv_mshoot.h
    \brief File defining \ref o2scl::ode_bv_mshoot 
*/

#include <o2scl/ode_bv_solve.h>

namespace o2scl {

  /** \brief Solve boundary-value ODE problems by multishooting
      with a generic nonlinear solver

      This class is experimental.

      Default template arguments
      - \c func_t - \ref ode_funct 
      - \c mat_t - \ref boost::numeric::ublas::matrix \< double \>
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
      - \c vec_int_t - \ref boost::numeric::ublas::vector \< int \>

      \future Make a class which performs an iterative linear
      solver which uses sparse matrices like ode_it_solve?
  */
  template<class func_t=ode_funct, 
    class mat_t=boost::numeric::ublas::matrix<double>,
    class vec_t=boost::numeric::ublas::vector<double>, 
    class vec_int_t=boost::numeric::ublas::vector<int> > 
    class ode_bv_mshoot : public ode_bv_solve {
    
    public:
    
    ode_bv_mshoot() {
      oisp=&def_ois;
      mrootp=&def_mroot;
      mem_size=0;
    }
    
    virtual ~ode_bv_mshoot() {
    }
    
    /** \brief Solve the boundary-value problem and store the 
	solution
    */
    int solve_final_value(double h, size_t n, size_t n_bound, 
			  vec_t &x_bound, mat_t &y_bound, 
			  vec_int_t &index, func_t &derivs) { 
      
      // Make copies of the inputs for later access
      this->l_index=&index;
      this->l_nbound=n_bound;
      this->l_xbound=&x_bound;
      this->l_ybound=&y_bound;
      this->l_h=h;
      this->l_derivs=&derivs;
      this->l_n=n;
  
      int lhs_unks=0, rhs_conts=0;
      this->l_lhs_unks=lhs_unks;

      // Count the number of variables we'll need to solve for
      for (size_t i=0;i<n;i++) {
	if (index[i]<2) lhs_unks++;
	if (index[i]%2==1) rhs_conts++;
      }

      // Make sure that the boundary conditions make sense
      if (lhs_unks!=rhs_conts) {
	O2SCL_ERR2("Incorrect boundary conditions in ",
		       "ode_bv_mshoot::solve()",gsl_einval);
      } 

      // The number of variables to solve for
      int n_solve=lhs_unks+n*(n_bound-2);

      // Make space for the solution
      ubvector tx(n_solve);
      
      // Copy initial guess from y_bound
      size_t j=0;
      for (size_t i=0;i<n;i++) {
	if (index[i]<2) {
	  tx[j]=y_bound[0][i];
	  j++;
	}
      }
      for(size_t k=1;k<n_bound-1;k++) {
	for(size_t i=0;i<n;i++) {
	  tx[j]=y_bound[k][i];
	  j++;
	}
      }

      // Allocate memory as needed
      if (n!=mem_size) {
	sy.resize(n);
	sy2.resize(n);
	syerr.resize(n);
	sdydx.resize(n);
	sdydx2.resize(n);
	mem_size=n;
      }
  
      // The function object
      mm_funct_mfptr<ode_bv_mshoot<func_t,mat_t,vec_t,vec_int_t> >
      mfm(this,&ode_bv_mshoot<func_t,mat_t,vec_t,vec_int_t>::solve_fun);
      
      // Solve
      int ret=this->mrootp->msolve(n_solve,tx,mfm);
      if (ret!=0) {
	O2SCL_ERR("Solver failed in ode_bv_mshoot::solve().",ret);
      }

      return ret;
    }
    
    /** \brief Solve the boundary-value problem and store the 
	solution
    */
    template<class mat_row_t> 
    int solve_store(double h, size_t n, size_t n_bound, vec_t &x_bound, 
		    mat_t &y_bound, vec_int_t &index, size_t &n_sol, 
		    vec_t &x_sol, mat_t &y_sol, mat_t &dydx_sol, 
		    mat_t &yerr_sol, func_t &derivs) {

      if (n_bound<2) {
	O2SCL_ERR2("Not enough boundaries (must be at least two) in ",
		       "ode_bv_mshoot::solve_store().",gsl_einval);
      }
      if (n_sol<n_bound) {
	O2SCL_ERR2("Not enough room to store boundaries in ",
		       "ode_bv_mshoot::solve_store().",gsl_einval);
      }

      // Solve to fix x_bound and y_bound
      solve_final_value(h,n,n_bound,x_bound,y_bound,index,derivs);

      // Find the indices which correspond to the boundaries
      ubvector_int inxs(n_bound);
      for(size_t k=0;k<n_bound;k++) {
	inxs[k]=((size_t)(((double)k)/((double)n_bound)*
			  (((double)(n_sol))-1.0+1.0e-12)));
	std::cout << k << " " << inxs[k] << " " << n_sol << std::endl;
      }
      // Double check that each interval has some space 
      for(size_t k=1;k<n_bound-1;k++) {
	if (inxs[k]==inxs[k-1] || inxs[k]==inxs[k+1]) {
	  O2SCL_ERR2("Not enough room to store boundaries in ",
			 "ode_bv_mshoot::solve_store().",gsl_einval);
	}
      }
	
      // Now create the table
      for(size_t k=0;k<n_bound-1;k++) {
	size_t n_sol_tmp=inxs[k+1];
	mat_row_t ystart(y_bound,k);

	std::cout << "Old boundaries: " << inxs << std::endl;

	oisp->template solve_store<mat_t,mat_row_t>
	(x_bound[k],x_bound[k+1],h,n,ystart,
	 n_sol_tmp,x_sol,y_sol,dydx_sol,yerr_sol,derivs,inxs[k]);

	std::cout << "New boundaries: " << n_sol_tmp << std::endl;

	// If it didn't use all the space, shift the indexes
	// accordingly
	if (((int)n_sol_tmp)<inxs[k+1]) {
	  for(size_t k2=k+1;k2<n_bound;k2++) {
	    inxs[k2]=((size_t)(((double)k2-k-1)/((double)(n_bound-k-1))*
			       (((double)(n_sol-n_sol_tmp)))))+n_sol_tmp-1;
	  }
	}

	std::cout << "New boundaries: " << inxs << std::endl;
      }
      
      n_sol=inxs[n_bound-1];

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
    vec_t *l_xbound;

    /// Storage for the ending vector
    mat_t *l_ybound;

    /// Storage for the stepsize
    double l_h;

    /// The functions to integrate
    func_t *l_derivs;

    /// The number of functions
    size_t l_n;

    /// The number of boundaries
    size_t l_nbound;

    /// The number of unknowns on the LHS
    size_t l_lhs_unks;

    /// \name Temporary storage for \ref solve_fun()
    //@{
    vec_t sy, sy2, syerr, sdydx, sdydx2;
    //@}
    
    /// Size of recent allocation
    size_t mem_size;

    /// The shooting function to be solved by the multidimensional solver
    int solve_fun(size_t nv, const vec_t &tx, vec_t &ty) {

      // Counters to index through function parameters tx and ty
      size_t j_x=0, j_y=0;
      
      // Shorthand for the boundary specifications
      vec_t &xb=*this->l_xbound;
      mat_t &yb=*this->l_ybound;

      // Set up the boundaries from the values in tx
      for(size_t k=0;k<l_nbound-1;k++) {
	if (k==0) {
	  for(size_t i=0;i<this->l_n;i++) {
	    if ((*this->l_index)[i]<2) {
	      yb[k][i]=tx[j_x];
	      j_x++;
	    }
	  }
	} else {
	  for(size_t i=0;i<this->l_n;i++) {
	    yb[k][i]=tx[j_x];
	    j_x++;
	  }
	}
      }

      // Integrate between all of the boundaries
      for(size_t k=0;k<l_nbound-1;k++) {

	double x0=xb[k];
	double x1=xb[k+1];

	// Setup the start vector sy. This extra copying from l_ybound
	// to sy might be avoided by using a mat_row_t object to send
	// l_ybound directly to the solver?

	for(size_t i=0;i<this->l_n;i++) {
	  sy[i]=yb[k][i];
	}

	// Shoot across. We specify memory for the derivatives and
	// errors to prevent ode_iv_solve from having to continuously
	// allocate and reallocate it.
	oisp->solve_final_value(x0,x1,this->l_h,this->l_n,sy,
				sy2,*this->l_derivs);

	if (k!=l_nbound-2) {
	  // Construct equations for the internal RHS boundaries
	  for(size_t i=0;i<this->l_n;i++) {
	    ty[j_y]=sy2[i]-yb[k+1][i];
	    j_y++;
	  }
	} else {
	  // Construct the equations for the rightmost boundary
	  for(size_t i=0;i<this->l_n;i++) {
	    if ((*this->l_index)[i]%2==1) {
	      double yright=yb[k+1][i];
	      if (yright==0.0) {
		ty[j_y]=sy2[i];
	      } else {
		ty[j_y]=(sy2[i]-yright)/yright;
	      }
	      j_y++;
	    }
	  }
	}
	
	// End of loop to integrate between all boundaries
      }

      return 0;
    }

#endif

  };

}

#endif
