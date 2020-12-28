/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#ifndef O2SCL_ODE_IV_TABLE_H
#define O2SCL_ODE_IV_TABLE_H

/** \file ode_iv_table.h
    \brief File defining \ref o2scl::ode_iv_table 
*/

#include <o2scl/astep.h>
#include <o2scl/astep_gsl.h>
#include <o2scl/ode_iv_solve.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Solve an initial-value ODE problem and store
      the result in a \ref table object

      This class is experimental.

      \future It would be nice not to have to copy the results from a
      matrix into a table, but this may require a nontrivial
      modification of the ODE solvers and/or the table class.

      \future One possible idea is to redo the name specification as
      a separate function, which allows one to either specify 
      prefixes or full column names. We also need to figure out
      how to handle column units.
  */
  template<class func_t=ode_funct<>, 
    class vec_t=ubvector, class alloc_vec_t=ubvector, 
    class alloc_t=ubvector_alloc> class ode_iv_table : 
  public ode_iv_solve<func_t,vec_t,alloc_vec_t,alloc_t> {
    
    public:

    /// Desc
    int solve_grid_table(size_t n, vec_t &ystart, table<> &t,
			 std::string x_col, std::string y_prefix,
			 std::string dydx_prefix, std::string yerr_prefix,
			 func_t &derivs) {
			 
      std::string cname;

      size_t n_sol=t.get_nlines();
      double x0=t.get(x_col,0);
      double x1=t.get(x_col,n_sol-1);
      double h=t.get(x_col,1)-x0;

      // Allocate vectors and matrices for solve_grid()
      alloc_vec_t x_sol;
      this->ao.allocate(x_sol,n_sol);
      ubmatrix ysol(n_sol,n), dydx_sol(n_sol,n), yerr_sol(n_sol,n);

      // Copy over x_grid from table
      ubvector &ob=t.get_column(x_col);
      vector_copy(n_sol,x_col,x_sol);

      // Perform solution
      this->template solve_grid<ubmatrix,ubmatrix_row>
      (x0,x1,h,n,ystart,n_sol,x_sol,ysol,dydx_sol,yerr_sol,derivs);

      // Pointers to columns in table
      std::vector<ubvector *> yt, dyt, errt;

      // Create new columns if necessary and get pointers
      for(size_t i=0;i<n;i++) {

	// TABLE FIXME

	/*
	cname=y_prefix+itos(i);
	if (!t.is_column(cname)) t.new_column(cname);
	yt.push_back(&t.get_column(cname));

	cname=dydx_prefix+itos(i);
	if (!t.is_column(cname)) t.new_column(cname);
	dyt.push_back(&t.get_column(cname));

	cname=yerr_prefix+itos(i);
	if (!t.is_column(cname)) t.new_column(cname);
	errt.push_back(&t.get_column(cname));

	// Now copy data over
	for(size_t j=0;j<n_sol;j++) {
	  (*yt[i])[j]=ysol[j][i];
	  (*dyt[i])[j]=dydx_sol[j][i];
	  (*errt[i])[j]=yerr_sol[j][i];
	}
	*/

      }

      return 0;
    }

    /// Desc
    int solve_store_table(double x0, double x1, double h, size_t n, 
			  vec_t &ystart, size_t &n_sol, table<> &t, 
			  std::string x_col, std::string y_prefix,
			  std::string dydx_prefix, std::string yerr_prefix, 
			  func_t &derivs) {

      std::string cname;
      if (t.get_nlines()<n_sol) t.set_nlines(n_sol);
      
      // Allocate vectors and matrices for solve_store()
      alloc_vec_t x_sol;
      this->ao.allocate(x_sol,n_sol);
      ubmatrix ysol(n_sol,n), dydx_sol(n_sol,n), yerr_sol(n_sol,n);

      // Perform solution
      this->template solve_store<ubmatrix,ubmatrix_row>
      (x0,x1,h,n,ystart,n_sol,x_sol,ysol,dydx_sol,yerr_sol,derivs);

      // Pointers to columns in table
      std::vector<ubvector *> yt, dyt, errt;

      if (!t.is_column(x_col)) t.new_column(x_col);
      ubvector &x_vec=t.get_column(x_col);

      // Create new columns if necessary and get pointers
      for(size_t i=0;i<n;i++) {

	// TABLE FIXME

	/*

	cname=y_prefix+itos(i);
	if (!t.is_column(cname)) t.new_column(cname);
	yt.push_back(&t.get_column(cname));

	cname=dydx_prefix+itos(i);
	if (!t.is_column(cname)) t.new_column(cname);
	dyt.push_back(&t.get_column(cname));

	cname=yerr_prefix+itos(i);
	if (!t.is_column(cname)) t.new_column(cname);
	errt.push_back(&t.get_column(cname));

	// Now copy data over
	for(size_t j=0;j<n_sol;j++) {
	  if (i==0) x_vec[j]=x_sol[j];
	  (*yt[i])[j]=ysol[j][i];
	  (*dyt[i])[j]=dydx_sol[j][i];
	  (*errt[i])[j]=yerr_sol[j][i];
	}

	*/

      }

      return 0;
    }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
