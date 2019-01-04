/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Jerry Gagelman and Andrew W. Steiner
  
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
#include <cmath>

// AWS, 11/11/15: cstdlib appears to be
// required for size_t in gsl_sf_legendre
#include <cstdlib>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>

#include <o2scl/test_mgr.h>
#include <o2scl/inte_qaws_gsl.h>

using namespace std;
using namespace o2scl;

double legendre_gsl(double x, void *params) {
  int deg=*((int *)params);
  return gsl_sf_legendre_Pl(deg,x);
}

double legendre(double x, int &deg) {
  return gsl_sf_legendre_Pl(deg,x);
}

/* Test adaptive integration \ref o2scl::inte_qaws_gsl against the smooth
   integrand with algebraic-logarithmic endpoint-singularities,
   \f[
   I=\int_0^1 P_n(x) x^\alpha (1-x)^\beta \log^\mu(x) log^\nu(1-x)~dx
   \f]
   for several configurations of paramters \f$ \alpha,beta,\mu,\nu \f$
   and various tolerance values. 
   
   The test manager compares the result of each integration against that 
   obtained by the GSL routine \c gsl_integration_qag() with success
   only if the values differ by \f$ \epsilon_\mathrm{mach}. \f$ 
*/
int main(void) {

  test_mgr t;
  t.set_output_level(1);

  /* The degree of the Legendre polynomial in the integrand should be
   * even to avoid "resolving" the weight's singularity at the origin. 
   */
  int deg=24;
  size_t limit=512;
	
  inte_qaws_gsl<funct> Q;
  funct f=std::bind(legendre,std::placeholders::_1,deg);
	
  double alpha=0.0, beta=0.0; 
  int mu=0, nu=0;

  // setup GSL data
  gsl_integration_workspace*
    work=gsl_integration_workspace_alloc(limit);
  gsl_function f_gsl={&legendre_gsl,&deg};
  gsl_integration_qaws_table *	
    table=gsl_integration_qaws_table_alloc(alpha,beta,mu,nu);
	
  double o2scl_res, o2scl_err, gsl_res, gsl_err;
	
  // absolute error bound
  Q.tol_abs=1e-8;	
  // relative error bound
  Q.tol_rel=0;		
	
  while (Q.tol_abs > 1.5e-13) {
    cout << "\nabsolute tolerance: " << Q.tol_abs << "...\n";
    cout.width(6); cout << "alpha";
    cout.width(6); cout << "beta";
    cout.width(6); cout << "mu";
    cout.width(6); cout << "nu";
    cout.width(15); cout << "err. estimate";
    cout.width(15); cout << "iterations";
    cout.width(15); cout << "|O2scl - GSL|";
    cout << endl;
		
    alpha=-0.5; 
    beta=0; 
    mu=0; 
    nu=0;
    
    for(int config=1;config<=4;++config) {
      
      switch (config) {
      case 2:
	beta=-0.5;
	break;
      case 3:
	mu=1;
	break;
      case 4:
	nu=1;
	break;
      default:
	alpha=-0.5;
	beta=0; 
	mu=0; 
	nu=0;
      }
			
      cout.width(6); cout << alpha;
      cout.width(6); cout << beta;
      cout.width(6); cout << mu;
      cout.width(6); cout << nu;
			
      gsl_integration_qaws_table_set(table,alpha,beta,mu,nu);
      gsl_integration_qaws(&f_gsl,0,1,table,Q.tol_abs,Q.tol_rel,
			   limit,work,&gsl_res,&gsl_err);

      Q.set_weight(alpha,beta,mu,nu);
      Q.integ_err(f,0,1,o2scl_res,o2scl_err);

      double dbl_eps=std::numeric_limits<double>::epsilon();
      t.test_abs(gsl_res,o2scl_res,dbl_eps,"QAWS: O2scl vs GSL");
			
      cout.width(15); cout << o2scl_err;
      cout.width(15); cout << Q.last_iter;
      cout.width(15); cout << fabs(o2scl_res - gsl_res);
      cout << endl;
    }

    Q.tol_abs/=10.0;
  }
	
  gsl_integration_workspace_free(work);
  gsl_integration_qaws_table_free(table);

  t.report();

  return 0;
}
