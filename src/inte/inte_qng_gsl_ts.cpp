/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Jerry Gagelman
  and Andrew W. Steiner
  
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
/** Test inte_qng_gsl
    
    Test QNG routine against standard polynomials:
    \f[
    \int_{-1}^1 P_n^2(x)~dx = 2/(2n+1)
    \f]
    for the Legendre polynomials \f$ P_n \f$ for 
    \f$ n \in [0,10,20,\ldots,60] \f$ and for tolerances
    between \f$ 10^{-8} \f$ and \f$ 10^{-13} \f$.
*/
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/inte_qng_gsl.h>

using namespace std;
using namespace o2scl;

double legendre_gsl(double x, void *params) {
  int deg = *(int *)params;
  double P = gsl_sf_legendre_Pl(deg, x);
  return P*P;
}

double legendre(double x, int &deg) {
  double P = gsl_sf_legendre_Pl(deg, x);
  return P*P;
}

double legendre_exact(int deg) {
  return 2.0/(2*deg+1.0);
}

int main(void) {
  test_mgr test;
  test.set_output_level(1);

  // order parameter
  int n;	

  inte_qng_gsl<funct11> Q;

  // setup GSL data
  size_t neval;
  gsl_function f_gsl = {&legendre_gsl, &n};
	
  double o2scl_res, o2scl_err, gsl_res, gsl_err;

  // absolute error bound
  Q.tol_abs = 1e-7;		
  // relative error bound
  Q.tol_rel = 0;		
  
  while (Q.tol_abs > 5.0e-14) {

    cout << "\nAbsolute tolerance: " << Q.tol_abs << "...\n";
    cout.width(10); cout << "degree";
    cout.width(15); cout << "err. estimate";
    cout.width(15); cout << "exact error";
    cout.width(15); cout << "GK-rule used" << endl;
		
    n = 0;
    
    for(size_t i=0;i<7;i++) {
      
      funct11 f=std::bind(legendre,std::placeholders::_1,n);
      //funct_fptr_param<int> f(legendre,n);
      
      Q.integ_err(f,-1.0,1.0,o2scl_res,o2scl_err);
    
      gsl_integration_qng(&f_gsl,-1,1,Q.tol_abs,Q.tol_rel,
			  &gsl_res,&gsl_err,&neval);
      
      test.test_abs(o2scl_res,gsl_res,GSL_DBL_MIN,
		    "QNG: O2scl vs GSL");
			
      cout.width(10); cout << 2*n;
      cout.width(15); cout << o2scl_err;
      cout.width(15); cout << fabs(o2scl_res-legendre_exact(n));
      cout.width(15); cout << Q.feval << endl;

      n+=5;
    }
		
    Q.tol_abs/=1.0e3;
  }
  cout << endl;
	
  test.report();
  return 0;
}
