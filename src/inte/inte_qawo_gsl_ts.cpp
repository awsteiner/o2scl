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
#include <gsl/gsl_integration.h>

#include <o2scl/funct.h>
#include <o2scl/inte_qawo_gsl.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double defn_odd_gsl(double theta, void* params) {
  return theta;
}

double defn_odd(double theta) {
  return theta;
}

double exact_odd(int n) {
  // sign = (-1)^{n+1}
  double sign = ((n+1) % 2) ? -1 : 1;
  return 2.0 * sign / n;
}

double defn_even_gsl(double theta, void* params) {
  return fabs(theta);
}

double defn_even(double theta) {
  return fabs(theta);
}

double exact_even(int n) {
  if (n == 0) {
    return 0.5 * M_PI;
  } else if (n % 2 == 0) {
    return 0;
  } 
  return -4.0/(M_PI*n*n);
}

int main(void) {

  cout.setf(ios::scientific);
  cout.precision(10);

  test_mgr t;
  t.set_output_level(2);

  cout << "inte_qawo_gsl:\n" << endl;

  double a=3.0, calc, exact, diff;
  inte_qawo_gsl_sin<funct> cs;
  inte_qawo_gsl_cos<funct> cc;
  
  funct_fptr f_odd(defn_odd);
  funct_fptr f_even(defn_even);

  calc=cs.integ(f_odd,0.0,1.0);
  t.test_rel(calc,0.301169,1.0e-5,"inte_qawo_gsl 2");
  
  cc.verbose=1;
  calc=cc.integ(f_odd,0.0,1.0);
  t.test_rel(calc,0.381773,1.0e-5,"inte_qawo_gsl 1");
  cc.verbose=0;

  /*
    Test adaptive integration \ref o2scl::inte_qawo_gsl against computation
    of the Fourier coefficients for non-smooth periodic functions:
    \f[
    \frac{1}{\pi} \int_{-\pi}^\pi \theta \sin(n\theta)~d\theta
    = \frac{2(-1)^{n+1}}{n}, \qquad n = 1, 2, \ldots
    \f]
    and
    \f[
    \frac{1}{\pi} \int_{-\pi}^\pi |\theta| \cos(n\theta)~d\theta
    = -\frac{4}{\pi n^2}, \quad n=1, 3, 5, \ldots.
    \f]
    for various tolerance values and levels \f$ n. \f$
    
    The test manager compares the result of each integration against that 
    obtained by the GSL routine \c gsl_integration_qawo() with success
    only if the values differ by \f$ \epsilon_\mathrm{mach}. \f$ 
  */
	
  cout.precision(6);
  t.set_output_level(1);

  size_t limit=512, levels=10;
  
  inte_qawo_gsl_cos<funct> Qcos;
  inte_qawo_gsl_sin<funct> Qsin;
	
  // setup GSL data
  gsl_integration_workspace *work=gsl_integration_workspace_alloc(limit);
  gsl_integration_qawo_table *table=
    gsl_integration_qawo_table_alloc(1,2*M_PI,GSL_INTEG_SINE,levels);
  gsl_function f_gsl;
	
  double o2scl_res,o2scl_err,gsl_res,gsl_err;

  // absolute error bound
  Qcos.tol_abs=1e-8;		
  Qsin.tol_abs=1e-8;		
  // relative error bound
  Qcos.tol_rel=0;		
  Qsin.tol_rel=0;		

  while (Qsin.tol_abs > 1.5e-13) {

    cout << "\nabsolute tolerance: " << Qsin.tol_abs << "...\n";

    cout.width(8); cout << "F-coef";
    cout.width(15); cout << "err. estimate";
    cout.width(15); cout << "exact error";
    cout.width(15); cout << "iterations";
    // cout.width(15); cout << "|O2scl - GSL|";
    cout << endl;
		
    f_gsl.function=&defn_odd_gsl;
		
    for (int n=1;n<=1000;n*=10) {
	
      gsl_integration_qawo_table_set(table,1.0*n,2*M_PI,GSL_INTEG_SINE);
      gsl_integration_qawo(&f_gsl,-M_PI,Qsin.tol_abs,Qsin.tol_rel,
			   limit,work,table,&gsl_res,&gsl_err);
      
      Qsin.omega=n;
      Qsin.integ_err(f_odd,-M_PI,M_PI,o2scl_res,o2scl_err);

      t.test_abs(gsl_res,o2scl_res,GSL_DBL_MIN,"QAWO: GSL vs O2scl");
      
      string coef="b_"+itos(n);
      cout.width(8); cout << coef;
      cout.width(15); cout << o2scl_err;
      cout.width(15); cout << fabs(o2scl_res - M_PI * exact_odd(n));
      cout.width(15); cout << Qsin.last_iter;
      cout.width(15); cout << fabs(o2scl_res - gsl_res);
      cout << endl;
    }

    f_gsl.function=&defn_even_gsl;
		
    for (int n=1;n<=1111;) {
			
      gsl_integration_qawo_table_set(table,1.0*n,2*M_PI,GSL_INTEG_COSINE);
      gsl_integration_qawo(&f_gsl,-M_PI,Qcos.tol_abs,Qcos.tol_rel,
			   limit,work,table,&gsl_res,&gsl_err);
			
      Qcos.omega=n;
      Qcos.integ_err(f_even,-M_PI,M_PI,o2scl_res,o2scl_err);
			
#ifdef O2SCL_CPP11
      double dbl_eps=std::numeric_limits<double>::epsilon();
#else 
      double dbl_eps=GSL_DBL_EPSILON;
#endif
      t.test_abs(gsl_res,o2scl_res,dbl_eps,"QAWO: GSL vs O2scl");
			
      string coef="a_"+itos(n);
      cout.width(8); cout << coef;
      cout.width(15); cout << o2scl_err;
      cout.width(15); cout << fabs(o2scl_res - M_PI * exact_even(n));
      cout.width(15); cout << Qcos.last_iter;
      cout.width(15); cout << fabs(o2scl_res - gsl_res);
      cout << endl;

      // Ensure n is always odd
      n*=10;
      n+=1;
    }
		
    Qsin.tol_abs /= 10.0;
    Qcos.tol_abs /= 10.0;
  }
	
  gsl_integration_workspace_free(work);
  gsl_integration_qawo_table_free(table);

  cout << endl;
  
  t.report();
  return 0;
}

