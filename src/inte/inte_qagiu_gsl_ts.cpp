 /*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Jerry Gagelman and Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_integration.h>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/inte_qagiu_gsl.h>

using namespace std;
using namespace o2scl;

double gaussian(double x) {
  return exp(-x*x);
}

double sin_recip(double x) {
  return sin(1.0/(x+0.01))*pow(x+0.01,-2.0);
}

double f(double x, void *params) {
  double f=sin(1.0/(x+0.01))*pow(x+0.01,-2.0);
  return f;
}

double defn_gamma_gsl(double x, void* params) {
  double alpha=*(double*)params;
  return pow(x, alpha)*exp(-x);
}

double defn_gamma(double x, double &alpha) {
  return pow(x, alpha)*exp(-x);
}

double exact_gamma(double alpha) {
  return gsl_sf_gamma(alpha+1);
}

double defn_laguerre_gsl(double x, void* params) {
  double alpha=*(double*)params;
  double L_n=gsl_sf_laguerre_n(1, alpha, x);
  return L_n*L_n*pow(x, alpha)*exp(-x);
}

double defn_laguerre(double x, double &alpha) {
  double L_n=gsl_sf_laguerre_n(1, alpha, x);
  return L_n*L_n*pow(x, alpha)*exp(-x);
}

double exact_laguerre(double alpha) {
  /* Orthogonality relationship:
     \f$ h_n=\Gamma(n+\alpha+1)/n! \f$
     See http://dlmf.nist.gov/18
  */ 
  return gsl_sf_gamma(alpha+2);
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  cout << "inte_qagil_gsl:" << endl;
  cout << endl;

  funct tf=gaussian;
  funct tf2=sin_recip;
  
  inte_qagiu_gsl<funct> it;

  double ans, exact;

  cout.setf(ios::scientific);
  cout.precision(8);
  
  exact=0.5*sqrt(acos(-1.0));
  it.verbose=1;
  ans=it.integ(tf,0.0,0.0);
  t.test_rel(ans,exact,1.0e-8,"qagiu test");
  it.verbose=0;

  exact=1.0-cos(100.0/101.0);
  ans=it.integ(tf2,1.0,0.0);
  t.test_rel(ans,exact,1.0e-8,"qagiu test 2");

  exact=1.0-cos(100.0);
  ans=it.integ(tf2,0.0,0.0);
  t.test_rel(ans,exact,1.0e-8,"qagiu test 2");

  {
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    
    double result, error;
    
    gsl_function F;
    F.function=&f;
    F.params=0;
    
    gsl_integration_qagiu(&F,0,0,1e-7,1000,w,&result,&error); 

    t.test_rel(ans,result,1.0e-11,"GSL vs. O2scl");
    
    gsl_integration_workspace_free(w);
  }

  /**
     Test adaptive integration \ref o2scl::inte_qagiu_gsl against the definition
     of the Gamma function,
     \f[
     \int_0^\infty x^\alpha e^{-x}~dx=\Gamma(\alpha+1),
     \f]
     for various tolerance values.
 
     The test manager compares the result of each integration against that 
     obtained by the GSL routine \c gsl_integration_qagiu() with success if
     the relative difference is at most \f$ \epsilon_\mathrm{mach}. \f$
  */

  cout.precision(6);
  t.set_output_level(1);

  size_t limit=512;
  double alpha=1.0;
	
  inte_qagiu_gsl<funct> Q;
  funct f=std::bind(defn_laguerre,std::placeholders::_1,alpha);
	
  // setup GSL data
  gsl_integration_workspace*
    work=gsl_integration_workspace_alloc(limit);
  gsl_function f_gsl={&defn_laguerre_gsl, &alpha};
	
  double o2scl_res, o2scl_err, gsl_res, gsl_err;
	
  // absolute error bound
  Q.tol_abs=1e-8;		
  // relative error bound
  Q.tol_rel=0.0;		
	
  while (Q.tol_abs > 1.5e-13) {
    cout << "\nabsolute tolerance: " << Q.tol_abs << "...\n";
    cout.width(15); cout << "err. estimate";
    cout.width(15); cout << "exact error";
    cout.width(15); cout << "iterations";
    cout.width(15); cout << "|O2scl-GSL|";
    cout << endl;
		
    gsl_integration_qagiu(&f_gsl, 0.0, Q.tol_abs, Q.tol_rel, limit, work, 
			  &gsl_res, &gsl_err);
		
    Q.integ_err(f, 0.0, 0, o2scl_res, o2scl_err);
		
    double dbl_eps=std::numeric_limits<double>::epsilon()*1.01;
    t.test_rel(gsl_res,o2scl_res,dbl_eps,"QAGIU: O2scl vs GSL");
		
    cout.width(15); cout << o2scl_err;
    cout.width(15); cout << fabs(o2scl_res-exact_laguerre(alpha));
    cout.width(15); cout << Q.last_iter;
    cout.width(15); cout << fabs(o2scl_res-gsl_res);
    cout << endl;
		
    Q.tol_abs /= 10;
  }
	
  gsl_integration_workspace_free(work);
  cout << endl;

  t.report();
  return 0;
}
