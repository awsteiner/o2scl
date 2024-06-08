 /*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Jerry Gagelman and Andrew W. Steiner
  
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
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/inte_qag_gsl.h>

using namespace std;
using namespace o2scl;

// This function oscillates quite rapidly near x=0
double test_func_1(double x) {
  return -sin(1.0/(x+0.01))*pow(x+0.01,-2.0);
}

// A GSL-version of the above function
double gsl_test_func_1(double x, void *pa) {
  return -sin(1.0/(x+0.01))*pow(x+0.01,-2.0);
}

// The QUADPACK test function
double quadpack_func(double x, double &alpha) {
  double sigma=pow(2.0,alpha);
  return cos(sigma*sin(x));
}

// GSL version of the QUADPACK test function
double quadpack_func_gsl(double x, void* params) {
  double alpha=*(double*)params;
  double sigma=pow(2.0,alpha);
  return cos(sigma*sin(x));
}

// Exact result for the QUADPACK test function
double quadpack_exact(double alpha) {
  double sigma=pow(2.0,alpha);
  return M_PI*gsl_sf_bessel_J0(sigma);
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);
  inte_qag_gsl<funct> it1;

  double ans, exact;
  int vp=0;

  cout.setf(ios::scientific);
  cout.precision(10);

  funct tf2=test_func_1;

  // Compare with the exact result
  ans=it1.integ(tf2,0.0,1.0);
  exact=cos(100.0)-cos(1/1.01);
  t.test_rel(ans,exact,1.0e-8,"qag test");

  // Compare with the GSL result for the same integral
  {
    double result, error;
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function=&gsl_test_func_1;
    gsl_integration_qag(&F,0,1,it1.tol_abs,it1.tol_rel,1000,
			GSL_INTEG_GAUSS15,w,&result,&error); 
    gsl_integration_workspace_free(w);

    cout << "O2scl            GSL              Exact" << endl;
    cout << ans << " " << result << " " << exact << endl;
    cout << it1.get_error() << " " << error << endl;
    t.test_rel(ans,result,1.0e-12,"GSL comparison 1");
    t.test_rel(it1.get_error(),error,1.0e-5,"GSL comparison 1");

  }
  cout << endl;

  cout << "Test using different integration rules for qag:" << endl;
  exact=cos(100.0)-cos(1/1.01);
  for(int h=1;h<=6;h++) {
    it1.set_rule(h);

    // For the last rule, show what verbose output looks like
    if (h==6) it1.verbose=1;
    cout.precision(6);
    ans=it1.integ(tf2,0.0,1.0);
    cout.precision(10);
    cout << h << " " << ans << " " << it1.get_error() << " " 
	 << exact << " " << fabs(ans-exact) << endl;
  }
  
  /* Test adaptive integration \ref o2scl::inte_qag_gsl the standard
     oscillatory integrand from the QUADPACK test-battery [\ref
     Piessens83],
     \f[
     \int_0^\pi cos(2^\alpha sin(x))~dx=\pi J_0(2^\alpha),
     \f]
     where \f$ J_0 \f$ is the Bessel function, for various tolerance
     values and all Gauss-Kronrod rules.
     
     The test manager compares the result of each integration against
     that obtained by the GSL routine \c gsl_integration_qag()
     with success only if the values are equal, i.e., their
     difference underflows.
  */
  cout.precision(6);
	
  double alpha=8.0;
  size_t limit=512;

  inte_qag_gsl<funct> Q;
  funct f=std::bind(quadpack_func,std::placeholders::_1,alpha);
	
  // setup GSL data
  gsl_integration_workspace *work=gsl_integration_workspace_alloc(limit);
  gsl_function f_gsl={&quadpack_func_gsl, &alpha};
	
  double o2scl_res, o2scl_err, gsl_res, gsl_err;
  
  // absolute error bound
  Q.tol_abs=1.0e-8;		
  // relative error bound
  Q.tol_rel=0.0;		

  while (Q.tol_abs > 1.5e-14) {
    cout << "\nabsolute tolerance: " << Q.tol_abs << "...\n";
    cout.width(15); cout << "GK-rule";
    cout.width(15); cout << "err. estimate";
    cout.width(15); cout << "exact error";
    cout.width(15); cout << "iterations";
    cout.width(15); cout << "|O2scl - GSL|";
    cout << endl;
		
    for (int key=1;key<=6;++key) {
      gsl_integration_qag(&f_gsl, 0, M_PI, Q.tol_abs, Q.tol_rel, 
			  limit, key, work, &gsl_res, &gsl_err);
      Q.set_rule(key);
      Q.integ_err(f, 0.0, M_PI, o2scl_res, o2scl_err);

      double dbl_eps=std::numeric_limits<double>::epsilon()*1.01;
      t.test_abs(o2scl_res,gsl_res,dbl_eps,"QAG: O2scl vs GSL");
			
      cout.width(15); cout << Q.get_rule();
      cout.width(15); cout << o2scl_err;
      cout.width(15); cout << fabs(o2scl_res-quadpack_exact(alpha));
      cout.width(15); cout << Q.last_iter;
      cout.width(15); cout << fabs(o2scl_res-gsl_res);
      cout << endl;
    }
		
    Q.tol_abs /= 10;
  }
	
  gsl_integration_workspace_free(work);

  t.report();
  return 0;
}

