/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
  along with O2scl. If not,see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#include <gsl/gsl_integration.h>

#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/inte_qags_gsl.h>

using namespace std;
using namespace o2scl;

double log_over_sqrt(double x,double &alpha) {
  double f=log(alpha*x)/sqrt(x);
  return f;
}

double log_over_sqrt_gsl(double x,void *params) {
  double alpha=*(double *)params;
  return log(alpha*x)/sqrt(x);
}

double log_recip_gsl(double x,void* params) {
  double alpha=*(double*)params;
  return pow(x,alpha)*log(1/x);
}

double log_recip(double x,double &alpha) {
  return pow(x,alpha)*log(1/x);
}

double log_recip_exact(double alpha) {
  return pow(alpha+1,-2.0);
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);
  inte_qags_gsl<funct11> it;

  double ans,exact;

  cout.setf(ios::scientific);
  cout.precision(10);

  cout << "inte_qags_gsl:" << endl;
  cout << endl;

  double alf=1.0, herr;
  funct11 tf3=std::bind(log_over_sqrt,std::placeholders::_1,alf);

  exact=-4.0;
  it.verbose=1;
  ans=it.integ(tf3,0.0,1.0);
  it.verbose=0;
  herr=it.get_error();
  t.test_rel(ans,exact,1.0e-8,"qags test");
  
  // Compare with GSL and make sure we get the same result
  {
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    
    double result,error;
    double alpha=1.0;
    
    gsl_function F;
    F.function=&log_over_sqrt_gsl;
    F.params=&alpha;
    
    gsl_integration_qags(&F,0,1,0,1e-7,1000,w,&result,&error); 
    
    t.test_rel(ans,result,1.0e-11,"GSL vs. O2scl");
    
    gsl_integration_workspace_free(w);
  }
  
  /*
    Test adaptive integration \ref o2scl::inte_qag_gsl against the
    integrand with algebraic-logarithmic singularity from the QUADPACK
    test-battery [\ref Piessens83],
    \f[
    \int_0^1 x^\alpha \log(1/x)~dx=\frac{1}{(\alpha+1)^2},
    \f]
    for various tolerance values.
    
    The test manager compares the result of each integration against
    that obtained by the GSL routine \c gsl_integration_qags() with
    success if the relative difference is at most \f$
    \epsilon_\mathrm{mach}. \f$
  */

  t.set_output_level(1);
  cout.precision(6);

  size_t limit=512;
  double alpha=-0.5;
	
  inte_qags_gsl<funct11> Q;
  funct11 f=std::bind(log_recip,std::placeholders::_1,alpha);
	
  // setup GSL data
  gsl_integration_workspace*
    work=gsl_integration_workspace_alloc(limit);
  gsl_function f_gsl={&log_recip_gsl,&alpha};
	
  double o2scl_res,o2scl_err,gsl_res,gsl_err;
	
  // absolute error bound
  Q.tol_abs=1e-9;		
  // relative error bound
  Q.tol_rel=0.0;		
	
  while (Q.tol_abs > 1.5e-14) {
    cout << "\nabsolute tolerance: " << Q.tol_abs << "...\n";
    cout.width(15); cout << "err. estimate";
    cout.width(15); cout << "exact error";
    cout.width(15); cout << "iterations";
    cout.width(15); cout << "|O2scl-GSL|";
    cout << endl;
		
    gsl_integration_qags(&f_gsl,0,1,Q.tol_abs,Q.tol_rel,
			 limit,work,&gsl_res,&gsl_err);

    Q.integ_err(f,0,1,o2scl_res,o2scl_err);
		
    double dbl_eps=std::numeric_limits<double>::epsilon();
    t.test_abs(gsl_res,o2scl_res,dbl_eps*1.01,"QAGS: GSL vs O2scl");
		
    cout.width(15); cout << o2scl_err;
    cout.width(15); cout << fabs(o2scl_res-log_recip_exact(alpha));
    cout.width(15); cout << Q.last_iter;
    cout.width(15); cout << fabs(o2scl_res-gsl_res);
    cout << endl;
    
    Q.tol_abs/=100.0;
  }
	
  gsl_integration_workspace_free(work);
  cout << endl;
  
  t.report();
  return 0;
}
