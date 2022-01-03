/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Jerry Gagelman and Andrew W. Steiner
  
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
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_integration.h>

#include <o2scl/constants.h>
#include <o2scl/test_mgr.h>
#include <o2scl/funct.h>
#include <o2scl/inte_qagi_gsl.h>

using namespace std;
using namespace o2scl;

double gaussian(double x) {
  return exp(-x*x);
}

double gsl_gaussian(double x, void *params) {
  return exp(-x*x);
}

double defn_symmetric_gsl(double x, void* params) {
  double alpha=*(double*)params;
  return pow(fabs(x),alpha)*exp(-x*x);
}

double defn_symmetric(double x,double &alpha) {
  return pow(fabs(x),alpha)*exp(-x*x);
}

double exact_symmetric(double alpha) {
  return gsl_sf_gamma(0.5*(alpha+1));
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  double ans, exact;

  cout << "inte_qagi_gsl:" << endl;
  cout << endl;

  cout.setf(ios::scientific);
  cout.precision(10);
  
  inte_qagi_gsl<funct> it;
  funct tf=gaussian;
  exact=sqrt(o2scl_const::pi);
  it.verbose=1;
  ans=it.integ(tf,0.0,0.0);
  t.test_rel(ans,exact,1.0e-8,"qagi test");

  /// Compare with GSL and make sure we get the same result
  {
    gsl_integration_workspace *w=gsl_integration_workspace_alloc(1000);
    
    double result,error;
    double alpha=0.01;
    
    gsl_function F;
    F.function=&gsl_gaussian;
    F.params=&alpha;
    
    gsl_integration_qagi(&F,0,1e-7,1000,w,&result,&error); 
    
    t.test_rel(ans,result,1.0e-14,"GSL vs. O2scl");
    
    gsl_integration_workspace_free(w);
  }

  /*
    Test adaptive integration \ref o2scl::inte_qagi_gsl against the 
    symmetric integral
    \f[
    \int_{-\infty}^{\infty} |x|^\alpha e^{-x^2}~dx 
    =\Gamma\left(\frac{\alpha+1}{2}\right),
    \f]
    for various tolerance values.
    
    The test manager compares the result of each integration against that 
    obtained by the GSL routine \c gsl_integration_qagi() with success if
    the relative difference is at most \f$ \epsilon_\mathrm{mach}. \f$
  */
	
  cout.precision(6);
  t.set_output_level(1);

  size_t limit=512;
  // optional value of alpha from command line
  double alpha=1.0;
	
  inte_qagi_gsl<funct> Q;
  funct f=std::bind(defn_symmetric,std::placeholders::_1,alpha);
  
  // setup GSL data
  gsl_integration_workspace *work=gsl_integration_workspace_alloc(limit);
  gsl_function f_gsl={&defn_symmetric_gsl,&alpha};
	
  double o2scl_res,o2scl_err,gsl_res,gsl_err;
	
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
		
    gsl_integration_qagi(&f_gsl,Q.tol_abs,Q.tol_rel,limit,work,
			 &gsl_res,&gsl_err);
    Q.integ_err(f,0.0,0.0,o2scl_res,o2scl_err);
		
      double dbl_eps=std::numeric_limits<double>::epsilon();
    t.test_rel(gsl_res,o2scl_res,dbl_eps,"QAGI: GSL vs O2scl");
		
    cout.width(15); cout << o2scl_err;
    cout.width(15); cout << fabs(o2scl_res-exact_symmetric(alpha));
    cout.width(15); cout << Q.last_iter;
    cout.width(15); cout << fabs(o2scl_res-gsl_res);
    cout << endl;
		
    Q.tol_abs /= 100.0;
  }
	
  gsl_integration_workspace_free(work);
  cout << endl;

  t.report();
  return 0;
}
