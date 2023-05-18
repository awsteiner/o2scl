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
#include <gsl/gsl_integration.h>

#include <o2scl/funct.h>
#include <o2scl/inte_qawc_gsl.h>
#include <o2scl/inte_cauchy_cern.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double exponential(double tx) {
  return exp(-tx);
}

double logp1_fun(double tx) {
  return log(tx+1.0);
}

double recip_x3_gsl(double x, void* params) {
  return 1.0/(5*x*x*x+6);
}

double recip_x3(double x) {
  return 1.0/(5*x*x*x+6);
}

double recip_x3_exact() {
  return log(125.0/631.0)/18.0;
}

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  
  double a=3.0, calc;
  inte_qawc_gsl<funct> cc;

  funct tf=exponential;
  funct tf2=logp1_fun;

  cout << "inte_qawc_gsl: " << endl;
  cout << endl;
  
  cc.s=2.0;
  calc=cc.integ(tf,1.0,3.0);
  t.test_rel(calc,-0.286166693342262,1.0e-14,"inte_qawc_gsl 1");
  
  cc.s=1.0;
  cc.verbose=1;
  calc=cc.integ(tf2,0.0,2.0);
  t.test_rel(calc,1.030654733388659,1.0e-14,"inte_qawc_gsl 2");
  cc.verbose=0;

  /**
     Test adaptive integration \ref o2scl::inte_qawc_gsl against the 
     integral from the QUADPACK test battery [\ref Piessens83]
     \f[
     \int_{-1}^5 \frac{dx}{x(5x^3+6)}=\frac{\log(125/631)}{18},
     \f]
     for various tolerance values. 
 
     The test manager compares the result of each integration against
     that obtained by the GSL routine \c gsl_integration_qawc() with
     success only if the values differ by \f$ \epsilon_\mathrm{mach}.
     \f$
  */
  t.set_output_level(1);

  size_t limit=512;
	
  inte_qawc_gsl<funct> Q;
  funct f=recip_x3;
	
  double aa=-1, b=5, c=0;
  Q.s=c;
	
  // setup GSL data
  gsl_integration_workspace *work=gsl_integration_workspace_alloc(limit);
  gsl_function f_gsl={&recip_x3_gsl,&aa};
	
  double o2scl_res, o2scl_err, gsl_res, gsl_err;
	
  // absolute error bound
  Q.tol_abs=1e-8;		
  // relative error bound
  Q.tol_rel=0;		
	
  while (Q.tol_abs > 1.5e-13) {
    cout << "\nabsolute tolerance: " << Q.tol_abs << "...\n";
    cout.width(15); cout << "err. estimate";
    cout.width(15); cout << "exact error";
    cout.width(15); cout << "iterations";
    cout.width(15); cout << "|O2scl-GSL|";
    cout << endl;
		
    gsl_integration_qawc(&f_gsl,aa,b,c,Q.tol_abs,Q.tol_rel,
			 limit,work,&gsl_res,&gsl_err);

    Q.integ_err(f,aa,b,o2scl_res,o2scl_err);
		
      double dbl_eps=std::numeric_limits<double>::epsilon();
    t.test_abs(gsl_res,o2scl_res,dbl_eps,"QAWC: GSL vs O2scl");
		
    cout.width(15); cout << o2scl_err;
    cout.width(15); cout << fabs(o2scl_res-recip_x3_exact());
    cout.width(15); cout << Q.last_iter;
    cout.width(15); cout << fabs(o2scl_res-gsl_res);
    cout << endl;

    Q.tol_abs/=10.0;
  }
	
  gsl_integration_workspace_free(work);
  cout << endl;
  
  t.report();
  return 0;
}

