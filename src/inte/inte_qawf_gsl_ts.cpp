/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2015, Jerry Gagelman
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

#include <cmath>
#include <gsl/gsl_integration.h>

#include <o2scl/funct.h>
#include <o2scl/inte_qawf_gsl.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double quad_exp(double x) {
  return x*x*exp(-x);
}

double exponential_gsl(double x, void* params) {
  double alpha=*(double*)params;
  return exp(-alpha*x);
}

double exponential(double x, double &alpha) {
  return exp(-alpha*x);
}

double exponential_exact(double alpha, double omega) {
  double r2=alpha*alpha+omega*omega;
  return alpha/r2;
}

int main(void) {
  double a=3.0,calc,exact,diff;
  inte_qawf_gsl_sin<funct> cs;
  inte_qawf_gsl_cos<funct> cc;
  test_mgr t;
  t.set_output_level(2);

  funct tf=quad_exp;

  cout.setf(ios::scientific);
  
  calc=cs.integ(tf,0.0,1.0);
  t.test_rel(calc,0.5,1.0e-14,"inte_qawf_gsl 2");
  
  cc.verbose=1;
  calc=cc.integ(tf,0.0,1.0);
  cc.verbose=0;
  
  t.test_rel(calc,-0.5,1.0e-14,"inte_qawf_gsl 1");

  /*
    The routine gsl_integration_qawf has trouble with the
    integral from the QUADPACK test battery [\ref Piessens83]
    \f[
    \int_0^\infty x^{-1/2}\cos(\pi x/2)~dx=1.
    \f]
    In its place,the (trivially \f$ L^1 \f$) integrand is used:
    \f[
    \int_0^\infty e^{-\alpha x}\cos(\omega x)~dx =
    \mathrm{Re} \int_0^\infty e^{(i\omega-\alpha)x}~dx =
    \frac{\alpha}{\alpha^2+\omega^2}.
    \f]
    
    The test manager compares the result of each integration against that 
    obtained by the GSL routine \c gsl_integration_qawf() with success
    only if the values differ by \f$ \epsilon_\mathrm{mach}. \f$ 
  */
  t.set_output_level(1);

  double alpha=1.0;
  double omega=M_PI;
	
  size_t limit=512, levels=10;
	
  inte_qawf_gsl_cos<funct> Q;
  Q.omega=omega;
  funct f=std::bind(exponential,std::placeholders::_1,alpha);
		
  // setup GSL data
  gsl_integration_workspace
    *ws_main=gsl_integration_workspace_alloc(limit),
    *ws_cycles=gsl_integration_workspace_alloc(limit);
  gsl_integration_qawo_table *table=gsl_integration_qawo_table_alloc
    (omega,1,GSL_INTEG_COSINE,levels);
	
  gsl_function f_gsl={&exponential_gsl,&alpha};
	
  double o2scl_res,o2scl_err,gsl_res,gsl_err;
	
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
		
    gsl_integration_qawf(&f_gsl,0,Q.tol_abs,limit,ws_main,
			 ws_cycles,table,&gsl_res,&gsl_err);
    Q.integ_err(f,0,0,o2scl_res,o2scl_err);
		
      double dbl_eps=std::numeric_limits<double>::epsilon();
    t.test_abs(gsl_res,o2scl_res,dbl_eps,"QAWF: GSL vs O2scl");
		
    cout.width(15); cout << o2scl_err;
    cout.width(15); cout << fabs(o2scl_res-exponential_exact(alpha,omega));
    cout.width(15); cout << Q.last_iter;
    cout.width(15); cout << fabs(o2scl_res-gsl_res);
    cout << endl;
		
    Q.tol_abs /= 10.0;
  }
	
  gsl_integration_workspace_free(ws_main);
  gsl_integration_workspace_free(ws_cycles);
  
  t.report();
  return 0;
}

