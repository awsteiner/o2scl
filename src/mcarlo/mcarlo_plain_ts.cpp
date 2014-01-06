/*
   -------------------------------------------------------------------
  
   Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#include <cstdlib>

/// For M_PI
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>

#include <o2scl/test_mgr.h>
#include <o2scl/multi_funct.h>
#include <o2scl/inte_multi_comp.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/mcarlo_plain.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double test_fun(size_t nv, const ubvector &x) {
  double y=1.0/(1.0-cos(x[0])*cos(x[1])*cos(x[2]))/M_PI/M_PI/M_PI;
  return y;
}

double g(double *k, size_t dim, void *params) {
  return 1.0/(1.0-cos(k[0])*cos(k[1])*cos(k[2]))/M_PI/M_PI/M_PI;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);
 
  cout.setf(ios::scientific);

  double exact = 1.3932039296;
  double res1, res2, res3;

  // GSL version
  {
    double err;

    gsl_monte_function G = {&g,3,0};

    double xl[3]={0,0,0};
    double xu[3]={M_PI,M_PI,M_PI};
    
    const gsl_rng_type *T;
    gsl_rng *r;
    
    gsl_rng_env_setup();
    
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
    
    gsl_monte_plain_state *s=gsl_monte_plain_alloc(3);

    gsl_monte_plain_integrate(&G,xl,xu,3,100000,r,s,&res1,&err);
    cout << "GSL   res,exact,err,rel: " 
	 << res1 << " " << exact << " " << err << " " 
	 << fabs(res1-exact)/err << endl;
    t.test_rel(res1,exact,err*10.0,"GSL");
    
    gsl_monte_plain_free(s);
  }

  // O2SCL version with GSL RNG
  {
    double err;

    mcarlo_plain<multi_funct<>,ubvector,int,rng_gsl> gm;
    ubvector a(3), b(3);
    a[0]=0.0;
    a[1]=0.0;
    a[2]=0.0;
    b[0]=M_PI;
    b[1]=M_PI;
    b[2]=M_PI;

    multi_funct_fptr<ubvector> tf(test_fun);

    //    gm.verbose=0;

    gm.n_points=100000;
    gm.minteg_err(tf,3,a,b,res2,err);

    cout << "O2scl res,exact,err,rel: " 
	 << res2 << " " << exact << " " << err << " " 
	 << fabs(res2-exact)/err << endl;
    t.test_rel(res2,exact,err*10.0,"O2SCL");
    cout << "O2scl-GSL: " << res2-res1 << endl;
    t.test_rel(res1,res2,1.0e-9,"O2SCL vs. GSL");
  }

#ifdef O2SCL_CPP11

  // O2SCL version with CPP11 RNG
  {
    double err;

    mcarlo_plain<multi_funct<>,ubvector,std::mt19937,
		 std::uniform_real_distribution<double> > gm;

#ifdef O2SCL_NEVER_DEFINED
}{
#endif

    ubvector a(3), b(3);
    a[0]=0.0;
    a[1]=0.0;
    a[2]=0.0;
    b[0]=M_PI;
    b[1]=M_PI;
    b[2]=M_PI;

    multi_funct_fptr<ubvector> tf(test_fun);

    //    gm.verbose=0;

    gm.n_points=100000;
    gm.minteg_err(tf,3,a,b,res2,err);

    cout << "O2scl res,exact,err,rel: " 
	 << res2 << " " << exact << " " << err << " " 
	 << fabs(res2-exact)/err << endl;
    t.test_rel(res2,exact,err*10.0,"O2SCL2");
  }

#endif

t.report();
 
return 0;
}


