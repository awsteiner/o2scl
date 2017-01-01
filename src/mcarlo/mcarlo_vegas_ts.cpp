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
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#include <o2scl/test_mgr.h>
#include <o2scl/multi_funct.h>
#include <o2scl/inte_qng_gsl.h>
#include <o2scl/mcarlo_vegas.h>

/// For M_PI
#include <gsl/gsl_math.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double test_fun(size_t nv, const ubvector &x) {
  // We make a small shift by 0.1 to avoid integrands which
  // are small everywhere
  double y=1.0/(1.0-cos(x[0])*cos(x[1])*cos(x[2]))/M_PI/M_PI/M_PI+0.1;
  return y;
}

double test_fun_gsl(double *k, size_t dim, void *params) {
  // We make a small shift by 0.1 to avoid integrands which
  // are small everywhere
  return 1.0/(1.0-cos(k[0])*cos(k[1])*cos(k[2]))/M_PI/M_PI/M_PI+0.1;
}

/*
  typedef struct {
  double alpha;
  size_t iterations;
  int stage;
  int mode;
  int verbose;
  FILE * ostream;
  } gsl_monte_vegas_params;
*/

int main(void) {
  test_mgr t;
  t.set_output_level(1);
 
  cout.setf(ios::scientific);

  double exact = 1.3932039296;
  double res1, res2, res3;

  for(int k=-1;k<=1;k+=2) {

    cout << "Mode: " << k << endl;

    // GSL version
    {
      double err;

      gsl_monte_function G = {&test_fun_gsl,3,0};

      double xl[3]={0,0,0};
      double xu[3]={M_PI,M_PI,M_PI};
    
      const gsl_rng_type *T;
      gsl_rng *r;
    
      gsl_rng_env_setup();
    
      T=gsl_rng_default;
      r=gsl_rng_alloc(T);
    
      gsl_monte_vegas_state *s=gsl_monte_vegas_alloc(3);

      gsl_monte_vegas_params params;
      gsl_monte_vegas_params_get(s,&params);
      params.mode=k;

      s->verbose=-1;
    
      gsl_monte_vegas_integrate(&G,xl,xu,3,100000,r,s,&res1,&err);
      res1-=0.1*pow(M_PI,3.0);
      cout << "GSL   res,exact,err,rel: " 
	   << res1 << " " << exact << " " << err << " " 
	   << fabs(res1-exact)/err << endl;
      t.test_rel(res1,exact,err*10.0,"GSL");
    
      gsl_monte_vegas_free (s);
    }

    // O2SCL version
    {
      double err;
      
      mcarlo_vegas<> gm;
      ubvector a(3), b(3);
      a[0]=0.0;
      a[1]=0.0;
      a[2]=0.0;
      b[0]=M_PI;
      b[1]=M_PI;
      b[2]=M_PI;

      gm.mode=k;

      multi_funct11 tf=test_fun;

      //gm.verbose=2;

      gm.n_points=100000;
      gm.minteg_err(tf,3,a,b,res2,err);
      res2-=0.1*pow(M_PI,3.0);

      cout << "O2scl res,exact,err,rel: " 
	   << res2 << " " << exact << " " << err << " " 
	   << fabs(res2-exact)/err << endl;
      t.test_rel(res2,exact,err*10.0,"O2SCL");
      cout << "O2scl-GSL: " << res2-res1 << endl;
      t.test_rel(res1,res2,1.0e-9,"O2SCL vs. GSL");
    }

    cout << endl;
  }

  for(int v=1;v<=3;v++) {

    cout << "Testing verbose output: " << v << endl;

    double res, err;
    
    ubvector a(3), b(3);
      a[0]=0.0;
      a[1]=0.0;
      a[2]=0.0;
      b[0]=M_PI;
      b[1]=M_PI;
      b[2]=M_PI;
    
    multi_funct11 tf=test_fun;
    
    mcarlo_vegas<> gm;
    gm.verbose=v;
    gm.n_points=100000;
    gm.minteg_err(tf,3,a,b,res,err);
  }
  
  t.report();
  
  return 0;
}


