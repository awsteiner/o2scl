/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#include <o2scl/funct.h>
#include <o2scl/min_quad_golden.h>
 
using namespace std;
using namespace o2scl;
 
double minfun(double x);
double minfun_gsl(double x, void *p);
double minfun2(double x);
double minfun2_gsl(double x, void *p);

double minfun_gsl(double x, void *p) {
  return -exp(-(x-0.5)*(x-0.5));
}

double minfun(double x) {
  return -exp(-(x-0.5)*(x-0.5));
}

// A more pathological function with a hidden sharp minimum
double minfun2_gsl(double x, void *p) {
  return pow(fabs(x),0.01)-1.0;
}

double minfun2(double x) {
  return pow(fabs(x),0.01)-1.0;
}

int main(void) {

  std::vector<double> c1, c2, c3;
  std::vector<double> c1gsl, c2gsl, c3gsl;

  cout.setf(ios::scientific);
  cout.setf(ios::showpos);
  
  test_mgr t;
  t.set_output_level(2);

  double x, min;
  funct mf=minfun;
  funct mf2=minfun2;
  min_quad_golden<funct> mb;
  
  for(size_t k=0;k<2;k++) {

    // ---------------------------------------------------
    // GSL version

    {
      int status;
      int iter=0, max_iter=100;
    
      double a,b;
      gsl_function F;
     
      if (k==0) F.function=&minfun_gsl;
      else F.function=&minfun2_gsl;
      F.params=0;
     
      const gsl_min_fminimizer_type *T=gsl_min_fminimizer_quad_golden;
      gsl_min_fminimizer *s=gsl_min_fminimizer_alloc(T);
      gsl_min_fminimizer_set(s,&F,0.2,-1.0,1.0);
     
      do {
	iter++;
	status=gsl_min_fminimizer_iterate(s);
     
	x=gsl_min_fminimizer_x_minimum(s);
	a=gsl_min_fminimizer_x_lower(s);
	b=gsl_min_fminimizer_x_upper(s);
     
	status=gsl_min_test_interval(a,b,0.001,0.0);
	c1gsl.push_back(a);     
	c2gsl.push_back(x);     
	c3gsl.push_back(b);     

      } while (status==gsl_continue && iter<max_iter);
    
      gsl_min_fminimizer_free(s);
    }

    {
    
      // ---------------------------------------------------
      // O2scl version with set() and iterate() to compare
      // with GSL
    
      x=0.2;
      if (k==0) mb.set(mf,x,-1.0,1.0);
      else mb.set(mf2,x,-1.0,1.0);
      bool done=false;
      int iter=0;
      while (!done) {
	mb.iterate();
	iter++;
	if (gsl_min_test_interval(mb.x_lower,mb.x_upper,
				  0.001,0.0)==success) {
	  done=true;
	}
	c1.push_back(mb.x_lower);     
	c2.push_back(mb.x_minimum);     
	c3.push_back(mb.x_upper);     
      }
      if (k==0) {
	t.test_rel(mb.x_minimum,0.5,1.0e-5,"val");
	t.test_rel(mb.f_minimum,-1.0,1.0e-5,"min");
      } else {
	t.test_rel(mb.x_minimum,0.0,1.0e-4,"val");
	t.test_rel(mb.f_minimum,-1.0,1.0,"min");
      }
    
      // ---------------------------------------------------
      // O2scl version with min_bkt()
    
      x=0.2;
      if (k==0) {
	mb.min_bkt(x,-1.0,1.0,min,mf);
	t.test_rel(x,0.5,1.0e-5,"val");
	t.test_rel(min,-1.0,1.0e-5,"min");
      } else {
	mb.min_bkt(x,-1.0,1.0,min,mf2);
	t.test_rel(x,0.0,1.0e-4,"val");
	t.test_rel(min,-1.0,1.0,"min");
      }
    }

    cout.precision(12);
    cout << c1.size() << " " << c1gsl.size() << endl;
    for(size_t i=0;i<c1.size() && i<c1gsl.size();i++) {
      cout.width(2);
      cout << i << " " << c2[i] << " " << c2gsl[i] << endl;
    }
    cout.precision(6);

  }

  t.report();
  return 0;
}

