/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2016, Andrew W. Steiner

  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  -------------------------------------------------------------------
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/deriv_cern.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/rng_gsl.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

double testfun(double x, double &pa) {
  return exp(pa*x)*sin(x)-cos(x)*exp(-pa*x);
}

double exact(double x, double &pa) {
  return 2.0*pa*cosh(pa*x)*(sin(x)+cos(x));
}

int tti(double x) {
  return ((int)(-log10(fabs(x))));
}


int main(void) {
  test_mgr t;
  t.set_output_level(2);
  double pa, x, gres, cres;
  size_t tmp, t1;
  rng_gsl gr;
  
  funct11 tf=std::bind(testfun,std::placeholders::_1,pa);
  deriv_cern<> cd;
  deriv_gsl<> gd;
  gd.h=1.0e-4;
  
  cout.setf(ios::scientific);

  const int max_index=15;
  double time[max_index];
  size_t n[max_index];
  double time2[max_index];
  size_t n2[max_index];
  for(int i=0;i<max_index;i++) {
    n[i]=0;
    n2[i]=0;
  }
  double acc;
  const size_t NX=100;
  const size_t N=10000;
  
  for(size_t j=0;j<20;j++) {
    
    gd.h=pow(10.0,-2-((double)j)/5.0);

    acc=0.0;
    tmp=clock();
    
    for(size_t ix=0;ix<NX;ix++) {
      pa=gr.random();
      x=gr.random();
      for(size_t k=0;k<N;k++) {
	gres=gd.deriv(x,tf);
      }
      double ex=exact(x,pa);
      acc+=fabs((gres-ex)/ex);
    }
    acc/=((double)NX);
    t1=(clock()-tmp)/((double)N);
    int ix=tti(acc);
    cout << gd.h << " " << acc << " " << ix << endl;
    if (ix>=0 && ix<max_index) {
      n[ix]++;
      time[ix]+=t1;
    }

  }
  cout << endl;

 cout << "Accuracy   deriv_gsl   deriv" << endl;
  cout << "--------------------------------------------------" << endl;
  for(int j=0;j<max_index;j++) {
    if (n[j]>0) {
      cout.width(2);
      cout << j << "         ";
      cout << time[j]/n[j] << "        " << endl;
      //<< time2[j]/n2[j] << " "
      //<< n[j] << " " << n2[j] << endl;
    }

  }


  t.report();
  return 0;
}

