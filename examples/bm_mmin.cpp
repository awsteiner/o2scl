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
#include <iostream>
#include <o2scl/constants.h>
#include <o2scl/multi_funct.h>
#include <o2scl/gsl_mmin_simp.h>
#include <o2scl/gsl_mmin_conp.h>
#include <o2scl/gsl_mmin_conf.h>
#include <o2scl/gsl_mmin_bfgs2.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

// A very annoying function to minimize
double spring(size_t nv, const ovector_view &x, void *&pa) {
  double theta=atan2(x[1],x[0]);
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  double z=x[2];
  while (z>pi) z-=2.0*pi;
  while (z<-pi) z+=2.0*pi;
  double tmz=theta-z;
  double rm1=r-1.0;
  double ret=exp(tmz*tmz+rm1*rm1)+fabs(x[2]/10.0);
  return ret;
}

ofstream fout;

double spring_print(size_t nv, const ovector_view &x, void *&pa) {
  double theta=atan2(x[1],x[0]);
  double r=sqrt(x[0]*x[0]+x[1]*x[1]);
  double z=x[2];
  while (z>pi) z-=2.0*pi;
  while (z<-pi) z+=2.0*pi;
  double tmz=theta-z;
  double rm1=r-1.0;
  double ret=exp(tmz*tmz+rm1*rm1)+fabs(x[2]/10.0);
  fout << x[0] << " " << x[1] << " " << x[2] << " " << ret << endl;
  return ret;
}

int main(void) {
  cout.setf(ios::scientific);
  cout.setf(ios::showpos);
  cout.precision(4);

  gsl_mmin_simp<void *,multi_funct<void *> > gs;
  gsl_mmin_conp<void *,multi_funct<void *> > cp;
  gsl_mmin_conf<void *,multi_funct<void *> > cf;
  gsl_mmin_bfgs2<void *,multi_funct<void *> > bf2;
  multi_funct_fptr_noerr<void *> mfm(spring);
  multi_funct_fptr_noerr<void *> mfm2(spring_print);

  ovector x(3);
  double fmin;
  void *vp;

  gs.ntrial*=4;
  x[0]=1.0;
  x[1]=0.0;
  x[2]=2*pi;
  gs.mmin(3,x,fmin,vp,mfm);
  cout << "gsl_mmin_simp:  " << " ";
  cout << gs.last_ntrial << " " << x << " " << fmin << endl;
  cout << err_hnd->get_str() << endl;

  cp.ntrial*=4;
  x[0]=1.0;
  x[1]=0.0;
  x[2]=2*pi;
  cp.mmin(3,x,fmin,vp,mfm);
  cout << "gsl_mmin_conp:  " << " ";
  cout << cp.last_ntrial << " " << x << " " << fmin << endl;
  cout << err_hnd->get_str() << endl;

  cf.ntrial*=4;
  x[0]=1.0;
  x[1]=0.0;
  x[2]=2*pi;
  cf.mmin(3,x,fmin,vp,mfm);
  cout << "gsl_mmin_conf:  " << " ";
  cout << cf.last_ntrial << " " << x << " " << fmin << endl;
  cout << err_hnd->get_str() << endl;

  bf2.ntrial*=4;
  x[0]=1.0;
  x[1]=0.0;
  x[2]=2*pi;
  bf2.verbose=2;
  bf2.mmin(3,x,fmin,vp,mfm);
  cout << "gsl_mmin_bfgs2: " << " ";
  cout << bf2.last_ntrial << " " << x << " " << fmin << endl;
  cout << err_hnd->get_str() << endl;

  x[0]=1.0;
  x[1]=0.0;
  x[2]=4*pi;
  gs.mmin(3,x,fmin,vp,mfm2);
  cout << "gsl_mmin_simp:  " << " ";
  cout << cp.last_ntrial << " " << x << " " << fmin << endl;
  cout << err_hnd->get_str() << endl;


  return 0;
}
