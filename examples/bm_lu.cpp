/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2014, Andrew W. Steiner

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
/*
  Typical output:
  
*/
#include <iostream>

#include <gsl/gsl_linalg.h>

#include <o2scl/lu.h>
#include <o2scl/test_mgr.h>
#include <o2scl/ovector_tlate.h>
#include <o2scl/omatrix_tlate.h>
#include <o2scl/columnify.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_linalg;

template<size_t MSIZE> int run() {

  omatrix gm1(MSIZE,MSIZE), om1(MSIZE,MSIZE);
  umatrix um1(MSIZE,MSIZE);
  permutation gp1(MSIZE), op1(MSIZE), mp1(MSIZE);
  double mm1[MSIZE][MSIZE], d;
  int sig;
  size_t t1, N=(size_t)(100000/MSIZE/MSIZE/sqrt(MSIZE)*1000);
  for(size_t i=0;i<MSIZE;i++) {
    for(size_t j=0;j<MSIZE;j++) {
      gm1[i][j]=1.0+sin(i)+tan(j);
      om1[i][j]=1.0+sin(i)+tan(j);
      mm1[i][j]=1.0+sin(i)+tan(j);
    }
  }

  bool show_res=false;

  cout << "size N GSL O2scl(omatrix) O2scl(umatrix) O2scl(2d-array)" << endl;
  cout << MSIZE << " ";
  cout.width(7);
  cout << N << " ";

  for(size_t ij=0;ij<2;ij++) {
    d=0.0;
    t1=clock();
    for(size_t k=0;k<N;k++) {
      for(size_t i=0;i<MSIZE;i++) {
	for(size_t j=0;j<MSIZE;j++) {
	  gm1[i][j]=1.0+sin(i)+tan(j);
	}
      }
      gsl_linalg_LU_decomp(&gm1,&gp1,&sig);
      d+=gm1[0][0];
    }
    t1=(clock()-t1)/10000;
    cout.width(4);
    cout << t1 << " ";
    if (show_res) {
      cout << gm1 << endl;
      cout << gp1 << endl;
    }

    d=0.0;
    t1=clock();
    for(size_t k=0;k<N;k++) {
      for(size_t i=0;i<MSIZE;i++) {
	for(size_t j=0;j<MSIZE;j++) {
	  om1[i][j]=1.0+sin(i)+tan(j);
	}
      }
      LU_decomp(MSIZE,om1,op1,sig);
      d+=om1[0][0];
    }
    t1=(clock()-t1)/10000;
    cout.width(4);
    cout << t1 << " ";
    if (show_res) {
      cout << om1 << endl;
      cout << op1 << endl;
    }

    d=0.0;
    t1=clock();
    for(size_t k=0;k<N;k++) {
      for(size_t i=0;i<MSIZE;i++) {
	for(size_t j=0;j<MSIZE;j++) {
	  um1[i][j]=1.0+sin(i)+tan(j);
	}
      }
      LU_decomp(MSIZE,um1,op1,sig);
      d+=um1[0][0];
    }
    t1=(clock()-t1)/10000;
    cout.width(4);
    cout << t1 << " ";
    if (show_res) {
      cout << um1 << endl;
      cout << op1 << endl;
    }

    d=0.0;
    t1=clock();
    for(size_t k=0;k<N;k++) {
      for(size_t i=0;i<MSIZE;i++) {
	for(size_t j=0;j<MSIZE;j++) {
	  mm1[i][j]=1.0+sin(i)+tan(j);
	}
      }
      LU_decomp(MSIZE,mm1,mp1,sig);
      d+=mm1[0][0];
    }
    t1=(clock()-t1)/10000;
    cout.width(4);
    cout << t1 << " ";
    if (show_res) {
      for(size_t i=0;i<MSIZE;i++) {
	for(size_t j=0;j<MSIZE;j++) {
	  if (mm1[i][j]>=0.0) cout << " " << mm1[i][j] << " ";
	  else cout << mm1[i][j] << " ";
	}
	cout << endl;
      }
      cout << op1 << " ";
    }
  }
  cout << endl;
  return 0;
}

int main(int argv, char *argc[]) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);
  
  if (true) {
    cout << endl;
    {
      // Show that GSL range checking is off
      gsl_vector *gv=gsl_vector_alloc(2);
      double x=gsl_vector_get(gv,3);
      cout << "Gsl   range checking: 0" << endl;
    }
    cout << "O2scl range checking: " << lib_settings.range_check() << endl;
    cout << endl;
  }
  
  run<10>();
  run<12>();
  run<15>();
  run<18>();
  run<20>();
  run<25>();
  run<60>();

  //  t.report();
  return 0;
}

