/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2018, Andrew W. Steiner

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

#include <gsl/gsl_linalg.h>

#include <o2scl/lu.h>
#include <o2scl/test_mgr.h>
#include <o2scl/columnify.h>

#define BOOST_UBLAS_NDEBUG
#define BOOST_UBLAS_INLINE
#define BOOST_UBLAS_USE_FAST_SAME 

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;
typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;

using namespace std;
using namespace o2scl;
using namespace o2scl_linalg;

template<size_t MSIZE> int run() {

  ubmatrix gm1(MSIZE,MSIZE), om1(MSIZE,MSIZE);
  ubmatrix ub1, ub2;
  pmatrix pm(MSIZE);
  ub1.resize(MSIZE,MSIZE);
  ub2.resize(MSIZE,MSIZE);
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
  
  cout << MSIZE << " ";
  cout.width(7);
  cout << N << " ";

  for(size_t ij=0;ij<2;ij++) {

    // GSL 
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

    // output
    cout.width(4);
    cout << t1 << " ";
    if (show_res) {
      cout << gm1 << endl;
      cout << gp1 << endl;
    }

    // O2scl
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

    // output
    cout.width(4);
    cout << t1 << " ";
    if (show_res) {
      cout << om1 << endl;
      cout << op1 << endl;
    }

    // O2scl (arrays)
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

    // output
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

    // ublas 1
    d=0.0;
    t1=clock();
    for(size_t k=0;k<N/10;k++) {
      for(size_t i=0;i<MSIZE;i++) {
	for(size_t j=0;j<MSIZE;j++) {
	  ub1(i,j)=1.0+sin(i)+tan(j);
	}
      }
      o2scl_linalg_paren::LU_decomp(MSIZE,ub1,op1,sig);
      d+=ub1(0,0);
    }
    t1=(clock()-t1)/10000;

    // output
    cout.width(4);
    cout << t1 << " ";

    // ublas 2
    d=0.0;
    t1=clock();
    for(size_t k=0;k<N/10;k++) {
      for(size_t i=0;i<MSIZE;i++) {
	for(size_t j=0;j<MSIZE;j++) {
	  ub1(i,j)=1.0+sin(i)+tan(j);
	}
      }
      int res = lu_factorize(ub1,pm);
      d+=ub1(0,0);
    }
    t1=(clock()-t1)/10000;

    // output
    cout.width(4);
    cout << t1 << " ";

    if (ij==0) cout << " | ";

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
    cout << "GSL range checking: " 
	 << gsl_check_range << endl;	  
    cout << "O2scl range checking: " 
	 << lib_settings.range_check() << endl;
    cout << endl;
  }
  
  cout << "sz  N       GSL  O2ve O2ar ub1  ub2   | GSL  O2ve O2ar ub1  ub2" 
       << endl;
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

