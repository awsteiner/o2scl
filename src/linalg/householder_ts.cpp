/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2025, Andrew W. Steiner

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
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_linalg.h>

#include <o2scl/householder.h>
#include <o2scl/test_mgr.h>
#include <o2scl/columnify.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  {
    using namespace o2scl_cblas;
    using namespace o2scl_linalg;

    gsl_vector *gv1=gsl_vector_alloc(5);
    gsl_vector *gv2=gsl_vector_alloc(5);
    gsl_matrix *gm1=gsl_matrix_alloc(5,5);
    gsl_matrix *gm2=gsl_matrix_alloc(5,5);
    ubmatrix om1(5,5), om2(5,5);
    ubvector ov1(5), ov2(5), gv3(5), ov3(5);
    permutation gp1(5), op1(5), mp1(5);
    int sig, ret;

    // Test householder_transform and householder_hm

    for(size_t i=0;i<5;i++) {
      gsl_vector_set(gv1,i,cos(((double)(i))));
      ov1[i]=cos(((double)(i)));
      for(size_t j=0;j<5;j++) {
	gsl_matrix_set(gm1,i,j,1.0+sin(((double)(i)))+tan(((double)(j))));
	om1(i,j)=1.0+sin(((double)(i)))+tan(((double)(j)));
	gsl_matrix_set(gm2,i,j,1.0+sin(((double)(i)))+tan(((double)(j))));
	om2(i,j)=1.0+sin(((double)(i)))+tan(((double)(j)));
      }
    }
  
    gsl_linalg_householder_transform(gv1);
    householder_transform<ubvector>(5,ov1);
    t.test_rel_vec(5,ov1,gsl_vector_wrap(gv1),1.0e-12,"h trans 1");

    gsl_linalg_householder_hm(1.5,gv1,gm1);
    householder_hm(5,5,1.5,ov1,om1);
    t.test_rel_vec(5,ov1,gsl_vector_wrap(gv1),1.0e-12,"h hm 1");

    // Test householder_hv
  
    for(size_t i=0;i<5;i++) {
      gsl_vector_set(gv1,i,cos(((double)(i))));
      ov1[i]=cos(((double)(i)));
      gsl_vector_set(gv2,i,cos(((double)(i))));
      ov2[i]=cos(((double)(i)));
    }

    gsl_linalg_householder_hv(1.5,gv2,gv1);
    householder_hv<ubvector>(5,1.5,ov2,ov1);
    t.test_rel_vec(5,ov1,gsl_vector_wrap(gv1),1.0e-12,"h hv 1");
    t.test_rel_vec(5,ov2,gsl_vector_wrap(gv2),1.0e-12,"h hv 2");

    // Test householder_mh

    for(size_t i=0;i<5;i++) {
      gsl_vector_set(gv1,i,cos(((double)(i))));
      ov1[i]=cos(((double)(i)));
      for(size_t j=0;j<5;j++) {
	gsl_matrix_set(gm1,i,j,1.0+sin(((double)(i)))+tan(((double)(j))));
	om1(i,j)=1.0+sin(((double)(i)))+tan(((double)(j)));
	gsl_matrix_set(gm2,i,j,1.0+sin(((double)(i)))+tan(((double)(j))));
	om2(i,j)=1.0+sin(((double)(i)))+tan(((double)(j)));
      }
    }
  
    gsl_linalg_householder_mh(1.5,gv1,gm1);
    householder_mh(5,5,1.5,ov1,om1);
    t.test_rel_vec(5,ov1,gsl_vector_wrap(gv1),1.0e-12,"h mh 1");
    t.test_rel_mat(5,5,om1,gsl_matrix_wrap(gm1),1.0e-12,"h mh 2");

    // Test householder_hm1

    for(size_t i=0;i<5;i++) {
      gsl_vector_set(gv1,i,cos(((double)(i))));
      ov1[i]=cos(((double)(i)));
      for(size_t j=0;j<5;j++) {
	gsl_matrix_set(gm1,i,j,1.0+sin(((double)(i)))+tan(((double)(j))));
	om1(i,j)=1.0+sin(((double)(i)))+tan(((double)(j)));
	gsl_matrix_set(gm2,i,j,1.0+sin(((double)(i)))+tan(((double)(j))));
	om2(i,j)=1.0+sin(((double)(i)))+tan(((double)(j)));
      }
    }
  
    gsl_linalg_householder_hm1(sqrt(2.0),gm1);
    householder_hm1(5,5,sqrt(2.0),om1);
    t.test_rel_mat(5,5,om1,gsl_matrix_wrap(gm1),1.0e-12,"h hm1 2");

  }

  t.report();
  return 0;
}

