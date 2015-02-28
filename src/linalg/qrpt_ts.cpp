/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2015, Andrew W. Steiner

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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <o2scl/qrpt.h>
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

    gsl_matrix *gm1=gsl_matrix_alloc(5,5);
    gsl_matrix *gm2=gsl_matrix_alloc(5,5);
    gsl_vector *gv1=gsl_vector_alloc(5);
    gsl_vector *gv2=gsl_vector_alloc(5);
    gsl_permutation *gp1=gsl_permutation_alloc(5);
    ubmatrix om1(5,5), om2(5,5);
    ubvector ov1(5), ov2(5);
    permutation op1(5);
    int gsig, osig, msig;

    for(size_t ik=0;ik<10;ik++) {

      // Test decomp

      for(size_t i=0;i<5;i++) {
	gsl_vector_set(gv1,i,cos(((double)(i))));
	gsl_vector_set(gv2,i,cos(((double)(i))));
	ov1[i]=cos(((double)(i)));
	ov2[i]=cos(((double)(i)));
	for(size_t j=0;j<5;j++) {
	  gsl_matrix_set(gm1,i,j,((double)(ik+1))/(1.0+i+j));
	  om1(i,j)=((double)(ik+1))/(1.0+i+j);
	  gsl_matrix_set(gm2,i,j,((double)(ik+1))/(1.0+i+j));
	  om2(i,j)=((double)(ik+1))/(1.0+i+j);
	}
      }
      
      gsl_linalg_QRPT_decomp(gm1,gv1,gp1,&gsig,gv2);
      QRPT_decomp(5,5,om1,ov1,op1,osig,ov2);
      t.test_rel_mat(5,5,om1,gsl_matrix_wrap(gm1),1.0e-11,"qrpt decomp 1");
      t.test_rel_vec(5,ov1,gsl_vector_wrap(gv1),1.0e-11,"qrpt decomp 2");
      // We decrease the accuracy to 10^{-6} here because the last
      // element is six orders of magnitude smaller than the first,
      // and thus easily susceptible to precision problems.
      t.test_rel_vec(5,ov2,gsl_vector_wrap(gv2),5.0e-6,"qrpt decomp 3");

    }

  }

  t.report();

  return 0;
}
