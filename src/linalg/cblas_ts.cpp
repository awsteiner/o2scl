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
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_blas.h>

#include <o2scl/test_mgr.h>
#include <o2scl/cblas.h>
#include <o2scl/columnify.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;
typedef boost::numeric::ublas::matrix<double> ubmatrix;

void reset_data(ubvector &v1, ubvector &v2, ubvector &v3,
		gsl_vector *g1, gsl_vector *g2, gsl_vector *g3,
		gsl_vector *h1, gsl_vector *h2, gsl_vector *h3,
		ubmatrix &m1, ubmatrix &m2, ubmatrix &m3,
		gsl_matrix *gm1, gsl_matrix *gm2, gsl_matrix *gm3,
		gsl_matrix *hm1, gsl_matrix *hm2, gsl_matrix *hm3) {
  
  for(size_t i=0;i<5;i++) {
    v1[i]=sin((double)(i+1));
    v2[i]=cos((double)(i+1));
    v3[i]=tan((double)(i+1));
    gsl_vector_set(g1,i,sin((double)(i+1)));
    gsl_vector_set(g2,i,cos((double)(i+1)));
    gsl_vector_set(g3,i,tan((double)(i+1)));
    gsl_vector_set(h1,i,sin((double)(i+1)));
    gsl_vector_set(h2,i,cos((double)(i+1)));
    gsl_vector_set(h3,i,tan((double)(i+1)));
    for(size_t j=0;j<5;j++) {
      m1(i,j)=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      m2(i,j)=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      m3(i,j)=0.0;
      gsl_matrix_set(gm1,i,j,sqrt(((double)(i+1.0)))*((double)(j+1.0)));
      gsl_matrix_set(hm1,i,j,sqrt(((double)(i+1.0)))*((double)(j+1.0)));
      gsl_matrix_set(gm2,i,j,sqrt(((double)(i+1.0)))*((double)(j+1.0)));
      gsl_matrix_set(hm2,i,j,sqrt(((double)(i+1.0)))*((double)(j+1.0)));
      gsl_matrix_set(gm3,i,j,0.0);
      gsl_matrix_set(hm3,i,j,0.0);
    }
  }
  
  return;
}

int main(void) {
  test_mgr t;
  t.set_output_level(1);

  ubvector v1(5), v2(5), v3(5);
  gsl_vector *g1, *g2, *g3;
  gsl_vector *h1, *h2, *h3;
  g1=gsl_vector_alloc(5);
  g2=gsl_vector_alloc(5);
  g3=gsl_vector_alloc(5);
  h1=gsl_vector_alloc(5);
  h2=gsl_vector_alloc(5);
  h3=gsl_vector_alloc(5);

  ubmatrix m1(5,5);
  ubmatrix m2(5,5);
  ubmatrix m3(5,5);
  gsl_matrix *gm1, *gm2, *gm3;
  gsl_matrix *hm1, *hm2, *hm3;
  gm1=gsl_matrix_alloc(5,5);
  gm2=gsl_matrix_alloc(5,5);
  gm3=gsl_matrix_alloc(5,5);
  hm1=gsl_matrix_alloc(5,5);
  hm2=gsl_matrix_alloc(5,5);
  hm3=gsl_matrix_alloc(5,5);

  double res1, res2, res3;

  // Test with operator[]
  {
    cout << "Test with operator[]: " << endl;
    
    // daxpy
    cout << "daxpy: " << endl;
    
    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);
  
    gsl_blas_daxpy(0.5,g1,g2);
    cblas_daxpy(5,0.5,h1->data,1,h2->data,1);
    o2scl_cblas::daxpy(0.5,5,v1,v2);
    t.test_rel_vec(5,gsl_vector_wrap(g2),
		   gsl_vector_wrap(h2),1.0e-12,"daxpy1 operator[]");
    t.test_rel_vec(5,v2,gsl_vector_wrap(h2),1.0e-12,"daxpy2 operator[] ");
    
    // ddot
    cout << "ddot: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_ddot(g1,g2,&res1);
    res2=cblas_ddot(5,h1->data,1,h2->data,1);
    res3=o2scl_cblas::ddot(5,v1,v2);
    t.test_rel(res1,res2,1.0e-12,"ddot1 operator[] ");
    t.test_rel(res1,res3,1.0e-12,"ddot2 operator[] ");
  
    // dnrm2
    cout << "dnrm2: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    res1=gsl_blas_dnrm2(g1);
    res2=cblas_dnrm2(5,h1->data,1);
    res3=o2scl_cblas::dnrm2(5,v1);
    t.test_rel(res1,res2,1.0e-12,"dnrm21 operator[] ");
    t.test_rel(res1,res3,1.0e-12,"dnrm22 operator[] ");
  
    // dscal
    cout << "dscal: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dscal(0.5,g3);
    cblas_dscal(5,0.5,h3->data,1);
    o2scl_cblas::dscal(0.5,5,v3);
    t.test_rel_vec(5,gsl_vector_wrap(g3),
		   gsl_vector_wrap(h3),1.0e-12,"dscal1 operator[] ");
    t.test_rel_vec(5,v3,gsl_vector_wrap(h3),1.0e-12,"dscal2 operator[] ");

    // dgemv, NoTrans
    cout << "dgemv, NoTrans: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dgemv(CblasNoTrans,0.1,gm1,g1,0.2,g2);
    cblas_dgemv(CblasRowMajor,CblasNoTrans,5,5,0.1,hm1->data,5,
		h1->data,1,0.2,h2->data,1);
    o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_NoTrans,
		       5,5,0.1,m1,v1,0.2,v2);
    t.test_rel_vec(5,gsl_vector_wrap(g2),
		   gsl_vector_wrap(h2),1.0e-12,"dgemv1 operator[] ");
    t.test_rel_vec(5,v2,gsl_vector_wrap(h2),1.0e-12,"dgemv2 operator[] ");

    // dgemv, Trans (Note that ConjTrans was not originally supported)
    cout << "dgemv, Trans: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dgemv(CblasTrans,0.1,gm1,g1,0.2,g2);
    cblas_dgemv(CblasRowMajor,CblasTrans,5,5,0.1,hm1->data,5,
		h1->data,1,0.2,h2->data,1);
    o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Trans,
		       5,5,0.1,m1,v1,0.2,v2);
    t.test_rel_vec(5,gsl_vector_wrap(g2),gsl_vector_wrap(h2),1.0e-12,
		   "dgemv1t operator[] ");
    t.test_rel_vec(5,v2,gsl_vector_wrap(h2),1.0e-12,"dgemv2t operator[] ");

    // dtrsv, Upper, NoTrans, NonUnit
    cout << "dtrsv, Upper, NoTrans, NonUnit: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dtrsv(CblasUpper,CblasNoTrans,CblasNonUnit,gm1,g3);
    cblas_dtrsv(CblasRowMajor,CblasUpper,CblasNoTrans,CblasNonUnit,
		5,hm1->data,5,h3->data,1);
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Upper,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,5,5,m1,v3);
    t.test_rel_vec(5,gsl_vector_wrap(g3),gsl_vector_wrap(h3),1.0e-12,
		   "dtrsv1unn operator[] ");
    t.test_rel_vec(5,v3,gsl_vector_wrap(h3),1.0e-12,"dtrsv2unn operator[] ");

    // dtrsv, Lower, NoTrans, NonUnit
    cout << "dtrsv, Lower, NoTrans, NonUnit: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dtrsv(CblasLower,CblasNoTrans,CblasNonUnit,gm1,g3);
    cblas_dtrsv(CblasRowMajor,CblasLower,CblasNoTrans,CblasNonUnit,
		5,hm1->data,5,h3->data,1);
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Lower,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,5,5,m1,v3);
    t.test_rel_vec(5,gsl_vector_wrap(g3),gsl_vector_wrap(h3),1.0e-12,
		   "dtrsv1lnn operator[] ");
    t.test_rel_vec(5,v3,gsl_vector_wrap(h3),1.0e-12,"dtrsv2lnn operator[] ");

    // dtrsv, Upper, Trans, NonUnit
    cout << "dtrsv, Upper, Trans, NonUnit: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dtrsv(CblasUpper,CblasTrans,CblasNonUnit,gm1,g3);
    cblas_dtrsv(CblasRowMajor,CblasUpper,CblasTrans,CblasNonUnit,
		5,hm1->data,5,h3->data,1);
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Upper,
		       o2scl_cblas::o2cblas_Trans,
		       o2scl_cblas::o2cblas_NonUnit,5,5,m1,v3);
    t.test_rel_vec(5,gsl_vector_wrap(g3),gsl_vector_wrap(h3),1.0e-12,
		   "dtrsv1utn operator[] ");
    t.test_rel_vec(5,v3,gsl_vector_wrap(h3),1.0e-12,"dtrsv2utn operator[] ");

    // dtrsv, Lower, Trans, NonUnit
    cout << "dtrsv, Lower, Trans, NonUnit: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dtrsv(CblasLower,CblasTrans,CblasNonUnit,gm1,g3);
    cblas_dtrsv(CblasRowMajor,CblasLower,CblasTrans,CblasNonUnit,
		5,hm1->data,5,h3->data,1);
    o2scl_cblas::dtrsv(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Lower,
		       o2scl_cblas::o2cblas_Trans,
		       o2scl_cblas::o2cblas_NonUnit,5,5,m1,v3);
    t.test_rel_vec(5,gsl_vector_wrap(g3),gsl_vector_wrap(h3),1.0e-12,
		   "dtrsv1ltn operator[] ");
    t.test_rel_vec(5,v3,gsl_vector_wrap(h3),1.0e-12,"dtrsv2ltn operator[] ");

    // dtrsm, Left, Upper, NoTrans, NonUnit
    cout << "dtrsm, Left, Upper, NoTrans, NonUnit: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dtrsm(CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,2.0,gm1,gm2);
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasUpper,CblasNoTrans,CblasNonUnit,
		5,5,2.0,hm1->data,5,hm2->data,5);
    o2scl_cblas::dtrsm(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Left,
		       o2scl_cblas::o2cblas_Upper,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,5,5,2.0,m1,m2);
    t.test_rel_nonzero_mat(5,5,gsl_matrix_wrap(gm2),
		   gsl_matrix_wrap(hm2),1.0e-12,1.0e-15,"dtrsm 1");
    t.test_rel_nonzero_mat(5,5,m2,gsl_matrix_wrap(hm2),
			   1.0e-12,1.0e-15,"dtrsm 2");

    // dtrsm, Right, Upper, NoTrans, NonUnit
    cout << "dtrsm, Right, Upper, NoTrans, NonUnit: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dtrsm(CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,2.0,gm1,gm2);
    cblas_dtrsm(CblasRowMajor,CblasRight,CblasUpper,CblasNoTrans,CblasNonUnit,
		5,5,2.0,hm1->data,5,hm2->data,5);
    o2scl_cblas::dtrsm(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Right,
		       o2scl_cblas::o2cblas_Upper,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,5,5,2.0,m1,m2);
    t.test_rel_nonzero_mat(5,5,gsl_matrix_wrap(gm2),
		   gsl_matrix_wrap(hm2),1.0e-12,1.0e-15,"dtrsm 3");
    t.test_rel_nonzero_mat(5,5,m2,gsl_matrix_wrap(hm2),
			   1.0e-12,1.0e-15,"dtrsm 4");

    // dtrsm, Left, Lower, NoTrans, NonUnit
    cout << "dtrsm, Left, Lower, NoTrans, NonUnit: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dtrsm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,2.0,gm1,gm2);
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,
		5,5,2.0,hm1->data,5,hm2->data,5);
    o2scl_cblas::dtrsm(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Left,
		       o2scl_cblas::o2cblas_Lower,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NonUnit,5,5,2.0,m1,m2);
    t.test_rel_nonzero_mat(5,5,gsl_matrix_wrap(gm2),
		   gsl_matrix_wrap(hm2),1.0e-12,1.0e-15,"dtrsm 5");
    t.test_rel_nonzero_mat(5,5,m2,gsl_matrix_wrap(hm2),
			   1.0e-12,1.0e-15,"dtrsm 6");

    // dtrsm, Left, Upper, Trans, NonUnit
    cout << "dtrsm, Left, Upper, Trans, NonUnit: " << endl;

    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    gsl_blas_dtrsm(CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,2.0,gm1,gm2);
    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasUpper,CblasTrans,CblasNonUnit,
		5,5,2.0,hm1->data,5,hm2->data,5);
    o2scl_cblas::dtrsm(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Left,
		       o2scl_cblas::o2cblas_Upper,
		       o2scl_cblas::o2cblas_Trans,
		       o2scl_cblas::o2cblas_NonUnit,5,5,2.0,m1,m2);
    t.test_rel_nonzero_mat(5,5,gsl_matrix_wrap(gm2),
		   gsl_matrix_wrap(hm2),1.0e-12,1.0e-15,"dtrsm 7");
    t.test_rel_nonzero_mat(5,5,m2,gsl_matrix_wrap(hm2),
			   1.0e-12,1.0e-15,"dtrsm 8");

    // dgemm, NoTrans, NoTrans
    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    cout << "dgemm, NoTrans, NoTrans: " << endl;
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,0.1,gm1,gm2,0.2,gm3);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,5,5,5,0.1,hm1->data,5,
		hm2->data,5,0.2,hm3->data,5);
    o2scl_cblas::dgemm(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NoTrans,
		       5,5,5,0.1,m1,m2,0.2,m3);
    t.test_rel_mat(5,5,gsl_matrix_wrap(gm3),
		   gsl_matrix_wrap(hm3),1.0e-12,
		   "dgemm1ntnt operator[] ");
    t.test_rel_mat(5,5,m3,gsl_matrix_wrap(hm3),1.0e-12,
		   "dgemm2ntnt operator[] ");
    
    // dgemm, Trans, NoTrans
    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    cout << "dgemm, Trans, NoTrans: " << endl;
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,0.1,gm1,gm2,0.2,gm3);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,5,5,5,0.1,hm1->data,5,
		hm2->data,5,0.2,hm3->data,5);
    o2scl_cblas::dgemm(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Trans,
		       o2scl_cblas::o2cblas_NoTrans,
		       5,5,5,0.1,m1,m2,0.2,m3);
    t.test_rel_mat(5,5,gsl_matrix_wrap(gm3),gsl_matrix_wrap(hm3),1.0e-12,
		   "dgemm1tnt operator[] ");
    t.test_rel_mat(5,5,m3,gsl_matrix_wrap(hm3),1.0e-12,"dgemm2tnt operator[] ");

    // dgemm, NoTrans, Trans
    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    cout << "dgemm, NoTrans, Trans: " << endl;
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,0.1,gm1,gm2,0.2,gm3);
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,5,5,5,0.1,hm1->data,5,
		hm2->data,5,0.2,hm3->data,5);
    o2scl_cblas::dgemm(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_Trans,
		       5,5,5,0.1,m1,m2,0.2,m3);
    t.test_rel_mat(5,5,gsl_matrix_wrap(gm3),gsl_matrix_wrap(hm3),1.0e-12,
		   "dgemm1ntt operator[] ");
    t.test_rel_mat(5,5,m3,gsl_matrix_wrap(hm3),1.0e-12,"dgemm2ntt operator[] ");

    // dgemm, Trans, Trans
    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    cout << "dgemm, Trans, Trans: " << endl;
    gsl_blas_dgemm(CblasTrans,CblasTrans,0.1,gm1,gm2,0.2,gm3);
    cblas_dgemm(CblasRowMajor,CblasTrans,CblasTrans,5,5,5,0.1,hm1->data,5,
		hm2->data,5,0.2,hm3->data,5);
    o2scl_cblas::dgemm(o2scl_cblas::o2cblas_RowMajor,
		       o2scl_cblas::o2cblas_Trans,
		       o2scl_cblas::o2cblas_Trans,
		       5,5,5,0.1,m1,m2,0.2,m3);
    t.test_rel_mat(5,5,gsl_matrix_wrap(gm3),gsl_matrix_wrap(hm3),1.0e-12,
		   "dgemm1tt operator[] ");
    t.test_rel_mat(5,5,m3,gsl_matrix_wrap(hm3),1.0e-12,"dgemm2tt operator[] ");

    // dgemm, ColMajor, NoTrans, NoTrans
    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    cout << "dgemm, ColMajor, NoTrans, NoTrans: " << endl;
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,5,5,5,0.1,hm1->data,5,
		hm2->data,5,0.2,hm3->data,5);
    o2scl_cblas::dgemm(o2scl_cblas::o2cblas_ColMajor,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_NoTrans,
		       5,5,5,0.1,m1,m2,0.2,m3);
    t.test_rel_mat(5,5,m3,gsl_matrix_wrap(hm3),1.0e-12,
		   "dgemmc2ntnt operator[] ");
    
    // dgemm, ColMajor, Trans, NoTrans
    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    cout << "dgemm, ColMajor, Trans, NoTrans: " << endl;
    cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,5,5,5,0.1,hm1->data,5,
		hm2->data,5,0.2,hm3->data,5);
    o2scl_cblas::dgemm(o2scl_cblas::o2cblas_ColMajor,
		       o2scl_cblas::o2cblas_Trans,
		       o2scl_cblas::o2cblas_NoTrans,
		       5,5,5,0.1,m1,m2,0.2,m3);
    t.test_rel_mat(5,5,m3,gsl_matrix_wrap(hm3),1.0e-12,
		   "dgemmc2tnt operator[] ");

    // dgemm, ColMajor, NoTrans, Trans
    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    cout << "dgemm, ColMajor, NoTrans, Trans: " << endl;
    cblas_dgemm(CblasColMajor,CblasNoTrans,CblasTrans,5,5,5,0.1,hm1->data,5,
		hm2->data,5,0.2,hm3->data,5);
    o2scl_cblas::dgemm(o2scl_cblas::o2cblas_ColMajor,
		       o2scl_cblas::o2cblas_NoTrans,
		       o2scl_cblas::o2cblas_Trans,
		       5,5,5,0.1,m1,m2,0.2,m3);
    t.test_rel_mat(5,5,m3,gsl_matrix_wrap(hm3),1.0e-12,
		   "dgemmc2ntt operator[] ");

    // dgemm, ColMajor, Trans, Trans
    reset_data(v1,v2,v3,g1,g2,g3,h1,h2,h3,
	       m1,m2,m3,gm1,gm2,gm3,hm1,hm2,hm3);

    cout << "dgemm, ColMajor, Trans, Trans: " << endl;
    cblas_dgemm(CblasColMajor,CblasTrans,CblasTrans,5,5,5,0.1,hm1->data,5,
		hm2->data,5,0.2,hm3->data,5);
    o2scl_cblas::dgemm(o2scl_cblas::o2cblas_ColMajor,
		       o2scl_cblas::o2cblas_Trans,
		       o2scl_cblas::o2cblas_Trans,
		       5,5,5,0.1,m1,m2,0.2,m3);
    t.test_rel_mat(5,5,m3,gsl_matrix_wrap(hm3),1.0e-12,
		   "dgemmc2tt operator[] ");

#ifdef O2SCL_NEVER_DEFINED

    // daxpy_subvec
    cout << "daxpy_subvec: " << endl;
    for(size_t i=0;i<5;i++) {
      v1[i]=sin((double)(i+1));
      v2[i]=cos((double)(i+1));
      gsl_vector_set(g1,i,sin((double)(i+1)));
      gsl_vector_set(g2,i,cos((double)(i+1)));
    }
    for(size_t k=0;k<5;k++) {
      ubvector_subvector vsub1(v1,k,5-k);
      ubvector_subvector vsub2(v2,k,5-k);
      o2scl_cblas::daxpy(0.5,5-k,vsub1,vsub2);
      o2scl_cblas::daxpy_subvec(0.5,5,g1,g2,k);
      t.test_rel_vec(5,v2,g2,1.0e-12,"daxpy_subvec1 operator[] ");
    }

    // ddot_subvec
    cout << "ddot_subvec: " << endl;
    for(size_t k=0;k<5;k++) {
      ubvector_subvector vsub1(v1,k,5-k);
      ubvector_subvector vsub2(v2,k,5-k);
      double tx1=o2scl_cblas::ddot(5-k,vsub1,vsub2);
      double tx2=o2scl_cblas::ddot_subvec(5,g1,g2,k);
      t.test_rel(tx1,tx2,1.0e-12,"ddot_subvec1 operator[] ");
    }

    // dscal_subvec
    cout << "dscal_subvec: " << endl;
    for(size_t k=0;k<5;k++) {
      ubvector_subvector vsub1(v1,k,5-k);
      o2scl_cblas::dscal(0.5,5-k,vsub1);
      o2scl_cblas::dscal_subvec(0.5,5,g1,k);
    }
  
    // dnrm2_subvec
    cout << "dnrm2_subvec: " << endl;
    for(size_t k=0;k<5;k++) {
      ubvector_subvector vsub1(v1,k,5-k);
      double tx1=o2scl_cblas::dnrm2(5-k,vsub1);
      double tx2=o2scl_cblas::dnrm2_subvec(5,g1,k);
      t.test_rel(tx1,tx2,1.0e-12,"dnrm2_subvec1 operator[] ");
    }

    // dnrm2_subcol
    cout << "dnrm2_subcol: " << endl;
    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
	m1[i][j]=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      }
    }
    for(size_t k=0;k<5;k++) {
      for(size_t ell=0;ell<5;ell++) {
	ubmatrix_col mcol1(m1,k);
	ubvector_subvector vsub1(mcol1,ell,5-ell);
	double tx1=o2scl_cblas::dnrm2(5-ell,vsub1);
	double tx2=o2scl_cblas::dnrm2_subcol(m1,ell,k,5);
	t.test_rel(tx1,tx2,1.0e-12,"dnrm2_subcol operator[] ");
      }
    }

    // dscal_subcol
    cout << "dscal_subcol: " << endl;
    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
	m1[i][j]=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      }
    }
    for(size_t k=0;k<5;k++) {
      for(size_t ell=0;ell<5;ell++) {
	ubmatrix_col mcol1(m1,k);
	ubvector_subvector vsub1(mcol1,ell,5-ell);
	o2scl_cblas::dscal(0.5,5-ell,vsub1);
	o2scl_cblas::dscal_subcol(m1,ell,k,5,0.5);
	for(size_t ij=0;ij<5-ell;ij++) {
	  t.test_rel(vsub1[ij],m1[ij+ell][k],1.0e-12,
		     "dscal_subcol operator[] ");
	}
      }
    }

    // daxpy_subcol
    cout << "daxpy_subcol: " << endl;
    for(size_t i=0;i<5;i++) {
      v1[i]=sin((double)(i+1));
      for(size_t j=0;j<5;j++) {
	m1[i][j]=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      }
    }
    for(size_t k=0;k<5;k++) {
      for(size_t ell=0;ell<5;ell++) {
	ubmatrix_col mcol1(m1,k);
	ubvector_subvector vsub1(mcol1,ell,5-ell);
	o2scl_cblas::daxpy(0.5,5-ell,vsub1,v1);
	o2scl_cblas::daxpy_subcol(0.5,5,m1,ell,k,v1);
	for(size_t ij=0;ij<5-ell;ij++) {
	  t.test_rel(vsub1[ij],m1[ij+ell][k],1.0e-12,
		     "daxpy_subcol operator[] ");
	}
      }
    }

    // ddot_subcol
    cout << "ddot_subcol: " << endl;
    for(size_t i=0;i<5;i++) {
      v1[i]=sin((double)(i+1));
      for(size_t j=0;j<5;j++) {
	m1[i][j]=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      }
    }
    for(size_t k=0;k<5;k++) {
      for(size_t ell=0;ell<5;ell++) {
	ubmatrix_col mcol1(m1,k);
	ubvector_subvector vsub1(mcol1,ell,5-ell);
	o2scl_cblas::ddot(5-ell,vsub1,v1);
	o2scl_cblas::ddot_subcol(5,m1,ell,k,v1);
	for(size_t ij=0;ij<5-ell;ij++) {
	  t.test_rel(vsub1[ij],m1[ij+ell][k],1.0e-12,
		     "ddot_subcol operator[] ");
	}
      }
    }

    // dnrm2_subrow
    cout << "dnrm2_subrow: " << endl;
    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
	m1[i][j]=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      }
    }
    for(size_t k=0;k<5;k++) {
      for(size_t ell=0;ell<5;ell++) {
	ubmatrix_row mrow1(m1,k);
	ubvector_subvector vsub1(mrow1,ell,5-ell);
	double tx1=o2scl_cblas::dnrm2(5-ell,vsub1);
	double tx2=o2scl_cblas::dnrm2_subrow(m1,k,ell,5);
	t.test_rel(tx1,tx2,1.0e-12,"dnrm2_subrow operator[] ");
      }
    }

    // dscal_subrow
    cout << "dscal_subrow: " << endl;
    for(size_t i=0;i<5;i++) {
      for(size_t j=0;j<5;j++) {
	m1[i][j]=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      }
    }
    for(size_t k=0;k<5;k++) {
      for(size_t ell=0;ell<5;ell++) {
	ubmatrix_row mrow1(m1,k);
	ubvector_subvector vsub1(mrow1,ell,5-ell);
	o2scl_cblas::dscal(0.5,5-ell,vsub1);
	o2scl_cblas::dscal_subrow(m1,k,ell,5,0.5);
	for(size_t ij=0;ij<5-ell;ij++) {
	  t.test_rel(vsub1[ij],m1[k][ij+ell],1.0e-12,
		     "dscal_subrow operator[] ");
	}
      }
    }

    // daxpy_subrow
    cout << "daxpy_subrow: " << endl;
    for(size_t i=0;i<5;i++) {
      v1[i]=sin((double)(i+1));
      for(size_t j=0;j<5;j++) {
	m1[i][j]=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      }
    }
    for(size_t k=0;k<5;k++) {
      for(size_t ell=0;ell<5;ell++) {
	ubmatrix_row mcol1(m1,k);
	ubvector_subvector vsub1(mcol1,ell,5-ell);
	o2scl_cblas::daxpy(0.5,5-ell,vsub1,v1);
	o2scl_cblas::daxpy_subrow(0.5,5,m1,k,ell,v1);
	for(size_t ij=0;ij<5-ell;ij++) {
	  t.test_rel(vsub1[ij],m1[k][ij+ell],1.0e-12,
		     "daxpy_subrow operator[] ");
	}
      }
    }

    // ddot_subrow
    cout << "ddot_subrow: " << endl;
    for(size_t i=0;i<5;i++) {
      v1[i]=sin((double)(i+1));
      for(size_t j=0;j<5;j++) {
	m1[i][j]=sqrt(((double)(i+1.0)))*((double)(j+1.0));
      }
    }
    for(size_t k=0;k<5;k++) {
      for(size_t ell=0;ell<5;ell++) {
	ubmatrix_row mcol1(m1,k);
	ubvector_subvector vsub1(mcol1,ell,5-ell);
	o2scl_cblas::ddot(5-ell,vsub1,v1);
	o2scl_cblas::ddot_subrow(5,m1,k,ell,v1);
	for(size_t ij=0;ij<5-ell;ij++) {
	  t.test_rel(vsub1[ij],m1[k][ij+ell],1.0e-12,
		     "ddot_subrow operator[] ");
	}
      }
    }

#endif

    cout << endl;
  }


  t.report();
  return 0;
}

