/*
  -------------------------------------------------------------------
  
  Copyright (C) 2009-2014, Marco Cammarata and Andrew W. Steiner
  
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
#include <o2scl/err_hnd.h>
#include <o2scl/smooth_gsl.h>

using namespace o2scl;
using namespace std;

void smooth_gsl::init_pointers_and_defs() {
  bw=0;
  B=0;
  c=0;
  mw=0;
  X=0;
  cov=0;
  ncoeffs=10; 
  norder=4;
  x_set=false;
}

int smooth_gsl::free(void) {
  if (bw!=0) gsl_bspline_free(bw);
  if (B!=0) gsl_vector_free(B);
  if (X!=0) gsl_matrix_free(X);
  if (c!=0) gsl_vector_free(c);
  if (cov!=0) gsl_matrix_free(cov);
  if (mw!=0) gsl_multifit_linear_free(mw);
  x_set=false;
  return success;
}

int smooth_gsl::init() {
  
  free(); 
  
  nbreak=ncoeffs+2-norder;
  
  size_t n=x->size;

  bw=gsl_bspline_alloc(norder, nbreak);
  B=gsl_vector_alloc(ncoeffs); 
  X=gsl_matrix_alloc(n,ncoeffs);
  c=gsl_vector_alloc(ncoeffs);
  cov=gsl_matrix_alloc(ncoeffs,ncoeffs);
  mw=gsl_multifit_linear_alloc(n,ncoeffs);
  
  gsl_bspline_knots_uniform(gsl_vector_get(x,0),gsl_vector_get(x,n-1),bw);
  
  // construct the matrix X
  for(size_t i=0;i<n;i++) {

    gsl_bspline_eval(gsl_vector_get(x,i),B,bw);

    for (size_t j=0;j<ncoeffs;++j) {
      double Bj=gsl_vector_get(B,j);
      gsl_matrix_set(X,i,j,Bj);
    }
  }

  x_set=true;

  return success;
}

int smooth_gsl::fit_errors(const gsl_vector *y, const gsl_vector *e) {

  double chisq;
  size_t n;
  n=x->size;

  gsl_vector *w;
  // generate weight vector (w= 1/e^2)
  w=gsl_vector_alloc(n);

  for (size_t i=0; i < n; i++) {
    gsl_vector_set(w,i,pow(gsl_vector_get(e,i),2));
  }

  gsl_multifit_wlinear(X,w,y,c,cov,&chisq,mw);

  gsl_vector_free(w);

  x_set=true;

  return success;
}

int smooth_gsl::fit(const gsl_vector *y) {
  double chisq;
  gsl_multifit_linear(X,y,c,cov,&chisq,mw);
  return success;
}

int smooth_gsl::smooth_data(const gsl_vector *y, 
			    gsl_vector *ys) {
  return smooth_data(y,0,ys);
}

int smooth_gsl::smooth_data(const gsl_vector *y, 
			    const gsl_vector *e, gsl_vector *ys) {
  
  
  if (!x_set) {
    O2SCL_ERR_RET("X values not set in smooth_gsl::smooth_data().",
		  exc_einval);
  }

  double xi, vs;
  size_t n;
  n=x->size;
  if (e == 0) {
    fit(y);
  } else {
    fit_errors(y,e);
  }

  for (size_t i=0; i < n; i++) {
    xi=gsl_vector_get(x,i);
    vs=calc_for_x(xi);
    gsl_vector_set(ys,i,vs);
  }

  return success;
}

double smooth_gsl::calc_for_x(double xi) {
  double ys,ye;
  gsl_bspline_eval(xi, B, bw);
  gsl_multifit_linear_est(B, c, cov, &ys, &ye);
  return ys;
}
