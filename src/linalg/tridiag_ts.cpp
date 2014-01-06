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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_linalg.h>

#include <o2scl/tridiag.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  void *vpx=0;
  size_t tmp;

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

  static const size_t arr_size=50;

  {
    using namespace o2scl_cblas;
    using namespace o2scl_linalg;

    // GSL
    gsl_vector *d=gsl_vector_alloc(arr_size);
    gsl_vector *od=gsl_vector_alloc(arr_size-1);
    gsl_vector *r=gsl_vector_alloc(arr_size);
    gsl_vector *s=gsl_vector_alloc(arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      gsl_vector_set(d,i,i+1.0);
      if (i!=arr_size-1) gsl_vector_set(od,i,1.0);
      gsl_vector_set(r,i,1.0);
      gsl_vector_set(s,i,0.0);
    }
    int ret=gsl_linalg_solve_symm_tridiag(d,od,r,s);

    // O2scl with ubvectors
    ubvector du(arr_size), odu(arr_size-1), ru(arr_size), su(arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      du[i]=i+1.0;
      if (i!=arr_size-1) odu[i]=1.0;
      ru[i]=1.0;
      su[i]=0.0;
    }
    solve_tridiag_sym(du,odu,ru,su,arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      t.test_rel(gsl_vector_get(s,i),su[i],1.0e-12,"o vs. u 5");
    }
  }

  {
    using namespace o2scl_cblas;
    using namespace o2scl_linalg;

    // GSL
    gsl_vector *d=gsl_vector_alloc(arr_size);
    gsl_vector *od=gsl_vector_alloc(arr_size-1);
    gsl_vector *od2=gsl_vector_alloc(arr_size-1);
    gsl_vector *r=gsl_vector_alloc(arr_size);
    gsl_vector *s=gsl_vector_alloc(arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      gsl_vector_set(d,i,i+1.0);
      if (i!=arr_size-1) gsl_vector_set(od,i,1.0);
      if (i!=arr_size-1) gsl_vector_set(od2,i,1.0);
      gsl_vector_set(r,i,1.0);
      gsl_vector_set(s,i,0.0);
    }
    int ret=gsl_linalg_solve_tridiag(d,od,od2,r,s);

    // O2scl with ubvectors
    ubvector du(arr_size), odu(arr_size-1), ru(arr_size), su(arr_size);
    ubvector odu2(arr_size-1);
  
    for(size_t i=0;i<arr_size;i++) {
      du[i]=i+1.0;
      if (i!=arr_size-1) odu[i]=1.0;
      if (i!=arr_size-1) odu2[i]=1.0;
      ru[i]=1.0;
      su[i]=0.0;
    }
    solve_tridiag_nonsym(du,odu,odu2,ru,su,arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      t.test_rel(gsl_vector_get(s,i),su[i],1.0e-12,"o vs. u 6");
    }
  }

  {
    using namespace o2scl_cblas;
    using namespace o2scl_linalg;

    // GSL
    gsl_vector *d=gsl_vector_alloc(arr_size);
    gsl_vector *od=gsl_vector_alloc(arr_size);
    gsl_vector *r=gsl_vector_alloc(arr_size);
    gsl_vector *s=gsl_vector_alloc(arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      gsl_vector_set(d,i,i+1.0);
      gsl_vector_set(od,i,1.0);
      gsl_vector_set(r,i,1.0);
      gsl_vector_set(s,i,0.0);
    }
    int ret=gsl_linalg_solve_symm_cyc_tridiag(d,od,r,s);

    // O2scl with ubvectors
    ubvector du(arr_size), odu(arr_size), ru(arr_size), su(arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      du[i]=i+1.0;
      odu[i]=1.0;
      ru[i]=1.0;
      su[i]=0.0;
    }
    solve_cyc_tridiag_sym(du,odu,ru,su,arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      t.test_rel(gsl_vector_get(s,i),su[i],1.0e-12,"o vs. u 7");
    }
  }

  {
    using namespace o2scl_cblas;
    using namespace o2scl_linalg;

    // GSL
    gsl_vector *d=gsl_vector_alloc(arr_size);
    gsl_vector *od=gsl_vector_alloc(arr_size);
    gsl_vector *od2=gsl_vector_alloc(arr_size);
    gsl_vector *r=gsl_vector_alloc(arr_size);
    gsl_vector *s=gsl_vector_alloc(arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      gsl_vector_set(d,i,i+1.0);
      gsl_vector_set(od,i,1.0);
      gsl_vector_set(od2,i,1.0);
      gsl_vector_set(r,i,1.0);
      gsl_vector_set(s,i,0.0);
    }
    int ret=gsl_linalg_solve_cyc_tridiag(d,od,od2,r,s);

    // O2scl with ubvectors
    ubvector du(arr_size), odu(arr_size), ru(arr_size), su(arr_size);
    ubvector odu2(arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      du[i]=i+1.0;
      odu[i]=1.0;
      odu2[i]=1.0;
      ru[i]=1.0;
      su[i]=0.0;
    }
    solve_cyc_tridiag_nonsym(du,odu,odu2,ru,su,arr_size);
  
    for(size_t i=0;i<arr_size;i++) {
      t.test_rel(gsl_vector_get(s,i),su[i],1.0e-12,"o vs. u 8");
    }
  }

  t.report();
  return 0;
}

