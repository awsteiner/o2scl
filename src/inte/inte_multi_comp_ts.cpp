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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <o2scl/test_mgr.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_multi_comp.h>
#include <o2scl/multi_funct.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

double test_fun(size_t nv, const ubvector &x) {
  double tmp;
  for(size_t i=0;i<nv;i++) tmp+=x[i];
  return 1.0;
}

int main(void) {
  test_mgr t;
  t.set_output_level(2);

  /// The integrator
  inte_multi_comp<multi_funct<> > ci4;
  
  /// The individual integration objects
  inte_qag_gsl<funct11> gl4a;
  inte_qag_gsl<funct11> gl4b;
  inte_qag_gsl<funct11> gl4c;
  inte_qag_gsl<funct11> gl4d;

  /// Set the 1-d objects in the integrator
  ci4.set_oned_inte(gl4a,0);
  ci4.set_oned_inte(gl4b,1);
  ci4.set_oned_inte(gl4c,2);
  ci4.set_oned_inte(gl4d,3);
  
  int vp=0;
  int ndim,i;
  double res;
  ubvector ll, ul;

  cout.setf(ios::scientific);
  
  multi_funct_fptr<> ffn(test_fun);

  // Calculate the volume of a 'ndim'-dimensional hypercube with
  // sides of length 2.
  for(ndim=1;ndim<=4;ndim++) {
    ll.resize(ndim);
    ul.resize(ndim);
    for(i=0;i<ndim;i++) {
      ll[i]=0.0;
      ul[i]=2.0;
    }
    res=ci4.minteg(ffn,ndim,ll,ul);

    cout << ndim << " " << res << endl;

    t.test_rel(res,pow(2.0,ndim),1.0e-6,"inte_multi_comp");
    
    ll.clear();
    ul.clear();
  }

  t.report();
  return 0;
}
