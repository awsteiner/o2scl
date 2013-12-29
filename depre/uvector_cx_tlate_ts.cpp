/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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

#include <complex>

#include <gsl/gsl_complex_math.h>

#include <o2scl/uvector_cx_tlate.h>
#include <o2scl/cx_arith.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {
  uvector_cx a(3);
 
  test_mgr t;
  t.set_output_level(1);
  
  a.set(0,1.0,2.0);
  a.set(1,2.0,3.0);
  a.set(2,-0.323529,-0.794118);
  
  gsl_complex ac=a.get(0), bc=a.get(1), cc=a.get(2), c;
  
  gsl_complex s=gsl_complex_div(gsl_complex_mul(ac,bc),
                                gsl_complex_add(ac,bc));
  c=gsl_complex_sub(s,ac);
  
  t.test_rel(GSL_REAL(c),GSL_REAL(cc),1.0e-5,"real");
  t.test_rel(GSL_IMAG(c),GSL_IMAG(cc),1.0e-5,"imag");
  
  std::complex<double> sa=a.get_stl(0), sb=a.get_stl(1), sc=a.get_stl(2), sd;
  
  sd=(sa*sb)/(sa+sb)-sa;
  
  t.test_rel(GSL_REAL(c),sd.real(),1.0e-5,"real stl");
  t.test_rel(GSL_IMAG(c),sd.imag(),1.0e-5,"imag stl");
  
  for(size_t j=0;j<3;j++) {
    gsl_complex x=a.get(j);
    t.test_rel(a.real(j),GSL_REAL(x),1.0e-12,"get_real");
    t.test_rel(a.imag(j),GSL_IMAG(x),1.0e-12,"get_imag");
  }
  
  uvector_cx_subvector sub(a,1,2);
  for(size_t k=0;k<2;k++) {
    gsl_complex g1=a.get(k+1), g2=sub.get(k);
    t.test_rel(GSL_REAL(g1),GSL_REAL(g2),1.0e-6,"sub real");
    t.test_rel(GSL_IMAG(g1),GSL_IMAG(g2),1.0e-6,"sub imag");
  }

  std::complex<double> ccx(1.0,2.0), ccx2;
  gsl_complex ggx=complex_to_gsl(ccx);
  ccx2=gsl_to_complex(ggx);
  t.test_rel(ccx.real(),ccx2.real(),1.0e-10,"convert");
  t.test_rel(ccx.imag(),ccx2.imag(),1.0e-10,"convert");
  
  t.report();
  return 0;
}
