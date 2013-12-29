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

#include <gsl/gsl_complex_math.h>

#include <o2scl/test_mgr.h>
#include <o2scl/ovector_cx_tlate.h>
// AWS 10/22/2012: It seems to be important that cx_arith is here, and
// after ovector_cx_tlate.h, otherwise, installation fails in unusual
// places. It's probably a header file dependency issue. 
#include <o2scl/cx_arith.h>

using namespace std;
using namespace o2scl;

int main(void) {

  test_mgr t;
  t.set_output_level(1);
  
  ovector_cx a(3);
  double dx=5.0;

  // Test simple set operations
  a.set(0,1.0,2.0);
  a.set(1,2.0,3.0);
  a.set(2,-0.323529,-0.794118);

  // Test the deep copy constructor
  ovector_cx acopy=a;

  cout << a << endl;
  cout << acopy << endl;

  // Test the shallow copy constructor
  ovector_cx_view acopy2=a;

  cout << "a: " << a << endl;
  cout << "acopy2: " << acopy2 << endl;

  // Test += complex
  acopy+=a;

  cout << "a: " << a << endl;
  cout << "acopy: " << acopy << endl;

  // Test += double
  acopy2+=dx;

  cout << "a: " << a << endl;
  cout << "dx: " << dx << endl;
  cout << "acopy: " << acopy << endl;
  cout << "acopy2: " << acopy2 << endl;

  gsl_complex ac=acopy.get(0), bc=acopy.get(1), cc=acopy.get(2), c;

  cout << GSL_REAL(ac) << " " << GSL_IMAG(ac) << endl;
  cout << GSL_REAL(bc) << " " << GSL_IMAG(bc) << endl;
  cout << GSL_REAL(cc) << " " << GSL_IMAG(cc) << endl;
  
  gsl_complex s=gsl_complex_div(gsl_complex_mul(ac,bc),
                                gsl_complex_add(ac,bc));
  c=gsl_complex_sub(s,ac);
  
  t.test_rel(GSL_REAL(c),GSL_REAL(cc),1.0e-5,"real");
  t.test_rel(GSL_IMAG(c),GSL_IMAG(cc),1.0e-5,"imag");
  
  std::complex<double> sa=acopy.get_stl(0), sb=acopy.get_stl(1), 
    sc=acopy.get_stl(2), sd;
  
  sd=(sa*sb)/(sa+sb)-sa;
  
  t.test_rel(GSL_REAL(c),sd.real(),1.0e-5,"real stl");
  t.test_rel(GSL_IMAG(c),sd.imag(),1.0e-5,"imag stl");
  
  gsl_vector_complex *gvc=(gsl_vector_complex *)(&acopy);
  
  for(size_t i=0;i<3;i++) {
    gsl_complex gc=gsl_vector_complex_get(gvc,i);
    if (i==0) {
      t.test_rel(GSL_REAL(gc),GSL_REAL(ac),1.0e-4,"gsl");
      t.test_rel(GSL_IMAG(gc),GSL_IMAG(ac),1.0e-4,"gsl");
    } else if (i==1) {
      t.test_rel(GSL_REAL(gc),GSL_REAL(bc),1.0e-4,"gsl");
      t.test_rel(GSL_IMAG(gc),GSL_IMAG(bc),1.0e-4,"gsl");
    } else {
      t.test_rel(GSL_REAL(gc),GSL_REAL(cc),1.0e-4,"gsl");
      t.test_rel(GSL_IMAG(gc),GSL_IMAG(cc),1.0e-4,"gsl");
    }
  }
  
  for(size_t j=0;j<3;j++) {
    gsl_complex x=a.get(j);
    t.test_rel(a.real(j),GSL_REAL(x),1.0e-12,"get_real");
    t.test_rel(a.imag(j),GSL_IMAG(x),1.0e-12,"get_imag");
  }
  
  ovector_cx_subvector sub(a,1,2);
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

  a.conjugate();
  cout << a.real(0) << " " << a.imag(0) << endl;
  cout << a.real(1) << " " << a.imag(1) << endl;
  cout << a.real(2) << " " << a.imag(2) << endl;
  
  t.report();
  return 0;
}
