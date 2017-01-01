/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2017, Andrew W. Steiner
  
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
#include <o2scl/vector_derint.h>
#include <o2scl/columnify.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;

int main(void) {
  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

  static const size_t n=50;

  typedef boost::numeric::ublas::vector<double> ubvector;
  ubvector sinx(n), cosx(n);
  for(size_t i=0;i<n;i++) {
    sinx[i]=sin(((double)i)/4.0);
    cosx[i]=cos(((double)i)/4.0);
  }

  {
    ubvector t1(n), t2(n), t3(n), t4(n), t5(n), t6(n);
    vector_deriv_threept(n,sinx,t1);
    vector_deriv_threept_tap(n,sinx,t2);
    vector_deriv_fivept(n,sinx,t3);
    vector_deriv_fivept_tap(n,sinx,t4);
    vector_deriv_interp(n,sinx,t5,itp_cspline);
    vector_deriv_interp(n,sinx,t6,itp_akima);

    cout << "Differentiation: " << endl;
    cout.setf(ios::showpos);
    cout.precision(4);
    for(size_t i=0;i<n;i++) {
      cout << fabs(t1[i]-0.25*cosx[i]) << " " 
	   << fabs(t2[i]-0.25*cosx[i]) << " " 
	   << fabs(t3[i]-0.25*cosx[i]) << " " 
	   << fabs(t4[i]-0.25*cosx[i]) << " " 
	   << fabs(t5[i]-0.25*cosx[i]) << " " 
	   << fabs(t6[i]-0.25*cosx[i]) << endl;
      t.test_abs(fabs(t1[i]-0.25*cosx[i]),0.0,7.0e-3,"threept");
      t.test_abs(fabs(t2[i]-0.25*cosx[i]),0.0,2.0e-2,"threept_tap");
      t.test_abs(fabs(t3[i]-0.25*cosx[i]),0.0,4.0e-4,"fivept");
      t.test_abs(fabs(t4[i]-0.25*cosx[i]),0.0,7.0e-3,"fivept_tap");
      t.test_abs(fabs(t5[i]-0.25*cosx[i]),0.0,8.0e-3,"cspline");
      t.test_abs(fabs(t6[i]-0.25*cosx[i]),0.0,7.0e-3,"akima");
    }
    cout << endl;
    cout.precision(6);
    cout.unsetf(ios::showpos);
  }
  
  cout << "Integration: " << endl;
  double exact=(1.0-cos(((double)(n-1))/4.0));
  double t1=fabs(exact-0.25*vector_integ_trap(n,sinx));
  t.test_abs(t1,0.0,4.0e-4,"trap");
  cout << t1 << endl;
  double t2=fabs(exact-0.25*vector_integ_threept(n,sinx));
  t.test_abs(t2,0.0,4.0e-4,"threept");
  cout << t2 << endl;
  double t3=fabs(exact-0.25*vector_integ_extended4(n,sinx));
  t.test_abs(t3,0.0,3.0e-4,"extended4");
  cout << t3 << endl;
  double t4=fabs(exact-0.25*vector_integ_durand(n,sinx));
  t.test_abs(t4,0.0,4.0e-4,"durand");
  cout << t4 << endl;
  double t5=fabs(exact-0.25*vector_integ_extended8(n,sinx));
  t.test_abs(t5,0.0,3.0e-6,"extended8");
  cout << t5 << endl;
  double t6=fabs(exact-0.25*vector_integ_interp(n,sinx,itp_cspline));
  t.test_abs(t6,0.0,3.0e-4,"cspline");
  cout << t6 << endl;
  double t7=fabs(exact-0.25*vector_integ_interp(n,sinx,itp_akima));
  t.test_abs(t7,0.0,3.0e-5,"akima");
  cout << t7 << endl;
  
  t.report();
  return 0;
}

