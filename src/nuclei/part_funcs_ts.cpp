/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <o2scl/test_mgr.h>
#include <o2scl/part_funcs.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;
using namespace o2scl_hdf;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);

  part_funcs pf;
  pf.load("../../data/o2scl/nucmass/");

  double xpf, TdpfdT;
  int Z, N;
  double T_K=1.15011964639613e+10;
  double T=5.02260981304584e-03;
  double v, vop;
  
  Z=2;
  N=1;
  pf.few78(Z,N,T_K,v,vop);
  vop/=T;
  t.test_rel(v,2.0,1.0e-6,"1");
  t.test_abs(vop,0.0,1.0e-6,"2");

  Z=41;
  N=74;
  pf.few78(Z,N,T_K,v,vop);
  vop/=T;
  t.test_rel(v,1.957695e+03,1.0e-6,"3");
  t.test_rel(vop,1.425872e+06,1.0e-6,"4");

  Z=78;
  N=109;
  pf.few78(Z,N,T_K,v,vop);
  vop/=T;
  t.test_rel(v,9.484188e+05,1.0e-6,"5");
  t.test_rel(vop,1.110043e+09,1.0e-6,"6");

  Z=99;
  N=152;
  pf.few78(Z,N,T_K,v,vop);
  vop/=T;
  t.test_rel(v,3.40619635e+07,1.0e-6,"7");
  t.test_rel(vop,4.07247895e+10,1.0e-6,"8");

  Z=111;
  N=201;
  pf.few78(Z,N,T_K,v,vop);
  vop/=T;
  t.test_rel(v,1.434830e+07,1.0e-6,"9");
  t.test_rel(vop,8.158737e+09,1.0e-6,"10");
  
  pf.rt00(28,28,5.8e9,xpf,TdpfdT);
  cout << xpf << " " << TdpfdT << endl;
  pf.r03(28,28,28.8e9,xpf,TdpfdT);
  cout << xpf << " " << TdpfdT << endl;

  t.report();
  return 0;
}

