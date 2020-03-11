/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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

#include <o2scl/test_mgr.h>
#include <o2scl/inte_adapt_cern.h>
#include <o2scl/polylog.h>

using namespace std;
using namespace o2scl;

int main(void) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(2);
  
#ifdef O2SCL_LD_TYPES

#ifdef O2SCL_NEW_BOOST_INTEGRATION
  
  fermi_dirac_integ_gsl f1;
  fermi_dirac_integ_direct f2;

  t.test_rel(f1.calc_1o2(0.5),f2.calc_1o2(0.5),4.0e-16,"fd 1");
  t.test_rel(f1.calc_m1o2(0.5),f2.calc_m1o2(0.5),4.0e-16,"fd 2");
  t.test_rel(f1.calc_3o2(0.5),f2.calc_3o2(0.5),4.0e-16,"fd 3");
  t.test_rel(f1.calc_2(0.5),f2.calc_2(0.5),4.0e-16,"fd 4");
  t.test_rel(f1.calc_3(0.5),f2.calc_3(0.5),4.0e-16,"fd 5");

  polylog p;
  t.test_rel(p.calc(2,-0.5),-0.448414206923646,4.0e-15,"pl 1");
  t.test_rel(p.calc(2,-2.0),-1.43674636688368,4.0e-15,"pl 2");
  t.test_rel(p.calc(3,-0.5),-0.472597844658897,4.0e-15,"pl 3");
  t.test_rel(p.calc(3,-2.0),-1.66828336396657,4.0e-15,"pl 4");
  t.test_rel(p.calc(2,0.5),0.5822405264650125,4.0e-15,"pl 5");
  t.test_rel(p.calc(3,0.5),0.5372131936080402,4.0e-15,"pl 6");

  bessel_K_exp_integ_tl<o2scl::inte_exp_sinh_boost
			<funct_ld,15,long double>,long double> be;
  long double res, err;
  be.calc_err(2,2.0,res,err);
  t.test_rel<long double>(res,1.875045062139460,2.0e-16,"be 1");
  be.calc_err(2,20.0,res,err);
  t.test_rel<long double>(res,0.3070874263512549,2.0e-16,"be 2");

  bessel_K_exp_integ_gsl beg;
  bessel_K_exp_integ_direct bed;

  t.test_rel(beg.K1exp(2.0),bed.K1exp(2.0),1.0e-15,"bed 1");
  t.test_rel(beg.K2exp(2.0),bed.K2exp(2.0),1.0e-15,"bed 2");
  t.test_rel(beg.K3exp(2.0),bed.K3exp(2.0),1.0e-15,"bed 3");
  
#endif
  
#endif
  
  t.report();
  return 0;
}
