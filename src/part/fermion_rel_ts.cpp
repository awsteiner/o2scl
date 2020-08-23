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
#include <o2scl/fermion_rel.h>
#include <o2scl/fermion_eff.h>
#include <o2scl/test_mgr.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/eos_sn.h>

#ifdef O2SCL_LD_TYPES
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

#ifdef O2SCL_LD_TYPES
typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
#endif

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);
  part_calibrate_class pcc;

  fermion f(1.0,2.0);
  fermion_rel fr;

  cout << "----------------------------------------------------" << endl;
  cout << "Function calibrate()." << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;
  
  double v1=pcc.part_calibrate<fermion,fermion_rel>
    (f,fr,true,"../../data/o2scl/fermion_deriv_cal.o2",false,1,true);
  t.test_rel(v1,0.0,4.0e-6,"calibrate");
  
  cout << "----------------------------------------------------" << endl;
  cout << "Function calibrate() with better limits." << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  // These seem to improve the accuracy. It's not clear that more
  // stringent tolerances will improve results.
  fr.upper_limit_fac=40.0;
  fr.dit->tol_abs=1.0e-13;
  fr.dit->tol_rel=1.0e-13;
  fr.nit->tol_abs=1.0e-13;
  fr.nit->tol_rel=1.0e-13;
  fr.density_root->tol_rel=1.0e-10;

  double v2=pcc.part_calibrate<fermion,fermion_rel>
    (f,fr,1,"../../data/o2scl/fermion_deriv_cal.o2",false,1,1);
  t.test_rel(v2,0.0,4.0e-10,"calibrate 2");

  // -----------------------------------------------------------------
  // Downcast the shared_ptr to the default integration type. This
  // shows how to get access the internal integration object that
  // fermion_rel is using.
  // 
  // From cppreference.com: "If the cast is successful, dynamic_cast
  // returns a value of type new_type. If the cast fails and new_type
  // is a pointer type, it returns a null pointer of that type. If the
  // cast fails and new_type is a reference type, it throws an
  // exception that matches a handler of type std::bad_cast."
  
  inte_qag_gsl<> *qag=dynamic_cast<inte_qag_gsl<> *>(fr.dit.get());
  inte_qag_gsl<> &qag2=dynamic_cast<inte_qag_gsl<> &>(*fr.dit.get());
  t.test_gen(qag->get_rule()==qag2.get_rule(),"downcast");
  
#ifdef O2SCL_LD_TYPES

  // These don't work yet. The next step is to make sure that
  // bessel_K_exp_integ_direct<long double,cpp_dec_float_50>
  // is tested in ../other/polylog_ts.cpp
  
  /*
    fermion_rel_tl<fermi_dirac_integ_direct<long double,cpp_dec_float_50>,
    bessel_K_exp_integ_direct<long double,cpp_dec_float_50>,
    long double> fermion_rel_ld;
  */
  //fermion_rel_tl<cpp_dec_float_50> fermion_rel_cdf;
  
#endif

  t.report();

  return 0;
}

