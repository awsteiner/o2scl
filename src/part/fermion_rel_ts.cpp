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
#include <o2scl/fermion_rel.h>
#include <o2scl/fermion_eff.h>
#include <o2scl/test_mgr.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/eos_sn.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
typedef boost::multiprecision::number<
  boost::multiprecision::cpp_dec_float<35> > cpp_dec_float_35;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);
  part_calibrate_class pcc;

  fermion f(1.0,2.0);
  fermion f2(1.0,2.0);
  fermion_rel fr;

  // AWS, 11/1/21: Note that massless fermions aren't necessary tested
  // by the particle calibrate classes so we test this separately here
  
  cout << "----------------------------------------------------" << endl;
  cout << "Test massless fermions." << endl;
  cout << "----------------------------------------------------" << endl;
  
  if (true) {
    
    double temper=3.0;
    f.m=1.0e-6;
    f2.m=1.0e-6;

    f.inc_rest_mass=true;
    f2.inc_rest_mass=true;
    cout << "Variable inc_rest_mass true.\n" << endl;
    
    cout << "mu=2.5: " << endl;
    f.mu=2.5;
    f2.mu=2.5;
    fr.calc_mu(f,temper);
    cout << f.n << " " << f.ed << " " << f.pr << " " << f.en << endl;
    fr.massless_calc_mu(f2,temper);
    cout << f2.n << " " << f2.ed << " " << f2.pr << " " << f2.en << endl;
    t.test_rel(f.n,f2.n,1.0e-6,"massless");
    t.test_rel(f.pr,f2.pr,1.0e-6,"massless");
    t.test_rel(f.ed,f2.ed,1.0e-6,"massless");
    fr.pair_mu(f,temper);
    cout << f.n << " " << f.ed << " " << f.pr << " " << f.en << endl;
    fr.massless_pair_mu(f2,temper);
    cout << f2.n << " " << f2.ed << " " << f2.pr << " " << f2.en << endl;
    t.test_rel(f.n,f2.n,1.0e-6,"massless");
    t.test_rel(f.pr,f2.pr,1.0e-6,"massless");
    t.test_rel(f.ed,f2.ed,1.0e-6,"massless");
    cout << endl;
    
    f.mu=0.0;
    f2.mu=0.0;
    cout << "mu=0.0: " << endl;
    fr.calc_mu(f,temper);
    cout << f.n << " " << f.ed << " " << f.pr << " " << f.en << endl;
    fr.massless_calc_mu(f2,temper);
    cout << f2.n << " " << f2.ed << " " << f2.pr << " " << f2.en << endl;
    t.test_rel(f.n,f2.n,1.0e-6,"massless");
    t.test_rel(f.pr,f2.pr,1.0e-6,"massless");
    t.test_rel(f.ed,f2.ed,1.0e-6,"massless");
    fr.pair_mu(f,temper);
    cout << f.n << " " << f.ed << " " << f.pr << " " << f.en << endl;
    fr.massless_pair_mu(f2,temper);
    cout << f2.n << " " << f2.ed << " " << f2.pr << " " << f2.en << endl;
    t.test_rel(f2.n,f2.n,1.0e-6,"massless");
    t.test_rel(f.pr,f2.pr,1.0e-6,"massless");
    t.test_rel(f.ed,f2.ed,1.0e-6,"massless");
    cout << endl;
    
    f.inc_rest_mass=false;
    f2.inc_rest_mass=false;
    cout << "Variable inc_rest_mass false.\n" << endl;
    
    cout << "mu=2.5: " << endl;
    f.mu=2.5;
    f2.mu=2.5;
    fr.calc_mu(f,temper);
    cout << f.n << " " << f.ed << " " << f.pr << " " << f.en << endl;
    fr.massless_calc_mu(f2,temper);
    cout << f2.n << " " << f2.ed << " " << f2.pr << " " << f2.en << endl;
    t.test_rel(f.n,f2.n,1.0e-6,"massless");
    t.test_rel(f.pr,f2.pr,1.0e-6,"massless");
    t.test_rel(f.ed,f2.ed,1.0e-6,"massless");
    fr.pair_mu(f,temper);
    cout << f.n << " " << f.ed << " " << f.pr << " " << f.en << endl;
    fr.massless_pair_mu(f2,temper);
    cout << f2.n << " " << f2.ed << " " << f2.pr << " " << f2.en << endl;
    t.test_rel(f.n,f2.n,1.0e-6,"massless");
    t.test_rel(f.pr,f2.pr,1.0e-6,"massless");
    t.test_rel(f.ed,f2.ed,1.0e-6,"massless");
    cout << endl;
    
    f.mu=0.0;
    f2.mu=0.0;
    cout << "mu=0.0: " << endl;
    fr.calc_mu(f,temper);
    cout << f.n << " " << f.ed << " " << f.pr << " " << f.en << endl;
    fr.massless_calc_mu(f2,temper);
    cout << f2.n << " " << f2.ed << " " << f2.pr << " " << f2.en << endl;
    t.test_rel(f.n,f2.n,1.0e-6,"massless");
    t.test_rel(f.pr,f2.pr,1.0e-6,"massless");
    t.test_rel(f.ed,f2.ed,1.0e-6,"massless");
    fr.pair_mu(f,temper);
    cout << f.n << " " << f.ed << " " << f.pr << " " << f.en << endl;
    fr.massless_pair_mu(f2,temper);
    cout << f2.n << " " << f2.ed << " " << f2.pr << " " << f2.en << endl;
    t.test_rel(f.n,f2.n,5.0e-6,"massless");
    t.test_rel(f.pr,f2.pr,1.0e-6,"massless");
    t.test_rel(f.ed,f2.ed,1.0e-6,"massless");
    cout << endl;
    
    f.m=1.0;
    f2.m=1.0;
  }
  
  cout << "----------------------------------------------------" << endl;
  cout << "Function calibrate()." << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  fr.verify_ti=true;
  double v1=pcc.part_calibrate<fermion,fermion_rel>
    (f,fr,true,"../../data/o2scl/fermion_deriv_cal.o2",false,0,true);
  t.test_rel(v1,0.0,4.0e-6,"calibrate");
  cout << endl;

  cout << "----------------------------------------------------" << endl;
  cout << "Function calibrate() with better limits." << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  // These seem to improve the accuracy. It's not clear that more
  // stringent tolerances will improve results.
  fr.upper_limit_fac=40.0;
  fr.fri.dit.tol_abs=1.0e-13;
  fr.fri.dit.tol_rel=1.0e-13;
  fr.fri.nit.tol_abs=1.0e-13;
  fr.fri.nit.tol_rel=1.0e-13;
  fr.density_root->tol_rel=1.0e-10;

  double v2=pcc.part_calibrate<fermion,fermion_rel>
    (f,fr,1,"../../data/o2scl/fermion_deriv_cal.o2",false,0,true);
  t.test_rel(v2,0.0,4.0e-10,"calibrate 2");

  cout << endl;
  
  fermion_ld fld;
  fermion_rel_ld frld;
  fermion_cdf25 f25;
  fermion_rel_cdf25 fr25;

#ifndef O2SCL_FAST_TEST
#ifdef O2SCL_NEVER_DEFINED
  
  cout << "----------------------------------------------------" << endl;
  cout << "Testing multiprecision" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  // AWS 6/17/22: this doesn't quite work, it fails on calc_density()

  frld.verify_ti=true;
  long double v3=pcc.part_calibrate<fermion_ld,fermion_rel_ld>
    (fld,frld,1,"../../data/o2scl/fermion_deriv_cal.o2",false,2,true);
  t.test_rel<long double>(v3,0.0,4.0e-10,"calibrate 3");

  //fr25.verify_ti=true;
  //cpp_dec_float_25 v4=pcc.part_calibrate<fermion_cdf25,fermion_rel_cdf25>
  //(f25,fr25,1,"../../data/o2scl/fermion_deriv_cal.o2",false,2,true);
  //t.test_rel<cpp_dec_float_25>(v3,0.0,4.0e-10,"calibrate 3");

#endif
#endif

#ifdef O2SCL_NEVER_DEFINED
  hdf_file hfx;
  hfx.open("../../data/o2scl/fermion_deriv_cal.o2");
  table_units<> tx;
  hdf_input(hfx,tx);
  hfx.close();
  f.non_interacting=true;
  f.inc_rest_mass=true;
  fld.non_interacting=true;
  fld.inc_rest_mass=true;
  f25.non_interacting=true;
  f25.inc_rest_mass=true;
  double T=1;
  long double Tld=1;
  long double T25=1;
  int ret;
  for(int lmot=-3;lmot<=3;lmot++) {
    for(int lpsi=-3;lpsi<=1;lpsi++) {
      
      double psi=pow(10.0,((double)lpsi));
      double mot=pow(10.0,((double)lmot));
      long double psi_ld=pow(10,((long double)lpsi));
      long double mot_ld=pow(10,((long double)lmot));
      cpp_dec_float_25 psi_25=pow(10,((cpp_dec_float_25)lpsi));
      cpp_dec_float_25 mot_25=pow(10,((cpp_dec_float_25)lmot));
      
      f.g=2;
      f.m=mot*T;
      f.mu=psi*T+f.m;
      ret=fr.calc_mu(f,T);
      cout << lmot << " " << lpsi << endl;
      fld.g=2;
      fld.m=mot_ld*Tld;
      fld.mu=psi_ld*Tld+fld.m;
      ret=frld.calc_mu(fld,Tld);
      f25.g=2;
      f25.m=mot_25*T25;
      f25.mu=psi_25*T25+f25.m;
      ret=fr25.calc_mu(f25,T25);
      cout << dtos(f.n,0) << " " << abs(f.n-fld.n)/abs(fld.n) << endl;
      cout << dtos(f.ed,0) << " " << abs(f.ed-fld.ed)/abs(fld.ed) << endl;
      cout << dtos(f.pr,0) << " " << abs(f.pr-fld.pr)/abs(fld.pr) << endl;
      cout << dtos(f.en,0) << " " << abs(f.en-fld.en)/abs(fld.en) << endl;
      cout << dtos(fld.n,0) << " " << abs(fld.n-f25.n)/abs(f25.n) << endl;
      cout << dtos(fld.ed,0) << " " << abs(fld.ed-f25.ed)/abs(f25.ed) << endl;
      cout << dtos(fld.pr,0) << " " << abs(fld.pr-f25.pr)/abs(f25.pr) << endl;
      cout << dtos(fld.en,0) << " " << abs(fld.en-f25.en)/abs(f25.en) << endl;
      cout << dtos(f25.n,0) << endl;
      cout << dtos(f25.ed,0) << endl;
      cout << dtos(f25.pr,0) << endl;
      cout << dtos(f25.en,0) << endl;
    }
  }
#endif
  
  t.report();

  return 0;
}

