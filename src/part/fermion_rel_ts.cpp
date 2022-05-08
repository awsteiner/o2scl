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

  /*
    AWS, 5/6/22: I'm commenting this out for now as it gives 
    identical results to the default
    
    typedef fermion_rel_tl<fermion_tl<double>,o2scl::fermi_dirac_integ_gsl,
    bessel_K_exp_integ_boost<double>,
    fermion_rel_integ<funct,double>,
    root_cern<>,root_cern<>,funct,double> fr2_t;
    
    fr2_t fr2;
    
    fr2.verify_ti=true;
    double v1x=pcc.part_calibrate<fermion,fr2_t>
    (f,fr2,true,"../../data/o2scl/fermion_deriv_cal.o2",false,0,true);
    t.test_rel(v1x,0.0,4.0e-6,"calibrate x");
  */

  /*
  cout << fr.upper_limit_fac << endl;
  cout << fr.fri.dit.tol_abs << endl;
  cout << fr.fri.dit.tol_rel << endl;
  cout << fr.fri.nit.tol_abs << endl;
  cout << fr.fri.nit.tol_rel << endl;
  cout << fr.density_root->tol_rel << endl;
  */
  
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
  cout << "----------------------------------------------------" << endl;
  cout << "Testing multiprecision" << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  fermion_ld fld;
  fermion_rel_ld frld;

  cout << "Non-degnerate test:" << endl;
  f.inc_rest_mass=true;
  f.non_interacting=true;
  f.mu=1.5;
  f.m=1.0;
  f.g=2.0;
  f.m=1.0;
  fr.calc_mu(f,0.3);
  cout << "double: " << endl;
  cout << f.n << " " << f.ed << endl;
  cout << f.pr << " " << f.en << endl;
  cout << endl;
  
  fld.m=1;
  fld.g=2;
  fld.mu=3;
  fld.mu/=2;
  long double Tld=3;
  Tld/=10;
  frld.verify_ti=true;
  frld.calc_mu(fld,Tld);
  cout << "long double: " << frld.last_method << endl;
  cout << "    n,ed: " << dtos(fld.n,0) << " " << dtos(fld.ed,0) << endl;
  cout << "   pr,en: " << dtos(fld.pr,0) << " " << dtos(fld.en,0) << endl;
  cout << "pr check: " << dtos(-fld.ed+fld.n*fld.mu+Tld*fld.en,0) << endl;
  cout << endl;

  cout << "Degnerate test:" << endl;
  f.mu=15.0;
  f.m=1.0;
  f.g=2.0;
  f.m=1.0;
  fr.calc_mu(f,0.3);
  cout << "double: " << endl;
  cout << f.n << " " << f.ed << endl;
  cout << f.pr << " " << f.en << endl;
  cout << endl;

  fld.mu=15;
  frld.fri.dit25.verbose=1;
  frld.fri.dit35.verbose=1;
  frld.fri.dit50.verbose=1;
  frld.calc_mu(fld,Tld);
  cout << "long double: " << frld.last_method << endl;
  cout << dtos(fld.n,0) << " " << dtos(fld.ed,0) << endl;
  cout << dtos(fld.pr,0) << " " << dtos(fld.en,0) << endl;
  cout << dtos(-fld.ed+fld.n*fld.mu+Tld*fld.en,0) << endl;
  cout << endl;

  // AWS 5/6/22: this doesn't quite work, it fails on one of
  // the degenerate entropy integrals

  if (1) {
    frld.verbose=2;
    long double v3=pcc.part_calibrate<fermion_ld,fermion_rel_ld>
      (fld,frld,1,"../../data/o2scl/fermion_deriv_cal.o2",false,1,true);
    t.test_rel<long double>(v3,0.0,4.0e-10,"calibrate 3");
  }

  t.report();

  return 0;
}

