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
  
  double v1=pcc.part_calibrate<fermion,fermion_rel>
    (f,fr,true,"../../data/o2scl/fermion_deriv_cal.o2",false,1,true);
  t.test_rel(v1,0.0,4.0e-6,"calibrate");

  fermion_rel_tl<fermion_tl<double>,o2scl::fermi_dirac_integ_gsl,
                 bessel_K_exp_integ_boost<double>,inte_qagiu_gsl<>,
                 inte_qag_gsl<>,fermion_rel_integ<funct,double>,
                 root_cern<>,root_cern<>,funct,double> fr2;
  
  double v1x=pcc.part_calibrate<fermion,fermion_rel>
    (f,fr,true,"../../data/o2scl/fermion_deriv_cal.o2",false,1,true);
  t.test_rel(v1x,0.0,4.0e-6,"calibrate x");

  cout << fr.upper_limit_fac << endl;
  cout << fr.fri.dit.tol_abs << endl;
  cout << fr.fri.dit.tol_rel << endl;
  cout << fr.fri.nit.tol_abs << endl;
  cout << fr.fri.nit.tol_rel << endl;
  cout << fr.fri.dit.tol_abs << endl;
  cout << fr.fri.dit.tol_rel << endl;
  cout << fr.fri.nit.tol_abs << endl;
  cout << fr.fri.nit.tol_rel << endl;
  cout << fr.density_root->tol_rel << endl;
  
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
  fr.fri.dit.tol_abs=1.0e-13;
  fr.fri.dit.tol_rel=1.0e-13;
  fr.fri.nit.tol_abs=1.0e-13;
  fr.fri.nit.tol_rel=1.0e-13;
  fr.density_root->tol_rel=1.0e-10;

  double v2=pcc.part_calibrate<fermion,fermion_rel>
    (f,fr,1,"../../data/o2scl/fermion_deriv_cal.o2",false,1,1);
  t.test_rel(v2,0.0,4.0e-10,"calibrate 2");

  // AWS, 11/1/21: taking the multiprecision types out while I
  // develop new fermion integrators.
#ifdef O2SCL_NEVER_DEFINED
  
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
  
  //inte_qag_gsl<> *qag=dynamic_cast<inte_qag_gsl<> *>(fr.dit.get());
  //inte_qag_gsl<> &qag2=dynamic_cast<inte_qag_gsl<> &>(*fr.dit.get());
  //t.test_gen(qag->get_rule()==qag2.get_rule(),"downcast");
  
  // The long double type isn't that much more precise than double, so
  // I'm concerned about precision loss in the integrands. Thus, the
  // plan is to have a particle type which operates in
  // cpp_dec_float_35 precision, with integrators which internally use
  // cpp_dec_float_50 precision. For the integrators, the template
  // parameter is the maximum number of refinements, we try 30.
  
  // One limitation for using multiprecision types is the lack of a
  // systematic expansion for massless fermions in fermion.h

  fermion_ld fld;
  fermion_rel_ld frld;
  fermion_rel_ld2 frld2;
  
  fermion_cdf35 f35;
  fermion_rel_cdf35 fr35;

  f.m=1;
  f.g=2;
  f.mu=3;
  f.mu/=2;
  double T=3;
  T/=10;
  fr.calc_mu(f,T);
  cout << "double:" << endl;
  cout << dtos(f.n,0) << " " << dtos(f.ed,0) << endl;
  cout << dtos(f.pr,0) << " " << dtos(f.en,0) << endl;
  cout << endl;

  frld.verbose=2;
  fld.m=1;
  fld.g=2;
  fld.mu=3;
  fld.mu/=2;
  long double Tld=3;
  Tld/=10;
  frld.calc_mu(fld,Tld);
  cout << "long double:" << endl;
  cout << dtos(fld.n,0) << " " << dtos(fld.ed,0) << endl;
  cout << dtos(fld.pr,0) << " " << dtos(fld.en,0) << endl;
  cout << endl;
       
  fld.m=1;
  fld.g=2;
  fld.mu=4;
  frld.calc_mu(fld,Tld);
  cout << "long double:" << endl;
  cout << dtos(fld.n,0) << " " << dtos(fld.ed,0) << endl;
  cout << dtos(fld.pr,0) << " " << dtos(fld.en,0) << endl;
  cout << endl;

  fr35.verbose=2;
  f35.m=1;
  f35.g=2;
  f35.mu=3;
  f35.mu/=2;
  cpp_dec_float_35 T35=3;
  T35/=10;
  fr35.calc_mu(f35,T35);
  cout << "cpp_dec_float_35:" << endl;
  cout << dtos(f35.n,0) << "\n" << dtos(f35.ed,0) << endl;
  cout << dtos(f35.pr,0) << "\n" << dtos(f35.en,0) << endl;
  cout << endl;

  f35.m=1;
  f35.g=2;
  f35.mu=4;
  fr35.calc_mu(f35,T35);
  cout << "cpp_dec_float_35:" << endl;
  cout << dtos(f35.n,0) << "\n" << dtos(f35.ed,0) << endl;
  cout << dtos(f35.pr,0) << "\n" << dtos(f35.en,0) << endl;
  cout << endl;

  // -------------

  // This doesn't work yet

  if (false) {
    f35.g=2;
    f35.m=1;
    f35.m/=100000;
    f35.mu=10001;
    f35.mu/=100000;
    T35=1;
    T35/=100;
    cout << dtos(f35.m,0) << endl;
    cout << dtos(f35.mu,0) << endl;
    cout << dtos(T35,0) << endl;
    fr35.calc_mu(f35,T35);
    cout << "cpp_dec_float_35:" << endl;
    cout << dtos(f35.n,0) << "\n" << dtos(f35.ed,0) << endl;
    cout << dtos(f35.pr,0) << "\n" << dtos(f35.en,0) << endl;
    cout << endl;
  }

  if (false) {
    Tld=0.01L;
    fld.m=1.0e-5L;
    fld.g=2;
    fld.mu=1.0001e-1L;
    frld2.pair_mu(fld,Tld);
    cout << "long double:" << endl;
    cout << dtos(fld.n,0) << " " << dtos(fld.ed,0) << endl;
    cout << dtos(fld.pr,0) << " " << dtos(fld.en,0) << endl;
    cout << endl;
  }

  // These don't work yet
  
  //part_calibrate_class_tl<long double> pcc_ld;
  //long double vx_ld=pcc_ld.part_calibrate<fermion_ld,fermion_rel_ld>
  //(fld,frld,1,"../../data/o2scl/fermion_deriv_cal.o2",false,2,true);
  //cout << vx_ld << endl;

  //part_calibrate_class_tl<long double> pcc_ld;
  //long double vx_ld2=pcc_ld.part_calibrate<fermion_ld,fermion_rel_ld2>
  //(fld,frld2,1,"../../data/o2scl/fermion_deriv_cal.o2",false,2,true);
  //cout << vx_ld2 << endl;

  //part_calibrate_class_tl<cpp_dec_float_35> pcc_cdf35;
  //cpp_dec_float_35 vx_cdf35=pcc_cdf35.part_calibrate<fermion_cdf35,
  //fermion_rel_cdf35>
  //(f35,fr35,1,"../../data/o2scl/fermion_deriv_cal.o2",false,3,true);
  //cout << vx_cdf35 << endl;
  
#endif
  
  t.report();

  return 0;
}

