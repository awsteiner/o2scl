/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
#include <o2scl/fermion_rel.h>
#include <o2scl/test_mgr.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/eos_sn.h>

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
  
  std::string arg;
  if (argc>=2) {
    arg=argv[1];
  }
  
  fermion f;
  fermion_ld fld;
  fermion_cdf25 f25;
  f.g=2;
  fld.g=2;
  f25.g=2;
  
  fermion_rel fr;
  fermion_rel_ld frld;
  fermion_rel_cdf25 fr25;
  
  fr.verify_ti=true;
  frld.verify_ti=true;
  fr25.verify_ti=true;

  fr.err_nonconv=false;
  fr.nit.err_nonconv=false;
  fr.dit.err_nonconv=false;
  fr.it_multip.err_nonconv=false;

  frld.err_nonconv=false;
  frld.nit.err_nonconv=false;
  frld.dit.err_nonconv=false;
  frld.it_multip.err_nonconv=false;
  
  fr25.err_nonconv=false;
  fr25.nit.err_nonconv=false;
  fr25.dit.err_nonconv=false;
  fr25.it_multip.err_nonconv=false;

  int first_test=0;
  double test_shift=1.0;
  
  // An exhaustive comparison of the fermion_rel class at various
  // levels of precision

  if (arg=="1") {

    // This tests higher-accuracy settings. Smaller tolerances for the
    // integrals cause some of the integrations to fail. Note that a
    // larger value of 'test_shift' is also required, but it's not
    // clear to me if this really represents an accuracy issue.
    
    fr.dit.tol_abs=1.0e-11;
    fr.dit.tol_rel=1.0e-11;
    fr.nit.tol_abs=1.0e-11;
    fr.nit.tol_rel=1.0e-11;
    fr.upper_limit_fac=40.0;
    fr.density_root.tol_rel=1.0e-10;
    fr.def_massless_root.tol_rel=1.0e-10;
    
    frld.dit.tol_abs=1.0e-15;
    frld.dit.tol_rel=1.0e-15;
    frld.nit.tol_abs=1.0e-15;
    frld.nit.tol_rel=1.0e-15;
    frld.upper_limit_fac=60.0;
    frld.density_root.tol_rel=1.0e-14;
    frld.def_massless_root.tol_rel=1.0e-14;
    
    fr25.dit.tol_abs=1.0e-18;
    fr25.dit.tol_rel=1.0e-18;
    fr25.nit.tol_abs=1.0e-18;
    fr25.nit.tol_rel=1.0e-18;
    fr25.upper_limit_fac=80.0;
    fr25.density_root.tol_rel=1.0e-18;
    fr25.def_massless_root.tol_rel=1.0e-18;
    
    test_shift=10.0;
  }
  
  if (arg=="2") {

    // I think this runs without throwing any exceptions, but it needs
    // some work, especially to ensure the thermodynamic identity is
    // satisfied at higher precision. Either way, it is quite slow.
    
    fr.multip=true;
    frld.multip=true;
    fr25.multip=true;
    
    fr.upper_limit_fac=40.0;
    fr.density_root.tol_rel=1.0e-10;
    fr.def_massless_root.tol_rel=1.0e-10;
    
    frld.upper_limit_fac=52.0;
    frld.deg_entropy_fac=52.0;
    frld.tol_expan=1.0e-18;
    frld.exp_limit=11400.0;
    frld.density_root.tol_rel=1.0e-18;
    frld.def_massless_root.tol_rel=1.0e-18;
    
    fr25.upper_limit_fac=62.0;
    fr25.deg_entropy_fac=62.0;
    fr25.tol_expan=1.0e-23;
    fr25.exp_limit=6.7e7;
    fr25.density_root.tol_rel=1.0e-23;
    fr25.def_massless_root.tol_rel=1.0e-23;
  }

  /*
  if (arg=="3") {
    
    fr.dit.tol_abs=1.0e-13;
    fr.dit.tol_rel=1.0e-13;
    fr.nit.tol_abs=1.0e-13;
    fr.nit.tol_rel=1.0e-13;
    fr.upper_limit_fac=40.0;
    fr.density_root.tol_rel=1.0e-10;
    
    frld.dit.tol_abs=1.0e-18;
    frld.dit.tol_rel=1.0e-18;
    frld.nit.tol_abs=1.0e-18;
    frld.nit.tol_rel=1.0e-18;
    frld.density_root.tol_rel=1.0e-18;
    frld.def_massless_root.tol_rel=1.0e-18;
    frld.upper_limit_fac=52.0;
    frld.deg_entropy_fac=52.0;
    frld.tol_expan=1.0e-18;
    frld.exp_limit=11400.0;
    
    fr25.dit.tol_abs=1.0e-25;
    fr25.dit.tol_rel=1.0e-25;
    fr25.nit.tol_abs=1.0e-25;
    fr25.nit.tol_rel=1.0e-25;
    fr25.density_root.tol_rel=1.0e-23;
    fr25.def_massless_root.tol_rel=1.0e-23;
    fr25.upper_limit_fac=62.0;
    fr25.deg_entropy_fac=62.0;
    fr25.tol_expan=1.0e-23;
    fr25.exp_limit=6.7e7;
    
  }
  */
  
  cout.precision(4);
  int count=0;
  // Sums of calc_mu() comparisons between fp types
  int cmu_n=0, cmu_en=0, cmu_ld_n=0, cmu_ld_en=0;
  // Sums of calc_mu() accuracy via. thermodynamic identity
  int cmu_ti=0, cmu_ld_ti=0, cmu_25_ti=0;
  // Sums of pair_mu() comparisons between fp types
  int pmu_n=0, pmu_en=0, pmu_ld_n=0, pmu_ld_en=0;
  // Sums of pair_mu() accuracy via. thermodynamic identity
  int pmu_ti=0, pmu_ld_ti=0, pmu_25_ti=0;
  // Sums of calc_density() comparisons between fp types
  int cd_mu=0, cd_en=0, cd_ld_mu=0, cd_ld_en=0;
  // Sums of calc_density() accuracy via. thermodynamic identity
  int cd_ti=0, cd_ld_ti=0, cd_25_ti=0;
  // Sums of pair_density() comparisons between fp types
  int pd_mu=0, pd_en=0, pd_ld_mu=0, pd_ld_en=0;
  // Sums of pair_density() accuracy via. thermodynamic identity
  int pd_ti=0, pd_ld_ti=0, pd_25_ti=0;

  part_cal_new<> pcn;
  if (false) {
    pcn.test_calc_mu(f,fld,f25,fr,frld,fr25,t,count,first_test,
                     1564,1544,2393,2176,1432,1906,2597);
  }
  //pcn.test_pair_mu(f,fld,f25,fr,frld,fr25,t,count,first_test,
  //1618,1604,2336,2218,1408,1642,1947);
  //exit(-1);

  cout << " cnt m          T           mu/n       "
       << "d-ld  ld-25 ti verify" << endl;
  for(int im=-2;im<=1;im++) {
    
    f.m=im;
    f.m=pow(10,f.m);
    fld.m=im;
    fld.m=pow(10,fld.m);
    f25.m=im;
    f25.m=pow(10,f25.m);
    
    for(int iT=-2;iT<=1;iT++) {
      
      double T=iT;
      T=pow(10,T);
      long double Tld=iT;
      Tld=pow(10,Tld);
      long double T25=iT;
      T25=pow(10,T25);
      
      for(int imu=-2;imu<=1;imu++) {
        
        f.mu=imu;
        f.mu=pow(10,f.mu);
        fld.mu=imu;
        fld.mu=pow(10,fld.mu);
        f25.mu=imu;
        f25.mu=pow(10,f25.mu);

        cout.width(4);
        cout << count << " ";
        
        if (count>=first_test) {
          
          cout << f.m << " " << T << " ";
          cout.setf(ios::showpos);
          cout << f.mu << " ";
          cout.unsetf(ios::showpos);

          int ret=fr.calc_mu(f,T);
          int retld=frld.calc_mu(fld,Tld);
          int ret25=fr25.calc_mu(f25,T25);
          
          if (ret==0) {

            // Test with inc_rest_mass=false
            fermion f2=f;
            f2.inc_rest_mass=false;
            f2.mu-=f2.m;
            fr.calc_mu(f2,T);
            t.test_rel(f.ed,f2.ed+f2.m*f2.n,1.0e-11,"irm false ed");
            t.test_rel(f.en,f2.en,1.0e-14,"irm false en");

            // Test with non_interacting=false
            fermion f3=f;
            f3.non_interacting=false;
            f3.nu=f3.mu;
            f3.ms=f3.m;
            f3.m*=sqrt(2.0);
            fr.calc_mu(f3,T);
            t.test_rel(f.n,f3.n,1.0e-14,"ni false n");
            t.test_rel(f.en,f3.en,1.0e-14,"ni false en");

            // Test with both equal to false
            fermion f4=f3;
            f4.inc_rest_mass=false;
            f4.non_interacting=false;
            f4.nu-=f4.m;
            fr.calc_mu(f4,T);
            t.test_rel(f3.ed,f4.ed+f4.n*f4.m,1.0e-11,"both false ed");
            t.test_rel(f3.en,f4.en,1.0e-13,"both false en");

            if (t.get_success()==false) {
              cout << "1" << endl;
              exit(-1);
            }
          }
          
          int idn=-2, iden=-2, ildn=-2, ilden=-2;
          if (ret==0 && retld==0) {
            idn=count_digits_same(f.n,fld.n);
            iden=count_digits_same(f.en,fld.en);
          }
          cmu_n+=idn;
          cmu_en+=iden;
          if (retld==0 && ret25==0) {
            ildn=count_digits_same(fld.n,f25.n);
            ilden=count_digits_same(fld.en,f25.en);
          }
          cmu_ld_en+=ilden;
          cmu_ld_n+=ildn;
          
          cout.width(2);
          cout << idn << " ";
          cout.width(2);
          cout << iden << " ";
          cout.width(2);
          cout << ildn << " ";
          cout.width(2);
          cout << ilden << " ";

          int x=-2, xld=-2, x25=-2;
          
          if (ret==0) {
            double pr2=-f.ed+f.n*f.mu+T*f.en;
            x=count_digits_same(f.pr,pr2);
          }
          cmu_ti+=x;
          cout.width(2);
          cout << x << " ";
          
          if (retld==0) {
            long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
            xld=count_digits_same(fld.pr,pr2ld);
          }
          cmu_ld_ti+=xld;
          cout.width(2);
          cout << xld << " ";

          if (ret25==0) {
            cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
            x25=count_digits_same(f25.pr,pr225);
          }
          cmu_25_ti+=x25;
          cout.width(2);
          cout << x25 << " cmu" << endl;

        } else {
          
          cout << endl;
          
        }

        count++;
        
        cout.width(4);
        cout << count << " ";
        
        if (count>=first_test) {
          
          cout << f.m << " " << T << " ";
          cout.setf(ios::showpos);
          cout << f.mu << " ";
          cout.unsetf(ios::showpos);
        
          int ret=fr.pair_mu(f,T);
          int retld=frld.pair_mu(fld,Tld);
          int ret25=fr25.pair_mu(f25,T25);

          if (ret==0) {

            // Test with inc_rest_mass=false
            fermion f2=f;
            f2.inc_rest_mass=false;
            f2.mu-=f2.m;
            fr.pair_mu(f2,T);
            t.test_rel(f.ed,f2.ed+f2.n*f2.m,1.0e-11,"irm false ed 2");
            t.test_rel(f.en,f2.en,1.0e-14,"irm false en 2");

            // Test with non_interacting=false
            fermion f3=f;
            f3.non_interacting=false;
            f3.nu=f3.mu;
            f3.ms=f3.m;
            f3.m*=sqrt(2.0);
            fr.pair_mu(f3,T);
            t.test_rel(f.n,f3.n,1.0e-14,"ni false n 2");
            t.test_rel(f.en,f3.en,1.0e-14,"ni false en 2");

            // Test with both equal to false
            fermion f4=f3;
            f4.inc_rest_mass=false;
            f4.non_interacting=false;
            f4.nu-=f4.m;
            fr.pair_mu(f4,T);
            t.test_rel(f3.ed,f4.ed+f4.n*f4.m,1.0e-11,"both false ed 2");
            t.test_rel(f3.en,f4.en,1.0e-13,"both false en 2");
            
            if (t.get_success()==false) {
              cout << "2" << endl;
              exit(-1);
            }
          }
          
          int kdn=-2, kden=-2, kldn=-2, klden=-2;
          if (ret==0 && retld==0) {
            kdn=count_digits_same(f.n,fld.n);
            kden=count_digits_same(f.en,fld.en);
          }
          pmu_n+=kdn;
          pmu_en+=kden;
          if (retld==0 && ret25==0) {
            kldn=count_digits_same(fld.n,f25.n);
            klden=count_digits_same(fld.en,f25.en);
          }
          pmu_ld_n+=kldn;
          pmu_ld_en+=klden;
        
          cout.width(2);
          cout << kdn << " ";
          cout.width(2);
          cout << kden << " ";
          cout.width(2);
          cout << kldn << " ";
          cout.width(2);
          cout << klden << " ";

          int x=-2, xld=-2, x25=-2;
          if (ret==0) {
            double pr2=-f.ed+f.n*f.mu+T*f.en;
            x=count_digits_same(f.pr,pr2);
          }
          pmu_ti+=x;
          cout.width(2);
          cout << x << " ";

          if (retld==0) {
            long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
            xld=count_digits_same(fld.pr,pr2ld);
          }
          pmu_ld_ti+=x;
          cout.width(2);
          cout << xld << " ";

          if (ret25==0) {
            cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
            x25=count_digits_same(f25.pr,pr225);
          }
          pmu_25_ti+=x;
          cout.width(2);
          cout << x25 << " pmu" << endl;
          
        } else {
          cout << endl;
        }
        
        count++;
        
      }
      
      for(int imu=-2;imu<=1;imu++) {
        
        f.mu=imu;
        f.mu=-pow(10,f.mu);
        fld.mu=imu;
        fld.mu=-pow(10,fld.mu);
        f25.mu=imu;
        f25.mu=-pow(10,f25.mu);
        
        cout.width(4);
        cout << count << " ";

        if (count>=first_test) {
        
          cout << f.m << " " << T << " ";
          cout.setf(ios::showpos);
          cout << f.mu << " ";
          cout.unsetf(ios::showpos);
        
          int ret=fr.calc_mu(f,T);
          int retld=frld.calc_mu(fld,Tld);
          int ret25=fr25.calc_mu(f25,T25);

          if (ret==0) {

            // Test with inc_rest_mass=false
            fermion f2=f;
            f2.inc_rest_mass=false;
            f2.mu-=f2.m;
            fr.calc_mu(f2,T);
            t.test_rel(f.ed,f2.ed+f2.n*f2.m,1.0e-11,"irm false ed 3");
            t.test_rel(f.en,f2.en,1.0e-14,"irm false en 3");

            // Test with non_interacting=false
            fermion f3=f;
            f3.non_interacting=false;
            f3.nu=f3.mu;
            f3.ms=f3.m;
            f3.m*=sqrt(2.0);
            fr.calc_mu(f3,T);
            t.test_rel(f.n,f3.n,1.0e-14,"ni false n 3");
            t.test_rel(f.en,f3.en,1.0e-14,"ni false en 3");

            // Test with both equal to false
            fermion f4=f3;
            f4.inc_rest_mass=false;
            f4.non_interacting=false;
            f4.nu-=f4.m;
            fr.calc_mu(f4,T);
            t.test_rel(f3.ed,f4.ed+f4.n*f4.m,1.0e-11,"both false ed 3");
            t.test_rel(f3.en,f4.en,1.0e-13,"both false en 3");
            
            if (t.get_success()==false) {
              cout << "3" << endl;
              exit(-1);
            }
          }
          
          int idn=-2, iden=-2, ildn=-2, ilden=-2;
          if (ret==0 && retld==0) {
            idn=count_digits_same(f.n,fld.n);
            iden=count_digits_same(f.en,fld.en);
          }
          if (retld==0 && ret25==0) {
            ildn=count_digits_same(fld.n,f25.n);
            ilden=count_digits_same(fld.en,f25.en);
          }
          cmu_n+=idn;
          cmu_en+=iden;
          cmu_ld_n+=ildn;
          cmu_ld_en+=ilden;
        
          cout.width(2);
          cout << idn << " ";
          cout.width(2);
          cout << iden << " ";
          cout.width(2);
          cout << ildn << " ";
          cout.width(2);
          cout << ilden << " ";
        
          double pr2=-f.ed+f.n*f.mu+T*f.en;
          int x=count_digits_same(f.pr,pr2);
          cmu_ti+=x;
          cout.width(2);
          cout << x << " ";
          long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
          int xld=count_digits_same(fld.pr,pr2ld);
          cmu_ld_ti+=xld;
          cout.width(2);
          cout << xld << " ";
          cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
          int x25=count_digits_same(f25.pr,pr225);
          cmu_25_ti+=x25;
          cout.width(2);
          cout << x25 << " cmu" << endl;
          
        } else {
          cout << endl;
        }
        
        count++;
        
        cout.width(4);
        cout << count << " ";

        if (count>=first_test) {
          
          cout << f.m << " " << T << " ";
          cout.setf(ios::showpos);
          cout << f.mu << " ";
          cout.unsetf(ios::showpos);
        
          int ret=fr.pair_mu(f,T);
          int retld=frld.pair_mu(fld,Tld);
          int ret25=fr25.pair_mu(f25,T25);
        
          if (ret==0) {

            // Test with inc_rest_mass=false
            fermion f2=f;
            f2.inc_rest_mass=false;
            f2.mu-=f2.m;
            fr.pair_mu(f2,T);
            t.test_rel(f.ed,f2.ed+f2.n*f2.m,1.0e-11,"irm false ed 4");
            t.test_rel(f.en,f2.en,1.0e-13,"irm false en 4");

            // Test with non_interacting=false
            fermion f3=f;
            f3.non_interacting=false;
            f3.nu=f3.mu;
            f3.ms=f3.m;
            f3.m*=sqrt(2.0);
            fr.pair_mu(f3,T);
            t.test_rel(f.n,f3.n,1.0e-14,"ni false n 4");
            t.test_rel(f.en,f3.en,1.0e-14,"ni false en 4");

            // Test with both equal to false
            fermion f4=f3;
            f4.inc_rest_mass=false;
            f4.non_interacting=false;
            f4.nu-=f4.m;
            fr.pair_mu(f4,T);
            t.test_rel(f3.ed,f4.ed+f4.n*f4.m,1.0e-11,"both false ed 4");
            t.test_rel(f3.en,f4.en,1.0e-13,"both false en 4");
            
            if (t.get_success()==false) {
              cout << "4" << endl;
              exit(-1);
            }
          }
          
          int kdn=-2, kden=-2, kldn=-2, klden=-2;
          if (ret==0 && retld==0) {
            kdn=count_digits_same(f.n,fld.n);
            kden=count_digits_same(f.en,fld.en);
          }
          pmu_n+=kdn;
          pmu_en+=kden;
          if (retld==0 && ret25==0) {
            kldn=count_digits_same(fld.n,f25.n);
            klden=count_digits_same(fld.en,f25.en);
          }
          pmu_ld_n+=kldn;
          pmu_ld_en+=klden;
        
          cout.width(2);
          cout << kdn << " ";
          cout.width(2);
          cout << kden << " ";
          cout.width(2);
          cout << kldn << " ";
          cout.width(2);
          cout << klden << " ";

          int x=-2, xld=-2, x25=-2;
          if (ret==0) {
            double pr2=-f.ed+f.n*f.mu+T*f.en;
            x=count_digits_same(f.pr,pr2);
          }
          pmu_ti+=x;
          cout.width(2);
          cout << x << " ";
          if (retld==0) {
            long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
            xld=count_digits_same(fld.pr,pr2ld);
          }
          pmu_ld_ti+=xld;
          cout.width(2);
          cout << xld << " ";
          if (ret25==0) {
            cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
            x25=count_digits_same(f25.pr,pr225);
          }
          pmu_25_ti+=x25;
          cout.width(2);
          cout << x25 << " pmu" << endl;
          
        } else {
          cout << endl;
        }
        
        count++;
        
      }
    }
  }

  for(int im=-2;im<=1;im++) {
    
    f.m=im;
    f.m=pow(10,f.m);
    fld.m=im;
    fld.m=pow(10,fld.m);
    f25.m=im;
    f25.m=pow(10,f25.m);
    
    for(int iT=-2;iT<=1;iT++) {
      
      double T=iT;
      T=pow(10,T);
      long double Tld=iT;
      Tld=pow(10,Tld);
      long double T25=iT;
      T25=pow(10,T25);
      
      for(int in=-2;in<=1;in++) {
        
        f.n=in;
        f.n=pow(10,f.n);
        fld.n=in;
        fld.n=pow(10,fld.n);
        f25.n=in;
        f25.n=pow(10,f25.n);
        
        cout.width(4);
        cout << count << " ";

        if (count>=first_test) {
        
          cout << f.m << " " << T << " ";
          cout.setf(ios::showpos);
          cout << f.n << " ";
          cout.unsetf(ios::showpos);
        
          if (false && count==322) {
            fr.verbose=2;
          }
        
          f.mu=f.m;
          int ret=fr.calc_density(f,T);
          fld.mu=fld.m;
          int retld=frld.calc_density(fld,Tld);
          f25.mu=f25.m;
          int ret25=fr25.calc_density(f25,T25);
        
          if (ret==0) {

            // Test with inc_rest_mass=false
            fermion f2=f;
            f2.inc_rest_mass=false;
            fr.calc_density(f2,T);
            t.test_rel(f.mu-f.m,f2.mu,1.0e-10,"irm false mu 5");
            t.test_rel(f.ed,f2.ed+f2.n*f2.m,1.0e-12,"irm false ed 5");
            t.test_rel(f.en,f2.en,1.0e-13,"irm false en 5");

            // Test with non_interacting=false
            fermion f3=f;
            f3.non_interacting=false;
            f3.ms=f3.m;
            f3.nu=f3.mu;
            f3.m*=sqrt(2.0);
            fr.calc_density(f3,T);
            t.test_rel(f.mu,f3.nu,1.0e-14,"ni false mu 5");
            t.test_rel(f.n,f3.n,1.0e-14,"ni false n 5");
            t.test_rel(f.en,f3.en,1.0e-14*test_shift,"ni false en 5");

            // Test with both equal to false
            fermion f4=f3;
            f4.inc_rest_mass=false;
            f4.non_interacting=false;
            fr.calc_density(f4,T);
            t.test_rel(f3.mu-f3.m,f4.nu,1.0e-14,"both false mu 5");
            t.test_rel(f3.ed,f4.ed+f4.n*f4.m,1.0e-12,"both false ed 5");
            t.test_rel(f3.en,f4.en,1.0e-13,"both false en 5");

            if (t.get_success()==false) {
              cout << "5" << endl;
              exit(-1);
            }
          }
          
          int idmu=-2, iden=-2, ildmu=-2, ilden=-2;
          if (ret==0 && retld==0) {
            idmu=count_digits_same(f.mu,fld.mu);
            iden=count_digits_same(f.en,fld.en);
          }
          cd_mu+=idmu;
          cd_en+=iden;
          if (retld==0 && ret25==0) {
            ildmu=count_digits_same(fld.mu,f25.mu);
            ilden=count_digits_same(fld.en,f25.en);
          }
          cd_ld_mu+=ildmu;
          cd_ld_en+=ilden;
        
          cout.width(2);
          cout << idmu << " ";
          cout.width(2);
          cout << iden << " ";
          cout.width(2);
          cout << ildmu << " ";
          cout.width(2);
          cout << ilden << " ";

          int x=-2, xld=-2, x25=-2;
          if (ret==0) {
            double pr2=-f.ed+f.n*f.mu+T*f.en;
            x=count_digits_same(f.pr,pr2);
          }
          cd_ti+=x;
          cout.width(2);
          cout << x << " ";
          if (retld==0) {
            long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
            xld=count_digits_same(fld.pr,pr2ld);
          }
          cd_ld_ti+=x;
          cout.width(2);
          cout << xld << " ";
          if (ret25==0) {
            cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
            x25=count_digits_same(f25.pr,pr225);
          }
          cd_25_ti+=x;
          cout.width(2);
          cout << x25 << " cde" << endl;
          
        } else {

          cout << endl;
        }
        
        count++;
        
        cout.width(4);
        cout << count << " ";

        if (count>=first_test) {
            
          cout << f.m << " " << T << " ";
          cout.setf(ios::showpos);
          cout << f.n << " ";
          cout.unsetf(ios::showpos);
        
          if (false && count==259) {
            fr.verbose=2;
          }
        
          f.mu=f.m;
          int ret=fr.pair_density(f,T);
          fld.mu=fld.m;
          int retld=frld.pair_density(fld,Tld);
          f25.mu=f25.m;
          int ret25=fr25.pair_density(f25,T25);
        
          if (ret==0) {

            // Test with inc_rest_mass=false
            fermion f2=f;
            f2.inc_rest_mass=false;
            fr.pair_density(f2,T);
            t.test_rel(f.mu-f.m,f2.mu,1.0e-10,"irm false mu 6");
            t.test_rel(f.ed,f2.ed+f2.n*f2.m,1.0e-10*test_shift,
                       "irm false ed 6");
            t.test_rel(f.en,f2.en,1.0e-12,"irm false en 6");

            // Test with non_interacting=false
            fermion f3=f;
            f3.non_interacting=false;
            f3.ms=f3.m;
            f3.nu=f3.mu;
            f3.m*=sqrt(2.0);
            fr.pair_density(f3,T);
            t.test_rel(f.mu,f3.nu,1.0e-10*test_shift,"ni false mu 6");
            t.test_rel(f.n,f3.n,1.0e-14,"ni false n 6");
            t.test_rel(f.en,f3.en,1.0e-14*test_shift,"ni false en 6");

            // Test with both equal to false
            fermion f4=f3;
            f4.inc_rest_mass=false;
            f4.non_interacting=false;
            fr.pair_density(f4,T);
            t.test_rel(f3.mu-f3.m,f4.nu,1.0e-12,"both false mu 6");
            t.test_rel(f3.ed,f4.ed+f4.n*f4.m,1.0e-9,"both false ed 6");
            t.test_rel(f3.en,f4.en,1.0e-13,"both false en 6");

            if (t.get_success()==false) {
              cout << "6" << endl;
              exit(-1);
            }
          }
          
          int kdmu=-2, kden=-2, kldmu=-2, klden=-2;
          if (ret==0 && retld==0) {
            kdmu=count_digits_same(f.mu,fld.mu);
            kden=count_digits_same(f.en,fld.en);
          }
          pd_mu+=kdmu;
          pd_en+=kden;
          if (retld==0 && ret25==0) {
            kldmu=count_digits_same(fld.mu,f25.mu);
            klden=count_digits_same(fld.en,f25.en);
          }
          pd_ld_mu+=kldmu;
          pd_ld_en+=klden;
        
          cout.width(2);
          cout << kdmu << " ";
          cout.width(2);
          cout << kden << " ";
          cout.width(2);
          cout << kldmu << " ";
          cout.width(2);
          cout << klden << " ";

          int x=-2, xld=-2, x25=-2;
          if (ret==0) {
            double pr2=-f.ed+f.n*f.mu+T*f.en;
            x=count_digits_same(f.pr,pr2);
          }
          pd_ti+=x;
          cout.width(2);
          cout << x << " ";
          if (retld==0) {
            long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
            xld=count_digits_same(fld.pr,pr2ld);
          }
          pd_ld_ti+=xld;
          cout.width(2);
          cout << xld << " ";
          if (ret25==0) {
            cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
            x25=count_digits_same(f25.pr,pr225);
          }
          pd_25_ti+=x25;
          cout.width(2);
          cout << x25 << " pde" << endl;
          
        } else {

          cout << endl;

        }
          
        count++;
        
      }
    }
  }

#endif

  cout << "calc_mu density (double <-> long double): "
       << cmu_n << endl;
  cout << "calc_mu entropy (double <-> long double): "
       << cmu_en << endl;
  cout << "calc_mu density (long double <-> cdf_25): "
       << cmu_ld_n << endl;
  cout << "calc_mu entropy (long double <-> cdf_25): "
       << cmu_ld_en << endl;
  cout << "calc_mu ti: " << cmu_ti << endl;
  cout << "calc_mu long double ti: " << cmu_ld_ti << endl;
  cout << "calc_mu cpp_dec_float_25 ti: " << cmu_25_ti << endl;
  cout << endl;
  
  cout << "pair_mu density (double <-> long double): "
       << pmu_n << endl;
  cout << "pair_mu entropy (double <-> long double): "
       << pmu_en << endl;
  cout << "pair_mu density (long double <-> cdf_25): "
       << pmu_ld_n << endl;
  cout << "pair_mu entropy (long double <-> cdf_25): "
       << pmu_ld_en << endl;
  cout << "pair_mu ti: " << pmu_ti << endl;
  cout << "pair_mu long double ti: " << pmu_ld_ti << endl;
  cout << "pair_mu cpp_dec_float_25 ti: " << pmu_25_ti << endl;
  cout << endl;
  
  cout << "calc_density density (double <-> long double): "
       << cd_mu << endl;
  cout << "calc_density entropy (double <-> long double): "
       << cd_en << endl;
  cout << "calc_density density (long double <-> cdf_25): "
       << cd_ld_mu << endl;
  cout << "calc_density entropy (long double <-> cdf_25): "
       << cd_ld_en << endl;
  cout << "calc_density ti: " << cd_ti << endl;
  cout << "calc_density long double ti: " << cd_ld_ti << endl;
  cout << "calc_density cpp_dec_float_25 ti: " << cd_25_ti << endl;
  cout << endl;
  
  cout << "pair_density density (double <-> long double): "
       << pd_mu << endl;
  cout << "pair_density entropy (double <-> long double): "
       << pd_en << endl;
  cout << "pair_density density (long double <-> cdf_25): "
       << pd_ld_mu << endl;
  cout << "pair_density entropy (long double <-> cdf_25): "
       << pd_ld_en << endl;
  cout << "pair_density ti: " << pd_ti << endl;
  cout << "pair_density long double ti: " << pd_ld_ti << endl;
  cout << "pair_density cpp_dec_float_25 ti: " << pd_25_ti << endl;
  cout << endl;

  if (argc<2) {

    // Normal case, used in 'make o2scl-test'
    t.test_gen(cmu_n>=1543,"cmu_n");
    t.test_gen(cmu_en>=1544,"cmu_en");
    t.test_gen(cmu_ld_n>=2407,"cmu_ld_n");
    t.test_gen(cmu_ld_en>=2179,"cmu_ld_en");
    t.test_gen(cmu_ti>=1400,"cmu_ti");
    t.test_gen(cmu_ld_ti>=1906,"cmu_ld_ti");
    t.test_gen(cmu_25_ti>=2597,"cmu_25_ti");
    
    t.test_gen(pmu_n>=1582,"pmu_n");
    t.test_gen(pmu_en>=1604,"pmu_en");
    t.test_gen(pmu_ld_n>=2336,"pmu_ld_n");
    t.test_gen(pmu_ld_en>=2218,"pmu_ld_en");
    t.test_gen(pmu_ti>=1376,"pmu_ti");
    t.test_gen(pmu_ld_ti>=1626,"pmu_ld_ti");
    t.test_gen(pmu_25_ti>=1931,"pmu_25_ti");
    
    t.test_gen(cd_mu>=891,"cd_mu");
    t.test_gen(cd_en>=837,"cd_en");
    t.test_gen(cd_ld_mu>=1253,"cd_ld_mu");
    t.test_gen(cd_ld_en>=1117,"cd_ld_en");
    t.test_gen(cd_ti>=967,"cd_ti");
    t.test_gen(cd_ld_ti>=967,"cd_ld_ti");
    t.test_gen(cd_25_ti>=967,"cd_25_ti");
    
    t.test_gen(pd_mu>=850,"pd_mu");
    t.test_gen(pd_en>=814,"pd_en");
    t.test_gen(pd_ld_mu>=1190,"pd_ld_mu");
    t.test_gen(pd_ld_en>=1101,"pd_ld_en");
    t.test_gen(pd_ti>=676,"pd_ti");
    t.test_gen(pd_ld_ti>=817,"pd_ld_ti");
    t.test_gen(pd_25_ti>=1025,"pd_25_ti");
    
  } else if (arg=="1") {
    
    // Improved accuracy case
    t.test_gen(cmu_n>=1738,"cmu_n");
    t.test_gen(cmu_en>=1612,"cmu_en");
    t.test_gen(cmu_ld_n>=2379,"cmu_ld_n");
    t.test_gen(cmu_ld_en>=2174,"cmu_ld_en");
    t.test_gen(cmu_ti>=1540,"cmu_ti");
    t.test_gen(cmu_ld_ti>=1992,"cmu_ld_ti");
    t.test_gen(cmu_25_ti>=2702,"cmu_25_ti");
    
    t.test_gen(pmu_n>=1756,"pmu_n");
    t.test_gen(pmu_en>=1670,"pmu_en");
    t.test_gen(pmu_ld_n>=2298,"pmu_ld_n");
    t.test_gen(pmu_ld_en>=2222,"pmu_ld_en");
    t.test_gen(pmu_ti>=1592,"pmu_ti");
    t.test_gen(pmu_ld_ti>=1821,"pmu_ld_ti");
    t.test_gen(pmu_25_ti>=2144,"pmu_25_ti");
    
    t.test_gen(cd_mu>=991,"cd_mu");
    t.test_gen(cd_en>=901,"cd_en");
    t.test_gen(cd_ld_mu>=1237,"cd_ld_mu");
    t.test_gen(cd_ld_en>=1104,"cd_ld_en");
    t.test_gen(cd_ti>=969,"cd_ti");
    t.test_gen(cd_ld_ti>=969,"cd_ld_ti");
    t.test_gen(cd_25_ti>=969,"cd_25_ti");
    
    t.test_gen(pd_mu>=939,"pd_mu");
    t.test_gen(pd_en>=888,"pd_en");
    t.test_gen(pd_ld_mu>=1190,"pd_ld_mu");
    t.test_gen(pd_ld_en>=1094,"pd_ld_en");
    t.test_gen(pd_ti>=819,"pd_ti");
    t.test_gen(pd_ld_ti>=997,"pd_ld_ti");
    t.test_gen(pd_25_ti>=1271,"pd_25_ti");
    
  }
  
  t.report();

  return 0;
}

