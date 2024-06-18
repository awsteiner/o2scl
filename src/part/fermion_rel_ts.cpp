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

template<class fp_t, class fp2_t>
int count_digits_same(fp_t &one, fp2_t &two, std::string msg="") {
  fp_t numer=one-static_cast<fp_t>(two);
  if (numer==0) return std::numeric_limits<fp_t>::max_digits10;
  int ret=((int)(-log10(fabs(numer)/fabs(static_cast<fp_t>(two)))));
  if (ret<-1000) return 99;
  if (ret<=0) {
    cout << "Problem: " << msg << endl;
    cout << dtos(one,-1) << " " << dtos(two,-1) << endl;
    cout << dtos(numer,-1) << " " << ret << endl;
    exit(-1);
  }
  return ret;
}

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(1);

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
  
  if (argc==1) {
    t.report();
    return 0;
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

  int first_test=0;
  
  // An exhaustive comparison of the two algorithms at
  // various levels of precision

  cout.precision(4);
  int count=0;

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

          fr.calc_mu(f,T);
          frld.calc_mu(fld,Tld);
          fr25.calc_mu(f25,T25);
          
          if (true) {

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
          
          int idn, iden, ildn, ilden;
          idn=count_digits_same(f.n,fld.n,"idn+ 1");
          iden=count_digits_same(f.en,fld.en,"iden+ 1");
          ildn=count_digits_same(fld.n,f25.n,"ildn+ 1");
          ilden=count_digits_same(fld.en,f25.en,"ilden+ 1");
          
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
          cout.width(2);
          cout << x << " ";
          long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
          int xld=count_digits_same(fld.pr,pr2ld);
          cout.width(2);
          cout << xld << " ";
          cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
          int x25=count_digits_same(f25.pr,pr225);
          cout.width(2);
          cout << x25 << endl;
          
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
        
          fr.pair_mu(f,T);
          frld.pair_mu(fld,Tld);
          fr25.pair_mu(f25,T25);

          if (true) {

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
          
          int kdn, kden, kldn, klden;
          kdn=count_digits_same(f.n,fld.n,"kdn+ 1");
          kden=count_digits_same(f.en,fld.en,"kden+ 1");
          kldn=count_digits_same(fld.n,f25.n,"kldn+ 1");
          klden=count_digits_same(fld.en,f25.en,"klden+ 1");
        
          cout.width(2);
          cout << kdn << " ";
          cout.width(2);
          cout << kden << " ";
          cout.width(2);
          cout << kldn << " ";
          cout.width(2);
          cout << klden << " ";
        
          double pr2=-f.ed+f.n*f.mu+T*f.en;
          int x=count_digits_same(f.pr,pr2);
          cout.width(2);
          cout << x << " ";
          long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
          int xld=count_digits_same(fld.pr,pr2ld);
          cout.width(2);
          cout << xld << " ";
          cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
          int x25=count_digits_same(f25.pr,pr225);
          cout.width(2);
          cout << x25 << endl;
          
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
        
          fr.calc_mu(f,T);
          frld.calc_mu(fld,Tld);
          fr25.calc_mu(f25,T25);

          if (true) {

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
          
          int idn, iden, ildn, ilden;
          idn=count_digits_same(f.n,fld.n,"idn- 1");
          iden=count_digits_same(f.en,fld.en,"iden- 1");
          ildn=count_digits_same(fld.n,f25.n,"ilnd- 1");
          ilden=count_digits_same(fld.en,f25.en,"ilden- 1");
        
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
          cout.width(2);
          cout << x << " ";
          long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
          int xld=count_digits_same(fld.pr,pr2ld);
          cout.width(2);
          cout << xld << " ";
          cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
          int x25=count_digits_same(f25.pr,pr225);
          cout.width(2);
          cout << x25 << endl;
          
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
        
          fr.pair_mu(f,T);
          frld.pair_mu(fld,Tld);
          fr25.pair_mu(f25,T25);
        
          if (true) {

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
          
          int kdn, kden, kldn, klden;
          kdn=count_digits_same(f.n,fld.n,"kdn- 1");
          kden=count_digits_same(f.en,fld.en,"kden- 1");
          kldn=count_digits_same(fld.n,f25.n,"kldn- 1");
          klden=count_digits_same(fld.en,f25.en,"klden- 1");
        
          cout.width(2);
          cout << kdn << " ";
          cout.width(2);
          cout << kden << " ";
          cout.width(2);
          cout << kldn << " ";
          cout.width(2);
          cout << klden << " ";
        
          double pr2=-f.ed+f.n*f.mu+T*f.en;
          int x=count_digits_same(f.pr,pr2);
          cout.width(2);
          cout << x << " ";
          long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
          int xld=count_digits_same(fld.pr,pr2ld);
          cout.width(2);
          cout << xld << " ";
          cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
          int x25=count_digits_same(f25.pr,pr225);
          cout.width(2);
          cout << x25 << endl;
          
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
          fr.calc_density(f,T);
          fld.mu=fld.m;
          frld.calc_density(fld,Tld);
          f25.mu=f25.m;
          fr25.calc_density(f25,T25);
        
          if (true) {

            // Test with inc_rest_mass=false
            fermion f2=f;
            f2.inc_rest_mass=false;
            fr.calc_density(f2,T);
            t.test_rel(f.mu-f.m,f2.mu,1.0e-10,"irm false mu 5");
            t.test_rel(f.ed,f2.ed+f2.n*f2.m,1.0e-13,"irm false ed 5");
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
            t.test_rel(f.en,f3.en,1.0e-14,"ni false en 5");

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
          
          int idmu, iden, ildmu, ilden;
          idmu=count_digits_same(f.mu,fld.mu,"idmu 1");
          iden=count_digits_same(f.en,fld.en,"iden 1");
          ildmu=count_digits_same(fld.mu,f25.mu,"ildmu 1");
          ilden=count_digits_same(fld.en,f25.en,"ilden 1");
        
          cout.width(2);
          cout << idmu << " ";
          cout.width(2);
          cout << iden << " ";
          cout.width(2);
          cout << ildmu << " ";
          cout.width(2);
          cout << ilden << " ";
        
          double pr2=-f.ed+f.n*f.mu+T*f.en;
          int x=count_digits_same(f.pr,pr2);
          cout.width(2);
          cout << x << " ";
          long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
          int xld=count_digits_same(fld.pr,pr2ld);
          cout.width(2);
          cout << xld << " ";
          cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
          int x25=count_digits_same(f25.pr,pr225);
          cout.width(2);
          cout << x25 << endl;
          
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
          fr.pair_density(f,T);
          fld.mu=fld.m;
          frld.pair_density(fld,Tld);
          f25.mu=f25.m;
          fr25.pair_density(f25,T25);
        
          if (true) {

            // Test with inc_rest_mass=false
            fermion f2=f;
            f2.inc_rest_mass=false;
            fr.pair_density(f2,T);
            t.test_rel(f.mu-f.m,f2.mu,1.0e-10,"irm false mu 6");
            t.test_rel(f.ed,f2.ed+f2.n*f2.m,1.0e-10,"irm false ed 6");
            t.test_rel(f.en,f2.en,1.0e-13,"irm false en 6");

            // Test with non_interacting=false
            fermion f3=f;
            f3.non_interacting=false;
            f3.ms=f3.m;
            f3.nu=f3.mu;
            f3.m*=sqrt(2.0);
            fr.pair_density(f3,T);
            t.test_rel(f.mu,f3.nu,1.0e-11,"ni false mu 6");
            t.test_rel(f.n,f3.n,1.0e-14,"ni false n 6");
            t.test_rel(f.en,f3.en,1.0e-14,"ni false en 6");

            // Test with both equal to false
            fermion f4=f3;
            f4.inc_rest_mass=false;
            f4.non_interacting=false;
            fr.pair_density(f4,T);
            t.test_rel(f3.mu-f3.m,f4.nu,1.0e-13,"both false mu 6");
            t.test_rel(f3.ed,f4.ed+f4.n*f4.m,1.0e-9,"both false ed 6");
            t.test_rel(f3.en,f4.en,1.0e-13,"both false en 6");

            if (t.get_success()==false) {
              cout << "6" << endl;
              exit(-1);
            }
          }
          
          int kdmu, kden, kldmu, klden;
          kdmu=count_digits_same(f.mu,fld.mu,"kdmu 1");
          kden=count_digits_same(f.en,fld.en,"kden 1");
          kldmu=count_digits_same(fld.mu,f25.mu,"kldmu 1");
          klden=count_digits_same(fld.en,f25.en,"klden 1");
        
          cout.width(2);
          cout << kdmu << " ";
          cout.width(2);
          cout << kden << " ";
          cout.width(2);
          cout << kldmu << " ";
          cout.width(2);
          cout << klden << " ";
        
          double pr2=-f.ed+f.n*f.mu+T*f.en;
          int x=count_digits_same(f.pr,pr2);
          cout.width(2);
          cout << x << " ";
          long double pr2ld=-fld.ed+fld.n*fld.mu+T*fld.en;
          int xld=count_digits_same(fld.pr,pr2ld);
          cout.width(2);
          cout << xld << " ";
          cpp_dec_float_25 pr225=-f25.ed+f25.n*f25.mu+T*f25.en;
          int x25=count_digits_same(f25.pr,pr225);
          cout.width(2);
          cout << x25 << endl;
          
        } else {

          cout << endl;

        }
          
        count++;
        
      }
    }
  }

#endif
  
  t.report();

  return 0;
}

