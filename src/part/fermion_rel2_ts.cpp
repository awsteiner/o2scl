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
#include <o2scl/fermion_rel2.h>
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

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  bool calibrate_ld=false;
  bool calibrate_cdf25=false;
  bool calibrate_new=false;
  if (argc>=2) {
    if (std::stoi(argv[1])==1) calibrate_ld=true;
    else if (std::stoi(argv[1])==2) calibrate_cdf25=true;
    else if (std::stoi(argv[1])==3) calibrate_new=true;
  }

  test_mgr t;
  t.set_output_level(2);
  part_calibrate_class pcc;
  part_calibrate_class_tl<long double> pcc_ld;
  part_calibrate_class_tl<cpp_dec_float_25> pcc_cdf25;

  fermion f(1.0,2.0);
  fermion f2(1.0,2.0);
  fermion_rel2 fr;

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
  double v1=pcc.part_calibrate<fermion,fermion_rel2>
    (f,fr,"../../data/o2scl/fermion_deriv_cal.o2",true,true,false,0,true);
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

  // Check calc_mu() alone
  double v2=pcc.part_calibrate<fermion,fermion_rel2>
    (f,fr,"../../data/o2scl/fermion_deriv_cal.o2",true,true,false,0,true);
  t.test_rel(v2,0.0,4.0e-10,"calibrate 2");
  cout << endl;
  
#ifndef O2SCL_NO_BOOST_MULTIPRECISION

  fermion_ld fld;
  fermion_rel2_ld frld;
  //fermion_rel2_ld_multip frld2;
  fermion_cdf25 f25;
  fermion_rel2_cdf25 fr25;
  //fermion_rel2_cdf252 fr25;

  if (calibrate_ld || calibrate_cdf25 || calibrate_new) {
    cout << "----------------------------------------------------" << endl;
    cout << "Testing multiprecision" << endl;
    cout << "----------------------------------------------------" << endl;
    cout << endl;
  }

  if (calibrate_ld) {
    frld.verify_ti=true;
    long double v3=pcc_ld.part_calibrate<fermion_ld,fermion_rel2_ld>
      (fld,frld,"../../data/o2scl/fermion_deriv_cal.o2",
       true,false,false,0,true);
    t.test_rel<long double>(v3,0.0,2.0e-12,"calibrate 3");
  }

  if (calibrate_cdf25) {
    fr25.verify_ti=true;
    cpp_dec_float_25 v4=pcc_cdf25.part_calibrate<fermion_cdf25,
                                                 fermion_rel2_cdf25>
      (f25,fr25,"../../data/o2scl/fermion_deriv_cal.o2",
       true,false,false,0,true);
    t.test_rel_boost<cpp_dec_float_25>(v4,0.0,8.0e-13,"calibrate 4");
  }
  
  if (calibrate_new) {

    t.set_output_level(1);
    
    f.non_interacting=true;
    f.inc_rest_mass=true;
    f.g=2;
    fld.non_interacting=true;
    fld.inc_rest_mass=true;
    fld.g=2;
    f25.non_interacting=true;
    f25.inc_rest_mass=true;
    f25.g=2;

    double T=1;
    long double Tld=1;
    long double T25=1;
    double max=0.0;
    cpp_dec_float_25 maxld=0.0;
    int ret;

    if (false) {
      
      fr.verbose=2;
      frld.verbose=2;
      fr25.verbose=2;
      fr.verify_ti=true;
      frld.verify_ti=true;
      fr25.verify_ti=true;

      int lmot=3;
      int lpsi=-1;
      
      double psi=pow(10.0,((double)lpsi));
      double mot=pow(10.0,((double)lmot));
      long double psi_ld=pow(10,((long double)lpsi));
      long double mot_ld=pow(10,((long double)lmot));
      cpp_dec_float_25 psi_25=pow(10,((cpp_dec_float_25)lpsi));
      cpp_dec_float_25 mot_25=pow(10,((cpp_dec_float_25)lmot));
      
      f.g=2;
      f.m=mot*T;
      f.mu=psi*T+f.m;
      
      fld.g=2;
      fld.m=mot_ld*Tld;
      fld.mu=psi_ld*Tld+fld.m;
      
      f25.g=2;
      f25.m=mot_25*T25;
      f25.mu=psi_25*T25+f25.m;
      
      fr.calc_mu(f,T);
      //frld.fri.verbose=2;
      //frld.fri.it.verbose=2;
      //frld.fri.it2.verbose=1;
      frld.calc_mu(fld,Tld);
      fr25.calc_mu(f25,T25);
      
      cout << mot << " " << psi << " mu,m: " << f.mu << " " << f.m << endl;
      cout.setf(ios::left);
      cout.width(32);
      cout << dtos(f.n,0) << " ";
      cout.width(32);
      cout << dtos(f.ed,0) << " ";
      cout.width(32);
      cout << dtos(f.pr,0) << " ";
      cout.width(32);
      cout << dtos(f.en,0) << endl;
      cout.width(32);
      cout << dtos(fld.n,0) << " ";
      cout.width(32);
      cout << dtos(fld.ed,0) << " ";
      cout.width(32);
      cout << dtos(fld.pr,0) << " ";
      cout.width(32);
      cout << dtos(fld.en,0) << endl;
      cout.width(32);
      cout << dtos(f25.n,0) << " ";
      cout.width(32);
      cout << dtos(f25.ed,0) << " ";
      cout.width(32);
      cout << dtos(f25.pr,0) << " ";
      cout.width(32);
      cout << dtos(f25.en,0) << endl;
      cout << endl;
      cout.unsetf(ios::left);
      cout << "  " << abs(f.n-fld.n)/abs(fld.n) << " ";
      cout << abs(f.ed-fld.ed)/abs(fld.ed) << " ";
      cout << abs(f.pr-fld.pr)/abs(fld.pr) << " ";
      cout << abs(f.en-fld.en)/abs(fld.en) << endl;
      cout << "  " << abs(fld.n-f25.n)/abs(f25.n) << " ";
      cout << abs(fld.ed-f25.ed)/abs(f25.ed) << " ";
      cout << abs(fld.pr-f25.pr)/abs(f25.pr) << " ";
      cout << abs(fld.en-f25.en)/abs(f25.en) << endl;
      cout << endl;
      fr.verbose=0;
      frld.verbose=0;
      fr25.verbose=0;
      
    }
    
    for(int lmot=-3;lmot<=3;lmot++) {
      for(int lpsi=-3;lpsi<=1;lpsi++) {
        
        cout << lmot << " " << lpsi << " ";
        
        double psi=pow(10.0,((double)lpsi));
        double mot=pow(10.0,((double)lmot));
        long double psi_ld=pow(10,((long double)lpsi));
        long double mot_ld=pow(10,((long double)lmot));
        cpp_dec_float_25 psi_25=pow(10,((cpp_dec_float_25)lpsi));
        cpp_dec_float_25 mot_25=pow(10,((cpp_dec_float_25)lmot));
        
        f.g=2;
        f.m=mot*T;
        f.mu=psi*T+f.m;
        cout << " mu,m: " << f.mu << " " << f.m << endl;
        ret=fr.calc_mu(f,T);
        t.test_gen(ret==0,"calibrate new ret 1");
        
        fld.g=2;
        fld.m=mot_ld*Tld;
        fld.mu=psi_ld*Tld+fld.m;
        ret=frld.calc_mu(fld,Tld);
        t.test_gen(ret==0,"calibrate new ret 2");
        
        f25.g=2;
        f25.m=mot_25*T25;
        f25.mu=psi_25*T25+f25.m;
        ret=fr25.calc_mu(f25,T25);
        t.test_gen(ret==0,"calibrate new ret 3");
        
        cout << "  " << abs(f.n-fld.n)/abs(fld.n) << " ";
        cout << abs(f.ed-fld.ed)/abs(fld.ed) << " ";
        cout << abs(f.pr-fld.pr)/abs(fld.pr) << " ";
        cout << abs(f.en-fld.en)/abs(fld.en) << " ";
        
        if (abs(f.n-fld.n)/abs(fld.n)>max) {
          max=abs(f.n-fld.n)/abs(fld.n);
        }
        if (abs(f.ed-fld.ed)/abs(fld.ed)>max) {
          max=abs(f.ed-fld.ed)/abs(fld.ed);
        }
        if (abs(f.pr-fld.pr)/abs(fld.pr)>max) {
          max=abs(f.pr-fld.pr)/abs(fld.pr);
        }
        if (abs(f.en-fld.en)/abs(fld.en)>max) {
          max=abs(f.en-fld.en)/abs(fld.en);
        }
        cout << max << endl;
        
        cout << "  " << abs(fld.n-f25.n)/abs(f25.n) << " ";
        cout << abs(fld.ed-f25.ed)/abs(f25.ed) << " ";
        cout << abs(fld.pr-f25.pr)/abs(f25.pr) << " ";
        cout << abs(fld.en-f25.en)/abs(f25.en) << " ";
        
        if (abs(fld.n-f25.n)/abs(f25.n)>maxld) {
          maxld=abs(fld.n-f25.n)/abs(f25.n);
        }
        if (abs(fld.ed-f25.ed)/abs(f25.ed)>maxld) {
          maxld=abs(fld.ed-f25.ed)/abs(f25.ed);
        }
        if (abs(fld.pr-f25.pr)/abs(f25.pr)>maxld) {
          maxld=abs(fld.pr-f25.pr)/abs(f25.pr);
        }
        if (abs(fld.en-f25.en)/abs(f25.en)>maxld) {
          maxld=abs(fld.en-f25.en)/abs(f25.en);
        }
        cout << maxld << endl;

        cout << "  " << fr.last_method << " " << frld.last_method << " "
             << fr25.last_method << endl;

      }
    }

    t.test_abs(max,0.0,8.0e-13,"calibrate_new");
    t.test_abs(static_cast<double>(maxld),0.0,1.0e-16,"calibrate_new ld");

  }

#endif
  
  t.report();

  return 0;
}

