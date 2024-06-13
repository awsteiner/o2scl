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
#include <o2scl/fermion_rel.h>
#include <o2scl/test_mgr.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/eos_sn.h>

#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_const;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

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

  int first_test=0;
  
  fermion_rel2 fr2;
  //fermion_rel2_ld fr2ld;
  //fermion_rel2_cdf25 fr225;

  // An exhaustive comparison of the two algorithms at
  // various levels of precision

  cout.precision(4);
  int count=0;

  cout << " cnt m          T           mu/n       "
       << "d-ld  ld-25 1 vs. 2" << endl;
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
          
          int idn, iden, ildn, ilden;
          if (f.n==0) {
            idn=99;
          } else if (f.n==fld.n) {
            idn=17;
          } else {
            idn=((int)(-log10(fabs(f.n-fld.n)/fabs(fld.n))));
          }
          if (f.en==0) {
            iden=99;
          } else if (f.en==fld.en) {
            iden=17;
          } else {
            iden=((int)(-log10(fabs(f.en-fld.en)/fabs(fld.en))));
          }
          if (fld.n==static_cast<long double>(f25.n)) {
            ildn=21;
          } else {
            ildn=((int)(-log10(fabs(fld.n-static_cast<long double>(f25.n))/
                               fabs(static_cast<long double>(f25.n)))));
          }
          if (fld.en==static_cast<long double>(f25.en)) {
            ilden=21;
          } else {
            ilden=((int)(-log10(fabs(fld.en-static_cast<long double>(f25.en))/
                                fabs(static_cast<long double>(f25.en)))));
          }
          if (idn<=0) {
            cout << "Problem in idn+ 1: " << endl;
            cout << dtos(f.n,-1) << " "
                 << dtos(fld.n,-1) << " " << idn << endl;
            exit(-1);
          }
          if (iden<=0) {
            cout << "Problem in iden+ 1: " << endl;
            cout << dtos(f.en,-1) << " "
                 << dtos(fld.en,-1) << " " << iden << endl;
            exit(-1);
          }
          if (ildn<=0) {
            cout << "Problem in ildn+ 1: " << endl;
            cout << dtos(fld.n,-1) << " "
                 << dtos(f25.n,-1) << " " << ildn << endl;
            cout << fabs(fld.n-static_cast<long double>(f25.n)) << endl;
            cout << fabs(fld.n-static_cast<long double>(f25.n))/
              fabs(static_cast<long double>(f25.n)) << endl;
            exit(-1);
          }
          if (ilden<=0) {
            cout << "Problem in ilden+ 1: " << endl;
            cout << dtos(fld.en,-1) << " "
                 << dtos(f25.en,-1) << " " << ilden << endl;
            cout << fabs(fld.en-static_cast<long double>(f25.en)) << endl;
            cout << fabs(fld.en-static_cast<long double>(f25.en))/
              fabs(static_cast<long double>(f25.en)) << endl;
            exit(-1);
          }
          
          cout.width(2);
          cout << idn << " ";
          cout.width(2);
          cout << iden << " ";
          cout.width(2);
          cout << ildn << " ";
          cout.width(2);
          cout << ilden << " ";
          
          double n1=f.n;
          double en1=f.en;
          long double nld1=fld.n;
          long double enld1=fld.en;
          cpp_dec_float_25 n251=f25.n;
          cpp_dec_float_25 en251=f25.en;
          
          fr2.calc_mu(f,T);
          //frld.calc_mu(fld,Tld);
          //fr25.calc_mu(f25,T25);
          
          int jdn, jden;
          if (f.n==n1) {
            jdn=17;
          } else {
            jdn=((int)(-log10(fabs(f.n-n1)/fabs(n1))));
          }
          if (f.en==en1) {
            jden=17;
          } else {
            jden=((int)(-log10(fabs(f.en-en1)/fabs(en1))));
          }
          
          cout.width(2);
          cout << jdn << " ";
          cout.width(2);
          cout << jden << " mu>0, calc_mu()." << endl;
          
          if (jdn<=0) {
            cout << "Problem in jdn+ 1: " << endl;
            cout << dtos(f.n,-1) << " "
                 << dtos(n1,-1) << " " << jdn << endl;
            cout << fabs(f.n-n1)/fabs(n1) << endl;
            cout << f.n << " " << n1 << endl;
            exit(-1);
          }
          if (jden<=9) {
            cout << "Problem in jden+ 1: " << endl;
            cout << dtos(f.en,-1) << " "
                 << dtos(en1,-1) << " " << jden << endl;
            cout << fabs(f.en-en1)/fabs(en1) << endl;
            cout << dtos(f.en,-1) << endl;
            cout << dtos(en1,-1) << endl;
            exit(-1);
          }
        
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

          int kdn, kden, kldn, klden;
          if (f.n==0) {
            kdn=99;
          } else if (f.n==fld.n) {
            kdn=17;
          } else {
            kdn=((int)(-log10(fabs(f.n-fld.n)/fabs(fld.n))));
          }
          if (f.en==0) {
            kden=99;
          } else if (f.en==fld.en) {
            kden=17;
          } else {
            kden=((int)(-log10(fabs(f.en-fld.en)/fabs(fld.en))));
          }
          if (fld.n==static_cast<long double>(f25.n)) {
            kldn=21;
          } else {
            kldn=((int)(-log10(fabs(fld.n-static_cast<long double>(f25.n))/
                               fabs(static_cast<long double>(f25.n)))));
          }
          if (fld.en==static_cast<long double>(f25.en)) {
            klden=21;
          } else {
            klden=((int)(-log10(fabs(fld.en-static_cast<long double>(f25.en))/
                                fabs(static_cast<long double>(f25.en)))));
          }
          if (kdn<=0) {
            cout << "Problem in kdn+ 1: " << endl;
            cout << dtos(f.n,-1) << " "
                 << dtos(fld.n,-1) << " " << kdn << endl;
            cout << f.n << " " << fld.n << endl;
            exit(-1);
          }
          if (kldn<=0) {
            cout << "Problem in kldn+ 1: " << endl;
            cout << dtos(fld.n,-1) << " "
                 << dtos(f25.n,-1) << " " << kldn << endl;
            cout << fabs(fld.n-static_cast<long double>(f25.n)) << endl;
            cout << fabs(fld.n-static_cast<long double>(f25.n))/
              fabs(static_cast<long double>(f25.n)) << endl;
            exit(-1);
          }
        
          cout.width(2);
          cout << kdn << " ";
          cout.width(2);
          cout << kden << " ";
          cout.width(2);
          cout << kldn << " ";
          cout.width(2);
          cout << klden << " ";
        
          double n1=f.n;
          double en1=f.en;
          long double nld1=fld.n;
          long double enld1=fld.en;
          cpp_dec_float_25 n251=f25.n;
          cpp_dec_float_25 en251=f25.en;
        
          fr2.pair_mu(f,T);

          int ldn, lden;
          if (f.n==n1) {
            ldn=17;
          } else {
            ldn=((int)(-log10(fabs(f.n-n1)/fabs(n1))));
          }
          if (f.en==en1) {
            lden=17;
          } else {
            lden=((int)(-log10(fabs(f.en-en1)/fabs(en1))));
          }

          cout.width(2);
          cout << ldn << " ";
          cout.width(2);
          cout << lden << " mu>0, pair_mu()." << endl;
        
          if (ldn<=0) {
            cout << "Problem in ldn+ 1: " << endl;
            cout << dtos(f.n,-1) << " "
                 << dtos(n1,-1) << " " << ldn << endl;
            exit(-1);
          }
          if (lden<=0) {
            cout << "Problem in lden+ 1: " << endl;
            cout << dtos(f.en,-1) << " "
                 << dtos(en1,-1) << " " << lden << endl;
            exit(-1);
          }

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

          int idn, iden, ildn, ilden;
          if (f.n==0) {
            idn=99;
          } else if (f.n==fld.n) {
            idn=17;
          } else {
            idn=((int)(-log10(fabs(f.n-fld.n)/fabs(fld.n))));
          }
          if (f.en==0) {
            iden=99;
          } else if (f.en==fld.en) {
            iden=17;
          } else {
            iden=((int)(-log10(fabs(f.en-fld.en)/fabs(fld.en))));
          }
          if (fld.n==static_cast<long double>(f25.n)) {
            ildn=21;
          } else {
            ildn=((int)(-log10(fabs(fld.n-static_cast<long double>(f25.n))/
                               fabs(static_cast<long double>(f25.n)))));
          }
          if (fld.en==static_cast<long double>(f25.en)) {
            ilden=21;
          } else {
            ilden=((int)(-log10(fabs(fld.en-static_cast<long double>(f25.en))/
                                fabs(static_cast<long double>(f25.en)))));
          }
          if (idn<=0) {
            cout << "Problem in idn- 1: " << endl;
            cout << dtos(f.n,-1) << " "
                 << dtos(fld.n,-1) << " " << idn << endl;
            cout << f.n << " " << fld.n << endl;
            cout << (f.n==0) << endl;
            exit(-1);
          }
          if (iden<=0) {
            cout << "Problem in iden- 1: " << endl;
            cout << dtos(f.en,-1) << " "
                 << dtos(fld.en,-1) << " " << iden << endl;
            exit(-1);
          }
          if (ildn<=0) {
            cout << "Problem in ildn- 1: " << endl;
            cout << dtos(fld.n,-1) << " "
                 << dtos(f25.n,-1) << " " << ildn << endl;
            cout << fabs(fld.n-static_cast<long double>(f25.n)) << endl;
            cout << fabs(fld.n-static_cast<long double>(f25.n))/
              fabs(static_cast<long double>(f25.n)) << endl;
            exit(-1);
          }
          if (ilden<=0) {
            cout << "Problem in ilden- 1: " << endl;
            cout << dtos(fld.en,-1) << " "
                 << dtos(f25.en,-1) << " " << ilden << endl;
            cout << fabs(fld.en-static_cast<long double>(f25.en)) << endl;
            cout << fabs(fld.en-static_cast<long double>(f25.en))/
              fabs(static_cast<long double>(f25.en)) << endl;
            exit(-1);
          }
        
          cout.width(2);
          cout << idn << " ";
          cout.width(2);
          cout << iden << " ";
          cout.width(2);
          cout << ildn << " ";
          cout.width(2);
          cout << ilden << " ";
        
          double n1=f.n;
          double en1=f.en;
          long double nld1=fld.n;
          long double enld1=fld.en;
          cpp_dec_float_25 n251=f25.n;
          cpp_dec_float_25 en251=f25.en;
        
          fr2.calc_mu(f,T);
          //frld.calc_mu(fld,Tld);
          //fr25.calc_mu(f25,T25);

          int jdn, jden;
          if (f.n==n1) {
            jdn=17;
          } else {
            jdn=((int)(-log10(fabs(f.n-n1)/fabs(n1))));
          }
          if (f.en==en1) {
            jden=17;
          } else {
            jden=((int)(-log10(fabs(f.en-en1)/fabs(en1))));
          }

          cout.width(2);
          cout << jdn << " ";
          cout.width(2);
          cout << jden << " mu<0, calc_mu()." << endl;
        
          if (jdn<=0) {
            cout << "Problem in jdn- 1: " << endl;
            cout << dtos(f.n,-1) << " "
                 << dtos(n1,-1) << " " << jdn << endl;
            cout << fabs(f.n-n1)/fabs(n1) << endl;
            cout << f.n << " " << n1 << endl;
            exit(-1);
          }
          if (jden<=0) {
            cout << "Problem in jden- 1: " << endl;
            cout << dtos(f.en,-1) << " "
                 << dtos(en1,-1) << " " << jden << endl;
            cout << fabs(f.en-en1)/fabs(en1) << endl;
            cout << f.en << " " << en1 << endl;
            exit(-1);
          }

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
        
          int kdn, kden, kldn, klden;
          if (f.n==0) {
            kdn=99;
          } else if (f.n==fld.n) {
            kdn=17;
          } else {
            kdn=((int)(-log10(fabs(f.n-fld.n)/fabs(fld.n))));
          }
          if (f.en==0) {
            kden=99;
          } else if (f.en==fld.en) {
            kden=17;
          } else {
            kden=((int)(-log10(fabs(f.en-fld.en)/fabs(fld.en))));
          }
          if (fld.n==static_cast<long double>(f25.n)) {
            kldn=21;
          } else {
            kldn=((int)(-log10(fabs(fld.n-static_cast<long double>(f25.n))/
                               fabs(static_cast<long double>(f25.n)))));
          }
          if (fld.en==static_cast<long double>(f25.en)) {
            klden=21;
          } else {
            klden=((int)(-log10(fabs(fld.en-static_cast<long double>(f25.en))/
                                fabs(static_cast<long double>(f25.en)))));
          }
          if (kdn<=0) {
            cout << "Problem in kdn- 1: " << endl;
            cout << dtos(f.n,-1) << " "
                 << dtos(fld.n,-1) << " " << kdn << endl;
            cout << f.n << " " << fld.n << endl;
            exit(-1);
          }
          if (kldn<=0) {
            cout << "Problem in kldn- 1: " << endl;
            cout << dtos(fld.n,-1) << " "
                 << dtos(f25.n,-1) << " " << kldn << endl;
            cout << fabs(fld.n-static_cast<long double>(f25.n)) << endl;
            cout << fabs(fld.n-static_cast<long double>(f25.n))/
              fabs(static_cast<long double>(f25.n)) << endl;
            exit(-1);
          }
        
          cout.width(2);
          cout << kdn << " ";
          cout.width(2);
          cout << kden << " ";
          cout.width(2);
          cout << kldn << " ";
          cout.width(2);
          cout << klden << " ";
        
          double n1=f.n;
          double en1=f.en;
          long double nld1=fld.n;
          long double enld1=fld.en;
          cpp_dec_float_25 n251=f25.n;
          cpp_dec_float_25 en251=f25.en;
        
          fr2.pair_mu(f,T);

          int ldn, lden;
          if (f.n==n1) {
            ldn=17;
          } else {
            ldn=((int)(-log10(fabs(f.n-n1)/fabs(n1))));
          }
          if (f.en==en1) {
            lden=17;
          } else {
            lden=((int)(-log10(fabs(f.en-en1)/fabs(en1))));
          }

          cout.width(2);
          cout << ldn << " ";
          cout.width(2);
          cout << lden << " mu<0, pair_mu()." << endl;
        
          if (ldn<=0) {
            cout << "Problem in ldn- 1: " << endl;
            cout << dtos(f.n,-1) << " "
                 << dtos(n1,-1) << " " << ldn << endl;
            exit(-1);
          }
          if (lden<=0) {
            cout << "Problem in lden- 1: " << endl;
            cout << dtos(f.en,-1) << " "
                 << dtos(en1,-1) << " " << lden << endl;
            exit(-1);
          }
        
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
            fr2.verbose=2;
          }
        
          f.mu=f.m;
          fr.calc_density(f,T);
          fld.mu=fld.m;
          frld.calc_density(fld,Tld);
          f25.mu=f25.m;
          fr25.calc_density(f25,T25);
        
          int idmu, iden, ildmu, ilden;
          if (f.mu==fld.mu) {
            idmu=17;
          } else {
            idmu=((int)(-log10(fabs(f.mu-fld.mu)/fabs(fld.mu))));
          }
          if (f.en==fld.en) {
            iden=17;
          } else {
            iden=((int)(-log10(fabs(f.en-fld.en)/fabs(fld.en))));
          }
          if (fld.mu==static_cast<long double>(f25.mu)) {
            ildmu=21;
          } else {
            ildmu=((int)(-log10(fabs(fld.mu-static_cast<long double>(f25.mu))/
                                fabs(static_cast<long double>(f25.mu)))));
          }
          if (fld.en==static_cast<long double>(f25.en)) {
            ilden=21;
          } else {
            ilden=((int)(-log10(fabs(fld.en-static_cast<long double>(f25.en))/
                                fabs(static_cast<long double>(f25.en)))));
          }
          if (idmu<=0) {
            cout << "Problem in idmu 1: " << endl;
            cout << dtos(f.mu,-1) << " "
                 << dtos(fld.mu,-1) << " " << idmu << endl;
            exit(-1);
          }
          if (iden<=0) {
            cout << "Problem in iden 1: " << endl;
            cout << dtos(f.en,-1) << " "
                 << dtos(fld.en,-1) << " " << iden << endl;
            exit(-1);
          }
          if (ildmu<=0) {
            cout << "Problem in ildmu 1: " << endl;
            cout << dtos(fld.mu,-1) << " "
                 << dtos(f25.mu,-1) << " " << ildmu << endl;
            cout << fabs(fld.mu-static_cast<long double>(f25.mu)) << endl;
            cout << fabs(fld.mu-static_cast<long double>(f25.mu))/
              fabs(static_cast<long double>(f25.mu)) << endl;
            exit(-1);
          }
          if (ilden<=0) {
            cout << "Problem in ilden 1: " << endl;
            cout << dtos(fld.en,-1) << " "
                 << dtos(f25.en,-1) << " " << ilden << endl;
            cout << fabs(fld.en-static_cast<long double>(f25.en)) << endl;
            cout << fabs(fld.en-static_cast<long double>(f25.en))/
              fabs(static_cast<long double>(f25.en)) << endl;
            exit(-1);
          }
        
          cout.width(2);
          cout << idmu << " ";
          cout.width(2);
          cout << iden << " ";
          cout.width(2);
          cout << ildmu << " ";
          cout.width(2);
          cout << ilden << " ";
        
          double mu1=f.mu;
          double en1=f.en;
          long double mu1ld=fld.mu;
          long double en1ld=fld.en;
          cpp_dec_float_25 mu125=f25.mu;
          cpp_dec_float_25 en125=f25.en;
        
          f.mu=f.m;
          fr2.calc_density(f,T);

          int jdmu, jden;
          if (f.mu==mu1) {
            jdmu=17;
          } else {
            jdmu=((int)(-log10(fabs(f.mu-mu1)/fabs(mu1))));
          }
          if (f.en==en1) {
            jden=17;
          } else {
            jden=((int)(-log10(fabs(f.en-en1)/fabs(en1))));
          }

          cout.width(2);
          cout << jdmu << " ";
          cout.width(2);
          cout << jden << " calc_density()." << endl;
        
          if (jdmu<=0) {
            cout << "Problem in jdmu 1: " << endl;
            cout << dtos(f.mu,-1) << " "
                 << dtos(mu1,-1) << " " << jdmu << endl;
            cout << fabs(f.mu-mu1)/fabs(mu1) << endl;
            cout << f.mu << " " << mu1 << endl;
            exit(-1);
          }
          if (jden<=0) {
            cout << "Problem in jden 1: " << endl;
            cout << dtos(f.en,-1) << " "
                 << dtos(en1,-1) << " " << jden << endl;
            cout << fabs(f.en-en1)/fabs(en1) << endl;
            cout << dtos(f.en,-1) << endl;
            cout << dtos(en1,-1) << endl;
            exit(-1);
          }

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
            fr2.verbose=2;
          }
        
          f.mu=f.m;
          fr.pair_density(f,T);
          fld.mu=fld.m;
          frld.pair_density(fld,Tld);
          f25.mu=f25.m;
          fr25.pair_density(f25,T25);
        
          int kdmu, kden, kldmu, klden;
          if (f.mu==0) {
            kdmu=99;
          } else if (f.mu==fld.mu) {
            kdmu=17;
          } else {
            kdmu=((int)(-log10(fabs(f.mu-fld.mu)/fabs(fld.mu))));
          }
          if (f.en==0) {
            kden=99;
          } else if (f.en==fld.en) {
            kden=17;
          } else {
            kden=((int)(-log10(fabs(f.en-fld.en)/fabs(fld.en))));
          }
          if (fld.mu==static_cast<long double>(f25.mu)) {
            kldmu=21;
          } else {
            kldmu=((int)(-log10(fabs(fld.mu-static_cast<long double>(f25.mu))/
                                fabs(static_cast<long double>(f25.mu)))));
          }
          if (fld.en==static_cast<long double>(f25.en)) {
            klden=21;
          } else {
            klden=((int)(-log10(fabs(fld.en-static_cast<long double>(f25.en))/
                                fabs(static_cast<long double>(f25.en)))));
          }
          if (kdmu<=0) {
            cout << "Problem in kdmu 1: " << endl;
            cout << dtos(f.mu,-1) << " "
                 << dtos(fld.mu,-1) << " " << kdmu << endl;
            cout << f.mu << " " << fld.mu << endl;
            exit(-1);
          }
          if (kldmu<=0) {
            cout << "Problem in kldmu 1: " << endl;
            cout << dtos(fld.mu,-1) << " "
                 << dtos(f25.mu,-1) << " " << kldmu << endl;
            cout << fabs(fld.mu-static_cast<long double>(f25.mu)) << endl;
            cout << fabs(fld.mu-static_cast<long double>(f25.mu))/
              fabs(static_cast<long double>(f25.mu)) << endl;
            exit(-1);
          }
        
          cout.width(2);
          cout << kdmu << " ";
          cout.width(2);
          cout << kden << " ";
          cout.width(2);
          cout << kldmu << " ";
          cout.width(2);
          cout << klden << " ";
        
          double mu1=f.mu;
          double en1=f.en;
          long double mu1ld=fld.mu;
          long double en1ld=fld.en;
          cpp_dec_float_25 mu125=f25.mu;
          cpp_dec_float_25 en125=f25.en;
        
          f.mu=f.m;
          fr2.pair_density(f,T);
        
          int ldmu, lden;
          if (f.mu==mu1) {
            ldmu=17;
          } else {
            ldmu=((int)(-log10(fabs(f.mu-mu1)/fabs(mu1))));
          }
          if (f.en==en1) {
            lden=17;
          } else {
            lden=((int)(-log10(fabs(f.en-en1)/fabs(en1))));
          }
        
          cout.width(2);
          cout << ldmu << " ";
          cout.width(2);
          cout << lden << " pair_density()." << endl;
        
          if (ldmu<=0) {
            cout << "Problem in ldmu 1: " << endl;
            cout << dtos(f.mu,-1) << " "
                 << dtos(mu1,-1) << " " << ldmu << endl;
            exit(-1);
          }
          if (lden<=0) {
            cout << "Problem in lden 1: " << endl;
            cout << dtos(f.en,-1) << " "
                 << dtos(en1,-1) << " " << lden << endl;
            exit(-1);
          }

        } else {

          cout << endl;

        }
          
        count++;
        
      }
    }
  }
  
  t.report();

  return 0;
}

