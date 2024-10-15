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
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/test_mgr.h>
// For access to global convert_units object
#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(int argc, char *argv[]) {

  cout.setf(ios::scientific);
  
  test_mgr t;
  t.set_output_level(1);
  
  if (true) {
    // Testing dndmu
    
    fermion_deriv fd(5.0,2.0);
    fermion_deriv_rel fdr;
    
    fd.inc_rest_mass=false;
    
    fd.n=0.5;
    fdr.calc_density(fd,0);
    double mu1=fd.mu;
    fd.n=0.5+1.0e-4;
    fdr.calc_density(fd,0);
    double mu2=fd.mu;
    fd.n=0.5;
    fdr.calc_density(fd,0);
    fdr.calc_deriv_zerot(fd);
    t.test_rel(1.0e-4/(mu2-mu1),fd.dndmu,1.0e-4,"dndmu at zerot.");
    
    fd.inc_rest_mass=true;
    
    fd.n=0.5;
    fdr.calc_density(fd,0);
    mu1=fd.mu;
    fd.n=0.5+1.0e-4;
    fdr.calc_density(fd,0);
    mu2=fd.mu;
    fd.n=0.5;
    fdr.calc_density(fd,0);
    fdr.calc_deriv_zerot(fd);
    t.test_rel(1.0e-4/(mu2-mu1),fd.dndmu,1.0e-4,"dndmu at zerot.");

  }
  
  if (true) {
    
    fermion_deriv fd(5.0,2.0);
    double T;

    fd.inc_rest_mass=false;

    fermion_deriv_rel fdr;

    // -----------------------------------------------------------------
    // Test the specific heat of degenerate fermions (As a reminder,
    // note that the degenerate specific heat differs for relativistic
    // and non-relativistic particles by a factor of two. Below is the
    // case for relativistic particles.) The code below demonstrates
    // that the computation of the specific heat (and that of the
    // entropy and the pressure) fails for sufficiently low
    // temperatures, i.e. in the extremely degenerate case. [AWS,
    // 12/20/09 - This now works better.]

    fdr.method=fermion_deriv_rel::automatic;

    cout << "Test degenerate specific heat:\n" << endl;
    cout << "err T(MeV)     pr           en           "
         << "C            C(deg appx)  C(next term)" << endl;
    fd.mu=0.0;
    for (T=3.0/hc_mev_fm;T>=0.001/hc_mev_fm;T/=3.0) {

      fd.init(o2scl_const::mass_electron_f<double>()*
              o2scl_settings.get_convert_units().convert("kg","1/fm",1.0),2.0);
      fd.non_interacting=true;
      fd.n=0.2;
      fdr.calc_density(fd,T);
      fd.kf=cbrt(6.0*pi2/fd.g*fd.n);
      cout << T*hc_mev_fm << " " << fd.pr << " " << fd.en << " "
           << T/fd.n*(fd.dsdT-fd.dndT*fd.dndT/fd.dndmu) << " "
           << o2scl_const::pi*o2scl_const::pi*T/fd.mu << " " 
           << pow(o2scl_const::pi/fd.kf,2.0)*T*fd.mu-
        pow(o2scl_const::pi*T/fd.kf,3.0)*o2scl_const::pi/15.0*
        (5.0*pow(fd.m,4.0)+4.0*fd.m*fd.m*fd.mu*fd.mu+14.0*pow(fd.mu,4.0)) 
           << endl;
    
      t.test_rel(T/fd.n*(fd.dsdT-fd.dndT*fd.dndT/fd.dndmu),
                 o2scl_const::pi*o2scl_const::pi*T/fd.mu,1.0e-2,"sh1");
      t.test_rel(T/fd.n*(fd.dsdT-fd.dndT*fd.dndT/fd.dndmu),
                 fd.en/fd.n,1.0e-3,"sh2");
    
      // Compare with direct differentiation
      fdr.calc_density(fd,T);
      double en1=fd.en;
      double h=1.0e-4;
      fdr.calc_density(fd,T+h);
      double en2=fd.en;
      if (true || T>0.1/hc_mev_fm) {
        cout << "\t" << (en2-en1)/h*T/fd.n << endl;
        t.test_rel(T/fd.n*(fd.dsdT-fd.dndT*fd.dndT/fd.dndmu),
                   (en2-en1)/h*T/fd.n,1.0e-2,"sh3");
      }
    }
    cout << endl;

  }

  // -----------------------------------------------------------------

  if (true) {
  
    fermion_deriv fd;
    fermion_deriv_rel fdr;
    
    fd.g=2.0;
    fd.inc_rest_mass=true;
    fd.non_interacting=true;
  
    for(int im=-2;im<=1;im++) {
    
      fd.m=im;
      fd.m=pow(10,fd.m);
    
      for(int iT=-2;iT<=1;iT++) {
      
        double T=iT;
        T=pow(10,T);
      
        for(int imu=-2;imu<=1;imu++) {
        
          fd.mu=imu;
          fd.mu=pow(10,fd.mu);
        
          cout.width(2);
          cout << im << " ";
          cout.width(2);
          cout << iT << " ";
          cout.width(2);
          cout << imu << " ";

          fdr.method=fermion_deriv_rel::direct;
          fdr.calc_mu(fd,T);
          double dndmu1=fd.dndmu;
          double dndT1=fd.dndT;
          double dsdT1=fd.dsdT;
          fdr.method=fermion_deriv_rel::by_parts;
          fdr.calc_mu(fd,T);
          cout << dndmu1 << " " << fd.dndmu << " ";
          cout.width(2);
          cout << count_digits_same(dndmu1,fd.dndmu) << " ";
          cout.width(2);
          cout << count_digits_same(dndT1,fd.dndT) << " ";
          cout.width(2);
          cout << count_digits_same(dsdT1,fd.dsdT) << " ";

          if (false) {
            
            fdr.multip=true;
            
            fdr.method=fermion_deriv_rel::direct;
            fdr.calc_mu(fd,T);
            dndmu1=fd.dndmu;
            dndT1=fd.dndT;
            dsdT1=fd.dsdT;
            fdr.method=fermion_deriv_rel::by_parts;
            fdr.calc_mu(fd,T);
            cout << dndmu1 << " " << fd.dndmu << " ";
            cout.width(2);
            cout << count_digits_same(dndmu1,fd.dndmu) << " ";
            cout.width(2);
            cout << count_digits_same(dndT1,fd.dndT) << " ";
            cout.width(2);
            cout << count_digits_same(dsdT1,fd.dsdT) << endl;
            
          } else {
            
            cout << endl;
            
          }

          fdr.multip=false;
        }
      }
    }

  }

  fermion_deriv_rel_ld fdrl;
  
  t.report();

  return 0;
}
