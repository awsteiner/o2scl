/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner

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
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/fermion_eff.h>
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
  
  fermion_deriv sf(5.0,2.0);
  fermion_deriv sf2(5.0,2.0);
  fermion ef(5.0,2.0);
  double T;

  sf.inc_rest_mass=false;

  fermion_deriv_rel snf;
  fermion_eff eff;
  
  // -----------------------------------------------------------------
  // Test the specific heat of degenerate fermions (As a reminder,
  // note that the degenerate specific heat differs for relativistic
  // and non-relativistic particles by a factor of two. Below is the
  // case for relativistic particles.) The code below demonstrates
  // that the computation of the specific heat (and that of the
  // entropy and the pressure) fails for sufficiently low temperatures,
  // i.e. in the extremely degenerate case. [12/20/09 - This now
  // works better]

  snf.method=fermion_deriv_rel::automatic;

  cout << "Test degenerate specific heat:\n" << endl;
  cout << "err T(MeV)     pr           en           "
       << "C            C(deg appx)  C(next term)" << endl;
  sf.mu=0.0;
  for (T=3.0/hc_mev_fm;T>=0.001/hc_mev_fm;T/=3.0) {

    sf.init(o2scl_mks::mass_electron*
	    o2scl_settings.get_convert_units().convert("kg","1/fm",1.0),2.0);
    sf.non_interacting=true;
    sf.n=0.2;
    snf.calc_density(sf,T);
    sf.kf=cbrt(6.0*pi2/sf.g*sf.n);
    cout << T*hc_mev_fm << " " << sf.pr << " " << sf.en << " "
	 << T/sf.n*(sf.dsdT-sf.dndT*sf.dndT/sf.dndmu) << " "
	 << o2scl_const::pi*o2scl_const::pi*T/sf.mu << " " 
	 << pow(o2scl_const::pi/sf.kf,2.0)*T*sf.mu-
      pow(o2scl_const::pi*T/sf.kf,3.0)*o2scl_const::pi/15.0*
      (5.0*pow(sf.m,4.0)+4.0*sf.m*sf.m*sf.mu*sf.mu+14.0*pow(sf.mu,4.0)) 
	 << endl;
    
    t.test_rel(T/sf.n*(sf.dsdT-sf.dndT*sf.dndT/sf.dndmu),
	       o2scl_const::pi*o2scl_const::pi*T/sf.mu,1.0e-2,"sh1");
    t.test_rel(T/sf.n*(sf.dsdT-sf.dndT*sf.dndT/sf.dndmu),
	       sf.en/sf.n,1.0e-3,"sh2");
    
    // Compare with direct differentiation
    snf.calc_density(sf,T);
    double en1=sf.en;
    double h=1.0e-4;
    snf.calc_density(sf,T+h);
    double en2=sf.en;
    if (true || T>0.1/hc_mev_fm) {
      cout << "\t" << (en2-en1)/h*T/sf.n << endl;
      t.test_rel(T/sf.n*(sf.dsdT-sf.dndT*sf.dndT/sf.dndmu),
		 (en2-en1)/h*T/sf.n,1.0e-2,"sh3");
    }
  }
  cout << endl;

  fermion_deriv sfx(1.0,2.0);
  sfx.inc_rest_mass=false;

  fermion_rel fr;
  double pc_fr=part_calibrate<fermion,fermion_rel>
    (ef,fr,true,"../../data/o2scl/fermion_cal2.o2",false,1,true);
  cout << pc_fr << endl;
  
  cout << "----------------------------------------------------" << endl;
  cout << "Function deriv_calibrate() method=direct." << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  snf.method=fermion_deriv_rel::direct;
  
  double pc_fdr_dir=part_calibrate<fermion_deriv,fermion_deriv_rel>
    (sfx,snf,true,"../../data/o2scl/fermion_cal2.o2",false,1,true);
  cout << pc_fdr_dir << endl;
  t.test_rel(pc_fr,pc_fdr_dir,1.0e-6,"nonderiv vs. deriv");
  
  double dc_dir=snf.deriv_calibrate(sfx,1,"../../data/o2scl/fermion_cal2.o2");
  cout << "Deriv_Calibrate: " << dc_dir << endl;
  t.test_rel(dc_dir,0.0,8.0e-6,"deriv_calibrate direct");

  cout << "----------------------------------------------------" << endl;
  cout << "Function deriv_calibrate() method=by_parts." << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  snf.method=fermion_deriv_rel::by_parts;
  
  double pc_fdr_byp=part_calibrate<fermion_deriv,fermion_deriv_rel>
    (sfx,snf,true,"../../data/o2scl/fermion_cal2.o2",false,1,true);
  cout << pc_fdr_byp << endl;
  t.test_rel(pc_fr,pc_fdr_byp,1.0e-6,"nonderiv vs. deriv");

  double dc_byp=snf.deriv_calibrate(sfx,1,"../../data/o2scl/fermion_cal2.o2");
  cout << "Deriv_Calibrate: " << dc_byp << endl;
  t.test_rel(dc_byp,0.0,1.0e-6,"deriv_calibrate by parts");
  cout << endl;

  cout << "----------------------------------------------------" << endl;
  cout << "Function deriv_calibrate() method=automatic." << endl;
  cout << "----------------------------------------------------" << endl;
  cout << endl;

  snf.method=fermion_deriv_rel::automatic;
  
  double pc_fdr_auto=part_calibrate<fermion_deriv,fermion_deriv_rel>
    (sfx,snf,true,"../../data/o2scl/fermion_cal2.o2",false,1,true);
  cout << pc_fdr_auto << endl;
  t.test_rel(pc_fr,pc_fdr_auto,1.0e-6,"nonderiv vs. deriv");

  double dc_auto=snf.deriv_calibrate(sfx,1,"../../data/o2scl/fermion_cal2.o2");
  cout << "Deriv_Calibrate: " << dc_auto << endl;
  t.test_rel(dc_auto,0.0,1.0e-6,"deriv_calibrate auto");
  cout << endl;

  cout << "------------------------------------------------------" << endl;

  // New function for derivative calibration
  double val2y=part_deriv_calibrate<fermion_deriv,fermion_deriv_rel>
    (sfx,snf,true,"../../data/o2scl/fermion_cal2.o2",1,true);
  cout << val2y << endl;

  cout << "------------------------------------------------------" << endl;
  {

    cout << "dndm test" << endl;
    snf.method=fermion_deriv_rel::direct;
    double d1, d2, eps=1.0e-4;
    double dndmu, dndT, dsdT, dndm;

    sf.inc_rest_mass=false;
    sf.non_interacting=true;

    T=1.0;
    sf.m=5.0;
    sf.ms=5.0;
    sf.mu=1.0*T-sf.m+sf.ms;

    sf.mu=1.0*T-sf.m+sf.ms+eps;
    snf.calc_mu(sf,T);
    d1=sf.n;
    sf.mu=1.0*T-sf.m+sf.ms;
    snf.calc_mu(sf,T);
    d2=sf.n;
    dndmu=(d1-d2)/eps;

    snf.calc_mu(sf,T+eps);
    d1=sf.n;
    snf.calc_mu(sf,T);
    d2=sf.n;
    dndT=(d1-d2)/eps;

    snf.calc_mu(sf,T+eps);
    d1=sf.en;
    snf.calc_mu(sf,T);
    d2=sf.en;
    dsdT=(d1-d2)/eps;

    sf.non_interacting=false;
    sf.nu=sf.mu;
    sf.ms=5.0+eps;
    snf.calc_mu(sf,T);
    d1=sf.n;

    sf.ms=5.0;
    snf.calc_mu(sf,T);
    d2=sf.n;
    dndm=(d1-d2)/eps;
    sf.non_interacting=true;

    snf.calc_mu(sf,T);
    cout << "sf: " << sf.dndmu << " " << sf.dndT << " " << sf.dsdT << endl;
    cout << "un: " << snf.unc.dndmu << " " << snf.unc.dndT << " " 
	 << snf.unc.dsdT << endl;
    cout << "nu: " << dndmu << " " << dndT << " " << dsdT << endl;
    double dndm2=3.0*sf.n/sf.m-(sf.dndT+sf.mu/T*sf.dndmu)*T/sf.m-sf.dndmu;
    cout << dndm2 << endl;
    cout << endl;

  }

  {

    cout << "dndm test" << endl;
    snf.method=fermion_deriv_rel::direct;
    double d1, d2, eps=1.0e-4;
    double dndmu, dndT, dsdT, dndm;

    sf.inc_rest_mass=true;
    sf.non_interacting=true;

    T=1.0;
    sf.m=5.0;
    sf.ms=5.0;
    sf.mu=1.0*T;

    sf.mu=1.0*T+eps;
    snf.calc_mu(sf,T);
    d1=sf.n;
    sf.mu=1.0*T;
    snf.calc_mu(sf,T);
    d2=sf.n;
    dndmu=(d1-d2)/eps;

    snf.calc_mu(sf,T+eps);
    d1=sf.n;
    snf.calc_mu(sf,T);
    d2=sf.n;
    dndT=(d1-d2)/eps;

    snf.calc_mu(sf,T+eps);
    d1=sf.en;
    snf.calc_mu(sf,T);
    d2=sf.en;
    dsdT=(d1-d2)/eps;

    sf.non_interacting=false;
    sf.nu=sf.mu;
    sf.ms=5.0+eps;
    snf.calc_mu(sf,T);
    d1=sf.n;
    
    sf.ms=5.0;
    snf.calc_mu(sf,T);
    d2=sf.n;
    dndm=(d1-d2)/eps;
    sf.non_interacting=true;

    snf.calc_mu(sf,T);
    cout << "sf: " << sf.dndmu << " " << sf.dndT << " " << sf.dsdT << endl;
    cout << "un: " << snf.unc.dndmu << " " << snf.unc.dndT << " " 
	 << snf.unc.dsdT << endl;
    cout << "nu: " << dndmu << " " << dndT << " " << dsdT << " "
	 << dndm << endl;
    double dndm2=3.0*sf.n/sf.ms-(sf.dndT*T/sf.ms+sf.nu/sf.ms*sf.dndmu);
    cout << dndm2 << endl;
    cout << endl;

    // Test the relation between specific heats
    double cp=snf.heat_cap_ppart_const_press(sf,T);
    double cv=snf.heat_cap_ppart_const_vol(sf,T);
    double alpha=snf.coeff_thermal_exp(sf,T);
    double beta=snf.compress_const_tptr(sf,T);
    t.test_rel(cp-cv,alpha*alpha/sf.n*T/beta,1.0e-10,"cp-cv");

  }

  t.report();

  return 0;
}
