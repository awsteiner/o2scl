/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
#include <o2scl/fermion_mag_zerot.h>
#include <o2scl/test_mgr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  test_mgr t;
  t.set_output_level(1);

  cout.setf(ios::scientific);

  // In "Gaussian units":
  cout << pow(4.414e13*o2scl_const::ec_gauss_fm2,2.0)*hc_mev_fm/8.0/pi/
    o2scl_const::fine_structure << endl;

  fermion_mag_zerot mfz;
  fermion e(0.5,2.0);
  e.mu=1.0;
  e.m=0.5;
  mfz.calc_mu_zerot(e);
  //cout << e.n << " " << e.mu << " " << e.ed << " " << e.pr << endl;

  double mag_field=1.0e16;

  mfz.calc_mu_zerot_mag(e,-mag_field*o2scl_const::ec_gauss_fm2,0.0);
  cout << e.n << " " << e.mu << " " << e.ed << " " << e.pr << endl;
  cout << mfz.nmax_up << " " << mfz.nmax_dn << endl;

  e.mu-=e.m;
  e.inc_rest_mass=false;
  mfz.calc_mu_zerot_mag(e,-mag_field*o2scl_const::ec_gauss_fm2,0.0);
  cout << e.n << " " << e.mu << " " << e.ed+e.n*e.m << " " << e.pr << endl;
  cout << mfz.nmax_up << " " << mfz.nmax_dn << endl;
  e.inc_rest_mass=true;
  e.mu+=e.m;

  e.mu*=0.9;
  mfz.calc_density_zerot_mag(e,-mag_field*o2scl_const::ec_gauss_fm2,0.0);
  cout << e.n << " " << e.mu << " " << e.ed << " " << e.pr << endl;
  cout << endl;

  e.m=0.5;
  for(e.n=1.0;e.n>0.9e-10;e.n/=10.0) {
    mfz.calc_density_zerot_mag(e,-mag_field*o2scl_const::ec_gauss_fm2,0.0);
    cout << e.n << " " << e.mu << " " << e.ed << " " << e.pr << endl;
  }
  cout << endl;

  e.n=pow(0.011,1.5)*sqrt(2.0)/2.0/pi2*(1.0-0.0001);
  cout << e.n << endl;
  mfz.calc_density_zerot_mag(e,-mag_field*o2scl_const::ec_gauss_fm2,0.0);
  cout << mfz.nmax_up << " " << mfz.nmax_dn << endl;
  // This point doesn't work unless err_nonconv is set to false in
  // the fermion_mag_zerot class. I don't know why.
  e.n=pow(0.011,1.5)*sqrt(2.0)/2.0/pi2*(1.0+0.0001);
  cout << e.n << endl;
  mfz.calc_density_zerot_mag(e,-mag_field*o2scl_const::ec_gauss_fm2,0.0);
  cout << mfz.nmax_up << " " << mfz.nmax_dn << endl;
  
  cout << endl;

  e.inc_rest_mass=false;
  e.m=0.5;
  for(e.n=1.0;e.n>0.9e-8;e.n/=10.0) {
    mfz.calc_density_zerot_mag(e,-mag_field*o2scl_const::ec_gauss_fm2,0.0);
    cout << e.n << " " << e.mu << " " << e.ed+e.n*e.m << " " << e.pr << endl;
  }

  t.report();

  return 0;
}
