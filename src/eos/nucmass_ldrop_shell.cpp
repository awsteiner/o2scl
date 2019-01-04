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

#include <o2scl/nucmass_ldrop_shell.h>

using namespace std;
using namespace o2scl;

nucmass_ldrop_shell::nucmass_ldrop_shell() {
  nfit=11;
  inc_shell=true;
}

double nucmass_ldrop_shell::drip_binding_energy_d
(double Z, double N, double npout, double nnout, double chi, double T) {
 
  double ret=nucmass_ldrop_pair::drip_binding_energy_d
    (Z,N,npout,nnout,chi,T);

  if (inc_shell==false) return ret;

  shell=shell_energy_interp(Z,N);

  return ret-shell;
}

int nucmass_ldrop_shell::fit_fun(size_t nv, const ubvector &x) {
  doi=x[0];
  surften=x[1];
  ss=x[2];
  coul_coeff=x[3];
  n1=x[4];
  n0=x[5];
  Epair=x[6];
  s_a1=x[7];
  s_a2=x[8];
  s_a3=x[9];
  s_anp=x[10];
  return 0;
}
    
int nucmass_ldrop_shell::guess_fun(size_t nv, ubvector &x) {
  x[0]=doi;
  x[1]=surften;
  x[2]=ss;
  x[3]=coul_coeff;
  x[4]=n1;
  x[5]=n0;
  x[6]=Epair;
  x[7]=s_a1;
  x[8]=s_a2;
  x[9]=s_a3;
  x[10]=s_anp;
  return 0;
}

nucmass_frdm_shell::nucmass_frdm_shell() {
  nfit=14;
}

double nucmass_frdm_shell::mass_excess_d(double Z, double N) {
  
  double ret=nucmass_frdm::mass_excess_d(Z,N);
  ret-=shell_energy_interp(Z,N);
      
  return ret;
}

int nucmass_frdm_shell::fit_fun(size_t nv, const ubvector &x) {
  K=x[0]*200.0;
  r0=x[1];
  W=x[2];
  a1=x[3];
  a2=x[4];
  J=x[5];
  Q=x[6];
  ca=x[7];
  C=x[8]*60.0;
  gamma=x[9];

  s_a1=x[10];
  s_a2=x[11];
  s_a3=x[12];
  s_anp=x[13];
  return 0;
}
    
int nucmass_frdm_shell::guess_fun(size_t nv, ubvector &x) {
  x[0]=K/200.0;
  x[1]=r0;
  x[2]=W;
  x[3]=a1;
  x[4]=a2;
  x[5]=J;
  x[6]=Q;
  x[7]=ca;
  x[8]=C/60.0;
  x[9]=gamma;

  x[10]=s_a1;
  x[11]=s_a2;
  x[12]=s_a3;
  x[13]=s_anp;
  return 0;
}
