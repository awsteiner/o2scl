/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/schematic_eos.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

schematic_eos::schematic_eos() {
  kprime=0.0;
  n0=0.16;
  eoa=-16.0/hc_mev_fm;
  comp=230.0/hc_mev_fm;
  kpp=0.0;
  a=17.0/hc_mev_fm;
  b=13.0/hc_mev_fm;
  gamma=1.0;
}

int schematic_eos::calc_e(fermion &ln, fermion &lp, thermo &lth) {
  double xp, barn, sym, symp;
  
  barn=ln.n+lp.n;
  
  if (barn<=0.0) {
    xp=0.0;
    lth.ed=0.0;
    ln.mu=0.0;
    lp.mu=0.0;
    lth.pr=0.0;
    return 0;
  } else {
    xp=lp.n/barn;
  }
  
  sym=a*pow(barn/n0,2.0/3.0)+b*pow(barn/n0,gamma);
  
  // The derivative of the symmetry energy w.r.t. density
  symp=a*2.0/3.0*pow(barn/n0,-1.0/3.0)+b*gamma*pow(barn/n0,gamma-1.0);
  
  lth.ed=ln.m*ln.n+lp.m*lp.n+barn*
    (eoa+comp/18.0*pow((barn/n0-1.0),2.0)+
     sym*pow((1-2*xp),2.0)+kprime/162.0*pow(barn/n0-1.0,3.0)+
     kpp/1944.0*pow(barn/n0-1.0,4.0));
  ln.mu=lth.ed/barn+barn/n0*
    (comp/9.0*(barn/n0-1)+kprime/54.0*pow(barn/n0-1.0,2.0)+
     kpp/486.0*pow(barn/n0-1.0,3.0)+
     pow(1-2.0*xp,2.0)*symp)+xp*4.0*(1-2.0*xp)*sym;
  lp.mu=ln.mu-4.0*(1.0-2.0*xp)*sym;
  lth.pr=-lth.ed+ln.mu*ln.n+lp.mu*lp.n;

  return 0;
}

