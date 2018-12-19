/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

#include <o2scl/classical.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

//--------------------------------------------
// classical_eval_thermo class

classical_eval_thermo::classical_eval_thermo() {
}

void classical_eval_thermo::calc_mu(part &p, double temper) {

  if (temper<0.0) {
    O2SCL_ERR2("Temperature less than zero in ",
	       "classical_eval_thermo::calc_mu().",exc_einval);
  }

  if (p.non_interacting==true) { p.nu=p.mu; p.ms=p.m; }

  // Handle zero temperature case
  if (temper==0.0) {
    if (p.inc_rest_mass) {
      p.n=0.0;
      p.ed=p.n*p.m;
    } else {
      p.n=0.0;
      p.ed=0.0;
    }
    p.pr=0.0;
    p.en=0.0;
    return;
  }

  if (p.inc_rest_mass) {
    if ((p.nu-p.m)/temper<-500.0) {
      p.n=0.0;
    } else {
      p.n=exp((p.nu-p.m)/temper)*p.g*pow(p.ms*temper/pi/2.0,1.5);
    }
    p.ed=1.5*temper*p.n+p.n*p.m;
  } else {
    if (p.nu/temper<-500.0) {
      p.n=0.0;
    } else {
      p.n=exp(p.nu/temper)*p.g*pow(p.ms*temper/pi/2.0,1.5);
    }
    p.ed=1.5*temper*p.n;
  }
  p.pr=p.n*temper;
  p.en=(p.ed+p.pr-p.n*p.nu)/temper;
  return;
}

void classical_eval_thermo::calc_density(part &p, double temper) {

  if (p.n<0.0 || temper<0.0) {
    O2SCL_ERR2("Density or temperature less than zero in ",
	       "classical_eval_thermo::calc_density().",exc_einval);
  }
    
  if (p.non_interacting==true) { p.ms=p.m; }

  // Handle zero density first
  if (p.n==0.0) {
    if (p.inc_rest_mass) {
      p.nu=p.m;
    } else {
      p.nu=0.0;
    }
    if (p.non_interacting==true) { p.mu=p.nu; }
    p.ed=0.0;
    p.pr=0.0;
    p.en=0.0;
    return;
  }

  // Handle the zero temperature case
  if (temper==0.0) {
    if (p.inc_rest_mass) {
      p.nu=p.m;
      p.ed=p.n*p.m;
    } else {
      p.nu=0.0;
      p.ed=0.0;
    }
    p.pr=0.0;
    p.en=0.0;
    return;
  }

  if (p.inc_rest_mass) {
    p.nu=p.m+temper*log(p.n/p.g*pow(2.0*pi/p.ms/temper,1.5));
    p.ed=1.5*temper*p.n+p.n*p.m;
  } else {
    p.nu=temper*log(p.n/p.g*pow(2.0*pi/p.ms/temper,1.5));
    p.ed=1.5*temper*p.n;
  }
  
  if (p.non_interacting==true) { p.mu=p.nu; }

  p.pr=p.n*temper;
  p.en=(p.ed+p.pr-p.n*p.nu)/temper;

  return;
}

