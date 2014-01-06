/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

#include <o2scl/sn_classical.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

//--------------------------------------------
// sn_classical class

sn_classical::sn_classical() {
}

sn_classical::~sn_classical() {
}

void sn_classical::calc_mu(part_deriv &p, double temper) {
  if (p.non_interacting==true) {p.nu=p.mu; p.ms=p.m;}
  if (p.nu/temper<-500.0) p.n=0.0;
  else p.n=exp(p.nu/temper)*p.g*pow(p.ms*temper/pi/2.0,1.5);
  p.ed=1.5*temper*p.n;
  p.pr=p.n*temper;
  p.en=(p.ed+p.pr-p.n*p.nu)/temper;
  p.dndT=-p.nu/temper/temper*p.n+1.5*p.n/temper;
  p.dndmu=p.n/temper;
  p.dsdT=2.5*p.dndT-p.nu*p.dndT/temper+p.n*p.nu/temper/temper;
  p.dndm=1.5*p.n/p.ms;
  return;
}

void sn_classical::calc_density(part_deriv &p, double temper) {
  if (p.non_interacting==true) {p.ms=p.m;}
  if (p.n<=0.0) {
    p.nu=p.ms;
    if (p.non_interacting==true) {p.mu=p.nu;}
    p.ed=0.0;
    p.pr=0.0;
    p.en=0.0;
    return;
  }
  p.nu=temper*log(p.n/p.g*pow(2.0*pi/p.ms/temper,1.5));
  if (p.non_interacting==true) {p.mu=p.nu;}
  p.ed=1.5*temper*p.n;
  p.pr=p.n*temper;
  p.en=(p.ed+p.pr-p.n*p.mu)/temper;
  p.dndT=-p.nu/temper/temper*p.n+1.5*p.n/temper;
  p.dndmu=p.n/temper;
  p.dsdT=2.5*p.dndT-p.nu*p.dndT/temper+p.n*p.nu/temper/temper;
  p.dndm=1.5*p.n/p.ms;
  return;
}

