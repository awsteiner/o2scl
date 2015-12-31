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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/classical_deriv.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

//--------------------------------------------
// classical_deriv class

classical_deriv::classical_deriv() {
}

classical_deriv::~classical_deriv() {
}

void classical_deriv::calc_mu(part_deriv &p, double temper) {

  cl.calc_mu(p,temper);

  if (temper==0.0) {
    p.dndT=0.0;
    p.dndmu=0.0;
    p.dsdT=0.0;
    return;
  }

  p.dndT=-p.nu/temper/temper*p.n+1.5*p.n/temper;
  p.dndmu=p.n/temper;
  p.dsdT=2.5*p.dndT-p.nu*p.dndT/temper+p.n*p.nu/temper/temper;

  return;
}

void classical_deriv::calc_density(part_deriv &p, double temper) {

  cl.calc_density(p,temper);

  // Handle zero density first
  if (p.n==0.0 || temper==0.0) {
    p.dndT=0.0;
    p.dndmu=0.0;
    p.dsdT=0.0;
    return;
  }

  p.dndT=-p.nu/temper/temper*p.n+1.5*p.n/temper;
  p.dndmu=p.n/temper;
  p.dsdT=2.5*p.dndT-p.nu*p.dndT/temper+p.n*p.nu/temper/temper;
  
  return;
}

