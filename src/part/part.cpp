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

#include <o2scl/part.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

thermo o2scl::operator+(const thermo &left, const thermo &right) {
  thermo th;
  th.ed=left.ed+right.ed;
  th.pr=left.pr+right.pr;
  th.en=left.en+right.en;
  return th;
}

thermo o2scl::operator+(const thermo &left, const part &right) {
  thermo th;
  th.ed=left.ed+right.ed;
  th.pr=left.pr+right.pr;
  th.en=left.en+right.en;
  return th;
}

thermo o2scl::operator-(const thermo &left, const thermo &right) {
  thermo th;
  th.ed=left.ed-right.ed;
  th.pr=left.pr-right.pr;
  th.en=left.en-right.en;
  return th;
}

thermo o2scl::operator-(const thermo &left, const part &right) {
  thermo th;
  th.ed=left.ed-right.ed;
  th.pr=left.pr-right.pr;
  th.en=left.en-right.en;
  return th;
}

//--------------------------------------------
// part class

part::part(double mass, double dof) {
  m=mass; 
  ms=mass; 
  g=dof;

  non_interacting=true;
  inc_rest_mass=true;
}

part::~part() {
}

void part::init(double mass, double dof) {
  m=mass; 
  ms=mass; 
  g=dof;
  return;
}

void part::anti(part& ax) {
  ax.g=g;
  ax.m=m;
  ax.ms=ms;
  ax.inc_rest_mass=inc_rest_mass;
  ax.non_interacting=non_interacting;
  if (inc_rest_mass) {
    ax.nu=-nu;
    ax.mu=-mu;
  } else {
    ax.nu=-nu-2.0*m;
    ax.mu=-mu-2.0*m;
  }
  return;
}

