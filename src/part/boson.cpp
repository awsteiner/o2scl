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

#include <o2scl/boson.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

//--------------------------------------------
// boson class

boson::boson(double mass, double dof) : part(mass,dof) {
  co=0.0;
}

void boson::massless_calc(double temper) {
  
  n=g*pow(temper,3.0)*zeta3/pi2;
  ed=g*pi2*pow(temper,4.0)/30.0;
  pr=ed/3.0;
  en=g*pi2*pow(temper,3.0)/22.5;

  return;
}

