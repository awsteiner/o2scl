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

#include <o2scl/quark_eos.h>

using namespace std;
using namespace o2scl;

quark_eos::quark_eos() {
  fet=&def_fet;
}

int quark_eos::calc_p(quark &u, quark &d, quark &s, thermo &th) {
  O2SCL_ERR("Tried to run missing function quark_eos::calc_p().",-1);
  return -1;
}

int quark_eos::calc_e(quark &u, quark &d, quark &s, thermo &th) {
  O2SCL_ERR("Tried to run missing function quark_eos::calc_e().",-1);
  return -1;
}

int quark_eos::calc_temp_p(quark &u, quark &d, quark &s,
			  double temper, thermo &th) {
  O2SCL_ERR("Tried to run missing function quark_eos::calc_temp_p().",-1);
  return -1;
}

int quark_eos::calc_temp_e(quark &u, quark &d, quark &s, 
			  double temper, thermo &th) {
  O2SCL_ERR("Tried to run missing function quark_eos::calc_temp_e().",-1);
  return -1;
}

