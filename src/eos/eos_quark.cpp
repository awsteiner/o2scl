/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/eos_quark.h>

using namespace std;
using namespace o2scl;

eos_quark::eos_quark() {
}

int eos_quark::calc_p(quark &u, quark &d, quark &s, thermo &th) {
  O2SCL_ERR("Tried to run missing function eos_quark::calc_p().",
	    o2scl::exc_eunimpl);
  return o2scl::exc_eunimpl;
}

int eos_quark::calc_e(quark &u, quark &d, quark &s, thermo &th) {
  O2SCL_ERR("Tried to run missing function eos_quark::calc_e().",
		o2scl::exc_eunimpl);
  return o2scl::exc_eunimpl;
}

int eos_quark::calc_temp_p(quark &u, quark &d, quark &s,
			   double temper, thermo &th) {
  O2SCL_ERR("Tried to run missing function eos_quark::calc_temp_p().",
		o2scl::exc_eunimpl);
  return o2scl::exc_eunimpl;
}

int eos_quark::calc_temp_e(quark &u, quark &d, quark &s, 
			   double temper, thermo &th) {
  O2SCL_ERR("Tried to run missing function eos_quark::calc_temp_e().",
		o2scl::exc_eunimpl);
  return o2scl::exc_eunimpl;
}

