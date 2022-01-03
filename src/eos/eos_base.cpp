/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

#include <o2scl/eos_base.h>

using namespace std;
using namespace o2scl;

eos_base::eos_base() {
  eos_thermo=&def_thermo;
}

void eos_base::set_thermo(thermo &th) {
  eos_thermo=&th;
}

const thermo &eos_base::get_thermo() {
  return *eos_thermo;
}

int eos_base::beta_eq_T0(ubvector &nB_grid, ubvector &guess,
			 fermion &e, bool include_muons,
			 fermion &mu, fermion_rel &frel,
			 std::shared_ptr<table_units<> > results) {
  O2SCL_ERR("Function beta_eq_T0() not implemented.",
	    o2scl::exc_eunimpl);
  return o2scl::exc_eunimpl;
}
