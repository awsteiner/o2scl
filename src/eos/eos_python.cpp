/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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

#include <o2scl/eos_python.h>
#include <o2scl/hdf_eos_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

void *o2scl_eos_had_strings(char *eos_type, char *model) {
  return eos_had_strings(eos_type,model);
}

void o2scl_eos_had_calc_density_zerot
(double nn, double np, double *ed, double *pr,
 double *mun, double *mup, double *nun, double *nup,
 double *msn, double *msp) {

  return;
}
