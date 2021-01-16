/*
  -------------------------------------------------------------------

  Copyright (C) 2020-2021, Andrew W. Steiner

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

#include <o2scl/eos_base.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/part.h>

extern "C" {

void *o2scl_create_eos_base();

void o2scl_free_eos_base(void *vp);

void o2scl_eos_base_get_def_thermo(void *vp, void *p_v);

void o2scl_eos_base_set_def_thermo(void *vp, void *p_v);

void o2scl_eos_base_set_thermo(void *vptr, void *ptr_th);

o2scl::thermo o2scl_eos_base_get_thermo(void *vptr);

double o2scl_eos_had_base_get_eoa(void *vp);

void o2scl_eos_had_base_set_eoa(void *vp, double v);

}
