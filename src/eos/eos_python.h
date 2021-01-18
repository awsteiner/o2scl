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
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/part.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>

extern "C" {

void *o2scl_create_eos_base();

void o2scl_free_eos_base(void *vp);

void o2scl_eos_base_get_def_thermo(void *vp, void *p_v);

void o2scl_eos_base_set_def_thermo(void *vp, void *p_v);

double o2scl_eos_had_base_get_eoa(void *vp);

void o2scl_eos_had_base_set_eoa(void *vp, double v);

double o2scl_eos_had_base_get_msom(void *vp);

void o2scl_eos_had_base_set_msom(void *vp, double v);

double o2scl_eos_had_base_get_comp(void *vp);

void o2scl_eos_had_base_set_comp(void *vp, double v);

double o2scl_eos_had_base_get_n0(void *vp);

void o2scl_eos_had_base_set_n0(void *vp, double v);

double o2scl_eos_had_base_get_esym(void *vp);

void o2scl_eos_had_base_set_esym(void *vp, double v);

double o2scl_eos_had_base_get_kprime(void *vp);

void o2scl_eos_had_base_set_kprime(void *vp, double v);

bool o2scl_eos_had_base_get_err_nonconv(void *vp);

void o2scl_eos_had_base_set_err_nonconv(void *vp, bool v);

void o2scl_eos_had_base_get_def_neutron(void *vp, void *p_v);

void o2scl_eos_had_base_set_def_neutron(void *vp, void *p_v);

void o2scl_eos_had_base_get_def_proton(void *vp, void *p_v);

void o2scl_eos_had_base_set_def_proton(void *vp, void *p_v);

int o2scl_eos_had_base_calc_e(void *vptr, void *ptr_n, void *ptr_p, void *ptr_th);

void *o2scl_create_eos_had_skyrme();

void o2scl_free_eos_had_skyrme(void *vp);

double o2scl_eos_had_skyrme_get_t0(void *vp);

void o2scl_eos_had_skyrme_set_t0(void *vp, double v);

double o2scl_eos_had_skyrme_get_t1(void *vp);

void o2scl_eos_had_skyrme_set_t1(void *vp, double v);

double o2scl_eos_had_skyrme_get_t2(void *vp);

void o2scl_eos_had_skyrme_set_t2(void *vp, double v);

double o2scl_eos_had_skyrme_get_t3(void *vp);

void o2scl_eos_had_skyrme_set_t3(void *vp, double v);

double o2scl_eos_had_skyrme_get_x0(void *vp);

void o2scl_eos_had_skyrme_set_x0(void *vp, double v);

double o2scl_eos_had_skyrme_get_x1(void *vp);

void o2scl_eos_had_skyrme_set_x1(void *vp, double v);

double o2scl_eos_had_skyrme_get_x2(void *vp);

void o2scl_eos_had_skyrme_set_x2(void *vp, double v);

double o2scl_eos_had_skyrme_get_x3(void *vp);

void o2scl_eos_had_skyrme_set_x3(void *vp, double v);

double o2scl_eos_had_skyrme_get_alpha(void *vp);

void o2scl_eos_had_skyrme_set_alpha(void *vp, double v);

double o2scl_eos_had_skyrme_get_a(void *vp);

void o2scl_eos_had_skyrme_set_a(void *vp, double v);

double o2scl_eos_had_skyrme_get_b(void *vp);

void o2scl_eos_had_skyrme_set_b(void *vp, double v);

double o2scl_eos_had_skyrme_get_W0(void *vp);

void o2scl_eos_had_skyrme_set_W0(void *vp, double v);

double o2scl_eos_had_skyrme_get_b4(void *vp);

void o2scl_eos_had_skyrme_set_b4(void *vp, double v);

double o2scl_eos_had_skyrme_get_b4p(void *vp);

void o2scl_eos_had_skyrme_set_b4p(void *vp, double v);

bool o2scl_eos_had_skyrme_get_parent_method(void *vp);

void o2scl_eos_had_skyrme_set_parent_method(void *vp, bool v);

void o2scl_eos_had_skyrme_get_reference(void *vp, void *p_v);

void o2scl_eos_had_skyrme_set_reference(void *vp, void *p_v);

void o2scl_eos_had_skyrme_get_nrfd(void *vp, void *p_v);

void o2scl_eos_had_skyrme_set_nrfd(void *vp, void *p_v);

}
