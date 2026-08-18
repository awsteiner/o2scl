/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2020-2026, Andrew W. Steiner

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

#include <o2scl/mmin.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/min.h>
#include <o2scl/min_brent_gsl.h>

extern "C" {

int o2scl_mmin_base__get_verbose(void *vptr);

void o2scl_mmin_base__set_verbose(void *vptr, int v);

int o2scl_mmin_base__get_ntrial(void *vptr);

void o2scl_mmin_base__set_ntrial(void *vptr, int v);

double o2scl_mmin_base__get_tol_rel(void *vptr);

void o2scl_mmin_base__set_tol_rel(void *vptr, double v);

double o2scl_mmin_base__get_tol_abs(void *vptr);

void o2scl_mmin_base__set_tol_abs(void *vptr, double v);

int o2scl_mmin_base__get_last_ntrial(void *vptr);

void o2scl_mmin_base__set_last_ntrial(void *vptr, int v);

bool o2scl_mmin_base__get_err_nonconv(void *vptr);

void o2scl_mmin_base__set_err_nonconv(void *vptr, bool v);

void *o2scl_create_mmin_simp2_();

void o2scl_free_mmin_simp2_(void *vptr);

double o2scl_mmin_simp2__get_size(void *vptr);

void o2scl_mmin_simp2__set_size(void *vptr, double v);

double o2scl_mmin_simp2__get_fval(void *vptr);

void o2scl_mmin_simp2__set_fval(void *vptr, double v);

int o2scl_mmin_simp2__get_print_simplex(void *vptr);

void o2scl_mmin_simp2__set_print_simplex(void *vptr, int v);

int o2scl_min_base__get_verbose(void *vptr);

void o2scl_min_base__set_verbose(void *vptr, int v);

int o2scl_min_base__get_ntrial(void *vptr);

void o2scl_min_base__set_ntrial(void *vptr, int v);

double o2scl_min_base__get_tol_rel(void *vptr);

void o2scl_min_base__set_tol_rel(void *vptr, double v);

double o2scl_min_base__get_tol_abs(void *vptr);

void o2scl_min_base__set_tol_abs(void *vptr, double v);

int o2scl_min_base__get_last_ntrial(void *vptr);

void o2scl_min_base__set_last_ntrial(void *vptr, int v);

bool o2scl_min_base__get_err_nonconv(void *vptr);

void o2scl_min_base__set_err_nonconv(void *vptr, bool v);

void *o2scl_create_min_brent_gsl_();

void o2scl_free_min_brent_gsl_(void *vptr);

double o2scl_min_brent_gsl__get_x_minimum(void *vptr);

void o2scl_min_brent_gsl__set_x_minimum(void *vptr, double v);

double o2scl_min_brent_gsl__get_x_lower(void *vptr);

void o2scl_min_brent_gsl__set_x_lower(void *vptr, double v);

double o2scl_min_brent_gsl__get_x_upper(void *vptr);

void o2scl_min_brent_gsl__set_x_upper(void *vptr, double v);

double o2scl_min_brent_gsl__get_f_minimum(void *vptr);

void o2scl_min_brent_gsl__set_f_minimum(void *vptr, double v);

double o2scl_min_brent_gsl__get_f_lower(void *vptr);

void o2scl_min_brent_gsl__set_f_lower(void *vptr, double v);

double o2scl_min_brent_gsl__get_f_upper(void *vptr);

void o2scl_min_brent_gsl__set_f_upper(void *vptr, double v);

}
