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

#include <o2scl/min_python.h>

using namespace std;
using namespace o2scl;

int o2scl_mmin_base__get_verbose(void *vptr) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  return ptr->verbose;
}

void o2scl_mmin_base__set_verbose(void *vptr, int v) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  ptr->verbose=v;
  return;
}

int o2scl_mmin_base__get_ntrial(void *vptr) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  return ptr->ntrial;
}

void o2scl_mmin_base__set_ntrial(void *vptr, int v) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  ptr->ntrial=v;
  return;
}

double o2scl_mmin_base__get_tol_rel(void *vptr) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  return ptr->tol_rel;
}

void o2scl_mmin_base__set_tol_rel(void *vptr, double v) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  ptr->tol_rel=v;
  return;
}

double o2scl_mmin_base__get_tol_abs(void *vptr) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  return ptr->tol_abs;
}

void o2scl_mmin_base__set_tol_abs(void *vptr, double v) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  ptr->tol_abs=v;
  return;
}

int o2scl_mmin_base__get_last_ntrial(void *vptr) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  return ptr->last_ntrial;
}

void o2scl_mmin_base__set_last_ntrial(void *vptr, int v) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  ptr->last_ntrial=v;
  return;
}

bool o2scl_mmin_base__get_err_nonconv(void *vptr) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  return ptr->err_nonconv;
}

void o2scl_mmin_base__set_err_nonconv(void *vptr, bool v) {
  mmin_base<> *ptr=(mmin_base<> *)vptr;
  ptr->err_nonconv=v;
  return;
}

void *o2scl_create_mmin_simp2_() {
  mmin_simp2<> *ptr=new mmin_simp2<>;
  return ptr;
}

void o2scl_free_mmin_simp2_(void *vptr) {
  mmin_simp2<> *ptr=(mmin_simp2<> *)vptr;
  delete ptr;
  return;
}

double o2scl_mmin_simp2__get_size(void *vptr) {
  mmin_simp2<> *ptr=(mmin_simp2<> *)vptr;
  return ptr->size;
}

void o2scl_mmin_simp2__set_size(void *vptr, double v) {
  mmin_simp2<> *ptr=(mmin_simp2<> *)vptr;
  ptr->size=v;
  return;
}

double o2scl_mmin_simp2__get_fval(void *vptr) {
  mmin_simp2<> *ptr=(mmin_simp2<> *)vptr;
  return ptr->fval;
}

void o2scl_mmin_simp2__set_fval(void *vptr, double v) {
  mmin_simp2<> *ptr=(mmin_simp2<> *)vptr;
  ptr->fval=v;
  return;
}

int o2scl_mmin_simp2__get_print_simplex(void *vptr) {
  mmin_simp2<> *ptr=(mmin_simp2<> *)vptr;
  return ptr->print_simplex;
}

void o2scl_mmin_simp2__set_print_simplex(void *vptr, int v) {
  mmin_simp2<> *ptr=(mmin_simp2<> *)vptr;
  ptr->print_simplex=v;
  return;
}

int o2scl_min_base__get_verbose(void *vptr) {
  min_base<> *ptr=(min_base<> *)vptr;
  return ptr->verbose;
}

void o2scl_min_base__set_verbose(void *vptr, int v) {
  min_base<> *ptr=(min_base<> *)vptr;
  ptr->verbose=v;
  return;
}

int o2scl_min_base__get_ntrial(void *vptr) {
  min_base<> *ptr=(min_base<> *)vptr;
  return ptr->ntrial;
}

void o2scl_min_base__set_ntrial(void *vptr, int v) {
  min_base<> *ptr=(min_base<> *)vptr;
  ptr->ntrial=v;
  return;
}

double o2scl_min_base__get_tol_rel(void *vptr) {
  min_base<> *ptr=(min_base<> *)vptr;
  return ptr->tol_rel;
}

void o2scl_min_base__set_tol_rel(void *vptr, double v) {
  min_base<> *ptr=(min_base<> *)vptr;
  ptr->tol_rel=v;
  return;
}

double o2scl_min_base__get_tol_abs(void *vptr) {
  min_base<> *ptr=(min_base<> *)vptr;
  return ptr->tol_abs;
}

void o2scl_min_base__set_tol_abs(void *vptr, double v) {
  min_base<> *ptr=(min_base<> *)vptr;
  ptr->tol_abs=v;
  return;
}

int o2scl_min_base__get_last_ntrial(void *vptr) {
  min_base<> *ptr=(min_base<> *)vptr;
  return ptr->last_ntrial;
}

void o2scl_min_base__set_last_ntrial(void *vptr, int v) {
  min_base<> *ptr=(min_base<> *)vptr;
  ptr->last_ntrial=v;
  return;
}

bool o2scl_min_base__get_err_nonconv(void *vptr) {
  min_base<> *ptr=(min_base<> *)vptr;
  return ptr->err_nonconv;
}

void o2scl_min_base__set_err_nonconv(void *vptr, bool v) {
  min_base<> *ptr=(min_base<> *)vptr;
  ptr->err_nonconv=v;
  return;
}

void *o2scl_create_min_brent_gsl_() {
  min_brent_gsl<> *ptr=new min_brent_gsl<>;
  return ptr;
}

void o2scl_free_min_brent_gsl_(void *vptr) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  delete ptr;
  return;
}

double o2scl_min_brent_gsl__get_x_minimum(void *vptr) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  return ptr->x_minimum;
}

void o2scl_min_brent_gsl__set_x_minimum(void *vptr, double v) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  ptr->x_minimum=v;
  return;
}

double o2scl_min_brent_gsl__get_x_lower(void *vptr) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  return ptr->x_lower;
}

void o2scl_min_brent_gsl__set_x_lower(void *vptr, double v) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  ptr->x_lower=v;
  return;
}

double o2scl_min_brent_gsl__get_x_upper(void *vptr) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  return ptr->x_upper;
}

void o2scl_min_brent_gsl__set_x_upper(void *vptr, double v) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  ptr->x_upper=v;
  return;
}

double o2scl_min_brent_gsl__get_f_minimum(void *vptr) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  return ptr->f_minimum;
}

void o2scl_min_brent_gsl__set_f_minimum(void *vptr, double v) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  ptr->f_minimum=v;
  return;
}

double o2scl_min_brent_gsl__get_f_lower(void *vptr) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  return ptr->f_lower;
}

void o2scl_min_brent_gsl__set_f_lower(void *vptr, double v) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  ptr->f_lower=v;
  return;
}

double o2scl_min_brent_gsl__get_f_upper(void *vptr) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  return ptr->f_upper;
}

void o2scl_min_brent_gsl__set_f_upper(void *vptr, double v) {
  min_brent_gsl<> *ptr=(min_brent_gsl<> *)vptr;
  ptr->f_upper=v;
  return;
}

