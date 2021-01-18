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

#include <o2scl/eos_python.h>
#include <o2scl/part.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>

using namespace std;
using namespace o2scl;

void *o2scl_create_eos_base() {
  eos_base *ptr=new eos_base;
  return ptr;
}

void o2scl_free_eos_base(void *vptr) {
  eos_base *ptr=(eos_base *)vptr;
  delete ptr;
}

void o2scl_eos_base_get_def_thermo(void *vptr, void *p_v) {
  eos_base *ptr=(eos_base *)vptr;
  o2scl::thermo *p_t=(o2scl::thermo *)p_v;
  *(p_t)=ptr->def_thermo;
  return;
}

void o2scl_eos_base_set_def_thermo(void *vptr, void *p_v) {
  eos_base *ptr=(eos_base *)vptr;
  o2scl::thermo *p_t=(o2scl::thermo *)p_v;
  ptr->def_thermo=*(p_t);
  return;
}

double o2scl_eos_had_base_get_eoa(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->eoa;
}

void o2scl_eos_had_base_set_eoa(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->eoa=v;
  return;
}

double o2scl_eos_had_base_get_msom(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->msom;
}

void o2scl_eos_had_base_set_msom(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->msom=v;
  return;
}

double o2scl_eos_had_base_get_comp(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->comp;
}

void o2scl_eos_had_base_set_comp(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->comp=v;
  return;
}

double o2scl_eos_had_base_get_n0(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->n0;
}

void o2scl_eos_had_base_set_n0(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->n0=v;
  return;
}

double o2scl_eos_had_base_get_esym(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->esym;
}

void o2scl_eos_had_base_set_esym(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->esym=v;
  return;
}

double o2scl_eos_had_base_get_kprime(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->kprime;
}

void o2scl_eos_had_base_set_kprime(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->kprime=v;
  return;
}

bool o2scl_eos_had_base_get_err_nonconv(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->err_nonconv;
}

void o2scl_eos_had_base_set_err_nonconv(void *vptr, bool v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->err_nonconv=v;
  return;
}

void o2scl_eos_had_base_get_def_neutron(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_t=(o2scl::fermion *)p_v;
  *(p_t)=ptr->def_neutron;
  return;
}

void o2scl_eos_had_base_set_def_neutron(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_t=(o2scl::fermion *)p_v;
  ptr->def_neutron=*(p_t);
  return;
}

void o2scl_eos_had_base_get_def_proton(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_t=(o2scl::fermion *)p_v;
  *(p_t)=ptr->def_proton;
  return;
}

void o2scl_eos_had_base_set_def_proton(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_t=(o2scl::fermion *)p_v;
  ptr->def_proton=*(p_t);
  return;
}

int o2scl_eos_had_base_calc_e(void *vptr, void *ptr_n, void *ptr_p, void *ptr_th) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *n=(o2scl::fermion *)ptr_n;
  o2scl::fermion *p=(o2scl::fermion *)ptr_p;
  o2scl::thermo *th=(o2scl::thermo *)ptr_th;
  int ret=ptr->calc_e(*n,*p,*th);
  return ret;
}

void *o2scl_create_eos_had_skyrme() {
  eos_had_skyrme *ptr=new eos_had_skyrme;
  return ptr;
}

void o2scl_free_eos_had_skyrme(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  delete ptr;
}

double o2scl_eos_had_skyrme_get_t0(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->t0;
}

void o2scl_eos_had_skyrme_set_t0(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->t0=v;
  return;
}

double o2scl_eos_had_skyrme_get_t1(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->t1;
}

void o2scl_eos_had_skyrme_set_t1(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->t1=v;
  return;
}

double o2scl_eos_had_skyrme_get_t2(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->t2;
}

void o2scl_eos_had_skyrme_set_t2(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->t2=v;
  return;
}

double o2scl_eos_had_skyrme_get_t3(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->t3;
}

void o2scl_eos_had_skyrme_set_t3(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->t3=v;
  return;
}

double o2scl_eos_had_skyrme_get_x0(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->x0;
}

void o2scl_eos_had_skyrme_set_x0(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->x0=v;
  return;
}

double o2scl_eos_had_skyrme_get_x1(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->x1;
}

void o2scl_eos_had_skyrme_set_x1(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->x1=v;
  return;
}

double o2scl_eos_had_skyrme_get_x2(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->x2;
}

void o2scl_eos_had_skyrme_set_x2(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->x2=v;
  return;
}

double o2scl_eos_had_skyrme_get_x3(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->x3;
}

void o2scl_eos_had_skyrme_set_x3(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->x3=v;
  return;
}

double o2scl_eos_had_skyrme_get_alpha(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->alpha;
}

void o2scl_eos_had_skyrme_set_alpha(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->alpha=v;
  return;
}

double o2scl_eos_had_skyrme_get_a(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->a;
}

void o2scl_eos_had_skyrme_set_a(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->a=v;
  return;
}

double o2scl_eos_had_skyrme_get_b(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->b;
}

void o2scl_eos_had_skyrme_set_b(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->b=v;
  return;
}

double o2scl_eos_had_skyrme_get_W0(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->W0;
}

void o2scl_eos_had_skyrme_set_W0(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->W0=v;
  return;
}

double o2scl_eos_had_skyrme_get_b4(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->b4;
}

void o2scl_eos_had_skyrme_set_b4(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->b4=v;
  return;
}

double o2scl_eos_had_skyrme_get_b4p(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->b4p;
}

void o2scl_eos_had_skyrme_set_b4p(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->b4p=v;
  return;
}

bool o2scl_eos_had_skyrme_get_parent_method(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->parent_method;
}

void o2scl_eos_had_skyrme_set_parent_method(void *vptr, bool v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->parent_method=v;
  return;
}

const char *o2scl_eos_had_skyrme_get_reference(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  python_temp_string=ptr->reference;
  return python_temp_string.c_str();
}

void o2scl_eos_had_skyrme_set_reference(void *vptr, void *p_v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  std::string *p_t=(std::string *)p_v;
  ptr->reference=*(p_t);
  return;
}

void o2scl_eos_had_skyrme_get_nrfd(void *vptr, void *p_v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  o2scl::fermion_deriv_nr *p_t=(o2scl::fermion_deriv_nr *)p_v;
  *(p_t)=ptr->nrfd;
  return;
}

void o2scl_eos_had_skyrme_set_nrfd(void *vptr, void *p_v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  o2scl::fermion_deriv_nr *p_t=(o2scl::fermion_deriv_nr *)p_v;
  ptr->nrfd=*(p_t);
  return;
}

