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

void o2scl_eos_base_set_thermo(void *vptr, void *ptr_th) {
  eos_base *ptr=(eos_base *)vptr;
  o2scl::thermo *th=(o2scl::thermo *)ptr_th;
  ptr->set_thermo(*th);
  return;
}

o2scl::thermo o2scl_eos_base_get_thermo(void *vptr) {
  eos_base *ptr=(eos_base *)vptr;
  o2scl::thermo ret=ptr->get_thermo();
  return ret;
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

