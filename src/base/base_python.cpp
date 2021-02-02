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

#include <o2scl/base_python.h>

using namespace std;
using namespace o2scl;

void *o2scl_create_table__() {
  table<> *ptr=new table<>;
  return ptr;
}

void o2scl_free_table__(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  delete ptr;
}

void o2scl_table___set(void *vptr, char *col, size_t row, double val) {
  table<> *ptr=(table<> *)vptr;
  ptr->set(col,row,val);
  return;
}

double o2scl_table___get(void *vptr, char *col, size_t row) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->get(col,row);
  return ret;
}

size_t o2scl_table___get_ncolumns(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->get_ncolumns();
  return ret;
}

size_t o2scl_table___get_nlines(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->get_nlines();
  return ret;
}

void o2scl_table___set_nlines(void *vptr, size_t lines) {
  table<> *ptr=(table<> *)vptr;
  ptr->set_nlines(lines);
  return;
}

void o2scl_table___new_column(void *vptr, char *col) {
  table<> *ptr=(table<> *)vptr;
  ptr->new_column(col);
  return;
}

const char *o2scl_table___get_column_name(void *vptr, size_t icol) {
  table<> *ptr=(table<> *)vptr;
  std::string ret=ptr->get_column_name(icol);
  python_temp_string=ret;
  return python_temp_string.c_str();
}

void o2scl_table___clear(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->clear();
  return;
}

void o2scl_table___clear_data(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->clear_data();
  return;
}

void o2scl_table___clear_table(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->clear_table();
  return;
}

void o2scl_table___clear_constants(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->clear_constants();
  return;
}

void o2scl_table___index_operator(void *vptr, char *col, double **dptr, int *n) {
  table<> *ptr=(table<> *)vptr;
  *dptr=(double *)(&(ptr->operator[](col)[0]));
  *n=ptr->get_nlines();
  return;
}

void *o2scl_create_table_units__() {
  table_units<> *ptr=new table_units<>;
  return ptr;
}

void o2scl_free_table_units__(void *vptr) {
  table_units<> *ptr=(table_units<> *)vptr;
  delete ptr;
}

const char *o2scl_table_units___get_unit(void *vptr, char *col) {
  table_units<> *ptr=(table_units<> *)vptr;
  std::string ret=ptr->get_unit(col);
  python_temp_string=ret;
  return python_temp_string.c_str();
}

void o2scl_table_units___set_unit(void *vptr, char *col, char *unit) {
  table_units<> *ptr=(table_units<> *)vptr;
  ptr->set_unit(col,unit);
  return;
}

int o2scl_table_units___convert_to_unit(void *vptr, char *col, char *unit, bool err_on_fail) {
  table_units<> *ptr=(table_units<> *)vptr;
  int ret=ptr->convert_to_unit(col,unit,err_on_fail);
  return ret;
}

void o2scl_table_units___clear_table(void *vptr) {
  table_units<> *ptr=(table_units<> *)vptr;
  ptr->clear_table();
  return;
}

void *o2scl_create_table3d() {
  table3d *ptr=new table3d;
  return ptr;
}

void o2scl_free_table3d(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  delete ptr;
}

void o2scl_table3d_set(void *vptr, size_t ix, size_t iy, char *name, double val) {
  table3d *ptr=(table3d *)vptr;
  ptr->set(ix,iy,name,val);
  return;
}

double o2scl_table3d_get(void *vptr, size_t ix, size_t iy, char *name) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->get(ix,iy,name);
  return ret;
}

void o2scl_table3d_new_slice(void *vptr, char *slice) {
  table3d *ptr=(table3d *)vptr;
  ptr->new_slice(slice);
  return;
}

size_t o2scl_table3d_get_nx(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  size_t ret=ptr->get_nx();
  return ret;
}

size_t o2scl_table3d_get_ny(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  size_t ret=ptr->get_ny();
  return ret;
}

size_t o2scl_table3d_get_nslices(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  size_t ret=ptr->get_nslices();
  return ret;
}

void *o2scl_create_tensor__() {
  tensor<> *ptr=new tensor<>;
  return ptr;
}

void o2scl_free_tensor__(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  delete ptr;
}

void o2scl_tensor___clear(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  ptr->clear();
  return;
}

void *o2scl_create_find_constants() {
  find_constants *ptr=new find_constants;
  return ptr;
}

void o2scl_free_find_constants(void *vptr) {
  find_constants *ptr=(find_constants *)vptr;
  delete ptr;
}

void o2scl_find_constants_find_print(void *vptr, char *name, char *unit, size_t prec, int verbose) {
  find_constants *ptr=(find_constants *)vptr;
  ptr->find_print(name,unit,prec,verbose);
  return;
}

double o2scl_find_constants_find_unique(void *vptr, char *name, char *unit) {
  find_constants *ptr=(find_constants *)vptr;
  double ret=ptr->find_unique(name,unit);
  return ret;
}

void *o2scl_create_convert_units__() {
  convert_units<> *ptr=new convert_units<>;
  return ptr;
}

void o2scl_free_convert_units__(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  delete ptr;
}

int o2scl_convert_units___get_verbose(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  return ptr->verbose;
}

void o2scl_convert_units___set_verbose(void *vptr, int v) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->verbose=v;
  return;
}

bool o2scl_convert_units___get_use_gnu_units(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  return ptr->use_gnu_units;
}

void o2scl_convert_units___set_use_gnu_units(void *vptr, bool v) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->use_gnu_units=v;
  return;
}

bool o2scl_convert_units___get_err_on_fail(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  return ptr->err_on_fail;
}

void o2scl_convert_units___set_err_on_fail(void *vptr, bool v) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->err_on_fail=v;
  return;
}

bool o2scl_convert_units___get_combine_two_conv(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  return ptr->combine_two_conv;
}

void o2scl_convert_units___set_combine_two_conv(void *vptr, bool v) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->combine_two_conv=v;
  return;
}

const char *o2scl_convert_units___get_units_cmd_string(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  python_temp_string=ptr->units_cmd_string;
  return python_temp_string.c_str();
}

void o2scl_convert_units___set_units_cmd_string(void *vptr, void *p_v) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  std::string *p_t=(std::string *)p_v;
  ptr->units_cmd_string=*(p_t);
  return;
}

double o2scl_convert_units___convert(void *vptr, char *frm, char *to, double val) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  double ret=ptr->convert(frm,to,val);
  return ret;
}

int o2scl_convert_units___convert_ret(void *vptr, char *frm, char *to, double val, double converted) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  int ret=ptr->convert_ret(frm,to,val,converted);
  return ret;
}

void o2scl_convert_units___print_cache(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->print_cache();
  return;
}

void *o2scl_create_shared_ptr_table_units__() {
  table_units<> *ptr=new table_units<>;
  return ptr;
}

void o2scl_free_shared_ptr_table_units__(void *vptr) {
  table_units<> *ptr=(table_units<> *)vptr;
  delete ptr;
}

void *o2scl_shared_ptr_table_units___ptr(void *vp) {
  std::shared_ptr<table_units<> > *p=(std::shared_ptr<table_units<> > *)vp;
  table_units<> *ref=p->get();
  return ref;
}

