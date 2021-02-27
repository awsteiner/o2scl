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

#include <o2scl/hdf_python.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

void *o2scl_hdf_create_hdf_file() {
  hdf_file *ptr=new hdf_file;
  return ptr;
}

void o2scl_hdf_free_hdf_file(void *vptr) {
  hdf_file *ptr=(hdf_file *)vptr;
  delete ptr;
  return;
}

 void o2scl_hdf_hdf_file_open(void *vptr, char *fname, bool write_access, bool err_on_fail) {
  hdf_file *ptr=(hdf_file *)vptr;
  ptr->open(fname,write_access,err_on_fail);
  return;
}

 void o2scl_hdf_hdf_file_open_or_create(void *vptr, char *fname) {
  hdf_file *ptr=(hdf_file *)vptr;
  ptr->open_or_create(fname);
  return;
}

 void o2scl_hdf_hdf_file_close(void *vptr) {
  hdf_file *ptr=(hdf_file *)vptr;
  ptr->close();
  return;
}

void o2scl_hdf_hdf_input_table_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table<> *t=(table<> *)ptr_t;
  hdf_input(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_output_table_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table<> *t=(table<> *)ptr_t;
  hdf_output(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_input_table_units_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table_units<> *t=(table_units<> *)ptr_t;
  hdf_input(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_output_table_units_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table_units<> *t=(table_units<> *)ptr_t;
  hdf_output(*hf,*t,name);
  return;
}

