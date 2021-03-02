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

#include <o2scl/table.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>

extern "C" {

void *o2scl_hdf_create_hdf_file();

void o2scl_hdf_free_hdf_file(void *vptr);

void o2scl_hdf_hdf_file_open(void *vptr, char *fname, bool write_access, bool err_on_fail);

void o2scl_hdf_hdf_file_open_or_create(void *vptr, char *fname);

void o2scl_hdf_hdf_file_close(void *vptr);

void o2scl_hdf_hdf_input_table_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_n_table_wrapper(void *ptr_hf, void *ptr_t, void *&ptr_name);

void o2scl_hdf_hdf_output_table_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_table_units_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_n_table_units_wrapper(void *ptr_hf, void *ptr_t, void *&ptr_name);

void o2scl_hdf_hdf_output_table_units_wrapper(void *ptr_hf, void *ptr_t, char *name);

}
