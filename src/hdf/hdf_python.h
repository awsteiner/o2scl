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

int o2scl_hdf_hdf_file_get_compr_type(void *vptr);

void o2scl_hdf_hdf_file_set_compr_type(void *vptr, int v);

size_t o2scl_hdf_hdf_file_get_min_compr_size(void *vptr);

void o2scl_hdf_hdf_file_set_min_compr_size(void *vptr, size_t v);

bool o2scl_hdf_hdf_file_has_write_access(void *vptr);

void o2scl_hdf_hdf_file_open(void *vptr, char *fname, bool write_access=false, bool err_on_fail=true);

void o2scl_hdf_hdf_file_open_or_create(void *vptr, char *fname);

void o2scl_hdf_hdf_file_close(void *vptr);

int o2scl_hdf_hdf_file_getc(void *vptr, char *name, char *c);

int o2scl_hdf_hdf_file_getd(void *vptr, char *name, double *d);

int o2scl_hdf_hdf_file_geti(void *vptr, char *name, int *i);

int o2scl_hdf_hdf_file_get_szt(void *vptr, char *name, size_t *u);

int o2scl_hdf_hdf_file_gets(void *vptr, char *name, void *ptr_s);

int o2scl_hdf_hdf_file_gets_var(void *vptr, char *name, void *ptr_s);

int o2scl_hdf_hdf_file_gets_fixed(void *vptr, char *name, void *ptr_s);

void o2scl_hdf_hdf_file_setc(void *vptr, char *name, char c);

void o2scl_hdf_hdf_file_setd(void *vptr, char *name, double d);

void o2scl_hdf_hdf_file_seti(void *vptr, char *name, int i);

void o2scl_hdf_hdf_file_set_szt(void *vptr, char *name, size_t u);

void o2scl_hdf_hdf_file_sets(void *vptr, char *name, char *s);

void o2scl_hdf_hdf_file_sets_fixed(void *vptr, char *name, char *s);

int o2scl_hdf_hdf_file_getd_vec(void *vptr, char *name, void *ptr_v);

int o2scl_hdf_hdf_file_geti_vec(void *vptr, char *name, void *ptr_v);

int o2scl_hdf_hdf_file_get_szt_vec(void *vptr, char *name, void *ptr_v);

int o2scl_hdf_hdf_file_gets_vec(void *vptr, char *name, void *ptr_s);

int o2scl_hdf_hdf_file_setd_vec(void *vptr, char *name, void *ptr_v);

int o2scl_hdf_hdf_file_seti_vec(void *vptr, char *name, void *ptr_v);

int o2scl_hdf_hdf_file_set_szt_vec(void *vptr, char *name, void *ptr_v);

int o2scl_hdf_hdf_file_sets_vec(void *vptr, char *name, void *ptr_s);

int o2scl_hdf_hdf_file_getd_ten(void *vptr, char *name, void *ptr_t);

int o2scl_hdf_hdf_file_geti_ten(void *vptr, char *name, void *ptr_t);

int o2scl_hdf_hdf_file_get_szt_ten(void *vptr, char *name, void *ptr_t);

int o2scl_hdf_hdf_file_setd_ten(void *vptr, char *name, void *ptr_t);

int o2scl_hdf_hdf_file_seti_ten(void *vptr, char *name, void *ptr_t);

int o2scl_hdf_hdf_file_set_szt_ten(void *vptr, char *name, void *ptr_t);

int o2scl_hdf_hdf_file_find_object_by_type(void *vptr, char *type, void *ptr_name, int verbose=0);

int o2scl_hdf_hdf_file_find_object_by_name(void *vptr, char *name, void *ptr_type, int verbose=0);

int o2scl_hdf_hdf_file_find_object_by_pattern(void *vptr, char *pattern, void *ptr_type, int verbose=0);

void o2scl_hdf_hdf_file_file_list(void *vptr, int verbose);

void o2scl_hdf_hdf_file_copy(void *vptr, int verbose, void *ptr_hf2);

void o2scl_hdf_hdf_input_table_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_n_table_wrapper(void *ptr_hf, void *ptr_t, void *&ptr_name);

void o2scl_hdf_hdf_output_table_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_table_units_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_n_table_units_wrapper(void *ptr_hf, void *ptr_t, void *&ptr_name);

void o2scl_hdf_hdf_output_table_units_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_table3d_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_n_table3d_wrapper(void *ptr_hf, void *ptr_t, void *&ptr_name);

void o2scl_hdf_hdf_output_table3d_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_uniform_grid_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_n_uniform_grid_wrapper(void *ptr_hf, void *ptr_t, void *&ptr_name);

void o2scl_hdf_hdf_output_uniform_grid_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_tensor_grid_wrapper(void *ptr_hf, void *ptr_t, char *name);

void o2scl_hdf_hdf_input_n_tensor_grid_wrapper(void *ptr_hf, void *ptr_t, void *&ptr_name);

void o2scl_hdf_hdf_output_tensor_grid_wrapper(void *ptr_hf, void *ptr_t, char *name);

int o2scl_hdf_value_spec_wrapper(char *spec, void *ptr_d, int verbose=0, bool err_on_fail=true);

int o2scl_hdf_vector_spec_std_vector_double__wrapper(char *spec, void *ptr_v, int verbose=0, bool err_on_fail=true);

int o2scl_hdf_strings_spec_std_vector_std_string__wrapper(char *spec, void *ptr_v, int verbose=0, bool err_on_fail=true);

}
