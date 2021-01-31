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

#include <o2scl/table_units.h>
#include <o2scl/table3d.h>
#include <o2scl/tensor.h>
#include <o2scl/tensor_grid.h>
#include <o2scl/find_constants.h>
#include <o2scl/convert_units.h>

extern "C" {

void *o2scl_create_table__();

void o2scl_free_table__(void *vp);

void o2scl_table___set(void *vptr, char *col, size_t row, double val);

double o2scl_table___get(void *vptr, char *col, size_t row);

size_t o2scl_table___get_ncolumns(void *vptr);

size_t o2scl_table___get_nlines(void *vptr);

void o2scl_table___set_nlines(void *vptr, size_t lines);

void o2scl_table___new_column(void *vptr, char *col);

const char *o2scl_table___get_column_name(void *vptr, size_t icol);

void o2scl_table___clear(void *vptr);

void o2scl_table___clear_data(void *vptr);

void o2scl_table___clear_table(void *vptr);

void o2scl_table___clear_constants(void *vptr);

void *o2scl_create_table_units__();

void o2scl_free_table_units__(void *vp);

const char *o2scl_table_units___get_unit(void *vptr, char *col);

void o2scl_table_units___set_unit(void *vptr, char *col, char *unit);

int o2scl_table_units___convert_to_unit(void *vptr, char *col, char *unit, bool err_on_fail);

void o2scl_table_units___clear_table(void *vptr);

void *o2scl_create_table3d();

void o2scl_free_table3d(void *vp);

void o2scl_table3d_set(void *vptr, size_t ix, size_t iy, char *name, double val);

double o2scl_table3d_get(void *vptr, size_t ix, size_t iy, char *name);

void o2scl_table3d_new_slice(void *vptr, char *slice);

size_t o2scl_table3d_get_nx(void *vptr);

size_t o2scl_table3d_get_ny(void *vptr);

size_t o2scl_table3d_get_nslices(void *vptr);

void *o2scl_create_tensor__();

void o2scl_free_tensor__(void *vp);

void o2scl_tensor___clear(void *vptr);

void *o2scl_create_find_constants();

void o2scl_free_find_constants(void *vp);

void o2scl_find_constants_find_print(void *vptr, char *name, char *unit, size_t prec, int verbose);

double o2scl_find_constants_find_unique(void *vptr, char *name, char *unit);

void *o2scl_create_convert_units__();

void o2scl_free_convert_units__(void *vp);

int o2scl_convert_units___get_verbose(void *vp);

void o2scl_convert_units___set_verbose(void *vp, int v);

bool o2scl_convert_units___get_use_gnu_units(void *vp);

void o2scl_convert_units___set_use_gnu_units(void *vp, bool v);

bool o2scl_convert_units___get_err_on_fail(void *vp);

void o2scl_convert_units___set_err_on_fail(void *vp, bool v);

bool o2scl_convert_units___get_combine_two_conv(void *vp);

void o2scl_convert_units___set_combine_two_conv(void *vp, bool v);

void o2scl_convert_units___get_units_cmd_string(void *vp, void *p_v);

void o2scl_convert_units___set_units_cmd_string(void *vp, void *p_v);

double o2scl_convert_units___convert(void *vptr, char *frm, char *to, double val);

int o2scl_convert_units___convert_ret(void *vptr, char *frm, char *to, double val, double converted);

void o2scl_convert_units___print_cache(void *vptr);

void *o2scl_create_shared_ptr_table_units__();

void o2scl_free_shared_ptr_table_units__(void *vp);

void *o2scl_shared_ptr_table_units___ptr(void *vp);

}
