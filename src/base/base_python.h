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
#include <o2scl/lib_settings.h>

extern "C" {

void *o2scl_create_std__string();

void o2scl_free_std__string(void *vptr);

 size_t o2scl_std__string_length(void *vptr);

 char o2scl_std__string_getitem(void *vptr, size_t n);

void o2scl_std__string_setitem(void *vptr, size_t n, char val);

 void o2scl_std__string_resize(void *vptr, size_t n);

const char *o2scl_std__string_c_str(void *vptr);

void *o2scl_create_std__vector_double_();

void o2scl_free_std__vector_double_(void *vptr);

 void o2scl_std__vector_double__resize(void *vptr, size_t n);

 size_t o2scl_std__vector_double__size(void *vptr);

 double o2scl_std__vector_double__getitem(void *vptr, size_t n);

void o2scl_std__vector_double__setitem(void *vptr, size_t n, double val);

void *o2scl_create_vector_int_();

void o2scl_free_vector_int_(void *vptr);

 void o2scl_vector_int__resize(void *vptr, size_t n);

 size_t o2scl_vector_int__size(void *vptr);

 int o2scl_vector_int__getitem(void *vptr, size_t n);

void o2scl_vector_int__setitem(void *vptr, size_t n, int val);

void *o2scl_create_vector_size_t_();

void o2scl_free_vector_size_t_(void *vptr);

 void o2scl_vector_size_t__resize(void *vptr, size_t n);

 size_t o2scl_vector_size_t__size(void *vptr);

 size_t o2scl_vector_size_t__getitem(void *vptr, size_t n);

void o2scl_vector_size_t__setitem(void *vptr, size_t n, size_t val);

void *o2scl_create_lib_settings_class();

void o2scl_free_lib_settings_class(void *vptr);

void *o2scl_lib_settings_class_get_data_dir(void *vptr);

 int o2scl_lib_settings_class_set_data_dir(void *vptr, char *dir);

void *o2scl_lib_settings_class_get_doc_dir(void *vptr);

 int o2scl_lib_settings_class_set_doc_dir(void *vptr, char *dir);

 bool o2scl_lib_settings_class_eos_installed(void *vptr);

 bool o2scl_lib_settings_class_part_installed(void *vptr);

 bool o2scl_lib_settings_class_hdf_support(void *vptr);

 bool o2scl_lib_settings_class_openmp_support(void *vptr);

 bool o2scl_lib_settings_class_readline_support(void *vptr);

 bool o2scl_lib_settings_class_ncurses_support(void *vptr);

 bool o2scl_lib_settings_class_gsl2_support(void *vptr);

 bool o2scl_lib_settings_class_armadillo_support(void *vptr);

 bool o2scl_lib_settings_class_eigen_support(void *vptr);

 bool o2scl_lib_settings_class_fftw_support(void *vptr);

 bool o2scl_lib_settings_class_python_support(void *vptr);

 bool o2scl_lib_settings_class_hdf5_compression_support(void *vptr);

void *o2scl_lib_settings_class_system_type(void *vptr);

 bool o2scl_lib_settings_class_range_check(void *vptr);

void *o2scl_lib_settings_class_time_compiled(void *vptr);

void *o2scl_lib_settings_class_date_compiled(void *vptr);

void *o2scl_lib_settings_class_o2scl_version(void *vptr);

 void o2scl_lib_settings_class_config_h_report(void *vptr);

void *o2scl_lib_settings_class_get_convert_units(void *vptr);

void *o2scl_create_table__();

void o2scl_free_table__(void *vptr);

void o2scl_copy_table__(void *vsrc, void *vdest);

void o2scl_table___getitem(void *vptr, char *col, double **dptr, int *n);

 void o2scl_table___set(void *vptr, char *col, size_t row, double val);

 double o2scl_table___get(void *vptr, char *col, size_t row);

 size_t o2scl_table___get_ncolumns(void *vptr);

 size_t o2scl_table___get_nlines(void *vptr);

 void o2scl_table___set_nlines(void *vptr, size_t lines);

 size_t o2scl_table___get_maxlines(void *vptr);

 void o2scl_table___set_maxlines(void *vptr, size_t llines);

 void o2scl_table___set_nlines_auto(void *vptr, size_t il);

 void o2scl_table___inc_maxlines(void *vptr, size_t llines);

 void o2scl_table___new_column(void *vptr, char *col);

void o2scl_table___get_column(void *vptr, char *col, double **dptr, int *n);

void *o2scl_table___get_column_name(void *vptr, size_t icol);

 void o2scl_table___rename_column(void *vptr, char *src, char *dest);

 void o2scl_table___delete_column(void *vptr, char *col);

void *o2scl_table___get_sorted_name(void *vptr, size_t icol);

 void o2scl_table___init_column(void *vptr, char *scol, double val);

 bool o2scl_table___is_column(void *vptr, char *scol);

 size_t o2scl_table___lookup_column(void *vptr, char *scol);

 void o2scl_table___copy_column(void *vptr, char *src, char *dest);

 void o2scl_table___add_col_from_table(void *vptr, void *ptr_source, char *src_index, char *src_col, char *dest_index, char *dest_col);

 void o2scl_table___insert_table(void *vptr, void *ptr_source, char *src_index, bool allow_extrap, char *dest_index);

 void o2scl_table___add_table(void *vptr, void *ptr_source);

 void o2scl_table___new_row(void *vptr, size_t n);

 void o2scl_table___copy_row(void *vptr, size_t src, size_t dest);

 void o2scl_table___delete_row(void *vptr, char *scol, double val);

 void o2scl_table___delete_rows_func(void *vptr, char *func);

 void o2scl_table___line_of_names(void *vptr, char *names);

 void o2scl_table___line_of_data(void *vptr, void *ptr_data);

 size_t o2scl_table___ordered_lookup(void *vptr, char *scol, double val);

 size_t o2scl_table___lookup(void *vptr, char *scol, double val);

 size_t o2scl_table___lookup_val(void *vptr, char *scol, double val, char *scol2);

 void o2scl_table___set_interp_type(void *vptr, size_t interp_type);

 size_t o2scl_table___get_interp_type(void *vptr);

 double o2scl_table___interp(void *vptr, char *sx, double x0, char *sy);

 double o2scl_table___interp_index(void *vptr, size_t ix, double x0, size_t iy);

 void o2scl_table___deriv_col(void *vptr, char *x, char *y, char *yp);

 double o2scl_table___deriv(void *vptr, char *sx, double x0, char *sy);

 double o2scl_table___deriv_index(void *vptr, size_t ix, double x0, size_t iy);

 void o2scl_table___deriv2_col(void *vptr, char *x, char *y, char *yp);

 double o2scl_table___deriv2(void *vptr, char *sx, double x0, char *sy);

 double o2scl_table___deriv2_index(void *vptr, size_t ix, double x0, size_t iy);

 double o2scl_table___integ(void *vptr, char *sx, double x1, double x2, char *sy);

 double o2scl_table___integ_index(void *vptr, size_t ix, double x1, double x2, size_t iy);

 void o2scl_table___integ_col(void *vptr, char *x, char *y, char *yi);

 double o2scl_table___max(void *vptr, char *max);

 double o2scl_table___min(void *vptr, char *min);

 void o2scl_table___zero_table(void *vptr);

 void o2scl_table___clear(void *vptr);

 void o2scl_table___clear_data(void *vptr);

 void o2scl_table___clear_table(void *vptr);

 void o2scl_table___clear_constants(void *vptr);

 void o2scl_table___sort_table(void *vptr, char *scol);

 void o2scl_table___sort_column(void *vptr, char *scol);

 void o2scl_table___average_col_roll(void *vptr, char *col_name, size_t window);

 void o2scl_table___average_rows(void *vptr, size_t window, bool rolling);

 void o2scl_table___is_valid(void *vptr);

 void o2scl_table___functions_columns(void *vptr, char *list);

 void o2scl_table___function_column(void *vptr, char *function, char *scol);

 double o2scl_table___row_function(void *vptr, char *scol, size_t row);

 size_t o2scl_table___function_find_row(void *vptr, char *function);

 void o2scl_table___summary(void *vptr);

void *o2scl_create_table_units__();

void o2scl_free_table_units__(void *vptr);

void o2scl_copy_table_units__(void *vsrc, void *vdest);

 void o2scl_table_units___set_unit(void *vptr, char *col, char *unit);

void *o2scl_table_units___get_unit(void *vptr, char *col);

 void o2scl_table_units___line_of_units(void *vptr, char *unit_line);

 void o2scl_table_units___remove_unit(void *vptr, char *col);

 int o2scl_table_units___convert_to_unit(void *vptr, char *col, char *unit, bool err_on_fail);

void *o2scl_create_uniform_grid__();

void o2scl_free_uniform_grid__(void *vptr);

 size_t o2scl_uniform_grid___get_nbins(void *vptr);

 size_t o2scl_uniform_grid___get_npoints(void *vptr);

 bool o2scl_uniform_grid___is_log(void *vptr);

 double o2scl_uniform_grid___get_start(void *vptr);

 double o2scl_uniform_grid___get_end(void *vptr);

 double o2scl_uniform_grid___get_width(void *vptr);

 double o2scl_uniform_grid___getitem(void *vptr, size_t n);

void o2scl_free_uniform_grid_end__(void *vptr);

void *o2scl_uniform_grid_end___init(double start, double end, size_t n_bins);

void o2scl_free_uniform_grid_width__(void *vptr);

void *o2scl_uniform_grid_width___init(double start, double width, size_t n_bins);

void o2scl_free_uniform_grid_end_width__(void *vptr);

void *o2scl_uniform_grid_end_width___init(double start, double end, double width);

void o2scl_free_uniform_grid_log_end__(void *vptr);

void *o2scl_uniform_grid_log_end___init(double start, double end, size_t n_bins);

void o2scl_free_uniform_grid_log_width__(void *vptr);

void *o2scl_uniform_grid_log_width___init(double start, double width, size_t n_bins);

void o2scl_free_uniform_grid_log_end_width__(void *vptr);

void *o2scl_uniform_grid_log_end_width___init(double start, double end, double width);

void *o2scl_create_table3d();

void o2scl_free_table3d(void *vptr);

void o2scl_copy_table3d(void *vsrc, void *vdest);

 void o2scl_table3d_set_xy(void *vptr, char *x_name, size_t nx, void *ptr_x, char *y_name, size_t ny, void *ptr_y);

 void o2scl_table3d_set(void *vptr, size_t ix, size_t iy, char *name, double val);

 double o2scl_table3d_get(void *vptr, size_t ix, size_t iy, char *name);

 void o2scl_table3d_new_slice(void *vptr, char *slice);

 size_t o2scl_table3d_get_nx(void *vptr);

 size_t o2scl_table3d_get_ny(void *vptr);

 size_t o2scl_table3d_get_nslices(void *vptr);

void *o2scl_create_tensor__();

void o2scl_free_tensor__(void *vptr);

void o2scl_copy_tensor__(void *vsrc, void *vdest);

 void o2scl_tensor___clear(void *vptr);

void *o2scl_create_find_constants();

void o2scl_free_find_constants(void *vptr);

 void o2scl_find_constants_find_print(void *vptr, char *name, char *unit, size_t prec, int verbose);

 double o2scl_find_constants_find_unique(void *vptr, char *name, char *unit);

void *o2scl_create_convert_units__();

void o2scl_free_convert_units__(void *vptr);

int o2scl_convert_units___get_verbose(void *vptr);

void o2scl_convert_units___set_verbose(void *vptr, int v);

bool o2scl_convert_units___get_use_gnu_units(void *vptr);

void o2scl_convert_units___set_use_gnu_units(void *vptr, bool v);

bool o2scl_convert_units___get_err_on_fail(void *vptr);

void o2scl_convert_units___set_err_on_fail(void *vptr, bool v);

bool o2scl_convert_units___get_combine_two_conv(void *vptr);

void o2scl_convert_units___set_combine_two_conv(void *vptr, bool v);

void *o2scl_convert_units___get_units_cmd_string(void *vptr);

void o2scl_convert_units___set_units_cmd_string(void *vptr, void *p_v);

 double o2scl_convert_units___convert(void *vptr, char *frm, char *to, double val);

 int o2scl_convert_units___convert_ret(void *vptr, char *frm, char *to, double val, double converted);

 void o2scl_convert_units___print_cache(void *vptr);

void *o2scl_create_shared_ptr_table_units__();

void o2scl_free_shared_ptr_table_units__(void *vptr);

void *o2scl_shared_ptr_table_units___ptr(void *vp);

}
