/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2020-2025, Andrew W. Steiner

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

#include <o2scl/table_units.h>
#include <o2scl/table3d.h>
#include <o2scl/tensor.h>
#include <o2scl/tensor_grid.h>
#include <o2scl/find_constants.h>
#include <o2scl/convert_units.h>
#include <o2scl/lib_settings.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <o2scl/format_float.h>
#include <o2scl/interp_krige.h>
#include <o2scl/cli.h>
#include <o2scl/funct.h>
#include <o2scl/cursesw.h>
#include <o2scl/funct_to_fp.h>

extern "C" {

void *o2scl_create_std_vector_double_();

void o2scl_free_std_vector_double_(void *vptr);

void o2scl_copy_std_vector_double_(void *vsrc, void *vdest);

void o2scl_std_vector_double__resize(void *vptr, size_t n);

size_t o2scl_std_vector_double__size(void *vptr);

double o2scl_std_vector_double__getitem(void *vptr, size_t n);

void o2scl_std_vector_double__setitem(void *vptr, size_t i, double val);

void o2scl_std_vector_double__push_back(void *vptr, double x);

void *o2scl_create_std_vector_int_();

void o2scl_free_std_vector_int_(void *vptr);

void o2scl_copy_std_vector_int_(void *vsrc, void *vdest);

void o2scl_std_vector_int__resize(void *vptr, size_t n);

size_t o2scl_std_vector_int__size(void *vptr);

int o2scl_std_vector_int__getitem(void *vptr, size_t n);

void o2scl_std_vector_int__setitem(void *vptr, size_t i, int val);

void *o2scl_create_std_vector_size_t_();

void o2scl_free_std_vector_size_t_(void *vptr);

void o2scl_copy_std_vector_size_t_(void *vsrc, void *vdest);

void o2scl_std_vector_size_t__resize(void *vptr, size_t n);

size_t o2scl_std_vector_size_t__size(void *vptr);

size_t o2scl_std_vector_size_t__getitem(void *vptr, size_t n);

void o2scl_std_vector_size_t__setitem(void *vptr, size_t i, size_t val);

void *o2scl_create_std_vector_std_string_();

void o2scl_free_std_vector_std_string_(void *vptr);

void o2scl_copy_std_vector_std_string_(void *vsrc, void *vdest);

void o2scl_std_vector_std_string__resize(void *vptr, size_t n);

size_t o2scl_std_vector_std_string__size(void *vptr);

void o2scl_std_vector_std_string__push_back(void *vptr, void *ptr_x);

void *o2scl_std_vector_std_string__getitem(void *vptr, size_t n);

void o2scl_std_vector_std_string__setitem(void *vptr, size_t i, std::string *val);

void *o2scl_create_boost_numeric_ublas_vector_double_();

void o2scl_free_boost_numeric_ublas_vector_double_(void *vptr);

void o2scl_copy_boost_numeric_ublas_vector_double_(void *vsrc, void *vdest);

size_t o2scl_boost_numeric_ublas_vector_double__size(void *vptr);

void o2scl_boost_numeric_ublas_vector_double__resize(void *vptr, size_t n);

double o2scl_boost_numeric_ublas_vector_double__getitem(void *vptr, size_t i);

void o2scl_boost_numeric_ublas_vector_double__setitem(void *vptr, size_t i, double val);

void *o2scl_create_boost_numeric_ublas_vector_int_();

void o2scl_free_boost_numeric_ublas_vector_int_(void *vptr);

void o2scl_copy_boost_numeric_ublas_vector_int_(void *vsrc, void *vdest);

size_t o2scl_boost_numeric_ublas_vector_int__size(void *vptr);

void o2scl_boost_numeric_ublas_vector_int__resize(void *vptr, size_t n);

int o2scl_boost_numeric_ublas_vector_int__getitem(void *vptr, size_t i);

void o2scl_boost_numeric_ublas_vector_int__setitem(void *vptr, size_t i, int val);

void *o2scl_create_boost_numeric_ublas_matrix_double_();

void o2scl_free_boost_numeric_ublas_matrix_double_(void *vptr);

void o2scl_copy_boost_numeric_ublas_matrix_double_(void *vsrc, void *vdest);

size_t o2scl_boost_numeric_ublas_matrix_double__size1(void *vptr);

size_t o2scl_boost_numeric_ublas_matrix_double__size2(void *vptr);

void o2scl_boost_numeric_ublas_matrix_double__resize(void *vptr, size_t m, size_t n);

double o2scl_boost_numeric_ublas_matrix_double__getitem(void *vptr, size_t m, size_t n);

void o2scl_boost_numeric_ublas_matrix_double__setitem(void *vptr, size_t i, size_t j, double val);

void *o2scl_create_boost_numeric_ublas_matrix_int_();

void o2scl_free_boost_numeric_ublas_matrix_int_(void *vptr);

void o2scl_copy_boost_numeric_ublas_matrix_int_(void *vsrc, void *vdest);

size_t o2scl_boost_numeric_ublas_matrix_int__size1(void *vptr);

size_t o2scl_boost_numeric_ublas_matrix_int__size2(void *vptr);

void o2scl_boost_numeric_ublas_matrix_int__resize(void *vptr, size_t m, size_t n);

int o2scl_boost_numeric_ublas_matrix_int__getitem(void *vptr, size_t m, size_t n);

void o2scl_boost_numeric_ublas_matrix_int__setitem(void *vptr, size_t i, size_t j, int val);

void *o2scl_create_std_vector_std_vector_double_();

void o2scl_free_std_vector_std_vector_double_(void *vptr);

void o2scl_copy_std_vector_std_vector_double_(void *vsrc, void *vdest);

void o2scl_std_vector_std_vector_double__resize(void *vptr, size_t n);

size_t o2scl_std_vector_std_vector_double__size(void *vptr);

void o2scl_std_vector_std_vector_double__getitem(void *vptr, size_t n, double **dptr, int *n_);

void o2scl_std_vector_std_vector_double__setitem(void *vptr, size_t i, void *valptr);

void *o2scl_create_std_vector_std_vector_std_string_();

void o2scl_free_std_vector_std_vector_std_string_(void *vptr);

void o2scl_copy_std_vector_std_vector_std_string_(void *vsrc, void *vdest);

void o2scl_std_vector_std_vector_std_string__resize(void *vptr, size_t n);

size_t o2scl_std_vector_std_vector_std_string__size(void *vptr);

void *o2scl_std_vector_std_vector_std_string__getitem(void *vptr, size_t n);

void o2scl_std_vector_std_vector_std_string__setitem(void *vptr, size_t i, std::vector<std::string> val);

void *o2scl_create_std_complex_double_();

void o2scl_free_std_complex_double_(void *vptr);

double o2scl_std_complex_double__real(void *vptr);

void o2scl_std_complex_double__real_set(void *vptr, double value);

double o2scl_std_complex_double__imag(void *vptr);

void o2scl_std_complex_double__imag_set(void *vptr, double value);

void *o2scl_std_complex_double__init(double re, double im);

void *o2scl_create_lib_settings_class();

void o2scl_free_lib_settings_class(void *vptr);

void *o2scl_lib_settings_class_get_data_dir(void *vptr);

int o2scl_lib_settings_class_set_data_dir(void *vptr, void *ptr_dir);

void *o2scl_lib_settings_class_get_doc_dir(void *vptr);

int o2scl_lib_settings_class_set_doc_dir(void *vptr, void *ptr_dir);

bool o2scl_lib_settings_class_openmp_support(void *vptr);

bool o2scl_lib_settings_class_readline_support(void *vptr);

bool o2scl_lib_settings_class_ncurses_support(void *vptr);

bool o2scl_lib_settings_class_armadillo_support(void *vptr);

bool o2scl_lib_settings_class_eigen_support(void *vptr);

bool o2scl_lib_settings_class_fftw_support(void *vptr);

bool o2scl_lib_settings_class_hdf5_compression_support(void *vptr);

void *o2scl_lib_settings_class_system_type(void *vptr);

bool o2scl_lib_settings_class_range_check(void *vptr);

void *o2scl_lib_settings_class_time_compiled(void *vptr);

void *o2scl_lib_settings_class_date_compiled(void *vptr);

void *o2scl_lib_settings_class_o2scl_version(void *vptr);

void o2scl_lib_settings_class_config_h_report(void *vptr);

void *o2scl_lib_settings_class_get_convert_units(void *vptr);

void *o2scl_create_table_();

void o2scl_free_table_(void *vptr);

void o2scl_copy_table_(void *vsrc, void *vdest);

void o2scl_table__getitem(void *vptr, void *ptr_col, double **dptr, int *n_);

void o2scl_table__set(void *vptr, void *ptr_col, size_t row, double val);

double o2scl_table__get(void *vptr, void *ptr_col, size_t row);

size_t o2scl_table__get_ncolumns(void *vptr);

size_t o2scl_table__get_nlines(void *vptr);

void o2scl_table__set_nlines(void *vptr, size_t lines);

size_t o2scl_table__get_maxlines(void *vptr);

void o2scl_table__set_maxlines(void *vptr, size_t llines);

void o2scl_table__set_nlines_auto(void *vptr, size_t il);

void o2scl_table__inc_maxlines(void *vptr, size_t llines);

void o2scl_table__new_column(void *vptr, void *ptr_col);

void *o2scl_table__get_column_name(void *vptr, size_t icol);

void o2scl_table__rename_column(void *vptr, void *ptr_src, void *ptr_dest);

void o2scl_table__delete_column(void *vptr, void *ptr_col);

void *o2scl_table__get_sorted_name(void *vptr, size_t icol);

void o2scl_table__init_column(void *vptr, void *ptr_scol, double val);

bool o2scl_table__is_column(void *vptr, void *ptr_scol);

size_t o2scl_table__lookup_column(void *vptr, void *ptr_scol);

void o2scl_table__copy_column(void *vptr, void *ptr_src, void *ptr_dest);

void o2scl_table__add_col_from_table(void *vptr, void *ptr_source, void *ptr_src_index, void *ptr_src_col, void *ptr_dest_index, void *ptr_dest_col);

void o2scl_table__insert_table(void *vptr, void *ptr_source, void *ptr_src_index, bool allow_extrap, void *ptr_dest_index);

void o2scl_table__add_table(void *vptr, void *ptr_source);

void o2scl_table__new_row(void *vptr, size_t n);

void o2scl_table__copy_row(void *vptr, size_t src, size_t dest);

void o2scl_table__delete_row(void *vptr, void *ptr_scol, double val);

void o2scl_table__delete_rows_func(void *vptr, void *ptr_func);

void o2scl_table__delete_rows_ends(void *vptr, size_t row_start, size_t row_end);

void o2scl_table__line_of_names(void *vptr, void *ptr_names);

void o2scl_table__line_of_data(void *vptr, void *ptr_data);

void o2scl_table__insert_row(void *vptr, size_t nv, void *ptr_data, size_t row);

size_t o2scl_table__ordered_lookup(void *vptr, void *ptr_scol, double val);

size_t o2scl_table__lookup(void *vptr, void *ptr_scol, double val);

size_t o2scl_table__lookup_val(void *vptr, void *ptr_scol, double val, void *ptr_scol2);

void o2scl_table__set_interp_type(void *vptr, size_t interp_type);

size_t o2scl_table__get_interp_type(void *vptr);

double o2scl_table__interp(void *vptr, void *ptr_sx, double x0, void *ptr_sy);

double o2scl_table__interp_index(void *vptr, size_t ix, double x0, size_t iy);

void o2scl_table__deriv_col(void *vptr, void *ptr_x, void *ptr_y, void *ptr_yp);

double o2scl_table__deriv(void *vptr, void *ptr_sx, double x0, void *ptr_sy);

double o2scl_table__deriv_index(void *vptr, size_t ix, double x0, size_t iy);

void o2scl_table__deriv2_col(void *vptr, void *ptr_x, void *ptr_y, void *ptr_yp);

double o2scl_table__deriv2(void *vptr, void *ptr_sx, double x0, void *ptr_sy);

double o2scl_table__deriv2_index(void *vptr, size_t ix, double x0, size_t iy);

double o2scl_table__integ(void *vptr, void *ptr_sx, double x1, double x2, void *ptr_sy);

double o2scl_table__integ_index(void *vptr, size_t ix, double x1, double x2, size_t iy);

void o2scl_table__integ_col(void *vptr, void *ptr_x, void *ptr_y, void *ptr_yi);

double o2scl_table__max(void *vptr, void *ptr_max);

double o2scl_table__min(void *vptr, void *ptr_min);

void o2scl_table__zero_table(void *vptr);

void o2scl_table__clear(void *vptr);

void o2scl_table__clear_data(void *vptr);

void o2scl_table__clear_table(void *vptr);

void o2scl_table__clear_constants(void *vptr);

void o2scl_table__sort_table(void *vptr, void *ptr_scol);

void o2scl_table__sort_column(void *vptr, void *ptr_scol);

void o2scl_table__average_col_roll(void *vptr, void *ptr_col_name, size_t window);

void o2scl_table__average_rows(void *vptr, size_t window, bool rolling);

void o2scl_table__is_valid(void *vptr);

void o2scl_table__functions_columns(void *vptr, void *ptr_list);

void o2scl_table__function_column(void *vptr, void *ptr_function, void *ptr_scol);

double o2scl_table__row_function(void *vptr, void *ptr_scol, size_t row);

size_t o2scl_table__function_find_row(void *vptr, void *ptr_function);

void o2scl_table__summary(void *vptr);

void *o2scl_create_table_units_();

void o2scl_free_table_units_(void *vptr);

void o2scl_copy_table_units_(void *vsrc, void *vdest);

void o2scl_table_units__set_unit(void *vptr, void *ptr_col, void *ptr_unit);

void *o2scl_table_units__get_unit(void *vptr, void *ptr_col);

void o2scl_table_units__line_of_units(void *vptr, void *ptr_unit_line);

void o2scl_table_units__remove_unit(void *vptr, void *ptr_col);

int o2scl_table_units__convert_to_unit(void *vptr, void *ptr_col, void *ptr_unit, bool err_on_fail=true);

void *o2scl_create_uniform_grid_();

void o2scl_free_uniform_grid_(void *vptr);

size_t o2scl_uniform_grid__get_nbins(void *vptr);

size_t o2scl_uniform_grid__get_npoints(void *vptr);

bool o2scl_uniform_grid__is_log(void *vptr);

double o2scl_uniform_grid__get_start(void *vptr);

double o2scl_uniform_grid__get_end(void *vptr);

double o2scl_uniform_grid__get_width(void *vptr);

double o2scl_uniform_grid__getitem(void *vptr, size_t n);

void o2scl_uniform_grid__vector(void *vptr, void *ptr_v);

void *o2scl_create_uniform_grid_end_();

void o2scl_free_uniform_grid_end_(void *vptr);

void *o2scl_uniform_grid_end__init(double start, double end, size_t n_bins);

void *o2scl_create_uniform_grid_width_();

void o2scl_free_uniform_grid_width_(void *vptr);

void *o2scl_uniform_grid_width__init(double start, double width, size_t n_bins);

void *o2scl_create_uniform_grid_end_width_();

void o2scl_free_uniform_grid_end_width_(void *vptr);

void *o2scl_uniform_grid_end_width__init(double start, double end, double width);

void *o2scl_create_uniform_grid_log_end_();

void o2scl_free_uniform_grid_log_end_(void *vptr);

void *o2scl_uniform_grid_log_end__init(double start, double end, size_t n_bins);

void *o2scl_create_uniform_grid_log_width_();

void o2scl_free_uniform_grid_log_width_(void *vptr);

void *o2scl_uniform_grid_log_width__init(double start, double width, size_t n_bins);

void *o2scl_create_uniform_grid_log_end_width_();

void o2scl_free_uniform_grid_log_end_width_(void *vptr);

void *o2scl_uniform_grid_log_end_width__init(double start, double end, double width);

void *o2scl_create_table3d();

void o2scl_free_table3d(void *vptr);

void o2scl_copy_table3d(void *vsrc, void *vdest);

void o2scl_table3d_set_size(void *vptr, size_t nx, size_t ny);

void o2scl_table3d_set_xy(void *vptr, void *ptr_x_name, size_t nx, void *ptr_x, void *ptr_y_name, size_t ny, void *ptr_y);

void o2scl_table3d_set_xy_grid(void *vptr, void *ptr_x_name, void *ptr_x_grid, void *ptr_y_name, void *ptr_y_grid);

void o2scl_table3d_set_interp_type(void *vptr, size_t interp_type);

void o2scl_table3d_line_of_names(void *vptr, void *ptr_names);

void o2scl_table3d_set(void *vptr, size_t ix, size_t iy, void *ptr_name, double val);

double o2scl_table3d_get(void *vptr, size_t ix, size_t iy, void *ptr_name);

double o2scl_table3d_get_i(void *vptr, size_t ix, size_t iy, size_t iz);

void o2scl_table3d_set_i(void *vptr, size_t ix, size_t iy, size_t iz, double val);

void o2scl_table3d_set_val(void *vptr, double x, double y, void *ptr_name, double val);

double o2scl_table3d_get_val(void *vptr, double x, double y, void *ptr_name);

void o2scl_table3d_set_grid_x(void *vptr, size_t ix, double val);

void o2scl_table3d_set_grid_y(void *vptr, size_t iy, double val);

double o2scl_table3d_get_grid_x(void *vptr, size_t ix);

double o2scl_table3d_get_grid_y(void *vptr, size_t iy);

void *o2scl_table3d_get_x_name(void *vptr);

void *o2scl_table3d_get_y_name(void *vptr);

void o2scl_table3d_set_x_name(void *vptr, void *ptr_name);

void o2scl_table3d_set_y_name(void *vptr, void *ptr_name);

void o2scl_table3d_get_size(void *vptr, size_t *nx, size_t *ny);

size_t o2scl_table3d_get_nx(void *vptr);

size_t o2scl_table3d_get_ny(void *vptr);

size_t o2scl_table3d_get_nslices(void *vptr);

bool o2scl_table3d_is_size_set(void *vptr);

bool o2scl_table3d_is_xy_set(void *vptr);

void *o2scl_table3d_get_slice_name(void *vptr, size_t i);

void o2scl_table3d_new_slice(void *vptr, void *ptr_slice);

void o2scl_table3d_set_slice_all(void *vptr, void *ptr_name, double val);

size_t o2scl_table3d_lookup_slice(void *vptr, void *ptr_name);

bool o2scl_table3d_is_slice(void *vptr, void *ptr_name, size_t *ix);

void o2scl_table3d_rename_slice(void *vptr, void *ptr_name1, void *ptr_name2);

void o2scl_table3d_copy_slice(void *vptr, void *ptr_name1, void *ptr_name2);

void *o2scl_table3d_get_slice(void *vptr, void *ptr_slice);

void *o2scl_table3d_get_slice_i(void *vptr, void *ptr_slice);

void o2scl_table3d_lookup_x(void *vptr, double val, size_t ix);

void o2scl_table3d_lookup_y(void *vptr, double val, size_t iy);

double o2scl_table3d_interp(void *vptr, double x, double y, void *ptr_name);

double o2scl_table3d_deriv_x(void *vptr, double x, double y, void *ptr_name);

double o2scl_table3d_deriv_y(void *vptr, double x, double y, void *ptr_name);

double o2scl_table3d_deriv_xy(void *vptr, double x, double y, void *ptr_name);

double o2scl_table3d_integ_x(void *vptr, double x1, double x2, double y, void *ptr_name);

double o2scl_table3d_integ_y(void *vptr, double x, double y1, double y2, void *ptr_name);

void o2scl_table3d_zero_table(void *vptr);

void o2scl_table3d_clear(void *vptr);

int o2scl_table3d_function_matrix(void *vptr, void *ptr_function, void *ptr_mat, bool throw_on_err);

void o2scl_table3d_function_slice(void *vptr, void *ptr_function, void *ptr_slice);

void o2scl_table3d_summary(void *vptr);

void *o2scl_create_index_spec();

void o2scl_free_index_spec(void *vptr);

size_t o2scl_index_spec_get_type(void *vptr);

void o2scl_index_spec_set_type(void *vptr, size_t v);

size_t o2scl_index_spec_get_ix1(void *vptr);

void o2scl_index_spec_set_ix1(void *vptr, size_t v);

size_t o2scl_index_spec_get_ix2(void *vptr);

void o2scl_index_spec_set_ix2(void *vptr, size_t v);

size_t o2scl_index_spec_get_ix3(void *vptr);

void o2scl_index_spec_set_ix3(void *vptr, size_t v);

double o2scl_index_spec_get_val1(void *vptr);

void o2scl_index_spec_set_val1(void *vptr, double v);

double o2scl_index_spec_get_val2(void *vptr);

void o2scl_index_spec_set_val2(void *vptr, double v);

double o2scl_index_spec_get_val3(void *vptr);

void o2scl_index_spec_set_val3(void *vptr, double v);

void o2scl_free_ix_index(void *vptr);

void *o2scl_ix_index_init(size_t ix);

void o2scl_free_ix_fixed(void *vptr);

void *o2scl_ix_fixed_init(size_t ix, size_t ix2);

void o2scl_free_ix_sum(void *vptr);

void *o2scl_ix_sum_init(size_t ix);

void o2scl_free_ix_trace(void *vptr);

void *o2scl_ix_trace_init(size_t ix, size_t ix2);

void o2scl_free_ix_reverse(void *vptr);

void *o2scl_ix_reverse_init(size_t ix);

void o2scl_free_ix_range(void *vptr);

void *o2scl_ix_range_init(size_t ix, size_t start, size_t end);

void o2scl_free_ix_interp(void *vptr);

void *o2scl_ix_interp_init(size_t ix, double v);

void o2scl_free_ix_grid(void *vptr);

void *o2scl_ix_grid_init(size_t ix, double start, double end, size_t n_bins, bool log);

void o2scl_free_ix_gridw(void *vptr);

void *o2scl_ix_gridw_init(size_t ix, double start, double end, double width, bool log);

void *o2scl_create_tensor_();

void o2scl_free_tensor_(void *vptr);

void o2scl_copy_tensor_(void *vsrc, void *vdest);

void o2scl_tensor__is_valid(void *vptr);

void o2scl_tensor__clear(void *vptr);

void o2scl_tensor__set(void *vptr, void *ptr_index, double val);

void o2scl_tensor__set_all(void *vptr, double x);

void o2scl_tensor__swap_data(void *vptr, void *ptr_data);

double o2scl_tensor__get(void *vptr, void *ptr_index);

void o2scl_tensor__resize(void *vptr, size_t n, void *ptr_index);

size_t o2scl_tensor__get_rank(void *vptr);

size_t o2scl_tensor__get_size(void *vptr, size_t i);

void *o2scl_tensor__get_size_arr(void *vptr);

void o2scl_tensor__get_data(void *vptr, double **dptr, int *n_);

size_t o2scl_tensor__total_size(void *vptr);

size_t o2scl_tensor__pack_indices(void *vptr, void *ptr_index);

void o2scl_tensor__unpack_index(void *vptr, size_t ix, void *ptr_index);

double o2scl_tensor__min_value(void *vptr);

void o2scl_tensor__min_index(void *vptr, void *ptr_index);

void o2scl_tensor__min(void *vptr, void *ptr_ix, double *value);

double o2scl_tensor__max_value(void *vptr);

void o2scl_tensor__max_index(void *vptr, void *ptr_index);

void o2scl_tensor__max(void *vptr, void *ptr_ix, double *value);

void o2scl_tensor__minmax_value(void *vptr, double *min, double *max);

void o2scl_tensor__minmax_index(void *vptr, void *ptr_min, void *ptr_max);

void o2scl_tensor__minmax(void *vptr, void *ptr_min_ix, double *min_value, void *ptr_max_ix, double *max_value);

double o2scl_tensor__total_sum(void *vptr);

void o2scl_tensor__copy_table3d_sum(void *vptr, size_t ix_x, size_t ix_y, void *ptr_tab, void *ptr_x_name, void *ptr_y_name, void *ptr_slice_name);

void o2scl_tensor__copy_table3d(void *vptr, size_t ix_x, size_t ix_y, void *ptr_tab, void *ptr_x_name, void *ptr_y_name, void *ptr_slice_name);

void *o2scl_tensor__create_size(size_t rank, void *ptr_sizes);

void *o2scl_create_tensor_int_std_vector_int_();

void o2scl_free_tensor_int_std_vector_int_(void *vptr);

void o2scl_copy_tensor_int_std_vector_int_(void *vsrc, void *vdest);

void o2scl_tensor_int_std_vector_int__is_valid(void *vptr);

void o2scl_tensor_int_std_vector_int__clear(void *vptr);

void o2scl_tensor_int_std_vector_int__set(void *vptr, void *ptr_index, int val);

void o2scl_tensor_int_std_vector_int__set_all(void *vptr, int x);

int o2scl_tensor_int_std_vector_int__get(void *vptr, void *ptr_index);

void o2scl_tensor_int_std_vector_int__resize(void *vptr, size_t n, void *ptr_index);

size_t o2scl_tensor_int_std_vector_int__get_rank(void *vptr);

size_t o2scl_tensor_int_std_vector_int__get_size(void *vptr, size_t i);

void *o2scl_tensor_int_std_vector_int__get_data(void *vptr);

size_t o2scl_tensor_int_std_vector_int__total_size(void *vptr);

size_t o2scl_tensor_int_std_vector_int__pack_indices(void *vptr, void *ptr_index);

void o2scl_tensor_int_std_vector_int__unpack_index(void *vptr, size_t ix, void *ptr_index);

int o2scl_tensor_int_std_vector_int__min_value(void *vptr);

void o2scl_tensor_int_std_vector_int__min_index(void *vptr, void *ptr_index);

void o2scl_tensor_int_std_vector_int__min(void *vptr, void *ptr_index, int *val);

int o2scl_tensor_int_std_vector_int__max_value(void *vptr);

void o2scl_tensor_int_std_vector_int__max_index(void *vptr, void *ptr_index);

void o2scl_tensor_int_std_vector_int__max(void *vptr, void *ptr_index, int *val);

void o2scl_tensor_int_std_vector_int__minmax_value(void *vptr, int *min, int *max);

void o2scl_tensor_int_std_vector_int__minmax_index(void *vptr, void *ptr_index_min, void *ptr_index_max);

void o2scl_tensor_int_std_vector_int__minmax(void *vptr, void *ptr_index_min, int *min, void *ptr_index_max, int *max);

int o2scl_tensor_int_std_vector_int__total_sum(void *vptr);

void o2scl_tensor_int_std_vector_int__copy_table3d_sum(void *vptr, size_t ix_x, size_t ix_y, void *ptr_tab, void *ptr_x_name, void *ptr_y_name, void *ptr_slice_name);

void o2scl_tensor_int_std_vector_int__copy_table3d(void *vptr, size_t ix_x, size_t ix_y, void *ptr_tab, void *ptr_x_name, void *ptr_y_name, void *ptr_slice_name);

void *o2scl_tensor_int_std_vector_int__create_size(size_t rank, void *ptr_sizes);

void *o2scl_create_tensor_size_t_std_vector_size_t_();

void o2scl_free_tensor_size_t_std_vector_size_t_(void *vptr);

void o2scl_copy_tensor_size_t_std_vector_size_t_(void *vsrc, void *vdest);

void o2scl_tensor_size_t_std_vector_size_t__is_valid(void *vptr);

void o2scl_tensor_size_t_std_vector_size_t__clear(void *vptr);

void o2scl_tensor_size_t_std_vector_size_t__set(void *vptr, void *ptr_index, size_t val);

void o2scl_tensor_size_t_std_vector_size_t__set_all(void *vptr, size_t x);

int o2scl_tensor_size_t_std_vector_size_t__get(void *vptr, void *ptr_index);

void o2scl_tensor_size_t_std_vector_size_t__resize(void *vptr, size_t n, void *ptr_index);

size_t o2scl_tensor_size_t_std_vector_size_t__get_rank(void *vptr);

size_t o2scl_tensor_size_t_std_vector_size_t__get_size(void *vptr, size_t i);

void *o2scl_tensor_size_t_std_vector_size_t__get_data(void *vptr);

size_t o2scl_tensor_size_t_std_vector_size_t__total_size(void *vptr);

size_t o2scl_tensor_size_t_std_vector_size_t__min_value(void *vptr);

void o2scl_tensor_size_t_std_vector_size_t__min_index(void *vptr, void *ptr_index);

void o2scl_tensor_size_t_std_vector_size_t__min(void *vptr, void *ptr_index, size_t *val);

size_t o2scl_tensor_size_t_std_vector_size_t__max_value(void *vptr);

void o2scl_tensor_size_t_std_vector_size_t__max_index(void *vptr, void *ptr_index);

void o2scl_tensor_size_t_std_vector_size_t__max(void *vptr, void *ptr_index, size_t *val);

void o2scl_tensor_size_t_std_vector_size_t__minmax_value(void *vptr, size_t *min, size_t *max);

void o2scl_tensor_size_t_std_vector_size_t__minmax_index(void *vptr, void *ptr_index_min, void *ptr_index_max);

void o2scl_tensor_size_t_std_vector_size_t__minmax(void *vptr, void *ptr_index_min, size_t *min, void *ptr_index_max, size_t *max);

size_t o2scl_tensor_size_t_std_vector_size_t__total_sum(void *vptr);

void o2scl_tensor_size_t_std_vector_size_t__copy_table3d_sum(void *vptr, size_t ix_x, size_t ix_y, void *ptr_tab, void *ptr_x_name, void *ptr_y_name, void *ptr_slice_name);

void o2scl_tensor_size_t_std_vector_size_t__copy_table3d(void *vptr, size_t ix_x, size_t ix_y, void *ptr_tab, void *ptr_x_name, void *ptr_y_name, void *ptr_slice_name);

void *o2scl_tensor_size_t_std_vector_size_t__create_size(size_t rank, void *ptr_sizes);

void *o2scl_create_tensor_grid_();

void o2scl_free_tensor_grid_(void *vptr);

void o2scl_copy_tensor_grid_(void *vsrc, void *vdest);

void o2scl_tensor_grid__is_valid(void *vptr);

void o2scl_tensor_grid__set_val(void *vptr, void *ptr_grid_point, double val);

double o2scl_tensor_grid__get_val(void *vptr, void *ptr_grid_point);

void o2scl_tensor_grid__resize(void *vptr, size_t rank, void *ptr_dim);

bool o2scl_tensor_grid__is_grid_set(void *vptr);

void o2scl_tensor_grid__set_grid_packed(void *vptr, void *ptr_grid);

void o2scl_tensor_grid__set_grid_vec_vec(void *vptr, void *ptr_grid_vecs);

void o2scl_tensor_grid__default_grid(void *vptr);

void o2scl_tensor_grid__set_grid_i_vec(void *vptr, size_t i, void *ptr_grid);

void o2scl_tensor_grid__set_grid_i_func(void *vptr, size_t ix, void *ptr_func);

double o2scl_tensor_grid__get_grid(void *vptr, size_t i, size_t j);

void o2scl_tensor_grid__get_grid_packed(void *vptr, double **dptr, int *n_);

void o2scl_tensor_grid__set_grid(void *vptr, size_t i, size_t j, double val);

size_t o2scl_tensor_grid__lookup_grid(void *vptr, size_t i, double val);

void *o2scl_tensor_grid__copy_slice_interp(void *vptr, void *ptr_ifix, void *ptr_vals);

void o2scl_tensor_grid__copy_table3d_align(void *vptr, size_t ix_x, size_t ix_y, void *ptr_index, void *ptr_tab, void *ptr_z_name);

void o2scl_tensor_grid__copy_table3d_align_setxy(void *vptr, size_t ix_x, size_t ix_y, void *ptr_index, void *ptr_tab, void *ptr_x_name, void *ptr_y_name, void *ptr_z_name);

void o2scl_tensor_grid__copy_table3d_interp(void *vptr, size_t ix_x, size_t ix_y, void *ptr_index, void *ptr_tab, void *ptr_slice_name);

void o2scl_tensor_grid__copy_table3d_interp_values(void *vptr, size_t ix_x, size_t ix_y, void *ptr_values, void *ptr_tab, void *ptr_slice_name, int verbose=0);

void o2scl_tensor_grid__copy_table3d_interp_values_setxy(void *vptr, size_t ix_x, size_t ix_y, void *ptr_values, void *ptr_tab, void *ptr_x_name, void *ptr_y_name, void *ptr_slice_name);

void o2scl_tensor_grid__clear(void *vptr);

void o2scl_tensor_grid__set_interp_type(void *vptr, size_t interp_type);

double o2scl_tensor_grid__interp_linear_partial(void *vptr, void *ptr_ix_to_interp, void *ptr_ix, void *ptr_val);

double o2scl_tensor_grid__interp_linear(void *vptr, void *ptr_v);

void o2scl_tensor_grid__from_table3d_fermi(void *vptr, void *ptr_t3d, void *ptr_slice, size_t n_points, double low=0.0, double high=0.0, double width=0.0);

void *o2scl_create_find_constants_const_entry();

void o2scl_free_find_constants_const_entry(void *vptr);

void *o2scl_find_constants_const_entry_get_names(void *vptr);

void o2scl_find_constants_const_entry_set_names(void *vptr, void *p_v);

void *o2scl_find_constants_const_entry_get_unit(void *vptr);

void o2scl_find_constants_const_entry_set_unit(void *vptr, void *p_v);

int o2scl_find_constants_const_entry_get_unit_flag(void *vptr);

void o2scl_find_constants_const_entry_set_unit_flag(void *vptr, int v);

double o2scl_find_constants_const_entry_get_val(void *vptr);

void o2scl_find_constants_const_entry_set_val(void *vptr, double v);

void *o2scl_find_constants_const_entry_get_source(void *vptr);

void o2scl_find_constants_const_entry_set_source(void *vptr, void *p_v);

int o2scl_find_constants_const_entry_get_m(void *vptr);

void o2scl_find_constants_const_entry_set_m(void *vptr, int v);

int o2scl_find_constants_const_entry_get_k(void *vptr);

void o2scl_find_constants_const_entry_set_k(void *vptr, int v);

int o2scl_find_constants_const_entry_get_s(void *vptr);

void o2scl_find_constants_const_entry_set_s(void *vptr, int v);

int o2scl_find_constants_const_entry_get_K(void *vptr);

void o2scl_find_constants_const_entry_set_K(void *vptr, int v);

int o2scl_find_constants_const_entry_get_A(void *vptr);

void o2scl_find_constants_const_entry_set_A(void *vptr, int v);

int o2scl_find_constants_const_entry_get_mol(void *vptr);

void o2scl_find_constants_const_entry_set_mol(void *vptr, int v);

int o2scl_find_constants_const_entry_get_cd(void *vptr);

void o2scl_find_constants_const_entry_set_cd(void *vptr, int v);

void *o2scl_create_find_constants_();

void o2scl_free_find_constants_(void *vptr);

void o2scl_find_constants__output_list_cout(void *vptr);

void o2scl_find_constants__add_constant(void *vptr, void *ptr_f, int verbose=0);

void o2scl_find_constants__del_constant(void *vptr, void *ptr_name, int verbose=0);

void *o2scl_create_convert_units_der_unit();

void o2scl_free_convert_units_der_unit(void *vptr);

void *o2scl_convert_units_der_unit_get_label(void *vptr);

void o2scl_convert_units_der_unit_set_label(void *vptr, void *p_v);

int o2scl_convert_units_der_unit_get_m(void *vptr);

void o2scl_convert_units_der_unit_set_m(void *vptr, int v);

int o2scl_convert_units_der_unit_get_k(void *vptr);

void o2scl_convert_units_der_unit_set_k(void *vptr, int v);

int o2scl_convert_units_der_unit_get_s(void *vptr);

void o2scl_convert_units_der_unit_set_s(void *vptr, int v);

int o2scl_convert_units_der_unit_get_K(void *vptr);

void o2scl_convert_units_der_unit_set_K(void *vptr, int v);

int o2scl_convert_units_der_unit_get_A(void *vptr);

void o2scl_convert_units_der_unit_set_A(void *vptr, int v);

int o2scl_convert_units_der_unit_get_mol(void *vptr);

void o2scl_convert_units_der_unit_set_mol(void *vptr, int v);

int o2scl_convert_units_der_unit_get_cd(void *vptr);

void o2scl_convert_units_der_unit_set_cd(void *vptr, int v);

double o2scl_convert_units_der_unit_get_val(void *vptr);

void o2scl_convert_units_der_unit_set_val(void *vptr, double v);

void *o2scl_convert_units_der_unit_get_name(void *vptr);

void o2scl_convert_units_der_unit_set_name(void *vptr, void *p_v);

void *o2scl_create_convert_units_();

void o2scl_free_convert_units_(void *vptr);

int o2scl_convert_units__get_verbose(void *vptr);

void o2scl_convert_units__set_verbose(void *vptr, int v);

bool o2scl_convert_units__get_err_on_fail(void *vptr);

void o2scl_convert_units__set_err_on_fail(void *vptr, bool v);

bool o2scl_convert_units__get_combine_two_conv(void *vptr);

void o2scl_convert_units__set_combine_two_conv(void *vptr, bool v);

double o2scl_convert_units__convert(void *vptr, void *ptr_frm, void *ptr_to, double val);

int o2scl_convert_units__convert_ret(void *vptr, void *ptr_frm, void *ptr_to, double val, double converted);

void o2scl_convert_units__del_unit(void *vptr, void *ptr_name);

void o2scl_convert_units__add_unit(void *vptr, void *ptr_d);

void o2scl_convert_units__set_natural_units(void *vptr, bool c_is_one=true, bool hbar_is_one=true, bool kb_is_one=true);

int o2scl_convert_units__is_in_cache(void *vptr, void *ptr_frm, void *ptr_to);

int o2scl_convert_units__remove_cache(void *vptr, void *ptr_frm, void *ptr_to);

void o2scl_convert_units__clear_cache(void *vptr);

void o2scl_convert_units__test_unique(void *vptr);

void o2scl_convert_units__print_cache(void *vptr);

void o2scl_convert_units__print_units_cout(void *vptr);

void o2scl_convert_units__find_print(void *vptr, void *ptr_name, void *ptr_unit, size_t prec, bool use_regex);

double o2scl_convert_units__find_unique(void *vptr, void *ptr_name, void *ptr_unit, bool use_regex=false);

void *o2scl_create_columnify();

void o2scl_free_columnify(void *vptr);

int o2scl_columnify_get_align_left(void *vptr);


int o2scl_columnify_get_align_right(void *vptr);


int o2scl_columnify_get_align_lmid(void *vptr);


int o2scl_columnify_get_align_rmid(void *vptr);


int o2scl_columnify_get_align_dp(void *vptr);


int o2scl_columnify_get_align_lnum(void *vptr);


void *o2scl_create_format_float();

void o2scl_free_format_float(void *vptr);

void o2scl_format_float_set_sig_figs(void *vptr, size_t sig_figs);

void o2scl_format_float_set_exp_limits(void *vptr, int min, int max);

void o2scl_format_float_set_pad_zeros(void *vptr, bool pad);

void o2scl_format_float_set_dec_point(void *vptr, void *ptr_dec_point);

void o2scl_format_float_set_exp_digits(void *vptr, int d);

void o2scl_format_float_html_mode(void *vptr);

void o2scl_format_float_latex_mode(void *vptr);

void o2scl_format_float_c_mode(void *vptr);

void *o2scl_format_float_convert(void *vptr, double x, bool debug=false);

void *o2scl_create_interp_vec_std_vector_double_();

void o2scl_free_interp_vec_std_vector_double_(void *vptr);

void o2scl_interp_vec_std_vector_double__set(void *vptr, size_t n, void *ptr_x, void *ptr_y, int interp_type);

void o2scl_interp_vec_std_vector_double__clear(void *vptr);

double o2scl_interp_vec_std_vector_double__eval(void *vptr, double x0);

double o2scl_interp_vec_std_vector_double__deriv(void *vptr, double x0);

double o2scl_interp_vec_std_vector_double__deriv2(void *vptr, double x0);

double o2scl_interp_vec_std_vector_double__integ(void *vptr, double x1, double x2);

void *o2scl_create_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise_();

void o2scl_free_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise_(void *vptr);

size_t o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__get_mode_loo_cv(void *vptr);


size_t o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__get_mode_loo_cv_bf(void *vptr);


size_t o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__get_mode_max_lml(void *vptr);


int o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__get_verbose(void *vptr);

void o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__set_verbose(void *vptr, int v);

size_t o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__get_mode(void *vptr);

void o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__set_mode(void *vptr, size_t v);

void o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__set(void *vptr, size_t size, void *ptr_x, void *ptr_y);

double o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__eval(void *vptr, double x0);

double o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__deriv(void *vptr, double x0);

double o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__deriv2(void *vptr, double x0);

double o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__sigma(void *vptr, double x0);

double o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__sample(void *vptr, double x0);

void o2scl_interp_krige_optim_std_vector_double_std_vector_double_covar_funct_rbf_noise__sample_vec(void *vptr, void *ptr_x, void *ptr_y);

void *o2scl_create_terminal();

void o2scl_free_terminal(void *vptr);

bool o2scl_terminal_is_redirected(void *vptr);

size_t o2scl_terminal_str_len(void *vptr, void *ptr_str);

void *o2scl_terminal_hrule(void *vptr, size_t n=78);

void *o2scl_terminal_cyan_fg(void *vptr);

void *o2scl_terminal_magenta_fg(void *vptr);

void *o2scl_terminal_yellow_fg(void *vptr);

void *o2scl_terminal_red_fg(void *vptr);

void *o2scl_terminal_green_fg(void *vptr);

void *o2scl_terminal_blue_fg(void *vptr);

void *o2scl_terminal_cyan_bg(void *vptr);

void *o2scl_terminal_magenta_bg(void *vptr);

void *o2scl_terminal_yellow_bg(void *vptr);

void *o2scl_terminal_red_bg(void *vptr);

void *o2scl_terminal_green_bg(void *vptr);

void *o2scl_terminal_blue_bg(void *vptr);

void *o2scl_terminal_default_fgbg(void *vptr);

void *o2scl_terminal_bold(void *vptr);

void *o2scl_terminal_eight_bit_fg(void *vptr, short col);

void *o2scl_terminal_eight_bit_bg(void *vptr, short col);

void *o2scl_terminal_lowint(void *vptr);

void *o2scl_terminal_underline(void *vptr);

void *o2scl_terminal_reverse(void *vptr);

void *o2scl_terminal_alt_font(void *vptr);

void *o2scl_terminal_normal_font(void *vptr);

void *o2scl_terminal_eight_bit_summ(void *vptr);

void *o2scl_terminal_three_byte_summ(void *vptr);

void *o2scl_terminal_three_byte_summ_long(void *vptr);

void *o2scl_create_gen_test_number_double_();

void o2scl_free_gen_test_number_double_(void *vptr);

void o2scl_gen_test_number_double__reset(void *vptr);

void o2scl_gen_test_number_double__set_radix(void *vptr, double r);

double o2scl_gen_test_number_double__gen(void *vptr);

void o2scl_free_funct_string_double_(void *vptr);

int o2scl_funct_string_double__set_parm(void *vptr, void *ptr_name, double val);

double o2scl_funct_string_double__getitem(void *vptr, double x);

void *o2scl_funct_string_double__init(char *expr, char *var);

void *o2scl_create_comm_option_s();

void o2scl_free_comm_option_s(void *vptr);

char o2scl_comm_option_s_get_shrt(void *vptr);

void o2scl_comm_option_s_set_shrt(void *vptr, char v);

void *o2scl_comm_option_s_get_lng(void *vptr);

void o2scl_comm_option_s_set_lng(void *vptr, void *p_v);

void *o2scl_comm_option_s_get_desc(void *vptr);

void o2scl_comm_option_s_set_desc(void *vptr, void *p_v);

int o2scl_comm_option_s_get_min_parms(void *vptr);

void o2scl_comm_option_s_set_min_parms(void *vptr, int v);

int o2scl_comm_option_s_get_max_parms(void *vptr);

void o2scl_comm_option_s_set_max_parms(void *vptr, int v);

void *o2scl_comm_option_s_get_parm_desc(void *vptr);

void o2scl_comm_option_s_set_parm_desc(void *vptr, void *p_v);

void *o2scl_comm_option_s_get_help(void *vptr);

void o2scl_comm_option_s_set_help(void *vptr, void *p_v);

int o2scl_comm_option_s_get_type(void *vptr);

void o2scl_comm_option_s_set_type(void *vptr, int v);

void *o2scl_create_cmd_line_arg();

void o2scl_free_cmd_line_arg(void *vptr);

void *o2scl_cmd_line_arg_get_arg(void *vptr);

void o2scl_cmd_line_arg_set_arg(void *vptr, void *p_v);

bool o2scl_cmd_line_arg_get_is_option(void *vptr);

void o2scl_cmd_line_arg_set_is_option(void *vptr, bool v);

bool o2scl_cmd_line_arg_get_is_valid(void *vptr);

void o2scl_cmd_line_arg_set_is_valid(void *vptr, bool v);

void *o2scl_cmd_line_arg_get_parms(void *vptr);

void o2scl_cmd_line_arg_set_parms(void *vptr, void *p_v);

void *o2scl_create_cli();

void o2scl_free_cli(void *vptr);

bool o2scl_cli_get_sync_verbose(void *vptr);

void o2scl_cli_set_sync_verbose(void *vptr, bool v);

bool o2scl_cli_get_gnu_intro(void *vptr);

void o2scl_cli_set_gnu_intro(void *vptr, bool v);

void *o2scl_cli_get_desc(void *vptr);

void o2scl_cli_set_desc(void *vptr, void *p_v);

void *o2scl_cli_get_cmd_name(void *vptr);

void o2scl_cli_set_cmd_name(void *vptr, void *p_v);

void *o2scl_cli_get_addl_help_cmd(void *vptr);

void o2scl_cli_set_addl_help_cmd(void *vptr, void *p_v);

void *o2scl_cli_get_addl_help_cli(void *vptr);

void o2scl_cli_set_addl_help_cli(void *vptr, void *p_v);

int o2scl_cli_set_verbose(void *vptr, int v);

int o2scl_cli_parse_for_aliases(void *vptr, void *ptr_sv, bool allow_undashed, bool debug=false);

int o2scl_cli_apply_aliases(void *vptr, void *ptr_sv, size_t istart, bool debug=false);

void *o2scl_cli_get_option_list(void *vptr);

void *o2scl_cli_parameter_desc(void *vptr, void *ptr_name);

void *o2scl_cli_option_short_desc(void *vptr, void *ptr_name);

void *o2scl_create_shared_ptr_table_units_();

void o2scl_free_shared_ptr_table_units_(void *vptr);

void *o2scl_shared_ptr_table_units__ptr(void *vp);

void *o2scl_rearrange_and_copy_tensor_double__wrapper(void *ptr_t, void *ptr_spec, int verbose=0, bool err_on_fail=true);

void *o2scl_rearrange_and_copy_tensor_int_int__wrapper(void *ptr_t, void *ptr_spec, int verbose=0, bool err_on_fail=true);

void *o2scl_rearrange_and_copy_tensor_size_t_size_t__wrapper(void *ptr_t, void *ptr_spec, int verbose=0, bool err_on_fail=true);

void *o2scl_grid_rearrange_and_copy_tensor_grid_double__wrapper(void *ptr_t, void *ptr_spec, int verbose=0, bool err_on_fail=true);

double o2scl_fermi_function_wrapper(double x);

double o2scl_bose_function_wrapper(double x);

double o2scl_quadratic_extremum_x_double__wrapper(double x1, double x2, double x3, double y1, double y2, double y3);

double o2scl_quadratic_extremum_y_double__wrapper(double x1, double x2, double x3, double y1, double y2, double y3);

void o2scl_screenify_vector_std_string__wrapper(size_t nin, void *ptr_in_cols, void *ptr_out_cols, size_t max_size=80);

bool o2scl_file_exists_wrapper(void *ptr_fname);

void o2scl_RGBtoHSV_wrapper(double r, double g, double b, void *ptr_h, void *ptr_s, void *ptr_v);

void o2scl_HSVtoRGB_wrapper(double h, double s, double v, void *ptr_r, void *ptr_g, void *ptr_b);

void o2scl_wordexp_single_file_wrapper(void *ptr_fname);

void o2scl_wordexp_wrapper_wrapper(void *ptr_word, void *ptr_matches);

double o2scl_function_to_double_wrapper(void *ptr_s, int verbose=0);

int o2scl_function_to_double_nothrow_wrapper(void *ptr_s, void *ptr_result, int verbose=0);

double o2scl_find_constant_wrapper(void *ptr_name, void *ptr_unit);

int o2scl_string_to_uint_list_vector_size_t__wrapper(void *ptr_x, void *ptr_list);

void o2scl_rewrap_keep_endlines_wrapper(void *ptr_str, void *ptr_sv, size_t ncol=79, int verbose=0, bool ignore_vt100=true);

size_t o2scl_vector_level_count_std_vector_double_std_vector_double__wrapper(double level, size_t n, void *ptr_x, void *ptr_y);

void o2scl_vector_deriv_interp_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_v, void *ptr_dv, size_t interp_type=2);

void o2scl_vector_deriv2_interp_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_v, void *ptr_dv, size_t interp_type=2);

void o2scl_vector_deriv_xy_interp_std_vector_double_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_vx, void *ptr_vy, void *ptr_dv, size_t interp_type=2);

void o2scl_vector_deriv2_xy_interp_std_vector_double_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_vx, void *ptr_vy, void *ptr_dv, size_t interp_type=2);

double o2scl_vector_integ_interp_std_vector_double__wrapper(size_t n, void *ptr_vx, size_t interp_type=2);

double o2scl_vector_integ_xy_interp_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_vx, void *ptr_vy, size_t interp_type=2);

double o2scl_vector_integ_ul_interp_std_vector_double__wrapper(size_t n, double x2, void *ptr_v, size_t interp_type=2);

double o2scl_vector_integ_ul_xy_interp_std_vector_double_std_vector_double__wrapper(size_t n, double x2, void *ptr_vx, void *ptr_vy, size_t interp_type=2);

void o2scl_vector_find_level_std_vector_double_std_vector_double__wrapper(double level, size_t n, void *ptr_x, void *ptr_y, void *ptr_locs);

void o2scl_vector_invert_enclosed_sum_std_vector_double_std_vector_double__wrapper(double sum, size_t n, void *ptr_x, void *ptr_y, void *ptr_lev, int boundaries=0, int verbose=0, bool err_on_fail=true);

int o2scl_vector_region_int_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_x, void *ptr_y, double intl, void *ptr_locs, int boundaries=0, int verbose=0, bool err_on_fail=true);

int o2scl_vector_region_fracint_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_x, void *ptr_y, double intl, void *ptr_locs, int boundaries=0, int verbose=0, bool err_on_fail=true);

int o2scl_vector_bound_fracint_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_x, void *ptr_y, double frac, void *ptr_low, void *ptr_high, int boundaries=0, int verbose=0, bool err_on_fail=true);

int o2scl_vector_bound_int_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_x, void *ptr_y, double frac, void *ptr_low, void *ptr_high, int boundaries=0, int verbose=0, bool err_on_fail=true);

void o2scl_rebin_xy_std_vector_double_std_vector_double_std_vector_double_std_vector_double__wrapper(void *ptr_x, void *ptr_y, void *ptr_x_out, void *ptr_y_out, size_t n_pts, size_t interp_type);

double o2scl_linear_or_log_chi2_std_vector_double_std_vector_double__wrapper(void *ptr_x, void *ptr_y);

void o2scl_linear_or_log_std_vector_double_std_vector_double__wrapper(void *ptr_x, void *ptr_y, void *ptr_log_x, void *ptr_log_y);

void o2scl_vector_refine_std_vector_double_std_vector_double_double__wrapper(size_t n, void *ptr_index, void *ptr_data, size_t factor, size_t interp_type=2);

void o2scl_linear_or_log_std_vector_double__wrapper(void *ptr_x, void *ptr_log_x);

int o2scl_get_screen_size_ioctl_wrapper(void *ptr_row, void *ptr_col);

}
