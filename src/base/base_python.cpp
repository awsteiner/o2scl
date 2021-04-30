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

void *o2scl_create_std__string() {
  std::string *ptr=new std::string;
  return ptr;
}

void o2scl_free_std__string(void *vptr) {
  std::string *ptr=(std::string *)vptr;
  delete ptr;
  return;
}

size_t o2scl_std__string_length(void *vptr) {
  std::string *ptr=(std::string *)vptr;
  size_t ret=ptr->length();
  return ret;
}

char o2scl_std__string_getitem(void *vptr, size_t n) {
  std::string *ptr=(std::string *)vptr;
  char ret=ptr->operator[](n);
  return ret;
}

void o2scl_std__string_setitem(void *vptr, size_t i, char val) {
  std::string *ptr=(std::string *)vptr;
  (*ptr)[i]=val;
  return;
}

void o2scl_std__string_resize(void *vptr, size_t n) {
  std::string *ptr=(std::string *)vptr;
  ptr->resize(n);
  return;
}

void *o2scl_create_std__vector_double_() {
  std::vector<double> *ptr=new std::vector<double>;
  return ptr;
}

void o2scl_free_std__vector_double_(void *vptr) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  delete ptr;
  return;
}

void o2scl_std__vector_double__resize(void *vptr, size_t n) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std__vector_double__size(void *vptr) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

double o2scl_std__vector_double__getitem(void *vptr, size_t n) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  double ret=ptr->operator[](n);
  return ret;
}

void o2scl_std__vector_double__setitem(void *vptr, size_t i, double val) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_std__vector_int_() {
  std::vector<int> *ptr=new std::vector<int>;
  return ptr;
}

void o2scl_free_std__vector_int_(void *vptr) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  delete ptr;
  return;
}

void o2scl_std__vector_int__resize(void *vptr, size_t n) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std__vector_int__size(void *vptr) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

int o2scl_std__vector_int__getitem(void *vptr, size_t n) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  int ret=ptr->operator[](n);
  return ret;
}

void o2scl_std__vector_int__setitem(void *vptr, size_t i, int val) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_std__vector_size_t_() {
  std::vector<size_t> *ptr=new std::vector<size_t>;
  return ptr;
}

void o2scl_free_std__vector_size_t_(void *vptr) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  delete ptr;
  return;
}

void o2scl_std__vector_size_t__resize(void *vptr, size_t n) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std__vector_size_t__size(void *vptr) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

size_t o2scl_std__vector_size_t__getitem(void *vptr, size_t n) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  size_t ret=ptr->operator[](n);
  return ret;
}

void o2scl_std__vector_size_t__setitem(void *vptr, size_t i, size_t val) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_std__vector_std__string_() {
  std::vector<std::string> *ptr=new std::vector<std::string>;
  return ptr;
}

void o2scl_free_std__vector_std__string_(void *vptr) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  delete ptr;
  return;
}

void o2scl_std__vector_std__string__resize(void *vptr, size_t n) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std__vector_std__string__size(void *vptr) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

void *o2scl_std__vector_std__string__getitem(void *vptr, size_t n) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->operator[](n);
  return sptr;
}

void o2scl_std__vector_std__string__setitem(void *vptr, size_t i, char *val) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_boost__numeric__ublas__vector_double_() {
  boost::numeric::ublas::vector<double> *ptr=new boost::numeric::ublas::vector<double>;
  return ptr;
}

void o2scl_free_boost__numeric__ublas__vector_double_(void *vptr) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  delete ptr;
  return;
}

size_t o2scl_boost__numeric__ublas__vector_double__size(void *vptr) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

void o2scl_boost__numeric__ublas__vector_double__resize(void *vptr, size_t n) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  ptr->resize(n);
  return;
}

double o2scl_boost__numeric__ublas__vector_double__getitem(void *vptr, size_t i) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  double ret=ptr->operator[](i);
  return ret;
}

void o2scl_boost__numeric__ublas__vector_double__setitem(void *vptr, size_t i, double val) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_boost__numeric__ublas__matrix_double_() {
  boost::numeric::ublas::matrix<double> *ptr=new boost::numeric::ublas::matrix<double>;
  return ptr;
}

void o2scl_free_boost__numeric__ublas__matrix_double_(void *vptr) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  delete ptr;
  return;
}

size_t o2scl_boost__numeric__ublas__matrix_double__size1(void *vptr) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  size_t ret=ptr->size1();
  return ret;
}

size_t o2scl_boost__numeric__ublas__matrix_double__size2(void *vptr) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  size_t ret=ptr->size2();
  return ret;
}

void o2scl_boost__numeric__ublas__matrix_double__resize(void *vptr, size_t m, size_t n) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  ptr->resize(m,n);
  return;
}

double o2scl_boost__numeric__ublas__matrix_double__getitem(void *vptr, size_t m, size_t n) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  double ret=ptr->operator()(m,n);
  return ret;
}

void o2scl_boost__numeric__ublas__matrix_double__setitem(void *vptr, size_t i, size_t j, double val) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  (*ptr)(i,j)=val;
  return;
}

void *o2scl_create_lib_settings_class() {
  lib_settings_class *ptr=new lib_settings_class;
  return ptr;
}

void o2scl_free_lib_settings_class(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  delete ptr;
  return;
}

void *o2scl_lib_settings_class_get_data_dir(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->get_data_dir();
  return sptr;
}

int o2scl_lib_settings_class_set_data_dir(void *vptr, char *dir) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  int ret=ptr->set_data_dir(dir);
  return ret;
}

void *o2scl_lib_settings_class_get_doc_dir(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->get_doc_dir();
  return sptr;
}

int o2scl_lib_settings_class_set_doc_dir(void *vptr, char *dir) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  int ret=ptr->set_doc_dir(dir);
  return ret;
}

bool o2scl_lib_settings_class_eos_installed(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->eos_installed();
  return ret;
}

bool o2scl_lib_settings_class_part_installed(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->part_installed();
  return ret;
}

bool o2scl_lib_settings_class_hdf_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->hdf_support();
  return ret;
}

bool o2scl_lib_settings_class_openmp_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->openmp_support();
  return ret;
}

bool o2scl_lib_settings_class_readline_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->readline_support();
  return ret;
}

bool o2scl_lib_settings_class_ncurses_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->ncurses_support();
  return ret;
}

bool o2scl_lib_settings_class_gsl2_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->gsl2_support();
  return ret;
}

bool o2scl_lib_settings_class_armadillo_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->armadillo_support();
  return ret;
}

bool o2scl_lib_settings_class_eigen_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->eigen_support();
  return ret;
}

bool o2scl_lib_settings_class_fftw_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->fftw_support();
  return ret;
}

bool o2scl_lib_settings_class_python_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->python_support();
  return ret;
}

bool o2scl_lib_settings_class_hdf5_compression_support(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->hdf5_compression_support();
  return ret;
}

void *o2scl_lib_settings_class_system_type(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->system_type();
  return sptr;
}

bool o2scl_lib_settings_class_range_check(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  bool ret=ptr->range_check();
  return ret;
}

void *o2scl_lib_settings_class_time_compiled(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->time_compiled();
  return sptr;
}

void *o2scl_lib_settings_class_date_compiled(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->date_compiled();
  return sptr;
}

void *o2scl_lib_settings_class_o2scl_version(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->o2scl_version();
  return sptr;
}

void o2scl_lib_settings_class_config_h_report(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  ptr->config_h_report();
  return;
}

void *o2scl_lib_settings_class_get_convert_units(void *vptr) {
  lib_settings_class *ptr=(lib_settings_class *)vptr;
  convert_units<> *ret=&ptr->get_convert_units();
  return ret;
}

void *o2scl_create_table__() {
  table<> *ptr=new table<>;
  return ptr;
}

void o2scl_free_table__(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_table__(void *vsrc, void *vdest) {
  table<> *src=(table<> *)vsrc;
  table<> *dest=(table<> *)vdest;
  *dest=*src;
}

void o2scl_table___getitem(void *vptr, char *col, double **dptr, int *n) {
  table<> *ptr=(table<> *)vptr;
  const std::vector<double> &r=ptr->operator[](col);
  *dptr=(double *)(&(r[0]));
  *n=r.size();
  return;
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

size_t o2scl_table___get_maxlines(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->get_maxlines();
  return ret;
}

void o2scl_table___set_maxlines(void *vptr, size_t llines) {
  table<> *ptr=(table<> *)vptr;
  ptr->set_maxlines(llines);
  return;
}

void o2scl_table___set_nlines_auto(void *vptr, size_t il) {
  table<> *ptr=(table<> *)vptr;
  ptr->set_nlines_auto(il);
  return;
}

void o2scl_table___inc_maxlines(void *vptr, size_t llines) {
  table<> *ptr=(table<> *)vptr;
  ptr->inc_maxlines(llines);
  return;
}

void o2scl_table___new_column(void *vptr, char *col) {
  table<> *ptr=(table<> *)vptr;
  ptr->new_column(col);
  return;
}

void *o2scl_table___get_column_name(void *vptr, size_t icol) {
  table<> *ptr=(table<> *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->get_column_name(icol);
  return sptr;
}

void o2scl_table___rename_column(void *vptr, char *src, char *dest) {
  table<> *ptr=(table<> *)vptr;
  ptr->rename_column(src,dest);
  return;
}

void o2scl_table___delete_column(void *vptr, char *col) {
  table<> *ptr=(table<> *)vptr;
  ptr->delete_column(col);
  return;
}

void *o2scl_table___get_sorted_name(void *vptr, size_t icol) {
  table<> *ptr=(table<> *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->get_sorted_name(icol);
  return sptr;
}

void o2scl_table___init_column(void *vptr, char *scol, double val) {
  table<> *ptr=(table<> *)vptr;
  ptr->init_column(scol,val);
  return;
}

bool o2scl_table___is_column(void *vptr, char *scol) {
  table<> *ptr=(table<> *)vptr;
  bool ret=ptr->is_column(scol);
  return ret;
}

size_t o2scl_table___lookup_column(void *vptr, char *scol) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->lookup_column(scol);
  return ret;
}

void o2scl_table___copy_column(void *vptr, char *src, char *dest) {
  table<> *ptr=(table<> *)vptr;
  ptr->copy_column(src,dest);
  return;
}

void o2scl_table___add_col_from_table(void *vptr, void *ptr_source, char *src_index, char *src_col, char *dest_index, char *dest_col) {
  table<> *ptr=(table<> *)vptr;
  table<> *source=(table<> *)ptr_source;
  ptr->add_col_from_table(*source,src_index,src_col,dest_index,dest_col);
  return;
}

void o2scl_table___insert_table(void *vptr, void *ptr_source, char *src_index, bool allow_extrap, char *dest_index) {
  table<> *ptr=(table<> *)vptr;
  table<> *source=(table<> *)ptr_source;
  ptr->insert_table(*source,src_index,allow_extrap,dest_index);
  return;
}

void o2scl_table___add_table(void *vptr, void *ptr_source) {
  table<> *ptr=(table<> *)vptr;
  table<> *source=(table<> *)ptr_source;
  ptr->add_table(*source);
  return;
}

void o2scl_table___new_row(void *vptr, size_t n) {
  table<> *ptr=(table<> *)vptr;
  ptr->new_row(n);
  return;
}

void o2scl_table___copy_row(void *vptr, size_t src, size_t dest) {
  table<> *ptr=(table<> *)vptr;
  ptr->copy_row(src,dest);
  return;
}

void o2scl_table___delete_row(void *vptr, char *scol, double val) {
  table<> *ptr=(table<> *)vptr;
  ptr->delete_row(scol,val);
  return;
}

void o2scl_table___delete_rows_func(void *vptr, char *func) {
  table<> *ptr=(table<> *)vptr;
  ptr->delete_rows_func(func);
  return;
}

void o2scl_table___line_of_names(void *vptr, char *names) {
  table<> *ptr=(table<> *)vptr;
  ptr->line_of_names(names);
  return;
}

void o2scl_table___line_of_data(void *vptr, void *ptr_data) {
  table<> *ptr=(table<> *)vptr;
  std::vector<double> *data=(std::vector<double> *)ptr_data;
  ptr->line_of_data(*data);
  return;
}

size_t o2scl_table___ordered_lookup(void *vptr, char *scol, double val) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->ordered_lookup(scol,val);
  return ret;
}

size_t o2scl_table___lookup(void *vptr, char *scol, double val) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->lookup(scol,val);
  return ret;
}

size_t o2scl_table___lookup_val(void *vptr, char *scol, double val, char *scol2) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->lookup_val(scol,val,scol2);
  return ret;
}

void o2scl_table___set_interp_type(void *vptr, size_t interp_type) {
  table<> *ptr=(table<> *)vptr;
  ptr->set_interp_type(interp_type);
  return;
}

size_t o2scl_table___get_interp_type(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->get_interp_type();
  return ret;
}

double o2scl_table___interp(void *vptr, char *sx, double x0, char *sy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->interp(sx,x0,sy);
  return ret;
}

double o2scl_table___interp_index(void *vptr, size_t ix, double x0, size_t iy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->interp(ix,x0,iy);
  return ret;
}

void o2scl_table___deriv_col(void *vptr, char *x, char *y, char *yp) {
  table<> *ptr=(table<> *)vptr;
  ptr->deriv(x,y,yp);
  return;
}

double o2scl_table___deriv(void *vptr, char *sx, double x0, char *sy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->deriv(sx,x0,sy);
  return ret;
}

double o2scl_table___deriv_index(void *vptr, size_t ix, double x0, size_t iy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->deriv(ix,x0,iy);
  return ret;
}

void o2scl_table___deriv2_col(void *vptr, char *x, char *y, char *yp) {
  table<> *ptr=(table<> *)vptr;
  ptr->deriv2(x,y,yp);
  return;
}

double o2scl_table___deriv2(void *vptr, char *sx, double x0, char *sy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->deriv2(sx,x0,sy);
  return ret;
}

double o2scl_table___deriv2_index(void *vptr, size_t ix, double x0, size_t iy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->deriv2(ix,x0,iy);
  return ret;
}

double o2scl_table___integ(void *vptr, char *sx, double x1, double x2, char *sy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->integ(sx,x1,x2,sy);
  return ret;
}

double o2scl_table___integ_index(void *vptr, size_t ix, double x1, double x2, size_t iy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->integ(ix,x1,x2,iy);
  return ret;
}

void o2scl_table___integ_col(void *vptr, char *x, char *y, char *yi) {
  table<> *ptr=(table<> *)vptr;
  ptr->integ(x,y,yi);
  return;
}

double o2scl_table___max(void *vptr, char *max) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->max(max);
  return ret;
}

double o2scl_table___min(void *vptr, char *min) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->min(min);
  return ret;
}

void o2scl_table___zero_table(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->zero_table();
  return;
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

void o2scl_table___sort_table(void *vptr, char *scol) {
  table<> *ptr=(table<> *)vptr;
  ptr->sort_table(scol);
  return;
}

void o2scl_table___sort_column(void *vptr, char *scol) {
  table<> *ptr=(table<> *)vptr;
  ptr->sort_column(scol);
  return;
}

void o2scl_table___average_col_roll(void *vptr, char *col_name, size_t window) {
  table<> *ptr=(table<> *)vptr;
  ptr->average_col_roll(col_name,window);
  return;
}

void o2scl_table___average_rows(void *vptr, size_t window, bool rolling) {
  table<> *ptr=(table<> *)vptr;
  ptr->average_rows(window,rolling);
  return;
}

void o2scl_table___is_valid(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->is_valid();
  return;
}

void o2scl_table___functions_columns(void *vptr, char *list) {
  table<> *ptr=(table<> *)vptr;
  ptr->functions_columns(list);
  return;
}

void o2scl_table___function_column(void *vptr, char *function, char *scol) {
  table<> *ptr=(table<> *)vptr;
  ptr->function_column(function,scol);
  return;
}

double o2scl_table___row_function(void *vptr, char *scol, size_t row) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->row_function(scol,row);
  return ret;
}

size_t o2scl_table___function_find_row(void *vptr, char *function) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->function_find_row(function);
  return ret;
}

void o2scl_table___summary(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->summary();
  return;
}

void *o2scl_create_table_units__() {
  table_units<> *ptr=new table_units<>;
  return ptr;
}

void o2scl_free_table_units__(void *vptr) {
  table_units<> *ptr=(table_units<> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_table_units__(void *vsrc, void *vdest) {
  table_units<> *src=(table_units<> *)vsrc;
  table_units<> *dest=(table_units<> *)vdest;
  *dest=*src;
}

void o2scl_table_units___set_unit(void *vptr, char *col, char *unit) {
  table_units<> *ptr=(table_units<> *)vptr;
  ptr->set_unit(col,unit);
  return;
}

void *o2scl_table_units___get_unit(void *vptr, char *col) {
  table_units<> *ptr=(table_units<> *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->get_unit(col);
  return sptr;
}

void o2scl_table_units___line_of_units(void *vptr, char *unit_line) {
  table_units<> *ptr=(table_units<> *)vptr;
  ptr->line_of_units(unit_line);
  return;
}

void o2scl_table_units___remove_unit(void *vptr, char *col) {
  table_units<> *ptr=(table_units<> *)vptr;
  ptr->remove_unit(col);
  return;
}

int o2scl_table_units___convert_to_unit(void *vptr, char *col, char *unit, bool err_on_fail) {
  table_units<> *ptr=(table_units<> *)vptr;
  int ret=ptr->convert_to_unit(col,unit,err_on_fail);
  return ret;
}

void *o2scl_create_uniform_grid__() {
  uniform_grid<> *ptr=new uniform_grid<>;
  return ptr;
}

void o2scl_free_uniform_grid__(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  delete ptr;
  return;
}

size_t o2scl_uniform_grid___get_nbins(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  size_t ret=ptr->get_nbins();
  return ret;
}

size_t o2scl_uniform_grid___get_npoints(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  size_t ret=ptr->get_npoints();
  return ret;
}

bool o2scl_uniform_grid___is_log(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  bool ret=ptr->is_log();
  return ret;
}

double o2scl_uniform_grid___get_start(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  double ret=ptr->get_start();
  return ret;
}

double o2scl_uniform_grid___get_end(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  double ret=ptr->get_end();
  return ret;
}

double o2scl_uniform_grid___get_width(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  double ret=ptr->get_width();
  return ret;
}

double o2scl_uniform_grid___getitem(void *vptr, size_t n) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  double ret=ptr->operator[](n);
  return ret;
}

void o2scl_uniform_grid___vector(void *vptr, void *ptr_v) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  std::vector<double> *v=(std::vector<double> *)ptr_v;
  ptr->vector(*v);
  return;
}

void *o2scl_create_uniform_grid_end__() {
  uniform_grid_end<> *ptr=new uniform_grid_end<>;
  return ptr;
}

void o2scl_free_uniform_grid_end__(void *vptr) {
  uniform_grid_end<> *ptr=(uniform_grid_end<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_end___init(double start, double end, size_t n_bins) {
  uniform_grid_end<> *ptr=new uniform_grid_end<>(start,end,n_bins);
  return ptr;
}

void *o2scl_create_uniform_grid_width__() {
  uniform_grid_width<> *ptr=new uniform_grid_width<>;
  return ptr;
}

void o2scl_free_uniform_grid_width__(void *vptr) {
  uniform_grid_width<> *ptr=(uniform_grid_width<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_width___init(double start, double width, size_t n_bins) {
  uniform_grid_width<> *ptr=new uniform_grid_width<>(start,width,n_bins);
  return ptr;
}

void *o2scl_create_uniform_grid_end_width__() {
  uniform_grid_end_width<> *ptr=new uniform_grid_end_width<>;
  return ptr;
}

void o2scl_free_uniform_grid_end_width__(void *vptr) {
  uniform_grid_end_width<> *ptr=(uniform_grid_end_width<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_end_width___init(double start, double end, double width) {
  uniform_grid_end_width<> *ptr=new uniform_grid_end_width<>(start,end,width);
  return ptr;
}

void *o2scl_create_uniform_grid_log_end__() {
  uniform_grid_log_end<> *ptr=new uniform_grid_log_end<>;
  return ptr;
}

void o2scl_free_uniform_grid_log_end__(void *vptr) {
  uniform_grid_log_end<> *ptr=(uniform_grid_log_end<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_log_end___init(double start, double end, size_t n_bins) {
  uniform_grid_log_end<> *ptr=new uniform_grid_log_end<>(start,end,n_bins);
  return ptr;
}

void *o2scl_create_uniform_grid_log_width__() {
  uniform_grid_log_width<> *ptr=new uniform_grid_log_width<>;
  return ptr;
}

void o2scl_free_uniform_grid_log_width__(void *vptr) {
  uniform_grid_log_width<> *ptr=(uniform_grid_log_width<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_log_width___init(double start, double width, size_t n_bins) {
  uniform_grid_log_width<> *ptr=new uniform_grid_log_width<>(start,width,n_bins);
  return ptr;
}

void *o2scl_create_uniform_grid_log_end_width__() {
  uniform_grid_log_end_width<> *ptr=new uniform_grid_log_end_width<>;
  return ptr;
}

void o2scl_free_uniform_grid_log_end_width__(void *vptr) {
  uniform_grid_log_end_width<> *ptr=(uniform_grid_log_end_width<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_log_end_width___init(double start, double end, double width) {
  uniform_grid_log_end_width<> *ptr=new uniform_grid_log_end_width<>(start,end,width);
  return ptr;
}

void *o2scl_create_table3d() {
  table3d *ptr=new table3d;
  return ptr;
}

void o2scl_free_table3d(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_table3d(void *vsrc, void *vdest) {
  table3d *src=(table3d *)vsrc;
  table3d *dest=(table3d *)vdest;
  *dest=*src;
}

void o2scl_table3d_set_size(void *vptr, size_t nx, size_t ny) {
  table3d *ptr=(table3d *)vptr;
  ptr->set_size(nx,ny);
  return;
}

void o2scl_table3d_set_xy(void *vptr, char *x_name, size_t nx, void *ptr_x, char *y_name, size_t ny, void *ptr_y) {
  table3d *ptr=(table3d *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  ptr->set_xy(x_name,nx,*x,y_name,ny,*y);
  return;
}

void o2scl_table3d_set_xy_grid(void *vptr, char *x_name, void *ptr_x_grid, char *y_name, void *ptr_y_grid) {
  table3d *ptr=(table3d *)vptr;
  uniform_grid<double> *x_grid=(uniform_grid<double> *)ptr_x_grid;
  uniform_grid<double> *y_grid=(uniform_grid<double> *)ptr_y_grid;
  ptr->set_xy(x_name,*x_grid,y_name,*y_grid);
  return;
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

double o2scl_table3d_get_i(void *vptr, size_t ix, size_t iy, size_t iz) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->get(ix,iy,iz);
  return ret;
}

void o2scl_table3d_set_i(void *vptr, size_t ix, size_t iy, size_t iz, double val) {
  table3d *ptr=(table3d *)vptr;
  ptr->set(ix,iy,iz,val);
  return;
}

void o2scl_table3d_set_val(void *vptr, double x, double y, char *name, double val) {
  table3d *ptr=(table3d *)vptr;
  ptr->set_val(x,y,name,val);
  return;
}

double o2scl_table3d_get_val(void *vptr, double x, double y, char *name) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->get_val(x,y,name);
  return ret;
}

void o2scl_table3d_set_grid_x(void *vptr, size_t ix, double val) {
  table3d *ptr=(table3d *)vptr;
  ptr->set_grid_x(ix,val);
  return;
}

void o2scl_table3d_set_grid_y(void *vptr, size_t iy, double val) {
  table3d *ptr=(table3d *)vptr;
  ptr->set_grid_y(iy,val);
  return;
}

double o2scl_table3d_get_grid_x(void *vptr, size_t ix) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->get_grid_x(ix);
  return ret;
}

double o2scl_table3d_get_grid_y(void *vptr, size_t iy) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->get_grid_y(iy);
  return ret;
}

void o2scl_table3d_get_size(void *vptr, size_t *nx, size_t *ny) {
  table3d *ptr=(table3d *)vptr;
  ptr->get_size(*nx,*ny);
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

bool o2scl_table3d_is_size_set(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  bool ret=ptr->is_size_set();
  return ret;
}

bool o2scl_table3d_is_xy_set(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  bool ret=ptr->is_xy_set();
  return ret;
}

void *o2scl_table3d_get_slice_name(void *vptr, size_t i) {
  table3d *ptr=(table3d *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->get_slice_name(i);
  return sptr;
}

void o2scl_table3d_new_slice(void *vptr, char *slice) {
  table3d *ptr=(table3d *)vptr;
  ptr->new_slice(slice);
  return;
}

void o2scl_table3d_set_slice_all(void *vptr, char *name, double val) {
  table3d *ptr=(table3d *)vptr;
  ptr->set_slice_all(name,val);
  return;
}

size_t o2scl_table3d_lookup_slice(void *vptr, char *name) {
  table3d *ptr=(table3d *)vptr;
  size_t ret=ptr->lookup_slice(name);
  return ret;
}

bool o2scl_table3d_is_slice(void *vptr, char *name, size_t *ix) {
  table3d *ptr=(table3d *)vptr;
  bool ret=ptr->is_slice(name,*ix);
  return ret;
}

void o2scl_table3d_rename_slice(void *vptr, char *name1, char *name2) {
  table3d *ptr=(table3d *)vptr;
  ptr->rename_slice(name1,name2);
  return;
}

void o2scl_table3d_copy_slice(void *vptr, char *name1, char *name2) {
  table3d *ptr=(table3d *)vptr;
  ptr->copy_slice(name1,name2);
  return;
}

void *o2scl_table3d_get_slice(void *vptr, char *slice) {
  table3d *ptr=(table3d *)vptr;
  boost::numeric::ublas::matrix<double> *ret=&ptr->get_slice(slice);
  return ret;
}

void *o2scl_table3d_get_slice_i(void *vptr, char *slice) {
  table3d *ptr=(table3d *)vptr;
  boost::numeric::ublas::matrix<double> *ret=&ptr->get_slice(slice);
  return ret;
}

void o2scl_table3d_lookup_x(void *vptr, double val, size_t ix) {
  table3d *ptr=(table3d *)vptr;
  ptr->lookup_x(val,ix);
  return;
}

void o2scl_table3d_lookup_y(void *vptr, double val, size_t iy) {
  table3d *ptr=(table3d *)vptr;
  ptr->lookup_y(val,iy);
  return;
}

double o2scl_table3d_interp(void *vptr, double x, double y, char *name) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->interp(x,y,name);
  return ret;
}

double o2scl_table3d_deriv_x(void *vptr, double x, double y, char *name) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->deriv_x(x,y,name);
  return ret;
}

double o2scl_table3d_deriv_y(void *vptr, double x, double y, char *name) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->deriv_y(x,y,name);
  return ret;
}

double o2scl_table3d_deriv_xy(void *vptr, double x, double y, char *name) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->deriv_xy(x,y,name);
  return ret;
}

double o2scl_table3d_integ_x(void *vptr, double x1, double x2, double y, char *name) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->integ_x(x1,x2,y,name);
  return ret;
}

double o2scl_table3d_integ_y(void *vptr, double x, double y1, double y2, char *name) {
  table3d *ptr=(table3d *)vptr;
  double ret=ptr->integ_y(x,y1,y2,name);
  return ret;
}

void o2scl_table3d_zero_table(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  ptr->zero_table();
  return;
}

void o2scl_table3d_clear(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  ptr->clear();
  return;
}

int o2scl_table3d_function_matrix(void *vptr, char *function, void *ptr_mat, bool throw_on_err) {
  table3d *ptr=(table3d *)vptr;
  boost::numeric::ublas::matrix<double> *mat=(boost::numeric::ublas::matrix<double> *)ptr_mat;
  int ret=ptr->function_matrix(function,*mat,throw_on_err);
  return ret;
}

void o2scl_table3d_function_slice(void *vptr, char *function, char *slice) {
  table3d *ptr=(table3d *)vptr;
  ptr->function_slice(function,slice);
  return;
}

void o2scl_table3d_summary(void *vptr) {
  table3d *ptr=(table3d *)vptr;
  ptr->summary();
  return;
}

void *o2scl_create_index_spec() {
  index_spec *ptr=new index_spec;
  return ptr;
}

void o2scl_free_index_spec(void *vptr) {
  index_spec *ptr=(index_spec *)vptr;
  delete ptr;
  return;
}

size_t o2scl_index_spec_get_type(void *vptr) {
  index_spec *ptr=(index_spec *)vptr;
  return ptr->type;
}

void o2scl_index_spec_set_type(void *vptr, size_t v) {
  index_spec *ptr=(index_spec *)vptr;
  ptr->type=v;
  return;
}

size_t o2scl_index_spec_get_ix1(void *vptr) {
  index_spec *ptr=(index_spec *)vptr;
  return ptr->ix1;
}

void o2scl_index_spec_set_ix1(void *vptr, size_t v) {
  index_spec *ptr=(index_spec *)vptr;
  ptr->ix1=v;
  return;
}

size_t o2scl_index_spec_get_ix2(void *vptr) {
  index_spec *ptr=(index_spec *)vptr;
  return ptr->ix2;
}

void o2scl_index_spec_set_ix2(void *vptr, size_t v) {
  index_spec *ptr=(index_spec *)vptr;
  ptr->ix2=v;
  return;
}

size_t o2scl_index_spec_get_ix3(void *vptr) {
  index_spec *ptr=(index_spec *)vptr;
  return ptr->ix3;
}

void o2scl_index_spec_set_ix3(void *vptr, size_t v) {
  index_spec *ptr=(index_spec *)vptr;
  ptr->ix3=v;
  return;
}

double o2scl_index_spec_get_val1(void *vptr) {
  index_spec *ptr=(index_spec *)vptr;
  return ptr->val1;
}

void o2scl_index_spec_set_val1(void *vptr, double v) {
  index_spec *ptr=(index_spec *)vptr;
  ptr->val1=v;
  return;
}

double o2scl_index_spec_get_val2(void *vptr) {
  index_spec *ptr=(index_spec *)vptr;
  return ptr->val2;
}

void o2scl_index_spec_set_val2(void *vptr, double v) {
  index_spec *ptr=(index_spec *)vptr;
  ptr->val2=v;
  return;
}

double o2scl_index_spec_get_val3(void *vptr) {
  index_spec *ptr=(index_spec *)vptr;
  return ptr->val3;
}

void o2scl_index_spec_set_val3(void *vptr, double v) {
  index_spec *ptr=(index_spec *)vptr;
  ptr->val3=v;
  return;
}

void *o2scl_create_tensor__() {
  tensor<> *ptr=new tensor<>;
  return ptr;
}

void o2scl_free_tensor__(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_tensor__(void *vsrc, void *vdest) {
  tensor<> *src=(tensor<> *)vsrc;
  tensor<> *dest=(tensor<> *)vdest;
  *dest=*src;
}

void o2scl_tensor___is_valid(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  ptr->is_valid();
  return;
}

void o2scl_tensor___clear(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  ptr->clear();
  return;
}

void o2scl_tensor___set(void *vptr, void *ptr_index, double val) {
  tensor<> *ptr=(tensor<> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  ptr->set(*index,val);
  return;
}

void o2scl_tensor___set_all(void *vptr, double x) {
  tensor<> *ptr=(tensor<> *)vptr;
  ptr->set_all(x);
  return;
}

double o2scl_tensor___get(void *vptr, void *ptr_index) {
  tensor<> *ptr=(tensor<> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  double ret=ptr->get(*index);
  return ret;
}

void o2scl_tensor___resize(void *vptr, size_t n, void *ptr_index) {
  tensor<> *ptr=(tensor<> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  ptr->resize(n,*index);
  return;
}

size_t o2scl_tensor___get_rank(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  size_t ret=ptr->get_rank();
  return ret;
}

size_t o2scl_tensor___get_size(void *vptr, size_t i) {
  tensor<> *ptr=(tensor<> *)vptr;
  size_t ret=ptr->get_size(i);
  return ret;
}

void o2scl_tensor___get_data(void *vptr, double **dptr, int *n) {
  tensor<> *ptr=(tensor<> *)vptr;
  const std::vector<double> &r=ptr->get_data();
  *dptr=(double *)(&(r[0]));
  *n=r.size();
  return;
}

size_t o2scl_tensor___total_size(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  size_t ret=ptr->total_size();
  return ret;
}

double o2scl_tensor___min_value(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  double ret=ptr->min_value();
  return ret;
}

double o2scl_tensor___max_value(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  double ret=ptr->max_value();
  return ret;
}

double o2scl_tensor___total_sum(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  double ret=ptr->total_sum();
  return ret;
}

void *o2scl_create_tensor_grid__() {
  tensor_grid<> *ptr=new tensor_grid<>;
  return ptr;
}

void o2scl_free_tensor_grid__(void *vptr) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_tensor_grid__(void *vsrc, void *vdest) {
  tensor_grid<> *src=(tensor_grid<> *)vsrc;
  tensor_grid<> *dest=(tensor_grid<> *)vdest;
  *dest=*src;
}

void o2scl_tensor_grid___is_valid(void *vptr) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  ptr->is_valid();
  return;
}

void o2scl_tensor_grid___set_val(void *vptr, void *ptr_grid_point, double val) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  vector<double> *grid_point=(vector<double> *)ptr_grid_point;
  ptr->set_val(*grid_point,val);
  return;
}

double o2scl_tensor_grid___get_val(void *vptr, void *ptr_grid_point) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  vector<double> *grid_point=(vector<double> *)ptr_grid_point;
  double ret=ptr->get_val(*grid_point);
  return ret;
}

bool o2scl_tensor_grid___is_grid_set(void *vptr) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  bool ret=ptr->is_grid_set();
  return ret;
}

void o2scl_tensor_grid___set_grid_packed(void *vptr, void *ptr_grid) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  vector<double> *grid=(vector<double> *)ptr_grid;
  ptr->set_grid_packed(*grid);
  return;
}

void o2scl_tensor_grid___default_grid(void *vptr) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  ptr->default_grid();
  return;
}

void o2scl_tensor_grid___set_grid_i_vec(void *vptr, size_t i, void *ptr_grid) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  vector<double> *grid=(vector<double> *)ptr_grid;
  ptr->set_grid_i_vec(i,*grid);
  return;
}

double o2scl_tensor_grid___get_grid(void *vptr, size_t i, size_t j) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  double ret=ptr->get_grid(i,j);
  return ret;
}

void o2scl_tensor_grid___set_grid(void *vptr, size_t i, size_t j, double val) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  ptr->set_grid(i,j,val);
  return;
}

void *o2scl_create_find_constants() {
  find_constants *ptr=new find_constants;
  return ptr;
}

void o2scl_free_find_constants(void *vptr) {
  find_constants *ptr=(find_constants *)vptr;
  delete ptr;
  return;
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
  return;
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

void *o2scl_convert_units___get_units_cmd_string(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->units_cmd_string;
  return sptr;
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
  std::shared_ptr<table_units<> > *ptr=new std::shared_ptr<table_units<> >(new table_units<>);
  return ptr;
}

void o2scl_free_shared_ptr_table_units__(void *vptr) {
  std::shared_ptr<table_units<> > *ptr=(std::shared_ptr<table_units<> > *)vptr;
  delete ptr;
}

void *o2scl_shared_ptr_table_units___ptr(void *vp) {
  std::shared_ptr<table_units<> > *p=(std::shared_ptr<table_units<> > *)vp;
  table_units<> *ref=p->get();
  return ref;
}

void *o2scl_ix_index_wrapper(size_t ix) {
  index_spec *ret=new index_spec;
  *ret=ix_index(ix);
  return ret;
}

void *o2scl_ix_fixed_wrapper(size_t ix, size_t ix2) {
  index_spec *ret=new index_spec;
  *ret=ix_fixed(ix,ix2);
  return ret;
}

void *o2scl_ix_sum_wrapper(size_t ix) {
  index_spec *ret=new index_spec;
  *ret=ix_sum(ix);
  return ret;
}

void *o2scl_ix_trace_wrapper(size_t ix, size_t ix2) {
  index_spec *ret=new index_spec;
  *ret=ix_trace(ix,ix2);
  return ret;
}

void *o2scl_ix_reverse_wrapper(size_t ix) {
  index_spec *ret=new index_spec;
  *ret=ix_reverse(ix);
  return ret;
}

void *o2scl_ix_range_wrapper(size_t ix, size_t start, size_t end) {
  index_spec *ret=new index_spec;
  *ret=ix_range(ix,start,end);
  return ret;
}

void *o2scl_ix_interp_wrapper(size_t ix, double v) {
  index_spec *ret=new index_spec;
  *ret=ix_interp(ix,v);
  return ret;
}

void *o2scl_ix_grid_wrapper(size_t ix, double start, double end, size_t n_bins, bool log) {
  index_spec *ret=new index_spec;
  *ret=ix_grid(ix,start,end,n_bins,log);
  return ret;
}

void *o2scl_ix_gridw_wrapper(size_t ix, double start, double end, double width, bool log) {
  index_spec *ret=new index_spec;
  *ret=ix_gridw(ix,start,end,width,log);
  return ret;
}

