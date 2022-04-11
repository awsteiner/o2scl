/*
  -------------------------------------------------------------------

  Copyright (C) 2020-2022, Andrew W. Steiner

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

void *o2scl_create_std_string() {
  std::string *ptr=new std::string;
  return ptr;
}

void o2scl_free_std_string(void *vptr) {
  std::string *ptr=(std::string *)vptr;
  delete ptr;
  return;
}

size_t o2scl_std_string_length(void *vptr) {
  std::string *ptr=(std::string *)vptr;
  size_t ret=ptr->length();
  return ret;
}

char o2scl_std_string_getitem(void *vptr, size_t n) {
  std::string *ptr=(std::string *)vptr;
  char ret=ptr->operator[](n);
  return ret;
}

void o2scl_std_string_setitem(void *vptr, size_t i, char val) {
  std::string *ptr=(std::string *)vptr;
  (*ptr)[i]=val;
  return;
}

void o2scl_std_string_resize(void *vptr, size_t n) {
  std::string *ptr=(std::string *)vptr;
  ptr->resize(n);
  return;
}

void *o2scl_create_std_vector_double_() {
  std::vector<double> *ptr=new std::vector<double>;
  return ptr;
}

void o2scl_free_std_vector_double_(void *vptr) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  delete ptr;
  return;
}

void o2scl_std_vector_double__resize(void *vptr, size_t n) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std_vector_double__size(void *vptr) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

double o2scl_std_vector_double__getitem(void *vptr, size_t n) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  double ret=ptr->operator[](n);
  return ret;
}

void o2scl_std_vector_double__setitem(void *vptr, size_t i, double val) {
  std::vector<double> *ptr=(std::vector<double> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_std_vector_int_() {
  std::vector<int> *ptr=new std::vector<int>;
  return ptr;
}

void o2scl_free_std_vector_int_(void *vptr) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  delete ptr;
  return;
}

void o2scl_std_vector_int__resize(void *vptr, size_t n) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std_vector_int__size(void *vptr) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

int o2scl_std_vector_int__getitem(void *vptr, size_t n) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  int ret=ptr->operator[](n);
  return ret;
}

void o2scl_std_vector_int__setitem(void *vptr, size_t i, int val) {
  std::vector<int> *ptr=(std::vector<int> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_std_vector_size_t_() {
  std::vector<size_t> *ptr=new std::vector<size_t>;
  return ptr;
}

void o2scl_free_std_vector_size_t_(void *vptr) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  delete ptr;
  return;
}

void o2scl_std_vector_size_t__resize(void *vptr, size_t n) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std_vector_size_t__size(void *vptr) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

size_t o2scl_std_vector_size_t__getitem(void *vptr, size_t n) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  size_t ret=ptr->operator[](n);
  return ret;
}

void o2scl_std_vector_size_t__setitem(void *vptr, size_t i, size_t val) {
  std::vector<size_t> *ptr=(std::vector<size_t> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_std_vector_std_string_() {
  std::vector<std::string> *ptr=new std::vector<std::string>;
  return ptr;
}

void o2scl_free_std_vector_std_string_(void *vptr) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  delete ptr;
  return;
}

void o2scl_std_vector_std_string__resize(void *vptr, size_t n) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std_vector_std_string__size(void *vptr) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

void *o2scl_std_vector_std_string__getitem(void *vptr, size_t n) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->operator[](n);
  return sptr;
}

void o2scl_std_vector_std_string__setitem(void *vptr, size_t i, char *val) {
  std::vector<std::string> *ptr=(std::vector<std::string> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_boost_numeric_ublas_vector_double_() {
  boost::numeric::ublas::vector<double> *ptr=new boost::numeric::ublas::vector<double>;
  return ptr;
}

void o2scl_free_boost_numeric_ublas_vector_double_(void *vptr) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  delete ptr;
  return;
}

size_t o2scl_boost_numeric_ublas_vector_double__size(void *vptr) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

void o2scl_boost_numeric_ublas_vector_double__resize(void *vptr, size_t n) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  ptr->resize(n);
  return;
}

double o2scl_boost_numeric_ublas_vector_double__getitem(void *vptr, size_t i) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  double ret=ptr->operator[](i);
  return ret;
}

void o2scl_boost_numeric_ublas_vector_double__setitem(void *vptr, size_t i, double val) {
  boost::numeric::ublas::vector<double> *ptr=(boost::numeric::ublas::vector<double> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_boost_numeric_ublas_matrix_double_() {
  boost::numeric::ublas::matrix<double> *ptr=new boost::numeric::ublas::matrix<double>;
  return ptr;
}

void o2scl_free_boost_numeric_ublas_matrix_double_(void *vptr) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  delete ptr;
  return;
}

size_t o2scl_boost_numeric_ublas_matrix_double__size1(void *vptr) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  size_t ret=ptr->size1();
  return ret;
}

size_t o2scl_boost_numeric_ublas_matrix_double__size2(void *vptr) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  size_t ret=ptr->size2();
  return ret;
}

void o2scl_boost_numeric_ublas_matrix_double__resize(void *vptr, size_t m, size_t n) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  ptr->resize(m,n);
  return;
}

double o2scl_boost_numeric_ublas_matrix_double__getitem(void *vptr, size_t m, size_t n) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  double ret=ptr->operator()(m,n);
  return ret;
}

void o2scl_boost_numeric_ublas_matrix_double__setitem(void *vptr, size_t i, size_t j, double val) {
  boost::numeric::ublas::matrix<double> *ptr=(boost::numeric::ublas::matrix<double> *)vptr;
  (*ptr)(i,j)=val;
  return;
}

void *o2scl_create_boost_numeric_ublas_matrix_int_() {
  boost::numeric::ublas::matrix<int> *ptr=new boost::numeric::ublas::matrix<int>;
  return ptr;
}

void o2scl_free_boost_numeric_ublas_matrix_int_(void *vptr) {
  boost::numeric::ublas::matrix<int> *ptr=(boost::numeric::ublas::matrix<int> *)vptr;
  delete ptr;
  return;
}

size_t o2scl_boost_numeric_ublas_matrix_int__size1(void *vptr) {
  boost::numeric::ublas::matrix<int> *ptr=(boost::numeric::ublas::matrix<int> *)vptr;
  size_t ret=ptr->size1();
  return ret;
}

size_t o2scl_boost_numeric_ublas_matrix_int__size2(void *vptr) {
  boost::numeric::ublas::matrix<int> *ptr=(boost::numeric::ublas::matrix<int> *)vptr;
  size_t ret=ptr->size2();
  return ret;
}

void o2scl_boost_numeric_ublas_matrix_int__resize(void *vptr, size_t m, size_t n) {
  boost::numeric::ublas::matrix<int> *ptr=(boost::numeric::ublas::matrix<int> *)vptr;
  ptr->resize(m,n);
  return;
}

int o2scl_boost_numeric_ublas_matrix_int__getitem(void *vptr, size_t m, size_t n) {
  boost::numeric::ublas::matrix<int> *ptr=(boost::numeric::ublas::matrix<int> *)vptr;
  int ret=ptr->operator()(m,n);
  return ret;
}

void o2scl_boost_numeric_ublas_matrix_int__setitem(void *vptr, size_t i, size_t j, int val) {
  boost::numeric::ublas::matrix<int> *ptr=(boost::numeric::ublas::matrix<int> *)vptr;
  (*ptr)(i,j)=val;
  return;
}

void *o2scl_create_std_vector_std_vector_double_() {
  std::vector<std::vector<double>> *ptr=new std::vector<std::vector<double>>;
  return ptr;
}

void o2scl_free_std_vector_std_vector_double_(void *vptr) {
  std::vector<std::vector<double>> *ptr=(std::vector<std::vector<double>> *)vptr;
  delete ptr;
  return;
}

void o2scl_std_vector_std_vector_double__resize(void *vptr, size_t n) {
  std::vector<std::vector<double>> *ptr=(std::vector<std::vector<double>> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std_vector_std_vector_double__size(void *vptr) {
  std::vector<std::vector<double>> *ptr=(std::vector<std::vector<double>> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

void o2scl_std_vector_std_vector_double__getitem(void *vptr, size_t n, double **dptr, int *n_) {
  std::vector<std::vector<double>> *ptr=(std::vector<std::vector<double>> *)vptr;
  const std::vector<double> &r=ptr->operator[](n);
  *dptr=(double *)(&(r[0]));
  *n_=r.size();
  return;
}

void o2scl_std_vector_std_vector_double__setitem(void *vptr, size_t i, void *valptr) {
  std::vector<std::vector<double>> *ptr=(std::vector<std::vector<double>> *)vptr;
  std::vector<double> *valptr2=(std::vector<double> *)valptr;
  (*ptr)[i]=*valptr2;
  return;
}

void *o2scl_create_std_vector_std_vector_std_string_() {
  std::vector<std::vector<std::string>> *ptr=new std::vector<std::vector<std::string>>;
  return ptr;
}

void o2scl_free_std_vector_std_vector_std_string_(void *vptr) {
  std::vector<std::vector<std::string>> *ptr=(std::vector<std::vector<std::string>> *)vptr;
  delete ptr;
  return;
}

void o2scl_std_vector_std_vector_std_string__resize(void *vptr, size_t n) {
  std::vector<std::vector<std::string>> *ptr=(std::vector<std::vector<std::string>> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std_vector_std_vector_std_string__size(void *vptr) {
  std::vector<std::vector<std::string>> *ptr=(std::vector<std::vector<std::string>> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

void *o2scl_std_vector_std_vector_std_string__getitem(void *vptr, size_t n) {
  std::vector<std::vector<std::string>> *ptr=(std::vector<std::vector<std::string>> *)vptr;
  std::vector<std::string> *vsptr=new std::vector<std::string>;
  *vsptr=ptr->operator[](n);
  return vsptr;
}

void o2scl_std_vector_std_vector_std_string__setitem(void *vptr, size_t i, std::vector<std::string> val) {
  std::vector<std::vector<std::string>> *ptr=(std::vector<std::vector<std::string>> *)vptr;
  (*ptr)[i]=val;
  return;
}

void *o2scl_create_std_complex_double_() {
  std::complex<double> *ptr=new std::complex<double>;
  return ptr;
}

void o2scl_free_std_complex_double_(void *vptr) {
  std::complex<double> *ptr=(std::complex<double> *)vptr;
  delete ptr;
  return;
}

double o2scl_std_complex_double__real(void *vptr) {
  std::complex<double> *ptr=(std::complex<double> *)vptr;
  double ret=ptr->real();
  return ret;
}

void o2scl_std_complex_double__real_set(void *vptr, double value) {
  std::complex<double> *ptr=(std::complex<double> *)vptr;
  ptr->real(value);
  return;
}

double o2scl_std_complex_double__imag(void *vptr) {
  std::complex<double> *ptr=(std::complex<double> *)vptr;
  double ret=ptr->imag();
  return ret;
}

void o2scl_std_complex_double__imag_set(void *vptr, double value) {
  std::complex<double> *ptr=(std::complex<double> *)vptr;
  ptr->imag(value);
  return;
}

void *o2scl_std_complex_double__init(double re, double im) {
  std::complex<double> *ptr=new std::complex<double>(re,im);
  return ptr;
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

void *o2scl_create_table_() {
  table<> *ptr=new table<>;
  return ptr;
}

void o2scl_free_table_(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_table_(void *vsrc, void *vdest) {
  table<> *src=(table<> *)vsrc;
  table<> *dest=(table<> *)vdest;
  *dest=*src;
}

void o2scl_table__getitem(void *vptr, char *col, double **dptr, int *n_) {
  table<> *ptr=(table<> *)vptr;
  const std::vector<double> &r=ptr->operator[](col);
  *dptr=(double *)(&(r[0]));
  *n_=r.size();
  return;
}

void o2scl_table__set(void *vptr, char *col, size_t row, double val) {
  table<> *ptr=(table<> *)vptr;
  ptr->set(col,row,val);
  return;
}

double o2scl_table__get(void *vptr, char *col, size_t row) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->get(col,row);
  return ret;
}

size_t o2scl_table__get_ncolumns(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->get_ncolumns();
  return ret;
}

size_t o2scl_table__get_nlines(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->get_nlines();
  return ret;
}

void o2scl_table__set_nlines(void *vptr, size_t lines) {
  table<> *ptr=(table<> *)vptr;
  ptr->set_nlines(lines);
  return;
}

size_t o2scl_table__get_maxlines(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->get_maxlines();
  return ret;
}

void o2scl_table__set_maxlines(void *vptr, size_t llines) {
  table<> *ptr=(table<> *)vptr;
  ptr->set_maxlines(llines);
  return;
}

void o2scl_table__set_nlines_auto(void *vptr, size_t il) {
  table<> *ptr=(table<> *)vptr;
  ptr->set_nlines_auto(il);
  return;
}

void o2scl_table__inc_maxlines(void *vptr, size_t llines) {
  table<> *ptr=(table<> *)vptr;
  ptr->inc_maxlines(llines);
  return;
}

void o2scl_table__new_column(void *vptr, char *col) {
  table<> *ptr=(table<> *)vptr;
  ptr->new_column(col);
  return;
}

void *o2scl_table__get_column_name(void *vptr, size_t icol) {
  table<> *ptr=(table<> *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->get_column_name(icol);
  return sptr;
}

void o2scl_table__rename_column(void *vptr, char *src, char *dest) {
  table<> *ptr=(table<> *)vptr;
  ptr->rename_column(src,dest);
  return;
}

void o2scl_table__delete_column(void *vptr, char *col) {
  table<> *ptr=(table<> *)vptr;
  ptr->delete_column(col);
  return;
}

void *o2scl_table__get_sorted_name(void *vptr, size_t icol) {
  table<> *ptr=(table<> *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->get_sorted_name(icol);
  return sptr;
}

void o2scl_table__init_column(void *vptr, char *scol, double val) {
  table<> *ptr=(table<> *)vptr;
  ptr->init_column(scol,val);
  return;
}

bool o2scl_table__is_column(void *vptr, char *scol) {
  table<> *ptr=(table<> *)vptr;
  bool ret=ptr->is_column(scol);
  return ret;
}

size_t o2scl_table__lookup_column(void *vptr, char *scol) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->lookup_column(scol);
  return ret;
}

void o2scl_table__copy_column(void *vptr, char *src, char *dest) {
  table<> *ptr=(table<> *)vptr;
  ptr->copy_column(src,dest);
  return;
}

void o2scl_table__add_col_from_table(void *vptr, void *ptr_source, char *src_index, char *src_col, char *dest_index, char *dest_col) {
  table<> *ptr=(table<> *)vptr;
  table<> *source=(table<> *)ptr_source;
  ptr->add_col_from_table(*source,src_index,src_col,dest_index,dest_col);
  return;
}

void o2scl_table__insert_table(void *vptr, void *ptr_source, char *src_index, bool allow_extrap, char *dest_index) {
  table<> *ptr=(table<> *)vptr;
  table<> *source=(table<> *)ptr_source;
  ptr->insert_table(*source,src_index,allow_extrap,dest_index);
  return;
}

void o2scl_table__add_table(void *vptr, void *ptr_source) {
  table<> *ptr=(table<> *)vptr;
  table<> *source=(table<> *)ptr_source;
  ptr->add_table(*source);
  return;
}

void o2scl_table__new_row(void *vptr, size_t n) {
  table<> *ptr=(table<> *)vptr;
  ptr->new_row(n);
  return;
}

void o2scl_table__copy_row(void *vptr, size_t src, size_t dest) {
  table<> *ptr=(table<> *)vptr;
  ptr->copy_row(src,dest);
  return;
}

void o2scl_table__delete_row(void *vptr, char *scol, double val) {
  table<> *ptr=(table<> *)vptr;
  ptr->delete_row(scol,val);
  return;
}

void o2scl_table__delete_rows_func(void *vptr, char *func) {
  table<> *ptr=(table<> *)vptr;
  ptr->delete_rows_func(func);
  return;
}

void o2scl_table__line_of_names(void *vptr, char *names) {
  table<> *ptr=(table<> *)vptr;
  ptr->line_of_names(names);
  return;
}

void o2scl_table__line_of_data(void *vptr, void *ptr_data) {
  table<> *ptr=(table<> *)vptr;
  std::vector<double> *data=(std::vector<double> *)ptr_data;
  ptr->line_of_data(*data);
  return;
}

size_t o2scl_table__ordered_lookup(void *vptr, char *scol, double val) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->ordered_lookup(scol,val);
  return ret;
}

size_t o2scl_table__lookup(void *vptr, char *scol, double val) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->lookup(scol,val);
  return ret;
}

size_t o2scl_table__lookup_val(void *vptr, char *scol, double val, char *scol2) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->lookup_val(scol,val,scol2);
  return ret;
}

void o2scl_table__set_interp_type(void *vptr, size_t interp_type) {
  table<> *ptr=(table<> *)vptr;
  ptr->set_interp_type(interp_type);
  return;
}

size_t o2scl_table__get_interp_type(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->get_interp_type();
  return ret;
}

double o2scl_table__interp(void *vptr, char *sx, double x0, char *sy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->interp(sx,x0,sy);
  return ret;
}

double o2scl_table__interp_index(void *vptr, size_t ix, double x0, size_t iy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->interp(ix,x0,iy);
  return ret;
}

void o2scl_table__deriv_col(void *vptr, char *x, char *y, char *yp) {
  table<> *ptr=(table<> *)vptr;
  ptr->deriv(x,y,yp);
  return;
}

double o2scl_table__deriv(void *vptr, char *sx, double x0, char *sy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->deriv(sx,x0,sy);
  return ret;
}

double o2scl_table__deriv_index(void *vptr, size_t ix, double x0, size_t iy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->deriv(ix,x0,iy);
  return ret;
}

void o2scl_table__deriv2_col(void *vptr, char *x, char *y, char *yp) {
  table<> *ptr=(table<> *)vptr;
  ptr->deriv2(x,y,yp);
  return;
}

double o2scl_table__deriv2(void *vptr, char *sx, double x0, char *sy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->deriv2(sx,x0,sy);
  return ret;
}

double o2scl_table__deriv2_index(void *vptr, size_t ix, double x0, size_t iy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->deriv2(ix,x0,iy);
  return ret;
}

double o2scl_table__integ(void *vptr, char *sx, double x1, double x2, char *sy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->integ(sx,x1,x2,sy);
  return ret;
}

double o2scl_table__integ_index(void *vptr, size_t ix, double x1, double x2, size_t iy) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->integ(ix,x1,x2,iy);
  return ret;
}

void o2scl_table__integ_col(void *vptr, char *x, char *y, char *yi) {
  table<> *ptr=(table<> *)vptr;
  ptr->integ(x,y,yi);
  return;
}

double o2scl_table__max(void *vptr, char *max) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->max(max);
  return ret;
}

double o2scl_table__min(void *vptr, char *min) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->min(min);
  return ret;
}

void o2scl_table__zero_table(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->zero_table();
  return;
}

void o2scl_table__clear(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->clear();
  return;
}

void o2scl_table__clear_data(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->clear_data();
  return;
}

void o2scl_table__clear_table(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->clear_table();
  return;
}

void o2scl_table__clear_constants(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->clear_constants();
  return;
}

void o2scl_table__sort_table(void *vptr, char *scol) {
  table<> *ptr=(table<> *)vptr;
  ptr->sort_table(scol);
  return;
}

void o2scl_table__sort_column(void *vptr, char *scol) {
  table<> *ptr=(table<> *)vptr;
  ptr->sort_column(scol);
  return;
}

void o2scl_table__average_col_roll(void *vptr, char *col_name, size_t window) {
  table<> *ptr=(table<> *)vptr;
  ptr->average_col_roll(col_name,window);
  return;
}

void o2scl_table__average_rows(void *vptr, size_t window, bool rolling) {
  table<> *ptr=(table<> *)vptr;
  ptr->average_rows(window,rolling);
  return;
}

void o2scl_table__is_valid(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->is_valid();
  return;
}

void o2scl_table__functions_columns(void *vptr, char *list) {
  table<> *ptr=(table<> *)vptr;
  ptr->functions_columns(list);
  return;
}

void o2scl_table__function_column(void *vptr, char *function, char *scol) {
  table<> *ptr=(table<> *)vptr;
  ptr->function_column(function,scol);
  return;
}

double o2scl_table__row_function(void *vptr, char *scol, size_t row) {
  table<> *ptr=(table<> *)vptr;
  double ret=ptr->row_function(scol,row);
  return ret;
}

size_t o2scl_table__function_find_row(void *vptr, char *function) {
  table<> *ptr=(table<> *)vptr;
  size_t ret=ptr->function_find_row(function);
  return ret;
}

void o2scl_table__summary(void *vptr) {
  table<> *ptr=(table<> *)vptr;
  ptr->summary();
  return;
}

void *o2scl_create_table_units_() {
  table_units<> *ptr=new table_units<>;
  return ptr;
}

void o2scl_free_table_units_(void *vptr) {
  table_units<> *ptr=(table_units<> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_table_units_(void *vsrc, void *vdest) {
  table_units<> *src=(table_units<> *)vsrc;
  table_units<> *dest=(table_units<> *)vdest;
  *dest=*src;
}

void o2scl_table_units__set_unit(void *vptr, char *col, char *unit) {
  table_units<> *ptr=(table_units<> *)vptr;
  ptr->set_unit(col,unit);
  return;
}

void *o2scl_table_units__get_unit(void *vptr, char *col) {
  table_units<> *ptr=(table_units<> *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->get_unit(col);
  return sptr;
}

void o2scl_table_units__line_of_units(void *vptr, char *unit_line) {
  table_units<> *ptr=(table_units<> *)vptr;
  ptr->line_of_units(unit_line);
  return;
}

void o2scl_table_units__remove_unit(void *vptr, char *col) {
  table_units<> *ptr=(table_units<> *)vptr;
  ptr->remove_unit(col);
  return;
}

int o2scl_table_units__convert_to_unit(void *vptr, char *col, char *unit, bool err_on_fail) {
  table_units<> *ptr=(table_units<> *)vptr;
  int ret=ptr->convert_to_unit(col,unit,err_on_fail);
  return ret;
}

void *o2scl_create_uniform_grid_() {
  uniform_grid<> *ptr=new uniform_grid<>;
  return ptr;
}

void o2scl_free_uniform_grid_(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  delete ptr;
  return;
}

size_t o2scl_uniform_grid__get_nbins(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  size_t ret=ptr->get_nbins();
  return ret;
}

size_t o2scl_uniform_grid__get_npoints(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  size_t ret=ptr->get_npoints();
  return ret;
}

bool o2scl_uniform_grid__is_log(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  bool ret=ptr->is_log();
  return ret;
}

double o2scl_uniform_grid__get_start(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  double ret=ptr->get_start();
  return ret;
}

double o2scl_uniform_grid__get_end(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  double ret=ptr->get_end();
  return ret;
}

double o2scl_uniform_grid__get_width(void *vptr) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  double ret=ptr->get_width();
  return ret;
}

double o2scl_uniform_grid__getitem(void *vptr, size_t n) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  double ret=ptr->operator[](n);
  return ret;
}

void o2scl_uniform_grid__vector(void *vptr, void *ptr_v) {
  uniform_grid<> *ptr=(uniform_grid<> *)vptr;
  std::vector<double> *v=(std::vector<double> *)ptr_v;
  ptr->vector(*v);
  return;
}

void *o2scl_create_uniform_grid_end_() {
  uniform_grid_end<> *ptr=new uniform_grid_end<>;
  return ptr;
}

void o2scl_free_uniform_grid_end_(void *vptr) {
  uniform_grid_end<> *ptr=(uniform_grid_end<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_end__init(double start, double end, size_t n_bins) {
  uniform_grid_end<> *ptr=new uniform_grid_end<>(start,end,n_bins);
  return ptr;
}

void *o2scl_create_uniform_grid_width_() {
  uniform_grid_width<> *ptr=new uniform_grid_width<>;
  return ptr;
}

void o2scl_free_uniform_grid_width_(void *vptr) {
  uniform_grid_width<> *ptr=(uniform_grid_width<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_width__init(double start, double width, size_t n_bins) {
  uniform_grid_width<> *ptr=new uniform_grid_width<>(start,width,n_bins);
  return ptr;
}

void *o2scl_create_uniform_grid_end_width_() {
  uniform_grid_end_width<> *ptr=new uniform_grid_end_width<>;
  return ptr;
}

void o2scl_free_uniform_grid_end_width_(void *vptr) {
  uniform_grid_end_width<> *ptr=(uniform_grid_end_width<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_end_width__init(double start, double end, double width) {
  uniform_grid_end_width<> *ptr=new uniform_grid_end_width<>(start,end,width);
  return ptr;
}

void *o2scl_create_uniform_grid_log_end_() {
  uniform_grid_log_end<> *ptr=new uniform_grid_log_end<>;
  return ptr;
}

void o2scl_free_uniform_grid_log_end_(void *vptr) {
  uniform_grid_log_end<> *ptr=(uniform_grid_log_end<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_log_end__init(double start, double end, size_t n_bins) {
  uniform_grid_log_end<> *ptr=new uniform_grid_log_end<>(start,end,n_bins);
  return ptr;
}

void *o2scl_create_uniform_grid_log_width_() {
  uniform_grid_log_width<> *ptr=new uniform_grid_log_width<>;
  return ptr;
}

void o2scl_free_uniform_grid_log_width_(void *vptr) {
  uniform_grid_log_width<> *ptr=(uniform_grid_log_width<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_log_width__init(double start, double width, size_t n_bins) {
  uniform_grid_log_width<> *ptr=new uniform_grid_log_width<>(start,width,n_bins);
  return ptr;
}

void *o2scl_create_uniform_grid_log_end_width_() {
  uniform_grid_log_end_width<> *ptr=new uniform_grid_log_end_width<>;
  return ptr;
}

void o2scl_free_uniform_grid_log_end_width_(void *vptr) {
  uniform_grid_log_end_width<> *ptr=(uniform_grid_log_end_width<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_uniform_grid_log_end_width__init(double start, double end, double width) {
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

void o2scl_free_ix_index(void *vptr) {
  ix_index *ptr=(ix_index *)vptr;
  delete ptr;
  return;
}

void *o2scl_ix_index_init(size_t ix) {
  ix_index *ptr=new ix_index(ix);
  return ptr;
}

void o2scl_free_ix_fixed(void *vptr) {
  ix_fixed *ptr=(ix_fixed *)vptr;
  delete ptr;
  return;
}

void *o2scl_ix_fixed_init(size_t ix, size_t ix2) {
  ix_fixed *ptr=new ix_fixed(ix,ix2);
  return ptr;
}

void o2scl_free_ix_sum(void *vptr) {
  ix_sum *ptr=(ix_sum *)vptr;
  delete ptr;
  return;
}

void *o2scl_ix_sum_init(size_t ix) {
  ix_sum *ptr=new ix_sum(ix);
  return ptr;
}

void o2scl_free_ix_trace(void *vptr) {
  ix_trace *ptr=(ix_trace *)vptr;
  delete ptr;
  return;
}

void *o2scl_ix_trace_init(size_t ix, size_t ix2) {
  ix_trace *ptr=new ix_trace(ix,ix2);
  return ptr;
}

void o2scl_free_ix_reverse(void *vptr) {
  ix_reverse *ptr=(ix_reverse *)vptr;
  delete ptr;
  return;
}

void *o2scl_ix_reverse_init(size_t ix) {
  ix_reverse *ptr=new ix_reverse(ix);
  return ptr;
}

void o2scl_free_ix_range(void *vptr) {
  ix_range *ptr=(ix_range *)vptr;
  delete ptr;
  return;
}

void *o2scl_ix_range_init(size_t ix, size_t start, size_t end) {
  ix_range *ptr=new ix_range(ix,start,end);
  return ptr;
}

void o2scl_free_ix_interp(void *vptr) {
  ix_interp *ptr=(ix_interp *)vptr;
  delete ptr;
  return;
}

void *o2scl_ix_interp_init(size_t ix, double v) {
  ix_interp *ptr=new ix_interp(ix,v);
  return ptr;
}

void o2scl_free_ix_grid(void *vptr) {
  ix_grid *ptr=(ix_grid *)vptr;
  delete ptr;
  return;
}

void *o2scl_ix_grid_init(size_t ix, double start, double end, size_t n_bins, bool log) {
  ix_grid *ptr=new ix_grid(ix,start,end,n_bins,log);
  return ptr;
}

void o2scl_free_ix_gridw(void *vptr) {
  ix_gridw *ptr=(ix_gridw *)vptr;
  delete ptr;
  return;
}

void *o2scl_ix_gridw_init(size_t ix, double start, double end, double width, bool log) {
  ix_gridw *ptr=new ix_gridw(ix,start,end,width,log);
  return ptr;
}

void *o2scl_create_tensor_() {
  tensor<> *ptr=new tensor<>;
  return ptr;
}

void o2scl_free_tensor_(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_tensor_(void *vsrc, void *vdest) {
  tensor<> *src=(tensor<> *)vsrc;
  tensor<> *dest=(tensor<> *)vdest;
  *dest=*src;
}

void o2scl_tensor__is_valid(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  ptr->is_valid();
  return;
}

void o2scl_tensor__clear(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  ptr->clear();
  return;
}

void o2scl_tensor__set(void *vptr, void *ptr_index, double val) {
  tensor<> *ptr=(tensor<> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  ptr->set(*index,val);
  return;
}

void o2scl_tensor__set_all(void *vptr, double x) {
  tensor<> *ptr=(tensor<> *)vptr;
  ptr->set_all(x);
  return;
}

void o2scl_tensor__swap_data(void *vptr, void *ptr_data) {
  tensor<> *ptr=(tensor<> *)vptr;
  std::vector<double> *data=(std::vector<double> *)ptr_data;
  ptr->swap_data(*data);
  return;
}

double o2scl_tensor__get(void *vptr, void *ptr_index) {
  tensor<> *ptr=(tensor<> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  double ret=ptr->get(*index);
  return ret;
}

void o2scl_tensor__resize(void *vptr, size_t n, void *ptr_index) {
  tensor<> *ptr=(tensor<> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  ptr->resize(n,*index);
  return;
}

size_t o2scl_tensor__get_rank(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  size_t ret=ptr->get_rank();
  return ret;
}

size_t o2scl_tensor__get_size(void *vptr, size_t i) {
  tensor<> *ptr=(tensor<> *)vptr;
  size_t ret=ptr->get_size(i);
  return ret;
}

void *o2scl_tensor__get_size_arr(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  const std::vector<size_t> *ret=&ptr->get_size_arr();
  return (void *)ret;
}

void o2scl_tensor__get_data(void *vptr, double **dptr, int *n_) {
  tensor<> *ptr=(tensor<> *)vptr;
  const std::vector<double> &r=ptr->get_data();
  *dptr=(double *)(&(r[0]));
  *n_=r.size();
  return;
}

size_t o2scl_tensor__total_size(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  size_t ret=ptr->total_size();
  return ret;
}

size_t o2scl_tensor__pack_indices(void *vptr, void *ptr_index) {
  tensor<> *ptr=(tensor<> *)vptr;
  std::vector<size_t> *index=(std::vector<size_t> *)ptr_index;
  size_t ret=ptr->pack_indices(*index);
  return ret;
}

void o2scl_tensor__unpack_index(void *vptr, size_t ix, void *ptr_index) {
  tensor<> *ptr=(tensor<> *)vptr;
  std::vector<size_t> *index=(std::vector<size_t> *)ptr_index;
  ptr->unpack_index(ix,*index);
  return;
}

double o2scl_tensor__min_value(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  double ret=ptr->min_value();
  return ret;
}

void o2scl_tensor__min_index(void *vptr, void *ptr_index) {
  tensor<> *ptr=(tensor<> *)vptr;
  std::vector<size_t> *index=(std::vector<size_t> *)ptr_index;
  ptr->min_index(*index);
  return;
}

void o2scl_tensor__min(void *vptr, void *ptr_ix, double *value) {
  tensor<> *ptr=(tensor<> *)vptr;
  std::vector<size_t> *ix=(std::vector<size_t> *)ptr_ix;
  ptr->min(*ix,*value);
  return;
}

double o2scl_tensor__max_value(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  double ret=ptr->max_value();
  return ret;
}

void o2scl_tensor__max_index(void *vptr, void *ptr_index) {
  tensor<> *ptr=(tensor<> *)vptr;
  std::vector<size_t> *index=(std::vector<size_t> *)ptr_index;
  ptr->max_index(*index);
  return;
}

void o2scl_tensor__max(void *vptr, void *ptr_ix, double *value) {
  tensor<> *ptr=(tensor<> *)vptr;
  std::vector<size_t> *ix=(std::vector<size_t> *)ptr_ix;
  ptr->max(*ix,*value);
  return;
}

void o2scl_tensor__minmax_value(void *vptr, double *min, double *max) {
  tensor<> *ptr=(tensor<> *)vptr;
  ptr->minmax_value(*min,*max);
  return;
}

void o2scl_tensor__minmax_index(void *vptr, void *ptr_min, void *ptr_max) {
  tensor<> *ptr=(tensor<> *)vptr;
  std::vector<size_t> *min=(std::vector<size_t> *)ptr_min;
  std::vector<size_t> *max=(std::vector<size_t> *)ptr_max;
  ptr->minmax_index(*min,*max);
  return;
}

void o2scl_tensor__minmax(void *vptr, void *ptr_min_ix, double *min_value, void *ptr_max_ix, double *max_value) {
  tensor<> *ptr=(tensor<> *)vptr;
  std::vector<size_t> *min_ix=(std::vector<size_t> *)ptr_min_ix;
  std::vector<size_t> *max_ix=(std::vector<size_t> *)ptr_max_ix;
  ptr->minmax(*min_ix,*min_value,*max_ix,*max_value);
  return;
}

double o2scl_tensor__total_sum(void *vptr) {
  tensor<> *ptr=(tensor<> *)vptr;
  double ret=ptr->total_sum();
  return ret;
}

void o2scl_tensor__convert_table3d_sum(void *vptr, size_t ix_x, size_t ix_y, void *ptr_tab, char *x_name, char *y_name, char *slice_name) {
  tensor<> *ptr=(tensor<> *)vptr;
  table3d *tab=(table3d *)ptr_tab;
  ptr->convert_table3d_sum(ix_x,ix_y,*tab,x_name,y_name,slice_name);
  return;
}

void *o2scl_tensor__rearrange_and_copy(void *vptr, char *spec, int verbose, bool err_on_fail) {
  tensor<> *ptr=(tensor<> *)vptr;
  tensor<> *ret=new tensor<>;
  *ret=ptr->rearrange_and_copy(spec,verbose,err_on_fail);
  return ret;
}

void *o2scl_tensor__create_size(size_t rank, void *ptr_sizes) {
  std::vector<size_t> *sizes=(std::vector<size_t> *)ptr_sizes;
  tensor<> *ptr=new tensor<>(rank,*sizes);
  return ptr;
}

void *o2scl_create_tensor_grid_() {
  tensor_grid<> *ptr=new tensor_grid<>;
  return ptr;
}

void o2scl_free_tensor_grid_(void *vptr) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_tensor_grid_(void *vsrc, void *vdest) {
  tensor_grid<> *src=(tensor_grid<> *)vsrc;
  tensor_grid<> *dest=(tensor_grid<> *)vdest;
  *dest=*src;
}

void o2scl_tensor_grid__is_valid(void *vptr) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  ptr->is_valid();
  return;
}

void o2scl_tensor_grid__set_val(void *vptr, void *ptr_grid_point, double val) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  vector<double> *grid_point=(vector<double> *)ptr_grid_point;
  ptr->set_val(*grid_point,val);
  return;
}

double o2scl_tensor_grid__get_val(void *vptr, void *ptr_grid_point) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  vector<double> *grid_point=(vector<double> *)ptr_grid_point;
  double ret=ptr->get_val(*grid_point);
  return ret;
}

bool o2scl_tensor_grid__is_grid_set(void *vptr) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  bool ret=ptr->is_grid_set();
  return ret;
}

void o2scl_tensor_grid__set_grid_packed(void *vptr, void *ptr_grid) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  vector<double> *grid=(vector<double> *)ptr_grid;
  ptr->set_grid_packed(*grid);
  return;
}

void o2scl_tensor_grid__default_grid(void *vptr) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  ptr->default_grid();
  return;
}

void o2scl_tensor_grid__set_grid_i_vec(void *vptr, size_t i, void *ptr_grid) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  vector<double> *grid=(vector<double> *)ptr_grid;
  ptr->set_grid_i_vec(i,*grid);
  return;
}

double o2scl_tensor_grid__get_grid(void *vptr, size_t i, size_t j) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  double ret=ptr->get_grid(i,j);
  return ret;
}

void o2scl_tensor_grid__set_grid(void *vptr, size_t i, size_t j, double val) {
  tensor_grid<> *ptr=(tensor_grid<> *)vptr;
  ptr->set_grid(i,j,val);
  return;
}

void *o2scl_create_tensor_int_std_vector_int_() {
  tensor<int,std::vector<int>> *ptr=new tensor<int,std::vector<int>>;
  return ptr;
}

void o2scl_free_tensor_int_std_vector_int_(void *vptr) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_tensor_int_std_vector_int_(void *vsrc, void *vdest) {
  tensor<int,std::vector<int>> *src=(tensor<int,std::vector<int>> *)vsrc;
  tensor<int,std::vector<int>> *dest=(tensor<int,std::vector<int>> *)vdest;
  *dest=*src;
}

void o2scl_tensor_int_std_vector_int__is_valid(void *vptr) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  ptr->is_valid();
  return;
}

void o2scl_tensor_int_std_vector_int__clear(void *vptr) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  ptr->clear();
  return;
}

void o2scl_tensor_int_std_vector_int__set(void *vptr, void *ptr_index, int val) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  ptr->set(*index,val);
  return;
}

void o2scl_tensor_int_std_vector_int__set_all(void *vptr, int x) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  ptr->set_all(x);
  return;
}

int o2scl_tensor_int_std_vector_int__get(void *vptr, void *ptr_index) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  int ret=ptr->get(*index);
  return ret;
}

void o2scl_tensor_int_std_vector_int__resize(void *vptr, size_t n, void *ptr_index) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  ptr->resize(n,*index);
  return;
}

size_t o2scl_tensor_int_std_vector_int__get_rank(void *vptr) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  size_t ret=ptr->get_rank();
  return ret;
}

size_t o2scl_tensor_int_std_vector_int__get_size(void *vptr, size_t i) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  size_t ret=ptr->get_size(i);
  return ret;
}

void *o2scl_tensor_int_std_vector_int__get_data(void *vptr) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  const std::vector<int> *ret=&ptr->get_data();
  return (void *)ret;
}

size_t o2scl_tensor_int_std_vector_int__total_size(void *vptr) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  size_t ret=ptr->total_size();
  return ret;
}

int o2scl_tensor_int_std_vector_int__min_value(void *vptr) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  int ret=ptr->min_value();
  return ret;
}

int o2scl_tensor_int_std_vector_int__max_value(void *vptr) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  int ret=ptr->max_value();
  return ret;
}

int o2scl_tensor_int_std_vector_int__total_sum(void *vptr) {
  tensor<int,std::vector<int>> *ptr=(tensor<int,std::vector<int>> *)vptr;
  int ret=ptr->total_sum();
  return ret;
}

void *o2scl_tensor_int_std_vector_int__create_size(size_t rank, void *ptr_sizes) {
  std::vector<size_t> *sizes=(std::vector<size_t> *)ptr_sizes;
  tensor<int,std::vector<int>> *ptr=new tensor<int,std::vector<int>>(rank,*sizes);
  return ptr;
}

void *o2scl_create_tensor_size_t_std_vector_size_t_() {
  tensor<size_t,std::vector<size_t>> *ptr=new tensor<size_t,std::vector<size_t>>;
  return ptr;
}

void o2scl_free_tensor_size_t_std_vector_size_t_(void *vptr) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_tensor_size_t_std_vector_size_t_(void *vsrc, void *vdest) {
  tensor<size_t,std::vector<size_t>> *src=(tensor<size_t,std::vector<size_t>> *)vsrc;
  tensor<size_t,std::vector<size_t>> *dest=(tensor<size_t,std::vector<size_t>> *)vdest;
  *dest=*src;
}

void o2scl_tensor_size_t_std_vector_size_t__is_valid(void *vptr) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  ptr->is_valid();
  return;
}

void o2scl_tensor_size_t_std_vector_size_t__clear(void *vptr) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  ptr->clear();
  return;
}

void o2scl_tensor_size_t_std_vector_size_t__set(void *vptr, void *ptr_index, size_t val) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  ptr->set(*index,val);
  return;
}

void o2scl_tensor_size_t_std_vector_size_t__set_all(void *vptr, size_t x) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  ptr->set_all(x);
  return;
}

int o2scl_tensor_size_t_std_vector_size_t__get(void *vptr, void *ptr_index) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  int ret=ptr->get(*index);
  return ret;
}

void o2scl_tensor_size_t_std_vector_size_t__resize(void *vptr, size_t n, void *ptr_index) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  vector<size_t> *index=(vector<size_t> *)ptr_index;
  ptr->resize(n,*index);
  return;
}

size_t o2scl_tensor_size_t_std_vector_size_t__get_rank(void *vptr) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  size_t ret=ptr->get_rank();
  return ret;
}

size_t o2scl_tensor_size_t_std_vector_size_t__get_size(void *vptr, size_t i) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  size_t ret=ptr->get_size(i);
  return ret;
}

void *o2scl_tensor_size_t_std_vector_size_t__get_data(void *vptr) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  const std::vector<size_t> *ret=&ptr->get_data();
  return (void *)ret;
}

size_t o2scl_tensor_size_t_std_vector_size_t__total_size(void *vptr) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  size_t ret=ptr->total_size();
  return ret;
}

size_t o2scl_tensor_size_t_std_vector_size_t__min_value(void *vptr) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  size_t ret=ptr->min_value();
  return ret;
}

size_t o2scl_tensor_size_t_std_vector_size_t__max_value(void *vptr) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  size_t ret=ptr->max_value();
  return ret;
}

size_t o2scl_tensor_size_t_std_vector_size_t__total_sum(void *vptr) {
  tensor<size_t,std::vector<size_t>> *ptr=(tensor<size_t,std::vector<size_t>> *)vptr;
  size_t ret=ptr->total_sum();
  return ret;
}

void *o2scl_tensor_size_t_std_vector_size_t__create_size(size_t rank, void *ptr_sizes) {
  std::vector<size_t> *sizes=(std::vector<size_t> *)ptr_sizes;
  tensor<size_t,std::vector<size_t>> *ptr=new tensor<size_t,std::vector<size_t>>(rank,*sizes);
  return ptr;
}

void *o2scl_create_find_constants_const_entry() {
  find_constants::const_entry *ptr=new find_constants::const_entry;
  return ptr;
}

void o2scl_free_find_constants_const_entry(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  delete ptr;
  return;
}

void *o2scl_find_constants_const_entry_get_names(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return (void *)(&(ptr->names));
}

void o2scl_find_constants_const_entry_set_names(void *vptr, void *p_v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  std::vector<std::string> *p_tsot=(std::vector<std::string> *)p_v;
  ptr->names=*(p_tsot);
  return;
}

void *o2scl_find_constants_const_entry_get_unit(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->unit;
  return sptr;
}

void o2scl_find_constants_const_entry_set_unit(void *vptr, void *p_v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->unit=*(p_tsot);
  return;
}

int o2scl_find_constants_const_entry_get_unit_flag(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return ptr->unit_flag;
}

void o2scl_find_constants_const_entry_set_unit_flag(void *vptr, int v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  ptr->unit_flag=v;
  return;
}

double o2scl_find_constants_const_entry_get_val(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return ptr->val;
}

void o2scl_find_constants_const_entry_set_val(void *vptr, double v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  ptr->val=v;
  return;
}

void *o2scl_find_constants_const_entry_get_source(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->source;
  return sptr;
}

void o2scl_find_constants_const_entry_set_source(void *vptr, void *p_v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->source=*(p_tsot);
  return;
}

int o2scl_find_constants_const_entry_get_m(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return ptr->m;
}

void o2scl_find_constants_const_entry_set_m(void *vptr, int v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  ptr->m=v;
  return;
}

int o2scl_find_constants_const_entry_get_k(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return ptr->k;
}

void o2scl_find_constants_const_entry_set_k(void *vptr, int v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  ptr->k=v;
  return;
}

int o2scl_find_constants_const_entry_get_s(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return ptr->s;
}

void o2scl_find_constants_const_entry_set_s(void *vptr, int v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  ptr->s=v;
  return;
}

int o2scl_find_constants_const_entry_get_K(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return ptr->K;
}

void o2scl_find_constants_const_entry_set_K(void *vptr, int v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  ptr->K=v;
  return;
}

int o2scl_find_constants_const_entry_get_A(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return ptr->A;
}

void o2scl_find_constants_const_entry_set_A(void *vptr, int v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  ptr->A=v;
  return;
}

int o2scl_find_constants_const_entry_get_mol(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return ptr->mol;
}

void o2scl_find_constants_const_entry_set_mol(void *vptr, int v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  ptr->mol=v;
  return;
}

int o2scl_find_constants_const_entry_get_cd(void *vptr) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  return ptr->cd;
}

void o2scl_find_constants_const_entry_set_cd(void *vptr, int v) {
  find_constants::const_entry *ptr=(find_constants::const_entry *)vptr;
  ptr->cd=v;
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

void o2scl_find_constants_output_list_cout(void *vptr) {
  find_constants *ptr=(find_constants *)vptr;
  ptr->output_list_cout();
  return;
}

void o2scl_find_constants_add_constant(void *vptr, void *ptr_f, int verbose) {
  find_constants *ptr=(find_constants *)vptr;
  find_constants::const_entry *f=(find_constants::const_entry *)ptr_f;
  ptr->add_constant(*f,verbose);
  return;
}

void o2scl_find_constants_del_constant(void *vptr, void *ptr_name, int verbose) {
  find_constants *ptr=(find_constants *)vptr;
  std::string *name=(std::string *)ptr_name;
  ptr->del_constant(*name,verbose);
  return;
}

void *o2scl_create_convert_units_der_unit() {
  convert_units<>::der_unit *ptr=new convert_units<>::der_unit;
  return ptr;
}

void o2scl_free_convert_units_der_unit(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  delete ptr;
  return;
}

void *o2scl_convert_units_der_unit_get_label(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->label;
  return sptr;
}

void o2scl_convert_units_der_unit_set_label(void *vptr, void *p_v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->label=*(p_tsot);
  return;
}

int o2scl_convert_units_der_unit_get_m(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  return ptr->m;
}

void o2scl_convert_units_der_unit_set_m(void *vptr, int v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  ptr->m=v;
  return;
}

int o2scl_convert_units_der_unit_get_k(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  return ptr->k;
}

void o2scl_convert_units_der_unit_set_k(void *vptr, int v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  ptr->k=v;
  return;
}

int o2scl_convert_units_der_unit_get_s(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  return ptr->s;
}

void o2scl_convert_units_der_unit_set_s(void *vptr, int v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  ptr->s=v;
  return;
}

int o2scl_convert_units_der_unit_get_K(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  return ptr->K;
}

void o2scl_convert_units_der_unit_set_K(void *vptr, int v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  ptr->K=v;
  return;
}

int o2scl_convert_units_der_unit_get_A(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  return ptr->A;
}

void o2scl_convert_units_der_unit_set_A(void *vptr, int v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  ptr->A=v;
  return;
}

int o2scl_convert_units_der_unit_get_mol(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  return ptr->mol;
}

void o2scl_convert_units_der_unit_set_mol(void *vptr, int v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  ptr->mol=v;
  return;
}

int o2scl_convert_units_der_unit_get_cd(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  return ptr->cd;
}

void o2scl_convert_units_der_unit_set_cd(void *vptr, int v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  ptr->cd=v;
  return;
}

double o2scl_convert_units_der_unit_get_val(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  return ptr->val;
}

void o2scl_convert_units_der_unit_set_val(void *vptr, double v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  ptr->val=v;
  return;
}

void *o2scl_convert_units_der_unit_get_name(void *vptr) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->name;
  return sptr;
}

void o2scl_convert_units_der_unit_set_name(void *vptr, void *p_v) {
  convert_units<>::der_unit *ptr=(convert_units<>::der_unit *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->name=*(p_tsot);
  return;
}

void *o2scl_create_convert_units_() {
  convert_units<> *ptr=new convert_units<>;
  return ptr;
}

void o2scl_free_convert_units_(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  delete ptr;
  return;
}

int o2scl_convert_units__get_verbose(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  return ptr->verbose;
}

void o2scl_convert_units__set_verbose(void *vptr, int v) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->verbose=v;
  return;
}

bool o2scl_convert_units__get_err_on_fail(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  return ptr->err_on_fail;
}

void o2scl_convert_units__set_err_on_fail(void *vptr, bool v) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->err_on_fail=v;
  return;
}

bool o2scl_convert_units__get_combine_two_conv(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  return ptr->combine_two_conv;
}

void o2scl_convert_units__set_combine_two_conv(void *vptr, bool v) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->combine_two_conv=v;
  return;
}

double o2scl_convert_units__convert(void *vptr, char *frm, char *to, double val) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  double ret=ptr->convert(frm,to,val);
  return ret;
}

int o2scl_convert_units__convert_ret(void *vptr, char *frm, char *to, double val, double converted) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  int ret=ptr->convert_ret(frm,to,val,converted);
  return ret;
}

void o2scl_convert_units__del_unit(void *vptr, void *ptr_name) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  std::string *name=(std::string *)ptr_name;
  ptr->del_unit(*name);
  return;
}

void o2scl_convert_units__add_unit(void *vptr, void *ptr_d) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  convert_units<>::der_unit *d=(convert_units<>::der_unit *)ptr_d;
  ptr->add_unit(*d);
  return;
}

void o2scl_convert_units__set_natural_units(void *vptr, bool c_is_one, bool hbar_is_one, bool kb_is_one) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->set_natural_units(c_is_one,hbar_is_one,kb_is_one);
  return;
}

int o2scl_convert_units__is_in_cache(void *vptr, char *frm, char *to) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  int ret=ptr->is_in_cache(frm,to);
  return ret;
}

int o2scl_convert_units__remove_cache(void *vptr, char *frm, char *to) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  int ret=ptr->remove_cache(frm,to);
  return ret;
}

void o2scl_convert_units__clear_cache(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->clear_cache();
  return;
}

void o2scl_convert_units__test_unique(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->test_unique();
  return;
}

void o2scl_convert_units__print_cache(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->print_cache();
  return;
}

void o2scl_convert_units__print_units_cout(void *vptr) {
  convert_units<> *ptr=(convert_units<> *)vptr;
  ptr->print_units_cout();
  return;
}

void *o2scl_create_columnify() {
  columnify *ptr=new columnify;
  return ptr;
}

void o2scl_free_columnify(void *vptr) {
  columnify *ptr=(columnify *)vptr;
  delete ptr;
  return;
}

int o2scl_columnify_get_align_left(void *vptr) {
  columnify *ptr=(columnify *)vptr;
  return ptr->align_left;
}


int o2scl_columnify_get_align_right(void *vptr) {
  columnify *ptr=(columnify *)vptr;
  return ptr->align_right;
}


int o2scl_columnify_get_align_lmid(void *vptr) {
  columnify *ptr=(columnify *)vptr;
  return ptr->align_lmid;
}


int o2scl_columnify_get_align_rmid(void *vptr) {
  columnify *ptr=(columnify *)vptr;
  return ptr->align_rmid;
}


int o2scl_columnify_get_align_dp(void *vptr) {
  columnify *ptr=(columnify *)vptr;
  return ptr->align_dp;
}


int o2scl_columnify_get_align_lnum(void *vptr) {
  columnify *ptr=(columnify *)vptr;
  return ptr->align_lnum;
}


void *o2scl_create_format_float() {
  format_float *ptr=new format_float;
  return ptr;
}

void o2scl_free_format_float(void *vptr) {
  format_float *ptr=(format_float *)vptr;
  delete ptr;
  return;
}

void o2scl_format_float_set_sig_figs(void *vptr, size_t sig_figs) {
  format_float *ptr=(format_float *)vptr;
  ptr->set_sig_figs(sig_figs);
  return;
}

void o2scl_format_float_set_exp_limits(void *vptr, int min, int max) {
  format_float *ptr=(format_float *)vptr;
  ptr->set_exp_limits(min,max);
  return;
}

void o2scl_format_float_set_pad_zeros(void *vptr, bool pad) {
  format_float *ptr=(format_float *)vptr;
  ptr->set_pad_zeros(pad);
  return;
}

void o2scl_format_float_set_dec_point(void *vptr, char *dec_point) {
  format_float *ptr=(format_float *)vptr;
  ptr->set_dec_point(dec_point);
  return;
}

void o2scl_format_float_set_exp_digits(void *vptr, int d) {
  format_float *ptr=(format_float *)vptr;
  ptr->set_exp_digits(d);
  return;
}

void o2scl_format_float_html_mode(void *vptr) {
  format_float *ptr=(format_float *)vptr;
  ptr->html_mode();
  return;
}

void o2scl_format_float_latex_mode(void *vptr) {
  format_float *ptr=(format_float *)vptr;
  ptr->latex_mode();
  return;
}

void o2scl_format_float_c_mode(void *vptr) {
  format_float *ptr=(format_float *)vptr;
  ptr->c_mode();
  return;
}

void *o2scl_format_float_convert(void *vptr, double x, bool debug) {
  format_float *ptr=(format_float *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->convert(x,debug);
  return sptr;
}

void *o2scl_create_interp_std_vector_double_() {
  interp<std::vector<double>> *ptr=new interp<std::vector<double>>;
  return ptr;
}

void o2scl_free_interp_std_vector_double_(void *vptr) {
  interp<std::vector<double>> *ptr=(interp<std::vector<double>> *)vptr;
  delete ptr;
  return;
}

double o2scl_interp_std_vector_double__eval(void *vptr, double x0, size_t n, void *ptr_x, void *ptr_y) {
  interp<std::vector<double>> *ptr=(interp<std::vector<double>> *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  double ret=ptr->eval(x0,n,*x,*y);
  return ret;
}

double o2scl_interp_std_vector_double__deriv(void *vptr, double x0, size_t n, void *ptr_x, void *ptr_y) {
  interp<std::vector<double>> *ptr=(interp<std::vector<double>> *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  double ret=ptr->deriv(x0,n,*x,*y);
  return ret;
}

double o2scl_interp_std_vector_double__deriv2(void *vptr, double x0, size_t n, void *ptr_x, void *ptr_y) {
  interp<std::vector<double>> *ptr=(interp<std::vector<double>> *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  double ret=ptr->deriv2(x0,n,*x,*y);
  return ret;
}

double o2scl_interp_std_vector_double__integ(void *vptr, double x1, double x2, size_t n, void *ptr_x, void *ptr_y) {
  interp<std::vector<double>> *ptr=(interp<std::vector<double>> *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  double ret=ptr->integ(x1,x2,n,*x,*y);
  return ret;
}

void o2scl_interp_std_vector_double__set_type(void *vptr, int interp_type) {
  interp<std::vector<double>> *ptr=(interp<std::vector<double>> *)vptr;
  ptr->set_type(interp_type);
  return;
}

void *o2scl_create_interp_vec_std_vector_double_() {
  interp_vec<std::vector<double>> *ptr=new interp_vec<std::vector<double>>;
  return ptr;
}

void o2scl_free_interp_vec_std_vector_double_(void *vptr) {
  interp_vec<std::vector<double>> *ptr=(interp_vec<std::vector<double>> *)vptr;
  delete ptr;
  return;
}

void o2scl_interp_vec_std_vector_double__set(void *vptr, size_t n, void *ptr_x, void *ptr_y, int interp_type) {
  interp_vec<std::vector<double>> *ptr=(interp_vec<std::vector<double>> *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  ptr->set(n,*x,*y,interp_type);
  return;
}

void o2scl_interp_vec_std_vector_double__clear(void *vptr) {
  interp_vec<std::vector<double>> *ptr=(interp_vec<std::vector<double>> *)vptr;
  ptr->clear();
  return;
}

double o2scl_interp_vec_std_vector_double__eval(void *vptr, double x0) {
  interp_vec<std::vector<double>> *ptr=(interp_vec<std::vector<double>> *)vptr;
  double ret=ptr->eval(x0);
  return ret;
}

double o2scl_interp_vec_std_vector_double__deriv(void *vptr, double x0) {
  interp_vec<std::vector<double>> *ptr=(interp_vec<std::vector<double>> *)vptr;
  double ret=ptr->deriv(x0);
  return ret;
}

double o2scl_interp_vec_std_vector_double__deriv2(void *vptr, double x0) {
  interp_vec<std::vector<double>> *ptr=(interp_vec<std::vector<double>> *)vptr;
  double ret=ptr->deriv2(x0);
  return ret;
}

double o2scl_interp_vec_std_vector_double__integ(void *vptr, double x1, double x2) {
  interp_vec<std::vector<double>> *ptr=(interp_vec<std::vector<double>> *)vptr;
  double ret=ptr->integ(x1,x2);
  return ret;
}

void *o2scl_create_interp_krige_optim_std_vector_double_() {
  interp_krige_optim<std::vector<double>> *ptr=new interp_krige_optim<std::vector<double>>;
  return ptr;
}

void o2scl_free_interp_krige_optim_std_vector_double_(void *vptr) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  delete ptr;
  return;
}

int o2scl_interp_krige_optim_std_vector_double__get_verbose(void *vptr) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  return ptr->verbose;
}

void o2scl_interp_krige_optim_std_vector_double__set_verbose(void *vptr, int v) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  ptr->verbose=v;
  return;
}

size_t o2scl_interp_krige_optim_std_vector_double__get_mode(void *vptr) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  return ptr->mode;
}

void o2scl_interp_krige_optim_std_vector_double__set_mode(void *vptr, size_t v) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  ptr->mode=v;
  return;
}

size_t o2scl_interp_krige_optim_std_vector_double__get_nlen(void *vptr) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  return ptr->nlen;
}

void o2scl_interp_krige_optim_std_vector_double__set_nlen(void *vptr, size_t v) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  ptr->nlen=v;
  return;
}

bool o2scl_interp_krige_optim_std_vector_double__get_full_min(void *vptr) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  return ptr->full_min;
}

void o2scl_interp_krige_optim_std_vector_double__set_full_min(void *vptr, bool v) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  ptr->full_min=v;
  return;
}

double o2scl_interp_krige_optim_std_vector_double__eval(void *vptr, double x0) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  double ret=ptr->eval(x0);
  return ret;
}

double o2scl_interp_krige_optim_std_vector_double__deriv(void *vptr, double x0) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  double ret=ptr->deriv(x0);
  return ret;
}

double o2scl_interp_krige_optim_std_vector_double__deriv2(void *vptr, double x0) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  double ret=ptr->deriv2(x0);
  return ret;
}

double o2scl_interp_krige_optim_std_vector_double__sigma(void *vptr, double x0) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  double ret=ptr->sigma(x0);
  return ret;
}

double o2scl_interp_krige_optim_std_vector_double__sample(void *vptr, double x0) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  double ret=ptr->sample(x0);
  return ret;
}

void o2scl_interp_krige_optim_std_vector_double__sample_vec(void *vptr, void *ptr_x, void *ptr_y) {
  interp_krige_optim<std::vector<double>> *ptr=(interp_krige_optim<std::vector<double>> *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  ptr->sample_vec(*x,*y);
  return;
}

void *o2scl_create_gen_test_number_double_() {
  gen_test_number<double> *ptr=new gen_test_number<double>;
  return ptr;
}

void o2scl_free_gen_test_number_double_(void *vptr) {
  gen_test_number<double> *ptr=(gen_test_number<double> *)vptr;
  delete ptr;
  return;
}

void o2scl_gen_test_number_double__reset(void *vptr) {
  gen_test_number<double> *ptr=(gen_test_number<double> *)vptr;
  ptr->reset();
  return;
}

void o2scl_gen_test_number_double__set_radix(void *vptr, double r) {
  gen_test_number<double> *ptr=(gen_test_number<double> *)vptr;
  ptr->set_radix(r);
  return;
}

double o2scl_gen_test_number_double__gen(void *vptr) {
  gen_test_number<double> *ptr=(gen_test_number<double> *)vptr;
  double ret=ptr->gen();
  return ret;
}

void o2scl_free_funct_string(void *vptr) {
  funct_string *ptr=(funct_string *)vptr;
  delete ptr;
  return;
}

int o2scl_funct_string_set_parm(void *vptr, char *name, double val) {
  funct_string *ptr=(funct_string *)vptr;
  int ret=ptr->set_parm(name,val);
  return ret;
}

double o2scl_funct_string_getitem(void *vptr, double x) {
  funct_string *ptr=(funct_string *)vptr;
  double ret=ptr->operator()(x);
  return ret;
}

void *o2scl_funct_string_init(char *expr, char *var) {
  funct_string *ptr=new funct_string(expr,var);
  return ptr;
}

void *o2scl_create_comm_option_s() {
  comm_option_s *ptr=new comm_option_s;
  return ptr;
}

void o2scl_free_comm_option_s(void *vptr) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  delete ptr;
  return;
}

char o2scl_comm_option_s_get_shrt(void *vptr) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  return ptr->shrt;
}

void o2scl_comm_option_s_set_shrt(void *vptr, char v) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  ptr->shrt=v;
  return;
}

void *o2scl_comm_option_s_get_lng(void *vptr) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->lng;
  return sptr;
}

void o2scl_comm_option_s_set_lng(void *vptr, void *p_v) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->lng=*(p_tsot);
  return;
}

void *o2scl_comm_option_s_get_desc(void *vptr) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->desc;
  return sptr;
}

void o2scl_comm_option_s_set_desc(void *vptr, void *p_v) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->desc=*(p_tsot);
  return;
}

int o2scl_comm_option_s_get_min_parms(void *vptr) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  return ptr->min_parms;
}

void o2scl_comm_option_s_set_min_parms(void *vptr, int v) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  ptr->min_parms=v;
  return;
}

int o2scl_comm_option_s_get_max_parms(void *vptr) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  return ptr->max_parms;
}

void o2scl_comm_option_s_set_max_parms(void *vptr, int v) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  ptr->max_parms=v;
  return;
}

void *o2scl_comm_option_s_get_parm_desc(void *vptr) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->parm_desc;
  return sptr;
}

void o2scl_comm_option_s_set_parm_desc(void *vptr, void *p_v) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->parm_desc=*(p_tsot);
  return;
}

void *o2scl_comm_option_s_get_help(void *vptr) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->help;
  return sptr;
}

void o2scl_comm_option_s_set_help(void *vptr, void *p_v) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->help=*(p_tsot);
  return;
}

int o2scl_comm_option_s_get_type(void *vptr) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  return ptr->type;
}

void o2scl_comm_option_s_set_type(void *vptr, int v) {
  comm_option_s *ptr=(comm_option_s *)vptr;
  ptr->type=v;
  return;
}

void *o2scl_create_cmd_line_arg() {
  cmd_line_arg *ptr=new cmd_line_arg;
  return ptr;
}

void o2scl_free_cmd_line_arg(void *vptr) {
  cmd_line_arg *ptr=(cmd_line_arg *)vptr;
  delete ptr;
  return;
}

void *o2scl_cmd_line_arg_get_arg(void *vptr) {
  cmd_line_arg *ptr=(cmd_line_arg *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->arg;
  return sptr;
}

void o2scl_cmd_line_arg_set_arg(void *vptr, void *p_v) {
  cmd_line_arg *ptr=(cmd_line_arg *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->arg=*(p_tsot);
  return;
}

bool o2scl_cmd_line_arg_get_is_option(void *vptr) {
  cmd_line_arg *ptr=(cmd_line_arg *)vptr;
  return ptr->is_option;
}

void o2scl_cmd_line_arg_set_is_option(void *vptr, bool v) {
  cmd_line_arg *ptr=(cmd_line_arg *)vptr;
  ptr->is_option=v;
  return;
}

bool o2scl_cmd_line_arg_get_is_valid(void *vptr) {
  cmd_line_arg *ptr=(cmd_line_arg *)vptr;
  return ptr->is_valid;
}

void o2scl_cmd_line_arg_set_is_valid(void *vptr, bool v) {
  cmd_line_arg *ptr=(cmd_line_arg *)vptr;
  ptr->is_valid=v;
  return;
}

void *o2scl_cmd_line_arg_get_parms(void *vptr) {
  cmd_line_arg *ptr=(cmd_line_arg *)vptr;
  return (void *)(&(ptr->parms));
}

void o2scl_cmd_line_arg_set_parms(void *vptr, void *p_v) {
  cmd_line_arg *ptr=(cmd_line_arg *)vptr;
  std::vector<std::string> *p_tsot=(std::vector<std::string> *)p_v;
  ptr->parms=*(p_tsot);
  return;
}

void *o2scl_create_cli() {
  cli *ptr=new cli;
  return ptr;
}

void o2scl_free_cli(void *vptr) {
  cli *ptr=(cli *)vptr;
  delete ptr;
  return;
}

bool o2scl_cli_get_sync_verbose(void *vptr) {
  cli *ptr=(cli *)vptr;
  return ptr->sync_verbose;
}

void o2scl_cli_set_sync_verbose(void *vptr, bool v) {
  cli *ptr=(cli *)vptr;
  ptr->sync_verbose=v;
  return;
}

bool o2scl_cli_get_gnu_intro(void *vptr) {
  cli *ptr=(cli *)vptr;
  return ptr->gnu_intro;
}

void o2scl_cli_set_gnu_intro(void *vptr, bool v) {
  cli *ptr=(cli *)vptr;
  ptr->gnu_intro=v;
  return;
}

void *o2scl_cli_get_desc(void *vptr) {
  cli *ptr=(cli *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->desc;
  return sptr;
}

void o2scl_cli_set_desc(void *vptr, void *p_v) {
  cli *ptr=(cli *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->desc=*(p_tsot);
  return;
}

void *o2scl_cli_get_cmd_name(void *vptr) {
  cli *ptr=(cli *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->cmd_name;
  return sptr;
}

void o2scl_cli_set_cmd_name(void *vptr, void *p_v) {
  cli *ptr=(cli *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->cmd_name=*(p_tsot);
  return;
}

void *o2scl_cli_get_addl_help_cmd(void *vptr) {
  cli *ptr=(cli *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->addl_help_cmd;
  return sptr;
}

void o2scl_cli_set_addl_help_cmd(void *vptr, void *p_v) {
  cli *ptr=(cli *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->addl_help_cmd=*(p_tsot);
  return;
}

void *o2scl_cli_get_addl_help_cli(void *vptr) {
  cli *ptr=(cli *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->addl_help_cli;
  return sptr;
}

void o2scl_cli_set_addl_help_cli(void *vptr, void *p_v) {
  cli *ptr=(cli *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->addl_help_cli=*(p_tsot);
  return;
}

int o2scl_cli_set_verbose(void *vptr, int v) {
  cli *ptr=(cli *)vptr;
  int ret=ptr->set_verbose(v);
  return ret;
}

void *o2scl_create_shared_ptr_table_units_() {
  std::shared_ptr<table_units<> > *ptr=new std::shared_ptr<table_units<> >(new table_units<>);
  return ptr;
}

void o2scl_free_shared_ptr_table_units_(void *vptr) {
  std::shared_ptr<table_units<> > *ptr=(std::shared_ptr<table_units<> > *)vptr;
  delete ptr;
}

void *o2scl_shared_ptr_table_units__ptr(void *vp) {
  std::shared_ptr<table_units<> > *p=(std::shared_ptr<table_units<> > *)vp;
  table_units<> *ref=p->get();
  return ref;
}

double o2scl_fermi_function_wrapper(double E, double mu, double T, double limit) {
  double ret=fermi_function(E,mu,T,limit);
  return ret;
}

double o2scl_bose_function_wrapper(double E, double mu, double T, double limit) {
  double ret=bose_function(E,mu,T,limit);
  return ret;
}

double o2scl_quadratic_extremum_x_double__wrapper(double x1, double x2, double x3, double y1, double y2, double y3) {
  double ret=quadratic_extremum_x<double>(x1,x2,x3,y1,y2,y3);
  return ret;
}

double o2scl_quadratic_extremum_y_double__wrapper(double x1, double x2, double x3, double y1, double y2, double y3) {
  double ret=quadratic_extremum_y<double>(x1,x2,x3,y1,y2,y3);
  return ret;
}

void o2scl_screenify_vector_std_string__wrapper(size_t nin, void *ptr_in_cols, void *ptr_out_cols, size_t max_size) {
  vector<std::string> *in_cols=(vector<std::string> *)ptr_in_cols;
  vector<std::string> *out_cols=(vector<std::string> *)ptr_out_cols;
  screenify<vector<std::string>>(nin,*in_cols,*out_cols,max_size);
  return;
}

bool o2scl_file_exists_wrapper(char *fname) {
  bool ret=file_exists(fname);
  return ret;
}

void o2scl_RGBtoHSV_wrapper(double r, double g, double b, void *ptr_h, void *ptr_s, void *ptr_v) {
  double *h=(double *)ptr_h;
  double *s=(double *)ptr_s;
  double *v=(double *)ptr_v;
  RGBtoHSV(r,g,b,*h,*s,*v);
  return;
}

void o2scl_HSVtoRGB_wrapper(double h, double s, double v, void *ptr_r, void *ptr_g, void *ptr_b) {
  double *r=(double *)ptr_r;
  double *g=(double *)ptr_g;
  double *b=(double *)ptr_b;
  HSVtoRGB(h,s,v,*r,*g,*b);
  return;
}

void o2scl_wordexp_single_file_wrapper(void *&ptr_fname) {
  std::string *fname=new std::string;
  wordexp_single_file(*fname);
  ptr_fname=(void *)fname;
  return;
}

void o2scl_wordexp_wrapper_wrapper(char *word, void *ptr_matches) {
  std::vector<std::string> *matches=(std::vector<std::string> *)ptr_matches;
  wordexp_wrapper(word,*matches);
  return;
}

double o2scl_function_to_double_wrapper(char *s, int verbose) {
  double ret=function_to_double(s,verbose);
  return ret;
}

int o2scl_function_to_double_nothrow_wrapper(char *s, void *ptr_result, int verbose) {
  double *result=(double *)ptr_result;
  int ret=function_to_double_nothrow(s,*result,verbose);
  return ret;
}

int o2scl_string_to_uint_list_vector_size_t__wrapper(void *&ptr_x, void *ptr_list) {
  std::string *x=new std::string;
  vector<size_t> *list=(vector<size_t> *)ptr_list;
  int ret=string_to_uint_list<vector<size_t>>(*x,*list);
  return ret;
}

size_t o2scl_vector_level_count_std_vector_double_std_vector_double__wrapper(double level, size_t n, void *ptr_x, void *ptr_y) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  size_t ret=vector_level_count<std::vector<double>,std::vector<double>>(level,n,*x,*y);
  return ret;
}

void o2scl_vector_deriv_interp_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_v, void *ptr_dv, size_t interp_type) {
  std::vector<double> *v=(std::vector<double> *)ptr_v;
  std::vector<double> *dv=(std::vector<double> *)ptr_dv;
  vector_deriv_interp<std::vector<double>,std::vector<double>>(n,*v,*dv,interp_type);
  return;
}

void o2scl_vector_deriv2_interp_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_v, void *ptr_dv, size_t interp_type) {
  std::vector<double> *v=(std::vector<double> *)ptr_v;
  std::vector<double> *dv=(std::vector<double> *)ptr_dv;
  vector_deriv2_interp<std::vector<double>,std::vector<double>>(n,*v,*dv,interp_type);
  return;
}

void o2scl_vector_deriv_xy_interp_std_vector_double_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_vx, void *ptr_vy, void *ptr_dv, size_t interp_type) {
  std::vector<double> *vx=(std::vector<double> *)ptr_vx;
  std::vector<double> *vy=(std::vector<double> *)ptr_vy;
  std::vector<double> *dv=(std::vector<double> *)ptr_dv;
  vector_deriv_xy_interp<std::vector<double>,std::vector<double>,std::vector<double>>(n,*vx,*vy,*dv,interp_type);
  return;
}

void o2scl_vector_deriv2_xy_interp_std_vector_double_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_vx, void *ptr_vy, void *ptr_dv, size_t interp_type) {
  std::vector<double> *vx=(std::vector<double> *)ptr_vx;
  std::vector<double> *vy=(std::vector<double> *)ptr_vy;
  std::vector<double> *dv=(std::vector<double> *)ptr_dv;
  vector_deriv2_xy_interp<std::vector<double>,std::vector<double>,std::vector<double>>(n,*vx,*vy,*dv,interp_type);
  return;
}

double o2scl_vector_integ_interp_std_vector_double__wrapper(size_t n, void *ptr_vx, size_t interp_type) {
  std::vector<double> *vx=(std::vector<double> *)ptr_vx;
  double ret=vector_integ_interp<std::vector<double>>(n,*vx,interp_type);
  return ret;
}

double o2scl_vector_integ_xy_interp_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_vx, void *ptr_vy, size_t interp_type) {
  std::vector<double> *vx=(std::vector<double> *)ptr_vx;
  std::vector<double> *vy=(std::vector<double> *)ptr_vy;
  double ret=vector_integ_xy_interp<std::vector<double>,std::vector<double>>(n,*vx,*vy,interp_type);
  return ret;
}

double o2scl_vector_integ_ul_interp_std_vector_double__wrapper(size_t n, double x2, void *ptr_v, size_t interp_type) {
  std::vector<double> *v=(std::vector<double> *)ptr_v;
  double ret=vector_integ_ul_interp<std::vector<double>>(n,x2,*v,interp_type);
  return ret;
}

double o2scl_vector_integ_ul_xy_interp_std_vector_double_std_vector_double__wrapper(size_t n, double x2, void *ptr_vx, void *ptr_vy, size_t interp_type) {
  std::vector<double> *vx=(std::vector<double> *)ptr_vx;
  std::vector<double> *vy=(std::vector<double> *)ptr_vy;
  double ret=vector_integ_ul_xy_interp<std::vector<double>,std::vector<double>>(n,x2,*vx,*vy,interp_type);
  return ret;
}

void o2scl_vector_find_level_std_vector_double_std_vector_double__wrapper(double level, size_t n, void *ptr_x, void *ptr_y, void *ptr_locs) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  std::vector<double> *locs=(std::vector<double> *)ptr_locs;
  vector_find_level<std::vector<double>,std::vector<double>>(level,n,*x,*y,*locs);
  return;
}

void o2scl_vector_invert_enclosed_sum_std_vector_double_std_vector_double__wrapper(double sum, size_t n, void *ptr_x, void *ptr_y, void *ptr_lev, int boundaries, int verbose, bool err_on_fail) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  double *lev=(double *)ptr_lev;
  vector_invert_enclosed_sum<std::vector<double>,std::vector<double>>(sum,n,*x,*y,*lev,boundaries,verbose,err_on_fail);
  return;
}

int o2scl_vector_region_int_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_x, void *ptr_y, double intl, void *ptr_locs, int boundaries, int verbose, bool err_on_fail) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  std::vector<double> *locs=(std::vector<double> *)ptr_locs;
  int ret=vector_region_int<std::vector<double>,std::vector<double>>(n,*x,*y,intl,*locs,boundaries,verbose,err_on_fail);
  return ret;
}

int o2scl_vector_region_fracint_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_x, void *ptr_y, double intl, void *ptr_locs, int boundaries, int verbose, bool err_on_fail) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  std::vector<double> *locs=(std::vector<double> *)ptr_locs;
  int ret=vector_region_fracint<std::vector<double>,std::vector<double>>(n,*x,*y,intl,*locs,boundaries,verbose,err_on_fail);
  return ret;
}

int o2scl_vector_bound_fracint_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_x, void *ptr_y, double frac, void *ptr_low, void *ptr_high, int boundaries, int verbose, bool err_on_fail) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  double *low=(double *)ptr_low;
  double *high=(double *)ptr_high;
  int ret=vector_bound_fracint<std::vector<double>,std::vector<double>>(n,*x,*y,frac,*low,*high,boundaries,verbose,err_on_fail);
  return ret;
}

int o2scl_vector_bound_int_std_vector_double_std_vector_double__wrapper(size_t n, void *ptr_x, void *ptr_y, double frac, void *ptr_low, void *ptr_high, int boundaries, int verbose, bool err_on_fail) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  double *low=(double *)ptr_low;
  double *high=(double *)ptr_high;
  int ret=vector_bound_int<std::vector<double>,std::vector<double>>(n,*x,*y,frac,*low,*high,boundaries,verbose,err_on_fail);
  return ret;
}

void o2scl_rebin_xy_std_vector_double_std_vector_double_std_vector_double_std_vector_double__wrapper(void *ptr_x, void *ptr_y, void *ptr_x_out, void *ptr_y_out, size_t n_pts, size_t interp_type) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  std::vector<double> *x_out=(std::vector<double> *)ptr_x_out;
  std::vector<double> *y_out=(std::vector<double> *)ptr_y_out;
  rebin_xy<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>>(*x,*y,*x_out,*y_out,n_pts,interp_type);
  return;
}

double o2scl_linear_or_log_chi2_std_vector_double_std_vector_double__wrapper(void *ptr_x, void *ptr_y) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  double ret=linear_or_log_chi2<std::vector<double>,std::vector<double>>(*x,*y);
  return ret;
}

void o2scl_linear_or_log_std_vector_double_std_vector_double__wrapper(void *ptr_x, void *ptr_y, void *ptr_log_x, void *ptr_log_y) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  std::vector<double> *y=(std::vector<double> *)ptr_y;
  bool *log_x=(bool *)ptr_log_x;
  bool *log_y=(bool *)ptr_log_y;
  linear_or_log<std::vector<double>,std::vector<double>>(*x,*y,*log_x,*log_y);
  return;
}

void o2scl_vector_refine_std_vector_double_std_vector_double_double__wrapper(size_t n, void *ptr_index, void *ptr_data, size_t factor, size_t interp_type) {
  std::vector<double> *index=(std::vector<double> *)ptr_index;
  std::vector<double> *data=(std::vector<double> *)ptr_data;
  vector_refine<std::vector<double>,std::vector<double>,double>(n,*index,*data,factor,interp_type);
  return;
}

void o2scl_linear_or_log_std_vector_double__wrapper(void *ptr_x, void *ptr_log_x) {
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  bool *log_x=(bool *)ptr_log_x;
  linear_or_log<std::vector<double>>(*x,*log_x);
  return;
}

