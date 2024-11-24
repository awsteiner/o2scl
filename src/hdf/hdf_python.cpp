/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2020-2024, Andrew W. Steiner

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

#include <o2scl/hdf_python.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;
using namespace o2scl_acol;

void *o2scl_hdf_create_hdf_file() {
  hdf_file *ptr=new hdf_file;
  return ptr;
}

void o2scl_hdf_free_hdf_file(void *vptr) {
  hdf_file *ptr=(hdf_file *)vptr;
  delete ptr;
  return;
}

int o2scl_hdf_hdf_file_get_compr_type(void *vptr) {
  hdf_file *ptr=(hdf_file *)vptr;
  return ptr->compr_type;
}

void o2scl_hdf_hdf_file_set_compr_type(void *vptr, int v) {
  hdf_file *ptr=(hdf_file *)vptr;
  ptr->compr_type=v;
  return;
}

size_t o2scl_hdf_hdf_file_get_min_compr_size(void *vptr) {
  hdf_file *ptr=(hdf_file *)vptr;
  return ptr->min_compr_size;
}

void o2scl_hdf_hdf_file_set_min_compr_size(void *vptr, size_t v) {
  hdf_file *ptr=(hdf_file *)vptr;
  ptr->min_compr_size=v;
  return;
}

bool o2scl_hdf_hdf_file_has_write_access(void *vptr) {
  hdf_file *ptr=(hdf_file *)vptr;
  bool ret=ptr->has_write_access();
  return ret;
}

void o2scl_hdf_hdf_file_open(void *vptr, void *ptr_fname, bool write_access, bool err_on_fail) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *fname=(std::string *)ptr_fname;
  ptr->open(*fname,write_access,err_on_fail);
  return;
}

void o2scl_hdf_hdf_file_open_or_create(void *vptr, void *ptr_fname) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *fname=(std::string *)ptr_fname;
  ptr->open_or_create(*fname);
  return;
}

void o2scl_hdf_hdf_file_close(void *vptr) {
  hdf_file *ptr=(hdf_file *)vptr;
  ptr->close();
  return;
}

int o2scl_hdf_hdf_file_getc(void *vptr, void *ptr_name, char *c) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  int ret=ptr->getc(*name,*c);
  return ret;
}

int o2scl_hdf_hdf_file_getd(void *vptr, void *ptr_name, double *d) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  int ret=ptr->getd(*name,*d);
  return ret;
}

int o2scl_hdf_hdf_file_geti(void *vptr, void *ptr_name, int *i) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  int ret=ptr->geti(*name,*i);
  return ret;
}

int o2scl_hdf_hdf_file_get_szt(void *vptr, void *ptr_name, size_t *u) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  int ret=ptr->get_szt(*name,*u);
  return ret;
}

int o2scl_hdf_hdf_file_gets(void *vptr, void *ptr_name, void *ptr_s) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::string *s=(std::string *)ptr_s;
  int ret=ptr->gets(*name,*s);
  return ret;
}

int o2scl_hdf_hdf_file_gets_var(void *vptr, void *ptr_name, void *ptr_s) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::string *s=(std::string *)ptr_s;
  int ret=ptr->gets_var(*name,*s);
  return ret;
}

int o2scl_hdf_hdf_file_gets_fixed(void *vptr, void *ptr_name, void *ptr_s) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::string *s=(std::string *)ptr_s;
  int ret=ptr->gets_fixed(*name,*s);
  return ret;
}

int o2scl_hdf_hdf_file_gets_def_fixed(void *vptr, void *ptr_name, void *ptr_deft, void *ptr_s) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::string *deft=(std::string *)ptr_deft;
  std::string *s=(std::string *)ptr_s;
  int ret=ptr->gets_def_fixed(*name,*deft,*s);
  return ret;
}

void o2scl_hdf_hdf_file_setc(void *vptr, void *ptr_name, char c) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  ptr->setc(*name,c);
  return;
}

void o2scl_hdf_hdf_file_setd(void *vptr, void *ptr_name, double d) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  ptr->setd(*name,d);
  return;
}

void o2scl_hdf_hdf_file_seti(void *vptr, void *ptr_name, int i) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  ptr->seti(*name,i);
  return;
}

void o2scl_hdf_hdf_file_set_szt(void *vptr, void *ptr_name, size_t u) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  ptr->set_szt(*name,u);
  return;
}

void o2scl_hdf_hdf_file_sets(void *vptr, void *ptr_name, void *ptr_s) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::string *s=(std::string *)ptr_s;
  ptr->sets(*name,*s);
  return;
}

void o2scl_hdf_hdf_file_sets_fixed(void *vptr, void *ptr_name, void *ptr_s) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::string *s=(std::string *)ptr_s;
  ptr->sets_fixed(*name,*s);
  return;
}

int o2scl_hdf_hdf_file_getd_vec(void *vptr, void *ptr_name, void *ptr_v) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::vector<double> *v=(std::vector<double> *)ptr_v;
  int ret=ptr->getd_vec(*name,*v);
  return ret;
}

int o2scl_hdf_hdf_file_geti_vec(void *vptr, void *ptr_name, void *ptr_v) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::vector<int> *v=(std::vector<int> *)ptr_v;
  int ret=ptr->geti_vec(*name,*v);
  return ret;
}

int o2scl_hdf_hdf_file_get_szt_vec(void *vptr, void *ptr_name, void *ptr_v) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::vector<size_t> *v=(std::vector<size_t> *)ptr_v;
  int ret=ptr->get_szt_vec(*name,*v);
  return ret;
}

int o2scl_hdf_hdf_file_gets_vec_copy(void *vptr, void *ptr_name, void *ptr_s) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::vector<std::string> *s=(std::vector<std::string> *)ptr_s;
  int ret=ptr->gets_vec_copy(*name,*s);
  return ret;
}

int o2scl_hdf_hdf_file_setd_vec(void *vptr, void *ptr_name, void *ptr_v) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::vector<double> *v=(std::vector<double> *)ptr_v;
  int ret=ptr->setd_vec(*name,*v);
  return ret;
}

int o2scl_hdf_hdf_file_seti_vec(void *vptr, void *ptr_name, void *ptr_v) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::vector<int> *v=(std::vector<int> *)ptr_v;
  int ret=ptr->seti_vec(*name,*v);
  return ret;
}

int o2scl_hdf_hdf_file_set_szt_vec(void *vptr, void *ptr_name, void *ptr_v) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::vector<size_t> *v=(std::vector<size_t> *)ptr_v;
  int ret=ptr->set_szt_vec(*name,*v);
  return ret;
}

int o2scl_hdf_hdf_file_sets_vec_copy(void *vptr, void *ptr_name, void *ptr_s) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::vector<std::string> *s=(std::vector<std::string> *)ptr_s;
  int ret=ptr->sets_vec_copy(*name,*s);
  return ret;
}

int o2scl_hdf_hdf_file_getd_mat_copy(void *vptr, void *ptr_name, void *ptr_m) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  boost::numeric::ublas::matrix<double> *m=(boost::numeric::ublas::matrix<double> *)ptr_m;
  int ret=ptr->getd_mat_copy(*name,*m);
  return ret;
}

int o2scl_hdf_hdf_file_geti_mat_copy(void *vptr, void *ptr_name, void *ptr_m) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  boost::numeric::ublas::matrix<int> *m=(boost::numeric::ublas::matrix<int> *)ptr_m;
  int ret=ptr->geti_mat_copy(*name,*m);
  return ret;
}

int o2scl_hdf_hdf_file_setd_mat_copy(void *vptr, void *ptr_name, void *ptr_m) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  boost::numeric::ublas::matrix<double> *m=(boost::numeric::ublas::matrix<double> *)ptr_m;
  int ret=ptr->setd_mat_copy(*name,*m);
  return ret;
}

int o2scl_hdf_hdf_file_seti_mat_copy(void *vptr, void *ptr_name, void *ptr_m) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  boost::numeric::ublas::matrix<int> *m=(boost::numeric::ublas::matrix<int> *)ptr_m;
  int ret=ptr->seti_mat_copy(*name,*m);
  return ret;
}

int o2scl_hdf_hdf_file_getd_ten(void *vptr, void *ptr_name, void *ptr_t) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  tensor<> *t=(tensor<> *)ptr_t;
  int ret=ptr->getd_ten(*name,*t);
  return ret;
}

int o2scl_hdf_hdf_file_geti_ten(void *vptr, void *ptr_name, void *ptr_t) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  tensor<int> *t=(tensor<int> *)ptr_t;
  int ret=ptr->geti_ten(*name,*t);
  return ret;
}

int o2scl_hdf_hdf_file_get_szt_ten(void *vptr, void *ptr_name, void *ptr_t) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  tensor<size_t> *t=(tensor<size_t> *)ptr_t;
  int ret=ptr->get_szt_ten(*name,*t);
  return ret;
}

int o2scl_hdf_hdf_file_setd_ten(void *vptr, void *ptr_name, void *ptr_t) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  tensor<> *t=(tensor<> *)ptr_t;
  int ret=ptr->setd_ten(*name,*t);
  return ret;
}

int o2scl_hdf_hdf_file_seti_ten(void *vptr, void *ptr_name, void *ptr_t) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  tensor<int> *t=(tensor<int> *)ptr_t;
  int ret=ptr->seti_ten(*name,*t);
  return ret;
}

int o2scl_hdf_hdf_file_set_szt_ten(void *vptr, void *ptr_name, void *ptr_t) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  tensor<size_t> *t=(tensor<size_t> *)ptr_t;
  int ret=ptr->set_szt_ten(*name,*t);
  return ret;
}

int o2scl_hdf_hdf_file_getc_def(void *vptr, void *ptr_name, char deft, char *c) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  int ret=ptr->getc_def(*name,deft,*c);
  return ret;
}

int o2scl_hdf_hdf_file_getd_def(void *vptr, void *ptr_name, double deft, double *d) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  int ret=ptr->getd_def(*name,deft,*d);
  return ret;
}

int o2scl_hdf_hdf_file_geti_def(void *vptr, void *ptr_name, int deft, int *i) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  int ret=ptr->geti_def(*name,deft,*i);
  return ret;
}

int o2scl_hdf_hdf_file_get_szt_def(void *vptr, void *ptr_name, size_t deft, size_t *u) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  int ret=ptr->get_szt_def(*name,deft,*u);
  return ret;
}

int o2scl_hdf_hdf_file_gets_def(void *vptr, void *ptr_name, void *ptr_deft, void *ptr_s) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::string *deft=(std::string *)ptr_deft;
  std::string *s=(std::string *)ptr_s;
  int ret=ptr->gets_def(*name,*deft,*s);
  return ret;
}

int o2scl_hdf_hdf_file_find_object_by_type(void *vptr, void *ptr_type, void *ptr_name, bool use_regex, int verbose) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *type=(std::string *)ptr_type;
  std::string *name=(std::string *)ptr_name;
  int ret=ptr->find_object_by_type(*type,*name,use_regex,verbose);
  return ret;
}

int o2scl_hdf_hdf_file_find_object_by_name(void *vptr, void *ptr_name, void *ptr_type, bool use_regex, int verbose) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *name=(std::string *)ptr_name;
  std::string *type=(std::string *)ptr_type;
  int ret=ptr->find_object_by_name(*name,*type,use_regex,verbose);
  return ret;
}

int o2scl_hdf_hdf_file_find_object_by_pattern(void *vptr, void *ptr_pattern, void *ptr_type, bool use_regex, int verbose) {
  hdf_file *ptr=(hdf_file *)vptr;
  std::string *pattern=(std::string *)ptr_pattern;
  std::string *type=(std::string *)ptr_type;
  int ret=ptr->find_object_by_pattern(*pattern,*type,use_regex,verbose);
  return ret;
}

void o2scl_hdf_hdf_file_file_list(void *vptr, int verbose) {
  hdf_file *ptr=(hdf_file *)vptr;
  ptr->file_list(verbose);
  return;
}

void o2scl_hdf_hdf_file_copy(void *vptr, int verbose, void *ptr_hf2) {
  hdf_file *ptr=(hdf_file *)vptr;
  hdf_file *hf2=(hdf_file *)ptr_hf2;
  ptr->copy(verbose,*hf2);
  return;
}

void *o2scl_hdf_create_acol_manager() {
  acol_manager *ptr=new acol_manager;
  return ptr;
}

void o2scl_hdf_free_acol_manager(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  delete ptr;
  return;
}

void *o2scl_hdf_acol_manager_get_env_var_name(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->env_var_name;
  return sptr;
}

void o2scl_hdf_acol_manager_set_env_var_name(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->env_var_name=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_cl(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)((ptr->cl));
}

void o2scl_hdf_acol_manager_set_cl(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  cli *p_tsptr=(cli *)p_v;
  ptr->cl=p_tsptr;
  return;
}

int o2scl_hdf_acol_manager_get_verbose(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return ptr->verbose;
}

void o2scl_hdf_acol_manager_set_verbose(void *vptr, int v) {
  acol_manager *ptr=(acol_manager *)vptr;
  ptr->verbose=v;
  return;
}

void *o2scl_hdf_acol_manager_get_def_args(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->def_args;
  return sptr;
}

void o2scl_hdf_acol_manager_set_def_args(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->def_args=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_type(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->type;
  return sptr;
}

void o2scl_hdf_acol_manager_set_type(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->type=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_table_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->table_obj));
}

void o2scl_hdf_acol_manager_set_table_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  table_units<> *p_tsot=(table_units<> *)p_v;
  ptr->table_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_table3d_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->table3d_obj));
}

void o2scl_hdf_acol_manager_set_table3d_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  table3d *p_tsot=(table3d *)p_v;
  ptr->table3d_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_hist_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->hist_obj));
}

void o2scl_hdf_acol_manager_set_hist_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  hist *p_tsot=(hist *)p_v;
  ptr->hist_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_hist_2d_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->hist_2d_obj));
}

void o2scl_hdf_acol_manager_set_hist_2d_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  hist_2d *p_tsot=(hist_2d *)p_v;
  ptr->hist_2d_obj=*(p_tsot);
  return;
}

int o2scl_hdf_acol_manager_get_int_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return ptr->int_obj;
}

void o2scl_hdf_acol_manager_set_int_obj(void *vptr, int v) {
  acol_manager *ptr=(acol_manager *)vptr;
  ptr->int_obj=v;
  return;
}

char o2scl_hdf_acol_manager_get_char_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return ptr->char_obj;
}

void o2scl_hdf_acol_manager_set_char_obj(void *vptr, char v) {
  acol_manager *ptr=(acol_manager *)vptr;
  ptr->char_obj=v;
  return;
}

double o2scl_hdf_acol_manager_get_double_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return ptr->double_obj;
}

void o2scl_hdf_acol_manager_set_double_obj(void *vptr, double v) {
  acol_manager *ptr=(acol_manager *)vptr;
  ptr->double_obj=v;
  return;
}

size_t o2scl_hdf_acol_manager_get_size_t_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return ptr->size_t_obj;
}

void o2scl_hdf_acol_manager_set_size_t_obj(void *vptr, size_t v) {
  acol_manager *ptr=(acol_manager *)vptr;
  ptr->size_t_obj=v;
  return;
}

void *o2scl_hdf_acol_manager_get_string_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->string_obj;
  return sptr;
}

void o2scl_hdf_acol_manager_set_string_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->string_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_cont_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->cont_obj));
}

void o2scl_hdf_acol_manager_set_cont_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::vector<contour_line> *p_tsot=(std::vector<contour_line> *)p_v;
  ptr->cont_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_ug_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->ug_obj));
}

void o2scl_hdf_acol_manager_set_ug_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  uniform_grid<double> *p_tsot=(uniform_grid<double> *)p_v;
  ptr->ug_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_intv_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->intv_obj));
}

void o2scl_hdf_acol_manager_set_intv_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::vector<int> *p_tsot=(std::vector<int> *)p_v;
  ptr->intv_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_doublev_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->doublev_obj));
}

void o2scl_hdf_acol_manager_set_doublev_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::vector<double> *p_tsot=(std::vector<double> *)p_v;
  ptr->doublev_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_size_tv_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->size_tv_obj));
}

void o2scl_hdf_acol_manager_set_size_tv_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::vector<size_t> *p_tsot=(std::vector<size_t> *)p_v;
  ptr->size_tv_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_stringv_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->stringv_obj));
}

void o2scl_hdf_acol_manager_set_stringv_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::vector<std::string> *p_tsot=(std::vector<std::string> *)p_v;
  ptr->stringv_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_vvstring_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->vvstring_obj));
}

void o2scl_hdf_acol_manager_set_vvstring_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::vector<std::vector<std::string>> *p_tsot=(std::vector<std::vector<std::string>> *)p_v;
  ptr->vvstring_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_vvdouble_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->vvdouble_obj));
}

void o2scl_hdf_acol_manager_set_vvdouble_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::vector<std::vector<double>> *p_tsot=(std::vector<std::vector<double>> *)p_v;
  ptr->vvdouble_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_tensor_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->tensor_obj));
}

void o2scl_hdf_acol_manager_set_tensor_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  tensor<> *p_tsot=(tensor<> *)p_v;
  ptr->tensor_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_tensor_int_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->tensor_int_obj));
}

void o2scl_hdf_acol_manager_set_tensor_int_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  tensor<int> *p_tsot=(tensor<int> *)p_v;
  ptr->tensor_int_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_tensor_size_t_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->tensor_size_t_obj));
}

void o2scl_hdf_acol_manager_set_tensor_size_t_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  tensor<size_t> *p_tsot=(tensor<size_t> *)p_v;
  ptr->tensor_size_t_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_tensor_grid_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->tensor_grid_obj));
}

void o2scl_hdf_acol_manager_set_tensor_grid_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  tensor_grid<> *p_tsot=(tensor_grid<> *)p_v;
  ptr->tensor_grid_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_pdma_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->pdma_obj));
}

void o2scl_hdf_acol_manager_set_pdma_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  prob_dens_mdim_amr<> *p_tsot=(prob_dens_mdim_amr<> *)p_v;
  ptr->pdma_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_pdmg_obj(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  return (void *)(&(ptr->pdmg_obj));
}

void o2scl_hdf_acol_manager_set_pdmg_obj(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  prob_dens_mdim_gaussian<> *p_tsot=(prob_dens_mdim_gaussian<> *)p_v;
  ptr->pdmg_obj=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_command_color(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->command_color;
  return sptr;
}

void o2scl_hdf_acol_manager_set_command_color(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->command_color=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_type_color(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->type_color;
  return sptr;
}

void o2scl_hdf_acol_manager_set_type_color(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->type_color=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_param_color(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->param_color;
  return sptr;
}

void o2scl_hdf_acol_manager_set_param_color(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->param_color=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_help_color(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->help_color;
  return sptr;
}

void o2scl_hdf_acol_manager_set_help_color(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->help_color=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_exec_color(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->exec_color;
  return sptr;
}

void o2scl_hdf_acol_manager_set_exec_color(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->exec_color=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_url_color(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->url_color;
  return sptr;
}

void o2scl_hdf_acol_manager_set_url_color(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->url_color=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_default_color(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->default_color;
  return sptr;
}

void o2scl_hdf_acol_manager_set_default_color(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->default_color=*(p_tsot);
  return;
}

void *o2scl_hdf_acol_manager_get_color_spec(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->color_spec;
  return sptr;
}

void o2scl_hdf_acol_manager_set_color_spec(void *vptr, void *p_v) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->color_spec=*(p_tsot);
  return;
}

bool o2scl_hdf_acol_manager_help_found(void *vptr, void *ptr_arg1, void *ptr_arg2) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *arg1=(std::string *)ptr_arg1;
  std::string *arg2=(std::string *)ptr_arg2;
  bool ret=ptr->help_found(*arg1,*arg2);
  return ret;
}

int o2scl_hdf_acol_manager_run_empty(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  int ret=ptr->run_empty();
  return ret;
}

int o2scl_hdf_acol_manager_validate_interp_type(void *vptr) {
  acol_manager *ptr=(acol_manager *)vptr;
  int ret=ptr->validate_interp_type();
  return ret;
}

void o2scl_hdf_acol_manager_parse_vec_string(void *vptr, void *ptr_args) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::vector<std::string> *args=(std::vector<std::string> *)ptr_args;
  ptr->parse_vec_string(*args);
  return;
}

void o2scl_hdf_acol_manager_command_add(void *vptr, void *ptr_new_type) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *new_type=(std::string *)ptr_new_type;
  ptr->command_add(*new_type);
  return;
}

void o2scl_hdf_acol_manager_command_del(void *vptr, void *ptr_ltype) {
  acol_manager *ptr=(acol_manager *)vptr;
  std::string *ltype=(std::string *)ptr_ltype;
  ptr->command_del(*ltype);
  return;
}

void *o2scl_hdf_create_cloud_file() {
  cloud_file *ptr=new cloud_file;
  return ptr;
}

void o2scl_hdf_free_cloud_file(void *vptr) {
  cloud_file *ptr=(cloud_file *)vptr;
  delete ptr;
  return;
}

int o2scl_hdf_cloud_file_get_hash_type(void *vptr) {
  cloud_file *ptr=(cloud_file *)vptr;
  return ptr->hash_type;
}

void o2scl_hdf_cloud_file_set_hash_type(void *vptr, int v) {
  cloud_file *ptr=(cloud_file *)vptr;
  ptr->hash_type=v;
  return;
}

int o2scl_hdf_cloud_file_get_verbose(void *vptr) {
  cloud_file *ptr=(cloud_file *)vptr;
  return ptr->verbose;
}

void o2scl_hdf_cloud_file_set_verbose(void *vptr, int v) {
  cloud_file *ptr=(cloud_file *)vptr;
  ptr->verbose=v;
  return;
}

bool o2scl_hdf_cloud_file_get_throw_on_fail(void *vptr) {
  cloud_file *ptr=(cloud_file *)vptr;
  return ptr->throw_on_fail;
}

void o2scl_hdf_cloud_file_set_throw_on_fail(void *vptr, bool v) {
  cloud_file *ptr=(cloud_file *)vptr;
  ptr->throw_on_fail=v;
  return;
}

bool o2scl_hdf_cloud_file_get_allow_wget(void *vptr) {
  cloud_file *ptr=(cloud_file *)vptr;
  return ptr->allow_wget;
}

void o2scl_hdf_cloud_file_set_allow_wget(void *vptr, bool v) {
  cloud_file *ptr=(cloud_file *)vptr;
  ptr->allow_wget=v;
  return;
}

bool o2scl_hdf_cloud_file_get_allow_curl(void *vptr) {
  cloud_file *ptr=(cloud_file *)vptr;
  return ptr->allow_curl;
}

void o2scl_hdf_cloud_file_set_allow_curl(void *vptr, bool v) {
  cloud_file *ptr=(cloud_file *)vptr;
  ptr->allow_curl=v;
  return;
}

int o2scl_hdf_cloud_file_get_file(void *vptr, void *ptr_file, void *ptr_url, void *ptr_dir) {
  cloud_file *ptr=(cloud_file *)vptr;
  std::string *file=(std::string *)ptr_file;
  std::string *url=(std::string *)ptr_url;
  std::string *dir=(std::string *)ptr_dir;
  int ret=ptr->get_file(*file,*url,*dir);
  return ret;
}

int o2scl_hdf_cloud_file_get_file_hash(void *vptr, void *ptr_file, void *ptr_url, void *ptr_hash, void *ptr_dir) {
  cloud_file *ptr=(cloud_file *)vptr;
  std::string *file=(std::string *)ptr_file;
  std::string *url=(std::string *)ptr_url;
  std::string *hash=(std::string *)ptr_hash;
  std::string *dir=(std::string *)ptr_dir;
  int ret=ptr->get_file_hash(*file,*url,*hash,*dir);
  return ret;
}

int o2scl_hdf_cloud_file_hdf5_open(void *vptr, void *ptr_hf, void *ptr_file, void *ptr_url, void *ptr_dir) {
  cloud_file *ptr=(cloud_file *)vptr;
  hdf_file *hf=(hdf_file *)ptr_hf;
  std::string *file=(std::string *)ptr_file;
  std::string *url=(std::string *)ptr_url;
  std::string *dir=(std::string *)ptr_dir;
  int ret=ptr->hdf5_open(*hf,*file,*url,*dir);
  return ret;
}

int o2scl_hdf_cloud_file_hdf5_open_hash(void *vptr, void *ptr_hf, void *ptr_file, void *ptr_url, void *ptr_hash, void *ptr_dir) {
  cloud_file *ptr=(cloud_file *)vptr;
  hdf_file *hf=(hdf_file *)ptr_hf;
  std::string *file=(std::string *)ptr_file;
  std::string *url=(std::string *)ptr_url;
  std::string *hash=(std::string *)ptr_hash;
  std::string *dir=(std::string *)ptr_dir;
  int ret=ptr->hdf5_open_hash(*hf,*file,*url,*hash,*dir);
  return ret;
}

void o2scl_hdf_hdf_input_table_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table<> *t=(table<> *)ptr_t;
  hdf_input(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_input_n_table_wrapper(void *ptr_hf, void *ptr_t, void *ptr_name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table<> *t=(table<> *)ptr_t;
  std::string *name=(std::string *)ptr_name;
  hdf_input_n(*hf,*t,*name);
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

void o2scl_hdf_hdf_input_n_table_units_wrapper(void *ptr_hf, void *ptr_t, void *ptr_name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table_units<> *t=(table_units<> *)ptr_t;
  std::string *name=(std::string *)ptr_name;
  hdf_input_n(*hf,*t,*name);
  return;
}

void o2scl_hdf_hdf_output_table_units_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table_units<> *t=(table_units<> *)ptr_t;
  hdf_output(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_input_table3d_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table3d *t=(table3d *)ptr_t;
  hdf_input(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_input_n_table3d_wrapper(void *ptr_hf, void *ptr_t, void *ptr_name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table3d *t=(table3d *)ptr_t;
  std::string *name=(std::string *)ptr_name;
  hdf_input_n(*hf,*t,*name);
  return;
}

void o2scl_hdf_hdf_output_table3d_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  table3d *t=(table3d *)ptr_t;
  hdf_output(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_input_uniform_grid_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  uniform_grid<> *t=(uniform_grid<> *)ptr_t;
  hdf_input(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_input_n_uniform_grid_wrapper(void *ptr_hf, void *ptr_t, void *ptr_name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  uniform_grid<> *t=(uniform_grid<> *)ptr_t;
  std::string *name=(std::string *)ptr_name;
  hdf_input_n(*hf,*t,*name);
  return;
}

void o2scl_hdf_hdf_output_uniform_grid_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  uniform_grid<> *t=(uniform_grid<> *)ptr_t;
  hdf_output(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_input_tensor_grid_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  tensor_grid<> *t=(tensor_grid<> *)ptr_t;
  hdf_input(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_input_n_tensor_grid_wrapper(void *ptr_hf, void *ptr_t, void *ptr_name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  tensor_grid<> *t=(tensor_grid<> *)ptr_t;
  std::string *name=(std::string *)ptr_name;
  hdf_input_n(*hf,*t,*name);
  return;
}

void o2scl_hdf_hdf_output_tensor_grid_wrapper(void *ptr_hf, void *ptr_t, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  tensor_grid<> *t=(tensor_grid<> *)ptr_t;
  hdf_output(*hf,*t,name);
  return;
}

void o2scl_hdf_hdf_input_vector_contour_line_wrapper(void *ptr_hf, void *ptr_v, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  std::vector<contour_line> *v=(std::vector<contour_line> *)ptr_v;
  hdf_input(*hf,*v,name);
  return;
}

void o2scl_hdf_hdf_input_n_vector_contour_line_wrapper(void *ptr_hf, void *ptr_v, void *ptr_name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  std::vector<contour_line> *v=(std::vector<contour_line> *)ptr_v;
  std::string *name=(std::string *)ptr_name;
  hdf_input_n(*hf,*v,*name);
  return;
}

void o2scl_hdf_hdf_output_vector_contour_line_wrapper(void *ptr_hf, void *ptr_v, char *name) {
  hdf_file *hf=(hdf_file *)ptr_hf;
  std::vector<contour_line> *v=(std::vector<contour_line> *)ptr_v;
  hdf_output(*hf,*v,name);
  return;
}

int o2scl_hdf_value_spec_wrapper(char *spec, void *ptr_d, int verbose, bool err_on_fail) {
  double *d=(double *)ptr_d;
  int func_ret=value_spec(spec,*d,verbose,err_on_fail);
  return func_ret;
}

int o2scl_hdf_vector_spec_std_vector_double__wrapper(char *spec, void *ptr_v, int verbose, bool err_on_fail) {
  std::vector<double> *v=(std::vector<double> *)ptr_v;
  int func_ret=vector_spec<std::vector<double>>(spec,*v,verbose,err_on_fail);
  return func_ret;
}

int o2scl_hdf_strings_spec_std_vector_std_string__wrapper(char *spec, void *ptr_v, int verbose, bool err_on_fail) {
  std::vector<std::string> *v=(std::vector<std::string> *)ptr_v;
  int func_ret=strings_spec<std::vector<std::string>>(spec,*v,verbose,err_on_fail);
  return func_ret;
}

void *o2scl_hdf_vector_spec_wrapper(char *spec) {
  std::vector<double> *func_ret=new std::vector<double>;
  *func_ret=vector_spec(spec);
  return func_ret;
}

int o2scl_hdf_mult_vector_spec_std_vector_double__wrapper(char *spec, void *ptr_v, bool use_regex, int verbose, bool err_on_fail) {
  std::vector<std::vector<double>> *v=(std::vector<std::vector<double>> *)ptr_v;
  int func_ret=mult_vector_spec<std::vector<double>>(spec,*v,use_regex,verbose,err_on_fail);
  return func_ret;
}

