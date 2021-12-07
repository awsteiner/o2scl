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

#include <o2scl/other_python.h>

using namespace std;
using namespace o2scl;

void *o2scl_create_slack_messenger() {
  slack_messenger *ptr=new slack_messenger;
  return ptr;
}

void o2scl_free_slack_messenger(void *vptr) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  delete ptr;
  return;
}

int o2scl_slack_messenger_get_verbose(void *vptr) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  return ptr->verbose;
}

void o2scl_slack_messenger_set_verbose(void *vptr, int v) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  ptr->verbose=v;
  return;
}

void *o2scl_slack_messenger_get_url(void *vptr) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->url;
  return sptr;
}

void o2scl_slack_messenger_set_url(void *vptr, void *p_v) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  std::string *p_t=(std::string *)p_v;
  ptr->url=*(p_t);
  return;
}

void *o2scl_slack_messenger_get_channel(void *vptr) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->channel;
  return sptr;
}

void o2scl_slack_messenger_set_channel(void *vptr, void *p_v) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  std::string *p_t=(std::string *)p_v;
  ptr->channel=*(p_t);
  return;
}

void *o2scl_slack_messenger_get_icon(void *vptr) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->icon;
  return sptr;
}

void o2scl_slack_messenger_set_icon(void *vptr, void *p_v) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  std::string *p_t=(std::string *)p_v;
  ptr->icon=*(p_t);
  return;
}

void *o2scl_slack_messenger_get_username(void *vptr) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->username;
  return sptr;
}

void o2scl_slack_messenger_set_username(void *vptr, void *p_v) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  std::string *p_t=(std::string *)p_v;
  ptr->username=*(p_t);
  return;
}

double o2scl_slack_messenger_get_min_time_between(void *vptr) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  return ptr->min_time_between;
}

void o2scl_slack_messenger_set_min_time_between(void *vptr, double v) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  ptr->min_time_between=v;
  return;
}

bool o2scl_slack_messenger_set_url_from_env(void *vptr, char *env_var) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  bool ret=ptr->set_url_from_env(env_var);
  return ret;
}

bool o2scl_slack_messenger_set_channel_from_env(void *vptr, char *env_var) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  bool ret=ptr->set_channel_from_env(env_var);
  return ret;
}

bool o2scl_slack_messenger_set_username_from_env(void *vptr, char *env_var) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  bool ret=ptr->set_username_from_env(env_var);
  return ret;
}

int o2scl_slack_messenger_send(void *vptr, char *message, bool err_on_fail) {
  slack_messenger *ptr=(slack_messenger *)vptr;
  int ret=ptr->send(message,err_on_fail);
  return ret;
}

void *o2scl_slack_messenger_init(char *p_channel, char *p_username, char *p_url, bool p_mpi_time) {
  slack_messenger *ptr=new slack_messenger(p_channel,p_username,p_url,p_mpi_time);
  return ptr;
}

void *o2scl_create_quadratic_real_coeff_gsl() {
  quadratic_real_coeff_gsl *ptr=new quadratic_real_coeff_gsl;
  return ptr;
}

void o2scl_free_quadratic_real_coeff_gsl(void *vptr) {
  quadratic_real_coeff_gsl *ptr=(quadratic_real_coeff_gsl *)vptr;
  delete ptr;
  return;
}

int o2scl_quadratic_real_coeff_gsl_solve_r(void *vptr, double a2, double b2, double c2, double *r1, double *r2) {
  quadratic_real_coeff_gsl *ptr=(quadratic_real_coeff_gsl *)vptr;
  int ret=ptr->solve_r(a2,b2,c2,*r1,*r2);
  return ret;
}

int o2scl_quadratic_real_coeff_gsl_solve_rc(void *vptr, double a2, double b2, double c2, void *ptr_r1, void *ptr_r2) {
  quadratic_real_coeff_gsl *ptr=(quadratic_real_coeff_gsl *)vptr;
  std::complex<double> *r1=(std::complex<double> *)ptr_r1;
  std::complex<double> *r2=(std::complex<double> *)ptr_r2;
  int ret=ptr->solve_rc(a2,b2,c2,*r1,*r2);
  return ret;
}

void *o2scl_create_quadratic_real_coeff_gsl2_() {
  quadratic_real_coeff_gsl2<> *ptr=new quadratic_real_coeff_gsl2<>;
  return ptr;
}

void o2scl_free_quadratic_real_coeff_gsl2_(void *vptr) {
  quadratic_real_coeff_gsl2<> *ptr=(quadratic_real_coeff_gsl2<> *)vptr;
  delete ptr;
  return;
}

int o2scl_quadratic_real_coeff_gsl2__solve_r(void *vptr, double a2, double b2, double c2, double *r1, double *r2) {
  quadratic_real_coeff_gsl2<> *ptr=(quadratic_real_coeff_gsl2<> *)vptr;
  int ret=ptr->solve_r(a2,b2,c2,*r1,*r2);
  return ret;
}

int o2scl_quadratic_real_coeff_gsl2__solve_rc(void *vptr, double a2, double b2, double c2, void *ptr_r1, void *ptr_r2) {
  quadratic_real_coeff_gsl2<> *ptr=(quadratic_real_coeff_gsl2<> *)vptr;
  std::complex<double> *r1=(std::complex<double> *)ptr_r1;
  std::complex<double> *r2=(std::complex<double> *)ptr_r2;
  int ret=ptr->solve_rc(a2,b2,c2,*r1,*r2);
  return ret;
}

void *o2scl_create_cubic_real_coeff_cern_() {
  cubic_real_coeff_cern<> *ptr=new cubic_real_coeff_cern<>;
  return ptr;
}

void o2scl_free_cubic_real_coeff_cern_(void *vptr) {
  cubic_real_coeff_cern<> *ptr=(cubic_real_coeff_cern<> *)vptr;
  delete ptr;
  return;
}

int o2scl_cubic_real_coeff_cern__solve_r(void *vptr, double a3, double b3, double c3, double d3, double *r1, double *r2, double *r3) {
  cubic_real_coeff_cern<> *ptr=(cubic_real_coeff_cern<> *)vptr;
  int ret=ptr->solve_r(a3,b3,c3,d3,*r1,*r2,*r3);
  return ret;
}

int o2scl_cubic_real_coeff_cern__solve_rc(void *vptr, double a3, double b3, double c3, double d3, double *r1, void *ptr_r2, void *ptr_r3) {
  cubic_real_coeff_cern<> *ptr=(cubic_real_coeff_cern<> *)vptr;
  std::complex<double> *r2=(std::complex<double> *)ptr_r2;
  std::complex<double> *r3=(std::complex<double> *)ptr_r3;
  int ret=ptr->solve_rc(a3,b3,c3,d3,*r1,*r2,*r3);
  return ret;
}

void *o2scl_create_cubic_real_coeff_gsl() {
  cubic_real_coeff_gsl *ptr=new cubic_real_coeff_gsl;
  return ptr;
}

void o2scl_free_cubic_real_coeff_gsl(void *vptr) {
  cubic_real_coeff_gsl *ptr=(cubic_real_coeff_gsl *)vptr;
  delete ptr;
  return;
}

int o2scl_cubic_real_coeff_gsl_solve_r(void *vptr, double a3, double b3, double c3, double d3, double *r1, double *r2, double *r3) {
  cubic_real_coeff_gsl *ptr=(cubic_real_coeff_gsl *)vptr;
  int ret=ptr->solve_r(a3,b3,c3,d3,*r1,*r2,*r3);
  return ret;
}

int o2scl_cubic_real_coeff_gsl_solve_rc(void *vptr, double a3, double b3, double c3, double d3, double *r1, void *ptr_r2, void *ptr_r3) {
  cubic_real_coeff_gsl *ptr=(cubic_real_coeff_gsl *)vptr;
  std::complex<double> *r2=(std::complex<double> *)ptr_r2;
  std::complex<double> *r3=(std::complex<double> *)ptr_r3;
  int ret=ptr->solve_rc(a3,b3,c3,d3,*r1,*r2,*r3);
  return ret;
}

void *o2scl_create_quartic_real_coeff_cern_() {
  quartic_real_coeff_cern<> *ptr=new quartic_real_coeff_cern<>;
  return ptr;
}

void o2scl_free_quartic_real_coeff_cern_(void *vptr) {
  quartic_real_coeff_cern<> *ptr=(quartic_real_coeff_cern<> *)vptr;
  delete ptr;
  return;
}

int o2scl_quartic_real_coeff_cern__solve_r(void *vptr, double a4, double b4, double c4, double d4, double e4, double *r1, double *r2, double *r3, double *r4) {
  quartic_real_coeff_cern<> *ptr=(quartic_real_coeff_cern<> *)vptr;
  int ret=ptr->solve_r(a4,b4,c4,d4,e4,*r1,*r2,*r3,*r4);
  return ret;
}

int o2scl_quartic_real_coeff_cern__solve_rc(void *vptr, double a4, double b4, double c4, double d4, double e4, void *ptr_r1, void *ptr_r2, void *ptr_r3, void *ptr_r4) {
  quartic_real_coeff_cern<> *ptr=(quartic_real_coeff_cern<> *)vptr;
  std::complex<double> *r1=(std::complex<double> *)ptr_r1;
  std::complex<double> *r2=(std::complex<double> *)ptr_r2;
  std::complex<double> *r3=(std::complex<double> *)ptr_r3;
  std::complex<double> *r4=(std::complex<double> *)ptr_r4;
  int ret=ptr->solve_rc(a4,b4,c4,d4,e4,*r1,*r2,*r3,*r4);
  return ret;
}

void *o2scl_create_fermi_dirac_integ_gsl() {
  fermi_dirac_integ_gsl *ptr=new fermi_dirac_integ_gsl;
  return ptr;
}

void o2scl_free_fermi_dirac_integ_gsl(void *vptr) {
  fermi_dirac_integ_gsl *ptr=(fermi_dirac_integ_gsl *)vptr;
  delete ptr;
  return;
}

double o2scl_fermi_dirac_integ_gsl_calc_m1o2(void *vptr, double x) {
  fermi_dirac_integ_gsl *ptr=(fermi_dirac_integ_gsl *)vptr;
  double ret=ptr->calc_m1o2(x);
  return ret;
}

double o2scl_fermi_dirac_integ_gsl_calc_1o2(void *vptr, double x) {
  fermi_dirac_integ_gsl *ptr=(fermi_dirac_integ_gsl *)vptr;
  double ret=ptr->calc_1o2(x);
  return ret;
}

double o2scl_fermi_dirac_integ_gsl_calc_3o2(void *vptr, double x) {
  fermi_dirac_integ_gsl *ptr=(fermi_dirac_integ_gsl *)vptr;
  double ret=ptr->calc_3o2(x);
  return ret;
}

double o2scl_fermi_dirac_integ_gsl_calc_2(void *vptr, double x) {
  fermi_dirac_integ_gsl *ptr=(fermi_dirac_integ_gsl *)vptr;
  double ret=ptr->calc_2(x);
  return ret;
}

double o2scl_fermi_dirac_integ_gsl_calc_3(void *vptr, double x) {
  fermi_dirac_integ_gsl *ptr=(fermi_dirac_integ_gsl *)vptr;
  double ret=ptr->calc_3(x);
  return ret;
}

void *o2scl_create_bessel_K_exp_integ_gsl() {
  bessel_K_exp_integ_gsl *ptr=new bessel_K_exp_integ_gsl;
  return ptr;
}

void o2scl_free_bessel_K_exp_integ_gsl(void *vptr) {
  bessel_K_exp_integ_gsl *ptr=(bessel_K_exp_integ_gsl *)vptr;
  delete ptr;
  return;
}

double o2scl_bessel_K_exp_integ_gsl_K1exp(void *vptr, double x) {
  bessel_K_exp_integ_gsl *ptr=(bessel_K_exp_integ_gsl *)vptr;
  double ret=ptr->K1exp(x);
  return ret;
}

double o2scl_bessel_K_exp_integ_gsl_K2exp(void *vptr, double x) {
  bessel_K_exp_integ_gsl *ptr=(bessel_K_exp_integ_gsl *)vptr;
  double ret=ptr->K2exp(x);
  return ret;
}

double o2scl_bessel_K_exp_integ_gsl_K3exp(void *vptr, double x) {
  bessel_K_exp_integ_gsl *ptr=(bessel_K_exp_integ_gsl *)vptr;
  double ret=ptr->K3exp(x);
  return ret;
}

void *o2scl_create_hist() {
  hist *ptr=new hist;
  return ptr;
}

void o2scl_free_hist(void *vptr) {
  hist *ptr=(hist *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_hist(void *vsrc, void *vdest) {
  hist *src=(hist *)vsrc;
  hist *dest=(hist *)vdest;
  *dest=*src;
}

bool o2scl_hist_get_extend_rhs(void *vptr) {
  hist *ptr=(hist *)vptr;
  return ptr->extend_rhs;
}

void o2scl_hist_set_extend_rhs(void *vptr, bool v) {
  hist *ptr=(hist *)vptr;
  ptr->extend_rhs=v;
  return;
}

bool o2scl_hist_get_extend_lhs(void *vptr) {
  hist *ptr=(hist *)vptr;
  return ptr->extend_lhs;
}

void o2scl_hist_set_extend_lhs(void *vptr, bool v) {
  hist *ptr=(hist *)vptr;
  ptr->extend_lhs=v;
  return;
}

size_t o2scl_hist_size(void *vptr) {
  hist *ptr=(hist *)vptr;
  size_t ret=ptr->size();
  return ret;
}

void o2scl_hist_set_bin_edges_grid(void *vptr, void *ptr_g) {
  hist *ptr=(hist *)vptr;
  uniform_grid<double> *g=(uniform_grid<double> *)ptr_g;
  ptr->set_bin_edges(*g);
  return;
}

void o2scl_hist_set_bin_edges_vec(void *vptr, size_t n, void *ptr_v) {
  hist *ptr=(hist *)vptr;
  vector<double> *v=(vector<double> *)ptr_v;
  ptr->set_bin_edges(n,*v);
  return;
}

void o2scl_hist_update(void *vptr, double x, double val) {
  hist *ptr=(hist *)vptr;
  ptr->update(x,val);
  return;
}

void o2scl_hist_update_i(void *vptr, size_t i, double val) {
  hist *ptr=(hist *)vptr;
  ptr->update_i(i,val);
  return;
}

double o2scl_hist_get_wgt_i(void *vptr, size_t i) {
  hist *ptr=(hist *)vptr;
  double ret=ptr->get_wgt_i(i);
  return ret;
}

double o2scl_hist_get_wgt(void *vptr, double x) {
  hist *ptr=(hist *)vptr;
  double ret=ptr->get_wgt(x);
  return ret;
}

void o2scl_hist_set_wgt_i(void *vptr, size_t i, double x) {
  hist *ptr=(hist *)vptr;
  ptr->set_wgt_i(i,x);
  return;
}

void o2scl_hist_set_wgt(void *vptr, double x, double val) {
  hist *ptr=(hist *)vptr;
  ptr->set_wgt(x,val);
  return;
}

const double o2scl_hist_getitem(void *vptr, size_t i) {
  hist *ptr=(hist *)vptr;
  double ret=ptr->operator[](i);
  return ret;
}

size_t o2scl_hist_get_bin_index(void *vptr, double x) {
  hist *ptr=(hist *)vptr;
  size_t ret=ptr->get_bin_index(x);
  return ret;
}

int o2scl_hist_function(void *vptr, char *func) {
  hist *ptr=(hist *)vptr;
  int ret=ptr->function(func);
  return ret;
}

void o2scl_hist_clear(void *vptr) {
  hist *ptr=(hist *)vptr;
  ptr->clear();
  return;
}

void *o2scl_create_fract() {
  fract *ptr=new fract;
  return ptr;
}

void o2scl_free_fract(void *vptr) {
  fract *ptr=(fract *)vptr;
  delete ptr;
  return;
}

void o2scl_fract_nrf_z4m1(void *vptr, void *ptr_gx, void *ptr_gy, size_t kmax, double rmax, void *ptr_t3d, void *ptr_roots_x, void *ptr_roots_y, void *ptr_min, void *ptr_max) {
  fract *ptr=(fract *)vptr;
  uniform_grid<> *gx=(uniform_grid<> *)ptr_gx;
  uniform_grid<> *gy=(uniform_grid<> *)ptr_gy;
  o2scl::table3d *t3d=(o2scl::table3d *)ptr_t3d;
  std::vector<double> *roots_x=(std::vector<double> *)ptr_roots_x;
  std::vector<double> *roots_y=(std::vector<double> *)ptr_roots_y;
  std::vector<double> *min=(std::vector<double> *)ptr_min;
  std::vector<double> *max=(std::vector<double> *)ptr_max;
  ptr->nrf_z4m1(*gx,*gy,kmax,rmax,*t3d,*roots_x,*roots_y,*min,*max);
  return;
}

int o2scl_fract_itf_mandel(void *vptr, void *ptr_gx, void *ptr_gy, size_t kmax, double rmax, void *ptr_t3d, size_t *min, size_t *max) {
  fract *ptr=(fract *)vptr;
  uniform_grid<> *gx=(uniform_grid<> *)ptr_gx;
  uniform_grid<> *gy=(uniform_grid<> *)ptr_gy;
  o2scl::table3d *t3d=(o2scl::table3d *)ptr_t3d;
  int ret=ptr->itf_mandel(*gx,*gy,kmax,rmax,*t3d,*min,*max);
  return ret;
}

