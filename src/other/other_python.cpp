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
  std::string *p_tsot=(std::string *)p_v;
  ptr->url=*(p_tsot);
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
  std::string *p_tsot=(std::string *)p_v;
  ptr->channel=*(p_tsot);
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
  std::string *p_tsot=(std::string *)p_v;
  ptr->icon=*(p_tsot);
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
  std::string *p_tsot=(std::string *)p_v;
  ptr->username=*(p_tsot);
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

void o2scl_hist_create_rep_vec(void *vptr, void *ptr_v) {
  hist *ptr=(hist *)vptr;
  std::vector<double> *v=(std::vector<double> *)ptr_v;
  ptr->create_rep_vec(*v);
  return;
}

void *o2scl_hist_get_wgts(void *vptr) {
  hist *ptr=(hist *)vptr;
  const boost::numeric::ublas::vector<double> *ret=&ptr->get_wgts();
  return (void *)ret;
}

void *o2scl_hist_get_bins(void *vptr) {
  hist *ptr=(hist *)vptr;
  const boost::numeric::ublas::vector<double> *ret=&ptr->get_bins();
  return (void *)ret;
}

void o2scl_hist_from_table(void *vptr, void *ptr_t, char *colx, size_t n_bins) {
  hist *ptr=(hist *)vptr;
  table<> *t=(table<> *)ptr_t;
  ptr->from_table(*t,colx,n_bins);
  return;
}

void o2scl_hist_from_table_twocol(void *vptr, void *ptr_t, char *colx, char *coly, size_t n_bins) {
  hist *ptr=(hist *)vptr;
  table<> *t=(table<> *)ptr_t;
  ptr->from_table(*t,colx,coly,n_bins);
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

void *o2scl_create_contour_line() {
  contour_line *ptr=new contour_line;
  return ptr;
}

void o2scl_free_contour_line(void *vptr) {
  contour_line *ptr=(contour_line *)vptr;
  delete ptr;
  return;
}

void o2scl_copy_contour_line(void *vsrc, void *vdest) {
  contour_line *src=(contour_line *)vsrc;
  contour_line *dest=(contour_line *)vdest;
  *dest=*src;
}

double o2scl_contour_line_get_level(void *vptr) {
  contour_line *ptr=(contour_line *)vptr;
  return ptr->level;
}

void o2scl_contour_line_set_level(void *vptr, double v) {
  contour_line *ptr=(contour_line *)vptr;
  ptr->level=v;
  return;
}

void *o2scl_contour_line_get_x(void *vptr) {
  contour_line *ptr=(contour_line *)vptr;
  return (void *)(&(ptr->x));
}

void o2scl_contour_line_set_x(void *vptr, void *p_v) {
  contour_line *ptr=(contour_line *)vptr;
  std::vector<double> *p_tsot=(std::vector<double> *)p_v;
  ptr->x=*(p_tsot);
  return;
}

void *o2scl_contour_line_get_y(void *vptr) {
  contour_line *ptr=(contour_line *)vptr;
  return (void *)(&(ptr->y));
}

void o2scl_contour_line_set_y(void *vptr, void *p_v) {
  contour_line *ptr=(contour_line *)vptr;
  std::vector<double> *p_tsot=(std::vector<double> *)p_v;
  ptr->y=*(p_tsot);
  return;
}

void *o2scl_create_std_vector_contour_line_() {
  std::vector<contour_line> *ptr=new std::vector<contour_line>;
  return ptr;
}

void o2scl_free_std_vector_contour_line_(void *vptr) {
  std::vector<contour_line> *ptr=(std::vector<contour_line> *)vptr;
  delete ptr;
  return;
}

void *o2scl_std_vector_contour_line__getitem(void *vptr, size_t n) {
  std::vector<contour_line> *ptr=(std::vector<contour_line> *)vptr;
  contour_line *ret=&(ptr->operator[](n));
  return ret;
}

void o2scl_std_vector_contour_line__setitem(void *vptr, size_t i, void *valptr) {
  std::vector<contour_line> *ptr=(std::vector<contour_line> *)vptr;
  contour_line *valptr2=(contour_line *)valptr;
  (*ptr)[i]=*valptr2;
  return;
}

void o2scl_std_vector_contour_line__resize(void *vptr, size_t n) {
  std::vector<contour_line> *ptr=(std::vector<contour_line> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std_vector_contour_line__size(void *vptr) {
  std::vector<contour_line> *ptr=(std::vector<contour_line> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

void *o2scl_create_contour() {
  contour *ptr=new contour;
  return ptr;
}

void o2scl_free_contour(void *vptr) {
  contour *ptr=(contour *)vptr;
  delete ptr;
  return;
}

int o2scl_contour_get_verbose(void *vptr) {
  contour *ptr=(contour *)vptr;
  return ptr->verbose;
}

void o2scl_contour_set_verbose(void *vptr, int v) {
  contour *ptr=(contour *)vptr;
  ptr->verbose=v;
  return;
}

double o2scl_contour_get_lev_adjust(void *vptr) {
  contour *ptr=(contour *)vptr;
  return ptr->lev_adjust;
}

void o2scl_contour_set_lev_adjust(void *vptr, double v) {
  contour *ptr=(contour *)vptr;
  ptr->lev_adjust=v;
  return;
}

bool o2scl_contour_get_debug_next_point(void *vptr) {
  contour *ptr=(contour *)vptr;
  return ptr->debug_next_point;
}

void o2scl_contour_set_debug_next_point(void *vptr, bool v) {
  contour *ptr=(contour *)vptr;
  ptr->debug_next_point=v;
  return;
}

void o2scl_contour_set_data(void *vptr, void *ptr_ugx, void *ptr_ugy, void *ptr_udata) {
  contour *ptr=(contour *)vptr;
  uniform_grid<double> *ugx=(uniform_grid<double> *)ptr_ugx;
  uniform_grid<double> *ugy=(uniform_grid<double> *)ptr_ugy;
  boost::numeric::ublas::matrix<double> *udata=(boost::numeric::ublas::matrix<double> *)ptr_udata;
  ptr->set_data(*ugx,*ugy,*udata);
  return;
}

void o2scl_contour_set_levels(void *vptr, size_t n_levels, void *ptr_levels) {
  contour *ptr=(contour *)vptr;
  vector<size_t> *levels=(vector<size_t> *)ptr_levels;
  ptr->set_levels(n_levels,*levels);
  return;
}

void o2scl_contour_calc_contours(void *vptr, void *ptr_clines) {
  contour *ptr=(contour *)vptr;
  vector<contour_line> *clines=(vector<contour_line> *)ptr_clines;
  ptr->calc_contours(*clines);
  return;
}

void *o2scl_create_prob_dens_mdim_std_vector_double_() {
  prob_dens_mdim<std::vector<double>> *ptr=new prob_dens_mdim<std::vector<double>>;
  return ptr;
}

void o2scl_free_prob_dens_mdim_std_vector_double_(void *vptr) {
  prob_dens_mdim<std::vector<double>> *ptr=(prob_dens_mdim<std::vector<double>> *)vptr;
  delete ptr;
  return;
}

double o2scl_prob_dens_mdim_std_vector_double__pdf(void *vptr, void *ptr_x) {
  prob_dens_mdim<std::vector<double>> *ptr=(prob_dens_mdim<std::vector<double>> *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  double ret=ptr->pdf(*x);
  return ret;
}

double o2scl_prob_dens_mdim_std_vector_double__log_pdf(void *vptr, void *ptr_x) {
  prob_dens_mdim<std::vector<double>> *ptr=(prob_dens_mdim<std::vector<double>> *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  double ret=ptr->log_pdf(*x);
  return ret;
}

size_t o2scl_prob_dens_mdim_std_vector_double__dim(void *vptr) {
  prob_dens_mdim<std::vector<double>> *ptr=(prob_dens_mdim<std::vector<double>> *)vptr;
  size_t ret=ptr->dim();
  return ret;
}

void o2scl_prob_dens_mdim_std_vector_double__getitem(void *vptr, void *ptr_x) {
  prob_dens_mdim<std::vector<double>> *ptr=(prob_dens_mdim<std::vector<double>> *)vptr;
  std::vector<double> *x=(std::vector<double> *)ptr_x;
  ptr->operator()(*x);
  return;
}

void *o2scl_create_prob_dens_mdim_biv_gaussian_std_vector_double_() {
  prob_dens_mdim_biv_gaussian<std::vector<double>> *ptr=new prob_dens_mdim_biv_gaussian<std::vector<double>>;
  return ptr;
}

void o2scl_free_prob_dens_mdim_biv_gaussian_std_vector_double_(void *vptr) {
  prob_dens_mdim_biv_gaussian<std::vector<double>> *ptr=(prob_dens_mdim_biv_gaussian<std::vector<double>> *)vptr;
  delete ptr;
  return;
}

void o2scl_prob_dens_mdim_biv_gaussian_std_vector_double__set(void *vptr, double x_cent, double y_cent, double x_std, double y_std, double covar) {
  prob_dens_mdim_biv_gaussian<std::vector<double>> *ptr=(prob_dens_mdim_biv_gaussian<std::vector<double>> *)vptr;
  ptr->set(x_cent,y_cent,x_std,y_std,covar);
  return;
}

void o2scl_prob_dens_mdim_biv_gaussian_std_vector_double__get(void *vptr, double *x_cent, double *y_cent, double *x_std, double *y_std, double *covar) {
  prob_dens_mdim_biv_gaussian<std::vector<double>> *ptr=(prob_dens_mdim_biv_gaussian<std::vector<double>> *)vptr;
  ptr->get(*x_cent,*y_cent,*x_std,*y_std,*covar);
  return;
}

double o2scl_prob_dens_mdim_biv_gaussian_std_vector_double__level_fixed_integral(void *vptr, double integral) {
  prob_dens_mdim_biv_gaussian<std::vector<double>> *ptr=(prob_dens_mdim_biv_gaussian<std::vector<double>> *)vptr;
  double ret=ptr->level_fixed_integral(integral);
  return ret;
}

void *o2scl_create_prob_dens_mdim_gaussian_() {
  prob_dens_mdim_gaussian<> *ptr=new prob_dens_mdim_gaussian<>;
  return ptr;
}

void o2scl_free_prob_dens_mdim_gaussian_(void *vptr) {
  prob_dens_mdim_gaussian<> *ptr=(prob_dens_mdim_gaussian<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_prob_dens_mdim_gaussian__make_biv(void *vptr) {
  prob_dens_mdim_gaussian<> *ptr=(prob_dens_mdim_gaussian<> *)vptr;
  prob_dens_mdim_biv_gaussian<> *ret=new prob_dens_mdim_biv_gaussian<>;
  *ret=ptr->make_biv();
  return ret;
}

void *o2scl_create_prob_dens_mdim_amr_hypercube() {
  prob_dens_mdim_amr<>::hypercube *ptr=new prob_dens_mdim_amr<>::hypercube;
  return ptr;
}

void o2scl_free_prob_dens_mdim_amr_hypercube(void *vptr) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  delete ptr;
  return;
}

size_t o2scl_prob_dens_mdim_amr_hypercube_get_n_dim(void *vptr) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  return ptr->n_dim;
}

void o2scl_prob_dens_mdim_amr_hypercube_set_n_dim(void *vptr, size_t v) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  ptr->n_dim=v;
  return;
}

void *o2scl_prob_dens_mdim_amr_hypercube_get_low(void *vptr) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  return (void *)(&(ptr->low));
}

void o2scl_prob_dens_mdim_amr_hypercube_set_low(void *vptr, void *p_v) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  std::vector<double> *p_tsot=(std::vector<double> *)p_v;
  ptr->low=*(p_tsot);
  return;
}

void *o2scl_prob_dens_mdim_amr_hypercube_get_high(void *vptr) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  return (void *)(&(ptr->high));
}

void o2scl_prob_dens_mdim_amr_hypercube_set_high(void *vptr, void *p_v) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  std::vector<double> *p_tsot=(std::vector<double> *)p_v;
  ptr->high=*(p_tsot);
  return;
}

void *o2scl_prob_dens_mdim_amr_hypercube_get_inside(void *vptr) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  return (void *)(&(ptr->inside));
}

void o2scl_prob_dens_mdim_amr_hypercube_set_inside(void *vptr, void *p_v) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  std::vector<size_t> *p_tsot=(std::vector<size_t> *)p_v;
  ptr->inside=*(p_tsot);
  return;
}

double o2scl_prob_dens_mdim_amr_hypercube_get_frac_vol(void *vptr) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  return ptr->frac_vol;
}

void o2scl_prob_dens_mdim_amr_hypercube_set_frac_vol(void *vptr, double v) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  ptr->frac_vol=v;
  return;
}

double o2scl_prob_dens_mdim_amr_hypercube_get_weight(void *vptr) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  return ptr->weight;
}

void o2scl_prob_dens_mdim_amr_hypercube_set_weight(void *vptr, double v) {
  prob_dens_mdim_amr<>::hypercube *ptr=(prob_dens_mdim_amr<>::hypercube *)vptr;
  ptr->weight=v;
  return;
}

void *o2scl_create_std_vector_prob_dens_mdim_amr_hypercube_() {
  std::vector<prob_dens_mdim_amr<>::hypercube> *ptr=new std::vector<prob_dens_mdim_amr<>::hypercube>;
  return ptr;
}

void o2scl_free_std_vector_prob_dens_mdim_amr_hypercube_(void *vptr) {
  std::vector<prob_dens_mdim_amr<>::hypercube> *ptr=(std::vector<prob_dens_mdim_amr<>::hypercube> *)vptr;
  delete ptr;
  return;
}

void o2scl_std_vector_prob_dens_mdim_amr_hypercube__resize(void *vptr, size_t n) {
  std::vector<prob_dens_mdim_amr<>::hypercube> *ptr=(std::vector<prob_dens_mdim_amr<>::hypercube> *)vptr;
  ptr->resize(n);
  return;
}

size_t o2scl_std_vector_prob_dens_mdim_amr_hypercube__size(void *vptr) {
  std::vector<prob_dens_mdim_amr<>::hypercube> *ptr=(std::vector<prob_dens_mdim_amr<>::hypercube> *)vptr;
  size_t ret=ptr->size();
  return ret;
}

void *o2scl_std_vector_prob_dens_mdim_amr_hypercube__getitem(void *vptr, size_t n) {
  std::vector<prob_dens_mdim_amr<>::hypercube> *ptr=(std::vector<prob_dens_mdim_amr<>::hypercube> *)vptr;
  prob_dens_mdim_amr<>::hypercube *ret=&(ptr->operator[](n));
  return ret;
}

void o2scl_std_vector_prob_dens_mdim_amr_hypercube__setitem(void *vptr, size_t i, void *valptr) {
  std::vector<prob_dens_mdim_amr<>::hypercube> *ptr=(std::vector<prob_dens_mdim_amr<>::hypercube> *)vptr;
  prob_dens_mdim_amr<>::hypercube *valptr2=(prob_dens_mdim_amr<>::hypercube *)valptr;
  (*ptr)[i]=*valptr2;
  return;
}

void *o2scl_create_prob_dens_mdim_amr_() {
  prob_dens_mdim_amr<> *ptr=new prob_dens_mdim_amr<>;
  return ptr;
}

void o2scl_free_prob_dens_mdim_amr_(void *vptr) {
  prob_dens_mdim_amr<> *ptr=(prob_dens_mdim_amr<> *)vptr;
  delete ptr;
  return;
}

void *o2scl_prob_dens_mdim_amr__get_mesh(void *vptr) {
  prob_dens_mdim_amr<> *ptr=(prob_dens_mdim_amr<> *)vptr;
  return (void *)(&(ptr->mesh));
}

void o2scl_prob_dens_mdim_amr__set_mesh(void *vptr, void *p_v) {
  prob_dens_mdim_amr<> *ptr=(prob_dens_mdim_amr<> *)vptr;
  std::vector<prob_dens_mdim_amr<>::hypercube> *p_tsot=(std::vector<prob_dens_mdim_amr<>::hypercube> *)p_v;
  ptr->mesh=*(p_tsot);
  return;
}

