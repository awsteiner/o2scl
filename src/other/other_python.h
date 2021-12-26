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

#include <o2scl/slack_messenger.h>
#include <o2scl/poly.h>
#include <o2scl/polylog.h>
#include <o2scl/hist.h>

extern "C" {

void *o2scl_create_slack_messenger();

void o2scl_free_slack_messenger(void *vptr);

int o2scl_slack_messenger_get_verbose(void *vptr);

void o2scl_slack_messenger_set_verbose(void *vptr, int v);

void *o2scl_slack_messenger_get_url(void *vptr);

void o2scl_slack_messenger_set_url(void *vptr, void *p_v);

void *o2scl_slack_messenger_get_channel(void *vptr);

void o2scl_slack_messenger_set_channel(void *vptr, void *p_v);

void *o2scl_slack_messenger_get_icon(void *vptr);

void o2scl_slack_messenger_set_icon(void *vptr, void *p_v);

void *o2scl_slack_messenger_get_username(void *vptr);

void o2scl_slack_messenger_set_username(void *vptr, void *p_v);

double o2scl_slack_messenger_get_min_time_between(void *vptr);

void o2scl_slack_messenger_set_min_time_between(void *vptr, double v);

bool o2scl_slack_messenger_set_url_from_env(void *vptr, char *env_var);

bool o2scl_slack_messenger_set_channel_from_env(void *vptr, char *env_var);

bool o2scl_slack_messenger_set_username_from_env(void *vptr, char *env_var);

int o2scl_slack_messenger_send(void *vptr, char *message, bool err_on_fail=true);

void *o2scl_slack_messenger_init(char *p_channel, char *p_username, char *p_url, bool p_mpi_time);

void *o2scl_create_quadratic_real_coeff_gsl();

void o2scl_free_quadratic_real_coeff_gsl(void *vptr);

int o2scl_quadratic_real_coeff_gsl_solve_r(void *vptr, double a2, double b2, double c2, double *r1, double *r2);

int o2scl_quadratic_real_coeff_gsl_solve_rc(void *vptr, double a2, double b2, double c2, void *ptr_r1, void *ptr_r2);

void *o2scl_create_quadratic_real_coeff_gsl2_();

void o2scl_free_quadratic_real_coeff_gsl2_(void *vptr);

int o2scl_quadratic_real_coeff_gsl2__solve_r(void *vptr, double a2, double b2, double c2, double *r1, double *r2);

int o2scl_quadratic_real_coeff_gsl2__solve_rc(void *vptr, double a2, double b2, double c2, void *ptr_r1, void *ptr_r2);

void *o2scl_create_cubic_real_coeff_cern_();

void o2scl_free_cubic_real_coeff_cern_(void *vptr);

int o2scl_cubic_real_coeff_cern__solve_r(void *vptr, double a3, double b3, double c3, double d3, double *r1, double *r2, double *r3);

int o2scl_cubic_real_coeff_cern__solve_rc(void *vptr, double a3, double b3, double c3, double d3, double *r1, void *ptr_r2, void *ptr_r3);

void *o2scl_create_cubic_real_coeff_gsl();

void o2scl_free_cubic_real_coeff_gsl(void *vptr);

int o2scl_cubic_real_coeff_gsl_solve_r(void *vptr, double a3, double b3, double c3, double d3, double *r1, double *r2, double *r3);

int o2scl_cubic_real_coeff_gsl_solve_rc(void *vptr, double a3, double b3, double c3, double d3, double *r1, void *ptr_r2, void *ptr_r3);

void *o2scl_create_quartic_real_coeff_cern_();

void o2scl_free_quartic_real_coeff_cern_(void *vptr);

int o2scl_quartic_real_coeff_cern__solve_r(void *vptr, double a4, double b4, double c4, double d4, double e4, double *r1, double *r2, double *r3, double *r4);

int o2scl_quartic_real_coeff_cern__solve_rc(void *vptr, double a4, double b4, double c4, double d4, double e4, void *ptr_r1, void *ptr_r2, void *ptr_r3, void *ptr_r4);

void *o2scl_create_fermi_dirac_integ_gsl();

void o2scl_free_fermi_dirac_integ_gsl(void *vptr);

double o2scl_fermi_dirac_integ_gsl_calc_m1o2(void *vptr, double x);

double o2scl_fermi_dirac_integ_gsl_calc_1o2(void *vptr, double x);

double o2scl_fermi_dirac_integ_gsl_calc_3o2(void *vptr, double x);

double o2scl_fermi_dirac_integ_gsl_calc_2(void *vptr, double x);

double o2scl_fermi_dirac_integ_gsl_calc_3(void *vptr, double x);

void *o2scl_create_bessel_K_exp_integ_gsl();

void o2scl_free_bessel_K_exp_integ_gsl(void *vptr);

double o2scl_bessel_K_exp_integ_gsl_K1exp(void *vptr, double x);

double o2scl_bessel_K_exp_integ_gsl_K2exp(void *vptr, double x);

double o2scl_bessel_K_exp_integ_gsl_K3exp(void *vptr, double x);

void *o2scl_create_hist();

void o2scl_free_hist(void *vptr);

void o2scl_copy_hist(void *vsrc, void *vdest);

bool o2scl_hist_get_extend_rhs(void *vptr);

void o2scl_hist_set_extend_rhs(void *vptr, bool v);

bool o2scl_hist_get_extend_lhs(void *vptr);

void o2scl_hist_set_extend_lhs(void *vptr, bool v);

size_t o2scl_hist_size(void *vptr);

void o2scl_hist_set_bin_edges_grid(void *vptr, void *ptr_g);

void o2scl_hist_set_bin_edges_vec(void *vptr, size_t n, void *ptr_v);

void o2scl_hist_update(void *vptr, double x, double val=1.0);

void o2scl_hist_update_i(void *vptr, size_t i, double val=1.0);

double o2scl_hist_get_wgt_i(void *vptr, size_t i);

double o2scl_hist_get_wgt(void *vptr, double x);

void o2scl_hist_set_wgt_i(void *vptr, size_t i, double x);

void o2scl_hist_set_wgt(void *vptr, double x, double val);

const double o2scl_hist_getitem(void *vptr, size_t i);

size_t o2scl_hist_get_bin_index(void *vptr, double x);

int o2scl_hist_function(void *vptr, char *func);

void o2scl_hist_clear(void *vptr);

}
