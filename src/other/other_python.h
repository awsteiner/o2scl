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

#include <o2scl/slack_messenger.h>
#include <o2scl/poly.h>
#include <o2scl/polylog.h>
#include <o2scl/hist.h>
#include <o2scl/hist_2d.h>
#include <o2scl/contour.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/prob_dens_mdim_amr.h>

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

void o2scl_hist_create_rep_vec(void *vptr, void *ptr_v);

void *o2scl_hist_get_wgts(void *vptr);

void *o2scl_hist_get_bins(void *vptr);

void o2scl_hist_from_table(void *vptr, void *ptr_t, char *colx, size_t n_bins);

void o2scl_hist_from_table_twocol(void *vptr, void *ptr_t, char *colx, char *coly, size_t n_bins);

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

void *o2scl_create_hist_2d();

void o2scl_free_hist_2d(void *vptr);

void o2scl_copy_hist_2d(void *vsrc, void *vdest);

bool o2scl_hist_2d_get_extend_rhs(void *vptr);

void o2scl_hist_2d_set_extend_rhs(void *vptr, bool v);

bool o2scl_hist_2d_get_extend_lhs(void *vptr);

void o2scl_hist_2d_set_extend_lhs(void *vptr, bool v);

void o2scl_hist_2d_create_x_rep_vec(void *vptr, void *ptr_v);

void o2scl_hist_2d_create_y_rep_vec(void *vptr, void *ptr_v);

void o2scl_hist_2d_from_table(void *vptr, void *ptr_t, char *colx, char *coly, size_t n_bins_x, size_t n_bins_y);

void o2scl_hist_2d_from_table_wgt(void *vptr, void *ptr_t, char *colx, char *coly, char *colz, size_t n_bins_x, size_t n_bins_y);

size_t o2scl_hist_2d_size_x(void *vptr);

size_t o2scl_hist_2d_size_y(void *vptr);

void o2scl_hist_2d_set_bin_edges_grid(void *vptr, void *ptr_gx, void *ptr_gy);

void o2scl_hist_2d_set_bin_edges_vec(void *vptr, size_t nx, void *ptr_vx, size_t ny, void *ptr_vy);

void o2scl_hist_2d_update(void *vptr, double x, double y, double val=1.0);

void o2scl_hist_2d_update_i(void *vptr, size_t i, size_t j, double val=1.0);

double o2scl_hist_2d_get_wgt_i(void *vptr, size_t i, size_t j);

double o2scl_hist_2d_get_wgt(void *vptr, double x, double y);

void o2scl_hist_2d_set_wgt_i(void *vptr, size_t i, size_t j, double val);

double o2scl_hist_2d_get_x_low_i(void *vptr, size_t i);

double o2scl_hist_2d_get_x_high_i(void *vptr, size_t i);

double o2scl_hist_2d_get_y_low_i(void *vptr, size_t i);

double o2scl_hist_2d_get_y_high_i(void *vptr, size_t i);

void o2scl_hist_2d_set_wgt(void *vptr, double x, double y, double val);

void *o2scl_hist_2d_get_wgts(void *vptr);

void o2scl_hist_2d_clear(void *vptr);

void o2scl_hist_2d_clear_wgts(void *vptr);

void *o2scl_create_contour_line();

void o2scl_free_contour_line(void *vptr);

void o2scl_copy_contour_line(void *vsrc, void *vdest);

double o2scl_contour_line_get_level(void *vptr);

void o2scl_contour_line_set_level(void *vptr, double v);

void *o2scl_contour_line_get_x(void *vptr);

void o2scl_contour_line_set_x(void *vptr, void *p_v);

void *o2scl_contour_line_get_y(void *vptr);

void o2scl_contour_line_set_y(void *vptr, void *p_v);

void *o2scl_create_std_vector_contour_line_();

void o2scl_free_std_vector_contour_line_(void *vptr);

void *o2scl_std_vector_contour_line__getitem(void *vptr, size_t n);

void o2scl_std_vector_contour_line__setitem(void *vptr, size_t i, void *valptr);

void o2scl_std_vector_contour_line__resize(void *vptr, size_t n);

size_t o2scl_std_vector_contour_line__size(void *vptr);

void *o2scl_create_contour();

void o2scl_free_contour(void *vptr);

int o2scl_contour_get_verbose(void *vptr);

void o2scl_contour_set_verbose(void *vptr, int v);

double o2scl_contour_get_lev_adjust(void *vptr);

void o2scl_contour_set_lev_adjust(void *vptr, double v);

bool o2scl_contour_get_debug_next_point(void *vptr);

void o2scl_contour_set_debug_next_point(void *vptr, bool v);

void o2scl_contour_set_data(void *vptr, void *ptr_ugx, void *ptr_ugy, void *ptr_udata);

void o2scl_contour_set_levels(void *vptr, size_t n_levels, void *ptr_levels);

void o2scl_contour_calc_contours(void *vptr, void *ptr_clines);

void *o2scl_create_prob_dens_mdim_std_vector_double_();

void o2scl_free_prob_dens_mdim_std_vector_double_(void *vptr);

double o2scl_prob_dens_mdim_std_vector_double__pdf(void *vptr, void *ptr_x);

double o2scl_prob_dens_mdim_std_vector_double__log_pdf(void *vptr, void *ptr_x);

size_t o2scl_prob_dens_mdim_std_vector_double__dim(void *vptr);

void o2scl_prob_dens_mdim_std_vector_double__getitem(void *vptr, void *ptr_x);

void *o2scl_create_prob_dens_mdim_biv_gaussian_std_vector_double_();

void o2scl_free_prob_dens_mdim_biv_gaussian_std_vector_double_(void *vptr);

void o2scl_prob_dens_mdim_biv_gaussian_std_vector_double__set(void *vptr, double x_cent, double y_cent, double x_std, double y_std, double covar);

void o2scl_prob_dens_mdim_biv_gaussian_std_vector_double__get(void *vptr, double *x_cent, double *y_cent, double *x_std, double *y_std, double *covar);

double o2scl_prob_dens_mdim_biv_gaussian_std_vector_double__level_fixed_integral(void *vptr, double integral);

void *o2scl_create_prob_dens_mdim_gaussian_();

void o2scl_free_prob_dens_mdim_gaussian_(void *vptr);

void *o2scl_prob_dens_mdim_gaussian__make_biv(void *vptr);

void *o2scl_create_prob_dens_mdim_amr_hypercube();

void o2scl_free_prob_dens_mdim_amr_hypercube(void *vptr);

size_t o2scl_prob_dens_mdim_amr_hypercube_get_n_dim(void *vptr);

void o2scl_prob_dens_mdim_amr_hypercube_set_n_dim(void *vptr, size_t v);

void *o2scl_prob_dens_mdim_amr_hypercube_get_low(void *vptr);

void o2scl_prob_dens_mdim_amr_hypercube_set_low(void *vptr, void *p_v);

void *o2scl_prob_dens_mdim_amr_hypercube_get_high(void *vptr);

void o2scl_prob_dens_mdim_amr_hypercube_set_high(void *vptr, void *p_v);

void *o2scl_prob_dens_mdim_amr_hypercube_get_inside(void *vptr);

void o2scl_prob_dens_mdim_amr_hypercube_set_inside(void *vptr, void *p_v);

double o2scl_prob_dens_mdim_amr_hypercube_get_frac_vol(void *vptr);

void o2scl_prob_dens_mdim_amr_hypercube_set_frac_vol(void *vptr, double v);

double o2scl_prob_dens_mdim_amr_hypercube_get_weight(void *vptr);

void o2scl_prob_dens_mdim_amr_hypercube_set_weight(void *vptr, double v);

void *o2scl_create_std_vector_prob_dens_mdim_amr_hypercube_();

void o2scl_free_std_vector_prob_dens_mdim_amr_hypercube_(void *vptr);

void o2scl_std_vector_prob_dens_mdim_amr_hypercube__resize(void *vptr, size_t n);

size_t o2scl_std_vector_prob_dens_mdim_amr_hypercube__size(void *vptr);

void *o2scl_std_vector_prob_dens_mdim_amr_hypercube__getitem(void *vptr, size_t n);

void o2scl_std_vector_prob_dens_mdim_amr_hypercube__setitem(void *vptr, size_t i, void *valptr);

void *o2scl_create_prob_dens_mdim_amr_();

void o2scl_free_prob_dens_mdim_amr_(void *vptr);

void *o2scl_prob_dens_mdim_amr__get_mesh(void *vptr);

void o2scl_prob_dens_mdim_amr__set_mesh(void *vptr, void *p_v);

size_t o2scl_prob_dens_mdim_amr__get_n_dim(void *vptr);

void o2scl_prob_dens_mdim_amr__set_n_dim(void *vptr, size_t v);

void *o2scl_prob_dens_mdim_amr__get_low(void *vptr);

void o2scl_prob_dens_mdim_amr__set_low(void *vptr, void *p_v);

void *o2scl_prob_dens_mdim_amr__get_high(void *vptr);

void o2scl_prob_dens_mdim_amr__set_high(void *vptr, void *p_v);

bool o2scl_prob_dens_mdim_amr__get_allow_resampling(void *vptr);

void o2scl_prob_dens_mdim_amr__set_allow_resampling(void *vptr, bool v);

void *o2scl_prob_dens_mdim_amr__get_scale(void *vptr);

void o2scl_prob_dens_mdim_amr__set_scale(void *vptr, void *p_v);

int o2scl_prob_dens_mdim_amr__get_verbose(void *vptr);

void o2scl_prob_dens_mdim_amr__set_verbose(void *vptr, int v);

void o2scl_prob_dens_mdim_amr__clear(void *vptr);

void o2scl_prob_dens_mdim_amr__clear_mesh(void *vptr);

double o2scl_prob_dens_mdim_amr__total_volume(void *vptr);

}
