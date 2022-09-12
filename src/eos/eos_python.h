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

#include <o2scl/eos_base.h>
#include <o2scl/eos_had_base.h>
#include <o2scl/eos_had_skyrme.h>
#include <o2scl/eos_had_apr.h>
#include <o2scl/eos_had_rmf.h>
#include <o2scl/eos_quark_bag.h>
#include <o2scl/eos_quark_njl.h>
#include <o2scl/part.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>
#include <o2scl/hdf_eos_io.h>
#include <o2scl/nstar_cold.h>
#include <o2scl/tov_love.h>
#include <o2scl/eos_tov.h>
#include <o2scl/nucleus_rmf.h>

extern "C" {

void *o2scl_create_eos_base();

void o2scl_free_eos_base(void *vptr);

void *o2scl_eos_base_get_def_thermo(void *vptr);

void o2scl_eos_base_set_def_thermo(void *vptr, void *p_v);

double o2scl_eos_had_base_get_eoa(void *vptr);

void o2scl_eos_had_base_set_eoa(void *vptr, double v);

double o2scl_eos_had_base_get_msom(void *vptr);

void o2scl_eos_had_base_set_msom(void *vptr, double v);

double o2scl_eos_had_base_get_comp(void *vptr);

void o2scl_eos_had_base_set_comp(void *vptr, double v);

double o2scl_eos_had_base_get_n0(void *vptr);

void o2scl_eos_had_base_set_n0(void *vptr, double v);

double o2scl_eos_had_base_get_esym(void *vptr);

void o2scl_eos_had_base_set_esym(void *vptr, double v);

double o2scl_eos_had_base_get_kprime(void *vptr);

void o2scl_eos_had_base_set_kprime(void *vptr, double v);

bool o2scl_eos_had_base_get_err_nonconv(void *vptr);

void o2scl_eos_had_base_set_err_nonconv(void *vptr, bool v);

void *o2scl_eos_had_base_get_def_neutron(void *vptr);

void o2scl_eos_had_base_set_def_neutron(void *vptr, void *p_v);

void *o2scl_eos_had_base_get_def_proton(void *vptr);

void o2scl_eos_had_base_set_def_proton(void *vptr, void *p_v);

int o2scl_eos_had_base_calc_e(void *vptr, void *ptr_n, void *ptr_p, void *ptr_th);

int o2scl_eos_had_base_calc_p(void *vptr, void *ptr_n, void *ptr_p, void *ptr_th);

double o2scl_eos_had_base_fcomp(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_fcomp_err(void *vptr, double nb, double delta, double *unc);

double o2scl_eos_had_base_feoa(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_fesym(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_fesym_err(void *vptr, double nb, double delta, double *unc);

double o2scl_eos_had_base_fesym_slope(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_fesym_curve(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_fesym_skew(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_fesym_diff(void *vptr, double nb);

double o2scl_eos_had_base_feta(void *vptr, double nb);

double o2scl_eos_had_base_feta_prime(void *vptr, double nb);

double o2scl_eos_had_base_fkprime(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_fmsom(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_f_effm_neut(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_f_effm_prot(void *vptr, double nb, double delta=0.0);

double o2scl_eos_had_base_f_effm_scalar(void *vptr, double nb, double delta);

double o2scl_eos_had_base_f_effm_vector(void *vptr, double nb, double delta);

int o2scl_eos_had_base_fn0(void *vptr, double delta, double *nb, double *leoa);

void o2scl_eos_had_base_f_number_suscept(void *vptr, double mun, double mup, double *dPdnn, double *dPdnp, double *dPdpp);

void o2scl_eos_had_base_f_inv_number_suscept(void *vptr, double mun, double mup, double *dednn, double *dednp, double *dedpp);

int o2scl_eos_had_base_saturation(void *vptr);

double o2scl_eos_had_base_calc_mun_e(void *vptr, double nn, double np);

double o2scl_eos_had_base_calc_mup_e(void *vptr, double nn, double np);

double o2scl_eos_had_base_calc_ed(void *vptr, double nn, double np);

double o2scl_eos_had_base_calc_pr(void *vptr, double nn, double np);

double o2scl_eos_had_base_calc_nn_p(void *vptr, double mun, double mup);

double o2scl_eos_had_base_calc_np_p(void *vptr, double nn, double mup);

double o2scl_eos_had_base_calc_dmu_delta(void *vptr, double nb, double delta);

double o2scl_eos_had_base_calc_musum_delta(void *vptr, double nb, double delta);

double o2scl_eos_had_base_calc_pressure_nb(void *vptr, double nb, double delta);

double o2scl_eos_had_base_calc_edensity_nb(void *vptr, double nb, double delta);

int o2scl_eos_had_temp_base_calc_temp_e(void *vptr, void *ptr_n, void *ptr_p, double T, void *ptr_th);

int o2scl_eos_had_temp_base_calc_temp_p(void *vptr, void *ptr_n, void *ptr_p, double T, void *ptr_th);

void *o2scl_create_eos_had_skyrme();

void o2scl_free_eos_had_skyrme(void *vptr);

double o2scl_eos_had_skyrme_get_t0(void *vptr);

void o2scl_eos_had_skyrme_set_t0(void *vptr, double v);

double o2scl_eos_had_skyrme_get_t1(void *vptr);

void o2scl_eos_had_skyrme_set_t1(void *vptr, double v);

double o2scl_eos_had_skyrme_get_t2(void *vptr);

void o2scl_eos_had_skyrme_set_t2(void *vptr, double v);

double o2scl_eos_had_skyrme_get_t3(void *vptr);

void o2scl_eos_had_skyrme_set_t3(void *vptr, double v);

double o2scl_eos_had_skyrme_get_x0(void *vptr);

void o2scl_eos_had_skyrme_set_x0(void *vptr, double v);

double o2scl_eos_had_skyrme_get_x1(void *vptr);

void o2scl_eos_had_skyrme_set_x1(void *vptr, double v);

double o2scl_eos_had_skyrme_get_x2(void *vptr);

void o2scl_eos_had_skyrme_set_x2(void *vptr, double v);

double o2scl_eos_had_skyrme_get_x3(void *vptr);

void o2scl_eos_had_skyrme_set_x3(void *vptr, double v);

double o2scl_eos_had_skyrme_get_alpha(void *vptr);

void o2scl_eos_had_skyrme_set_alpha(void *vptr, double v);

double o2scl_eos_had_skyrme_get_a(void *vptr);

void o2scl_eos_had_skyrme_set_a(void *vptr, double v);

double o2scl_eos_had_skyrme_get_b(void *vptr);

void o2scl_eos_had_skyrme_set_b(void *vptr, double v);

double o2scl_eos_had_skyrme_get_W0(void *vptr);

void o2scl_eos_had_skyrme_set_W0(void *vptr, double v);

double o2scl_eos_had_skyrme_get_b4(void *vptr);

void o2scl_eos_had_skyrme_set_b4(void *vptr, double v);

double o2scl_eos_had_skyrme_get_b4p(void *vptr);

void o2scl_eos_had_skyrme_set_b4p(void *vptr, double v);

bool o2scl_eos_had_skyrme_get_parent_method(void *vptr);

void o2scl_eos_had_skyrme_set_parent_method(void *vptr, bool v);

void *o2scl_eos_had_skyrme_get_reference(void *vptr);

void o2scl_eos_had_skyrme_set_reference(void *vptr, void *p_v);

void *o2scl_eos_had_skyrme_get_nrfd(void *vptr);

void o2scl_eos_had_skyrme_set_nrfd(void *vptr, void *p_v);

void *o2scl_create_eos_had_apr();

void o2scl_free_eos_had_apr(void *vptr);

int o2scl_eos_had_apr_get_pion(void *vptr);

void o2scl_eos_had_apr_set_pion(void *vptr, int v);

bool o2scl_eos_had_apr_get_parent_method(void *vptr);

void o2scl_eos_had_apr_set_parent_method(void *vptr, bool v);

void *o2scl_create_eos_had_rmf();

void o2scl_free_eos_had_rmf(void *vptr);

size_t o2scl_eos_had_rmf_get_calc_e_steps(void *vptr);

void o2scl_eos_had_rmf_set_calc_e_steps(void *vptr, size_t v);

bool o2scl_eos_had_rmf_get_calc_e_relative(void *vptr);

void o2scl_eos_had_rmf_set_calc_e_relative(void *vptr, bool v);

bool o2scl_eos_had_rmf_get_zm_mode(void *vptr);

void o2scl_eos_had_rmf_set_zm_mode(void *vptr, bool v);

int o2scl_eos_had_rmf_get_verbose(void *vptr);

void o2scl_eos_had_rmf_set_verbose(void *vptr, int v);

double o2scl_eos_had_rmf_get_mnuc(void *vptr);

void o2scl_eos_had_rmf_set_mnuc(void *vptr, double v);

double o2scl_eos_had_rmf_get_ms(void *vptr);

void o2scl_eos_had_rmf_set_ms(void *vptr, double v);

double o2scl_eos_had_rmf_get_mw(void *vptr);

void o2scl_eos_had_rmf_set_mw(void *vptr, double v);

double o2scl_eos_had_rmf_get_mr(void *vptr);

void o2scl_eos_had_rmf_set_mr(void *vptr, double v);

double o2scl_eos_had_rmf_get_cs(void *vptr);

void o2scl_eos_had_rmf_set_cs(void *vptr, double v);

double o2scl_eos_had_rmf_get_cw(void *vptr);

void o2scl_eos_had_rmf_set_cw(void *vptr, double v);

double o2scl_eos_had_rmf_get_cr(void *vptr);

void o2scl_eos_had_rmf_set_cr(void *vptr, double v);

double o2scl_eos_had_rmf_get_b(void *vptr);

void o2scl_eos_had_rmf_set_b(void *vptr, double v);

double o2scl_eos_had_rmf_get_c(void *vptr);

void o2scl_eos_had_rmf_set_c(void *vptr, double v);

double o2scl_eos_had_rmf_get_zeta(void *vptr);

void o2scl_eos_had_rmf_set_zeta(void *vptr, double v);

double o2scl_eos_had_rmf_get_xi(void *vptr);

void o2scl_eos_had_rmf_set_xi(void *vptr, double v);

double o2scl_eos_had_rmf_get_a1(void *vptr);

void o2scl_eos_had_rmf_set_a1(void *vptr, double v);

double o2scl_eos_had_rmf_get_a2(void *vptr);

void o2scl_eos_had_rmf_set_a2(void *vptr, double v);

double o2scl_eos_had_rmf_get_a3(void *vptr);

void o2scl_eos_had_rmf_set_a3(void *vptr, double v);

double o2scl_eos_had_rmf_get_a4(void *vptr);

void o2scl_eos_had_rmf_set_a4(void *vptr, double v);

double o2scl_eos_had_rmf_get_a5(void *vptr);

void o2scl_eos_had_rmf_set_a5(void *vptr, double v);

double o2scl_eos_had_rmf_get_a6(void *vptr);

void o2scl_eos_had_rmf_set_a6(void *vptr, double v);

double o2scl_eos_had_rmf_get_b1(void *vptr);

void o2scl_eos_had_rmf_set_b1(void *vptr, double v);

double o2scl_eos_had_rmf_get_b2(void *vptr);

void o2scl_eos_had_rmf_set_b2(void *vptr, double v);

double o2scl_eos_had_rmf_get_b3(void *vptr);

void o2scl_eos_had_rmf_set_b3(void *vptr, double v);

int o2scl_eos_had_rmf_get_fields(void *vptr, double *sig, double *ome, double *rho);

int o2scl_eos_had_rmf_set_fields(void *vptr, double *sig, double *ome, double *rho);

void *o2scl_create_eos_quark();

void o2scl_free_eos_quark(void *vptr);

void *o2scl_create_eos_quark_bag();

void o2scl_free_eos_quark_bag(void *vptr);

double o2scl_eos_quark_bag_get_bag_constant(void *vptr);

void o2scl_eos_quark_bag_set_bag_constant(void *vptr, double v);

void *o2scl_create_eos_quark_njl();

void o2scl_free_eos_quark_njl(void *vptr);

double o2scl_eos_quark_njl_get_B0(void *vptr);

void o2scl_eos_quark_njl_set_B0(void *vptr, double v);

double o2scl_eos_quark_njl_get_L(void *vptr);

void o2scl_eos_quark_njl_set_L(void *vptr, double v);

double o2scl_eos_quark_njl_get_G(void *vptr);

void o2scl_eos_quark_njl_set_G(void *vptr, double v);

double o2scl_eos_quark_njl_get_K(void *vptr);

void o2scl_eos_quark_njl_set_K(void *vptr, double v);

double o2scl_eos_quark_njl_get_limit(void *vptr);

void o2scl_eos_quark_njl_set_limit(void *vptr, double v);

bool o2scl_eos_quark_njl_get_from_qq(void *vptr);

void o2scl_eos_quark_njl_set_from_qq(void *vptr, bool v);

int o2scl_eos_tov_get_verbose(void *vptr);

void o2scl_eos_tov_set_verbose(void *vptr, int v);

bool o2scl_eos_tov_has_baryons(void *vptr);

double o2scl_eos_tov_ed_from_pr(void *vptr, double pr);

double o2scl_eos_tov_pr_from_ed(void *vptr, double ed);

double o2scl_eos_tov_nb_from_ed(void *vptr, double ed);

double o2scl_eos_tov_nb_from_pr(void *vptr, double pr);

double o2scl_eos_tov_ed_from_nb(void *vptr, double nb);

double o2scl_eos_tov_pr_from_nb(void *vptr, double nb);

void o2scl_eos_tov_ed_nb_from_pr(void *vptr, double pr, double *ed, double *nb);

void *o2scl_create_eos_tov_buchdahl();

void o2scl_free_eos_tov_buchdahl(void *vptr);

double o2scl_eos_tov_buchdahl_get_Pstar(void *vptr);

void o2scl_eos_tov_buchdahl_set_Pstar(void *vptr, double v);

double o2scl_eos_tov_buchdahl_get_G_km_Msun(void *vptr);

void o2scl_eos_tov_buchdahl_set_G_km_Msun(void *vptr, double v);

void o2scl_eos_tov_buchdahl_set_baryon_density(void *vptr, double nb, double ed);

double o2scl_eos_tov_buchdahl_rad_from_gm(void *vptr, double gm);

double o2scl_eos_tov_buchdahl_ed_from_r_gm(void *vptr, double r, double beta);

double o2scl_eos_tov_buchdahl_pr_from_r_gm(void *vptr, double r, double beta);

double o2scl_eos_tov_buchdahl_exp2lam_from_r_gm(void *vptr, double r, double beta);

double o2scl_eos_tov_buchdahl_exp2phi_from_r_gm(void *vptr, double r, double beta);

void *o2scl_create_eos_tov_polytrope();

void o2scl_free_eos_tov_polytrope(void *vptr);

void o2scl_eos_tov_polytrope_set_coeff_index(void *vptr, double coeff, double index);

void *o2scl_create_eos_tov_linear();

void o2scl_free_eos_tov_linear(void *vptr);

void o2scl_eos_tov_linear_set_cs2_eps0(void *vptr, double cs2, double eps0);

void *o2scl_create_eos_tov_interp();

void o2scl_free_eos_tov_interp(void *vptr);

bool o2scl_eos_tov_interp_get_err_nonconv(void *vptr);

void o2scl_eos_tov_interp_set_err_nonconv(void *vptr, bool v);

void o2scl_eos_tov_interp_read_table(void *vptr, void *ptr_eos, char *s_cole, char *s_colp, char *s_colnb);

void o2scl_eos_tov_interp_default_low_dens_eos(void *vptr);

void o2scl_eos_tov_interp_sho11_low_dens_eos(void *vptr);

void o2scl_eos_tov_interp_s12_low_dens_eos(void *vptr, char *model, bool external=false);

void o2scl_eos_tov_interp_gcp10_low_dens_eos(void *vptr, char *model, bool external=false);

void o2scl_eos_tov_interp_ngl13_low_dens_eos(void *vptr, double L, char *model, bool external=false);

void o2scl_eos_tov_interp_ngl13_low_dens_eos2(void *vptr, double S, double L, double nt, char *fname);

void o2scl_eos_tov_interp_no_low_dens_eos(void *vptr);

void *o2scl_create_tov_solve();

void o2scl_free_tov_solve(void *vptr);

size_t o2scl_tov_solve_get_buffer_size(void *vptr);

void o2scl_tov_solve_set_buffer_size(void *vptr, size_t v);

size_t o2scl_tov_solve_get_max_table_size(void *vptr);

void o2scl_tov_solve_set_max_table_size(void *vptr, size_t v);

double o2scl_tov_solve_get_mass(void *vptr);

void o2scl_tov_solve_set_mass(void *vptr, double v);

double o2scl_tov_solve_get_rad(void *vptr);

void o2scl_tov_solve_set_rad(void *vptr, double v);

double o2scl_tov_solve_get_bmass(void *vptr);

void o2scl_tov_solve_set_bmass(void *vptr, double v);

double o2scl_tov_solve_get_gpot(void *vptr);

void o2scl_tov_solve_set_gpot(void *vptr, double v);

double o2scl_tov_solve_get_last_rjw(void *vptr);

void o2scl_tov_solve_set_last_rjw(void *vptr, double v);

double o2scl_tov_solve_get_last_f(void *vptr);

void o2scl_tov_solve_set_last_f(void *vptr, double v);

double o2scl_tov_solve_get_domega_rat(void *vptr);

void o2scl_tov_solve_set_domega_rat(void *vptr, double v);

double o2scl_tov_solve_get_pcent_max(void *vptr);

void o2scl_tov_solve_set_pcent_max(void *vptr, double v);

bool o2scl_tov_solve_get_reformat_results(void *vptr);

void o2scl_tov_solve_set_reformat_results(void *vptr, bool v);

double o2scl_tov_solve_get_baryon_mass(void *vptr);

void o2scl_tov_solve_set_baryon_mass(void *vptr, double v);

bool o2scl_tov_solve_get_ang_vel(void *vptr);

void o2scl_tov_solve_set_ang_vel(void *vptr, bool v);

bool o2scl_tov_solve_get_gen_rel(void *vptr);

void o2scl_tov_solve_set_gen_rel(void *vptr, bool v);

bool o2scl_tov_solve_get_calc_gpot(void *vptr);

void o2scl_tov_solve_set_calc_gpot(void *vptr, bool v);

double o2scl_tov_solve_get_step_min(void *vptr);

void o2scl_tov_solve_set_step_min(void *vptr, double v);

double o2scl_tov_solve_get_step_max(void *vptr);

void o2scl_tov_solve_set_step_max(void *vptr, double v);

double o2scl_tov_solve_get_step_start(void *vptr);

void o2scl_tov_solve_set_step_start(void *vptr, double v);

int o2scl_tov_solve_get_verbose(void *vptr);

void o2scl_tov_solve_set_verbose(void *vptr, int v);

size_t o2scl_tov_solve_get_max_integ_steps(void *vptr);

void o2scl_tov_solve_set_max_integ_steps(void *vptr, size_t v);

bool o2scl_tov_solve_get_err_nonconv(void *vptr);

void o2scl_tov_solve_set_err_nonconv(void *vptr, bool v);

double o2scl_tov_solve_get_pmax_default(void *vptr);

void o2scl_tov_solve_set_pmax_default(void *vptr, double v);

double o2scl_tov_solve_get_prbegin(void *vptr);

void o2scl_tov_solve_set_prbegin(void *vptr, double v);

double o2scl_tov_solve_get_prend(void *vptr);

void o2scl_tov_solve_set_prend(void *vptr, double v);

double o2scl_tov_solve_get_princ(void *vptr);

void o2scl_tov_solve_set_princ(void *vptr, double v);

double o2scl_tov_solve_get_fixed_pr_guess(void *vptr);

void o2scl_tov_solve_set_fixed_pr_guess(void *vptr, double v);

double o2scl_tov_solve_get_max_begin(void *vptr);

void o2scl_tov_solve_set_max_begin(void *vptr, double v);

double o2scl_tov_solve_get_max_end(void *vptr);

void o2scl_tov_solve_set_max_end(void *vptr, double v);

double o2scl_tov_solve_get_max_inc(void *vptr);

void o2scl_tov_solve_set_max_inc(void *vptr, double v);

void o2scl_tov_solve_set_eos(void *vptr, void *ptr_eos);

int o2scl_tov_solve_mvsr(void *vptr);

int o2scl_tov_solve_fixed(void *vptr, double mass, double pmax=1.0e20);

int o2scl_tov_solve_fixed_pr(void *vptr, double pcent, double pmax=1.0e20);

int o2scl_tov_solve_max(void *vptr);

void *o2scl_tov_solve_get_results(void *vptr);

void *o2scl_create_tov_love();

void o2scl_free_tov_love(void *vptr);

int o2scl_tov_love_get_show_ode(void *vptr);

void o2scl_tov_love_set_show_ode(void *vptr, int v);

bool o2scl_tov_love_get_addl_testing(void *vptr);

void o2scl_tov_love_set_addl_testing(void *vptr, bool v);

bool o2scl_tov_love_get_err_nonconv(void *vptr);

void o2scl_tov_love_set_err_nonconv(void *vptr, bool v);

void *o2scl_tov_love_get_results(void *vptr);

void o2scl_tov_love_set_results(void *vptr, void *p_v);

double o2scl_tov_love_get_delta(void *vptr);

void o2scl_tov_love_set_delta(void *vptr, double v);

double o2scl_tov_love_get_eps(void *vptr);

void o2scl_tov_love_set_eps(void *vptr, double v);

void o2scl_tov_love_get_tab(void *vptr, void *p_v);

void o2scl_tov_love_set_tab(void *vptr, void *p_v);

int o2scl_tov_love_calc_y(void *vptr, double *yR, double *beta, double *k2, double *lambda_km5, double *lambda_cgs, bool tabulate);

void o2scl_tov_love_add_disc(void *vptr, double rd);

void o2scl_tov_love_clear_discs(void *vptr);

int o2scl_tov_love_calc_H(void *vptr, double *yR, double *beta, double *k2, double *lambda_km5, double *lambda_cgs);

void *o2scl_create_nstar_cold();

void o2scl_free_nstar_cold(void *vptr);

double o2scl_nstar_cold_get_pressure_dec_nb(void *vptr);

void o2scl_nstar_cold_set_pressure_dec_nb(void *vptr, double v);

double o2scl_nstar_cold_get_allow_urca_nb(void *vptr);

void o2scl_nstar_cold_set_allow_urca_nb(void *vptr, double v);

double o2scl_nstar_cold_get_deny_urca_nb(void *vptr);

void o2scl_nstar_cold_set_deny_urca_nb(void *vptr, double v);

double o2scl_nstar_cold_get_acausal_nb(void *vptr);

void o2scl_nstar_cold_set_acausal_nb(void *vptr, double v);

double o2scl_nstar_cold_get_acausal_ed(void *vptr);

void o2scl_nstar_cold_set_acausal_ed(void *vptr, double v);

double o2scl_nstar_cold_get_acausal_pr(void *vptr);

void o2scl_nstar_cold_set_acausal_pr(void *vptr, double v);

void *o2scl_nstar_cold_get_def_tov(void *vptr);


bool o2scl_nstar_cold_get_eos_neg(void *vptr);

void o2scl_nstar_cold_set_eos_neg(void *vptr, bool v);

int o2scl_nstar_cold_get_verbose(void *vptr);

void o2scl_nstar_cold_set_verbose(void *vptr, int v);

double o2scl_nstar_cold_get_nb_start(void *vptr);

void o2scl_nstar_cold_set_nb_start(void *vptr, double v);

double o2scl_nstar_cold_get_nb_end(void *vptr);

void o2scl_nstar_cold_set_nb_end(void *vptr, double v);

double o2scl_nstar_cold_get_dnb(void *vptr);

void o2scl_nstar_cold_set_dnb(void *vptr, double v);

size_t o2scl_nstar_cold_get_max_row(void *vptr);

void o2scl_nstar_cold_set_max_row(void *vptr, size_t v);

bool o2scl_nstar_cold_get_remove_rows(void *vptr);

void o2scl_nstar_cold_set_remove_rows(void *vptr, bool v);

bool o2scl_nstar_cold_get_include_muons(void *vptr);

void o2scl_nstar_cold_set_include_muons(void *vptr, bool v);

bool o2scl_nstar_cold_get_err_nonconv(void *vptr);

void o2scl_nstar_cold_set_err_nonconv(void *vptr, bool v);

void o2scl_nstar_cold_set_eos(void *vptr, void *ptr_eos);

int o2scl_nstar_cold_calc_eos(void *vptr, double np_0=0.0);

int o2scl_nstar_cold_calc_nstar(void *vptr);

int o2scl_nstar_cold_fixed(void *vptr, double target_mass);

void *o2scl_nstar_cold_get_eos_results(void *vptr);

void *o2scl_nstar_cold_get_tov_results(void *vptr);

void *o2scl_create_nucleus_rmf();

void o2scl_free_nucleus_rmf(void *vptr);

double o2scl_nucleus_rmf_get_stens(void *vptr);

void o2scl_nucleus_rmf_set_stens(void *vptr, double v);

double o2scl_nucleus_rmf_get_rnrp(void *vptr);

void o2scl_nucleus_rmf_set_rnrp(void *vptr, double v);

double o2scl_nucleus_rmf_get_rnrms(void *vptr);

void o2scl_nucleus_rmf_set_rnrms(void *vptr, double v);

double o2scl_nucleus_rmf_get_rprms(void *vptr);

void o2scl_nucleus_rmf_set_rprms(void *vptr, double v);

double o2scl_nucleus_rmf_get_etot(void *vptr);

void o2scl_nucleus_rmf_set_etot(void *vptr, double v);

double o2scl_nucleus_rmf_get_r_charge(void *vptr);

void o2scl_nucleus_rmf_set_r_charge(void *vptr, double v);

double o2scl_nucleus_rmf_get_r_charge_cm(void *vptr);

void o2scl_nucleus_rmf_set_r_charge_cm(void *vptr, double v);

int o2scl_nucleus_rmf_run_nucleus(void *vptr, int nucleus_Z, int nucleus_N, int unocc_Z, int unocc_N);

void *o2scl_nucleus_rmf_get_profiles(void *vptr);

void *o2scl_nucleus_rmf_get_chden(void *vptr);

void o2scl_skyrme_load_wrapper(void *ptr_sk, char *model, bool external=false, int verbose=0);

void o2scl_rmf_load_wrapper(void *ptr_rmf, char *model, bool external=false);

}
