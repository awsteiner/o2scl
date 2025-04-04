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

#include <o2scl/eos_base.h>
#include <o2scl/eos_leptons.h>
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
#include <o2scl/nstar_rot.h>
#include <o2scl/tov_love.h>
#include <o2scl/eos_tov.h>
#include <o2scl/nucleus_rmf.h>
#include <o2scl/nucmass_ldrop.h>
#include <o2scl/nucmass_ldrop_shell.h>
#include <o2scl/nucleus_bin.h>

extern "C" {

void *o2scl_eos_base_get_def_thermo(void *vptr);

void o2scl_eos_base_set_def_thermo(void *vptr, void *p_v);

void *o2scl_create_eos_leptons();

void o2scl_free_eos_leptons(void *vptr);

void *o2scl_eos_leptons_get_th(void *vptr);

void o2scl_eos_leptons_set_th(void *vptr, void *p_v);

void *o2scl_eos_leptons_get_e(void *vptr);

void o2scl_eos_leptons_set_e(void *vptr, void *p_v);

void *o2scl_eos_leptons_get_mu(void *vptr);

void o2scl_eos_leptons_set_mu(void *vptr, void *p_v);

void *o2scl_eos_leptons_get_ph(void *vptr);

void o2scl_eos_leptons_set_ph(void *vptr, void *p_v);

void *o2scl_eos_leptons_get_ed(void *vptr);

void o2scl_eos_leptons_set_ed(void *vptr, void *p_v);

void *o2scl_eos_leptons_get_mud(void *vptr);

void o2scl_eos_leptons_set_mud(void *vptr, void *p_v);

void *o2scl_eos_leptons_get_phd(void *vptr);

void o2scl_eos_leptons_set_phd(void *vptr, void *p_v);

bool o2scl_eos_leptons_get_include_muons(void *vptr);

void o2scl_eos_leptons_set_include_muons(void *vptr, bool v);

bool o2scl_eos_leptons_get_include_deriv(void *vptr);

void o2scl_eos_leptons_set_include_deriv(void *vptr, bool v);

bool o2scl_eos_leptons_get_pde_from_density(void *vptr);

void o2scl_eos_leptons_set_pde_from_density(void *vptr, bool v);

int o2scl_eos_leptons_get_verbose(void *vptr);

void o2scl_eos_leptons_set_verbose(void *vptr, int v);

bool o2scl_eos_leptons_get_err_nonconv(void *vptr);

void o2scl_eos_leptons_set_err_nonconv(void *vptr, bool v);

void *o2scl_eos_leptons_get_frel(void *vptr);

void o2scl_eos_leptons_set_frel(void *vptr, void *p_v);

void o2scl_eos_leptons_default_acc(void *vptr);

void o2scl_eos_leptons_improved_acc(void *vptr);

void o2scl_eos_leptons_ld_acc(void *vptr);

void o2scl_eos_leptons_fp_25_acc(void *vptr);

int o2scl_eos_leptons_pair_mu(void *vptr, double T);

int o2scl_eos_leptons_pair_mu_eq(void *vptr, double T);

int o2scl_eos_leptons_pair_density(void *vptr, double T);

int o2scl_eos_leptons_pair_density_eq(void *vptr, double nq, double T);

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

void *o2scl_eos_tov_interp_get_full_vece(void *vptr);

void o2scl_eos_tov_interp_set_full_vece(void *vptr, void *p_v);

void *o2scl_eos_tov_interp_get_full_vecp(void *vptr);

void o2scl_eos_tov_interp_set_full_vecp(void *vptr, void *p_v);

void *o2scl_eos_tov_interp_get_full_vecnb(void *vptr);

void o2scl_eos_tov_interp_set_full_vecnb(void *vptr, void *p_v);

void o2scl_eos_tov_interp_read_table(void *vptr, void *ptr_eos, void *ptr_s_cole, void *ptr_s_colp, void *ptr_s_colnb);

void o2scl_eos_tov_interp_default_low_dens_eos(void *vptr);

void o2scl_eos_tov_interp_sho11_low_dens_eos(void *vptr);

void o2scl_eos_tov_interp_s12_low_dens_eos(void *vptr, void *ptr_model, bool external=false);

void o2scl_eos_tov_interp_gcp10_low_dens_eos(void *vptr, void *ptr_model, bool external=false);

void o2scl_eos_tov_interp_ngl13_low_dens_eos(void *vptr, double L, void *ptr_model, bool external=false);

void o2scl_eos_tov_interp_ngl13_low_dens_eos2(void *vptr, double S, double L, double nt, void *ptr_fname);

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

double o2scl_eos_nstar_rot_enth_from_pr(void *vptr, double pr);

double o2scl_eos_nstar_rot_pr_from_enth(void *vptr, double pr);

double o2scl_eos_nstar_rot_enth_from_nb(void *vptr, double pr);

void *o2scl_create_eos_nstar_rot_interp();

void o2scl_free_eos_nstar_rot_interp(void *vptr);

void o2scl_eos_nstar_rot_interp_set_eos_fm(void *vptr, size_t n, void *ptr_eden, void *ptr_pres, void *ptr_nb);

void *o2scl_create_eos_nstar_rot_C();

void o2scl_free_eos_nstar_rot_C(void *vptr);

void o2scl_eos_nstar_rot_C_set(void *vptr, bool rns_constants=false);

void *o2scl_create_eos_nstar_rot_L();

void o2scl_free_eos_nstar_rot_L(void *vptr);

void o2scl_eos_nstar_rot_L_set(void *vptr, bool rns_constants=false);

void *o2scl_create_nstar_rot();

void o2scl_free_nstar_rot(void *vptr);

int o2scl_nstar_rot_get_verbose(void *vptr);

void o2scl_nstar_rot_set_verbose(void *vptr, int v);

double o2scl_nstar_rot_get_Mass(void *vptr);

void o2scl_nstar_rot_set_Mass(void *vptr, double v);

double o2scl_nstar_rot_get_MSUN(void *vptr);

void o2scl_nstar_rot_set_MSUN(void *vptr, double v);

double o2scl_nstar_rot_get_R_e(void *vptr);

void o2scl_nstar_rot_set_R_e(void *vptr, double v);

double o2scl_nstar_rot_get_Omega(void *vptr);

void o2scl_nstar_rot_set_Omega(void *vptr, double v);

double o2scl_nstar_rot_get_Omega_K(void *vptr);

void o2scl_nstar_rot_set_Omega_K(void *vptr, double v);

double o2scl_nstar_rot_get_T(void *vptr);

void o2scl_nstar_rot_set_T(void *vptr, double v);

double o2scl_nstar_rot_get_W(void *vptr);

void o2scl_nstar_rot_set_W(void *vptr, double v);

double o2scl_nstar_rot_get_J(void *vptr);

void o2scl_nstar_rot_set_J(void *vptr, double v);

double o2scl_nstar_rot_get_C(void *vptr);

void o2scl_nstar_rot_set_C(void *vptr, double v);

double o2scl_nstar_rot_get_G(void *vptr);

void o2scl_nstar_rot_set_G(void *vptr, double v);

double o2scl_nstar_rot_get_I(void *vptr);

void o2scl_nstar_rot_set_I(void *vptr, double v);

double o2scl_nstar_rot_get_mass_quadrupole(void *vptr);

void o2scl_nstar_rot_set_mass_quadrupole(void *vptr, double v);

double o2scl_nstar_rot_get_KAPPA(void *vptr);

void o2scl_nstar_rot_set_KAPPA(void *vptr, double v);

double o2scl_nstar_rot_get_h_plus(void *vptr);

void o2scl_nstar_rot_set_h_plus(void *vptr, double v);

double o2scl_nstar_rot_get_h_minus(void *vptr);

void o2scl_nstar_rot_set_h_minus(void *vptr, double v);

double o2scl_nstar_rot_get_Z_p(void *vptr);

void o2scl_nstar_rot_set_Z_p(void *vptr, double v);

double o2scl_nstar_rot_get_Z_f(void *vptr);

void o2scl_nstar_rot_set_Z_f(void *vptr, double v);

double o2scl_nstar_rot_get_Z_b(void *vptr);

void o2scl_nstar_rot_set_Z_b(void *vptr, double v);

double o2scl_nstar_rot_get_om_over_Om(void *vptr);

void o2scl_nstar_rot_set_om_over_Om(void *vptr, double v);

double o2scl_nstar_rot_get_r_e(void *vptr);

void o2scl_nstar_rot_set_r_e(void *vptr, double v);

double o2scl_nstar_rot_get_r_ratio(void *vptr);

void o2scl_nstar_rot_set_r_ratio(void *vptr, double v);

void o2scl_nstar_rot_output_table(void *vptr, void *ptr_t);

void o2scl_nstar_rot_set_eos(void *vptr, void *ptr_eos);

void o2scl_nstar_rot_polytrope_eos(void *vptr, double index);

void o2scl_nstar_rot_fix_cent_eden_axis_rat(void *vptr, double cent_eden, double axis_rat, bool use_guess=false);

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

void *o2scl_create_nucmass_ldrop();

void o2scl_free_nucmass_ldrop(void *vptr);

double o2scl_nucmass_ldrop_get_n1(void *vptr);

void o2scl_nucmass_ldrop_set_n1(void *vptr, double v);

double o2scl_nucmass_ldrop_get_n0(void *vptr);

void o2scl_nucmass_ldrop_set_n0(void *vptr, double v);

double o2scl_nucmass_ldrop_get_surften(void *vptr);

void o2scl_nucmass_ldrop_set_surften(void *vptr, double v);

double o2scl_nucmass_ldrop_get_coul_coeff(void *vptr);

void o2scl_nucmass_ldrop_set_coul_coeff(void *vptr, double v);

double o2scl_nucmass_ldrop_get_nn(void *vptr);

void o2scl_nucmass_ldrop_set_nn(void *vptr, double v);

double o2scl_nucmass_ldrop_get_np(void *vptr);

void o2scl_nucmass_ldrop_set_np(void *vptr, double v);

double o2scl_nucmass_ldrop_get_Rn(void *vptr);

void o2scl_nucmass_ldrop_set_Rn(void *vptr, double v);

double o2scl_nucmass_ldrop_get_Rp(void *vptr);

void o2scl_nucmass_ldrop_set_Rp(void *vptr, double v);

double o2scl_nucmass_ldrop_get_surf(void *vptr);

void o2scl_nucmass_ldrop_set_surf(void *vptr, double v);

double o2scl_nucmass_ldrop_get_bulk(void *vptr);

void o2scl_nucmass_ldrop_set_bulk(void *vptr, double v);

double o2scl_nucmass_ldrop_get_coul(void *vptr);

void o2scl_nucmass_ldrop_set_coul(void *vptr, double v);

bool o2scl_nucmass_ldrop_get_large_vals_unphys(void *vptr);

void o2scl_nucmass_ldrop_set_large_vals_unphys(void *vptr, bool v);

void *o2scl_nucmass_ldrop_get_def_had_eos(void *vptr);

void o2scl_nucmass_ldrop_set_def_had_eos(void *vptr, void *p_v);

void *o2scl_nucmass_ldrop_get_def_neutron(void *vptr);

void o2scl_nucmass_ldrop_set_def_neutron(void *vptr, void *p_v);

void *o2scl_nucmass_ldrop_get_def_proton(void *vptr);

void o2scl_nucmass_ldrop_set_def_proton(void *vptr, void *p_v);

void *o2scl_nucmass_ldrop_get_th(void *vptr);

void o2scl_nucmass_ldrop_set_th(void *vptr, void *p_v);

double o2scl_nucmass_ldrop_mass_excess_d(void *vptr, double Z, double N);

double o2scl_nucmass_ldrop_mass_excess(void *vptr, int Z, int N);

double o2scl_nucmass_ldrop_binding_energy_densmat(void *vptr, double Z, double N, double npout, double nnout, double neout, double T);

void o2scl_nucmass_ldrop_set_n_and_p(void *vptr, void *ptr_un, void *ptr_up);

int o2scl_nucmass_ldrop_set_eos_had_temp_base(void *vptr, void *ptr_uhe);

void *o2scl_create_nucmass_ldrop_skin();

void o2scl_free_nucmass_ldrop_skin(void *vptr);

bool o2scl_nucmass_ldrop_skin_get_full_surface(void *vptr);

void o2scl_nucmass_ldrop_skin_set_full_surface(void *vptr, bool v);

bool o2scl_nucmass_ldrop_skin_get_new_skin_mode(void *vptr);

void o2scl_nucmass_ldrop_skin_set_new_skin_mode(void *vptr, bool v);

double o2scl_nucmass_ldrop_skin_get_doi(void *vptr);

void o2scl_nucmass_ldrop_skin_set_doi(void *vptr, double v);

double o2scl_nucmass_ldrop_skin_get_ss(void *vptr);

void o2scl_nucmass_ldrop_skin_set_ss(void *vptr, double v);

double o2scl_nucmass_ldrop_skin_get_pp(void *vptr);

void o2scl_nucmass_ldrop_skin_set_pp(void *vptr, double v);

double o2scl_nucmass_ldrop_skin_get_a0(void *vptr);

void o2scl_nucmass_ldrop_skin_set_a0(void *vptr, double v);

double o2scl_nucmass_ldrop_skin_get_a2(void *vptr);

void o2scl_nucmass_ldrop_skin_set_a2(void *vptr, double v);

double o2scl_nucmass_ldrop_skin_get_a4(void *vptr);

void o2scl_nucmass_ldrop_skin_set_a4(void *vptr, double v);

bool o2scl_nucmass_ldrop_skin_get_rel_vacuum(void *vptr);

void o2scl_nucmass_ldrop_skin_set_rel_vacuum(void *vptr, bool v);

double o2scl_nucmass_ldrop_skin_get_Tchalf(void *vptr);

void o2scl_nucmass_ldrop_skin_set_Tchalf(void *vptr, double v);

void *o2scl_create_nucmass_ldrop_pair();

void o2scl_free_nucmass_ldrop_pair(void *vptr);

double o2scl_nucmass_ldrop_pair_get_Epair(void *vptr);

void o2scl_nucmass_ldrop_pair_set_Epair(void *vptr, double v);

double o2scl_nucmass_ldrop_pair_get_pair(void *vptr);

void o2scl_nucmass_ldrop_pair_set_pair(void *vptr, double v);

void *o2scl_create_nucmass_ldrop_shell();

void o2scl_free_nucmass_ldrop_shell(void *vptr);

void *o2scl_create_nucmass_frdm_shell();

void o2scl_free_nucmass_frdm_shell(void *vptr);

void *o2scl_create_nucleus_bin();

void o2scl_free_nucleus_bin(void *vptr);

void *o2scl_nucleus_bin_get_ame16(void *vptr);

void o2scl_nucleus_bin_set_ame16(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ame20exp(void *vptr);

void o2scl_nucleus_bin_set_ame20exp(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ame20round(void *vptr);

void o2scl_nucleus_bin_set_ame20round(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ame95rmd(void *vptr);

void o2scl_nucleus_bin_set_ame95rmd(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ame03round(void *vptr);

void o2scl_nucleus_bin_set_ame03round(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ame03(void *vptr);

void o2scl_nucleus_bin_set_ame03(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ame95exp(void *vptr);

void o2scl_nucleus_bin_set_ame95exp(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ame12(void *vptr);

void o2scl_nucleus_bin_set_ame12(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ddme2(void *vptr);

void o2scl_nucleus_bin_set_ddme2(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ddmed(void *vptr);

void o2scl_nucleus_bin_set_ddmed(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_ddpc1(void *vptr);

void o2scl_nucleus_bin_set_ddpc1(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_nl3s(void *vptr);

void o2scl_nucleus_bin_set_nl3s(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_sly4(void *vptr);

void o2scl_nucleus_bin_set_sly4(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_skms(void *vptr);

void o2scl_nucleus_bin_set_skms(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_skp(void *vptr);

void o2scl_nucleus_bin_set_skp(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_sv_min(void *vptr);

void o2scl_nucleus_bin_set_sv_min(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_unedf0(void *vptr);

void o2scl_nucleus_bin_set_unedf0(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_unedf1(void *vptr);

void o2scl_nucleus_bin_set_unedf1(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_m95(void *vptr);

void o2scl_nucleus_bin_set_m95(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_m16(void *vptr);

void o2scl_nucleus_bin_set_m16(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_kt(void *vptr);

void o2scl_nucleus_bin_set_kt(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_kt2(void *vptr);

void o2scl_nucleus_bin_set_kt2(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_wlw1(void *vptr);

void o2scl_nucleus_bin_set_wlw1(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_wlw2(void *vptr);

void o2scl_nucleus_bin_set_wlw2(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_wlw3(void *vptr);

void o2scl_nucleus_bin_set_wlw3(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_wlw4(void *vptr);

void o2scl_nucleus_bin_set_wlw4(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_wlw5(void *vptr);

void o2scl_nucleus_bin_set_wlw5(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_sdnp1(void *vptr);

void o2scl_nucleus_bin_set_sdnp1(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_sdnp2(void *vptr);

void o2scl_nucleus_bin_set_sdnp2(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_sdnp3(void *vptr);

void o2scl_nucleus_bin_set_sdnp3(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_dz(void *vptr);

void o2scl_nucleus_bin_set_dz(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb2(void *vptr);

void o2scl_nucleus_bin_set_hfb2(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb8(void *vptr);

void o2scl_nucleus_bin_set_hfb8(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb14(void *vptr);

void o2scl_nucleus_bin_set_hfb14(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb14_v0(void *vptr);

void o2scl_nucleus_bin_set_hfb14_v0(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb17(void *vptr);

void o2scl_nucleus_bin_set_hfb17(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb21(void *vptr);

void o2scl_nucleus_bin_set_hfb21(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb22(void *vptr);

void o2scl_nucleus_bin_set_hfb22(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb23(void *vptr);

void o2scl_nucleus_bin_set_hfb23(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb24(void *vptr);

void o2scl_nucleus_bin_set_hfb24(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb25(void *vptr);

void o2scl_nucleus_bin_set_hfb25(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb26(void *vptr);

void o2scl_nucleus_bin_set_hfb26(void *vptr, void *p_v);

void *o2scl_nucleus_bin_get_hfb27(void *vptr);

void o2scl_nucleus_bin_set_hfb27(void *vptr, void *p_v);

void o2scl_skyrme_load_wrapper(void *ptr_sk, void *ptr_model, bool external=false, int verbose=0);

void o2scl_rmf_load_wrapper(void *ptr_rmf, void *ptr_model, bool external=false);

}
