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

#include <o2scl/part.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/boson_rel.h>
#include <o2scl/classical.h>
#include <o2scl/classical_deriv.h>
#include <o2scl/fermion_mag_zerot.h>
#include <o2scl/quark.h>

extern "C" {

void *o2scl_create_thermo();

void o2scl_free_thermo(void *vptr);

double o2scl_thermo_get_ed(void *vptr);

void o2scl_thermo_set_ed(void *vptr, double v);

double o2scl_thermo_get_pr(void *vptr);

void o2scl_thermo_set_pr(void *vptr, double v);

double o2scl_thermo_get_en(void *vptr);

void o2scl_thermo_set_en(void *vptr, double v);

void *o2scl_create_part();

void o2scl_free_part(void *vptr);

double o2scl_part_get_g(void *vptr);

void o2scl_part_set_g(void *vptr, double v);

double o2scl_part_get_m(void *vptr);

void o2scl_part_set_m(void *vptr, double v);

double o2scl_part_get_ms(void *vptr);

void o2scl_part_set_ms(void *vptr, double v);

double o2scl_part_get_mu(void *vptr);

void o2scl_part_set_mu(void *vptr, double v);

double o2scl_part_get_nu(void *vptr);

void o2scl_part_set_nu(void *vptr, double v);

double o2scl_part_get_n(void *vptr);

void o2scl_part_set_n(void *vptr, double v);

double o2scl_part_get_ed(void *vptr);

void o2scl_part_set_ed(void *vptr, double v);

double o2scl_part_get_pr(void *vptr);

void o2scl_part_set_pr(void *vptr, double v);

double o2scl_part_get_en(void *vptr);

void o2scl_part_set_en(void *vptr, double v);

bool o2scl_part_get_inc_rest_mass(void *vptr);

void o2scl_part_set_inc_rest_mass(void *vptr, bool v);

bool o2scl_part_get_non_interacting(void *vptr);

void o2scl_part_set_non_interacting(void *vptr, bool v);

void o2scl_part_init(void *vptr, double mass, double dof);

void o2scl_part_anti(void *vptr, void *ptr_ax);

void *o2scl_create_fermion();

void o2scl_free_fermion(void *vptr);

double o2scl_fermion_get_kf(void *vptr);

void o2scl_fermion_set_kf(void *vptr, double v);

double o2scl_fermion_get_del(void *vptr);

void o2scl_fermion_set_del(void *vptr, double v);

void *o2scl_create_quark();

void o2scl_free_quark(void *vptr);

double o2scl_quark_get_B(void *vptr);

void o2scl_quark_set_B(void *vptr, double v);

double o2scl_quark_get_qq(void *vptr);

void o2scl_quark_set_qq(void *vptr, double v);

void *o2scl_create_fermion_zerot();

void o2scl_free_fermion_zerot(void *vptr);

void o2scl_fermion_zerot_kf_from_density(void *vptr, void *ptr_f);

void o2scl_fermion_zerot_energy_density_zerot(void *vptr, void *ptr_f);

void o2scl_fermion_zerot_pressure_zerot(void *vptr, void *ptr_f);

void o2scl_fermion_zerot_calc_mu_zerot(void *vptr, void *ptr_f);

void o2scl_fermion_zerot_calc_density_zerot(void *vptr, void *ptr_f);

bool o2scl_fermion_thermo_calc_mu_deg(void *vptr, void *ptr_f, double T, double prec);

bool o2scl_fermion_thermo_calc_mu_ndeg(void *vptr, void *ptr_f, double T, double prec, bool inc_antip);

void o2scl_fermion_thermo_massless_calc_mu(void *vptr, void *ptr_f, double T);

void o2scl_fermion_thermo_massless_pair_mu(void *vptr, void *ptr_f, double T);

void o2scl_fermion_thermo_massless_calc_density(void *vptr, void *ptr_f, double T);

void o2scl_fermion_thermo_massless_pair_density(void *vptr, void *ptr_f, double T);

void *o2scl_create_fermion_rel();

void o2scl_free_fermion_rel(void *vptr);

bool o2scl_fermion_rel_get_err_nonconv(void *vptr);

void o2scl_fermion_rel_set_err_nonconv(void *vptr, bool v);

double o2scl_fermion_rel_get_min_psi(void *vptr);

void o2scl_fermion_rel_set_min_psi(void *vptr, double v);

double o2scl_fermion_rel_get_deg_limit(void *vptr);

void o2scl_fermion_rel_set_deg_limit(void *vptr, double v);

double o2scl_fermion_rel_get_upper_limit_fac(void *vptr);

void o2scl_fermion_rel_set_upper_limit_fac(void *vptr, double v);

int o2scl_fermion_rel_get_verbose(void *vptr);

void o2scl_fermion_rel_set_verbose(void *vptr, int v);

bool o2scl_fermion_rel_get_use_expansions(void *vptr);

void o2scl_fermion_rel_set_use_expansions(void *vptr, bool v);

double o2scl_fermion_rel_get_tol_expan(void *vptr);

void o2scl_fermion_rel_set_tol_expan(void *vptr, double v);

bool o2scl_fermion_rel_get_verify_ti(void *vptr);

void o2scl_fermion_rel_set_verify_ti(void *vptr, bool v);

double o2scl_fermion_rel_get_therm_ident(void *vptr);

void o2scl_fermion_rel_set_therm_ident(void *vptr, double v);

void *o2scl_fermion_rel_get_unc(void *vptr);

void o2scl_fermion_rel_set_unc(void *vptr, void *p_v);

int o2scl_fermion_rel_nu_from_n(void *vptr, void *ptr_f, double T);

int o2scl_fermion_rel_calc_density(void *vptr, void *ptr_f, double T);

int o2scl_fermion_rel_pair_density(void *vptr, void *ptr_f, double T);

void o2scl_fermion_rel_calc_mu(void *vptr, void *ptr_f, double T);

void o2scl_fermion_rel_pair_mu(void *vptr, void *ptr_f, double T);

void *o2scl_create_fermion_nonrel();

void o2scl_free_fermion_nonrel(void *vptr);

int o2scl_fermion_nonrel_calc_density(void *vptr, void *ptr_f, double T);

void o2scl_fermion_nonrel_calc_mu(void *vptr, void *ptr_f, double T);

void o2scl_fermion_nonrel_nu_from_n(void *vptr, void *ptr_f, double T);

void *o2scl_create_boson();

void o2scl_free_boson(void *vptr);

double o2scl_boson_get_co(void *vptr);

void o2scl_boson_set_co(void *vptr, double v);

void *o2scl_create_boson_rel();

void o2scl_free_boson_rel(void *vptr);

void o2scl_boson_rel_calc_density(void *vptr, void *ptr_b, double T);

void o2scl_boson_rel_calc_mu(void *vptr, void *ptr_b, double T);

void o2scl_boson_rel_nu_from_n(void *vptr, void *ptr_b, double T);

void o2scl_boson_rel_pair_density(void *vptr, void *ptr_b, double T);

void o2scl_boson_rel_pair_mu(void *vptr, void *ptr_b, double T);

void *o2scl_create_classical_thermo();

void o2scl_free_classical_thermo(void *vptr);

void o2scl_classical_thermo_calc_density(void *vptr, void *ptr_p, double T);

void o2scl_classical_thermo_calc_mu(void *vptr, void *ptr_p, double T);

void *o2scl_create_thermo_np_deriv_press();

void o2scl_free_thermo_np_deriv_press(void *vptr);

double o2scl_thermo_np_deriv_press_get_dsdT(void *vptr);

void o2scl_thermo_np_deriv_press_set_dsdT(void *vptr, double v);

double o2scl_thermo_np_deriv_press_get_dnndT(void *vptr);

void o2scl_thermo_np_deriv_press_set_dnndT(void *vptr, double v);

double o2scl_thermo_np_deriv_press_get_dnpdT(void *vptr);

void o2scl_thermo_np_deriv_press_set_dnpdT(void *vptr, double v);

double o2scl_thermo_np_deriv_press_get_dnndmun(void *vptr);

void o2scl_thermo_np_deriv_press_set_dnndmun(void *vptr, double v);

double o2scl_thermo_np_deriv_press_get_dndmu_mixed(void *vptr);

void o2scl_thermo_np_deriv_press_set_dndmu_mixed(void *vptr, double v);

double o2scl_thermo_np_deriv_press_get_dnpdmup(void *vptr);

void o2scl_thermo_np_deriv_press_set_dnpdmup(void *vptr, double v);

void *o2scl_create_thermo_np_deriv_helm();

void o2scl_free_thermo_np_deriv_helm(void *vptr);

double o2scl_thermo_np_deriv_helm_get_dsdT(void *vptr);

void o2scl_thermo_np_deriv_helm_set_dsdT(void *vptr, double v);

double o2scl_thermo_np_deriv_helm_get_dmundT(void *vptr);

void o2scl_thermo_np_deriv_helm_set_dmundT(void *vptr, double v);

double o2scl_thermo_np_deriv_helm_get_dmupdT(void *vptr);

void o2scl_thermo_np_deriv_helm_set_dmupdT(void *vptr, double v);

double o2scl_thermo_np_deriv_helm_get_dmundnn(void *vptr);

void o2scl_thermo_np_deriv_helm_set_dmundnn(void *vptr, double v);

double o2scl_thermo_np_deriv_helm_get_dmudn_mixed(void *vptr);

void o2scl_thermo_np_deriv_helm_set_dmudn_mixed(void *vptr, double v);

double o2scl_thermo_np_deriv_helm_get_dmupdnp(void *vptr);

void o2scl_thermo_np_deriv_helm_set_dmupdnp(void *vptr, double v);

void *o2scl_create_part_deriv_press();

void o2scl_free_part_deriv_press(void *vptr);

double o2scl_part_deriv_press_get_dndmu(void *vptr);

void o2scl_part_deriv_press_set_dndmu(void *vptr, double v);

double o2scl_part_deriv_press_get_dndT(void *vptr);

void o2scl_part_deriv_press_set_dndT(void *vptr, double v);

double o2scl_part_deriv_press_get_dsdT(void *vptr);

void o2scl_part_deriv_press_set_dsdT(void *vptr, double v);

void o2scl_part_deriv_press_deriv_f(void *vptr, double *dmudn, double *dmudT, double *dsdT_n);

void *o2scl_create_part_deriv();

void o2scl_free_part_deriv(void *vptr);

void *o2scl_create_fermion_deriv();

void o2scl_free_fermion_deriv(void *vptr);

void *o2scl_create_deriv_thermo_base();

void o2scl_free_deriv_thermo_base(void *vptr);

double o2scl_deriv_thermo_base_heat_cap_ppart_const_vol(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_heat_cap_ppart_const_press(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_compress_adiabatic(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_compress_const_tptr(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_coeff_thermal_exp(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_squared_sound_speed(void *vptr, void *ptr_p, double T);

void *o2scl_create_classical_deriv_thermo();

void o2scl_free_classical_deriv_thermo(void *vptr);

void o2scl_classical_deriv_thermo_calc_density(void *vptr, void *ptr_p, double T);

void o2scl_classical_deriv_thermo_calc_mu(void *vptr, void *ptr_p, double T);

void *o2scl_create_fermion_mag_zerot();

void o2scl_free_fermion_mag_zerot(void *vptr);

int o2scl_fermion_mag_zerot_get_nmax_up(void *vptr);

void o2scl_fermion_mag_zerot_set_nmax_up(void *vptr, int v);

int o2scl_fermion_mag_zerot_get_nmax_dn(void *vptr);

void o2scl_fermion_mag_zerot_set_nmax_dn(void *vptr, int v);

int o2scl_fermion_mag_zerot_get_sum_limit(void *vptr);

void o2scl_fermion_mag_zerot_set_sum_limit(void *vptr, int v);

void o2scl_fermion_mag_zerot_calc_mu_zerot_mag(void *vptr, void *ptr_f, double qB, double kappa);

void o2scl_fermion_mag_zerot_calc_density_zerot_mag(void *vptr, void *ptr_f, double qB, double kappa);

}
