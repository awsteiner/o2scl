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

#include <o2scl/part.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/boson_rel.h>
#include <o2scl/classical.h>
#include <o2scl/classical_deriv.h>
#include <o2scl/fermion_mag_zerot.h>

extern "C" {

void *o2scl_create_thermo();

void o2scl_free_thermo(void *vp);

double o2scl_thermo_get_ed(void *vp);

void o2scl_thermo_set_ed(void *vp, double v);

double o2scl_thermo_get_pr(void *vp);

void o2scl_thermo_set_pr(void *vp, double v);

double o2scl_thermo_get_en(void *vp);

void o2scl_thermo_set_en(void *vp, double v);

void *o2scl_create_part();

void o2scl_free_part(void *vp);

double o2scl_part_get_g(void *vp);

void o2scl_part_set_g(void *vp, double v);

double o2scl_part_get_m(void *vp);

void o2scl_part_set_m(void *vp, double v);

double o2scl_part_get_ms(void *vp);

void o2scl_part_set_ms(void *vp, double v);

double o2scl_part_get_mu(void *vp);

void o2scl_part_set_mu(void *vp, double v);

double o2scl_part_get_nu(void *vp);

void o2scl_part_set_nu(void *vp, double v);

double o2scl_part_get_ed(void *vp);

void o2scl_part_set_ed(void *vp, double v);

double o2scl_part_get_pr(void *vp);

void o2scl_part_set_pr(void *vp, double v);

double o2scl_part_get_en(void *vp);

void o2scl_part_set_en(void *vp, double v);

bool o2scl_part_get_inc_rest_mass(void *vp);

void o2scl_part_set_inc_rest_mass(void *vp, bool v);

bool o2scl_part_get_non_interacting(void *vp);

void o2scl_part_set_non_interacting(void *vp, bool v);

void o2scl_part_init(void *vptr, double mass, double dof);

void o2scl_part_anti(void *vptr, void *ptr_ax);

void *o2scl_create_fermion();

void o2scl_free_fermion(void *vp);

double o2scl_fermion_get_kf(void *vp);

void o2scl_fermion_set_kf(void *vp, double v);

double o2scl_fermion_get_del(void *vp);

void o2scl_fermion_set_del(void *vp, double v);

void *o2scl_create_fermion_zerot();

void o2scl_free_fermion_zerot(void *vp);

void o2scl_fermion_zerot_kf_from_density(void *vptr, void *ptr_f);

void o2scl_fermion_zerot_energy_density_zerot(void *vptr, void *ptr_f);

void o2scl_fermion_zerot_pressure_zerot(void *vptr, void *ptr_f);

void o2scl_fermion_zerot_calc_mu_zerot(void *vptr, void *ptr_f);

void o2scl_fermion_zerot_calc_density_zerot(void *vptr, void *ptr_f);

void *o2scl_create_fermion_rel();

void o2scl_free_fermion_rel(void *vp);

bool o2scl_fermion_rel_get_err_nonconv(void *vp);

void o2scl_fermion_rel_set_err_nonconv(void *vp, bool v);

double o2scl_fermion_rel_get_min_psi(void *vp);

void o2scl_fermion_rel_set_min_psi(void *vp, double v);

double o2scl_fermion_rel_get_deg_limit(void *vp);

void o2scl_fermion_rel_set_deg_limit(void *vp, double v);

double o2scl_fermion_rel_get_exp_limit(void *vp);

void o2scl_fermion_rel_set_exp_limit(void *vp, double v);

double o2scl_fermion_rel_get_upper_limit_fac(void *vp);

void o2scl_fermion_rel_set_upper_limit_fac(void *vp, double v);

double o2scl_fermion_rel_get_deg_entropy_fac(void *vp);

void o2scl_fermion_rel_set_deg_entropy_fac(void *vp, double v);

int o2scl_fermion_rel_get_verbose(void *vp);

void o2scl_fermion_rel_set_verbose(void *vp, int v);

bool o2scl_fermion_rel_get_use_expansions(void *vp);

void o2scl_fermion_rel_set_use_expansions(void *vp, bool v);

double o2scl_fermion_rel_get_tol_expan(void *vp);

void o2scl_fermion_rel_set_tol_expan(void *vp, double v);

bool o2scl_fermion_rel_get_verify_ti(void *vp);

void o2scl_fermion_rel_set_verify_ti(void *vp, bool v);

double o2scl_fermion_rel_get_therm_ident(void *vp);

void o2scl_fermion_rel_set_therm_ident(void *vp, double v);

void *o2scl_fermion_rel_get_unc(void *vp);

void o2scl_fermion_rel_set_unc(void *vp, void *p_v);

bool o2scl_fermion_rel_calc_mu_deg(void *vptr, void *ptr_f, double T, double prec);

bool o2scl_fermion_rel_calc_mu_ndeg(void *vptr, void *ptr_f, double T, double prec, bool inc_antip);

void o2scl_fermion_rel_massless_calc_mu(void *vptr, void *ptr_f, double T);

void o2scl_fermion_rel_massless_pair_mu(void *vptr, void *ptr_f, double T);

void o2scl_fermion_rel_massless_calc_density(void *vptr, void *ptr_f, double T);

void o2scl_fermion_rel_massless_pair_density(void *vptr, void *ptr_f, double T);

int o2scl_fermion_rel_nu_from_n(void *vptr, void *ptr_f, double T);

int o2scl_fermion_rel_calc_density(void *vptr, void *ptr_f, double T);

int o2scl_fermion_rel_pair_density(void *vptr, void *ptr_f, double T);

void o2scl_fermion_rel_calc_mu(void *vptr, void *ptr_f, double T);

void o2scl_fermion_rel_pair_mu(void *vptr, void *ptr_f, double T);

void *o2scl_create_fermion_nonrel();

void o2scl_free_fermion_nonrel(void *vp);

int o2scl_fermion_nonrel_calc_density(void *vptr, void *ptr_f, double T);

void o2scl_fermion_nonrel_calc_mu(void *vptr, void *ptr_f, double T);

void o2scl_fermion_nonrel_nu_from_n(void *vptr, void *ptr_f, double T);

void *o2scl_create_boson();

void o2scl_free_boson(void *vp);

double o2scl_boson_get_co(void *vp);

void o2scl_boson_set_co(void *vp, double v);

void *o2scl_create_boson_rel();

void o2scl_free_boson_rel(void *vp);

void o2scl_boson_rel_calc_density(void *vptr, void *ptr_b, double T);

void o2scl_boson_rel_calc_mu(void *vptr, void *ptr_b, double T);

void o2scl_boson_rel_nu_from_n(void *vptr, void *ptr_b, double T);

void o2scl_boson_rel_pair_density(void *vptr, void *ptr_b, double T);

void o2scl_boson_rel_pair_mu(void *vptr, void *ptr_b, double T);

void *o2scl_create_classical_thermo();

void o2scl_free_classical_thermo(void *vp);

void o2scl_classical_thermo_calc_density(void *vptr, void *ptr_p, double T);

void o2scl_classical_thermo_calc_mu(void *vptr, void *ptr_p, double T);

void *o2scl_create_thermo_np_deriv_press();

void o2scl_free_thermo_np_deriv_press(void *vp);

double o2scl_thermo_np_deriv_press_get_dsdT(void *vp);

void o2scl_thermo_np_deriv_press_set_dsdT(void *vp, double v);

double o2scl_thermo_np_deriv_press_get_dnndT(void *vp);

void o2scl_thermo_np_deriv_press_set_dnndT(void *vp, double v);

double o2scl_thermo_np_deriv_press_get_dnpdT(void *vp);

void o2scl_thermo_np_deriv_press_set_dnpdT(void *vp, double v);

double o2scl_thermo_np_deriv_press_get_dnndmun(void *vp);

void o2scl_thermo_np_deriv_press_set_dnndmun(void *vp, double v);

double o2scl_thermo_np_deriv_press_get_dndmu_mixed(void *vp);

void o2scl_thermo_np_deriv_press_set_dndmu_mixed(void *vp, double v);

double o2scl_thermo_np_deriv_press_get_dnpdmup(void *vp);

void o2scl_thermo_np_deriv_press_set_dnpdmup(void *vp, double v);

void *o2scl_create_thermo_np_deriv_helm();

void o2scl_free_thermo_np_deriv_helm(void *vp);

double o2scl_thermo_np_deriv_helm_get_dsdT(void *vp);

void o2scl_thermo_np_deriv_helm_set_dsdT(void *vp, double v);

double o2scl_thermo_np_deriv_helm_get_dmundT(void *vp);

void o2scl_thermo_np_deriv_helm_set_dmundT(void *vp, double v);

double o2scl_thermo_np_deriv_helm_get_dmupdT(void *vp);

void o2scl_thermo_np_deriv_helm_set_dmupdT(void *vp, double v);

double o2scl_thermo_np_deriv_helm_get_dmundnn(void *vp);

void o2scl_thermo_np_deriv_helm_set_dmundnn(void *vp, double v);

double o2scl_thermo_np_deriv_helm_get_dmudn_mixed(void *vp);

void o2scl_thermo_np_deriv_helm_set_dmudn_mixed(void *vp, double v);

double o2scl_thermo_np_deriv_helm_get_dmupdnp(void *vp);

void o2scl_thermo_np_deriv_helm_set_dmupdnp(void *vp, double v);

void *o2scl_create_part_deriv_press();

void o2scl_free_part_deriv_press(void *vp);

double o2scl_part_deriv_press_get_dndmu(void *vp);

void o2scl_part_deriv_press_set_dndmu(void *vp, double v);

double o2scl_part_deriv_press_get_dndT(void *vp);

void o2scl_part_deriv_press_set_dndT(void *vp, double v);

double o2scl_part_deriv_press_get_dsdT(void *vp);

void o2scl_part_deriv_press_set_dsdT(void *vp, double v);

void o2scl_part_deriv_press_deriv_f(void *vptr, void *ptr_dmudn, void *ptr_dmudT, void *ptr_dsdT_n);

void *o2scl_create_part_deriv();

void o2scl_free_part_deriv(void *vp);

void *o2scl_create_fermion_deriv();

void o2scl_free_fermion_deriv(void *vp);

void *o2scl_create_deriv_thermo_base();

void o2scl_free_deriv_thermo_base(void *vp);

double o2scl_deriv_thermo_base_heat_cap_ppart_const_vol(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_heat_cap_ppart_const_press(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_compress_adiabatic(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_compress_const_tptr(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_coeff_thermal_exp(void *vptr, void *ptr_p, double T);

double o2scl_deriv_thermo_base_squared_sound_speed(void *vptr, void *ptr_p, double T);

void *o2scl_create_fermion_deriv_rel();

void o2scl_free_fermion_deriv_rel(void *vp);

double o2scl_fermion_deriv_rel_get_exp_limit(void *vp);

void o2scl_fermion_deriv_rel_set_exp_limit(void *vp, double v);

double o2scl_fermion_deriv_rel_get_deg_limit(void *vp);

void o2scl_fermion_deriv_rel_set_deg_limit(void *vp, double v);

double o2scl_fermion_deriv_rel_get_upper_limit_fac(void *vp);

void o2scl_fermion_deriv_rel_set_upper_limit_fac(void *vp, double v);

void *o2scl_fermion_deriv_rel_get_unc(void *vp);

void o2scl_fermion_deriv_rel_set_unc(void *vp, void *p_v);

int o2scl_fermion_deriv_rel_get_method(void *vp);

void o2scl_fermion_deriv_rel_set_method(void *vp, int v);

int o2scl_fermion_deriv_rel_get_last_method(void *vp);

void o2scl_fermion_deriv_rel_set_last_method(void *vp, int v);

bool o2scl_fermion_deriv_rel_get_err_nonconv(void *vp);

void o2scl_fermion_deriv_rel_set_err_nonconv(void *vp, bool v);

int o2scl_fermion_deriv_rel_nu_from_n(void *vptr, void *ptr_f, double T);

int o2scl_fermion_deriv_rel_calc_density(void *vptr, void *ptr_f, double T);

int o2scl_fermion_deriv_rel_pair_density(void *vptr, void *ptr_f, double T);

int o2scl_fermion_deriv_rel_calc_mu(void *vptr, void *ptr_f, double T);

int o2scl_fermion_deriv_rel_pair_mu(void *vptr, void *ptr_f, double T);

void *o2scl_create_fermion_deriv_nr();

void o2scl_free_fermion_deriv_nr(void *vp);

double o2scl_fermion_deriv_nr_get_flimit(void *vp);

void o2scl_fermion_deriv_nr_set_flimit(void *vp, double v);

void *o2scl_fermion_deriv_nr_get_unc(void *vp);

void o2scl_fermion_deriv_nr_set_unc(void *vp, void *p_v);

void o2scl_fermion_deriv_nr_calc_density_zerot(void *vptr, void *ptr_f);

void o2scl_fermion_deriv_nr_calc_mu_zerot(void *vptr, void *ptr_f);

int o2scl_fermion_deriv_nr_nu_from_n(void *vptr, void *ptr_f, double T);

int o2scl_fermion_deriv_nr_calc_density(void *vptr, void *ptr_f, double T);

int o2scl_fermion_deriv_nr_calc_mu(void *vptr, void *ptr_f, double T);

void *o2scl_create_classical_deriv_thermo();

void o2scl_free_classical_deriv_thermo(void *vp);

void o2scl_classical_deriv_thermo_calc_density(void *vptr, void *ptr_p, double T);

void o2scl_classical_deriv_thermo_calc_mu(void *vptr, void *ptr_p, double T);

}
