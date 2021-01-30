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

extern "C" {

void *o2scl_create_eos_base();

void o2scl_free_eos_base(void *vp);

void o2scl_eos_base_get_def_thermo(void *vp, void *p_v);

void o2scl_eos_base_set_def_thermo(void *vp, void *p_v);

double o2scl_eos_had_base_get_eoa(void *vp);

void o2scl_eos_had_base_set_eoa(void *vp, double v);

double o2scl_eos_had_base_get_msom(void *vp);

void o2scl_eos_had_base_set_msom(void *vp, double v);

double o2scl_eos_had_base_get_comp(void *vp);

void o2scl_eos_had_base_set_comp(void *vp, double v);

double o2scl_eos_had_base_get_n0(void *vp);

void o2scl_eos_had_base_set_n0(void *vp, double v);

double o2scl_eos_had_base_get_esym(void *vp);

void o2scl_eos_had_base_set_esym(void *vp, double v);

double o2scl_eos_had_base_get_kprime(void *vp);

void o2scl_eos_had_base_set_kprime(void *vp, double v);

bool o2scl_eos_had_base_get_err_nonconv(void *vp);

void o2scl_eos_had_base_set_err_nonconv(void *vp, bool v);

void o2scl_eos_had_base_get_def_neutron(void *vp, void *p_v);

void o2scl_eos_had_base_set_def_neutron(void *vp, void *p_v);

void o2scl_eos_had_base_get_def_proton(void *vp, void *p_v);

void o2scl_eos_had_base_set_def_proton(void *vp, void *p_v);

int o2scl_eos_had_base_calc_e(void *vptr, void *ptr_n, void *ptr_p, void *ptr_th);

int o2scl_eos_had_base_calc_p(void *vptr, void *ptr_n, void *ptr_p, void *ptr_th);

double o2scl_eos_had_base_fcomp(void *vptr, double nb, double delta);

double o2scl_eos_had_base_fcomp_err(void *vptr, double nb, double delta, void *ptr_unc);

double o2scl_eos_had_base_feoa(void *vptr, double nb, double delta);

double o2scl_eos_had_base_fesym(void *vptr, double nb, double delta);

double o2scl_eos_had_base_fesym_err(void *vptr, double nb, double delta, void *ptr_unc);

double o2scl_eos_had_base_fesym_slope(void *vptr, double nb, double delta);

double o2scl_eos_had_base_fesym_curve(void *vptr, double nb, double delta);

double o2scl_eos_had_base_fesym_skew(void *vptr, double nb, double delta);

double o2scl_eos_had_base_fesym_diff(void *vptr, double nb);

double o2scl_eos_had_base_feta(void *vptr, double nb);

double o2scl_eos_had_base_feta_prime(void *vptr, double nb);

double o2scl_eos_had_base_fkprime(void *vptr, double nb, double delta);

double o2scl_eos_had_base_fmsom(void *vptr, double nb, double delta);

double o2scl_eos_had_base_f_effm_neut(void *vptr, double nb, double delta);

double o2scl_eos_had_base_f_effm_prot(void *vptr, double nb, double delta);

double o2scl_eos_had_base_f_effm_scalar(void *vptr, double nb, double delta);

double o2scl_eos_had_base_f_effm_vector(void *vptr, double nb, double delta);

double o2scl_eos_had_base_fn0(void *vptr, double delta, void *ptr_leoa);

void o2scl_eos_had_base_f_number_suscept(void *vptr, double mun, double mup, void *ptr_dPdnn, void *ptr_dPdnp, void *ptr_dPdpp);

void o2scl_eos_had_base_f_inv_number_suscept(void *vptr, double mun, double mup, void *ptr_dednn, void *ptr_dednp, void *ptr_dedpp);

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

void *o2scl_create_eos_had_skyrme();

void o2scl_free_eos_had_skyrme(void *vp);

double o2scl_eos_had_skyrme_get_t0(void *vp);

void o2scl_eos_had_skyrme_set_t0(void *vp, double v);

double o2scl_eos_had_skyrme_get_t1(void *vp);

void o2scl_eos_had_skyrme_set_t1(void *vp, double v);

double o2scl_eos_had_skyrme_get_t2(void *vp);

void o2scl_eos_had_skyrme_set_t2(void *vp, double v);

double o2scl_eos_had_skyrme_get_t3(void *vp);

void o2scl_eos_had_skyrme_set_t3(void *vp, double v);

double o2scl_eos_had_skyrme_get_x0(void *vp);

void o2scl_eos_had_skyrme_set_x0(void *vp, double v);

double o2scl_eos_had_skyrme_get_x1(void *vp);

void o2scl_eos_had_skyrme_set_x1(void *vp, double v);

double o2scl_eos_had_skyrme_get_x2(void *vp);

void o2scl_eos_had_skyrme_set_x2(void *vp, double v);

double o2scl_eos_had_skyrme_get_x3(void *vp);

void o2scl_eos_had_skyrme_set_x3(void *vp, double v);

double o2scl_eos_had_skyrme_get_alpha(void *vp);

void o2scl_eos_had_skyrme_set_alpha(void *vp, double v);

double o2scl_eos_had_skyrme_get_a(void *vp);

void o2scl_eos_had_skyrme_set_a(void *vp, double v);

double o2scl_eos_had_skyrme_get_b(void *vp);

void o2scl_eos_had_skyrme_set_b(void *vp, double v);

double o2scl_eos_had_skyrme_get_W0(void *vp);

void o2scl_eos_had_skyrme_set_W0(void *vp, double v);

double o2scl_eos_had_skyrme_get_b4(void *vp);

void o2scl_eos_had_skyrme_set_b4(void *vp, double v);

double o2scl_eos_had_skyrme_get_b4p(void *vp);

void o2scl_eos_had_skyrme_set_b4p(void *vp, double v);

bool o2scl_eos_had_skyrme_get_parent_method(void *vp);

void o2scl_eos_had_skyrme_set_parent_method(void *vp, bool v);

void o2scl_eos_had_skyrme_get_reference(void *vp, void *p_v);

void o2scl_eos_had_skyrme_set_reference(void *vp, void *p_v);

void o2scl_eos_had_skyrme_get_nrfd(void *vp, void *p_v);

void o2scl_eos_had_skyrme_set_nrfd(void *vp, void *p_v);

void *o2scl_create_eos_had_apr();

void o2scl_free_eos_had_apr(void *vp);

int o2scl_eos_had_apr_get_pion(void *vp);

void o2scl_eos_had_apr_set_pion(void *vp, int v);

bool o2scl_eos_had_apr_get_parent_method(void *vp);

void o2scl_eos_had_apr_set_parent_method(void *vp, bool v);

void *o2scl_create_eos_had_rmf();

void o2scl_free_eos_had_rmf(void *vp);

size_t o2scl_eos_had_rmf_get_calc_e_steps(void *vp);

void o2scl_eos_had_rmf_set_calc_e_steps(void *vp, size_t v);

bool o2scl_eos_had_rmf_get_calc_e_relative(void *vp);

void o2scl_eos_had_rmf_set_calc_e_relative(void *vp, bool v);

bool o2scl_eos_had_rmf_get_zm_mode(void *vp);

void o2scl_eos_had_rmf_set_zm_mode(void *vp, bool v);

int o2scl_eos_had_rmf_get_verbose(void *vp);

void o2scl_eos_had_rmf_set_verbose(void *vp, int v);

bool o2scl_eos_had_rmf_get_err_nonconv(void *vp);

void o2scl_eos_had_rmf_set_err_nonconv(void *vp, bool v);

double o2scl_eos_had_rmf_get_mnuc(void *vp);

void o2scl_eos_had_rmf_set_mnuc(void *vp, double v);

double o2scl_eos_had_rmf_get_ms(void *vp);

void o2scl_eos_had_rmf_set_ms(void *vp, double v);

double o2scl_eos_had_rmf_get_mw(void *vp);

void o2scl_eos_had_rmf_set_mw(void *vp, double v);

double o2scl_eos_had_rmf_get_mr(void *vp);

void o2scl_eos_had_rmf_set_mr(void *vp, double v);

double o2scl_eos_had_rmf_get_cs(void *vp);

void o2scl_eos_had_rmf_set_cs(void *vp, double v);

double o2scl_eos_had_rmf_get_cw(void *vp);

void o2scl_eos_had_rmf_set_cw(void *vp, double v);

double o2scl_eos_had_rmf_get_cr(void *vp);

void o2scl_eos_had_rmf_set_cr(void *vp, double v);

double o2scl_eos_had_rmf_get_b(void *vp);

void o2scl_eos_had_rmf_set_b(void *vp, double v);

double o2scl_eos_had_rmf_get_c(void *vp);

void o2scl_eos_had_rmf_set_c(void *vp, double v);

double o2scl_eos_had_rmf_get_zeta(void *vp);

void o2scl_eos_had_rmf_set_zeta(void *vp, double v);

double o2scl_eos_had_rmf_get_xi(void *vp);

void o2scl_eos_had_rmf_set_xi(void *vp, double v);

double o2scl_eos_had_rmf_get_a1(void *vp);

void o2scl_eos_had_rmf_set_a1(void *vp, double v);

double o2scl_eos_had_rmf_get_a2(void *vp);

void o2scl_eos_had_rmf_set_a2(void *vp, double v);

double o2scl_eos_had_rmf_get_a3(void *vp);

void o2scl_eos_had_rmf_set_a3(void *vp, double v);

double o2scl_eos_had_rmf_get_a4(void *vp);

void o2scl_eos_had_rmf_set_a4(void *vp, double v);

double o2scl_eos_had_rmf_get_a5(void *vp);

void o2scl_eos_had_rmf_set_a5(void *vp, double v);

double o2scl_eos_had_rmf_get_a6(void *vp);

void o2scl_eos_had_rmf_set_a6(void *vp, double v);

double o2scl_eos_had_rmf_get_b1(void *vp);

void o2scl_eos_had_rmf_set_b1(void *vp, double v);

double o2scl_eos_had_rmf_get_b2(void *vp);

void o2scl_eos_had_rmf_set_b2(void *vp, double v);

double o2scl_eos_had_rmf_get_b3(void *vp);

void o2scl_eos_had_rmf_set_b3(void *vp, double v);

void *o2scl_create_eos_quark();

void o2scl_free_eos_quark(void *vp);

void *o2scl_create_eos_quark_bag();

void o2scl_free_eos_quark_bag(void *vp);

double o2scl_eos_quark_bag_get_bag_constant(void *vp);

void o2scl_eos_quark_bag_set_bag_constant(void *vp, double v);

void *o2scl_create_eos_quark_njl();

void o2scl_free_eos_quark_njl(void *vp);

double o2scl_eos_quark_njl_get_B0(void *vp);

void o2scl_eos_quark_njl_set_B0(void *vp, double v);

double o2scl_eos_quark_njl_get_L(void *vp);

void o2scl_eos_quark_njl_set_L(void *vp, double v);

double o2scl_eos_quark_njl_get_G(void *vp);

void o2scl_eos_quark_njl_set_G(void *vp, double v);

double o2scl_eos_quark_njl_get_K(void *vp);

void o2scl_eos_quark_njl_set_K(void *vp, double v);

double o2scl_eos_quark_njl_get_limit(void *vp);

void o2scl_eos_quark_njl_set_limit(void *vp, double v);

bool o2scl_eos_quark_njl_get_fromqq(void *vp);

void o2scl_eos_quark_njl_set_fromqq(void *vp, bool v);

void o2scl_skyrme_load_wrapper(void *ptr_sk, char *model, bool external, int verbose);

void o2scl_rmf_load_wrapper(void *ptr_rmf, char *model, bool external);

}
