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

#include <o2scl/nucleus.h>
#include <o2scl/nucmass.h>
#include <o2scl/nucmass_ame.h>
#include <o2scl/nucmass_dz.h>
#include <o2scl/nucmass_frdm.h>
#include <o2scl/nucmass_dglg.h>
#include <o2scl/nucmass_sdnp.h>
#include <o2scl/nucmass_wlw.h>
#include <o2scl/nucmass_hfb.h>
#include <o2scl/nucmass_ktuy.h>
#include <o2scl/nucmass_fit.h>
#include <o2scl/nucmass_gen.h>
#include <o2scl/hdf_nucmass_io.h>

extern "C" {

void *o2scl_create_nucleus();

void o2scl_free_nucleus(void *vp);

int o2scl_nucleus_get_Z(void *vp);

void o2scl_nucleus_set_Z(void *vp, int v);

int o2scl_nucleus_get_N(void *vp);

void o2scl_nucleus_set_N(void *vp, int v);

int o2scl_nucleus_get_A(void *vp);

void o2scl_nucleus_set_A(void *vp, int v);

double o2scl_nucleus_get_mex(void *vp);

void o2scl_nucleus_set_mex(void *vp, double v);

double o2scl_nucleus_get_be(void *vp);

void o2scl_nucleus_set_be(void *vp, double v);

void *o2scl_create_nucmass_info();

void o2scl_free_nucmass_info(void *vp);

int o2scl_nucmass_info_parse_elstring(void *vptr, char *ela, void *ptr_Z, void *ptr_N, void *ptr_A);

int o2scl_nucmass_info_eltoZ(void *vptr, char *el);

const char *o2scl_nucmass_info_Ztoel(void *vptr, size_t Z);

const char *o2scl_nucmass_info_Ztoname(void *vptr, size_t Z);

const char *o2scl_nucmass_info_tostring(void *vptr, size_t Z, size_t N);

const char *o2scl_nucmass_info_int_to_spinp(void *vptr, int g);

int o2scl_nucmass_info_spinp_to_int(void *vptr, char *s);

double o2scl_nucmass_get_m_prot(void *vp);

void o2scl_nucmass_set_m_prot(void *vp, double v);

double o2scl_nucmass_get_m_neut(void *vp);

void o2scl_nucmass_set_m_neut(void *vp, double v);

double o2scl_nucmass_get_m_elec(void *vp);

void o2scl_nucmass_set_m_elec(void *vp, double v);

double o2scl_nucmass_get_m_amu(void *vp);

void o2scl_nucmass_set_m_amu(void *vp, double v);

bool o2scl_nucmass_is_included(void *vptr, int Z, int N);

int o2scl_nucmass_get_nucleus(void *vptr, int Z, int N, void *ptr_n);

double o2scl_nucmass_mass_excess(void *vptr, int Z, int N);

double o2scl_nucmass_mass_excess_d(void *vptr, double Z, double N);

double o2scl_nucmass_electron_binding(void *vptr, double Z);

double o2scl_nucmass_binding_energy(void *vptr, int Z, int N);

double o2scl_nucmass_binding_energy_d(void *vptr, double Z, double N);

double o2scl_nucmass_total_mass(void *vptr, int Z, int N);

double o2scl_nucmass_total_mass_d(void *vptr, double Z, double N);

double o2scl_nucmass_neutron_sep(void *vptr, int Z, int N);

double o2scl_nucmass_two_neutron_sep(void *vptr, int Z, int N);

double o2scl_nucmass_proton_sep(void *vptr, int Z, int N);

double o2scl_nucmass_two_proton_sep(void *vptr, int Z, int N);

double o2scl_nucmass_atomic_mass(void *vptr, int Z, int N);

double o2scl_nucmass_atomic_mass_d(void *vptr, double Z, double N);

size_t o2scl_nucmass_table_get_n(void *vp);

void o2scl_nucmass_table_set_n(void *vp, size_t v);

void o2scl_nucmass_table_get_reference(void *vp, void *p_v);

void o2scl_nucmass_table_set_reference(void *vp, void *p_v);

bool o2scl_nucmass_table_is_loaded(void *vptr);

size_t o2scl_nucmass_table_get_nentries(void *vptr);

size_t o2scl_nucmass_fit_base_get_nfit(void *vp);

void o2scl_nucmass_fit_base_set_nfit(void *vp, size_t v);

void *o2scl_create_nucmass_semi_empirical();

void o2scl_free_nucmass_semi_empirical(void *vp);

double o2scl_nucmass_semi_empirical_get_B(void *vp);

void o2scl_nucmass_semi_empirical_set_B(void *vp, double v);

double o2scl_nucmass_semi_empirical_get_Sv(void *vp);

void o2scl_nucmass_semi_empirical_set_Sv(void *vp, double v);

double o2scl_nucmass_semi_empirical_get_Ss(void *vp);

void o2scl_nucmass_semi_empirical_set_Ss(void *vp, double v);

double o2scl_nucmass_semi_empirical_get_Ec(void *vp);

void o2scl_nucmass_semi_empirical_set_Ec(void *vp, double v);

double o2scl_nucmass_semi_empirical_get_Epair(void *vp);

void o2scl_nucmass_semi_empirical_set_Epair(void *vp, double v);

double o2scl_nucmass_semi_empirical_mass_excess(void *vptr, int Z, int N);

double o2scl_nucmass_semi_empirical_mass_excess_d(void *vptr, double Z, double N);

void *o2scl_create_nucmass_ame();

void o2scl_free_nucmass_ame(void *vp);

void *o2scl_create_nucmass_dz_table();

void o2scl_free_nucmass_dz_table(void *vp);

void *o2scl_create_nucmass_dz_fit();

void o2scl_free_nucmass_dz_fit(void *vp);

void *o2scl_create_nucmass_dz_fit_33();

void o2scl_free_nucmass_dz_fit_33(void *vp);

void *o2scl_create_nucmass_frdm();

void o2scl_free_nucmass_frdm(void *vp);

double o2scl_nucmass_frdm_get_a1(void *vp);

void o2scl_nucmass_frdm_set_a1(void *vp, double v);

double o2scl_nucmass_frdm_get_J(void *vp);

void o2scl_nucmass_frdm_set_J(void *vp, double v);

double o2scl_nucmass_frdm_get_K(void *vp);

void o2scl_nucmass_frdm_set_K(void *vp, double v);

double o2scl_nucmass_frdm_get_a2(void *vp);

void o2scl_nucmass_frdm_set_a2(void *vp, double v);

double o2scl_nucmass_frdm_get_Q(void *vp);

void o2scl_nucmass_frdm_set_Q(void *vp, double v);

double o2scl_nucmass_frdm_get_a3(void *vp);

void o2scl_nucmass_frdm_set_a3(void *vp, double v);

double o2scl_nucmass_frdm_get_ca(void *vp);

void o2scl_nucmass_frdm_set_ca(void *vp, double v);

double o2scl_nucmass_frdm_get_W(void *vp);

void o2scl_nucmass_frdm_set_W(void *vp, double v);

double o2scl_nucmass_frdm_get_ael(void *vp);

void o2scl_nucmass_frdm_set_ael(void *vp, double v);

double o2scl_nucmass_frdm_get_rp(void *vp);

void o2scl_nucmass_frdm_set_rp(void *vp, double v);

double o2scl_nucmass_frdm_get_r0(void *vp);

void o2scl_nucmass_frdm_set_r0(void *vp, double v);

double o2scl_nucmass_frdm_get_MH(void *vp);

void o2scl_nucmass_frdm_set_MH(void *vp, double v);

double o2scl_nucmass_frdm_get_Mn(void *vp);

void o2scl_nucmass_frdm_set_Mn(void *vp, double v);

double o2scl_nucmass_frdm_get_e2(void *vp);

void o2scl_nucmass_frdm_set_e2(void *vp, double v);

double o2scl_nucmass_frdm_get_a(void *vp);

void o2scl_nucmass_frdm_set_a(void *vp, double v);

double o2scl_nucmass_frdm_get_aden(void *vp);

void o2scl_nucmass_frdm_set_aden(void *vp, double v);

double o2scl_nucmass_frdm_get_rmac(void *vp);

void o2scl_nucmass_frdm_set_rmac(void *vp, double v);

double o2scl_nucmass_frdm_get_h(void *vp);

void o2scl_nucmass_frdm_set_h(void *vp, double v);

double o2scl_nucmass_frdm_get_L(void *vp);

void o2scl_nucmass_frdm_set_L(void *vp, double v);

double o2scl_nucmass_frdm_get_C(void *vp);

void o2scl_nucmass_frdm_set_C(void *vp, double v);

double o2scl_nucmass_frdm_get_gamma(void *vp);

void o2scl_nucmass_frdm_set_gamma(void *vp, double v);

double o2scl_nucmass_frdm_get_amu(void *vp);

void o2scl_nucmass_frdm_set_amu(void *vp, double v);

double o2scl_nucmass_frdm_get_nn(void *vp);

void o2scl_nucmass_frdm_set_nn(void *vp, double v);

double o2scl_nucmass_frdm_get_np(void *vp);

void o2scl_nucmass_frdm_set_np(void *vp, double v);

double o2scl_nucmass_frdm_get_Rn(void *vp);

void o2scl_nucmass_frdm_set_Rn(void *vp, double v);

double o2scl_nucmass_frdm_get_Rp(void *vp);

void o2scl_nucmass_frdm_set_Rp(void *vp, double v);

void *o2scl_create_nucmass_mnmsk();

void o2scl_free_nucmass_mnmsk(void *vp);

void *o2scl_create_nucmass_mnmsk_exp();

void o2scl_free_nucmass_mnmsk_exp(void *vp);

void *o2scl_create_nucmass_gen();

void o2scl_free_nucmass_gen(void *vp);

void *o2scl_create_nucmass_dglg();

void o2scl_free_nucmass_dglg(void *vp);

void *o2scl_create_nucmass_hfb();

void o2scl_free_nucmass_hfb(void *vp);

void *o2scl_create_nucmass_hfb_sp();

void o2scl_free_nucmass_hfb_sp(void *vp);

void *o2scl_create_nucmass_ktuy();

void o2scl_free_nucmass_ktuy(void *vp);

void *o2scl_create_nucmass_sdnp();

void o2scl_free_nucmass_sdnp(void *vp);

void *o2scl_create_nucmass_wlw();

void o2scl_free_nucmass_wlw(void *vp);

void o2scl_ame_load_wrapper(void *ptr_ame, char *name, bool exp_only);

void o2scl_ame_load_ext_wrapper(void *ptr_ame, char *file_name, char *table_name, bool exp_only);

void o2scl_mnmsk_load_wrapper(void *ptr_mnmsk, char *model, char *filename);

void o2scl_hfb_load_wrapper(void *ptr_hfb, size_t model, char *filename);

void o2scl_hfb_sp_load_wrapper(void *ptr_hfb, size_t model, char *filename);

}
