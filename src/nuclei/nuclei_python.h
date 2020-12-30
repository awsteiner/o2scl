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

int o2scl_nucmass_info_parse_elstring(void *vptr, std::string ela, void *ptr_Z, void *ptr_N, void *ptr_A);

int o2scl_nucmass_info_eltoZ(void *vptr, std::string el);

std::string o2scl_nucmass_info_Ztoel(void *vptr, size_t Z);

std::string o2scl_nucmass_info_Ztoname(void *vptr, size_t Z);

std::string o2scl_nucmass_info_tostring(void *vptr, size_t Z, size_t N);

std::string o2scl_nucmass_info_int_to_spinp(void *vptr, int g);

int o2scl_nucmass_info_spinp_to_int(void *vptr, std::string s);

void *o2scl_create_nucmass();

void o2scl_free_nucmass(void *vp);

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

void *o2scl_create_nucmass_table();

void o2scl_free_nucmass_table(void *vp);

size_t o2scl_nucmass_table_get_n(void *vp);

void o2scl_nucmass_table_set_n(void *vp, size_t v);

void o2scl_nucmass_table_get_reference(void *vp, void *p_v);

void o2scl_nucmass_table_set_reference(void *vp, void *p_v);

bool o2scl_nucmass_table_is_loaded(void *vptr, );

size_t o2scl_nucmass_table_get_nentries(void *vptr, );

void *o2scl_create_nucmass_fit_base();

void o2scl_free_nucmass_fit_base(void *vp);

size_t o2scl_nucmass_fit_base_get_nfit(void *vp);

void o2scl_nucmass_fit_base_set_nfit(void *vp, size_t v);

int o2scl_nucmass_fit_base_fit_fun(void *vptr, void *ptr_x);

int o2scl_nucmass_fit_base_guess_fun(void *vptr, void *ptr_x);

}
