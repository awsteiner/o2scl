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

#include <o2scl/nuclei_python.h>

using namespace std;
using namespace o2scl;

void *o2scl_create_nucleus() {
  nucleus *ptr=new nucleus;
  return ptr;
}

void o2scl_free_nucleus(void *vptr) {
  nucleus *ptr=(nucleus *)vptr;
  delete ptr;
}

int o2scl_nucleus_get_Z(void *vptr) {
  nucleus *ptr=(nucleus *)vptr;
  return ptr->Z;
}

void o2scl_nucleus_set_Z(void *vptr, int v) {
  nucleus *ptr=(nucleus *)vptr;
  ptr->Z=v;
  return;
}

int o2scl_nucleus_get_N(void *vptr) {
  nucleus *ptr=(nucleus *)vptr;
  return ptr->N;
}

void o2scl_nucleus_set_N(void *vptr, int v) {
  nucleus *ptr=(nucleus *)vptr;
  ptr->N=v;
  return;
}

int o2scl_nucleus_get_A(void *vptr) {
  nucleus *ptr=(nucleus *)vptr;
  return ptr->A;
}

void o2scl_nucleus_set_A(void *vptr, int v) {
  nucleus *ptr=(nucleus *)vptr;
  ptr->A=v;
  return;
}

double o2scl_nucleus_get_mex(void *vptr) {
  nucleus *ptr=(nucleus *)vptr;
  return ptr->mex;
}

void o2scl_nucleus_set_mex(void *vptr, double v) {
  nucleus *ptr=(nucleus *)vptr;
  ptr->mex=v;
  return;
}

double o2scl_nucleus_get_be(void *vptr) {
  nucleus *ptr=(nucleus *)vptr;
  return ptr->be;
}

void o2scl_nucleus_set_be(void *vptr, double v) {
  nucleus *ptr=(nucleus *)vptr;
  ptr->be=v;
  return;
}

void *o2scl_create_nucmass_info() {
  nucmass_info *ptr=new nucmass_info;
  return ptr;
}

void o2scl_free_nucmass_info(void *vptr) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  delete ptr;
}

int o2scl_nucmass_info_parse_elstring(void *vptr, char *ela, void *ptr_Z, void *ptr_N, void *ptr_A) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  int *Z=(int *)ptr_Z;
  int *N=(int *)ptr_N;
  int *A=(int *)ptr_A;
  int ret=ptr->parse_elstring(ela,*Z,*N,*A);
  return ret;
}

int o2scl_nucmass_info_eltoZ(void *vptr, char *el) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  int ret=ptr->eltoZ(el);
  return ret;
}

std::string o2scl_nucmass_info_Ztoel(void *vptr, size_t Z) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  std::string ret=ptr->Ztoel(Z);
  return ret;
}

std::string o2scl_nucmass_info_Ztoname(void *vptr, size_t Z) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  std::string ret=ptr->Ztoname(Z);
  return ret;
}

std::string o2scl_nucmass_info_tostring(void *vptr, size_t Z, size_t N) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  std::string ret=ptr->tostring(Z,N);
  return ret;
}

std::string o2scl_nucmass_info_int_to_spinp(void *vptr, int g) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  std::string ret=ptr->int_to_spinp(g);
  return ret;
}

int o2scl_nucmass_info_spinp_to_int(void *vptr, char *s) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  int ret=ptr->spinp_to_int(s);
  return ret;
}

double o2scl_nucmass_get_m_prot(void *vptr) {
  nucmass *ptr=(nucmass *)vptr;
  return ptr->m_prot;
}

void o2scl_nucmass_set_m_prot(void *vptr, double v) {
  nucmass *ptr=(nucmass *)vptr;
  ptr->m_prot=v;
  return;
}

double o2scl_nucmass_get_m_neut(void *vptr) {
  nucmass *ptr=(nucmass *)vptr;
  return ptr->m_neut;
}

void o2scl_nucmass_set_m_neut(void *vptr, double v) {
  nucmass *ptr=(nucmass *)vptr;
  ptr->m_neut=v;
  return;
}

double o2scl_nucmass_get_m_elec(void *vptr) {
  nucmass *ptr=(nucmass *)vptr;
  return ptr->m_elec;
}

void o2scl_nucmass_set_m_elec(void *vptr, double v) {
  nucmass *ptr=(nucmass *)vptr;
  ptr->m_elec=v;
  return;
}

double o2scl_nucmass_get_m_amu(void *vptr) {
  nucmass *ptr=(nucmass *)vptr;
  return ptr->m_amu;
}

void o2scl_nucmass_set_m_amu(void *vptr, double v) {
  nucmass *ptr=(nucmass *)vptr;
  ptr->m_amu=v;
  return;
}

bool o2scl_nucmass_is_included(void *vptr, int Z, int N) {
  nucmass *ptr=(nucmass *)vptr;
  bool ret=ptr->is_included(Z,N);
  return ret;
}

int o2scl_nucmass_get_nucleus(void *vptr, int Z, int N, void *ptr_n) {
  nucmass *ptr=(nucmass *)vptr;
  nucleus *n=(nucleus *)ptr_n;
  int ret=ptr->get_nucleus(Z,N,*n);
  return ret;
}

double o2scl_nucmass_mass_excess(void *vptr, int Z, int N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->mass_excess(Z,N);
  return ret;
}

double o2scl_nucmass_mass_excess_d(void *vptr, double Z, double N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->mass_excess_d(Z,N);
  return ret;
}

double o2scl_nucmass_electron_binding(void *vptr, double Z) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->electron_binding(Z);
  return ret;
}

double o2scl_nucmass_binding_energy(void *vptr, int Z, int N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->binding_energy(Z,N);
  return ret;
}

double o2scl_nucmass_binding_energy_d(void *vptr, double Z, double N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->binding_energy_d(Z,N);
  return ret;
}

double o2scl_nucmass_total_mass(void *vptr, int Z, int N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->total_mass(Z,N);
  return ret;
}

double o2scl_nucmass_total_mass_d(void *vptr, double Z, double N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->total_mass_d(Z,N);
  return ret;
}

double o2scl_nucmass_neutron_sep(void *vptr, int Z, int N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->neutron_sep(Z,N);
  return ret;
}

double o2scl_nucmass_two_neutron_sep(void *vptr, int Z, int N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->two_neutron_sep(Z,N);
  return ret;
}

double o2scl_nucmass_proton_sep(void *vptr, int Z, int N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->proton_sep(Z,N);
  return ret;
}

double o2scl_nucmass_two_proton_sep(void *vptr, int Z, int N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->two_proton_sep(Z,N);
  return ret;
}

double o2scl_nucmass_atomic_mass(void *vptr, int Z, int N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->atomic_mass(Z,N);
  return ret;
}

double o2scl_nucmass_atomic_mass_d(void *vptr, double Z, double N) {
  nucmass *ptr=(nucmass *)vptr;
  double ret=ptr->atomic_mass_d(Z,N);
  return ret;
}

size_t o2scl_nucmass_table_get_n(void *vptr) {
  nucmass_table *ptr=(nucmass_table *)vptr;
  return ptr->n;
}

void o2scl_nucmass_table_set_n(void *vptr, size_t v) {
  nucmass_table *ptr=(nucmass_table *)vptr;
  ptr->n=v;
  return;
}

void o2scl_nucmass_table_get_reference(void *vptr, void *p_v) {
  nucmass_table *ptr=(nucmass_table *)vptr;
  std::string *p_t=(std::string *)p_v;
  *(p_t)=ptr->reference;
  return;
}

void o2scl_nucmass_table_set_reference(void *vptr, void *p_v) {
  nucmass_table *ptr=(nucmass_table *)vptr;
  std::string *p_t=(std::string *)p_v;
  ptr->reference=*(p_t);
  return;
}

bool o2scl_nucmass_table_is_loaded(void *vptr) {
  nucmass_table *ptr=(nucmass_table *)vptr;
  bool ret=ptr->is_loaded();
  return ret;
}

size_t o2scl_nucmass_table_get_nentries(void *vptr) {
  nucmass_table *ptr=(nucmass_table *)vptr;
  size_t ret=ptr->get_nentries();
  return ret;
}

size_t o2scl_nucmass_fit_base_get_nfit(void *vptr) {
  nucmass_fit_base *ptr=(nucmass_fit_base *)vptr;
  return ptr->nfit;
}

void o2scl_nucmass_fit_base_set_nfit(void *vptr, size_t v) {
  nucmass_fit_base *ptr=(nucmass_fit_base *)vptr;
  ptr->nfit=v;
  return;
}

