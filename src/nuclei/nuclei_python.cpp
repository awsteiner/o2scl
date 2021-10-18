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
using namespace o2scl_hdf;

void *o2scl_create_nucleus() {
  nucleus *ptr=new nucleus;
  return ptr;
}

void o2scl_free_nucleus(void *vptr) {
  nucleus *ptr=(nucleus *)vptr;
  delete ptr;
  return;
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
  return;
}

int o2scl_nucmass_info_parse_elstring(void *vptr, char *ela, int *Z, int *N, int *A) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  int ret=ptr->parse_elstring(ela,*Z,*N,*A);
  return ret;
}

int o2scl_nucmass_info_eltoZ(void *vptr, char *el) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  int ret=ptr->eltoZ(el);
  return ret;
}

void *o2scl_nucmass_info_Ztoel(void *vptr, size_t Z) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->Ztoel(Z);
  return sptr;
}

void *o2scl_nucmass_info_Ztoname(void *vptr, size_t Z) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->Ztoname(Z);
  return sptr;
}

void *o2scl_nucmass_info_tostring(void *vptr, size_t Z, size_t N) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->tostring(Z,N);
  return sptr;
}

void *o2scl_nucmass_info_int_to_spinp(void *vptr, int g) {
  nucmass_info *ptr=(nucmass_info *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->int_to_spinp(g);
  return sptr;
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

void *o2scl_nucmass_table_get_reference(void *vptr) {
  nucmass_table *ptr=(nucmass_table *)vptr;
  std::string *sptr=new std::string;
  *sptr=ptr->reference;
  return sptr;
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

void *o2scl_create_nucmass_semi_empirical() {
  nucmass_semi_empirical *ptr=new nucmass_semi_empirical;
  return ptr;
}

void o2scl_free_nucmass_semi_empirical(void *vptr) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  delete ptr;
  return;
}

double o2scl_nucmass_semi_empirical_get_B(void *vptr) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  return ptr->B;
}

void o2scl_nucmass_semi_empirical_set_B(void *vptr, double v) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  ptr->B=v;
  return;
}

double o2scl_nucmass_semi_empirical_get_Sv(void *vptr) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  return ptr->Sv;
}

void o2scl_nucmass_semi_empirical_set_Sv(void *vptr, double v) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  ptr->Sv=v;
  return;
}

double o2scl_nucmass_semi_empirical_get_Ss(void *vptr) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  return ptr->Ss;
}

void o2scl_nucmass_semi_empirical_set_Ss(void *vptr, double v) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  ptr->Ss=v;
  return;
}

double o2scl_nucmass_semi_empirical_get_Ec(void *vptr) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  return ptr->Ec;
}

void o2scl_nucmass_semi_empirical_set_Ec(void *vptr, double v) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  ptr->Ec=v;
  return;
}

double o2scl_nucmass_semi_empirical_get_Epair(void *vptr) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  return ptr->Epair;
}

void o2scl_nucmass_semi_empirical_set_Epair(void *vptr, double v) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  ptr->Epair=v;
  return;
}

double o2scl_nucmass_semi_empirical_mass_excess(void *vptr, int Z, int N) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  double ret=ptr->mass_excess(Z,N);
  return ret;
}

double o2scl_nucmass_semi_empirical_mass_excess_d(void *vptr, double Z, double N) {
  nucmass_semi_empirical *ptr=(nucmass_semi_empirical *)vptr;
  double ret=ptr->mass_excess_d(Z,N);
  return ret;
}

void *o2scl_create_nucmass_ame() {
  nucmass_ame *ptr=new nucmass_ame;
  return ptr;
}

void o2scl_free_nucmass_ame(void *vptr) {
  nucmass_ame *ptr=(nucmass_ame *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_dz_table() {
  nucmass_dz_table *ptr=new nucmass_dz_table;
  return ptr;
}

void o2scl_free_nucmass_dz_table(void *vptr) {
  nucmass_dz_table *ptr=(nucmass_dz_table *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_dz_fit() {
  nucmass_dz_fit *ptr=new nucmass_dz_fit;
  return ptr;
}

void o2scl_free_nucmass_dz_fit(void *vptr) {
  nucmass_dz_fit *ptr=(nucmass_dz_fit *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_dz_fit_33() {
  nucmass_dz_fit_33 *ptr=new nucmass_dz_fit_33;
  return ptr;
}

void o2scl_free_nucmass_dz_fit_33(void *vptr) {
  nucmass_dz_fit_33 *ptr=(nucmass_dz_fit_33 *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_frdm() {
  nucmass_frdm *ptr=new nucmass_frdm;
  return ptr;
}

void o2scl_free_nucmass_frdm(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  delete ptr;
  return;
}

double o2scl_nucmass_frdm_get_a1(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->a1;
}

void o2scl_nucmass_frdm_set_a1(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->a1=v;
  return;
}

double o2scl_nucmass_frdm_get_J(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->J;
}

void o2scl_nucmass_frdm_set_J(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->J=v;
  return;
}

double o2scl_nucmass_frdm_get_K(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->K;
}

void o2scl_nucmass_frdm_set_K(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->K=v;
  return;
}

double o2scl_nucmass_frdm_get_a2(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->a2;
}

void o2scl_nucmass_frdm_set_a2(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->a2=v;
  return;
}

double o2scl_nucmass_frdm_get_Q(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->Q;
}

void o2scl_nucmass_frdm_set_Q(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->Q=v;
  return;
}

double o2scl_nucmass_frdm_get_a3(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->a3;
}

void o2scl_nucmass_frdm_set_a3(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->a3=v;
  return;
}

double o2scl_nucmass_frdm_get_ca(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->ca;
}

void o2scl_nucmass_frdm_set_ca(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->ca=v;
  return;
}

double o2scl_nucmass_frdm_get_W(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->W;
}

void o2scl_nucmass_frdm_set_W(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->W=v;
  return;
}

double o2scl_nucmass_frdm_get_ael(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->ael;
}

void o2scl_nucmass_frdm_set_ael(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->ael=v;
  return;
}

double o2scl_nucmass_frdm_get_rp(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->rp;
}

void o2scl_nucmass_frdm_set_rp(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->rp=v;
  return;
}

double o2scl_nucmass_frdm_get_r0(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->r0;
}

void o2scl_nucmass_frdm_set_r0(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->r0=v;
  return;
}

double o2scl_nucmass_frdm_get_MH(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->MH;
}

void o2scl_nucmass_frdm_set_MH(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->MH=v;
  return;
}

double o2scl_nucmass_frdm_get_Mn(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->Mn;
}

void o2scl_nucmass_frdm_set_Mn(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->Mn=v;
  return;
}

double o2scl_nucmass_frdm_get_e2(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->e2;
}

void o2scl_nucmass_frdm_set_e2(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->e2=v;
  return;
}

double o2scl_nucmass_frdm_get_a(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->a;
}

void o2scl_nucmass_frdm_set_a(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->a=v;
  return;
}

double o2scl_nucmass_frdm_get_aden(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->aden;
}

void o2scl_nucmass_frdm_set_aden(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->aden=v;
  return;
}

double o2scl_nucmass_frdm_get_rmac(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->rmac;
}

void o2scl_nucmass_frdm_set_rmac(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->rmac=v;
  return;
}

double o2scl_nucmass_frdm_get_h(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->h;
}

void o2scl_nucmass_frdm_set_h(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->h=v;
  return;
}

double o2scl_nucmass_frdm_get_L(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->L;
}

void o2scl_nucmass_frdm_set_L(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->L=v;
  return;
}

double o2scl_nucmass_frdm_get_C(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->C;
}

void o2scl_nucmass_frdm_set_C(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->C=v;
  return;
}

double o2scl_nucmass_frdm_get_gamma(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->gamma;
}

void o2scl_nucmass_frdm_set_gamma(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->gamma=v;
  return;
}

double o2scl_nucmass_frdm_get_amu(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->amu;
}

void o2scl_nucmass_frdm_set_amu(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->amu=v;
  return;
}

double o2scl_nucmass_frdm_get_nn(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->nn;
}

void o2scl_nucmass_frdm_set_nn(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->nn=v;
  return;
}

double o2scl_nucmass_frdm_get_np(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->np;
}

void o2scl_nucmass_frdm_set_np(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->np=v;
  return;
}

double o2scl_nucmass_frdm_get_Rn(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->Rn;
}

void o2scl_nucmass_frdm_set_Rn(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->Rn=v;
  return;
}

double o2scl_nucmass_frdm_get_Rp(void *vptr) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  return ptr->Rp;
}

void o2scl_nucmass_frdm_set_Rp(void *vptr, double v) {
  nucmass_frdm *ptr=(nucmass_frdm *)vptr;
  ptr->Rp=v;
  return;
}

void *o2scl_create_nucmass_mnmsk() {
  nucmass_mnmsk *ptr=new nucmass_mnmsk;
  return ptr;
}

void o2scl_free_nucmass_mnmsk(void *vptr) {
  nucmass_mnmsk *ptr=(nucmass_mnmsk *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_mnmsk_exp() {
  nucmass_mnmsk_exp *ptr=new nucmass_mnmsk_exp;
  return ptr;
}

void o2scl_free_nucmass_mnmsk_exp(void *vptr) {
  nucmass_mnmsk_exp *ptr=(nucmass_mnmsk_exp *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_gen() {
  nucmass_gen *ptr=new nucmass_gen;
  return ptr;
}

void o2scl_free_nucmass_gen(void *vptr) {
  nucmass_gen *ptr=(nucmass_gen *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_dglg() {
  nucmass_dglg *ptr=new nucmass_dglg;
  return ptr;
}

void o2scl_free_nucmass_dglg(void *vptr) {
  nucmass_dglg *ptr=(nucmass_dglg *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_hfb() {
  nucmass_hfb *ptr=new nucmass_hfb;
  return ptr;
}

void o2scl_free_nucmass_hfb(void *vptr) {
  nucmass_hfb *ptr=(nucmass_hfb *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_hfb_sp() {
  nucmass_hfb_sp *ptr=new nucmass_hfb_sp;
  return ptr;
}

void o2scl_free_nucmass_hfb_sp(void *vptr) {
  nucmass_hfb_sp *ptr=(nucmass_hfb_sp *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_ktuy() {
  nucmass_ktuy *ptr=new nucmass_ktuy;
  return ptr;
}

void o2scl_free_nucmass_ktuy(void *vptr) {
  nucmass_ktuy *ptr=(nucmass_ktuy *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_sdnp() {
  nucmass_sdnp *ptr=new nucmass_sdnp;
  return ptr;
}

void o2scl_free_nucmass_sdnp(void *vptr) {
  nucmass_sdnp *ptr=(nucmass_sdnp *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_nucmass_wlw() {
  nucmass_wlw *ptr=new nucmass_wlw;
  return ptr;
}

void o2scl_free_nucmass_wlw(void *vptr) {
  nucmass_wlw *ptr=(nucmass_wlw *)vptr;
  delete ptr;
  return;
}

void o2scl_ame_load_wrapper(void *ptr_ame, char *name, bool exp_only) {
  nucmass_ame *ame=(nucmass_ame *)ptr_ame;
  ame_load(*ame,name,exp_only);
  return;
}

void o2scl_ame_load_ext_wrapper(void *ptr_ame, char *file_name, char *table_name, bool exp_only) {
  nucmass_ame *ame=(nucmass_ame *)ptr_ame;
  ame_load_ext(*ame,file_name,table_name,exp_only);
  return;
}

void o2scl_mnmsk_load_wrapper(void *ptr_mnmsk, char *model, char *filename) {
  nucmass_mnmsk *mnmsk=(nucmass_mnmsk *)ptr_mnmsk;
  mnmsk_load(*mnmsk,model,filename);
  return;
}

void o2scl_hfb_load_wrapper(void *ptr_hfb, size_t model, char *filename) {
  nucmass_hfb *hfb=(nucmass_hfb *)ptr_hfb;
  hfb_load(*hfb,model,filename);
  return;
}

void o2scl_hfb_sp_load_wrapper(void *ptr_hfb, size_t model, char *filename) {
  nucmass_hfb_sp *hfb=(nucmass_hfb_sp *)ptr_hfb;
  hfb_sp_load(*hfb,model,filename);
  return;
}

