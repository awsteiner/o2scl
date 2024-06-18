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

#include <o2scl/part_python.h>

using namespace std;
using namespace o2scl;

void *o2scl_create_thermo() {
  thermo *ptr=new thermo;
  return ptr;
}

void o2scl_free_thermo(void *vptr) {
  thermo *ptr=(thermo *)vptr;
  delete ptr;
  return;
}

double o2scl_thermo_get_ed(void *vptr) {
  thermo *ptr=(thermo *)vptr;
  return ptr->ed;
}

void o2scl_thermo_set_ed(void *vptr, double v) {
  thermo *ptr=(thermo *)vptr;
  ptr->ed=v;
  return;
}

double o2scl_thermo_get_pr(void *vptr) {
  thermo *ptr=(thermo *)vptr;
  return ptr->pr;
}

void o2scl_thermo_set_pr(void *vptr, double v) {
  thermo *ptr=(thermo *)vptr;
  ptr->pr=v;
  return;
}

double o2scl_thermo_get_en(void *vptr) {
  thermo *ptr=(thermo *)vptr;
  return ptr->en;
}

void o2scl_thermo_set_en(void *vptr, double v) {
  thermo *ptr=(thermo *)vptr;
  ptr->en=v;
  return;
}

void *o2scl_create_part() {
  part *ptr=new part;
  return ptr;
}

void o2scl_free_part(void *vptr) {
  part *ptr=(part *)vptr;
  delete ptr;
  return;
}

double o2scl_part_get_g(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->g;
}

void o2scl_part_set_g(void *vptr, double v) {
  part *ptr=(part *)vptr;
  ptr->g=v;
  return;
}

double o2scl_part_get_m(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->m;
}

void o2scl_part_set_m(void *vptr, double v) {
  part *ptr=(part *)vptr;
  ptr->m=v;
  return;
}

double o2scl_part_get_ms(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->ms;
}

void o2scl_part_set_ms(void *vptr, double v) {
  part *ptr=(part *)vptr;
  ptr->ms=v;
  return;
}

double o2scl_part_get_mu(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->mu;
}

void o2scl_part_set_mu(void *vptr, double v) {
  part *ptr=(part *)vptr;
  ptr->mu=v;
  return;
}

double o2scl_part_get_nu(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->nu;
}

void o2scl_part_set_nu(void *vptr, double v) {
  part *ptr=(part *)vptr;
  ptr->nu=v;
  return;
}

double o2scl_part_get_n(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->n;
}

void o2scl_part_set_n(void *vptr, double v) {
  part *ptr=(part *)vptr;
  ptr->n=v;
  return;
}

double o2scl_part_get_ed(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->ed;
}

void o2scl_part_set_ed(void *vptr, double v) {
  part *ptr=(part *)vptr;
  ptr->ed=v;
  return;
}

double o2scl_part_get_pr(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->pr;
}

void o2scl_part_set_pr(void *vptr, double v) {
  part *ptr=(part *)vptr;
  ptr->pr=v;
  return;
}

double o2scl_part_get_en(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->en;
}

void o2scl_part_set_en(void *vptr, double v) {
  part *ptr=(part *)vptr;
  ptr->en=v;
  return;
}

bool o2scl_part_get_inc_rest_mass(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->inc_rest_mass;
}

void o2scl_part_set_inc_rest_mass(void *vptr, bool v) {
  part *ptr=(part *)vptr;
  ptr->inc_rest_mass=v;
  return;
}

bool o2scl_part_get_non_interacting(void *vptr) {
  part *ptr=(part *)vptr;
  return ptr->non_interacting;
}

void o2scl_part_set_non_interacting(void *vptr, bool v) {
  part *ptr=(part *)vptr;
  ptr->non_interacting=v;
  return;
}

void o2scl_part_init(void *vptr, double mass, double dof) {
  part *ptr=(part *)vptr;
  ptr->init(mass,dof);
  return;
}

void o2scl_part_anti(void *vptr, void *ptr_ax) {
  part *ptr=(part *)vptr;
  part *ax=(part *)ptr_ax;
  ptr->anti(*ax);
  return;
}

void *o2scl_create_fermion() {
  fermion *ptr=new fermion;
  return ptr;
}

void o2scl_free_fermion(void *vptr) {
  fermion *ptr=(fermion *)vptr;
  delete ptr;
  return;
}

double o2scl_fermion_get_kf(void *vptr) {
  fermion *ptr=(fermion *)vptr;
  return ptr->kf;
}

void o2scl_fermion_set_kf(void *vptr, double v) {
  fermion *ptr=(fermion *)vptr;
  ptr->kf=v;
  return;
}

double o2scl_fermion_get_del(void *vptr) {
  fermion *ptr=(fermion *)vptr;
  return ptr->del;
}

void o2scl_fermion_set_del(void *vptr, double v) {
  fermion *ptr=(fermion *)vptr;
  ptr->del=v;
  return;
}

void *o2scl_create_quark() {
  quark *ptr=new quark;
  return ptr;
}

void o2scl_free_quark(void *vptr) {
  quark *ptr=(quark *)vptr;
  delete ptr;
  return;
}

double o2scl_quark_get_B(void *vptr) {
  quark *ptr=(quark *)vptr;
  return ptr->B;
}

void o2scl_quark_set_B(void *vptr, double v) {
  quark *ptr=(quark *)vptr;
  ptr->B=v;
  return;
}

double o2scl_quark_get_qq(void *vptr) {
  quark *ptr=(quark *)vptr;
  return ptr->qq;
}

void o2scl_quark_set_qq(void *vptr, double v) {
  quark *ptr=(quark *)vptr;
  ptr->qq=v;
  return;
}

void *o2scl_create_fermion_zerot() {
  fermion_zerot *ptr=new fermion_zerot;
  return ptr;
}

void o2scl_free_fermion_zerot(void *vptr) {
  fermion_zerot *ptr=(fermion_zerot *)vptr;
  delete ptr;
  return;
}

void o2scl_fermion_zerot_kf_from_density(void *vptr, void *ptr_f) {
  fermion_zerot *ptr=(fermion_zerot *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->kf_from_density(*f);
  return;
}

void o2scl_fermion_zerot_energy_density_zerot(void *vptr, void *ptr_f) {
  fermion_zerot *ptr=(fermion_zerot *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->energy_density_zerot(*f);
  return;
}

void o2scl_fermion_zerot_pressure_zerot(void *vptr, void *ptr_f) {
  fermion_zerot *ptr=(fermion_zerot *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->pressure_zerot(*f);
  return;
}

void o2scl_fermion_zerot_calc_mu_zerot(void *vptr, void *ptr_f) {
  fermion_zerot *ptr=(fermion_zerot *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->calc_mu_zerot(*f);
  return;
}

void o2scl_fermion_zerot_calc_density_zerot(void *vptr, void *ptr_f) {
  fermion_zerot *ptr=(fermion_zerot *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->calc_density_zerot(*f);
  return;
}

bool o2scl_fermion_thermo_calc_mu_deg(void *vptr, void *ptr_f, double T, double prec) {
  fermion_thermo *ptr=(fermion_thermo *)vptr;
  fermion *f=(fermion *)ptr_f;
  bool ret=ptr->calc_mu_deg(*f,T,prec);
  return ret;
}

bool o2scl_fermion_thermo_calc_mu_ndeg(void *vptr, void *ptr_f, double T, double prec, bool inc_antip) {
  fermion_thermo *ptr=(fermion_thermo *)vptr;
  fermion *f=(fermion *)ptr_f;
  bool ret=ptr->calc_mu_ndeg(*f,T,prec,inc_antip);
  return ret;
}

void o2scl_fermion_thermo_massless_calc_mu(void *vptr, void *ptr_f, double T) {
  fermion_thermo *ptr=(fermion_thermo *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->massless_calc_mu(*f,T);
  return;
}

void o2scl_fermion_thermo_massless_pair_mu(void *vptr, void *ptr_f, double T) {
  fermion_thermo *ptr=(fermion_thermo *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->massless_pair_mu(*f,T);
  return;
}

void o2scl_fermion_thermo_massless_calc_density(void *vptr, void *ptr_f, double T) {
  fermion_thermo *ptr=(fermion_thermo *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->massless_calc_density(*f,T);
  return;
}

void o2scl_fermion_thermo_massless_pair_density(void *vptr, void *ptr_f, double T) {
  fermion_thermo *ptr=(fermion_thermo *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->massless_pair_density(*f,T);
  return;
}

void *o2scl_create_fermion_rel() {
  fermion_rel *ptr=new fermion_rel;
  return ptr;
}

void o2scl_free_fermion_rel(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  delete ptr;
  return;
}

bool o2scl_fermion_rel_get_err_nonconv(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return ptr->err_nonconv;
}

void o2scl_fermion_rel_set_err_nonconv(void *vptr, bool v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  ptr->err_nonconv=v;
  return;
}

double o2scl_fermion_rel_get_min_psi(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return ptr->min_psi;
}

void o2scl_fermion_rel_set_min_psi(void *vptr, double v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  ptr->min_psi=v;
  return;
}

double o2scl_fermion_rel_get_deg_limit(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return ptr->deg_limit;
}

void o2scl_fermion_rel_set_deg_limit(void *vptr, double v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  ptr->deg_limit=v;
  return;
}

double o2scl_fermion_rel_get_upper_limit_fac(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return ptr->upper_limit_fac;
}

void o2scl_fermion_rel_set_upper_limit_fac(void *vptr, double v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  ptr->upper_limit_fac=v;
  return;
}

int o2scl_fermion_rel_get_verbose(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return ptr->verbose;
}

void o2scl_fermion_rel_set_verbose(void *vptr, int v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  ptr->verbose=v;
  return;
}

bool o2scl_fermion_rel_get_use_expansions(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return ptr->use_expansions;
}

void o2scl_fermion_rel_set_use_expansions(void *vptr, bool v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  ptr->use_expansions=v;
  return;
}

double o2scl_fermion_rel_get_tol_expan(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return ptr->tol_expan;
}

void o2scl_fermion_rel_set_tol_expan(void *vptr, double v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  ptr->tol_expan=v;
  return;
}

bool o2scl_fermion_rel_get_verify_ti(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return ptr->verify_ti;
}

void o2scl_fermion_rel_set_verify_ti(void *vptr, bool v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  ptr->verify_ti=v;
  return;
}

double o2scl_fermion_rel_get_therm_ident(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return ptr->therm_ident;
}

void o2scl_fermion_rel_set_therm_ident(void *vptr, double v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  ptr->therm_ident=v;
  return;
}

void *o2scl_fermion_rel_get_unc(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  return (void *)(&(ptr->unc));
}

void o2scl_fermion_rel_set_unc(void *vptr, void *p_v) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  fermion *p_tsot=(fermion *)p_v;
  ptr->unc=*(p_tsot);
  return;
}

int o2scl_fermion_rel_nu_from_n(void *vptr, void *ptr_f, double T) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  fermion *f=(fermion *)ptr_f;
  int ret=ptr->nu_from_n(*f,T);
  return ret;
}

int o2scl_fermion_rel_calc_density(void *vptr, void *ptr_f, double T) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  fermion *f=(fermion *)ptr_f;
  int ret=ptr->calc_density(*f,T);
  return ret;
}

int o2scl_fermion_rel_pair_density(void *vptr, void *ptr_f, double T) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  fermion *f=(fermion *)ptr_f;
  int ret=ptr->pair_density(*f,T);
  return ret;
}

void o2scl_fermion_rel_calc_mu(void *vptr, void *ptr_f, double T) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->calc_mu(*f,T);
  return;
}

void o2scl_fermion_rel_pair_mu(void *vptr, void *ptr_f, double T) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->pair_mu(*f,T);
  return;
}

void *o2scl_create_fermion_nonrel() {
  fermion_nonrel *ptr=new fermion_nonrel;
  return ptr;
}

void o2scl_free_fermion_nonrel(void *vptr) {
  fermion_nonrel *ptr=(fermion_nonrel *)vptr;
  delete ptr;
  return;
}

int o2scl_fermion_nonrel_calc_density(void *vptr, void *ptr_f, double T) {
  fermion_nonrel *ptr=(fermion_nonrel *)vptr;
  fermion *f=(fermion *)ptr_f;
  int ret=ptr->calc_density(*f,T);
  return ret;
}

void o2scl_fermion_nonrel_calc_mu(void *vptr, void *ptr_f, double T) {
  fermion_nonrel *ptr=(fermion_nonrel *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->calc_mu(*f,T);
  return;
}

void o2scl_fermion_nonrel_nu_from_n(void *vptr, void *ptr_f, double T) {
  fermion_nonrel *ptr=(fermion_nonrel *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->nu_from_n(*f,T);
  return;
}

void *o2scl_create_boson() {
  boson *ptr=new boson;
  return ptr;
}

void o2scl_free_boson(void *vptr) {
  boson *ptr=(boson *)vptr;
  delete ptr;
  return;
}

double o2scl_boson_get_co(void *vptr) {
  boson *ptr=(boson *)vptr;
  return ptr->co;
}

void o2scl_boson_set_co(void *vptr, double v) {
  boson *ptr=(boson *)vptr;
  ptr->co=v;
  return;
}

void *o2scl_create_boson_rel() {
  boson_rel *ptr=new boson_rel;
  return ptr;
}

void o2scl_free_boson_rel(void *vptr) {
  boson_rel *ptr=(boson_rel *)vptr;
  delete ptr;
  return;
}

void o2scl_boson_rel_calc_density(void *vptr, void *ptr_b, double T) {
  boson_rel *ptr=(boson_rel *)vptr;
  boson *b=(boson *)ptr_b;
  ptr->calc_density(*b,T);
  return;
}

void o2scl_boson_rel_calc_mu(void *vptr, void *ptr_b, double T) {
  boson_rel *ptr=(boson_rel *)vptr;
  boson *b=(boson *)ptr_b;
  ptr->calc_mu(*b,T);
  return;
}

void o2scl_boson_rel_nu_from_n(void *vptr, void *ptr_b, double T) {
  boson_rel *ptr=(boson_rel *)vptr;
  boson *b=(boson *)ptr_b;
  ptr->nu_from_n(*b,T);
  return;
}

void o2scl_boson_rel_pair_density(void *vptr, void *ptr_b, double T) {
  boson_rel *ptr=(boson_rel *)vptr;
  boson *b=(boson *)ptr_b;
  ptr->pair_density(*b,T);
  return;
}

void o2scl_boson_rel_pair_mu(void *vptr, void *ptr_b, double T) {
  boson_rel *ptr=(boson_rel *)vptr;
  boson *b=(boson *)ptr_b;
  ptr->pair_mu(*b,T);
  return;
}

void *o2scl_create_classical_thermo() {
  classical_thermo *ptr=new classical_thermo;
  return ptr;
}

void o2scl_free_classical_thermo(void *vptr) {
  classical_thermo *ptr=(classical_thermo *)vptr;
  delete ptr;
  return;
}

void o2scl_classical_thermo_calc_density(void *vptr, void *ptr_p, double T) {
  classical_thermo *ptr=(classical_thermo *)vptr;
  part *p=(part *)ptr_p;
  ptr->calc_density(*p,T);
  return;
}

void o2scl_classical_thermo_calc_mu(void *vptr, void *ptr_p, double T) {
  classical_thermo *ptr=(classical_thermo *)vptr;
  part *p=(part *)ptr_p;
  ptr->calc_mu(*p,T);
  return;
}

void *o2scl_create_thermo_np_deriv_press() {
  thermo_np_deriv_press *ptr=new thermo_np_deriv_press;
  return ptr;
}

void o2scl_free_thermo_np_deriv_press(void *vptr) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  delete ptr;
  return;
}

double o2scl_thermo_np_deriv_press_get_dsdT(void *vptr) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  return ptr->dsdT;
}

void o2scl_thermo_np_deriv_press_set_dsdT(void *vptr, double v) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  ptr->dsdT=v;
  return;
}

double o2scl_thermo_np_deriv_press_get_dnndT(void *vptr) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  return ptr->dnndT;
}

void o2scl_thermo_np_deriv_press_set_dnndT(void *vptr, double v) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  ptr->dnndT=v;
  return;
}

double o2scl_thermo_np_deriv_press_get_dnpdT(void *vptr) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  return ptr->dnpdT;
}

void o2scl_thermo_np_deriv_press_set_dnpdT(void *vptr, double v) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  ptr->dnpdT=v;
  return;
}

double o2scl_thermo_np_deriv_press_get_dnndmun(void *vptr) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  return ptr->dnndmun;
}

void o2scl_thermo_np_deriv_press_set_dnndmun(void *vptr, double v) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  ptr->dnndmun=v;
  return;
}

double o2scl_thermo_np_deriv_press_get_dndmu_mixed(void *vptr) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  return ptr->dndmu_mixed;
}

void o2scl_thermo_np_deriv_press_set_dndmu_mixed(void *vptr, double v) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  ptr->dndmu_mixed=v;
  return;
}

double o2scl_thermo_np_deriv_press_get_dnpdmup(void *vptr) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  return ptr->dnpdmup;
}

void o2scl_thermo_np_deriv_press_set_dnpdmup(void *vptr, double v) {
  thermo_np_deriv_press *ptr=(thermo_np_deriv_press *)vptr;
  ptr->dnpdmup=v;
  return;
}

void *o2scl_create_thermo_np_deriv_helm() {
  thermo_np_deriv_helm *ptr=new thermo_np_deriv_helm;
  return ptr;
}

void o2scl_free_thermo_np_deriv_helm(void *vptr) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  delete ptr;
  return;
}

double o2scl_thermo_np_deriv_helm_get_dsdT(void *vptr) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  return ptr->dsdT;
}

void o2scl_thermo_np_deriv_helm_set_dsdT(void *vptr, double v) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  ptr->dsdT=v;
  return;
}

double o2scl_thermo_np_deriv_helm_get_dmundT(void *vptr) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  return ptr->dmundT;
}

void o2scl_thermo_np_deriv_helm_set_dmundT(void *vptr, double v) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  ptr->dmundT=v;
  return;
}

double o2scl_thermo_np_deriv_helm_get_dmupdT(void *vptr) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  return ptr->dmupdT;
}

void o2scl_thermo_np_deriv_helm_set_dmupdT(void *vptr, double v) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  ptr->dmupdT=v;
  return;
}

double o2scl_thermo_np_deriv_helm_get_dmundnn(void *vptr) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  return ptr->dmundnn;
}

void o2scl_thermo_np_deriv_helm_set_dmundnn(void *vptr, double v) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  ptr->dmundnn=v;
  return;
}

double o2scl_thermo_np_deriv_helm_get_dmudn_mixed(void *vptr) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  return ptr->dmudn_mixed;
}

void o2scl_thermo_np_deriv_helm_set_dmudn_mixed(void *vptr, double v) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  ptr->dmudn_mixed=v;
  return;
}

double o2scl_thermo_np_deriv_helm_get_dmupdnp(void *vptr) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  return ptr->dmupdnp;
}

void o2scl_thermo_np_deriv_helm_set_dmupdnp(void *vptr, double v) {
  thermo_np_deriv_helm *ptr=(thermo_np_deriv_helm *)vptr;
  ptr->dmupdnp=v;
  return;
}

void *o2scl_create_part_deriv_press() {
  part_deriv_press *ptr=new part_deriv_press;
  return ptr;
}

void o2scl_free_part_deriv_press(void *vptr) {
  part_deriv_press *ptr=(part_deriv_press *)vptr;
  delete ptr;
  return;
}

double o2scl_part_deriv_press_get_dndmu(void *vptr) {
  part_deriv_press *ptr=(part_deriv_press *)vptr;
  return ptr->dndmu;
}

void o2scl_part_deriv_press_set_dndmu(void *vptr, double v) {
  part_deriv_press *ptr=(part_deriv_press *)vptr;
  ptr->dndmu=v;
  return;
}

double o2scl_part_deriv_press_get_dndT(void *vptr) {
  part_deriv_press *ptr=(part_deriv_press *)vptr;
  return ptr->dndT;
}

void o2scl_part_deriv_press_set_dndT(void *vptr, double v) {
  part_deriv_press *ptr=(part_deriv_press *)vptr;
  ptr->dndT=v;
  return;
}

double o2scl_part_deriv_press_get_dsdT(void *vptr) {
  part_deriv_press *ptr=(part_deriv_press *)vptr;
  return ptr->dsdT;
}

void o2scl_part_deriv_press_set_dsdT(void *vptr, double v) {
  part_deriv_press *ptr=(part_deriv_press *)vptr;
  ptr->dsdT=v;
  return;
}

void o2scl_part_deriv_press_deriv_f(void *vptr, double *dmudn, double *dmudT, double *dsdT_n) {
  part_deriv_press *ptr=(part_deriv_press *)vptr;
  ptr->deriv_f(*dmudn,*dmudT,*dsdT_n);
  return;
}

void *o2scl_create_part_deriv() {
  part_deriv *ptr=new part_deriv;
  return ptr;
}

void o2scl_free_part_deriv(void *vptr) {
  part_deriv *ptr=(part_deriv *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_fermion_deriv() {
  fermion_deriv *ptr=new fermion_deriv;
  return ptr;
}

void o2scl_free_fermion_deriv(void *vptr) {
  fermion_deriv *ptr=(fermion_deriv *)vptr;
  delete ptr;
  return;
}

void *o2scl_create_deriv_thermo_base() {
  deriv_thermo_base *ptr=new deriv_thermo_base;
  return ptr;
}

void o2scl_free_deriv_thermo_base(void *vptr) {
  deriv_thermo_base *ptr=(deriv_thermo_base *)vptr;
  delete ptr;
  return;
}

double o2scl_deriv_thermo_base_heat_cap_ppart_const_vol(void *vptr, void *ptr_p, double T) {
  deriv_thermo_base *ptr=(deriv_thermo_base *)vptr;
  part_deriv *p=(part_deriv *)ptr_p;
  double ret=ptr->heat_cap_ppart_const_vol(*p,T);
  return ret;
}

double o2scl_deriv_thermo_base_heat_cap_ppart_const_press(void *vptr, void *ptr_p, double T) {
  deriv_thermo_base *ptr=(deriv_thermo_base *)vptr;
  part_deriv *p=(part_deriv *)ptr_p;
  double ret=ptr->heat_cap_ppart_const_press(*p,T);
  return ret;
}

double o2scl_deriv_thermo_base_compress_adiabatic(void *vptr, void *ptr_p, double T) {
  deriv_thermo_base *ptr=(deriv_thermo_base *)vptr;
  part_deriv *p=(part_deriv *)ptr_p;
  double ret=ptr->compress_adiabatic(*p,T);
  return ret;
}

double o2scl_deriv_thermo_base_compress_const_tptr(void *vptr, void *ptr_p, double T) {
  deriv_thermo_base *ptr=(deriv_thermo_base *)vptr;
  part_deriv *p=(part_deriv *)ptr_p;
  double ret=ptr->compress_const_tptr(*p,T);
  return ret;
}

double o2scl_deriv_thermo_base_coeff_thermal_exp(void *vptr, void *ptr_p, double T) {
  deriv_thermo_base *ptr=(deriv_thermo_base *)vptr;
  part_deriv *p=(part_deriv *)ptr_p;
  double ret=ptr->coeff_thermal_exp(*p,T);
  return ret;
}

double o2scl_deriv_thermo_base_squared_sound_speed(void *vptr, void *ptr_p, double T) {
  deriv_thermo_base *ptr=(deriv_thermo_base *)vptr;
  part_deriv *p=(part_deriv *)ptr_p;
  double ret=ptr->squared_sound_speed(*p,T);
  return ret;
}

void *o2scl_create_classical_deriv_thermo() {
  classical_deriv_thermo *ptr=new classical_deriv_thermo;
  return ptr;
}

void o2scl_free_classical_deriv_thermo(void *vptr) {
  classical_deriv_thermo *ptr=(classical_deriv_thermo *)vptr;
  delete ptr;
  return;
}

void o2scl_classical_deriv_thermo_calc_density(void *vptr, void *ptr_p, double T) {
  classical_deriv_thermo *ptr=(classical_deriv_thermo *)vptr;
  part_deriv *p=(part_deriv *)ptr_p;
  ptr->calc_density(*p,T);
  return;
}

void o2scl_classical_deriv_thermo_calc_mu(void *vptr, void *ptr_p, double T) {
  classical_deriv_thermo *ptr=(classical_deriv_thermo *)vptr;
  part_deriv *p=(part_deriv *)ptr_p;
  ptr->calc_mu(*p,T);
  return;
}

void *o2scl_create_fermion_mag_zerot() {
  fermion_mag_zerot *ptr=new fermion_mag_zerot;
  return ptr;
}

void o2scl_free_fermion_mag_zerot(void *vptr) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  delete ptr;
  return;
}

int o2scl_fermion_mag_zerot_get_nmax_up(void *vptr) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  return ptr->nmax_up;
}

void o2scl_fermion_mag_zerot_set_nmax_up(void *vptr, int v) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  ptr->nmax_up=v;
  return;
}

int o2scl_fermion_mag_zerot_get_nmax_dn(void *vptr) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  return ptr->nmax_dn;
}

void o2scl_fermion_mag_zerot_set_nmax_dn(void *vptr, int v) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  ptr->nmax_dn=v;
  return;
}

int o2scl_fermion_mag_zerot_get_sum_limit(void *vptr) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  return ptr->sum_limit;
}

void o2scl_fermion_mag_zerot_set_sum_limit(void *vptr, int v) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  ptr->sum_limit=v;
  return;
}

void o2scl_fermion_mag_zerot_calc_mu_zerot_mag(void *vptr, void *ptr_f, double qB, double kappa) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->calc_mu_zerot_mag(*f,qB,kappa);
  return;
}

void o2scl_fermion_mag_zerot_calc_density_zerot_mag(void *vptr, void *ptr_f, double qB, double kappa) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->calc_density_zerot_mag(*f,qB,kappa);
  return;
}

