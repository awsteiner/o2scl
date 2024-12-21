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

#include <o2scl/eos_python.h>
#include <o2scl/part.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

void *o2scl_eos_base_get_def_thermo(void *vptr) {
  eos_base *ptr=(eos_base *)vptr;
  return (void *)(&(ptr->def_thermo));
}

void o2scl_eos_base_set_def_thermo(void *vptr, void *p_v) {
  eos_base *ptr=(eos_base *)vptr;
  o2scl::thermo *p_tsot=(o2scl::thermo *)p_v;
  ptr->def_thermo=*(p_tsot);
  return;
}

void *o2scl_create_eos_leptons() {
  eos_leptons *ptr=new eos_leptons;
  return ptr;
}

void o2scl_free_eos_leptons(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  delete ptr;
  return;
}

void *o2scl_eos_leptons_get_th(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return (void *)(&(ptr->th));
}

void o2scl_eos_leptons_set_th(void *vptr, void *p_v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  o2scl::thermo *p_tsot=(o2scl::thermo *)p_v;
  ptr->th=*(p_tsot);
  return;
}

void *o2scl_eos_leptons_get_e(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return (void *)(&(ptr->e));
}

void o2scl_eos_leptons_set_e(void *vptr, void *p_v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  o2scl::fermion *p_tsot=(o2scl::fermion *)p_v;
  ptr->e=*(p_tsot);
  return;
}

void *o2scl_eos_leptons_get_mu(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return (void *)(&(ptr->mu));
}

void o2scl_eos_leptons_set_mu(void *vptr, void *p_v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  o2scl::fermion *p_tsot=(o2scl::fermion *)p_v;
  ptr->mu=*(p_tsot);
  return;
}

void *o2scl_eos_leptons_get_ph(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return (void *)(&(ptr->ph));
}

void o2scl_eos_leptons_set_ph(void *vptr, void *p_v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  o2scl::boson *p_tsot=(o2scl::boson *)p_v;
  ptr->ph=*(p_tsot);
  return;
}

void *o2scl_eos_leptons_get_ed(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return (void *)(&(ptr->ed));
}

void o2scl_eos_leptons_set_ed(void *vptr, void *p_v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  o2scl::part_deriv_press *p_tsot=(o2scl::part_deriv_press *)p_v;
  ptr->ed=*(p_tsot);
  return;
}

void *o2scl_eos_leptons_get_mud(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return (void *)(&(ptr->mud));
}

void o2scl_eos_leptons_set_mud(void *vptr, void *p_v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  o2scl::part_deriv_press *p_tsot=(o2scl::part_deriv_press *)p_v;
  ptr->mud=*(p_tsot);
  return;
}

void *o2scl_eos_leptons_get_phd(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return (void *)(&(ptr->phd));
}

void o2scl_eos_leptons_set_phd(void *vptr, void *p_v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  o2scl::part_deriv_press *p_tsot=(o2scl::part_deriv_press *)p_v;
  ptr->phd=*(p_tsot);
  return;
}

bool o2scl_eos_leptons_get_include_muons(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return ptr->include_muons;
}

void o2scl_eos_leptons_set_include_muons(void *vptr, bool v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  ptr->include_muons=v;
  return;
}

bool o2scl_eos_leptons_get_include_deriv(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return ptr->include_deriv;
}

void o2scl_eos_leptons_set_include_deriv(void *vptr, bool v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  ptr->include_deriv=v;
  return;
}

bool o2scl_eos_leptons_get_pde_from_density(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return ptr->pde_from_density;
}

void o2scl_eos_leptons_set_pde_from_density(void *vptr, bool v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  ptr->pde_from_density=v;
  return;
}

int o2scl_eos_leptons_get_verbose(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return ptr->verbose;
}

void o2scl_eos_leptons_set_verbose(void *vptr, int v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  ptr->verbose=v;
  return;
}

bool o2scl_eos_leptons_get_err_nonconv(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return ptr->err_nonconv;
}

void o2scl_eos_leptons_set_err_nonconv(void *vptr, bool v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  ptr->err_nonconv=v;
  return;
}

void *o2scl_eos_leptons_get_frel(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  return (void *)(&(ptr->frel));
}

void o2scl_eos_leptons_set_frel(void *vptr, void *p_v) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  o2scl::fermion_rel *p_tsot=(o2scl::fermion_rel *)p_v;
  ptr->frel=*(p_tsot);
  return;
}

void o2scl_eos_leptons_default_acc(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  ptr->default_acc();
  return;
}

void o2scl_eos_leptons_improved_acc(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  ptr->improved_acc();
  return;
}

void o2scl_eos_leptons_ld_acc(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  ptr->ld_acc();
  return;
}

void o2scl_eos_leptons_fp_25_acc(void *vptr) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  ptr->fp_25_acc();
  return;
}

int o2scl_eos_leptons_pair_mu(void *vptr, double T) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  int ret=ptr->pair_mu(T);
  return ret;
}

int o2scl_eos_leptons_pair_mu_eq(void *vptr, double T) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  int ret=ptr->pair_mu_eq(T);
  return ret;
}

int o2scl_eos_leptons_pair_density(void *vptr, double T) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  int ret=ptr->pair_density(T);
  return ret;
}

int o2scl_eos_leptons_pair_density_eq(void *vptr, double nq, double T) {
  eos_leptons *ptr=(eos_leptons *)vptr;
  int ret=ptr->pair_density_eq(nq,T);
  return ret;
}

double o2scl_eos_had_base_get_eoa(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->eoa;
}

void o2scl_eos_had_base_set_eoa(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->eoa=v;
  return;
}

double o2scl_eos_had_base_get_msom(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->msom;
}

void o2scl_eos_had_base_set_msom(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->msom=v;
  return;
}

double o2scl_eos_had_base_get_comp(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->comp;
}

void o2scl_eos_had_base_set_comp(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->comp=v;
  return;
}

double o2scl_eos_had_base_get_n0(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->n0;
}

void o2scl_eos_had_base_set_n0(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->n0=v;
  return;
}

double o2scl_eos_had_base_get_esym(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->esym;
}

void o2scl_eos_had_base_set_esym(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->esym=v;
  return;
}

double o2scl_eos_had_base_get_kprime(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->kprime;
}

void o2scl_eos_had_base_set_kprime(void *vptr, double v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->kprime=v;
  return;
}

bool o2scl_eos_had_base_get_err_nonconv(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return ptr->err_nonconv;
}

void o2scl_eos_had_base_set_err_nonconv(void *vptr, bool v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->err_nonconv=v;
  return;
}

void *o2scl_eos_had_base_get_def_neutron(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return (void *)(&(ptr->def_neutron));
}

void o2scl_eos_had_base_set_def_neutron(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_tsot=(o2scl::fermion *)p_v;
  ptr->def_neutron=*(p_tsot);
  return;
}

void *o2scl_eos_had_base_get_def_proton(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  return (void *)(&(ptr->def_proton));
}

void o2scl_eos_had_base_set_def_proton(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_tsot=(o2scl::fermion *)p_v;
  ptr->def_proton=*(p_tsot);
  return;
}

int o2scl_eos_had_base_calc_e(void *vptr, void *ptr_n, void *ptr_p, void *ptr_th) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *n=(o2scl::fermion *)ptr_n;
  o2scl::fermion *p=(o2scl::fermion *)ptr_p;
  o2scl::thermo *th=(o2scl::thermo *)ptr_th;
  int ret=ptr->calc_e(*n,*p,*th);
  return ret;
}

int o2scl_eos_had_base_calc_p(void *vptr, void *ptr_n, void *ptr_p, void *ptr_th) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *n=(o2scl::fermion *)ptr_n;
  o2scl::fermion *p=(o2scl::fermion *)ptr_p;
  o2scl::thermo *th=(o2scl::thermo *)ptr_th;
  int ret=ptr->calc_p(*n,*p,*th);
  return ret;
}

double o2scl_eos_had_base_fcomp(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fcomp(nb,delta);
  return ret;
}

double o2scl_eos_had_base_fcomp_err(void *vptr, double nb, double delta, double *unc) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fcomp_err(nb,delta,*unc);
  return ret;
}

double o2scl_eos_had_base_feoa(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->feoa(nb,delta);
  return ret;
}

double o2scl_eos_had_base_fesym(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fesym(nb,delta);
  return ret;
}

double o2scl_eos_had_base_fesym_err(void *vptr, double nb, double delta, double *unc) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fesym_err(nb,delta,*unc);
  return ret;
}

double o2scl_eos_had_base_fesym_slope(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fesym_slope(nb,delta);
  return ret;
}

double o2scl_eos_had_base_fesym_curve(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fesym_curve(nb,delta);
  return ret;
}

double o2scl_eos_had_base_fesym_skew(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fesym_skew(nb,delta);
  return ret;
}

double o2scl_eos_had_base_fesym_diff(void *vptr, double nb) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fesym_diff(nb);
  return ret;
}

double o2scl_eos_had_base_feta(void *vptr, double nb) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->feta(nb);
  return ret;
}

double o2scl_eos_had_base_feta_prime(void *vptr, double nb) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->feta_prime(nb);
  return ret;
}

double o2scl_eos_had_base_fkprime(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fkprime(nb,delta);
  return ret;
}

double o2scl_eos_had_base_fmsom(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->fmsom(nb,delta);
  return ret;
}

double o2scl_eos_had_base_f_effm_neut(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->f_effm_neut(nb,delta);
  return ret;
}

double o2scl_eos_had_base_f_effm_prot(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->f_effm_prot(nb,delta);
  return ret;
}

double o2scl_eos_had_base_f_effm_scalar(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->f_effm_scalar(nb,delta);
  return ret;
}

double o2scl_eos_had_base_f_effm_vector(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->f_effm_vector(nb,delta);
  return ret;
}

int o2scl_eos_had_base_fn0(void *vptr, double delta, double *nb, double *leoa) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  int ret=ptr->fn0(delta,*nb,*leoa);
  return ret;
}

void o2scl_eos_had_base_f_number_suscept(void *vptr, double mun, double mup, double *dPdnn, double *dPdnp, double *dPdpp) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->f_number_suscept(mun,mup,*dPdnn,*dPdnp,*dPdpp);
  return;
}

void o2scl_eos_had_base_f_inv_number_suscept(void *vptr, double mun, double mup, double *dednn, double *dednp, double *dedpp) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  ptr->f_inv_number_suscept(mun,mup,*dednn,*dednp,*dedpp);
  return;
}

int o2scl_eos_had_base_saturation(void *vptr) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  int ret=ptr->saturation();
  return ret;
}

double o2scl_eos_had_base_calc_mun_e(void *vptr, double nn, double np) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_mun_e(nn,np);
  return ret;
}

double o2scl_eos_had_base_calc_mup_e(void *vptr, double nn, double np) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_mup_e(nn,np);
  return ret;
}

double o2scl_eos_had_base_calc_ed(void *vptr, double nn, double np) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_ed(nn,np);
  return ret;
}

double o2scl_eos_had_base_calc_pr(void *vptr, double nn, double np) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_pr(nn,np);
  return ret;
}

double o2scl_eos_had_base_calc_nn_p(void *vptr, double mun, double mup) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_nn_p(mun,mup);
  return ret;
}

double o2scl_eos_had_base_calc_np_p(void *vptr, double nn, double mup) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_np_p(nn,mup);
  return ret;
}

double o2scl_eos_had_base_calc_dmu_delta(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_dmu_delta(nb,delta);
  return ret;
}

double o2scl_eos_had_base_calc_musum_delta(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_musum_delta(nb,delta);
  return ret;
}

double o2scl_eos_had_base_calc_pressure_nb(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_pressure_nb(nb,delta);
  return ret;
}

double o2scl_eos_had_base_calc_edensity_nb(void *vptr, double nb, double delta) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double ret=ptr->calc_edensity_nb(nb,delta);
  return ret;
}

int o2scl_eos_had_temp_base_calc_temp_e(void *vptr, void *ptr_n, void *ptr_p, double T, void *ptr_th) {
  eos_had_temp_base *ptr=(eos_had_temp_base *)vptr;
  fermion *n=(fermion *)ptr_n;
  fermion *p=(fermion *)ptr_p;
  thermo *th=(thermo *)ptr_th;
  int ret=ptr->calc_temp_e(*n,*p,T,*th);
  return ret;
}

int o2scl_eos_had_temp_base_calc_temp_p(void *vptr, void *ptr_n, void *ptr_p, double T, void *ptr_th) {
  eos_had_temp_base *ptr=(eos_had_temp_base *)vptr;
  fermion *n=(fermion *)ptr_n;
  fermion *p=(fermion *)ptr_p;
  thermo *th=(thermo *)ptr_th;
  int ret=ptr->calc_temp_p(*n,*p,T,*th);
  return ret;
}

void *o2scl_create_eos_had_skyrme() {
  eos_had_skyrme *ptr=new eos_had_skyrme;
  return ptr;
}

void o2scl_free_eos_had_skyrme(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  delete ptr;
  return;
}

double o2scl_eos_had_skyrme_get_t0(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->t0;
}

void o2scl_eos_had_skyrme_set_t0(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->t0=v;
  return;
}

double o2scl_eos_had_skyrme_get_t1(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->t1;
}

void o2scl_eos_had_skyrme_set_t1(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->t1=v;
  return;
}

double o2scl_eos_had_skyrme_get_t2(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->t2;
}

void o2scl_eos_had_skyrme_set_t2(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->t2=v;
  return;
}

double o2scl_eos_had_skyrme_get_t3(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->t3;
}

void o2scl_eos_had_skyrme_set_t3(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->t3=v;
  return;
}

double o2scl_eos_had_skyrme_get_x0(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->x0;
}

void o2scl_eos_had_skyrme_set_x0(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->x0=v;
  return;
}

double o2scl_eos_had_skyrme_get_x1(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->x1;
}

void o2scl_eos_had_skyrme_set_x1(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->x1=v;
  return;
}

double o2scl_eos_had_skyrme_get_x2(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->x2;
}

void o2scl_eos_had_skyrme_set_x2(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->x2=v;
  return;
}

double o2scl_eos_had_skyrme_get_x3(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->x3;
}

void o2scl_eos_had_skyrme_set_x3(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->x3=v;
  return;
}

double o2scl_eos_had_skyrme_get_alpha(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->alpha;
}

void o2scl_eos_had_skyrme_set_alpha(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->alpha=v;
  return;
}

double o2scl_eos_had_skyrme_get_a(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->a;
}

void o2scl_eos_had_skyrme_set_a(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->a=v;
  return;
}

double o2scl_eos_had_skyrme_get_b(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->b;
}

void o2scl_eos_had_skyrme_set_b(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->b=v;
  return;
}

double o2scl_eos_had_skyrme_get_W0(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->W0;
}

void o2scl_eos_had_skyrme_set_W0(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->W0=v;
  return;
}

double o2scl_eos_had_skyrme_get_b4(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->b4;
}

void o2scl_eos_had_skyrme_set_b4(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->b4=v;
  return;
}

double o2scl_eos_had_skyrme_get_b4p(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->b4p;
}

void o2scl_eos_had_skyrme_set_b4p(void *vptr, double v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->b4p=v;
  return;
}

bool o2scl_eos_had_skyrme_get_parent_method(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return ptr->parent_method;
}

void o2scl_eos_had_skyrme_set_parent_method(void *vptr, bool v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  ptr->parent_method=v;
  return;
}

void *o2scl_eos_had_skyrme_get_reference(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  // The ownership of the string pointer is passed to the Python class
  // and the memory is freed later.
  std::string *sptr=new std::string;
  *sptr=ptr->reference;
  return sptr;
}

void o2scl_eos_had_skyrme_set_reference(void *vptr, void *p_v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  std::string *p_tsot=(std::string *)p_v;
  ptr->reference=*(p_tsot);
  return;
}

void *o2scl_eos_had_skyrme_get_nrfd(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  return (void *)(&(ptr->nrfd));
}

void o2scl_eos_had_skyrme_set_nrfd(void *vptr, void *p_v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  o2scl::fermion_deriv_nr *p_tsot=(o2scl::fermion_deriv_nr *)p_v;
  ptr->nrfd=*(p_tsot);
  return;
}

void *o2scl_create_eos_had_apr() {
  eos_had_apr *ptr=new eos_had_apr;
  return ptr;
}

void o2scl_free_eos_had_apr(void *vptr) {
  eos_had_apr *ptr=(eos_had_apr *)vptr;
  delete ptr;
  return;
}

int o2scl_eos_had_apr_get_pion(void *vptr) {
  eos_had_apr *ptr=(eos_had_apr *)vptr;
  return ptr->pion;
}

void o2scl_eos_had_apr_set_pion(void *vptr, int v) {
  eos_had_apr *ptr=(eos_had_apr *)vptr;
  ptr->pion=v;
  return;
}

bool o2scl_eos_had_apr_get_parent_method(void *vptr) {
  eos_had_apr *ptr=(eos_had_apr *)vptr;
  return ptr->parent_method;
}

void o2scl_eos_had_apr_set_parent_method(void *vptr, bool v) {
  eos_had_apr *ptr=(eos_had_apr *)vptr;
  ptr->parent_method=v;
  return;
}

void *o2scl_create_eos_had_rmf() {
  eos_had_rmf *ptr=new eos_had_rmf;
  return ptr;
}

void o2scl_free_eos_had_rmf(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  delete ptr;
  return;
}

size_t o2scl_eos_had_rmf_get_calc_e_steps(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->calc_e_steps;
}

void o2scl_eos_had_rmf_set_calc_e_steps(void *vptr, size_t v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->calc_e_steps=v;
  return;
}

bool o2scl_eos_had_rmf_get_calc_e_relative(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->calc_e_relative;
}

void o2scl_eos_had_rmf_set_calc_e_relative(void *vptr, bool v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->calc_e_relative=v;
  return;
}

bool o2scl_eos_had_rmf_get_zm_mode(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->zm_mode;
}

void o2scl_eos_had_rmf_set_zm_mode(void *vptr, bool v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->zm_mode=v;
  return;
}

int o2scl_eos_had_rmf_get_verbose(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->verbose;
}

void o2scl_eos_had_rmf_set_verbose(void *vptr, int v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->verbose=v;
  return;
}

double o2scl_eos_had_rmf_get_mnuc(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->mnuc;
}

void o2scl_eos_had_rmf_set_mnuc(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->mnuc=v;
  return;
}

double o2scl_eos_had_rmf_get_ms(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->ms;
}

void o2scl_eos_had_rmf_set_ms(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->ms=v;
  return;
}

double o2scl_eos_had_rmf_get_mw(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->mw;
}

void o2scl_eos_had_rmf_set_mw(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->mw=v;
  return;
}

double o2scl_eos_had_rmf_get_mr(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->mr;
}

void o2scl_eos_had_rmf_set_mr(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->mr=v;
  return;
}

double o2scl_eos_had_rmf_get_cs(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->cs;
}

void o2scl_eos_had_rmf_set_cs(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->cs=v;
  return;
}

double o2scl_eos_had_rmf_get_cw(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->cw;
}

void o2scl_eos_had_rmf_set_cw(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->cw=v;
  return;
}

double o2scl_eos_had_rmf_get_cr(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->cr;
}

void o2scl_eos_had_rmf_set_cr(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->cr=v;
  return;
}

double o2scl_eos_had_rmf_get_b(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->b;
}

void o2scl_eos_had_rmf_set_b(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->b=v;
  return;
}

double o2scl_eos_had_rmf_get_c(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->c;
}

void o2scl_eos_had_rmf_set_c(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->c=v;
  return;
}

double o2scl_eos_had_rmf_get_zeta(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->zeta;
}

void o2scl_eos_had_rmf_set_zeta(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->zeta=v;
  return;
}

double o2scl_eos_had_rmf_get_xi(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->xi;
}

void o2scl_eos_had_rmf_set_xi(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->xi=v;
  return;
}

double o2scl_eos_had_rmf_get_a1(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->a1;
}

void o2scl_eos_had_rmf_set_a1(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->a1=v;
  return;
}

double o2scl_eos_had_rmf_get_a2(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->a2;
}

void o2scl_eos_had_rmf_set_a2(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->a2=v;
  return;
}

double o2scl_eos_had_rmf_get_a3(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->a3;
}

void o2scl_eos_had_rmf_set_a3(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->a3=v;
  return;
}

double o2scl_eos_had_rmf_get_a4(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->a4;
}

void o2scl_eos_had_rmf_set_a4(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->a4=v;
  return;
}

double o2scl_eos_had_rmf_get_a5(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->a5;
}

void o2scl_eos_had_rmf_set_a5(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->a5=v;
  return;
}

double o2scl_eos_had_rmf_get_a6(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->a6;
}

void o2scl_eos_had_rmf_set_a6(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->a6=v;
  return;
}

double o2scl_eos_had_rmf_get_b1(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->b1;
}

void o2scl_eos_had_rmf_set_b1(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->b1=v;
  return;
}

double o2scl_eos_had_rmf_get_b2(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->b2;
}

void o2scl_eos_had_rmf_set_b2(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->b2=v;
  return;
}

double o2scl_eos_had_rmf_get_b3(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->b3;
}

void o2scl_eos_had_rmf_set_b3(void *vptr, double v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->b3=v;
  return;
}

int o2scl_eos_had_rmf_get_fields(void *vptr, double *sig, double *ome, double *rho) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  int ret=ptr->get_fields(*sig,*ome,*rho);
  return ret;
}

int o2scl_eos_had_rmf_set_fields(void *vptr, double *sig, double *ome, double *rho) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  int ret=ptr->set_fields(*sig,*ome,*rho);
  return ret;
}

void *o2scl_create_eos_quark_bag() {
  eos_quark_bag *ptr=new eos_quark_bag;
  return ptr;
}

void o2scl_free_eos_quark_bag(void *vptr) {
  eos_quark_bag *ptr=(eos_quark_bag *)vptr;
  delete ptr;
  return;
}

double o2scl_eos_quark_bag_get_bag_constant(void *vptr) {
  eos_quark_bag *ptr=(eos_quark_bag *)vptr;
  return ptr->bag_constant;
}

void o2scl_eos_quark_bag_set_bag_constant(void *vptr, double v) {
  eos_quark_bag *ptr=(eos_quark_bag *)vptr;
  ptr->bag_constant=v;
  return;
}

void *o2scl_create_eos_quark_njl() {
  eos_quark_njl *ptr=new eos_quark_njl;
  return ptr;
}

void o2scl_free_eos_quark_njl(void *vptr) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  delete ptr;
  return;
}

double o2scl_eos_quark_njl_get_B0(void *vptr) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  return ptr->B0;
}

void o2scl_eos_quark_njl_set_B0(void *vptr, double v) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  ptr->B0=v;
  return;
}

double o2scl_eos_quark_njl_get_L(void *vptr) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  return ptr->L;
}

void o2scl_eos_quark_njl_set_L(void *vptr, double v) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  ptr->L=v;
  return;
}

double o2scl_eos_quark_njl_get_G(void *vptr) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  return ptr->G;
}

void o2scl_eos_quark_njl_set_G(void *vptr, double v) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  ptr->G=v;
  return;
}

double o2scl_eos_quark_njl_get_K(void *vptr) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  return ptr->K;
}

void o2scl_eos_quark_njl_set_K(void *vptr, double v) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  ptr->K=v;
  return;
}

double o2scl_eos_quark_njl_get_limit(void *vptr) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  return ptr->limit;
}

void o2scl_eos_quark_njl_set_limit(void *vptr, double v) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  ptr->limit=v;
  return;
}

bool o2scl_eos_quark_njl_get_from_qq(void *vptr) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  return ptr->from_qq;
}

void o2scl_eos_quark_njl_set_from_qq(void *vptr, bool v) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  ptr->from_qq=v;
  return;
}

int o2scl_eos_tov_get_verbose(void *vptr) {
  eos_tov *ptr=(eos_tov *)vptr;
  return ptr->verbose;
}

void o2scl_eos_tov_set_verbose(void *vptr, int v) {
  eos_tov *ptr=(eos_tov *)vptr;
  ptr->verbose=v;
  return;
}

bool o2scl_eos_tov_has_baryons(void *vptr) {
  eos_tov *ptr=(eos_tov *)vptr;
  bool ret=ptr->has_baryons();
  return ret;
}

double o2scl_eos_tov_ed_from_pr(void *vptr, double pr) {
  eos_tov *ptr=(eos_tov *)vptr;
  double ret=ptr->ed_from_pr(pr);
  return ret;
}

double o2scl_eos_tov_pr_from_ed(void *vptr, double ed) {
  eos_tov *ptr=(eos_tov *)vptr;
  double ret=ptr->pr_from_ed(ed);
  return ret;
}

double o2scl_eos_tov_nb_from_ed(void *vptr, double ed) {
  eos_tov *ptr=(eos_tov *)vptr;
  double ret=ptr->nb_from_ed(ed);
  return ret;
}

double o2scl_eos_tov_nb_from_pr(void *vptr, double pr) {
  eos_tov *ptr=(eos_tov *)vptr;
  double ret=ptr->nb_from_pr(pr);
  return ret;
}

double o2scl_eos_tov_ed_from_nb(void *vptr, double nb) {
  eos_tov *ptr=(eos_tov *)vptr;
  double ret=ptr->ed_from_nb(nb);
  return ret;
}

double o2scl_eos_tov_pr_from_nb(void *vptr, double nb) {
  eos_tov *ptr=(eos_tov *)vptr;
  double ret=ptr->pr_from_nb(nb);
  return ret;
}

void o2scl_eos_tov_ed_nb_from_pr(void *vptr, double pr, double *ed, double *nb) {
  eos_tov *ptr=(eos_tov *)vptr;
  ptr->ed_nb_from_pr(pr,*ed,*nb);
  return;
}

void *o2scl_create_eos_tov_buchdahl() {
  eos_tov_buchdahl *ptr=new eos_tov_buchdahl;
  return ptr;
}

void o2scl_free_eos_tov_buchdahl(void *vptr) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  delete ptr;
  return;
}

double o2scl_eos_tov_buchdahl_get_Pstar(void *vptr) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  return ptr->Pstar;
}

void o2scl_eos_tov_buchdahl_set_Pstar(void *vptr, double v) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  ptr->Pstar=v;
  return;
}

double o2scl_eos_tov_buchdahl_get_G_km_Msun(void *vptr) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  return ptr->G_km_Msun;
}

void o2scl_eos_tov_buchdahl_set_G_km_Msun(void *vptr, double v) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  ptr->G_km_Msun=v;
  return;
}

void o2scl_eos_tov_buchdahl_set_baryon_density(void *vptr, double nb, double ed) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  ptr->set_baryon_density(nb,ed);
  return;
}

double o2scl_eos_tov_buchdahl_rad_from_gm(void *vptr, double gm) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  double ret=ptr->rad_from_gm(gm);
  return ret;
}

double o2scl_eos_tov_buchdahl_ed_from_r_gm(void *vptr, double r, double beta) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  double ret=ptr->ed_from_r_gm(r,beta);
  return ret;
}

double o2scl_eos_tov_buchdahl_pr_from_r_gm(void *vptr, double r, double beta) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  double ret=ptr->pr_from_r_gm(r,beta);
  return ret;
}

double o2scl_eos_tov_buchdahl_exp2lam_from_r_gm(void *vptr, double r, double beta) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  double ret=ptr->exp2lam_from_r_gm(r,beta);
  return ret;
}

double o2scl_eos_tov_buchdahl_exp2phi_from_r_gm(void *vptr, double r, double beta) {
  eos_tov_buchdahl *ptr=(eos_tov_buchdahl *)vptr;
  double ret=ptr->exp2phi_from_r_gm(r,beta);
  return ret;
}

void *o2scl_create_eos_tov_polytrope() {
  eos_tov_polytrope *ptr=new eos_tov_polytrope;
  return ptr;
}

void o2scl_free_eos_tov_polytrope(void *vptr) {
  eos_tov_polytrope *ptr=(eos_tov_polytrope *)vptr;
  delete ptr;
  return;
}

void o2scl_eos_tov_polytrope_set_coeff_index(void *vptr, double coeff, double index) {
  eos_tov_polytrope *ptr=(eos_tov_polytrope *)vptr;
  ptr->set_coeff_index(coeff,index);
  return;
}

void *o2scl_create_eos_tov_linear() {
  eos_tov_linear *ptr=new eos_tov_linear;
  return ptr;
}

void o2scl_free_eos_tov_linear(void *vptr) {
  eos_tov_linear *ptr=(eos_tov_linear *)vptr;
  delete ptr;
  return;
}

void o2scl_eos_tov_linear_set_cs2_eps0(void *vptr, double cs2, double eps0) {
  eos_tov_linear *ptr=(eos_tov_linear *)vptr;
  ptr->set_cs2_eps0(cs2,eps0);
  return;
}

void *o2scl_create_eos_tov_interp() {
  eos_tov_interp *ptr=new eos_tov_interp;
  return ptr;
}

void o2scl_free_eos_tov_interp(void *vptr) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  delete ptr;
  return;
}

bool o2scl_eos_tov_interp_get_err_nonconv(void *vptr) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  return ptr->err_nonconv;
}

void o2scl_eos_tov_interp_set_err_nonconv(void *vptr, bool v) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  ptr->err_nonconv=v;
  return;
}

void *o2scl_eos_tov_interp_get_full_vece(void *vptr) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  return (void *)(&(ptr->full_vece));
}

void o2scl_eos_tov_interp_set_full_vece(void *vptr, void *p_v) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  std::vector<double> *p_tsot=(std::vector<double> *)p_v;
  ptr->full_vece=*(p_tsot);
  return;
}

void *o2scl_eos_tov_interp_get_full_vecp(void *vptr) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  return (void *)(&(ptr->full_vecp));
}

void o2scl_eos_tov_interp_set_full_vecp(void *vptr, void *p_v) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  std::vector<double> *p_tsot=(std::vector<double> *)p_v;
  ptr->full_vecp=*(p_tsot);
  return;
}

void *o2scl_eos_tov_interp_get_full_vecnb(void *vptr) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  return (void *)(&(ptr->full_vecnb));
}

void o2scl_eos_tov_interp_set_full_vecnb(void *vptr, void *p_v) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  std::vector<double> *p_tsot=(std::vector<double> *)p_v;
  ptr->full_vecnb=*(p_tsot);
  return;
}

void o2scl_eos_tov_interp_read_table(void *vptr, void *ptr_eos, void *ptr_s_cole, void *ptr_s_colp, void *ptr_s_colnb) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  table_units<> *eos=(table_units<> *)ptr_eos;
  std::string *s_cole=(std::string *)ptr_s_cole;
  std::string *s_colp=(std::string *)ptr_s_colp;
  std::string *s_colnb=(std::string *)ptr_s_colnb;
  ptr->read_table(*eos,*s_cole,*s_colp,*s_colnb);
  return;
}

void o2scl_eos_tov_interp_default_low_dens_eos(void *vptr) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  ptr->default_low_dens_eos();
  return;
}

void o2scl_eos_tov_interp_sho11_low_dens_eos(void *vptr) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  ptr->sho11_low_dens_eos();
  return;
}

void o2scl_eos_tov_interp_s12_low_dens_eos(void *vptr, void *ptr_model, bool external) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  std::string *model=(std::string *)ptr_model;
  ptr->s12_low_dens_eos(*model,external);
  return;
}

void o2scl_eos_tov_interp_gcp10_low_dens_eos(void *vptr, void *ptr_model, bool external) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  std::string *model=(std::string *)ptr_model;
  ptr->gcp10_low_dens_eos(*model,external);
  return;
}

void o2scl_eos_tov_interp_ngl13_low_dens_eos(void *vptr, double L, void *ptr_model, bool external) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  std::string *model=(std::string *)ptr_model;
  ptr->ngl13_low_dens_eos(L,*model,external);
  return;
}

void o2scl_eos_tov_interp_ngl13_low_dens_eos2(void *vptr, double S, double L, double nt, void *ptr_fname) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  std::string *fname=(std::string *)ptr_fname;
  ptr->ngl13_low_dens_eos2(S,L,nt,*fname);
  return;
}

void o2scl_eos_tov_interp_no_low_dens_eos(void *vptr) {
  eos_tov_interp *ptr=(eos_tov_interp *)vptr;
  ptr->no_low_dens_eos();
  return;
}

void *o2scl_create_tov_solve() {
  tov_solve *ptr=new tov_solve;
  return ptr;
}

void o2scl_free_tov_solve(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  delete ptr;
  return;
}

size_t o2scl_tov_solve_get_buffer_size(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->buffer_size;
}

void o2scl_tov_solve_set_buffer_size(void *vptr, size_t v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->buffer_size=v;
  return;
}

size_t o2scl_tov_solve_get_max_table_size(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->max_table_size;
}

void o2scl_tov_solve_set_max_table_size(void *vptr, size_t v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->max_table_size=v;
  return;
}

double o2scl_tov_solve_get_mass(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->mass;
}

void o2scl_tov_solve_set_mass(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->mass=v;
  return;
}

double o2scl_tov_solve_get_rad(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->rad;
}

void o2scl_tov_solve_set_rad(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->rad=v;
  return;
}

double o2scl_tov_solve_get_bmass(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->bmass;
}

void o2scl_tov_solve_set_bmass(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->bmass=v;
  return;
}

double o2scl_tov_solve_get_gpot(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->gpot;
}

void o2scl_tov_solve_set_gpot(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->gpot=v;
  return;
}

double o2scl_tov_solve_get_last_rjw(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->last_rjw;
}

void o2scl_tov_solve_set_last_rjw(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->last_rjw=v;
  return;
}

double o2scl_tov_solve_get_last_f(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->last_f;
}

void o2scl_tov_solve_set_last_f(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->last_f=v;
  return;
}

double o2scl_tov_solve_get_domega_rat(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->domega_rat;
}

void o2scl_tov_solve_set_domega_rat(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->domega_rat=v;
  return;
}

double o2scl_tov_solve_get_pcent_max(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->pcent_max;
}

void o2scl_tov_solve_set_pcent_max(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->pcent_max=v;
  return;
}

bool o2scl_tov_solve_get_reformat_results(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->reformat_results;
}

void o2scl_tov_solve_set_reformat_results(void *vptr, bool v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->reformat_results=v;
  return;
}

double o2scl_tov_solve_get_baryon_mass(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->baryon_mass;
}

void o2scl_tov_solve_set_baryon_mass(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->baryon_mass=v;
  return;
}

bool o2scl_tov_solve_get_ang_vel(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->ang_vel;
}

void o2scl_tov_solve_set_ang_vel(void *vptr, bool v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->ang_vel=v;
  return;
}

bool o2scl_tov_solve_get_gen_rel(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->gen_rel;
}

void o2scl_tov_solve_set_gen_rel(void *vptr, bool v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->gen_rel=v;
  return;
}

bool o2scl_tov_solve_get_calc_gpot(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->calc_gpot;
}

void o2scl_tov_solve_set_calc_gpot(void *vptr, bool v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->calc_gpot=v;
  return;
}

double o2scl_tov_solve_get_step_min(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->step_min;
}

void o2scl_tov_solve_set_step_min(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->step_min=v;
  return;
}

double o2scl_tov_solve_get_step_max(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->step_max;
}

void o2scl_tov_solve_set_step_max(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->step_max=v;
  return;
}

double o2scl_tov_solve_get_step_start(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->step_start;
}

void o2scl_tov_solve_set_step_start(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->step_start=v;
  return;
}

int o2scl_tov_solve_get_verbose(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->verbose;
}

void o2scl_tov_solve_set_verbose(void *vptr, int v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->verbose=v;
  return;
}

size_t o2scl_tov_solve_get_max_integ_steps(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->max_integ_steps;
}

void o2scl_tov_solve_set_max_integ_steps(void *vptr, size_t v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->max_integ_steps=v;
  return;
}

bool o2scl_tov_solve_get_err_nonconv(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->err_nonconv;
}

void o2scl_tov_solve_set_err_nonconv(void *vptr, bool v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->err_nonconv=v;
  return;
}

double o2scl_tov_solve_get_pmax_default(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->pmax_default;
}

void o2scl_tov_solve_set_pmax_default(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->pmax_default=v;
  return;
}

double o2scl_tov_solve_get_prbegin(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->prbegin;
}

void o2scl_tov_solve_set_prbegin(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->prbegin=v;
  return;
}

double o2scl_tov_solve_get_prend(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->prend;
}

void o2scl_tov_solve_set_prend(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->prend=v;
  return;
}

double o2scl_tov_solve_get_princ(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->princ;
}

void o2scl_tov_solve_set_princ(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->princ=v;
  return;
}

double o2scl_tov_solve_get_fixed_pr_guess(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->fixed_pr_guess;
}

void o2scl_tov_solve_set_fixed_pr_guess(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->fixed_pr_guess=v;
  return;
}

double o2scl_tov_solve_get_max_begin(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->max_begin;
}

void o2scl_tov_solve_set_max_begin(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->max_begin=v;
  return;
}

double o2scl_tov_solve_get_max_end(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->max_end;
}

void o2scl_tov_solve_set_max_end(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->max_end=v;
  return;
}

double o2scl_tov_solve_get_max_inc(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  return ptr->max_inc;
}

void o2scl_tov_solve_set_max_inc(void *vptr, double v) {
  tov_solve *ptr=(tov_solve *)vptr;
  ptr->max_inc=v;
  return;
}

void o2scl_tov_solve_set_eos(void *vptr, void *ptr_eos) {
  tov_solve *ptr=(tov_solve *)vptr;
  eos_tov *eos=(eos_tov *)ptr_eos;
  ptr->set_eos(*eos);
  return;
}

int o2scl_tov_solve_mvsr(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  int ret=ptr->mvsr();
  return ret;
}

int o2scl_tov_solve_fixed(void *vptr, double mass, double pmax) {
  tov_solve *ptr=(tov_solve *)vptr;
  int ret=ptr->fixed(mass,pmax);
  return ret;
}

int o2scl_tov_solve_fixed_pr(void *vptr, double pcent, double pmax) {
  tov_solve *ptr=(tov_solve *)vptr;
  int ret=ptr->fixed_pr(pcent,pmax);
  return ret;
}

int o2scl_tov_solve_max(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  int ret=ptr->max();
  return ret;
}

void *o2scl_tov_solve_get_results(void *vptr) {
  tov_solve *ptr=(tov_solve *)vptr;
  std::shared_ptr<table_units<> > *ret=new std::shared_ptr<table_units<> >;
  *ret=ptr->get_results();
  return ret;
}

void *o2scl_create_tov_love() {
  tov_love *ptr=new tov_love;
  return ptr;
}

void o2scl_free_tov_love(void *vptr) {
  tov_love *ptr=(tov_love *)vptr;
  delete ptr;
  return;
}

int o2scl_tov_love_get_show_ode(void *vptr) {
  tov_love *ptr=(tov_love *)vptr;
  return ptr->show_ode;
}

void o2scl_tov_love_set_show_ode(void *vptr, int v) {
  tov_love *ptr=(tov_love *)vptr;
  ptr->show_ode=v;
  return;
}

bool o2scl_tov_love_get_addl_testing(void *vptr) {
  tov_love *ptr=(tov_love *)vptr;
  return ptr->addl_testing;
}

void o2scl_tov_love_set_addl_testing(void *vptr, bool v) {
  tov_love *ptr=(tov_love *)vptr;
  ptr->addl_testing=v;
  return;
}

bool o2scl_tov_love_get_err_nonconv(void *vptr) {
  tov_love *ptr=(tov_love *)vptr;
  return ptr->err_nonconv;
}

void o2scl_tov_love_set_err_nonconv(void *vptr, bool v) {
  tov_love *ptr=(tov_love *)vptr;
  ptr->err_nonconv=v;
  return;
}

void *o2scl_tov_love_get_results(void *vptr) {
  tov_love *ptr=(tov_love *)vptr;
  return (void *)(&(ptr->results));
}

void o2scl_tov_love_set_results(void *vptr, void *p_v) {
  tov_love *ptr=(tov_love *)vptr;
  table_units<> *p_tsot=(table_units<> *)p_v;
  ptr->results=*(p_tsot);
  return;
}

double o2scl_tov_love_get_delta(void *vptr) {
  tov_love *ptr=(tov_love *)vptr;
  return ptr->delta;
}

void o2scl_tov_love_set_delta(void *vptr, double v) {
  tov_love *ptr=(tov_love *)vptr;
  ptr->delta=v;
  return;
}

double o2scl_tov_love_get_eps(void *vptr) {
  tov_love *ptr=(tov_love *)vptr;
  return ptr->eps;
}

void o2scl_tov_love_set_eps(void *vptr, double v) {
  tov_love *ptr=(tov_love *)vptr;
  ptr->eps=v;
  return;
}

void o2scl_tov_love_get_tab(void *vptr, void *p_v) {
  tov_love *ptr=(tov_love *)vptr;
  std::shared_ptr<table_units<> > *p_tgsp=(std::shared_ptr<table_units<> > *)p_v;
  *(p_tgsp)=ptr->tab;
  return;
}

void o2scl_tov_love_set_tab(void *vptr, void *p_v) {
  tov_love *ptr=(tov_love *)vptr;
  std::shared_ptr<table_units<> > *p_tssp=(std::shared_ptr<table_units<> > *)p_v;
  ptr->tab=*(p_tssp);
  return;
}

int o2scl_tov_love_calc_y(void *vptr, double *yR, double *beta, double *k2, double *lambda_km5, double *lambda_cgs, bool tabulate) {
  tov_love *ptr=(tov_love *)vptr;
  int ret=ptr->calc_y(*yR,*beta,*k2,*lambda_km5,*lambda_cgs,tabulate);
  return ret;
}

void o2scl_tov_love_add_disc(void *vptr, double rd) {
  tov_love *ptr=(tov_love *)vptr;
  ptr->add_disc(rd);
  return;
}

void o2scl_tov_love_clear_discs(void *vptr) {
  tov_love *ptr=(tov_love *)vptr;
  ptr->clear_discs();
  return;
}

int o2scl_tov_love_calc_H(void *vptr, double *yR, double *beta, double *k2, double *lambda_km5, double *lambda_cgs) {
  tov_love *ptr=(tov_love *)vptr;
  int ret=ptr->calc_H(*yR,*beta,*k2,*lambda_km5,*lambda_cgs);
  return ret;
}

void *o2scl_create_nstar_cold() {
  nstar_cold *ptr=new nstar_cold;
  return ptr;
}

void o2scl_free_nstar_cold(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  delete ptr;
  return;
}

double o2scl_nstar_cold_get_pressure_dec_nb(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->pressure_dec_nb;
}

void o2scl_nstar_cold_set_pressure_dec_nb(void *vptr, double v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->pressure_dec_nb=v;
  return;
}

double o2scl_nstar_cold_get_allow_urca_nb(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->allow_urca_nb;
}

void o2scl_nstar_cold_set_allow_urca_nb(void *vptr, double v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->allow_urca_nb=v;
  return;
}

double o2scl_nstar_cold_get_deny_urca_nb(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->deny_urca_nb;
}

void o2scl_nstar_cold_set_deny_urca_nb(void *vptr, double v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->deny_urca_nb=v;
  return;
}

double o2scl_nstar_cold_get_acausal_nb(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->acausal_nb;
}

void o2scl_nstar_cold_set_acausal_nb(void *vptr, double v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->acausal_nb=v;
  return;
}

double o2scl_nstar_cold_get_acausal_ed(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->acausal_ed;
}

void o2scl_nstar_cold_set_acausal_ed(void *vptr, double v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->acausal_ed=v;
  return;
}

double o2scl_nstar_cold_get_acausal_pr(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->acausal_pr;
}

void o2scl_nstar_cold_set_acausal_pr(void *vptr, double v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->acausal_pr=v;
  return;
}

void *o2scl_nstar_cold_get_def_tov(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return (void *)(&(ptr->def_tov));
}


bool o2scl_nstar_cold_get_eos_neg(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->eos_neg;
}

void o2scl_nstar_cold_set_eos_neg(void *vptr, bool v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->eos_neg=v;
  return;
}

int o2scl_nstar_cold_get_verbose(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->verbose;
}

void o2scl_nstar_cold_set_verbose(void *vptr, int v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->verbose=v;
  return;
}

double o2scl_nstar_cold_get_nb_start(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->nb_start;
}

void o2scl_nstar_cold_set_nb_start(void *vptr, double v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->nb_start=v;
  return;
}

double o2scl_nstar_cold_get_nb_end(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->nb_end;
}

void o2scl_nstar_cold_set_nb_end(void *vptr, double v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->nb_end=v;
  return;
}

double o2scl_nstar_cold_get_dnb(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->dnb;
}

void o2scl_nstar_cold_set_dnb(void *vptr, double v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->dnb=v;
  return;
}

size_t o2scl_nstar_cold_get_max_row(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->max_row;
}

void o2scl_nstar_cold_set_max_row(void *vptr, size_t v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->max_row=v;
  return;
}

bool o2scl_nstar_cold_get_remove_rows(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->remove_rows;
}

void o2scl_nstar_cold_set_remove_rows(void *vptr, bool v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->remove_rows=v;
  return;
}

bool o2scl_nstar_cold_get_include_muons(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->include_muons;
}

void o2scl_nstar_cold_set_include_muons(void *vptr, bool v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->include_muons=v;
  return;
}

bool o2scl_nstar_cold_get_err_nonconv(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  return ptr->err_nonconv;
}

void o2scl_nstar_cold_set_err_nonconv(void *vptr, bool v) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  ptr->err_nonconv=v;
  return;
}

void o2scl_nstar_cold_set_eos(void *vptr, void *ptr_eos) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  eos_had_base *eos=(eos_had_base *)ptr_eos;
  ptr->set_eos(*eos);
  return;
}

int o2scl_nstar_cold_calc_eos(void *vptr, double np_0) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  int ret=ptr->calc_eos(np_0);
  return ret;
}

int o2scl_nstar_cold_calc_nstar(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  int ret=ptr->calc_nstar();
  return ret;
}

int o2scl_nstar_cold_fixed(void *vptr, double target_mass) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  int ret=ptr->fixed(target_mass);
  return ret;
}

void *o2scl_nstar_cold_get_eos_results(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  std::shared_ptr<table_units<> > *ret=new std::shared_ptr<table_units<> >;
  *ret=ptr->get_eos_results();
  return ret;
}

void *o2scl_nstar_cold_get_tov_results(void *vptr) {
  nstar_cold *ptr=(nstar_cold *)vptr;
  std::shared_ptr<table_units<> > *ret=new std::shared_ptr<table_units<> >;
  *ret=ptr->get_tov_results();
  return ret;
}

void *o2scl_create_nucleus_rmf() {
  nucleus_rmf *ptr=new nucleus_rmf;
  return ptr;
}

void o2scl_free_nucleus_rmf(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  delete ptr;
  return;
}

double o2scl_nucleus_rmf_get_stens(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  return ptr->stens;
}

void o2scl_nucleus_rmf_set_stens(void *vptr, double v) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  ptr->stens=v;
  return;
}

double o2scl_nucleus_rmf_get_rnrp(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  return ptr->rnrp;
}

void o2scl_nucleus_rmf_set_rnrp(void *vptr, double v) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  ptr->rnrp=v;
  return;
}

double o2scl_nucleus_rmf_get_rnrms(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  return ptr->rnrms;
}

void o2scl_nucleus_rmf_set_rnrms(void *vptr, double v) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  ptr->rnrms=v;
  return;
}

double o2scl_nucleus_rmf_get_rprms(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  return ptr->rprms;
}

void o2scl_nucleus_rmf_set_rprms(void *vptr, double v) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  ptr->rprms=v;
  return;
}

double o2scl_nucleus_rmf_get_etot(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  return ptr->etot;
}

void o2scl_nucleus_rmf_set_etot(void *vptr, double v) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  ptr->etot=v;
  return;
}

double o2scl_nucleus_rmf_get_r_charge(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  return ptr->r_charge;
}

void o2scl_nucleus_rmf_set_r_charge(void *vptr, double v) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  ptr->r_charge=v;
  return;
}

double o2scl_nucleus_rmf_get_r_charge_cm(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  return ptr->r_charge_cm;
}

void o2scl_nucleus_rmf_set_r_charge_cm(void *vptr, double v) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  ptr->r_charge_cm=v;
  return;
}

int o2scl_nucleus_rmf_run_nucleus(void *vptr, int nucleus_Z, int nucleus_N, int unocc_Z, int unocc_N) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  int ret=ptr->run_nucleus(nucleus_Z,nucleus_N,unocc_Z,unocc_N);
  return ret;
}

void *o2scl_nucleus_rmf_get_profiles(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  std::shared_ptr<table_units<> > *ret=new std::shared_ptr<table_units<> >;
  *ret=ptr->get_profiles();
  return ret;
}

void *o2scl_nucleus_rmf_get_chden(void *vptr) {
  nucleus_rmf *ptr=(nucleus_rmf *)vptr;
  std::shared_ptr<table_units<> > *ret=new std::shared_ptr<table_units<> >;
  *ret=ptr->get_chden();
  return ret;
}

void *o2scl_create_nucmass_ldrop() {
  nucmass_ldrop *ptr=new nucmass_ldrop;
  return ptr;
}

void o2scl_free_nucmass_ldrop(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  delete ptr;
  return;
}

double o2scl_nucmass_ldrop_get_n1(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->n1;
}

void o2scl_nucmass_ldrop_set_n1(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->n1=v;
  return;
}

double o2scl_nucmass_ldrop_get_n0(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->n0;
}

void o2scl_nucmass_ldrop_set_n0(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->n0=v;
  return;
}

double o2scl_nucmass_ldrop_get_surften(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->surften;
}

void o2scl_nucmass_ldrop_set_surften(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->surften=v;
  return;
}

double o2scl_nucmass_ldrop_get_coul_coeff(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->coul_coeff;
}

void o2scl_nucmass_ldrop_set_coul_coeff(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->coul_coeff=v;
  return;
}

double o2scl_nucmass_ldrop_get_nn(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->nn;
}

void o2scl_nucmass_ldrop_set_nn(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->nn=v;
  return;
}

double o2scl_nucmass_ldrop_get_np(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->np;
}

void o2scl_nucmass_ldrop_set_np(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->np=v;
  return;
}

double o2scl_nucmass_ldrop_get_Rn(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->Rn;
}

void o2scl_nucmass_ldrop_set_Rn(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->Rn=v;
  return;
}

double o2scl_nucmass_ldrop_get_Rp(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->Rp;
}

void o2scl_nucmass_ldrop_set_Rp(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->Rp=v;
  return;
}

double o2scl_nucmass_ldrop_get_surf(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->surf;
}

void o2scl_nucmass_ldrop_set_surf(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->surf=v;
  return;
}

double o2scl_nucmass_ldrop_get_bulk(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->bulk;
}

void o2scl_nucmass_ldrop_set_bulk(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->bulk=v;
  return;
}

double o2scl_nucmass_ldrop_get_coul(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->coul;
}

void o2scl_nucmass_ldrop_set_coul(void *vptr, double v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->coul=v;
  return;
}

bool o2scl_nucmass_ldrop_get_large_vals_unphys(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return ptr->large_vals_unphys;
}

void o2scl_nucmass_ldrop_set_large_vals_unphys(void *vptr, bool v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  ptr->large_vals_unphys=v;
  return;
}

void *o2scl_nucmass_ldrop_get_def_had_eos(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return (void *)(&(ptr->def_had_eos));
}

void o2scl_nucmass_ldrop_set_def_had_eos(void *vptr, void *p_v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  eos_had_rmf *p_tsot=(eos_had_rmf *)p_v;
  ptr->def_had_eos=*(p_tsot);
  return;
}

void *o2scl_nucmass_ldrop_get_def_neutron(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return (void *)(&(ptr->def_neutron));
}

void o2scl_nucmass_ldrop_set_def_neutron(void *vptr, void *p_v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  fermion *p_tsot=(fermion *)p_v;
  ptr->def_neutron=*(p_tsot);
  return;
}

void *o2scl_nucmass_ldrop_get_def_proton(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return (void *)(&(ptr->def_proton));
}

void o2scl_nucmass_ldrop_set_def_proton(void *vptr, void *p_v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  fermion *p_tsot=(fermion *)p_v;
  ptr->def_proton=*(p_tsot);
  return;
}

void *o2scl_nucmass_ldrop_get_th(void *vptr) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  return (void *)(&(ptr->th));
}

void o2scl_nucmass_ldrop_set_th(void *vptr, void *p_v) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  thermo *p_tsot=(thermo *)p_v;
  ptr->th=*(p_tsot);
  return;
}

double o2scl_nucmass_ldrop_mass_excess_d(void *vptr, double Z, double N) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  double ret=ptr->mass_excess_d(Z,N);
  return ret;
}

double o2scl_nucmass_ldrop_mass_excess(void *vptr, int Z, int N) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  double ret=ptr->mass_excess(Z,N);
  return ret;
}

double o2scl_nucmass_ldrop_drip_binding_energy_d(void *vptr, double Z, double N, double npout, double nnout, double chi, double dim, double T) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  double ret=ptr->drip_binding_energy_d(Z,N,npout,nnout,chi,dim,T);
  return ret;
}

void o2scl_nucmass_ldrop_set_n_and_p(void *vptr, void *ptr_un, void *ptr_up) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  fermion *un=(fermion *)ptr_un;
  fermion *up=(fermion *)ptr_up;
  ptr->set_n_and_p(*un,*up);
  return;
}

int o2scl_nucmass_ldrop_set_eos_had_temp_base(void *vptr, void *ptr_uhe) {
  nucmass_ldrop *ptr=(nucmass_ldrop *)vptr;
  eos_had_temp_base *uhe=(eos_had_temp_base *)ptr_uhe;
  int ret=ptr->set_eos_had_temp_base(*uhe);
  return ret;
}

void *o2scl_create_nucmass_ldrop_skin() {
  nucmass_ldrop_skin *ptr=new nucmass_ldrop_skin;
  return ptr;
}

void o2scl_free_nucmass_ldrop_skin(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  delete ptr;
  return;
}

bool o2scl_nucmass_ldrop_skin_get_full_surface(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->full_surface;
}

void o2scl_nucmass_ldrop_skin_set_full_surface(void *vptr, bool v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->full_surface=v;
  return;
}

bool o2scl_nucmass_ldrop_skin_get_new_skin_mode(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->new_skin_mode;
}

void o2scl_nucmass_ldrop_skin_set_new_skin_mode(void *vptr, bool v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->new_skin_mode=v;
  return;
}

double o2scl_nucmass_ldrop_skin_get_doi(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->doi;
}

void o2scl_nucmass_ldrop_skin_set_doi(void *vptr, double v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->doi=v;
  return;
}

double o2scl_nucmass_ldrop_skin_get_ss(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->ss;
}

void o2scl_nucmass_ldrop_skin_set_ss(void *vptr, double v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->ss=v;
  return;
}

double o2scl_nucmass_ldrop_skin_get_pp(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->pp;
}

void o2scl_nucmass_ldrop_skin_set_pp(void *vptr, double v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->pp=v;
  return;
}

double o2scl_nucmass_ldrop_skin_get_a0(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->a0;
}

void o2scl_nucmass_ldrop_skin_set_a0(void *vptr, double v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->a0=v;
  return;
}

double o2scl_nucmass_ldrop_skin_get_a2(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->a2;
}

void o2scl_nucmass_ldrop_skin_set_a2(void *vptr, double v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->a2=v;
  return;
}

double o2scl_nucmass_ldrop_skin_get_a4(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->a4;
}

void o2scl_nucmass_ldrop_skin_set_a4(void *vptr, double v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->a4=v;
  return;
}

bool o2scl_nucmass_ldrop_skin_get_rel_vacuum(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->rel_vacuum;
}

void o2scl_nucmass_ldrop_skin_set_rel_vacuum(void *vptr, bool v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->rel_vacuum=v;
  return;
}

double o2scl_nucmass_ldrop_skin_get_Tchalf(void *vptr) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  return ptr->Tchalf;
}

void o2scl_nucmass_ldrop_skin_set_Tchalf(void *vptr, double v) {
  nucmass_ldrop_skin *ptr=(nucmass_ldrop_skin *)vptr;
  ptr->Tchalf=v;
  return;
}

void *o2scl_create_nucmass_ldrop_pair() {
  nucmass_ldrop_pair *ptr=new nucmass_ldrop_pair;
  return ptr;
}

void o2scl_free_nucmass_ldrop_pair(void *vptr) {
  nucmass_ldrop_pair *ptr=(nucmass_ldrop_pair *)vptr;
  delete ptr;
  return;
}

double o2scl_nucmass_ldrop_pair_get_Epair(void *vptr) {
  nucmass_ldrop_pair *ptr=(nucmass_ldrop_pair *)vptr;
  return ptr->Epair;
}

void o2scl_nucmass_ldrop_pair_set_Epair(void *vptr, double v) {
  nucmass_ldrop_pair *ptr=(nucmass_ldrop_pair *)vptr;
  ptr->Epair=v;
  return;
}

double o2scl_nucmass_ldrop_pair_get_pair(void *vptr) {
  nucmass_ldrop_pair *ptr=(nucmass_ldrop_pair *)vptr;
  return ptr->pair;
}

void o2scl_nucmass_ldrop_pair_set_pair(void *vptr, double v) {
  nucmass_ldrop_pair *ptr=(nucmass_ldrop_pair *)vptr;
  ptr->pair=v;
  return;
}

void *o2scl_create_nucleus_bin() {
  nucleus_bin *ptr=new nucleus_bin;
  return ptr;
}

void o2scl_free_nucleus_bin(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  delete ptr;
  return;
}

void *o2scl_nucleus_bin_get_ame16(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ame16));
}

void o2scl_nucleus_bin_set_ame16(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ame *p_tsot=(nucmass_ame *)p_v;
  ptr->ame16=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ame20exp(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ame20exp));
}

void o2scl_nucleus_bin_set_ame20exp(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ame *p_tsot=(nucmass_ame *)p_v;
  ptr->ame20exp=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ame20round(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ame20round));
}

void o2scl_nucleus_bin_set_ame20round(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ame *p_tsot=(nucmass_ame *)p_v;
  ptr->ame20round=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ame95rmd(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ame95rmd));
}

void o2scl_nucleus_bin_set_ame95rmd(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ame *p_tsot=(nucmass_ame *)p_v;
  ptr->ame95rmd=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ame03round(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ame03round));
}

void o2scl_nucleus_bin_set_ame03round(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ame *p_tsot=(nucmass_ame *)p_v;
  ptr->ame03round=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ame03(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ame03));
}

void o2scl_nucleus_bin_set_ame03(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ame *p_tsot=(nucmass_ame *)p_v;
  ptr->ame03=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ame95exp(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ame95exp));
}

void o2scl_nucleus_bin_set_ame95exp(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ame *p_tsot=(nucmass_ame *)p_v;
  ptr->ame95exp=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ame12(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ame12));
}

void o2scl_nucleus_bin_set_ame12(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ame *p_tsot=(nucmass_ame *)p_v;
  ptr->ame12=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ddme2(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ddme2));
}

void o2scl_nucleus_bin_set_ddme2(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->ddme2=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ddmed(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ddmed));
}

void o2scl_nucleus_bin_set_ddmed(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->ddmed=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_ddpc1(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->ddpc1));
}

void o2scl_nucleus_bin_set_ddpc1(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->ddpc1=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_nl3s(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->nl3s));
}

void o2scl_nucleus_bin_set_nl3s(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->nl3s=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_sly4(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->sly4));
}

void o2scl_nucleus_bin_set_sly4(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->sly4=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_skms(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->skms));
}

void o2scl_nucleus_bin_set_skms(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->skms=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_skp(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->skp));
}

void o2scl_nucleus_bin_set_skp(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->skp=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_sv_min(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->sv_min));
}

void o2scl_nucleus_bin_set_sv_min(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->sv_min=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_unedf0(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->unedf0));
}

void o2scl_nucleus_bin_set_unedf0(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->unedf0=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_unedf1(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->unedf1));
}

void o2scl_nucleus_bin_set_unedf1(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_gen *p_tsot=(nucmass_gen *)p_v;
  ptr->unedf1=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_m95(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->m95));
}

void o2scl_nucleus_bin_set_m95(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_mnmsk *p_tsot=(nucmass_mnmsk *)p_v;
  ptr->m95=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_m16(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->m16));
}

void o2scl_nucleus_bin_set_m16(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_mnmsk *p_tsot=(nucmass_mnmsk *)p_v;
  ptr->m16=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_kt(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->kt));
}

void o2scl_nucleus_bin_set_kt(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ktuy *p_tsot=(nucmass_ktuy *)p_v;
  ptr->kt=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_kt2(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->kt2));
}

void o2scl_nucleus_bin_set_kt2(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_ktuy *p_tsot=(nucmass_ktuy *)p_v;
  ptr->kt2=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_wlw1(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->wlw1));
}

void o2scl_nucleus_bin_set_wlw1(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_wlw *p_tsot=(nucmass_wlw *)p_v;
  ptr->wlw1=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_wlw2(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->wlw2));
}

void o2scl_nucleus_bin_set_wlw2(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_wlw *p_tsot=(nucmass_wlw *)p_v;
  ptr->wlw2=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_wlw3(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->wlw3));
}

void o2scl_nucleus_bin_set_wlw3(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_wlw *p_tsot=(nucmass_wlw *)p_v;
  ptr->wlw3=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_wlw4(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->wlw4));
}

void o2scl_nucleus_bin_set_wlw4(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_wlw *p_tsot=(nucmass_wlw *)p_v;
  ptr->wlw4=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_wlw5(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->wlw5));
}

void o2scl_nucleus_bin_set_wlw5(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_wlw *p_tsot=(nucmass_wlw *)p_v;
  ptr->wlw5=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_sdnp1(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->sdnp1));
}

void o2scl_nucleus_bin_set_sdnp1(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_sdnp *p_tsot=(nucmass_sdnp *)p_v;
  ptr->sdnp1=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_sdnp2(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->sdnp2));
}

void o2scl_nucleus_bin_set_sdnp2(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_sdnp *p_tsot=(nucmass_sdnp *)p_v;
  ptr->sdnp2=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_sdnp3(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->sdnp3));
}

void o2scl_nucleus_bin_set_sdnp3(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_sdnp *p_tsot=(nucmass_sdnp *)p_v;
  ptr->sdnp3=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_dz(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->dz));
}

void o2scl_nucleus_bin_set_dz(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_dz_table *p_tsot=(nucmass_dz_table *)p_v;
  ptr->dz=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb2(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb2));
}

void o2scl_nucleus_bin_set_hfb2(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb *p_tsot=(nucmass_hfb *)p_v;
  ptr->hfb2=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb8(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb8));
}

void o2scl_nucleus_bin_set_hfb8(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb *p_tsot=(nucmass_hfb *)p_v;
  ptr->hfb8=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb14(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb14));
}

void o2scl_nucleus_bin_set_hfb14(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb *p_tsot=(nucmass_hfb *)p_v;
  ptr->hfb14=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb14_v0(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb14_v0));
}

void o2scl_nucleus_bin_set_hfb14_v0(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb *p_tsot=(nucmass_hfb *)p_v;
  ptr->hfb14_v0=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb17(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb17));
}

void o2scl_nucleus_bin_set_hfb17(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb_sp *p_tsot=(nucmass_hfb_sp *)p_v;
  ptr->hfb17=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb21(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb21));
}

void o2scl_nucleus_bin_set_hfb21(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb_sp *p_tsot=(nucmass_hfb_sp *)p_v;
  ptr->hfb21=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb22(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb22));
}

void o2scl_nucleus_bin_set_hfb22(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb_sp *p_tsot=(nucmass_hfb_sp *)p_v;
  ptr->hfb22=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb23(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb23));
}

void o2scl_nucleus_bin_set_hfb23(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb_sp *p_tsot=(nucmass_hfb_sp *)p_v;
  ptr->hfb23=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb24(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb24));
}

void o2scl_nucleus_bin_set_hfb24(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb_sp *p_tsot=(nucmass_hfb_sp *)p_v;
  ptr->hfb24=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb25(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb25));
}

void o2scl_nucleus_bin_set_hfb25(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb_sp *p_tsot=(nucmass_hfb_sp *)p_v;
  ptr->hfb25=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb26(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb26));
}

void o2scl_nucleus_bin_set_hfb26(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb_sp *p_tsot=(nucmass_hfb_sp *)p_v;
  ptr->hfb26=*(p_tsot);
  return;
}

void *o2scl_nucleus_bin_get_hfb27(void *vptr) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  return (void *)(&(ptr->hfb27));
}

void o2scl_nucleus_bin_set_hfb27(void *vptr, void *p_v) {
  nucleus_bin *ptr=(nucleus_bin *)vptr;
  nucmass_hfb_sp *p_tsot=(nucmass_hfb_sp *)p_v;
  ptr->hfb27=*(p_tsot);
  return;
}

void o2scl_skyrme_load_wrapper(void *ptr_sk, void *ptr_model, bool external, int verbose) {
  eos_had_skyrme *sk=(eos_had_skyrme *)ptr_sk;
  std::string *model=(std::string *)ptr_model;
  skyrme_load(*sk,*model,external,verbose);
  return;
}

void o2scl_rmf_load_wrapper(void *ptr_rmf, void *ptr_model, bool external) {
  eos_had_rmf *rmf=(eos_had_rmf *)ptr_rmf;
  std::string *model=(std::string *)ptr_model;
  rmf_load(*rmf,*model,external);
  return;
}

