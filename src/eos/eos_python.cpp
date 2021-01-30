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

#include <o2scl/eos_python.h>
#include <o2scl/part.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

void *o2scl_create_eos_base() {
  eos_base *ptr=new eos_base;
  return ptr;
}

void o2scl_free_eos_base(void *vptr) {
  eos_base *ptr=(eos_base *)vptr;
  delete ptr;
}

void o2scl_eos_base_get_def_thermo(void *vptr, void *p_v) {
  eos_base *ptr=(eos_base *)vptr;
  o2scl::thermo *p_t=(o2scl::thermo *)p_v;
  *(p_t)=ptr->def_thermo;
  return;
}

void o2scl_eos_base_set_def_thermo(void *vptr, void *p_v) {
  eos_base *ptr=(eos_base *)vptr;
  o2scl::thermo *p_t=(o2scl::thermo *)p_v;
  ptr->def_thermo=*(p_t);
  return;
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

void o2scl_eos_had_base_get_def_neutron(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_t=(o2scl::fermion *)p_v;
  *(p_t)=ptr->def_neutron;
  return;
}

void o2scl_eos_had_base_set_def_neutron(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_t=(o2scl::fermion *)p_v;
  ptr->def_neutron=*(p_t);
  return;
}

void o2scl_eos_had_base_get_def_proton(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_t=(o2scl::fermion *)p_v;
  *(p_t)=ptr->def_proton;
  return;
}

void o2scl_eos_had_base_set_def_proton(void *vptr, void *p_v) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  o2scl::fermion *p_t=(o2scl::fermion *)p_v;
  ptr->def_proton=*(p_t);
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

double o2scl_eos_had_base_fcomp_err(void *vptr, double nb, double delta, void *ptr_unc) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double *unc=(double *)ptr_unc;
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

double o2scl_eos_had_base_fesym_err(void *vptr, double nb, double delta, void *ptr_unc) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double *unc=(double *)ptr_unc;
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

double o2scl_eos_had_base_fn0(void *vptr, double delta, void *ptr_leoa) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double *leoa=(double *)ptr_leoa;
  double ret=ptr->fn0(delta,*leoa);
  return ret;
}

void o2scl_eos_had_base_f_number_suscept(void *vptr, double mun, double mup, void *ptr_dPdnn, void *ptr_dPdnp, void *ptr_dPdpp) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double *dPdnn=(double *)ptr_dPdnn;
  double *dPdnp=(double *)ptr_dPdnp;
  double *dPdpp=(double *)ptr_dPdpp;
  ptr->f_number_suscept(mun,mup,*dPdnn,*dPdnp,*dPdpp);
  return;
}

void o2scl_eos_had_base_f_inv_number_suscept(void *vptr, double mun, double mup, void *ptr_dednn, void *ptr_dednp, void *ptr_dedpp) {
  eos_had_base *ptr=(eos_had_base *)vptr;
  double *dednn=(double *)ptr_dednn;
  double *dednp=(double *)ptr_dednp;
  double *dedpp=(double *)ptr_dedpp;
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

void *o2scl_create_eos_had_skyrme() {
  eos_had_skyrme *ptr=new eos_had_skyrme;
  return ptr;
}

void o2scl_free_eos_had_skyrme(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  delete ptr;
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

const char *o2scl_eos_had_skyrme_get_reference(void *vptr) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  python_temp_string=ptr->reference;
  return python_temp_string.c_str();
}

void o2scl_eos_had_skyrme_set_reference(void *vptr, void *p_v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  std::string *p_t=(std::string *)p_v;
  ptr->reference=*(p_t);
  return;
}

void o2scl_eos_had_skyrme_get_nrfd(void *vptr, void *p_v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  o2scl::fermion_deriv_nr *p_t=(o2scl::fermion_deriv_nr *)p_v;
  *(p_t)=ptr->nrfd;
  return;
}

void o2scl_eos_had_skyrme_set_nrfd(void *vptr, void *p_v) {
  eos_had_skyrme *ptr=(eos_had_skyrme *)vptr;
  o2scl::fermion_deriv_nr *p_t=(o2scl::fermion_deriv_nr *)p_v;
  ptr->nrfd=*(p_t);
  return;
}

void *o2scl_create_eos_had_apr() {
  eos_had_apr *ptr=new eos_had_apr;
  return ptr;
}

void o2scl_free_eos_had_apr(void *vptr) {
  eos_had_apr *ptr=(eos_had_apr *)vptr;
  delete ptr;
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

bool o2scl_eos_had_rmf_get_err_nonconv(void *vptr) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  return ptr->err_nonconv;
}

void o2scl_eos_had_rmf_set_err_nonconv(void *vptr, bool v) {
  eos_had_rmf *ptr=(eos_had_rmf *)vptr;
  ptr->err_nonconv=v;
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

void *o2scl_create_eos_quark() {
  eos_quark *ptr=new eos_quark;
  return ptr;
}

void o2scl_free_eos_quark(void *vptr) {
  eos_quark *ptr=(eos_quark *)vptr;
  delete ptr;
}

void *o2scl_create_eos_quark_bag() {
  eos_quark_bag *ptr=new eos_quark_bag;
  return ptr;
}

void o2scl_free_eos_quark_bag(void *vptr) {
  eos_quark_bag *ptr=(eos_quark_bag *)vptr;
  delete ptr;
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

bool o2scl_eos_quark_njl_get_fromqq(void *vptr) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  return ptr->fromqq;
}

void o2scl_eos_quark_njl_set_fromqq(void *vptr, bool v) {
  eos_quark_njl *ptr=(eos_quark_njl *)vptr;
  ptr->fromqq=v;
  return;
}

void o2scl_skyrme_load_wrapper(void *ptr_sk, char *model, bool external, int verbose) {
  eos_had_skyrme *sk=(eos_had_skyrme *)ptr_sk;
  skyrme_load(*sk,model,external,verbose);
  return;
}

void o2scl_rmf_load_wrapper(void *ptr_rmf, char *model, bool external) {
  eos_had_rmf *rmf=(eos_had_rmf *)ptr_rmf;
  rmf_load(*rmf,model,external);
  return;
}

