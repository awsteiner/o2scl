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

#include <o2scl/part_python.h>

using namespace std;
using namespace o2scl;

void *o2scl_create_part() {
  part *ptr=new part;
  return ptr;
}

void o2scl_free_part(void *vptr) {
  part *ptr=(part *)vptr;
  delete ptr;
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

void *o2scl_create_fermion() {
  fermion *ptr=new fermion;
  return ptr;
}

void o2scl_free_fermion(void *vptr) {
  fermion *ptr=(fermion *)vptr;
  delete ptr;
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

void *o2scl_create_fermion_rel() {
  fermion_rel *ptr=new fermion_rel;
  return ptr;
}

void o2scl_free_fermion_rel(void *vptr) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  delete ptr;
}

int o2scl_fermion_rel_calc_density(void *vptr, void *ptr_f, double T) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  fermion *f=(fermion *)ptr_f;
  int ret=ptr->calc_density(*f,T);
  return ret;
}

void o2scl_fermion_rel_calc_mu(void *vptr, void *ptr_f, double T) {
  fermion_rel *ptr=(fermion_rel *)vptr;
  fermion *f=(fermion *)ptr_f;
  ptr->calc_mu(*f,T);
  return;
}

void *o2scl_create_fermion_nonrel() {
  fermion_nonrel *ptr=new fermion_nonrel;
  return ptr;
}

void o2scl_free_fermion_nonrel(void *vptr) {
  fermion_nonrel *ptr=(fermion_nonrel *)vptr;
  delete ptr;
}

void *o2scl_create_fermion_deriv_nr() {
  fermion_deriv_nr *ptr=new fermion_deriv_nr;
  return ptr;
}

void o2scl_free_fermion_deriv_nr(void *vptr) {
  fermion_deriv_nr *ptr=(fermion_deriv_nr *)vptr;
  delete ptr;
}

void *o2scl_create_fermion_deriv_rel() {
  fermion_deriv_rel *ptr=new fermion_deriv_rel;
  return ptr;
}

void o2scl_free_fermion_deriv_rel(void *vptr) {
  fermion_deriv_rel *ptr=(fermion_deriv_rel *)vptr;
  delete ptr;
}

void *o2scl_create_boson_rel() {
  boson_rel *ptr=new boson_rel;
  return ptr;
}

void o2scl_free_boson_rel(void *vptr) {
  boson_rel *ptr=(boson_rel *)vptr;
  delete ptr;
}

void *o2scl_create_classical_thermo() {
  classical_thermo *ptr=new classical_thermo;
  return ptr;
}

void o2scl_free_classical_thermo(void *vptr) {
  classical_thermo *ptr=(classical_thermo *)vptr;
  delete ptr;
}

void *o2scl_create_classical_deriv_thermo() {
  classical_deriv_thermo *ptr=new classical_deriv_thermo;
  return ptr;
}

void o2scl_free_classical_deriv_thermo(void *vptr) {
  classical_deriv_thermo *ptr=(classical_deriv_thermo *)vptr;
  delete ptr;
}

void *o2scl_create_fermion_mag_zerot() {
  fermion_mag_zerot *ptr=new fermion_mag_zerot;
  return ptr;
}

void o2scl_free_fermion_mag_zerot(void *vptr) {
  fermion_mag_zerot *ptr=(fermion_mag_zerot *)vptr;
  delete ptr;
}

