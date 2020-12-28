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

void *o2scl_create_fermion();

void o2scl_free_fermion(void *vp);

double o2scl_fermion_get_kf(void *vp);

void o2scl_fermion_set_kf(void *vp, double v);

double o2scl_fermion_get_del(void *vp);

void o2scl_fermion_set_del(void *vp, double v);

void *o2scl_create_fermion_rel();

void o2scl_free_fermion_rel(void *vp);

int o2scl_fermion_rel_calc_density(void *vptr, void *ptr_f, double T);

void o2scl_fermion_rel_calc_mu(void *vptr, void *ptr_f, double T);

void *o2scl_create_fermion_nonrel();

void o2scl_free_fermion_nonrel(void *vp);

void *o2scl_create_fermion_deriv_nr();

void o2scl_free_fermion_deriv_nr(void *vp);

void *o2scl_create_fermion_deriv_rel();

void o2scl_free_fermion_deriv_rel(void *vp);

void *o2scl_create_boson_rel();

void o2scl_free_boson_rel(void *vp);

void *o2scl_create_classical_thermo();

void o2scl_free_classical_thermo(void *vp);

void *o2scl_create_classical_deriv_thermo();

void o2scl_free_classical_deriv_thermo(void *vp);

void *o2scl_create_fermion_mag_zerot();

void o2scl_free_fermion_mag_zerot(void *vp);

}
