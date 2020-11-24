/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/part_python.h>
#include <o2scl/fermion_rel.h>
#include <o2scl/fermion_nonrel.h>
#include <o2scl/fermion_deriv_nr.h>
#include <o2scl/fermion_deriv_rel.h>
#include <o2scl/boson_rel.h>
#include <o2scl/classical.h>
#include <o2scl/classical_deriv.h>
#include <o2scl/fermion_mag_zerot.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

void *o2scl_create_part(double *&g, double *&m, double *&ms,
			double *&mu, double *&nu,
			double *&ed, double *&pr,
			double *&en, bool *&inc_rest_mass,
			bool *&non_interacting) {
  part *pp=new part;
  g=&pp->g;
  ms=&pp->ms;
  m=&pp->m;
  mu=&pp->mu;
  nu=&pp->nu;
  ed=&pp->ed;
  pr=&pp->pr;
  en=&pp->en;
  inc_rest_mass=&pp->inc_rest_mass;
  non_interacting=&pp->non_interacting;
  return pp;
}

void *o2scl_create_fermion() {
  fermion *fr=new fermion;
  return fr;
}

void *o2scl_create_fermion_rel() {
  fermion_rel *fr=new fermion_rel;
  return fr;
}

void *o2scl_create_fermion_nonrel() {
  fermion_nonrel *fr=new fermion_nonrel;
  return fr;
}

void *o2scl_create_fermion_deriv_nr() {
  fermion_deriv_nr *fr=new fermion_deriv_nr;
  return fr;
}

void *o2scl_create_fermion_deriv_rel() {
  fermion_deriv_rel *fr=new fermion_deriv_rel;
  return fr;
}

void *o2scl_create_boson_rel() {
  boson_rel *fr=new boson_rel;
  return fr;
}

void *o2scl_create_classical_thermo() {
  classical_thermo *fr=new classical_thermo;
  return fr;
}

void *o2scl_create_classical_deriv_thermo() {
  classical_deriv_thermo *fr=new classical_deriv_thermo;
  return fr;
}

void *o2scl_create_fermion_mag_zerot() {
  fermion_mag_zerot *fr=new fermion_mag_zerot;
  return fr;
}

void o2scl_free_part(void *vp) {
  part *fr=(part *)vp;
  delete fr;
  return;
}

void o2scl_free_fermion(void *vp) {
  fermion *fr=(fermion *)vp;
  delete fr;
  return;
}

void o2scl_free_fermion_rel(void *vp) {
  fermion_rel *fr=(fermion_rel *)vp;
  delete fr;
  return;
}

void o2scl_free_fermion_nonrel(void *vp) {
  fermion_nonrel *fr=(fermion_nonrel *)vp;
  delete fr;
}

void o2scl_free_fermion_deriv_nr(void *vp) {
  fermion_deriv_nr *fr=(fermion_deriv_nr *)vp;
  delete fr;
}

void o2scl_free_fermion_deriv_rel(void *vp) {
  fermion_deriv_rel *fr=(fermion_deriv_rel *)vp;
  delete fr;
}

void o2scl_free_boson_rel(void *vp) {
  boson_rel *fr=(boson_rel *)vp;
  delete fr;
}

void o2scl_free_classical_thermo(void *vp) {
  classical_thermo *fr=(classical_thermo *)vp;
  delete fr;
}

void o2scl_free_classical_deriv_thermo(void *vp) {
  classical_deriv_thermo *fr=(classical_deriv_thermo *)vp;
  delete fr;
}

void o2scl_free_fermion_mag_zerot(void *vp) {
  fermion_mag_zerot *fr=(fermion_mag_zerot *)vp;
  delete fr;
}

void o2scl_fermion_rel_calc_density(void *frp, void *fp, double T) {

  fermion_rel *fr=(fermion_rel *)frp;
  fermion *f=(fermion *)fp;
  fr->calc_density(*f,T);

  return;
}

void o2scl_fermion_rel_calc_mu(void *frp, void *fp, double T) {

  fermion_rel *fr=(fermion_rel *)frp;
  fermion *f=(fermion *)fp;
  fr->calc_mu(*f,T);

  return;
}

void o2scl_classical_calc_density(void *frp, void *fp, double T) {

  classical_thermo *fr=(classical_thermo *)frp;
  part *f=(part *)fp;
  fr->calc_density(*f,T);

  return;
}

void o2scl_classical_calc_mu(void *frp, void *fp, double T) {

  classical_thermo *fr=(classical_thermo *)frp;
  part *f=(part *)fp;
  fr->calc_mu(*f,T);

  return;
}

void o2scl_fermion_nonrel_calc_density(void *frp, void *fp, double T) {

  fermion_nonrel *fr=(fermion_nonrel *)frp;
  fermion *f=(fermion *)fp;
  fr->calc_density(*f,T);

  return;
}

void o2scl_fermion_nonrel_calc_mu(void *frp, void *fp, double T) {

  fermion_nonrel *fr=(fermion_nonrel *)frp;
  fermion *f=(fermion *)fp;
  fr->calc_mu(*f,T);

  return;
}

void o2scl_fermion_deriv_rel_calc_density(void *frp, void *fp, double T) {

  fermion_deriv_rel *fr=(fermion_deriv_rel *)frp;
  fermion_deriv *f=(fermion_deriv *)fp;
  fr->calc_density(*f,T);

  return;
}

void o2scl_fermion_deriv_rel_calc_mu(void *frp, void *fp, double T) {

  fermion_deriv_rel *fr=(fermion_deriv_rel *)frp;
  fermion_deriv *f=(fermion_deriv *)fp;
  fr->calc_mu(*f,T);

  return;
}

