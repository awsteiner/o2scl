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

void o2scl_fermion_density(void *vp, double m, double g,
			   double T, double n,
			   double *mu, double *ed, double *pr,
			   double *en) {
  fermion_rel *fr=(fermion_rel *)vp;
  fermion f(m,g);
  f.n=n;
  fr->calc_density(f,T);
  *mu=f.mu;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_fermion_int_density(void *vp, double m, double ms, double g,
			       double T, double n,
			       double *nu, double *ed, double *pr,
			       double *en) {
  fermion_rel *fr=(fermion_rel *)vp;
  fermion f(m,g);
  f.non_interacting=false;
  f.n=n;
  fr->calc_density(f,T);
  *nu=f.nu;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_fermion_mu(void *vp, double m, double g,
		      double T, double mu,
		      double *n, double *ed, double *pr,
		      double *en) {
  fermion_rel *fr=(fermion_rel *)vp;
  fermion f(m,g);
  f.mu=mu;
  fr->calc_mu(f,T);
  *n=f.n;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_fermion_int_mu(void *vp, double m, double ms, double g,
			  double T, double nu,
			  double *n, double *ed, double *pr,
			  double *en) {
  fermion_rel *fr=(fermion_rel *)vp;
  fermion f(m,g);
  f.non_interacting=false;
  f.nu=nu;
  fr->calc_mu(f,T);
  *n=f.n;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_fermion_nonrel_density(void *vp, double m, double g,
			   double T, double n,
			   double *mu, double *ed, double *pr,
			   double *en) {
  fermion_nonrel *fr=(fermion_nonrel *)vp;
  fermion f(m,g);
  f.n=n;
  fr->calc_density(f,T);
  *mu=f.mu;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_fermion_nonrel_int_density(void *vp, double m, double ms, double g,
			       double T, double n,
			       double *nu, double *ed, double *pr,
			       double *en) {
  fermion_nonrel *fr=(fermion_nonrel *)vp;
  fermion f(m,g);
  f.non_interacting=false;
  f.n=n;
  fr->calc_density(f,T);
  *nu=f.nu;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_fermion_nonrel_mu(void *vp, double m, double g,
		      double T, double mu,
		      double *n, double *ed, double *pr,
		      double *en) {
  fermion_nonrel *fr=(fermion_nonrel *)vp;
  fermion f(m,g);
  f.mu=mu;
  fr->calc_mu(f,T);
  *n=f.n;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_fermion_nonrel_int_mu(void *vp, double m, double ms, double g,
			  double T, double nu,
			  double *n, double *ed, double *pr,
			  double *en) {
  fermion_nonrel *fr=(fermion_nonrel *)vp;
  fermion f(m,g);
  f.non_interacting=false;
  f.nu=nu;
  fr->calc_mu(f,T);
  *n=f.n;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_classical_density(void *vp, double m, double g,
			   double T, double n,
			   double *mu, double *ed, double *pr,
			   double *en) {
  classical_thermo *cl=(classical_thermo *)vp;
  fermion f(m,g);
  f.n=n;
  cl->calc_density(f,T);
  *mu=f.mu;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_classical_mu(void *vp, double m, double g,
		      double T, double mu,
		      double *n, double *ed, double *pr,
		      double *en) {
  classical_thermo *cl=(classical_thermo *)vp;
  fermion f(m,g);
  f.mu=mu;
  cl->calc_mu(f,T);
  *n=f.n;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;

  return;
}

void o2scl_fermion_deriv_mu(void *vp, double m, double g,
			    double T, double mu,
			    double *n, double *ed, double *pr,
			    double *en, double *dndT,
			    double *dsdT, double *dndmu) {
  fermion_deriv_rel *fr=(fermion_deriv_rel *)vp;
  fermion_deriv f(m,g);
  f.mu=mu;
  fr->calc_mu(f,T);
  *n=f.n;
  *ed=f.ed;
  *pr=f.pr;
  *en=f.en;
  *dndT=f.dndT;
  *dsdT=f.dsdT;
  *dndmu=f.dndmu;

  return;
}
