/* -------------------------------------------------------------------
  
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
#ifndef O2SCL_PART_PYTHON_H
#define O2SCL_PART_PYTHON_H

/** \file part_python.h
    \brief File for python interface
*/

extern "C" {

  void o2scl_set_err_hnd_gsl();
  
  void *o2scl_create_part(double *&g, double *&m, double *&ms,
			  double *&mu, double *&nu,
			  double *&ed, double *&pr,
			  double *&en, bool *&inc_rest_mass,
			  bool *&non_interacting);
  void *o2scl_create_fermion(double *&g, double *&m, double *&ms,
			     double *&mu, double *&nu,
			     double *&ed, double *&pr,
			     double *&en, bool *&inc_rest_mass,
			     bool *&non_interacting, double *&kf,
			     double *&del);

  void o2scl_free_part(void *vp);
  void o2scl_free_fermion(void *vp);

  void *o2scl_create_fermion_rel();
  void *o2scl_create_fermion_nonrel();
  void *o2scl_create_fermion_deriv_rel();
  void *o2scl_create_fermion_deriv_nr();
  void *o2scl_create_boson_rel();
  void *o2scl_create_classical_thermo();
  void *o2scl_create_classical_deriv_thermo();
  void *o2scl_create_fermion_mag_zerot();
  
  void o2scl_free_fermion_rel(void *vp);
  void o2scl_free_fermion_nonrel(void *vp);
  void o2scl_free_fermion_deriv_rel(void *vp);
  void o2scl_free_fermion_deriv_nr(void *vp);
  void o2scl_free_boson_rel(void *vp);
  void o2scl_free_classical_thermo(void *vp);
  void o2scl_free_classical_deriv_thermo(void *vp);
  void o2scl_free_fermion_mag_zerot(void *vp);
  
  void o2scl_fermion_rel_calc_density(void *frp, void *fp, double T);
  void o2scl_fermion_rel_calc_mu(void *frp, void *fp, double T);
  void o2scl_classical_calc_density(void *frp, void *fp, double T);
  void o2scl_classical_calc_mu(void *frp, void *fp, double T);
  void o2scl_fermion_nonrel_calc_density(void *frp, void *fp, double T);
  void o2scl_fermion_nonrel_calc_mu(void *frp, void *fp, double T);
  void o2scl_fermion_deriv_rel_calc_density(void *frp, void *fp, double T);
  void o2scl_fermion_deriv_rel_calc_mu(void *frp, void *fp, double T);
  
}

#endif
