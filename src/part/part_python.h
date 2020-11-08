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
  
  void *o2scl_create_fermion_rel();
  void *o2scl_create_fermion_nonrel();
  void *o2scl_create_fermion_deriv_rel();
  void *o2scl_create_fermion_deriv_nr();
  void *o2scl_create_boson_rel();
  void *o2scl_create_classical_thermo();
  void *o2scl_create_classical_deriv_thermo();
  void *o2scl_create_fermion_mag_zerot();
  void o2scl_free_fermion_rel(void *vp);
  void o2scl_fermion_density
  (void *vp, double m, double g, double T, double n,
   double *mu, double *ed, double *pr, double *en);
  void o2scl_fermion_int_density
  (void *vp, double m, double g, double T, double n,
   double *nu, double *ed, double *pr, double *en);
  void o2scl_fermion_mu
  (void *vp, double m, double ms, double g, double T, double mu,
   double *n, double *ed, double *pr, double *en);
  void o2scl_fermion_int_mu
  (void *vp, double m, double ms, double g, double T, double nu,
   double *n, double *ed, double *pr, double *en);
  
}

#endif
