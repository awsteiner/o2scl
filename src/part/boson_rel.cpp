/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/boson_rel.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

//--------------------------------------------
// boson_rel class

boson_rel::boson_rel() {
}

boson_rel::~boson_rel() {
}

void boson_rel::set_inte(inte<> &l_nit, inte<> &l_dit) {
}

void boson_rel::calc_mu(boson &b, double temper) {
}

void boson_rel::calc_max_density(boson &b, double temper) {
}

void boson_rel::nu_from_n(boson &b, double temper) {
}

void boson_rel::calc_density(boson &b, double temper) {
}

double boson_rel::deg_density_fun(double k, boson &b, double T) {
}
  
double boson_rel::deg_energy_fun(double k, boson &b, double T) {
  
double boson_rel::deg_entropy_fun(double k, boson &b, double T) {
  
double boson_rel::density_fun(double u, boson &b, double T) {

double boson_rel::energy_fun(double u, boson &b, double T) {

double boson_rel::entropy_fun(double u, boson &b, double T) {
int boson_rel::solve_fun(size_t nv, const ubvector &x, ubvector &y,
                         double density, boson &b, double T) {

void boson_rel::pair_mu(boson &b, double temper) {
void boson_rel::pair_density(boson &b, double temper) {
int boson_rel::pair_density_fun(size_t nv,
                                const ubvector &x, ubvector &y,
                                double density, boson &b, double T) {
