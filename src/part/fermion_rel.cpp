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

#include <o2scl/fermion_rel.h>
#include <o2scl/root_cern.h>
#include <o2scl/root_bkt_cern.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/inte_qag_gsl.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

void *o2scl_create_fermion_rel() {
  fermion_rel *fr=new fermion_rel;
  return fr;
}

void o2scl_free_fermion_rel(void *vp) {
  fermion_rel *fr=(fermion_rel *)vp;
  delete fr;
  return;
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
