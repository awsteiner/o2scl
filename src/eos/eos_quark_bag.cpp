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

#include <o2scl/eos_quark_bag.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

eos_quark_bag::eos_quark_bag() {
  bag_constant=200.0/hc_mev_fm;
}

int eos_quark_bag::calc_p(quark &u, quark &d, quark &s, thermo &th) {

  if (u.mu>u.m) u.kf=sqrt(u.mu*u.mu-u.m*u.m);
  else u.kf=0.0;
  if (d.mu>d.m) d.kf=sqrt(d.mu*d.mu-d.m*d.m);
  else d.kf=0.0;
  if (s.mu>s.m) s.kf=sqrt(s.mu*s.mu-s.m*s.m);
  else s.kf=0.0;

  u.ms=u.m;
  d.ms=d.m;
  s.ms=s.m;
  
  fet.calc_mu_zerot(u);
  fet.calc_mu_zerot(d);
  fet.calc_mu_zerot(s);
  
  th.ed=u.ed+d.ed+s.ed+bag_constant;
  th.pr=u.pr+d.pr+s.pr-bag_constant;
  th.en=0.0;
  
  return 0;
}

int eos_quark_bag::calc_e(quark &u, quark &d, quark &s, thermo &th) {

  fet.kf_from_density(u);
  fet.kf_from_density(d);
  fet.kf_from_density(s);

  u.mu=sqrt(u.kf*u.kf+u.m*u.m);
  d.mu=sqrt(d.kf*d.kf+d.m*d.m);
  s.mu=sqrt(s.kf*s.kf+s.m*s.m);

  u.ms=u.m;
  d.ms=d.m;
  s.ms=s.m;

  fet.calc_density_zerot(u);
  fet.calc_density_zerot(d);
  fet.calc_density_zerot(s);

  th.ed=u.ed+d.ed+s.ed+bag_constant;
  th.pr=u.pr+d.pr+s.pr-bag_constant;
  th.en=0.0;
  
  return 0;
}

int eos_quark_bag::calc_temp_p(quark& u, quark& d, quark& s, 
			      const double temper, thermo &th) {

  if (temper<=0.0) {
    return calc_p(u,d,s,th);
  }

  u.nu=u.mu;
  d.nu=d.mu;
  s.nu=s.mu;
  u.ms=u.m;
  d.ms=d.m;
  s.ms=s.m;
  
  fet.pair_mu(u,temper);
  fet.pair_mu(d,temper);
  fet.pair_mu(s,temper);
  
  fet.kf_from_density(u);
  fet.kf_from_density(d);
  fet.kf_from_density(s);
  
  th.ed=u.ed+d.ed+s.ed+bag_constant;
  th.pr=u.pr+d.pr+s.pr-bag_constant;
  th.en=(th.ed+th.pr-u.n*u.mu-d.n*d.mu-s.n*s.mu)/temper;

  return 0;
}

int eos_quark_bag::calc_temp_e(quark& u, quark& d, quark& s, 
			      const double temper, thermo &th) {

  u.ms=u.m;
  d.ms=d.m;
  s.ms=s.m;
  
  fet.pair_density(u,temper);
  fet.pair_density(d,temper);
  fet.pair_density(s,temper);

  fet.kf_from_density(u);
  fet.kf_from_density(d);
  fet.kf_from_density(s);
  
  th.ed=u.ed+d.ed+s.ed+bag_constant;
  th.pr=u.pr+d.pr+s.pr-bag_constant;
  th.en=(th.ed+th.pr-u.n*u.mu-d.n*d.mu-s.n*s.mu)/temper;

  return 0;
}

