/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015, Andrew W. Steiner
  
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

#include <o2scl/test_mgr.h>
#include <o2scl/tov_love.h>
#include <o2scl/tov_solve.h>
#include <o2scl/eos_had_apr.h>
#include <o2scl/nstar_cold.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

int main(void) {

  cout.setf(ios::scientific);

  test_mgr t;
  t.set_output_level(2);

  double schwarz_km=o2scl_cgs::schwarzchild_radius/1.0e5;
  
  eos_had_apr apr;
  apr.pion=1;
  
  nstar_cold nst;
  nst.set_eos(apr);
  nst.calc_eos();
  
  o2_shared_ptr<table_units<> >::type te=nst.get_eos_results();

  eos_tov_interp eti;
  eti.default_low_dens_eos();
  eti.read_table(*te,"ed","pr","nb");

  tov_solve ts;
  ts.verbose=0;
  ts.set_eos(eti);
  ts.calc_gpot=true;
  ts.ang_vel=true;
  ts.fixed(1.4);

  o2_shared_ptr<table_units<> >::type profile=ts.get_results();

  // I'm not sure why cs2 isn't computed properly, so
  // we have to recompute it here
  profile->delete_column("cs2");
  profile->deriv("ed","pr","cs2");

  double I=ts.domega_rat*pow(ts.rad,4.0)/3.0/schwarz_km;
  cout << "Radius: " << ts.rad << " km" << endl;
  cout << "Moment of inertia: " << I << " Msun*km^2" << endl;

  // Compute dimensionless tidal deformability from
  // Yagi and Yunes (2013) correlation
  double l_ibar=log(I/pow(schwarz_km/2.0,2.0)/pow(1.4,3.0));
  double lbar=exp(-30.5395+38.3931*l_ibar-16.3071*pow(l_ibar,2.0)+
		  3.36972*pow(l_ibar,3.0)-0.26105*pow(l_ibar,4.0));

  // Direct calculation
  tov_love tl;
  tl.tab=profile;
  double yR, beta, k2, lambda_km5, lambda_cgs;
  tl.calc_y(yR,beta,k2,lambda_km5,lambda_cgs);
  
  cout << "Dimensionless tidal deformability (YY13 correlation): "
       << lbar << endl;
  cout << "Dimensionless tidal deformability (direct calculation): " 
       << lambda_km5/pow(schwarz_km/2.0,5.0)/pow(1.4,5.0) << endl;

  t.report();
  
  return 0;
}


