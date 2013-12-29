/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#include <o2scl/ddc_eos.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

ddc_eos::ddc_eos() {
  mnuc=939.0/hc_mev_fm;
  ms=550.0/197.33;
  mw=783.0/197.33;
  mr=763.0/197.33;
  Gs=10.72854;
  Gw=13.29015;
  Gr=7.32196;
  as=1.365469;
  aw=1.402488;
  ar=0.515;
  bs=0.226061;
  bw=0.172577;
  cs=0.409704;
  cw=0.344293;
  ds=0.901995;
  dw=0.983955;
  rho0=0.153;

}

int ddc_eos::calc_eq_e(fermion &ne, fermion &pr, double sig, double ome, 
		       double lrho, double &f1, double &f2, double &f3, 
		       thermo &lth) {

  ne.non_interacting=false;
  pr.non_interacting=false;

  ne.ms=ne.m-Gs*sig;
  pr.ms=pr.m-Gs*sig;
  
  double rho=ne.n+pr.n;
  double rhoI=(pr.n-ne.n)/2.0;
  double x=rho/rho0;
  double xpds2=(x+ds)*(x+ds);
  double Gsr=Gs*as*(1.0+bs*xpds2)/(1.0+cs*xpds2);
  double xpdw2=(x+dw)*(x+dw);
  double Gwr=Gw*aw*(1.0+bw*xpdw2)/(1.0+cs*xpdw2);
  double Grr=Gr*exp(-ar*(x-1.0));
  
  double dgsdr=Gs*(as*bs*2.0*(x+ds)/(1.0+cs*xpds2)-
		as*cs*2.0*(x+ds)*(1.0+bs*xpds2)/pow((1.0+cs*xpds2),2.0))/rho0;
  double dgwdr=Gw*(aw*bw*2.0*(x+dw)/(1.0+cw*xpdw2)-
		aw*cw*2.0*(x+dw)*(1.0+bw*xpdw2)/pow((1.0+cw*xpdw2),2.0))/rho0;
  double dgrdr=-ar*Gr*exp(-ar*(x-1.0))/rho0;
  
  fzt.kf_from_density(ne);
  fzt.kf_from_density(pr);

  ne.nu=sqrt(ne.kf*ne.kf+ne.ms*ne.ms);
  pr.nu=sqrt(ne.kf*ne.kf+ne.ms*ne.ms);
  
  // We don't record error values, since these functions usually
  // always succeed. The Fermi momentum is computed automatically as
  // well.
  fzt.calc_density_zerot(ne);
  fzt.calc_density_zerot(pr);
  double nsn=1.0/ne.ms*(ne.ed-3.0*ne.pr);
  double nsp=1.0/pr.ms*(pr.ed-3.0*pr.pr);
  
  ne.mu=ne.nu+(Gw+dgwdr*rho)*ome+(-Gr/2.0+dgrdr*rhoI)*lrho-
    dgsdr*(nsn+nsp)*sig;
  pr.mu=pr.nu+(Gw+dgwdr*rho)*ome+(Gr/2.0+dgrdr*rhoI)*lrho-
    dgsdr*(nsn+nsp)*sig;
  
  f1=ms*ms*sig-Gs*(nsn+nsp);
  f2=mw*mw*ome-Gw*rho;
  f3=mr*mr*lrho-Gr*rhoI;

  lth.ed=ne.ed+pr.ed+0.5*ms*ms*sig*sig+0.5*mw*mw*ome*ome+0.5*mr*mr*lrho*lrho;
  lth.pr=ne.pr+pr.pr-0.5*ms*ms*sig*sig*(1.0+2.0*rho/Gs*dgsdr)+
    0.5*(mw*mw*ome*ome*(1.0+2.0*rho/Gw*dgwdr)+mr*mr*lrho*lrho*
	 (1.0+2.0*rho/Gr*dgrdr));

  return success;
}
