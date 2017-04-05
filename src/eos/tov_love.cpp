/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2017, Andrew W. Steiner
  
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
#include "tov_love.h"

#include <o2scl/lib_settings.h>

using namespace std;
using namespace o2scl;

typedef boost::numeric::ublas::vector<double> ubvector;

tov_love::tov_love() {
  oisp=&def_ois;
  schwarz_km=o2scl_cgs::schwarzchild_radius/1.0e5;
  eps=0.02;
  /*
    oisp->gsl_astp.con.a_dydt=1.0;
    oisp->gsl_astp.con.eps_rel=1.0e-10;
    oisp->gsl_astp.con.eps_abs=1.0e-10;
    oisp->ntrial*=10;
  */
}

double tov_love::eval_k2(double beta, double yR) {
  return (8.0/5.0*pow(beta,5.0)*pow(1.0-2.0*beta,2.0)*
	  (2.0-yR+2.0*beta*(yR-1.0))/
	  (2.0*beta*(6.0-3.0*yR+3.0*beta*(5.0*yR-8.0))+
	   4.0*beta*beta*beta*(13.0-11.0*yR+beta*(3.0*yR-2.0)+2.0*
			       beta*beta*(1.0+yR))+
	   3.0*pow(1.0-2.0*beta,2.0)*(2.0-yR+2.0*beta*(yR-1.0))*
	   log(1.0-2.0*beta)));
}

int tov_love::y_derivs(double r, size_t nv, const ubvector &vals,
			  ubvector &ders) {
    
  double ed=tab->interp("r",r,"ed");
  double pr=tab->interp("r",r,"pr");
  double cs2=tab->interp("r",r,"cs2");
  double gm=tab->interp("r",r,"gm");
  double pi=o2scl_const::pi;

  if (r==0.0) {
    ders[0]=0.0;
  } else {
    double elam=1.0/(1.0-schwarz_km*gm/r);
    double nup=schwarz_km*elam*(gm+4.0*pi*pr*r*r*r)/r/r;
      
    double Q=2.0*pi*schwarz_km*elam*(5.0*ed+9.0*pr+(ed+pr)/cs2)-
      6.0*elam/r/r-nup*nup;
      
    double y=vals[0];
      
    ders[0]=(-r*r*Q-y*elam*(1.0+2.0*pi*schwarz_km*r*r*(pr-ed))-y*y)/r;

  } 

  if (!gsl_finite(ders[0])) {
    return 1;
    std::cout << "Derivative not finite." << std::endl;
    std::cout << "r,ed,pr,cs2,gm: " 
	      << r << " " << ed << " " << pr << " " << cs2 << " " 
	      << gm << std::endl;
    std::cout << "y': " << ders[0] << std::endl;
    exit(-1);
  }
  return 0;
}

int tov_love::H_derivs(double r, size_t nv, const ubvector &vals,
			  ubvector &ders) {
    
  double ed=tab->interp("r",r,"ed");
  double pr=tab->interp("r",r,"pr");
  double cs2=tab->interp("r",r,"cs2");
  double gm=tab->interp("r",r,"gm");
  double pi=o2scl_const::pi;

  double dHdr=vals[1];
  double H=vals[0];

  if (r==0.0) {
    ders[0]=2.0*r;
    ders[1]=2.0;
    return 0;
  }

  double fact=1.0-schwarz_km*gm/r;
  ders[0]=dHdr;
  /*
    Units of G (schwarz_km/2.0) is [km/Msun].
    Units of d^2H/dr^2 are H/r/r, so to convert pressure to 1/r/r
    we have to multiply by G. Same to convert radius*pressure to 1/r.
  */
  ders[1]=2.0/fact*H*(-schwarz_km*pi*(5.0*ed+9.0*pr+(ed+pr)/cs2)+3.0/r/r
		      +2.0/fact*pow(schwarz_km/2.0*gm/r/r+
				    2.0*pi*pr*r*schwarz_km,2.0))+
    2.0*dHdr/r/fact*(-1.0+schwarz_km/2.0*gm/r+pi*r*r*(ed-pr)*schwarz_km);

  return 0;
}

void tov_love::calc_y(double &yR, double &beta, double &k2, 
		      double &lambda_km5, double &lambda_cgs,
		      bool tabulate) {
  
  size_t count;

  double R=tab->max("r");
  double gm=tab->max("gm");

  double r=eps, h=1.0e-1;
  ubvector y(1), dydx_out(1), yerr(1);
  y[0]=2.0;

  if (tabulate) {
    results.clear();
    results.line_of_names("r y dydr");
    results.set_unit("r","km");
    results.set_unit("dydr","1/km");
  }

  ode_funct od=std::bind
    (std::mem_fn<int(double,size_t,const ubvector &,ubvector &)>
     (&tov_love::y_derivs),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::placeholders::_4);

    ubvector yout(1);
    oisp->solve_final_value(eps,R,h,1,y,yout,od);
    y[0]=yout[0];

  yR=y[0];
  beta=schwarz_km/2.0*gm/R;

  k2=eval_k2(beta,yR);
  
  lambda_km5=2.0/3.0*k2*pow(R,5.0);

  // First compute in Msun*km^2*s^2
  lambda_cgs=2.0/3.0*k2*pow(R,5.0)/
    pow(o2scl_mks::speed_of_light/1.0e3,2.0)/schwarz_km;
  // Convert to g*cm^2*s^2
  lambda_cgs*=o2scl_cgs::solar_mass*1.0e10;
  
  if (tabulate) {
    results.add_constant("k2",k2);
    results.add_constant("lambda_cgs",lambda_cgs);
    results.add_constant("lambda_km5",lambda_km5);
  }

  //std::cout << "Love number=" << lambda_cgs << std::endl;

  return;
}

void tov_love::calc_H(double &yR, double &beta, double &k2, 
			 double &lambda_km5, double &lambda_cgs) {

  // Try again, with Hinderer eqn.

  double R=tab->max("r");
  double gm=tab->max("gm");

  double r=eps, h=1.0e-1;
  ubvector y(2), dydx_out(2), yerr(2);

  y[0]=r*r;
  y[1]=2.0*r;

  ode_funct od2=std::bind
    (std::mem_fn<int(double,size_t,const ubvector &,ubvector &)>
     (&tov_love::H_derivs),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::placeholders::_4);
  
  ubvector yout(2);

  size_t N=10;
  for(size_t i=0;i<N;i++) {
    double r2=eps+(R-eps)/((double)N)*((double)(i+1));
    oisp->solve_final_value(r,r2,h,2,y,yout,od2);
    y[0]=yout[0];
    y[1]=yout[1];
    r=r2;
  }
    
  yR=R*y[1]/y[0];

  beta=schwarz_km/2.0*gm/R;

  k2=eval_k2(beta,yR);
    
  lambda_km5=2.0/3.0*k2*pow(R,5.0);
  
  // First compute in Msun*km^2*s^2
  lambda_cgs=2.0/3.0*k2*pow(R,5.0)/
    pow(o2scl_mks::speed_of_light/1.0e3,2.0)/schwarz_km;
  // Convert to g*cm^2*s^2
  lambda_cgs*=o2scl_cgs::solar_mass*1.0e10;

  return;
}
