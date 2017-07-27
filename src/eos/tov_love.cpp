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
  delta=1.0e-4;
}

double tov_love::eval_k2(double beta, double yR) {
  /*
    This is a slightly reformatted but equivalent expression:
    
    double k2_new=8.0/5.0*pow(beta,5.0)*pow(1.0-2.0*beta,2.0)*
    (2.0-yR+2.0*beta*(yR-1.0))/
    (2.0*beta*(6.0-3.0*yR+3.0*beta*(5.0*yR-8.0)+
    2.0*beta*beta*(13.0-11.0*yR+beta*(3.0*yR-2.0)+
    2.0*beta*beta*(1.0+yR)))+
    3.0*pow(1.0-2.0*beta,2.0)*(2.0-yR+2.0*beta*(yR-1.0))*
    log(1.0-2.0*beta));
  */
  double k2=8.0/5.0*pow(beta,5.0)*pow(1.0-2.0*beta,2.0)*
    (2.0-yR+2.0*beta*(yR-1.0))/
    (2.0*beta*(6.0-3.0*yR+3.0*beta*(5.0*yR-8.0))+
     4.0*beta*beta*beta*(13.0-11.0*yR+beta*(3.0*yR-2.0)+2.0*
			 beta*beta*(1.0+yR))+
     3.0*pow(1.0-2.0*beta,2.0)*(2.0-yR+2.0*beta*(yR-1.0))*
     log(1.0-2.0*beta));
  return k2;
}

int tov_love::y_derivs(double r, size_t nv, const std::vector<double> &vals,
			  std::vector<double> &ders) {
    
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

int tov_love::H_derivs(double r, size_t nv, const std::vector<double> &vals,
			  std::vector<double> &ders) {
    
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
    Units of G (schwarz_km/2.0) is [km/Msun]. Units of d^2H/dr^2 are
    H/r/r, so to convert pressure to 1/r/r we have to multiply by G.
    The same method is used to convert radius*pressure to 1/r.
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
  
  if (disc.size()>0 && tabulate) {
    O2SCL_ERR2("Function tov_love::calc_y() does not yet handle ",
	       "discontinuities when tabulate is true.",
	       o2scl::exc_eunimpl);
  }

  size_t count;

  double R=tab->get_constant("rad");
  double gm=tab->get_constant("mass");

  double r=eps, h=1.0e-1;
  std::vector<double> y(1), dydx_out(1), yerr(1);
  y[0]=2.0;

  // Storage for 'tabulate' option
  std::vector<double> rt;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  ubmatrix yt, dy, ye;
  size_t n_sol=10000;

  // If 'tabulate' is true, allocate space
  if (tabulate) {
    rt.resize(n_sol);
    yt.resize(n_sol,1);
    dy.resize(n_sol,1);
    ye.resize(n_sol,1);
  }

  ode_funct2 od=std::bind
    (std::mem_fn<int(double,size_t,const std::vector<double> &,
		     std::vector<double> &)>
     (&tov_love::y_derivs),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::placeholders::_4);
  
  if (tabulate) {
    
    yt(0,0)=2.0;
    oisp->solve_store(eps,R,h,1,n_sol,rt,yt,ye,dy,od,0);
    yR=yt(n_sol-1,0);
    
  } else {

    if (disc.size()>0 && eps>disc[0]) {
      O2SCL_ERR2("Discontinuity is at smaller radius than eps ",
		 "in tov_love::calc_y().",o2scl::exc_einval);
    }
    std::vector<double> yout(1);

    // Loop over intervals between discontinuities and the
    // r=0 and r=R boundaries
    double x0=eps, x1;

    for(size_t j=0;j<disc.size()+1;j++) {

      // Set the inner and outer radial boundaries for this interval
      if (j!=0) {
	x0=disc[j-1]+delta;
      } 
      if (j==disc.size()) {
	x1=R;
      } else {
	// Make sure the discontinuity is not at the neutron
	// star radius
	if (disc[j]>=R) {
	  O2SCL_ERR2("Discontinuity is at or past full radius ",
		     "in tov_love::calc_y().",o2scl::exc_einval);
	}
	x1=disc[j]-delta;
      }

      // Solve the ODE in this interval
      oisp->solve_final_value(x0,x1,h,1,y,yout,od);

      // Add the correction at the discontinuity
      if (j!=disc.size()) {
	y[0]=yout[0]+(tab->interp("r",disc[j]+delta,"ed")-
		      tab->interp("r",disc[j]-delta,"ed"))/
	  (tab->interp("r",disc[j],"gm")+4.0*o2scl_const::pi+
	   disc[j]*disc[j]*disc[j]*tab->interp("r",disc[j],"pr"));
      } else {
	y[0]=yout[0];
      }

      // Proceed to the next interval
    }

    // Set the final value of yR at the star's radius
    yR=y[0];
  }

  beta=schwarz_km/2.0*gm/R;

  k2=eval_k2(beta,yR);
  
  lambda_km5=2.0/3.0*k2*pow(R,5.0);

  // First compute in Msun*km^2*s^2
  lambda_cgs=2.0/3.0*k2*pow(R,5.0)/
    pow(o2scl_mks::speed_of_light/1.0e3,2.0)/schwarz_km;
  // Convert to g*cm^2*s^2
  lambda_cgs*=o2scl_cgs::solar_mass*1.0e10;
  
  if (tabulate) {
    results.clear();

    results.add_constant("k2",k2);
    results.add_constant("lambda_cgs",lambda_cgs);
    results.add_constant("lambda_km5",lambda_km5);

    results.line_of_names("r y dydr ye ed pr cs2 gm");
    results.set_unit("r","km");
    results.set_unit("dydr","1/km");
    results.set_unit("ye","km");
    results.set_unit("ed","Msun/km^3");
    results.set_unit("pr","Msun/km^3");
    results.set_unit("gm","Msun");
    for(size_t j=0;j<n_sol;j++) {
      double line[8]={rt[j],yt(j,0),dy(j,0),ye(j,0),
		      tab->interp("r",rt[j],"ed"),
		      tab->interp("r",rt[j],"pr"),
		      tab->interp("r",rt[j],"cs2"),
		      tab->interp("r",rt[j],"gm")};
      results.line_of_data(8,line);
    }
  }

  return;
}

void tov_love::calc_H(double &yR, double &beta, double &k2, 
			 double &lambda_km5, double &lambda_cgs) {

  if (disc.size()>0) {
    O2SCL_ERR2("Function tov_love::calc_H() does not yet handle ",
	       "discontinuities.",o2scl::exc_eunimpl);
  }
  
  double R=tab->max("r");
  double gm=tab->max("gm");

  double r=eps, h=1.0e-1;
  std::vector<double> y(2), dydx_out(2), yerr(2);

  y[0]=r*r;
  y[1]=2.0*r;

  ode_funct2 od2=std::bind
    (std::mem_fn<int(double,size_t,const std::vector<double> &,
		     std::vector<double> &)>
     (&tov_love::H_derivs),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::placeholders::_4);
  
  std::vector<double> yout(2);

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
