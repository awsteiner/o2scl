/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2021, Andrew W. Steiner
  
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
  err_nonconv=true;
  addl_testing=false;
  show_ode=0;
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

  if (show_ode>0) {
    std::cout << "r,ed,pr,cs2,gm: " 
	      << r << " " << ed << " " << pr << " " << cs2 << " " 
	      << gm << std::endl;
    std::cout << "y, y': " << vals[0] << " " << ders[0] << std::endl;
  }
  
  if (!gsl_finite(ders[0])) {
    return 1;
    /*
      std::cout << "Derivative not finite." << std::endl;
      std::cout << "r,ed,pr,cs2,gm: " 
      << r << " " << ed << " " << pr << " " << cs2 << " " 
      << gm << std::endl;
      std::cout << "y': " << ders[0] << std::endl;
      //exit(-1);
      */
  }
  return 0;
}

int tov_love::H_derivs(double r, size_t nv, const std::vector<double> &vals,
			  std::vector<double> &ders) {

  tab->is_valid();
  
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

int tov_love::calc_y(double &yR, double &beta, double &k2, 
		      double &lambda_km5, double &lambda_cgs,
		      bool tabulate) {
  
  tab->is_valid();

  if (disc.size()>0 && eps>disc[0]) {
    O2SCL_ERR2("Discontinuity is at smaller radius than eps ",
	       "in tov_love::calc_y().",o2scl::exc_einval);
  }
  
  size_t count;

  double R=tab->get_constant("rad");
  double gm=tab->get_constant("mass");

  double r=eps, h=1.0e-1;
  std::vector<double> y(1), dydx_out(1), yerr(1);
  y[0]=2.0;

  ode_funct2 od=std::bind
    (std::mem_fn<int(double,size_t,const std::vector<double> &,
		     std::vector<double> &)>
     (&tov_love::y_derivs),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::placeholders::_4);

  if (tabulate) {
    
    // Storage
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    std::vector< std::vector<double> > rt(disc.size()+1);
    std::vector<ubmatrix> yt(disc.size()+1);
    std::vector<ubmatrix> dy(disc.size()+1);
    std::vector<ubmatrix> ye(disc.size()+1);
    std::vector<size_t> n_sol(disc.size()+1);
    
    // Loop over intervals between discontinuities and the
    // r=0 and r=R boundaries
    double x0=eps, x1;

    for(size_t j=0;j<disc.size()+1;j++) {

      // This will be changed by solve_store() below to reflect
      // the size of the final solution
      n_sol[j]=10000;
      
      // Allocate space ahead of time, so we can fill the next section
      // with the previous value from this section
      if (j==0) {
	rt[j].resize(10000);
	yt[j].resize(10000,1);
	dy[j].resize(10000,1);
	ye[j].resize(10000,1);
      }
      if (j<disc.size()) {
	n_sol[j+1]=10000;
	rt[j+1].resize(10000);
	yt[j+1].resize(10000,1);
	dy[j+1].resize(10000,1);
	ye[j+1].resize(10000,1);
      }

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
      
      // Set up initial condition near r=0
      if (j==0) {
	yt[0](0,0)=2.0;
      }

      // Solve the ODE in this interval
      if (x0>x1) {
	O2SCL_CONV2_RET("Discontinuities too close to resolve ",
			" in tov_love::calc_y().",
			o2scl::exc_etol,err_nonconv);
      }
      
      int ois_ret=oisp->solve_store(x0,x1,h,1,n_sol[j],rt[j],yt[j],
				    ye[j],dy[j],od,0);
      if (ois_ret!=0) {
	O2SCL_CONV2_RET("ODE function solve_store() failed ",
			" in tov_love::calc_y().",
			o2scl::exc_efailed,err_nonconv);
      }
      
      // Add the correction at the discontinuity
      if (j!=disc.size()) {
	yt[j+1](0,0)=yt[j](n_sol[j]-1,0)+
	  (tab->interp("r",disc[j]+delta,"ed")-
	   tab->interp("r",disc[j]-delta,"ed"))/
	  (tab->interp("r",disc[j],"gm")/
	   (4.0*o2scl_const::pi*disc[j]*disc[j]*disc[j])+
	   tab->interp("r",disc[j],"pr"));
	/*
	  std::cout << endl;
	  std::cout << "Here3: "
	  << tab->interp("r",disc[j]+delta,"ed") << " "
	  << tab->interp("r",disc[j]-delta,"ed") << " "
	  << tab->interp("r",disc[j]+delta,"ed") << " "
	  << disc[j] << std::endl;
	  std::cout << (tab->interp("r",disc[j]+delta,"ed")-
	  tab->interp("r",disc[j]-delta,"ed"))/
	  (tab->interp("r",disc[j],"gm")/(4.0*o2scl_const::pi*
	  disc[j]*disc[j]*disc[j])+
	  tab->interp("r",disc[j],"pr")) << std::endl;
	  std::cout << tab->get_unit("ed") << " "
	  << tab->get_unit("pr") << std::endl;
	*/
      }
      
    }

    // Final value of y at r=R
    yR=yt[disc.size()](n_sol[disc.size()]-1,0);

    /*
      for(size_t j=0;j<disc.size()+1;j++) {
      std::cout << "Here: " << j << " "
      << yt[j](0,0) << " " << yt[j](n_sol[j]-1,0) << std::endl;
      }
      std::cout << "Here2: " << yR << std::endl;
    */
    
    results.clear();
    results.line_of_names("r y dydr ye ed pr cs2 gm");
    results.set_unit("r","km");
    results.set_unit("dydr","1/km");
    results.set_unit("ye","km");
    results.set_unit("ed","Msun/km^3");
    results.set_unit("pr","Msun/km^3");
    results.set_unit("gm","Msun");
    for(size_t j=0;j<disc.size()+1;j++) {
      for(size_t k=0;k<n_sol[j];k++) {
	double line[8]={rt[j][k],yt[j](k,0),dy[j](k,0),ye[j](k,0),
			tab->interp("r",rt[j][k],"ed"),
			tab->interp("r",rt[j][k],"pr"),
			tab->interp("r",rt[j][k],"cs2"),
			tab->interp("r",rt[j][k],"gm")};
	results.line_of_data(8,line);
      }
    }
      
  } else {

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
      if (x0>x1) {
	O2SCL_CONV2_RET("Discontinuities too close to resolve ",
			" in tov_love::calc_y().",
			o2scl::exc_etol,err_nonconv);
      }
      int ois_ret=oisp->solve_final_value(x0,x1,h,1,y,yout,od);
      if (ois_ret!=0 && addl_testing) {
	cout << "H0: " << x0 << " " << x1 << " " << h << " "
	     << y[0] << " " << yout[0] << endl;
	cout << "H1: ";
	vector_out(std::cout,disc,true);
	cout << "H2: " << tab->get_nlines() << " "
	     << tab->get("r",0) << " "
	     << tab->get("r",tab->get_nlines()-1) << endl;
	cout << "H3: " << tab->get_constant("rad") << " "
	     << tab->get_constant("mass") << endl;
	show_ode=1;
	oisp->verbose=1;
	oisp->solve_final_value(x0,x1,h,1,y,yout,od);
	O2SCL_CONV2_RET("ODE function solve_final_value() failed ",
			" in tov_love::calc_y().",
			o2scl::exc_efailed,err_nonconv);
      }

      // Add the correction at the discontinuity
      if (j!=disc.size()) {
	y[0]=yout[0]+(tab->interp("r",disc[j]+delta,"ed")-
		      tab->interp("r",disc[j]-delta,"ed"))/
	  (tab->interp("r",disc[j],"gm")/
	   (4.0*o2scl_const::pi*disc[j]*disc[j]*disc[j])+
	   tab->interp("r",disc[j],"pr"));
      } else {
	y[0]=yout[0];
      }

      // Proceed to the next interval
    }

    // Set the final value of y at r=R
    yR=y[0];
  }
  
  beta=schwarz_km/2.0*gm/R;

  k2=eval_k2(beta,yR);

  // lambda in km^5
  lambda_km5=2.0/3.0*k2*pow(R,5.0);

  // First compute in Msun*km^2*s^2
  lambda_cgs=2.0/3.0*k2*pow(R,5.0)/
    pow(o2scl_mks::speed_of_light/1.0e3,2.0)/schwarz_km;

  // Convert to g*cm^2*s^2
  lambda_cgs*=o2scl_cgs::solar_mass*1.0e10;

  return 0;
}

int tov_love::calc_H(double &yR, double &beta, double &k2, 
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

  return 0;
}
