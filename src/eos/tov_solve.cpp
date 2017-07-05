/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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

#include <cstdlib>

#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/tov_solve.h>
#include <o2scl/root_cern.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

tov_solve::tov_solve() : out_table(new table_units<>) {

  // EOS
  eos_set=false;
  te=0;

  // Verbosity parameter
  verbose=1;

  // ODE solver
  as_ptr=&def_stepper;
  step_min=1.0e-4;
  step_max=0.05;
  step_start=4.0e-3;

  // Minimum log(pressure)
  min_log_pres=-1.0e2;

  // Allocation sizes
  buffer_size=1e5;
  max_table_size=400;

  // Units
  efactor=1.0;
  pfactor=1.0;
  nfactor=1.0;
  eunits="Msun/km^3";
  punits="Msun/km^3";
  nunits="1/fm^3";

  // Pressure loop for mvsr()
  prbegin=7.0e-7;
  prend=8.0e-3;
  princ=1.1;

  // Guess for pressure for fixed()
  fixed_pr_guess=5.2e-5;

  // Maximum mass initial guess
  max_begin=7.0e-5;
  max_end=5.0e-3;
  max_inc=1.3;

  // Solver
  mroot_ptr=&def_solver;
  def_solver.tol_abs=1.0e-12;
  def_solver.tol_rel=1.0e-7;
  def_solver.verbose=0;
  def_solver.ntrial=100;

  // Minimizer for max()
  min_ptr=&def_min;

  // Constants
  schwarz_km=o2scl_cgs::schwarzchild_radius/1.0e5;

  out_table->add_constant("schwarz",schwarz_km);
  out_table->add_constant("Msun",o2scl_mks::solar_mass);
  out_table->add_constant("pi",pi);
  out_table->add_constant("mproton",o2scl_mks::mass_proton);
  
  // Stellar properties
  mass=0.0;
  rad=0.0;
  bmass=0.0;
  gpot=0.0;
  domega_rat=0.0;

  // Other options
  gen_rel=true;
  ang_vel=false;
  calc_gpot=false;
  err_nonconv=true;
  
  // Initial value for target mass
  tmass=0.0;

  max_integ_steps=100000;
  pmax_default=1.0e20;
  pcent_max=pmax_default;
  reformat_results=true;
  
  // Mass of a proton in kg
  baryon_mass=o2scl_mks::mass_proton;
  
  // Initialize ODE function pointer
  ofm=std::bind
    (std::mem_fn<int(double,size_t,const ubvector &,ubvector &)>
     (&tov_solve::derivs),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3,std::placeholders::_4);

}

tov_solve::~tov_solve() {
}

void tov_solve::set_units(double s_efactor, double s_pfactor, 
			  double s_nfactor) {
  efactor=s_efactor;
  pfactor=s_pfactor;
  nfactor=s_nfactor;
  return;
}

void tov_solve::set_units(std::string leunits, 
			      std::string lpunits, std::string lnunits) {
  eunits=leunits;
  punits=lpunits;
  nunits=lnunits;

  if (eunits.length()==0 || eunits=="Msun/km^3" ||
      eunits=="solarmass/km^3") {
    efactor=1.0;
  } else {
    efactor=o2scl_settings.get_convert_units().convert
      (eunits,"Msun/km^3",1.0);
  }

  if (punits.length()==0 || punits=="Msun/km^3" ||
      punits=="solarmass/km^3") {
    pfactor=1.0;
  } else {
    pfactor=o2scl_settings.get_convert_units().convert
      (punits,"Msun/km^3",1.0);
  }

  if (nunits=="1/cm^3") {
    nfactor=1.0e-39;
  } else if (nunits=="1/m^3") {
    nfactor=1.0e-42;
  } else if (nunits=="1/fm^3") {
    nfactor=1.0;
  } else if (nunits=="") {
    nfactor=1.0;
  } else {
    nfactor=o2scl_settings.get_convert_units().convert
      (nunits,"1/fm^3",1.0);
  }

  out_table->add_constant("efactor",efactor);
  out_table->add_constant("pfactor",pfactor);
  out_table->add_constant("nfactor",nfactor);

  return;
}
 
int tov_solve::derivs(double r, size_t nv, const ubvector &y, 
		      ubvector &dydx) {
  
  // We don't throw the error handler for non-finite values here 
  // in case the ODE solver can recover from a bad step

  if (y[1]<=min_log_pres || y[1]>0.0) {
    return exc_efailed;
  }

  // Index counter
  size_t ix;

  // Handle zero radius case
  if (r==0.0) {
    dydx[0]=0.0;
    dydx[1]=0.0;
    ix=2;
    if (calc_gpot) {
      dydx[ix]=0.0;
      ix++;
      if (ang_vel) {
	dydx[ix]=0.0;
	ix++;
	dydx[ix]=0.0;
	ix++;
      }
    }
    if (te->has_baryons()) {
      dydx[ix]=0.0;
      ix++;
    }
    return success;
  }

  double pres=exp(y[1]);
  double gm=y[0];

  if (!std::isfinite(pres)) {
    return exc_efailed;
  }

  // The function get_eden() is now already in the proper units,
  // so there's no need for unit conversion here
  double ed, nb;
  te->ed_nb_from_pr(pres,ed,nb);
  
  if (!std::isfinite(ed)) {
    return exc_efailed;
  }

  double term1, term2, term3;
  dydx[0]=4.0*pi*r*r*ed;
  if (gen_rel) {
    term1=ed+pres;
    if (gm<0.0) term2=4.0*pi*pow(r,3.0)*pres;
    else term2=gm+4.0*pi*pow(r,3.0)*pres;
    term3=r-schwarz_km*gm;
  } else {
    term1=ed;
    if (gm<0.0) term2=0.0;
    else term2=gm;
    term3=r;
  }
  dydx[1]=-schwarz_km/2.0/r*term1*term2/term3/pres;

  ix=2;
  double gp=0.0;
  if (calc_gpot) {
    gp=y[ix];
    double dgpdr=-1.0/term1*dydx[1]*pres;
    dydx[ix]=dgpdr;
    ix++;
    if (ang_vel) {
      double j=sqrt(1.0-schwarz_km*gm/r)*exp(-gp);
      double djdr=-sqrt(1.0-schwarz_km*gm/r)*exp(-gp)*dgpdr+
	exp(-gp)/2.0/sqrt(1.0-schwarz_km*gm/r)*schwarz_km/r*(gm/r-dydx[0]);
      dydx[ix]=-4.0*r*r*r*djdr*y[ix+1];
      dydx[ix+1]=y[ix]/pow(r,4.0)/j;
      ix+=2;
    }
  }
  if (te->has_baryons()) {
    if (gm<0.0) {
      dydx[ix]=4.0*pi*r*r*nb*baryon_mass/o2scl_mks::solar_mass*
	1.0e54;
    } else {
      dydx[ix]=4.0*pi*r*r*nb*baryon_mass/o2scl_mks::solar_mass*
	1.0e54/sqrt(1.0-schwarz_km*gm/r);
    }
    ix++;
  }

  return success;
}

void tov_solve::make_unique_name(string &col, std::vector<string> &cnames) {
  bool done;
  do {
    // Exhaustively go through 'cnames' looking for 'col'
    done=true;
    for(size_t i=0;i<cnames.size();i++) {
      if (col==cnames[i]) {
	done=false;
	i=cnames.size();
      }
    }
    // If found, append an underscore
    if (done==false) {
      col+='_';
    }
  } while (done==false);
  return;
}

void tov_solve::column_setup(bool mvsr_mode) {

  // ---------------------------------------------------------------
  // Add internal columns

  vector<string> inames, iunits;
  string line2;
  inames.push_back("gm");
  iunits.push_back("Msun");
  inames.push_back("r");
  iunits.push_back("km");
  if (calc_gpot) {
    inames.push_back("gp");
    iunits.push_back("");
    if (ang_vel) {
      inames.push_back("rjw");
      iunits.push_back("km^4");
      inames.push_back("omega_rat");
      iunits.push_back("");
    }
  }
  if (te->has_baryons()) {
    inames.push_back("bm");
    iunits.push_back("Msun");
  }
  inames.push_back("pr");
  iunits.push_back(punits);
  inames.push_back("ed");
  iunits.push_back(eunits);
  if (te->has_baryons()) {
    inames.push_back("nb");
    iunits.push_back(nunits);
  }
  inames.push_back("sg");
  iunits.push_back("1/km");
  inames.push_back("rs");
  iunits.push_back("");
  inames.push_back("dmdr");
  iunits.push_back("Msun/km");
  inames.push_back("dlogpdr");
  iunits.push_back("1/km");
  if (calc_gpot) {
    inames.push_back("dgpdr");
    iunits.push_back("1/km");
  }
  if (te->has_baryons()) {
    inames.push_back("dbmdr");
    iunits.push_back("Msun/km");
  }
  if (mvsr_mode && pr_list.size()>0) {
    for(size_t i=0;i<pr_list.size();i++) {
      inames.push_back(((string)"r")+szttos(i));
      iunits.push_back("km");
      inames.push_back(((string)"gm")+szttos(i));
      iunits.push_back("Msun");
      if (te->has_baryons()) {
	inames.push_back(((string)"bm")+szttos(i));
	iunits.push_back("Msun");
      }
    }
  }
  for(size_t il=0;il<inames.size();il++) {
    line2+=inames[il];
    line2+=' ';
  }
  
  out_table->line_of_names(line2);
  for(size_t il=0;il<inames.size();il++) {
    out_table->set_unit(inames[il],iunits[il]);
  }
  
  return;
}

void tov_solve::make_table() {

  // ---------------------------------------------------------------
  // Clear previously stored data and set up table

  out_table->clear();
  column_setup();

  // ---------------------------------------------------------------
  // Add constants for mass and radius
  
  out_table->add_constant("mass",mass);
  out_table->add_constant("rad",rad);

  // ---------------------------------------------------------------
  // Create output file of profile for star

  if (ix_last+1<=max_table_size) {
    out_table->set_nlines(ix_last+1);
  } else {
    out_table->set_nlines(max_table_size);
  }

  // We scroll through all the buffer entries, and at the end
  // of the loop, add an additional two iterations to double
  // check that the first and last row are correct
  for(size_t j=0;j<ix_last+3;j++) {

    // Set buffer index and table index
    size_t bix=j;
    size_t tix=((size_t)(((double)out_table->get_nlines())*((double)j)/
			 ((double)(ix_last+1))));
    if (j==ix_last+1) {
      bix=0;
      tix=0;
    } else if (j==ix_last+2) {
      bix=ix_last;
      tix=max_table_size-1;
    }

    // Double check that we're not out of range
    if (bix>=ix_last+1) bix=ix_last;
    if (tix>=out_table->get_nlines()) tix=out_table->get_nlines()-1;
    
    ubvector line;
    size_t iv=0;

    // Output enclosed mass and radius
    out_table->set("gm",tix,rky[bix][0]);
    out_table->set("r",tix,rkx[bix]);

    iv=2;
    // Gravitational potential
    if (calc_gpot) {
      out_table->set("gp",tix,rky[bix][iv]);
      iv++;
      if (ang_vel) {
	out_table->set("rjw",tix,rky[bix][iv]);
	iv++;
	out_table->set("omega_rat",tix,rky[bix][iv]);
	iv++;
      }
    }

    // Enclosed baryon mass
    if (te->has_baryons()) {
      out_table->set("bm",tix,rky[bix][iv]);
      iv++;
    }

    // Energy density, pressure, and baryon density
    if (rky[bix][1]>min_log_pres) {
      double ed, nb;
      if (!std::isfinite(exp(rky[bix][1]))) {
	O2SCL_ERR2("Pressure not finite in ",
		   "tov_solve::make_table().",exc_efailed);
      }
      te->ed_nb_from_pr(exp(rky[bix][1]),ed,nb);
      // Convert pressure, energy density, and baryon density to user 
      // units by dividing by their factors
      out_table->set("pr",tix,exp(rky[bix][1])/pfactor);
      out_table->set("ed",tix,ed/efactor);
      if (te->has_baryons()) {
	out_table->set("nb",tix,nb/nfactor);
      }
    } else {
      out_table->set("pr",tix,0.0);
      out_table->set("ed",tix,0.0);
      if (te->has_baryons()) {
	out_table->set("nb",tix,0.0);
      }
    }

    // surface gravity and redshift
    if (rkx[bix]!=0.0) {
      out_table->set("sg",tix,(schwarz_km/2.0*rky[bix][0]/rkx[bix]/rkx[bix]/
			       sqrt(1.0-schwarz_km*rky[bix][0]/rkx[bix])));
      out_table->set("rs",tix,1.0/sqrt(1.0-rky[bix][0]*schwarz_km/
				       rkx[bix])-1.0);
    } else {
      out_table->set("sg",tix,0.0);
      out_table->set("rs",tix,0.0);
    }

    // derivatives
    out_table->set("dmdr",tix,rkdydx[bix][0]);
    out_table->set("dlogpdr",tix,rkdydx[bix][1]);
    iv=2;
    if (calc_gpot) {
      out_table->set("dgpdr",tix,rkdydx[bix][iv]);
      iv++;
    }
    if (te->has_baryons()) {
      out_table->set("dbmdr",tix,rkdydx[bix][iv]);
      iv++;
    }
    
    // Check for non-finite values
    for(size_t ik=0;ik<out_table->get_ncolumns();ik++) {
      if (!std::isfinite(out_table->get(ik,tix))) {
	O2SCL_ERR((((string)"Non-finite value for column '")+
		   out_table->get_column_name(ik)+
		   "' in tov_solve::make_table().").c_str(),exc_efailed);
      }
    }
      
  }

  return;
}

double tov_solve::max_fun(double maxx) {
  ubvector x(1), y(1);
  double retval;
  x[0]=maxx;

  int ret=integ_star(1,x,y);

  if (maxx<0.0) retval=1000.0*fabs(maxx);
  else if (ret!=0) retval=1000.0;
  else retval=-mass;

  return retval;
}

int tov_solve::integ_star(size_t ndvar, const ubvector &ndx, 
			  ubvector &ndy) {
  
  // ---------------------------------------------------------------
  // Initial failure conditions
  
  if (ndx[0]>pcent_max) {
    return cent_press_large;
  }

  if (ndx[0]<0.0) {
    // We don't call the error handler here, because the default
    // solver might be able to recover from a negative pressure. See
    // the documentation in tov_solve.h for a discussion of this.
    return cent_press_neg;
  }

  // ---------------------------------------------------------------
  // Count number of diff eqs. to solve

  size_t nvar=2;
  if (calc_gpot) {
    nvar++;
    if (ang_vel) nvar+=2;
  }
  if (te->has_baryons()) nvar++;

  // ---------------------------------------------------------------
  // Resize and allocate memory if necessary

  if (rkx.size()!=buffer_size || rky.size()==0 || rky[0].size()!=nvar) {
    rkx.resize(buffer_size);
    rky.resize(buffer_size);
    rkdydx.resize(buffer_size);
    for(size_t i=0;i<buffer_size;i++) {
      rky[i].resize(nvar);
      rkdydx[i].resize(nvar);
    }
  }
  ubvector yerr(nvar);
  
  // ---------------------------------------------------------------
  // Set up inital point

  rkx[0]=0.0;
  rky[0][0]=0.0;
  rky[0][1]=log(ndx[0]);
  size_t iv=2;
  if (calc_gpot) {
    rky[0][iv]=0.0;
    iv++;
    if (ang_vel) {
      rky[0][iv]=0.0;
      iv++;
      rky[0][iv]=0.5;
      iv++;
    }
  }
  if (te->has_baryons()) {
    rky[0][iv]=0.0;
    iv++;
  }

  // ---------------------------------------------------------------
    
  double outrad=0.0;
  if (verbose>=2) {
    cout << "Central pressure: " << ndx[0] << endl;
    cout << "   it Radius       Enc. mass    X=Log(P)          P"
	 << endl;
  }

  size_t ix=0, ix_next=1;
  int test;

  // ---------------------------------------------------------------
  // Main loop

  bool done=false;
  size_t it;
  for(it=0;it<max_integ_steps && done==false;it++) {
    
    // ---------------------------------------------------------------
    // Fix step size if too large or too small
    
    double h=step_start;
    if (h>step_max) h=step_max;
    if (h<step_min) h=step_min;

    // ---------------------------------------------------------------
    // Take an adaptive step 
    
    test=as_ptr->astep_full(rkx[ix],rkx[ix]+step_max,rkx[ix_next],h,nvar,
			    rky[ix],rkdydx[ix],rky[ix_next],yerr,
			    rkdydx[ix_next],ofm);
    if (test!=0) {
      done=true;
    }

    if (done==false) {

      ix++;
      ix_next++;

      // ---------------------------------------------------------------
      // Verbose output
    
      if ((verbose>=3 && it%20==0) || (verbose>1 && rkx[ix]>outrad)) {
	cout.width(5);
	cout << it << " " << rkx[ix] << " " << rky[ix][0] << " " 
	     << rky[ix][1] << " " << exp(rky[ix][1]) << endl;
	outrad+=1.0;
      }

      // ---------------------------------------------------------------
      // We've run out of room, so thin the data accordingly. Always
      // leave one empty row at the end for the final calculations
      // of mass and radius done outside of the loop below. 
    
      if (ix_next>=buffer_size-1) {

	//cout << "Rearrangement. ix=" << ix << " ix_next=" << ix_next 
	//<< " buffer_size=" << buffer_size << endl;

	// Leave initial point fixed, and thin everything except the
	// last point
	for(size_t j=2;j<=ix-1;j++) {
	  //cout << "Copying " << j << " to " << j/2 << endl;
	  rkx[j/2]=rkx[j];
	  for(size_t k=0;k<nvar;k++) {
	    rky[j/2][k]=rky[j][k];
	    rkdydx[j/2][k]=rkdydx[j][k];
	  }
	}

	// Copy the last point to its new location, ix_new
	size_t ix_new=(ix-1)/2;
	//cout << "Copying (last) " << ix << " to " << ix_new << endl;
	rkx[ix_new]=rkx[ix];
	for(size_t k=0;k<nvar;k++) {
	  rky[ix_new][k]=rky[ix][k];
	  rkdydx[ix_new][k]=rkdydx[ix][k];
	}

	// Set ix and ix_next
	ix=ix_new;
	ix_next=ix+1;
	//cout << "New values ix=" << ix << " ix_next=" << ix_next << endl;
      }

      // End of loop 'if (done==false)'
    }

    // End of main loop
  }

  if (it>=max_integ_steps) {
    O2SCL_CONV((((string)"Integration of star with central pressure ")+
		o2scl::dtos(ndx[0])+" exceeded "+
		o2scl::szttos(max_integ_steps)+
		" steps in integ_star().").c_str(),
	       exc_efailed,err_nonconv);
    return over_max_steps;
  }
  
  // --------------------------------------------------------------
  // Final calculations of mass and radius. This simple linear
  // interpolation to zero pressure could be improved, possibly by
  // using an analytical solution

  rad=-1.0/rkdydx[ix][1]+rkx[ix];
  mass=rky[ix][0]-rkdydx[ix][0]*(rkx[ix]-rad);

  // If extrapolation to zero pressure is more than 10 percent,
  // then register as a failure
  if (rad>1.1*rkx[ix]) {
    O2SCL_CONV("Last radial step too large in integ_star().",
	       exc_efailed,err_nonconv);
    return last_step_large;
  }

  // Verbose output
  if (verbose>=3) {
    cout << "Final surface interpolation: " << endl;
    cout.setf(ios::showpos);
    cout << "r, m(r), log(P): " 
	 << rkx[ix] << " " << rky[ix][0] << " " << rky[ix][1] << endl;
    cout << "r, m(r), log(P): " 
	 << rad << " " << mass << " " << 0.0 << endl;
    cout << "r, dm/dr, dlog(P)/dr: " 
	 << rkx[ix] << " " << rkdydx[ix][0] << " " << rkdydx[ix][1] << endl;
    cout.unsetf(ios::showpos);
  }

  // Extrapolate final gravitational potential and baryon mass
  iv=2;
  double lastgpot=0.0;
  if (calc_gpot) {
    lastgpot=rky[ix][iv]-rkdydx[ix][iv]*(rkx[ix]-rad);
    iv++;
    if (ang_vel) iv+=2;
  }
  if (te->has_baryons()) {
    bmass=rky[ix][iv]-rkdydx[ix][iv]*(rkx[ix]-rad);
  }
  
  // --------------------------------------------------------------
  // Output if verbose>0

  if (verbose>0) {
    cout.precision(4);
    cout << "Central P: " << ndx[0] << " (Msun/km^3), M: " 
	 <<  mass << " (Msun), R: " <<  rad << " (km)" << endl;
    cout.precision(6);
    if (verbose>=3) {
      cout << "Press a key and 'Enter' to continue." << endl;
      char ch;
      cin >> ch;
    }
  }

  // --------------------------------------------------------------
  // Final value

  ndy[0]=mass-tmass;
  if (!integ_star_final) {
    return 0;
  }
  
  // --------------------------------------------------------------
  // Store the last point for radius, mass, and pressure

  ix_last=ix_next;
  rkx[ix_last]=rad;
  rky[ix_last][0]=mass;
  rky[ix_last][1]=min_log_pres;

  // --------------------------------------------------------------

  iv=2;
  if (calc_gpot) {

    // Add a constant to the grav. pot to satisfy the boundary
    // conditions
    
    rky[ix_last][iv]=lastgpot;
    double phi_shift=0.5*log(1-schwarz_km*mass/rad)-lastgpot;
    for(size_t k=0;k<=ix_last;k++) {
      rky[k][iv]+=phi_shift;
    }
    gpot=rky[ix_last][iv];
    iv++;
    
    // --------------------------------------------------------------
    // Adjustments for rotation
    
    if (ang_vel) {

      last_rjw=rky[ix][iv]*exp(-phi_shift);
      last_f=rky[ix][iv+1];
      double f_corr=last_f+last_rjw/pow(rad,3.0)/3.0;

      // Correction for rjw
      for(size_t k=0;k<ix_last;k++) {
	rky[k][iv]*=exp(-phi_shift)/f_corr;
      }
      rky[ix_last][iv]=last_rjw/f_corr;
      iv++;

      // Correction for omega_rat
      rky[ix_last][iv]=last_f;
      for(size_t k=0;k<=ix_last;k++) {
	rky[k][iv]/=f_corr;
      }
      iv++;

      last_rjw/=f_corr;
      last_f/=f_corr;
      
      // Compute d ( (omega bar)/Omega ) / dr
      domega_rat=last_rjw/pow(rad,4.0);
    }
  }
  
  // --------------------------------------------------------------
  // Store the last point for baryonic mass

  if (te->has_baryons()) {
    rky[ix_last][iv]=bmass;
    iv++;
  }

  // --------------------------------------------------------------
  // Last row of derivatives
  
  for(size_t k=0;k<nvar;k++) {
    rkdydx[ix_last][k]=rkdydx[ix][k];
  }

  return 0;
}

int tov_solve::mvsr() {

  int info=0;
  pcent_max=pmax_default;

  if (eos_set==false) {
    O2SCL_ERR
      ("EOS not specified tov_solve::mvsr().",exc_efailed);
  }

  if (verbose>0) cout << "Mass versus radius mode." << endl;

  // ---------------------------------------------------------------
  // Clear previously stored data and setup table
  
  out_table->clear();
  column_setup(true);

  // ---------------------------------------------------------------
  // Main loop

  ubvector x(1), y(1);
  for (x[0]=prbegin;((prend>prbegin && x[0]<=prend) ||
		     (prend<prbegin && x[0]>=prend));) {
    
    integ_star_final=true;
    int ret=integ_star(1,x,y);
    if (ret!=0 && info==0) {
      O2SCL_CONV((((string)"Integration of star with central pressure ")
		  +dtos(x[0])+" failed in mvsr().").c_str(),exc_efailed,
		 err_nonconv);
      info+=mvsr_integ_star_failed+ret;
    }
      
    // --------------------------------------------------------------
    // Fill line of data for table

    std::vector<double> line;

    // output mass and radius
    line.push_back(mass);
    line.push_back(rad);

    // Gravitational potential and angular velocity columns
    if (calc_gpot) {
      line.push_back(gpot);
      if (ang_vel) {
	line.push_back(last_rjw);
	line.push_back(last_f);
      }
    }

    // output baryon mass
    if (te->has_baryons()) line.push_back(bmass);
    
    // output central pressure, energy density, and baryon density

    double ed, nb;
    if (!std::isfinite(x[0])) {
      O2SCL_ERR2("Central pressure not finite in ",
		 "tov_solve::mvsr().",exc_efailed);
    }
    te->ed_nb_from_pr(x[0],ed,nb);
    
    // Convert pressure, energy density, and baryon density to user 
    // units by dividing by their factors
    line.push_back(x[0]/pfactor);
    line.push_back(ed/efactor);
    if (te->has_baryons()) {
      line.push_back(nb/nfactor);
    }

    // output surface gravity and redshift

    if (rad!=0.0) {
      line.push_back(schwarz_km/2.0*mass/rad/rad/
		     sqrt(1.0-schwarz_km*mass/rad));
      line.push_back(1.0/sqrt(1.0-mass*schwarz_km/rad)-1.0);
    } else {
      line.push_back(0.0);
      line.push_back(0.0);
    }
    
    // output derivatives

    line.push_back(0.0);
    line.push_back(0.0);
    if (calc_gpot) line.push_back(0.0);
    if (te->has_baryons()) line.push_back(0.0);

    // Radius interpolation
    if (pr_list.size()>0) {
      iop.set_type(itp_linear);
      ubvector lpr_col(rky.size()), gm_col(rky.size()), bm_col(rky.size());
      for(size_t ii=0;ii<rky.size();ii++) {
	lpr_col[ii]=rky[ii][1];
	gm_col[ii]=rky[ii][0];
	if (te->has_baryons()) {
	  size_t index=2;
	  if (calc_gpot) {
	    index++;
	    if (ang_vel) index+=2;
	  }
	  bm_col[ii]=rky[ii][index];
	}
      }
      for(size_t ii=0;ii<pr_list.size();ii++) {
	double thisr=iop.eval(log(pr_list[ii]*pfactor),
				ix_last-1,lpr_col,rkx);
	double thisgm=iop.eval(log(pr_list[ii]*pfactor),
				 ix_last-1,lpr_col,gm_col);
	if (!std::isfinite(thisr)) {
	  string str=((string)"Obtained non-finite value when ")+
	    "interpolating radius for pressure "+dtos(pr_list[ii])+
	    " in tov_solve::mvsr().";
	  O2SCL_ERR(str.c_str(),exc_efailed);
	}
	line.push_back(thisr);
	if (!std::isfinite(thisgm)) {
	  string str=((string)"Obtained non-finite value when ")+
	    "interpolating gravitational mass for pressure "+dtos(pr_list[ii])+
	    " in tov_solve::mvsr().";
	  O2SCL_ERR(str.c_str(),exc_efailed);
	}
	line.push_back(thisgm);
	if (te->has_baryons()) {
	  double thisbm=iop.eval(log(pr_list[ii]*pfactor),
				 ix_last-1,lpr_col,bm_col);
	  if (!std::isfinite(thisbm)) {
	    string str=((string)"Obtained non-finite value when ")+
	      "interpolating baryon mass for pressure "+dtos(pr_list[ii])+
	      " in tov_solve::mvsr().";
	    O2SCL_ERR(str.c_str(),exc_efailed);
	  }
	  line.push_back(thisbm);
	}
      }
    }

    // --------------------------------------------------------------
    // Copy line of data to table
    
    out_table->line_of_data(line.size(),&(line[0]));
    if (line.size()!=out_table->get_ncolumns()) {
      O2SCL_ERR("Table size problem in tov_solve::mvsr().",
		exc_esanity);
    }
    
    // --------------------------------------------------------------
    // Get next central pressure

    x[0]*=princ;

  }

  // Find the row that refers to the maximum mass star
  size_t ix=out_table->lookup("gm",out_table->max("gm"));
  pcent_max=out_table->get("pr",ix);

  return info;
}

int tov_solve::max() {
  
  int info=0;

  pcent_max=pmax_default;

  if (verbose>0) cout << "Maximum mass mode." << endl;

  // --------------------------------------------------------------
  // Handle basic exceptions

  if (!gen_rel) {
    O2SCL_ERR2("Maximum mass not supported for non-relativistic ",
	       "stars in tov_solve::max().",exc_efailed);
  }
  if (ang_vel && !calc_gpot) {
    O2SCL_ERR2("Requested rotation without potential in ",
	       "tov_solve::max().",exc_efailed);
  }
  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified and/or memory not allocated in ",
	       "tov_solve::max().",exc_efailed);
  }

  // --------------------------------------------------------------

  integ_star_final=false;

  // --------------------------------------------------------------
  // Create a good initial guess for the minimizer using a
  // naive search over several central pressures

  double initmin=max_begin;
  double initval=max_fun(initmin);
  
  ubvector x(1), y(1);
  for(x[0]=max_begin*max_inc;x[0]<=max_end;x[0]*=max_inc) {
    y[0]=max_fun(x[0]);
    if (y[0]<initval) {
      initval=y[0];
      initmin=x[0];
    }
  }
  x[0]=initmin;

  // --------------------------------------------------------------
  // Full minimization

  funct mm=std::bind(std::mem_fn<double(double)>
		       (&tov_solve::max_fun),
		       this,std::placeholders::_1);
  if (min_ptr->min(x[0],y[0],mm)!=0) {
    info+=max_minimizer_failed;
  }

  // --------------------------------------------------------------
  // Final call to integ_star()

  integ_star_final=true;
  int ret=integ_star(1,x,y);
  if (ret!=0) {
    O2SCL_CONV("Last call to integ_star() failed in max().",exc_efailed,
	       err_nonconv);
    info+=max_integ_star_failed+ret;
  }
  
  // --------------------------------------------------------------
  // Final call to integ_star()

  if (verbose>0) {
    cout << "Maximum gravitational mass is: " << mass << endl;
    if (te->has_baryons()) {
      cout << "Corresponding baryon mass is: " << bmass << endl;
    }
    cout << "Corresponding radius is: " << rad << endl;
  }

  // --------------------------------------------------------------
  // Set central pressure of maximum mass star

  pcent_max=x[0];

  // --------------------------------------------------------------
  // Output stellar profile

  if (reformat_results) {
    make_table();
  }

  return info;
}

int tov_solve::fixed_pr(double pcent) {

  int info=0;

  // --------------------------------------------------------------
  // Handle basic exceptions

  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified in ",
	       "tov_solve::fixed_pr().",exc_efailed);
  }
  if (ang_vel && !calc_gpot) {
    O2SCL_ERR2("Requested rotation without potential in ",
	       "tov_solve::fixed_pr().",exc_efailed);
  }

  // --------------------------------------------------------------
  
  if (verbose>0) {
    cout << "Central pressure " << pcent << " (" << punits << ") "
	 << pcent*pfactor << " (Msun/km^3) mode." << endl;
  }
  pcent*=pfactor;
  
  // --------------------------------------------------------------

  // Stellar integration

  ubvector x(1), y(1);
  x[0]=pcent;
  integ_star_final=true;
  int ret=integ_star(1,x,y);
  if (ret!=0) {
    info+=fixed_integ_star_failed+ret;
    O2SCL_CONV("Failed to integrate star in tov_solve::fixed_pr().",
	       exc_efailed,err_nonconv);
  }

  if (verbose>0) {
    cout << "Gravitational mass is: " << mass << endl;
    cout << "Radius is: " << rad << endl;
    if (te->has_baryons()) {
      cout << "Baryon mass is: " << bmass << endl;
    }
  }

  if (reformat_results) {
    // Output stellar profile
    make_table();
  }
  
  return info;
}

int tov_solve::fixed(double target_mass, double pmax) {

  int info=0;

  // --------------------------------------------------------------
  // Handle basic exceptions

  if (eos_set==false) {
    O2SCL_ERR2("EOS not specified in ",
	       "tov_solve::fixed().",exc_efailed);
  }
  if (ang_vel && !calc_gpot) {
    O2SCL_ERR2("Requested rotation without potential in ",
	       "tov_solve::fixed().",exc_efailed);
  }

  // --------------------------------------------------------------
    
  tmass=target_mass;
  if (verbose>0) cout << tmass << " solar mass mode." << endl;

  // --------------------------------------------------------------

  // Compute maximum mass star if necessary, and set pmax_cent
  if (pmax>=pmax_default || tmass<0.0) {
    max();
  } else if (pmax<pmax_default) {
    pcent_max=pmax;
  }

  // Negative value given, so adjust target mass accordingly.
  if (tmass<0.0) tmass+=mass;
  
  // Initial guess

  ubvector x(1), y(1);
  x[0]=fixed_pr_guess;
  
  mm_funct fmf=std::bind
    (std::mem_fn<int(size_t,const ubvector &,ubvector &)>
     (&tov_solve::integ_star),
     this,std::placeholders::_1,std::placeholders::_2,
     std::placeholders::_3);

  // Set to false to avoid uninit'ed variable warnings
  bool save1=false, save2=false;
  if (err_nonconv==false) {
    save1=mroot_ptr->err_nonconv;
    save2=def_solver.def_jac.err_nonconv;
    mroot_ptr->err_nonconv=false;
    def_solver.def_jac.err_nonconv=false;
  }
  integ_star_final=false;

  // Before trying to solve, evaluate the initial guess
  if (integ_star(1,x,y)!=o2scl::success) {
    // If it failed, we must return (otherwise the
    // solver might throw an exception)
    info+=fixed_solver_failed;
    O2SCL_CONV("Initial guess failed in tov_solve::fixed().",
	       exc_efailed,err_nonconv);
    return info;
  }

  int ret=mroot_ptr->msolve(1,x,fmf);
  if (ret!=0) {
    info+=fixed_solver_failed;
  }

  if (err_nonconv==false) {
    mroot_ptr->err_nonconv=save1;
    def_solver.def_jac.err_nonconv=save2;
  }
  
  // Calculate gravitational potential and enclosed baryon mass
  integ_star_final=true;
  ret=integ_star(1,x,y);
  if (ret!=0) {
    info+=fixed_integ_star_failed+ret;
    O2SCL_CONV("Failed to integrate star in tov_solve::fixed().",
	       exc_efailed,err_nonconv);
  }

  if (verbose>0) {
    cout << "Gravitational mass is: " << mass << endl;
    cout << "Radius is: " << rad << endl;
    if (te->has_baryons()) {
      cout << "Baryon mass is: " << bmass << endl;
    }
  }

  if (fabs(y[0])>1.0e-3) {
    if (verbose>0) {
      cout << "Solution failed for target mass=" << tmass 
	   << " and gave mass=" << mass << "." << endl;
    }
    info+=fixed_wrong_mass;
    O2SCL_CONV("Failed to return correct mass in tov_solve::fixed().",
	       exc_efailed,err_nonconv);
  }

  if (reformat_results) {
    // Output stellar profile
    make_table();
  }
  
  return info;
}

