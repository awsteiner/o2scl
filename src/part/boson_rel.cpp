/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

#include <o2scl/boson_rel.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_const;

//--------------------------------------------
// boson_rel class

boson_rel::boson_rel() {
  density_mroot=&def_density_mroot;
  nit=&def_nit;
  dit=&def_dit;
  verify_ti=false;
  use_expansions=true;
  deg_limit=-0.5;
}

boson_rel::~boson_rel() {
}

void boson_rel::set_inte(inte<> &l_nit, inte<> &l_dit) {
  nit=&l_nit;
  dit=&l_dit;
  return;
}

void boson_rel::calc_mu(boson &b, double temper) {

  if (temper<=0.0) {
    O2SCL_ERR2("Temperature less than or equal to zero in ",
	       "boson_rel::calc_mu().",exc_einval);
  }
  if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }

  double psi;
  if (b.inc_rest_mass) {
    psi=(b.nu-b.ms)/temper;
  } else {
    psi=(b.nu+(b.m-b.ms))/temper;
  }

  if (psi>0.0) {
    O2SCL_ERR2("Chemical potential must be smaller than mass in ",
               "boson_rel::calc_mu().",o2scl::exc_einval);
  }

  // Try the non-degenerate expansion if psi is small enough
  if (use_expansions) {
    bool acc=this->calc_mu_ndeg(b,temper,1.0e-14);
    if (verbose>1) {
      std::cout << "calc_mu(): non-deg expan " << acc
                << std::endl;
    }
    if (acc) {
      /*
      unc.n=f.n*tol_expan;
      unc.ed=f.ed*tol_expan;
      unc.pr=f.pr*tol_expan;
      unc.en=f.en*tol_expan;
      */
      //last_method=4;
      return;
    }
  }
  
  bool deg=true;
  if (psi<deg_limit) deg=false;

  if (verbose>1) {
    cout << "boson_rel::calc_mu() psi: " << psi << endl;
  }
  
  if (deg) {
    
    // Compute the upper limit for degenerate integrals

    double arg, upper_limit_fac=30.0;
    if (b.inc_rest_mass) {
      arg=pow(upper_limit_fac*temper+b.nu,2.0)-b.ms*b.ms;
    } else {
      arg=pow(upper_limit_fac*temper+b.nu+b.m,2.0)-b.ms*b.ms;
    }
    double ul;
    if (arg>0.0) {
      ul=sqrt(arg);
    } else {
      O2SCL_ERR2("Zero density in degenerate limit in boson_rel::",
		 "calc_mu().",exc_efailed);
      return;
    }

    cout << "boson_rel::calc_mu() arg, ulf, ul: " << arg << " "
         << upper_limit_fac << " " << ul << endl;
    
    funct fd=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_density_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
    funct fe=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_energy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
    funct fs=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_entropy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);

    if (verbose>1) {
      cout << "boson_rel::calc_mu() deg density integral" << endl;
    }

    int iret;
    double err;

    dit->err_nonconv=false;
    iret=dit->integ_err(fd,0.0,ul,b.n,err);
    if (iret!=0) {
      cout << "Problem 1." << endl;
      for(double xx=0.0;xx<ul*(1.01);xx+=ul/20.0) {
        cout << xx << " " << fd(xx) << endl;
      }
      exit(-1);
    }
    b.n*=b.g/2.0/pi2;
    
    if (verbose>1) {
      cout << "boson_rel::calc_mu() deg energy density integral" << endl;
    }
    iret=dit->integ_err(fe,0.0,ul,b.ed,err);
    if (iret!=0) {
      cout << "Problem 2." << endl;
      for(double xx=0.0;xx<ul*(1.01);xx+=ul/20.0) {
        cout << xx << " " << fe(xx) << endl;
      }
      exit(-1);
    }
    b.ed*=b.g/2.0/pi2;
    
    if (verbose>1) {
      cout << "boson_rel::calc_mu() deg entropy density integral" << endl;
    }
    iret=dit->integ_err(fs,0.0,ul,b.en,err);
    if (iret!=0) {
      cout << "Problem 3." << endl;
      for(double xx=0.0;xx<ul*(1.01);xx+=ul/20.0) {
        cout << xx << " " << fs(xx) << endl;
      }
      exit(-1);
    }
    b.en*=b.g/2.0/pi2;

    dit->err_nonconv=true;
    
  } else {
    
    // If the temperature is large enough, perform the full integral
    
    funct mfd=std::bind(std::mem_fn<double(double,boson &,double)>
			  (&boson_rel::density_fun),
			  this,std::placeholders::_1,std::ref(b),temper);
    funct mfe=std::bind(std::mem_fn<double(double,boson &,double)>
			  (&boson_rel::energy_fun),
			  this,std::placeholders::_1,std::ref(b),temper);
    funct mfs=std::bind(std::mem_fn<double(double,boson &,double)>
			  (&boson_rel::entropy_fun),
			  this,std::placeholders::_1,std::ref(b),temper);
      
    double prefac=b.g*pow(temper,3.0)/2.0/pi2;

    // Compute the number density

    double err;
    int iret;

    nit->err_nonconv=false;
    
    if (verbose>1) {
      cout << "boson_rel::calc_mu() ndeg density integral" << endl;
    }
    iret=nit->integ_err(mfd,0.0,0.0,b.n,err);
    if (iret!=0 || b.n==0.0) {
      cout << "Problem 4." << endl;
      exit(-1);
    }
    b.n*=prefac;

    // Compute the energy density

    if (verbose>1) {
      cout << "boson_rel::calc_mu() ndeg energy density integral" << endl;
    }
    iret=nit->integ_err(mfe,0.0,0.0,b.ed,err);
    if (iret!=0) {
      cout << "Problem 5." << endl;
      exit(-1);
    }
    b.ed*=prefac*temper;
    if (!b.inc_rest_mass) b.ed-=b.n*b.m;
    
    // Compute the entropy

    if (verbose>1) {
      cout << "boson_rel::calc_mu() ndeg entropy density integral" << endl;
    }
    iret=nit->integ_err(mfs,0.0,0.0,b.en,err);
    if (iret!=0) {
      cout << "Problem 6." << endl;
      exit(-1);
    }
    b.en*=prefac;

    nit->err_nonconv=true;
    
  }
  
  b.pr=-b.ed+temper*b.en+b.mu*b.n;

  return;
}

void boson_rel::nu_from_n(boson &b, double temper) {
  
  ubvector x(1);

  x[0]=b.nu/temper;
  
  mm_funct mf=std::bind(std::mem_fn<int(size_t nv, const ubvector &,
                                        ubvector &,double,boson &,double)>
                        (&boson_rel::solve_fun),
                        this,std::placeholders::_1,std::placeholders::_2,
                        std::placeholders::_3,b.n,std::ref(b),temper);
  
  bool ec=density_mroot->err_nonconv;
  density_mroot->err_nonconv=false;
  int ret1=density_mroot->msolve(1,x,mf);
  density_mroot->err_nonconv=ec;

  if (ret1!=0) {
    density_mroot->verbose=2;
    int ret1=density_mroot->msolve(1,x,mf);
  }

  /*
  if (ret1!=0) {

    root_brent_gsl<> rbg;
    rbg.err_nonconv=false;
    int ret2=rbg.solve(nex,mf);
  */

  if (ret1!=0) {
    O2SCL_ERR("Solvers failed in boson_rel::nu_from_n().",
              o2scl::exc_efailed);
  }
  //}
  
  b.nu=x[0]*temper;
  
  return;
}

void boson_rel::calc_density(boson &b, double temper) {

  if (temper<=0.0) {
    O2SCL_ERR2("Temperature less than or equal to zero in ",
	       "boson_rel::calc_density().",exc_einval);
  }
  if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }

  nu_from_n(b,temper);

  calc_mu(b,temper);

  /*
  funct fe=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_energy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);
  funct fs=std::bind(std::mem_fn<double(double,boson &,double)>
		       (&boson_rel::deg_entropy_fun),
		       this,std::placeholders::_1,std::ref(b),temper);

  b.ed=dit->integ(fe,0.0,sqrt(pow(20.0*temper+b.nu,2.0)-b.ms*b.ms));
  b.ed*=b.g/2.0/pi2;
  b.en=dit->integ(fs,0.0,sqrt(pow(20.0*temper+b.nu,2.0)-b.ms*b.ms));
  b.en*=b.g/2.0/pi2;

  b.pr=-b.ed+temper*b.en+b.mu*b.n;
  */

  return;
}

double boson_rel::deg_density_fun(double k, boson &b, double T) {

  double nx;
  if (b.ms==b.nu) {
    nx=1.0/(exp(k/T)-1.0);
  } else {
    double E=o2hypot(k,b.ms);
    nx=o2scl::bose_function(E,b.nu,T);
  }
  double ret=k*k*nx;

  if (nx<0.0) {
    long double m2=b.ms;
    long double mu2=b.nu;
    long double T2=T;
    long double k2=k;
    long double E2=o2hypot(k2,m2);
    long double nx2=1.0L/(exp((E2-mu2)/T2)-1.0L);
    cout << b.ms << " " << b.nu << " " << k << " " << T << endl;
    cout << E2-mu2 << " " << b.ms-b.nu << " " << nx << " "
         << nx2 << endl;
    cout << "Problem 17." << endl;
    exit(-1);
  }
  if (!std::isfinite(ret)) {
    return 0.0;
    /*
      cout << "1: " << k << " " << b.ms << " " << b.nu << " " << T << endl;
      cout << exp(E/T-b.nu/T)-1.0 << " " << E/T-b.nu/T << endl;
      cout << b.nu-b.ms << endl;
      exit(-1);
    */
  }
  
  return ret;
}
  
double boson_rel::deg_energy_fun(double k, boson &b, double T) {

  double E=o2hypot(k,b.ms);
  double nx=o2scl::bose_function(E,b.nu,T);
  double ret=k*k*E*nx;
  
  if (!std::isfinite(ret)) {
    return 0.0;
    //cout << "2: " << k << " " << b.ms << " " << b.nu << " " << T << endl;
    //exit(-1);
  }

  return ret;
}
  
double boson_rel::deg_entropy_fun(double k, boson &b, double T) {

  double E=o2hypot(k,b.ms);
  double nx=o2scl::bose_function(E,b.nu,T);
  double ret;
  ret=-k*k*(nx*log(nx)-(1.0+nx)*log(1.0+nx));
  
  if (!std::isfinite(ret)) {
    return 0.0;
    /*
    double psi;
    if (b.inc_rest_mass) {
      psi=(b.nu-b.ms)/T;
    } else {
      psi=(b.nu+(b.m-b.ms))/T;
    }
    cout << "3: " << k << " " << b.ms << " " << b.nu << " " << T << endl;
    cout << "psi: " << psi << endl;
    cout << exp(E/T-b.nu/T)-1.0 << " " << E/T-b.nu/T << endl;
    cout << b.nu-b.ms << " " << nx << endl;
    cout << ret << endl;
    exit(-1);
    */
  }

  return ret;
}
  
double boson_rel::density_fun(double u, boson &b, double T) {

  double y;
  if (b.inc_rest_mass) {
    y=b.nu/T;
  } else {
    y=(b.nu+b.m)/T;
  }
  double eta=b.ms/T;

  double ret;
  if (y>700.0 && eta+u>700.0) {
    if (eta+u-y>700.0) {
      ret=0.0;
    } else {
      ret=(eta+u)*sqrt(u*u+2.0*eta*u)/(exp(eta+u-y)-1.0);
    }
  } else {
    ret=(eta+u)*sqrt(u*u+2.0*eta*u)*exp(y)/(exp(eta+u)-exp(y));
  }

  if (!std::isfinite(ret)) {
    cout << "4: " << u << " " << y << " " << eta << " " 
	 << b.ms << " " << b.nu << " " << T << endl;
    exit(-1);
  }

  return ret;
}

double boson_rel::energy_fun(double u, boson &b, double T) {

  double y;
  if (b.inc_rest_mass) {
    y=b.nu/T;
  } else {
    y=(b.nu+b.m)/T;
  }
  double eta=b.ms/T;

  double ret;
  if (y-u>200.0 && eta-u>200.0) {
    if (eta+u+y>100.0) {
      ret=0.0;
    } else {
      ret=(eta+u)*(eta+u)*sqrt(u*u+2.0*eta*u)/(exp(eta+u-y)-1.0);
    }
  } else {
    ret=(eta+u)*(eta+u)*sqrt(u*u+2.0*eta*u)*exp(y)/(exp(eta+u)-exp(y));
  }
  
  if (!std::isfinite(ret)) {
    cout << "5: " << u << " " << b.ms << " " << b.nu << " " << T << endl;
    exit(-1);
  }

  return ret;
}

double boson_rel::entropy_fun(double u, boson &b, double T) {

  double y;
  if (b.inc_rest_mass) {
    y=b.nu/T;
  } else {
    y=(b.nu+b.m)/T;
  }
  double eta=b.ms/T;

  double arg1=u*u+2*eta*u;
  double arg2=eta+u-y;
  double arg3=eta+u;

  double fb=1.0/(-1.0+exp(arg2));
  double ret=arg3*sqrt(arg1)*((1.0+fb)*log(1.0+fb)-fb*log(fb));

  if (!std::isfinite(ret)) {
    return 0.0;
  }
  /*
  double arg4=y-eta-u;
  double arg5=1+exp(arg4);
  double arg6=1+exp(arg2);
  double term1=log(arg5)/arg5;
  double term2=log(arg6)/arg6;
  double ret=arg3*sqrt(arg1)*(term1+term2);
  return ret;

  double ret;
  if (u-eta>200.0 && u-y>200.0) {
    ret=0.0;
  } else {
    double term1=exp(eta+u)*log(1.0/1.0-exp(y-eta-u));
    double term2=exp(y)*log(1.0/(exp(eta+u-y)-1.0));
    ret=(eta+u)*sqrt(u*u+2.0*eta*u)*(term1+term2)/
      (exp(eta+u)-exp(y));
  }

  */

  /*
  if (false) {
    cout << "6: " << u << " " << eta << " " << y << endl;
    cout << b.ms << " " << b.nu << " " << T << endl;

    u=200.0;
    term1=exp(eta+u)*log(1.0/1.0-exp(y-eta-u));
    term2=exp(y)*log(1.0/(exp(eta+u-y)-1.0));
    ret=(eta+u)*sqrt(u*u+2.0*eta*u)*(term1+term2)/
      (exp(eta+u)-exp(y));
    cout << ret << endl;
    
    exit(-1);
  }
  */

  return ret;
}

int boson_rel::solve_fun(size_t nv, const ubvector &x, ubvector &y,
                         double density, boson &b, double T) {

  double nden;
  
  b.nu=x[0]*T;
  if (b.non_interacting) b.mu=b.nu;
  
  double psi;
  if (b.inc_rest_mass) {
    psi=(b.nu-b.ms)/T;
  } else {
    psi=(b.nu+(b.m-b.ms))/T;
  }

  if (b.nu>b.ms) return 1;

  bool deg=true;
  if (psi<deg_limit) deg=false;
  
  if (deg) {

      // Compute the upper limit for degenerate integrals

    double arg, upper_limit_fac=30.0;
    if (b.inc_rest_mass) {
      arg=pow(upper_limit_fac*T+b.nu,2.0)-b.ms*b.ms;
    } else {
      arg=pow(upper_limit_fac*T+b.nu+b.m,2.0)-b.ms*b.ms;
    }
    double ul=sqrt(arg);
    
    funct fd=std::bind(std::mem_fn<double(double,boson &,double)>
                       (&boson_rel::deg_density_fun),
                       this,std::placeholders::_1,std::ref(b),T);
    double err;
    dit->err_nonconv=false;
    int iret=dit->integ_err(fd,0.0,ul,nden,err);
    dit->err_nonconv=true;
    if (iret!=0) {
      /*
        table_units<> t;
        def_dit.get_workspace().make_table(t);
        o2scl_hdf::hdf_file hf;
        hf.open_or_create("br.o2");
        hdf_output(hf,t,"br");
        hf.close();
        cout << b.nu << " " << b.ms << endl;
      */
      cout << "Problem 7b." << endl;
      exit(-1);
    }
    nden*=b.g/2.0/pi2;

  } else {

    // If the temperature is large enough, perform the full integral
    
    funct mfd=std::bind(std::mem_fn<double(double,boson &,double)>
			(&boson_rel::density_fun),
			this,std::placeholders::_1,std::ref(b),T);
    
    double prefac=b.g*pow(T,3.0)/2.0/pi2;
    
    // Compute the number density

    double err;
    nit->err_nonconv=false;
    int iret=nit->integ_err(mfd,0.0,0.0,nden,err);
    nit->err_nonconv=true;
    if (iret!=0) {
      cout << "Problem 8." << endl;
      exit(-1);
    }
    nden*=prefac;

  }

  y[0]=nden/density-1.0;
  cout.precision(12);
  cout << "A: " << x[0] << " " << y[0] << endl;

  return 0;
}

void boson_rel::pair_mu(boson &b, double temper) {
  
  if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }
  calc_mu(b,temper);
  
  boson antip(b.ms,b.g);
  b.anti(antip);
  calc_mu(antip,temper);
  b.n-=antip.n;
  b.pr+=antip.pr;
  b.ed+=antip.ed;
  b.en+=antip.en;

  return;
}

void boson_rel::pair_density(boson &b, double temper) {
  
  if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }

  ubvector x(1);
  x[0]=b.nu/temper;

  mm_funct mf=std::bind(std::mem_fn<int(size_t nv, const ubvector &,
                                        ubvector &,double,boson &,double)>
                        (&boson_rel::pair_density_fun),
                        this,std::placeholders::_1,std::placeholders::_2,
                        std::placeholders::_3,b.n,std::ref(b),temper);
  bool ec=density_mroot->err_nonconv;
  density_mroot->err_nonconv=false;
  int ret1=density_mroot->msolve(1,x,mf);
  density_mroot->err_nonconv=ec;

  /*
  if (ret1!=0) {

    root_brent_gsl<> rbg;
    rbg.err_nonconv=false;
    int ret2=rbg.solve(x,mf);
  */

  if (ret1!=0) {
    O2SCL_ERR("Solvers failed in boson_rel::nu_from_n().",
              o2scl::exc_efailed);
  }
  //}
  
  b.nu=x[0]*temper;

  if (b.non_interacting==true) { b.mu=b.nu; }
  
  pair_mu(b,temper);

  return;
}

int boson_rel::pair_density_fun(size_t nv,
                                const ubvector &x, ubvector &y,
                                double density, boson &b, double T) {

  b.nu=x[0]*T;
  if (b.non_interacting) {
    b.mu=b.nu;
  }

  pair_mu(b,T);
  
  y[0]=(b.n-density)/density;
  
  cout << "H: " << x[0] << " " << y[0] << " " << b.nu << " " << b.ms
       << endl;
  cout << "\t: " << b.n << " " << density << endl;
  
  return 0;
}
