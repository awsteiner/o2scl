/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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
#ifndef O2SCL_BOSON_REL_H
#define O2SCL_BOSON_REL_H

/** \file boson_rel.h
    \brief File defining \ref o2scl::boson_rel_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/inte.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>

#include <o2scl/boson.h>

namespace o2scl {

  /** \brief Equation of state for a relativistic boson
      
      \verbatim embed:rst

      .. todo:: 

         - In class boson_rel: Testing not completely finished.
         
      \endverbatim
  */
  template<class be_inte_t=bessel_K_exp_integ_gsl, class fp_t=double>
  class boson_rel_tl : public boson_thermo_tl<be_inte_t,fp_t> {

  public:

    /// Internal value of pi squared for convenience
    fp_t pi2;
    
    /// Create a boson with mass \c m and degeneracy \c g
    boson_rel_tl() {
      pi2=boost::math::constants::pi_sqr<fp_t>();
      density_mroot=&def_density_mroot;
      nit=&def_nit;
      dit=&def_dit;
      verify_ti=false;
      use_expansions=true;
      deg_limit=-0.5;
      upper_limit_fac=30.0;
      verbose=0;
    }
    
    virtual ~boson_rel_tl() {
    }
    
    /** \brief Calculate properties as function of chemical potential
     */
    virtual void calc_mu(boson &b, fp_t temper) {

      if (temper<=0) {
        O2SCL_ERR2("Temperature less than or equal to zero in ",
                   "boson_rel_tl<be_inte_t,fp_t>::calc_mu().",exc_einval);
      }
      if (b.non_interacting==true) {
        b.nu=b.mu;
        b.ms=b.m;
      }

      fp_t psi;
      if (b.inc_rest_mass) {
        psi=(b.nu-b.ms)/temper;
        if (b.nu>b.ms) {
          O2SCL_ERR2("Chemical potential must be smaller than mass in ",
                     "boson_rel_tl<be_inte_t,fp_t>::calc_mu().",
                     o2scl::exc_einval);
        }
      } else {
        psi=(b.nu+(b.m-b.ms))/temper;
        if (b.nu+b.m>b.ms) {
          std::cout.precision(12);
          std::cout << "Here: " << b.nu << " " << b.m << " " << b.ms << " "
                    << b.ms/b.m << std::endl;
          O2SCL_ERR2("Chemical potential must be smaller than mass in ",
                     "boson_rel_tl<be_inte_t,fp_t>::calc_mu().",
                     o2scl::exc_einval);
        }
      }

      // Try the non-degenerate expansion if psi is small enough
      if (use_expansions) {
        bool acc=this->calc_mu_ndeg(b,temper,1.0e-14);
        if (verbose>1) {
          std::cout << "calc_mu(): non-deg expan (bosons) " << acc
                    << " " << verbose << std::endl;
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
        std::cout << "boson_rel_tl<be_inte_t,fp_t>::calc_mu() psi: "
                  << psi << std::endl;
      }
  
      if (deg) {
    
        // Compute the upper limit for degenerate integrals

        fp_t arg;
        if (b.inc_rest_mass) {
          arg=pow(upper_limit_fac*temper+b.nu,2.0)-b.ms*b.ms;
        } else {
          arg=pow(upper_limit_fac*temper+b.nu+b.m,2.0)-b.ms*b.ms;
        }
        fp_t ul;
        if (arg>0) {
          ul=sqrt(arg);
        } else {
          O2SCL_ERR2("Zero density in degenerate limit in ",
                     "boson_rel_tl<be_inte_t,fp_t>::calc_mu().",
                     exc_efailed);
          return;
        }

        funct fd=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
                           (&boson_rel_tl<be_inte_t,fp_t>::deg_density_fun),
                           this,std::placeholders::_1,std::ref(b),temper);
        funct fe=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
                           (&boson_rel_tl<be_inte_t,fp_t>::deg_energy_fun),
                           this,std::placeholders::_1,std::ref(b),temper);
        funct fs=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
                           (&boson_rel_tl<be_inte_t,fp_t>::deg_entropy_fun),
                           this,std::placeholders::_1,std::ref(b),temper);

        if (verbose>1) {
          std::cout << "boson_rel_tl<be_inte_t,fp_t>::calc_mu() "
                    << "deg density integral" << std::endl;
        }

        int iret;
        fp_t err;

        dit->err_nonconv=false;
        iret=dit->integ_err(fd,0,ul,b.n,err);
        if (iret!=0) {
          std::cout << "Problem 1." << std::endl;
          for(fp_t xx=0;xx<ul*(1.01);xx+=ul/20) {
            std::cout << xx << " " << fd(xx) << std::endl;
          }
          exit(-1);
        }
        b.n*=b.g/2.0/pi2;
    
        if (verbose>1) {
          std::cout << "boson_rel_tl<be_inte_t,fp_t>::calc_mu() deg "
                    << "energy density integral" << std::endl;
        }
        iret=dit->integ_err(fe,0,ul,b.ed,err);
        if (iret!=0) {
          std::cout << "Problem 2." << std::endl;
          for(fp_t xx=0;xx<ul*(1.01);xx+=ul/20) {
            std::cout << xx << " " << fe(xx) << std::endl;
          }
          exit(-1);
        }
        b.ed*=b.g/2.0/pi2;
    
        if (verbose>1) {
          std::cout << "boson_rel_tl<be_inte_t,fp_t>::calc_mu() deg "
                    << "entropy density integral" << std::endl;
        }
        iret=dit->integ_err(fs,0,ul,b.en,err);
        if (iret!=0) {
          std::cout << "Problem 3." << std::endl;
          for(fp_t xx=0;xx<ul*(1.01);xx+=ul/20) {
            std::cout << xx << " " << fs(xx) << std::endl;
          }
          exit(-1);
        }
        b.en*=b.g/2.0/pi2;

        dit->err_nonconv=true;
    
      } else {
    
        // If the temperature is large enough, perform the full integral
    
        funct mfd=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
                            (&boson_rel_tl<be_inte_t,fp_t>::density_fun),
                            this,std::placeholders::_1,std::ref(b),temper);
        funct mfe=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
                            (&boson_rel_tl<be_inte_t,fp_t>::energy_fun),
                            this,std::placeholders::_1,std::ref(b),temper);
        funct mfs=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
                            (&boson_rel_tl<be_inte_t,fp_t>::entropy_fun),
                            this,std::placeholders::_1,std::ref(b),temper);
      
        fp_t prefac=b.g*pow(temper,3.0)/2.0/pi2;

        // Compute the number density

        fp_t err;
        int iret;

        nit->err_nonconv=false;
    
        if (verbose>1) {
          std::cout << "boson_rel_tl<be_inte_t,fp_t>::calc_mu() ndeg "
                    << "density integral" << std::endl;
        }
        iret=nit->integ_err(mfd,0,0,b.n,err);
        if (iret!=0 || b.n==0) {
          std::cout << "Problem 4." << std::endl;
          exit(-1);
        }
        b.n*=prefac;

        // Compute the energy density

        if (verbose>1) {
          std::cout << "boson_rel_tl<be_inte_t,fp_t>::calc_mu() ndeg "
                    << "energy density integral" << std::endl;
        }
        iret=nit->integ_err(mfe,0,0,b.ed,err);
        if (iret!=0) {
          std::cout << "Problem 5." << std::endl;
          exit(-1);
        }
        b.ed*=prefac*temper;
        if (!b.inc_rest_mass) b.ed-=b.n*b.m;
    
        // Compute the entropy

        if (verbose>1) {
          std::cout << "boson_rel_tl<be_inte_t,fp_t>::calc_mu() ndeg "
                    << "entropy density integral" << std::endl;
        }
        iret=nit->integ_err(mfs,0,0,b.en,err);
        if (iret!=0) {
          std::cout << "Problem 6." << std::endl;
          exit(-1);
        }
        b.en*=prefac;

        nit->err_nonconv=true;
    
      }
  
      b.pr=-b.ed+temper*b.en+b.mu*b.n;

      return;
    }
    
    /** \brief Calculate properties as function of density
     */
    virtual void calc_density(boson &b, fp_t temper) {

      if (temper<=0) {
        O2SCL_ERR2("Temperature less than or equal to zero in ",
                   "boson_rel_tl<be_inte_t,fp_t>::calc_density().",
                   exc_einval);
      }
      if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }

      if (b.non_interacting) {
        fp_t mu_temp=b.mu;
        fp_t n_temp=b.n;
        calc_max_density(b,temper);
        if (b.n<n_temp) {
          std::cout << "m: " << b.m << " mu: " << mu_temp << std::endl;
          std::cout.precision(12);
          std::cout << "requested density: " << n_temp << " max density: "
               << b.n << std::endl;
          O2SCL_ERR2("Density larger than max in ",
                     "boson_rel_tl<be_inte_t,fp_t>::calc_density().",
                     o2scl::exc_einval);
        }
        b.mu=mu_temp;
        b.n=n_temp;
      } else {
        fp_t mu_temp=b.nu;
        fp_t n_temp=b.n;
        calc_max_density(b,temper);
        if (b.n<n_temp) {
          O2SCL_ERR2("Density larger than max in ",
                     "boson_rel_tl<be_inte_t,fp_t>::calc_density().",
                     o2scl::exc_einval);
        }
        b.nu=mu_temp;
        b.n=n_temp;
      }
  
      nu_from_n(b,temper);

      calc_mu(b,temper);

      //cout << "Here: " << b.n << " " << b.nu << std::endl;
      //exit(-1);
      /*
        funct fe=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
        (&boson_rel_tl<be_inte_t,fp_t>::deg_energy_fun),
        this,std::placeholders::_1,std::ref(b),temper);
        funct fs=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
        (&boson_rel_tl<be_inte_t,fp_t>::deg_entropy_fun),
        this,std::placeholders::_1,std::ref(b),temper);

        b.ed=dit->integ(fe,0,sqrt(pow(20*temper+b.nu,2.0)-b.ms*b.ms));
        b.ed*=b.g/2.0/pi2;
        b.en=dit->integ(fs,0,sqrt(pow(20*temper+b.nu,2.0)-b.ms*b.ms));
        b.en*=b.g/2.0/pi2;

        b.pr=-b.ed+temper*b.en+b.mu*b.n;
      */

      return;
    }
    
    /** \brief Calculate the maximum density as a function of temperature
     */
    virtual void calc_max_density(boson &b, fp_t temper) {

      if (temper<=0) {
        O2SCL_ERR2("Temperature less than or equal to zero in ",
                   "boson_rel_tl<be_inte_t,fp_t>::calc_mu().",exc_einval);
      }
  
      if (b.non_interacting==true) {
        b.mu=b.m;
      } else {
        b.nu=b.ms;
      }

      return calc_mu(b,temper);
    }      
    
    /** \brief Calculate properties with antiparticles as function of
        chemical potential
    */
    virtual void pair_mu(boson &b, fp_t temper) {
  
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


    /** \brief Calculate properties with antiparticles as function of
	density
    */
    virtual void pair_density(boson &b, fp_t temper) {
  
      if (b.non_interacting==true) { b.nu=b.mu; b.ms=b.m; }

      ubvector x(1);
      x[0]=b.nu/temper;

      mm_funct mf=std::bind
        (std::mem_fn<int(size_t nv, const ubvector &,
                         ubvector &,fp_t,boson &,fp_t)>
         (&boson_rel_tl<be_inte_t,fp_t>::pair_density_fun),
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
        O2SCL_ERR2("Solvers failed in boson_rel_tl<be_inte_t,fp_t>",
                   "::nu_from_n().",o2scl::exc_efailed);
      }
      //}
  
      b.nu=x[0]*temper;

      if (b.non_interacting==true) { b.mu=b.nu; }
  
      pair_mu(b,temper);

      return;
    }


    /// Calculate effective chemical potential from density
    virtual void nu_from_n(boson &b, fp_t temper) {
  
      ubvector x(1);

      if (b.non_interacting==true) {
        b.nu=b.mu;
        b.ms=b.m;
      }

      // If the chemical potential is too large, create a valid
      // initial guess
      if (b.nu>=b.ms) {
        b.nu=b.ms-1.0e-4;
      }
  
      x[0]=b.nu/temper;
  
      mm_funct mf=std::bind(std::mem_fn<int(size_t nv, const ubvector &,
                                            ubvector &,fp_t,boson &,fp_t)>
                            (&boson_rel_tl<be_inte_t,fp_t>::solve_fun),
                            this,std::placeholders::_1,std::placeholders::_2,
                            std::placeholders::_3,b.n,std::ref(b),temper);

      bool ec=density_mroot->err_nonconv;
      density_mroot->err_nonconv=false;
      int ret1=density_mroot->msolve(1,x,mf);
      density_mroot->err_nonconv=ec;

      if (ret1!=0) {
        std::cout << "Problem A." << std::endl;
        std::cout << b.nu << " " << b.ms << std::endl;
        ubvector y(1);
        for(x[0]=9.88e-2;x[0]<9.91e-2;x[0]+=1.0e-6) {
          mf(1,x,y);
          std::cout << x[0] << " " << y[0] << std::endl;
        }
        exit(-1);
        density_mroot->verbose=2;
        x[0]*=0.99;
        ret1=density_mroot->msolve(1,x,mf);
      }

      if (ret1!=0) {
        O2SCL_ERR2("Solvers failed in ",
                   "boson_rel_tl<be_inte_t,fp_t>::nu_from_n().",
                  o2scl::exc_efailed);
      }
      //}
  
      b.nu=x[0]*temper;
  
      return;
    }
    
    /// Set degenerate and nondegenerate integrators
    void set_inte(inte<> &l_nit, inte<> &l_dit) {
      nit=&l_nit;
      dit=&l_dit;
      return;
    }      

    /** \brief Set the solver for use in calculating the chemical
	potential from the density */
    void set_density_mroot(mroot<> &rp) {
      density_mroot=&rp;
      return;
    }

    /// The default solver for calc_density().
    mroot_hybrids<> def_density_mroot;
    
    /// Default nondegenerate integrator
    inte_qagiu_gsl<> def_nit;

    /// Default degenerate integrator
    inte_qag_gsl<> def_dit;

    /// Return string denoting type ("boson_rel")
    virtual const char *type() { return "boson_rel"; }

    /** \brief If true, verify the thermodynamic identity
     */
    bool verify_ti;

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief Verbosity parameter
     */
    bool use_expansions;

    fp_t deg_limit;

    fp_t upper_limit_fac;
    
  protected:

    /// The non-degenerate integrator
    inte<> *nit;
    /// The degenerate integrator
    inte<> *dit;
    /// The solver for calc_density()
    mroot<> *density_mroot;

    /// Non-degenerate density integral
    fp_t density_fun(fp_t u, boson &b, fp_t T) {

      fp_t y;
      if (b.inc_rest_mass) {
        y=b.nu/T;
      } else {
        y=(b.nu+b.m)/T;
      }
      fp_t eta=b.ms/T;

      fp_t ret;
      if (y>700 && eta+u>700) {
        if (eta+u-y>700) {
          ret=0;
        } else {
          ret=(eta+u)*sqrt(u*u+2.0*eta*u)/(exp(eta+u-y)-1);
        }
      } else {
        ret=(eta+u)*sqrt(u*u+2.0*eta*u)*exp(y)/(exp(eta+u)-exp(y));
      }

      if (!std::isfinite(ret)) {
        std::cout << "4: " << u << " " << y << " " << eta << " " 
             << b.ms << " " << b.nu << " " << T << std::endl;
        exit(-1);
      }

      return ret;
    }

    /// Non-degenerate energy density integral
    fp_t energy_fun(fp_t u, boson &b, fp_t T) {

      fp_t y;
      if (b.inc_rest_mass) {
        y=b.nu/T;
      } else {
        y=(b.nu+b.m)/T;
      }
      fp_t eta=b.ms/T;

      fp_t ret;
      if (y-u>200 && eta-u>200) {
        if (eta+u+y>100) {
          ret=0;
        } else {
          ret=(eta+u)*(eta+u)*sqrt(u*u+2.0*eta*u)/(exp(eta+u-y)-1);
        }
      } else {
        ret=(eta+u)*(eta+u)*sqrt(u*u+2.0*eta*u)*exp(y)/(exp(eta+u)-exp(y));
      }
  
      if (!std::isfinite(ret)) {
        std::cout << "5: " << u << " " << b.ms << " " << b.nu
                  << " " << T << std::endl;
        exit(-1);
      }

      return ret;
    }

    /// Non-degenerate entropy integral
    fp_t entropy_fun(fp_t u, boson &b, fp_t T) {

      fp_t y;
      if (b.inc_rest_mass) {
        y=b.nu/T;
      } else {
        y=(b.nu+b.m)/T;
      }
      fp_t eta=b.ms/T;

      fp_t arg1=u*u+2*eta*u;
      fp_t arg2=eta+u-y;
      fp_t arg3=eta+u;

      fp_t fb=1/(-1+exp(arg2));
      fp_t ret=arg3*sqrt(arg1)*((1+fb)*log(1+fb)-fb*log(fb));

      if (!std::isfinite(ret)) {
        return 0;
      }
      /*
        fp_t arg4=y-eta-u;
        fp_t arg5=1+exp(arg4);
        fp_t arg6=1+exp(arg2);
        fp_t term1=log(arg5)/arg5;
        fp_t term2=log(arg6)/arg6;
        fp_t ret=arg3*sqrt(arg1)*(term1+term2);
        return ret;

        fp_t ret;
        if (u-eta>200 && u-y>200) {
        ret=0;
        } else {
        fp_t term1=exp(eta+u)*log(1/1-exp(y-eta-u));
        fp_t term2=exp(y)*log(1/(exp(eta+u-y)-1));
        ret=(eta+u)*sqrt(u*u+2.0*eta*u)*(term1+term2)/
        (exp(eta+u)-exp(y));
        }

      */

      /*
        if (false) {
        std::cout << "6: " << u << " " << eta << " " << y << std::endl;
        std::cout << b.ms << " " << b.nu << " " << T << std::endl;

        u=200;
        term1=exp(eta+u)*log(1/1-exp(y-eta-u));
        term2=exp(y)*log(1/(exp(eta+u-y)-1));
        ret=(eta+u)*sqrt(u*u+2.0*eta*u)*(term1+term2)/
        (exp(eta+u)-exp(y));
        std::cout << ret << std::endl;
    
        exit(-1);
        }
      */

      return ret;
    }

      
    /// Degenerate density integral
    fp_t deg_density_fun(fp_t k, boson &b, fp_t T) {

      fp_t E=hypot(k,b.ms);
      fp_t nx=o2scl::bose_function((E-b.nu)/T);
      fp_t ret=k*k*nx;

      if (!std::isfinite(ret)) {
        return 0;
        /*
          std::cout << "1: " << k << " " << b.ms << " " 
          << b.nu << " " << T << std::endl;
          std::cout << exp(E/T-b.nu/T)-1 << " " << E/T-b.nu/T << std::endl;
          std::cout << b.nu-b.ms << std::endl;
          exit(-1);
        */
      }
  
      return ret;
    }
    
    /// Degenerate energy density integral
    fp_t deg_energy_fun(fp_t k, boson &b, fp_t T) {

      fp_t E=hypot(k,b.ms);
      fp_t nx=o2scl::bose_function((E-b.nu)/T);
      fp_t ret=k*k*E*nx;
  
      if (!std::isfinite(ret)) {
        return 0;
        //std::cout << "2: " << k << " " << b.ms << " " << b.nu
        //<< " " << T << std::endl;
        //exit(-1);
      }

      return ret;
    }
      
    /// Degenerate entropy integral
    fp_t deg_entropy_fun(fp_t k, boson &b, fp_t T) {

      fp_t E=hypot(k,b.ms);
      fp_t nx=o2scl::bose_function((E-b.nu)/T);
      fp_t ret;
      ret=-k*k*(nx*log(nx)-(1+nx)*log(1+nx));
  
      if (!std::isfinite(ret)) {
        return 0;
        /*
          fp_t psi;
          if (b.inc_rest_mass) {
          psi=(b.nu-b.ms)/T;
          } else {
          psi=(b.nu+(b.m-b.ms))/T;
          }
          std::cout << "3: " << k << " " << b.ms << " " << b.nu 
          << " " << T << std::endl;
          std::cout << "psi: " << psi << std::endl;
          std::cout << exp(E/T-b.nu/T)-1 << " " << E/T-b.nu/T << std::endl;
          std::cout << b.nu-b.ms << " " << nx << std::endl;
          std::cout << ret << std::endl;
          exit(-1);
        */
      }

      return ret;
    }

    /// Solve for the density in calc_density()
    int solve_fun(size_t nv, const ubvector &x, ubvector &y,
                  fp_t density, boson &b, fp_t T) {

      fp_t nden;
  
      b.nu=x[0]*T;
      if (b.non_interacting) b.mu=b.nu;
  
      fp_t psi;
      if (b.inc_rest_mass) {
        psi=(b.nu-b.ms)/T;
      } else {
        psi=(b.nu+(b.m-b.ms))/T;
      }

      if (b.nu>b.ms) {
        return 1;
      }

      bool deg=true;
      if (psi<deg_limit) deg=false;

      if (deg) {

        // Compute the upper limit for degenerate integrals

        fp_t arg;
        if (b.inc_rest_mass) {
          arg=pow(upper_limit_fac*T+b.nu,2.0)-b.ms*b.ms;
        } else {
          arg=pow(upper_limit_fac*T+b.nu+b.m,2.0)-b.ms*b.ms;
        }
        fp_t ul=sqrt(arg);
    
        funct fd=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
                           (&boson_rel_tl<be_inte_t,fp_t>::deg_density_fun),
                           this,std::placeholders::_1,std::ref(b),T);
        fp_t err;
        dit->err_nonconv=false;
    
        int iret=dit->integ_err(fd,0,ul,nden,err);
        if (fabs(err)/fabs(nden)>1.0e-4) {
          def_dit.tol_rel=1.0e-10;
          def_dit.tol_abs=1.0e-10;
          iret=dit->integ_err(fd,0,ul,nden,err);
          def_dit.tol_rel=1.0e-8;
          def_dit.tol_abs=1.0e-8;
        }
        
        dit->err_nonconv=true;

        if (iret!=0) {
          /*
            table_units<> t;
            def_dit.get_workspace().make_table(t);
            o2scl_hdf::hdf_file hf;
            hf.open_or_create("br.o2");
            hdf_output(hf,t,"br");
            hf.close();
            std::cout << b.nu << " " << b.ms << std::endl;
          */
          std::cout << "Problem 7b." << std::endl;
          exit(-1);
        }
        nden*=b.g/2.0/pi2;

      } else {

        // If the temperature is large enough, perform the full integral
    
        funct mfd=std::bind(std::mem_fn<fp_t(fp_t,boson &,fp_t)>
                            (&boson_rel_tl<be_inte_t,fp_t>::density_fun),
                            this,std::placeholders::_1,std::ref(b),T);
    
        fp_t prefac=b.g*pow(T,3.0)/2.0/pi2;
    
        // Compute the number density

        fp_t err;
        nit->err_nonconv=false;
        int iret=nit->integ_err(mfd,0,0,nden,err);
        nit->err_nonconv=true;
        if (iret!=0) {
          std::cout << "Problem 8." << std::endl;
          exit(-1);
        }
        nden*=prefac;

      }

      y[0]=nden/density-1;

      return 0;
    }
      
    /// Solve for the density in pair_density()
    int pair_density_fun(size_t nv, const ubvector &x, ubvector &y,
                         fp_t density, boson &b, fp_t T) {

      b.nu=x[0]*T;
      if (b.non_interacting) {
        b.mu=b.nu;
      }

      pair_mu(b,T);
  
      y[0]=(b.n-density)/density;
  
      std::cout << "H: " << x[0] << " " << y[0] << " " << b.nu << " "
                << b.ms << std::endl;
      std::cout << "\t: " << b.n << " " << density << std::endl;
  
      return 0;
    }

  };

  /** \brief Double precision version of \ref o2scl::boson_rel_tl
  */
  typedef boson_rel_tl<> boson_rel;

  /** \brief Long double version of 
      \ref o2scl::boson_rel_tl 
  */
  typedef boson_rel_tl
  <bessel_K_exp_integ_boost<long double,
                            cpp_dec_float_25>,
   long double> boson_rel_ld;
  
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
  
  /** \brief 25-digit version of 
      \ref o2scl::boson_rel_tl 
  */
  typedef boson_rel_tl
  <bessel_K_exp_integ_boost<cpp_dec_float_25,
                            cpp_dec_float_35>,
   cpp_dec_float_25> boson_rel_cdf25;
  
#endif
  
}

#endif
