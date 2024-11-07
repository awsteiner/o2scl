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
#ifndef O2SCL_NONREL_FERMION_H
#define O2SCL_NONREL_FERMION_H

/** \file fermion_nonrel.h
    \brief File defining \ref o2scl::fermion_nonrel_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/root_cern.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/classical.h>
#include <o2scl/inte_qagiu_gsl.h>

#include <o2scl/fermion.h>
#include <o2scl/polylog.h>

namespace o2scl {

  /** \brief Nonrelativistic fermion class

      The effective mass computed by this class and stored in \ref
      part_tl::ms is the Landau mass, not the Dirac mass, as computed by
      \ref o2scl::fermion_rel_tl .

      This class works with both true and false values for either \ref
      part_tl::non_interacting or \ref part_tl::inc_rest_mass.

      Pressure is computed with
      \f[
      P = 2 \varepsilon/3
      \f]
      and entropy density with
      \f[
      s = \frac{5 \varepsilon}{3 T} - \frac{n \mu}{T}
      \f]
      \verbatim embed:rst
      These relations can be verified with an integration by
      parts. See, e.g. [Callen85]_ pg. 403 or [Landau80]_.
      \endverbatim

      The functions \ref pair_density() and \ref pair_mu() have not 
      been implemented and just call the error handler.
      
      \note The function \ref calc_density() calls the error handler
      at zero density and finite temperature, because the correct
      answer implies \f$ \mu = - \infty \f$ . At zero density and zero
      temperature the function \ref calc_density() calls \ref
      calc_density_zerot() which gives the proper chemical potential
      of \f$ mu = m \f$ without calling the error handler.

      \verbatim embed:rst

      .. todo::

      In class fermion_nonrel_tl

      - Future: Implement \ref
      o2scl::fermion_nonrel_tl::pair_density() and \ref
      o2scl::fermion_nonrel_tl::pair_mu(). [AWS, 1/23/19: it is not
      entirely clear to me that antiparticles will be useful.]

      - Future: This could be improved by performing a Chebyshev
      approximation (for example) to invert the density integral so
      that we don't need to use a solver.
         
      \endverbatim
  */
  template<class fermion_t=fermion_tl<double>,
	   class fd_inte_t=fermi_dirac_integ_gsl,
	   class be_inte_t=bessel_K_exp_integ_gsl,
	   class root_t=root_cern<>,
	   class func_t=funct, class fp_t=double>
  class fermion_nonrel_tl :
    public fermion_thermo_tl<fermion_t,fd_inte_t,be_inte_t,root_t,
			     func_t,fp_t> {
    
  public:
    
    /// Create a nonrelativistic fermion with mass 'm' and degeneracy 'g'
    fermion_nonrel_tl() {
      density_root=&def_density_root;
      verify_ti=false;
    }
    
    virtual ~fermion_nonrel_tl() {
    }

    /** \brief Zero temperature fermions
     */
    virtual void calc_mu_zerot(fermion_t &f) {
      if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
      if (f.inc_rest_mass) {
        f.kf=sqrt(2.0*f.ms*(f.nu-f.m));
      } else {
        f.kf=sqrt(2.0*f.ms*f.nu);
      }
      f.n=f.kf*f.kf*f.kf*f.g/6.0/o2scl_const::pi2;
      f.ed=f.g*pow(f.kf,5.0)/20.0/o2scl_const::pi2/f.ms;
      if (f.inc_rest_mass) f.ed+=f.n*f.m;
      f.pr=-f.ed+f.n*f.nu;
      f.en=0.0;
      return;
    }


    /** \brief Zero temperature fermions
     */
    virtual void calc_density_zerot(fermion_t &f) {
      if (f.non_interacting) { f.ms=f.m; }
      this->kf_from_density(f);
      f.nu=f.kf*f.kf/2.0/f.ms;
      f.ed=f.g*pow(f.kf,5.0)/20.0/o2scl_const::pi2/f.ms;
      if (f.inc_rest_mass) {
        f.ed+=f.n*f.m;
        f.nu+=f.m;
      }
      f.pr=-f.ed+f.n*f.nu;
      f.en=0.0;
  
      if (f.non_interacting) { f.mu=f.nu; }
      return;
    }
    
    /** \brief Calculate properties as function of chemical potential
     */
    virtual int calc_mu(fermion_t &f, fp_t temper) {

      if (temper<0.0) {
        O2SCL_ERR("Temperature less than zero in fermion_nonrel::calc_mu().",
                  exc_einval);
      }
      if (temper==0.0) {
        calc_mu_zerot(f);
        return 0;
      }

      fp_t y, sy, spi, ey, int1, int2;

      if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }

      if (f.ms<0.0) {
        O2SCL_ERR2("Mass negative in ",
                   "fermion_nonrel::calc_mu().",exc_einval);
      }

      if (temper<=0.0) {
        calc_mu_zerot(f);
        return 0;
      }
  
      if (f.inc_rest_mass) {
        y=(f.nu-f.m)/temper;
      } else {
        y=f.nu/temper;
      }

      // Number density
    
      f.n=this->fd_integ.calc_1o2(y);
      f.n*=f.g*pow(2.0*f.ms*temper,1.5)/4.0/o2scl_const::pi2;
    
      // Energy density:
    
      f.ed=this->fd_integ.calc_3o2(y);
      f.ed*=f.g*pow(2.0*f.ms*temper,2.5)/8.0/o2scl_const::pi2/f.ms;
    
      if (f.inc_rest_mass) {
      
        // Finish energy density
        f.ed+=f.n*f.m;
      
        // entropy density
        f.en=(5.0*(f.ed-f.n*f.m)/3.0-(f.nu-f.m)*f.n)/temper;
      
        // pressure
        f.pr=2.0*(f.ed-f.n*f.m)/3.0;
      
      } else {
      
        // entropy density
        f.en=(5.0*f.ed/3.0-f.nu*f.n)/temper;
      
        // pressure
        f.pr=2.0*f.ed/3.0;
      
      }

      if (!isfinite(f.nu) || !isfinite(f.n)) {
        O2SCL_ERR2("Chemical potential or density in ",
                   "fermion_nonrel::calc_mu().",exc_efailed);
      }
  
      return 0;
    }


    /** \brief Calculate properties as function of density

        If the density is zero, this function just sets part::mu,
        part::nu, part::ed, part::pr, and part::en to zero and returns
        without calling the error handler (even though at 
        zero density and finite	temperature, the chemical potentials
        formally are equal to \f$ -\infty \f$). 
    */
    virtual int calc_density(fermion_t &f, fp_t temper) {

      if (f.m<0.0 || (f.non_interacting==false && f.ms<0.0)) {
        O2SCL_ERR2("Mass negative in ",
                   "fermion_nonrel::calc_density().",exc_einval);
      }
      if (temper<0.0) {
        O2SCL_ERR2("Temperature less than zero in ",
                   "fermion_nonrel::calc_density().",exc_einval);
      }
      if (temper==0.0) {
        calc_density_zerot(f);
        return 0;
      }

      // Note it is important to throw if n=0 because the correct chemical
      // potential in that case is mu=-infty and we don't want encourage
      // the user to use this code in that case.
      if (f.n<=0.0) {
        std::string str=((std::string)"Density, ")+o2scl::dtos(f.n)+
          ", less than zero when temperature is positive in "+
          "fermion_nonrel::calc_density().";
        O2SCL_ERR(str.c_str(),exc_einval);
      }
  
      if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
      nu_from_n(f,temper);

      if (f.non_interacting) { f.mu=f.nu; }

      fp_t y, spi, ey, sy;
      if (f.inc_rest_mass) {
        y=-(f.nu-f.m)/temper;
      } else {
        y=-f.nu/temper;
      }

      // energy density
      f.ed=this->fd_integ.calc_3o2(-y);

      if (f.inc_rest_mass) {
    
        // Finish energy density
        f.ed*=f.g*pow(2.0*f.ms*temper,2.5)/8.0/o2scl_const::pi2/f.ms;
        f.ed+=f.n*f.m;
    
        // entropy density
        f.en=(5.0*(f.ed-f.n*f.m)/3.0-(f.nu-f.m)*f.n)/temper;
    
        // pressure
        f.pr=2.0*(f.ed-f.n*f.m)/3.0;
    
      } else {

        // Finish energy density
        f.ed*=f.g*pow(2.0*f.ms*temper,2.5)/8.0/o2scl_const::pi2/f.ms;

        // entropy density
        f.en=(5.0*f.ed/3.0-f.nu*f.n)/temper;
    
        // pressure
        f.pr=2.0*f.ed/3.0;
    
      }
  
      return 0;
    }

    /** \brief Compute thermodynamics with antiparticles at fixed
        chemical potential (unimplemented)
    */
    virtual int pair_mu(fermion_t &f, fp_t temper) {
      O2SCL_ERR2("Function fermion_nonrel::pair_mu() not ",
                 "implemented.",exc_eunimpl);
      return 0;
    }

    /** \brief Compute thermodynamics with antiparticles at fixed
        density (unimplemented)
    */
    virtual int pair_density(fermion_t &f, fp_t temper) {
      O2SCL_ERR2("Function fermion_nonrel::pair_density() not ",
                 "implemented.",exc_eunimpl);
      return 0;
    }
    
    /// Calculate effective chemical potential from density
    virtual void nu_from_n(fermion_t &f, fp_t temper) {

      fp_t init_n=f.n, init_m=f.m, init_ms=f.ms, init_nu=f.nu;
  
      // Use initial value of nu for initial guess
      fp_t nex;
      if (f.inc_rest_mass) {
        nex=-(f.nu-f.m)/temper;
      } else {
        nex=-f.nu/temper;
      } 

      // Choose a new initial guess if nex is too large
      // 
      // 0.4343 is approximately log10(exp(1)) so this value
      // is approximately -706 in double precision. The
      // value 0.9 is chosen as a buffer to the maximum value.
      if (nex>std::numeric_limits<fp_t>::max_exponent10/0.4343*0.9) {
        nex=std::numeric_limits<fp_t>::max_exponent10/0.4343/2.0;
      }
  
      func_t mf=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t)>
                          (&fermion_nonrel_tl<fermion_t,fd_inte_t,
                           be_inte_t,root_t,func_t,fp_t>::solve_fun),
                          this,std::placeholders::_1,f.n/f.g,f.ms*temper);
      
      // Turn off convergence errors temporarily, since we'll
      // try again if it fails
      bool enc=density_root->err_nonconv;
      density_root->err_nonconv=false;
      int ret=density_root->solve(nex,mf);

      // The root_cern class has a hard time when nex is near zero,
      // so we try a bracketing solver
      if (ret!=0) {
        fp_t blow=fabs(nex), bhigh=-blow;
        fp_t yhigh=mf(bhigh), ylow=mf(blow);
        for(size_t j=0;j<10 && yhigh*ylow>0.0;j++) {
          fp_t delta=fabs(blow);
          blow+=delta;
          bhigh-=delta;
          yhigh=mf(bhigh);
          ylow=mf(blow);
        }
        if ((yhigh<0.0 && ylow>0.0) || (yhigh>0.0 && ylow<0.0)) {
          o2scl::root_brent_gsl<func_t,fp_t> rbg;
          rbg.err_nonconv=false;
          ret=rbg.solve_bkt(blow,bhigh,mf);
          if (ret==0) nex=blow;
        }
      }
      
      if (ret!=0) {

        // If it failed, try to get a guess from classical_thermo particle
    
        classical_thermo_tl<fp_t> cl;
        cl.calc_density(f,temper);
        if (f.inc_rest_mass) {
          nex=-(f.nu-f.m)/temper;
        } else {
          nex=-f.nu/temper;
        } 
        ret=density_root->solve(nex,mf);
    
        // If it failed again, add error information
        if (ret!=0) {
          std::cout << "Function fermion_nonrel::nu_from_n() failed."
                    << std::endl;
          std::cout.precision(14);
          std::cout << "  n,m,ms,T,nu: " << init_n << " " << init_m << " "
                    << init_ms << " " << temper << " " << init_nu << std::endl;
          std::cout << "  ni,irm: " << f.non_interacting << " "
                    << f.inc_rest_mass << std::endl;
          O2SCL_ERR("Solver failed in fermion_nonrel::nu_from_n().",ret);
        }
      }

      // Restore the value of density_root->err_nonconv
      density_root->err_nonconv=enc;
  
      if (f.inc_rest_mass) {
        f.nu=-nex*temper+f.m;
      } else {
        f.nu=-nex*temper;
      }

      return;
    }


    /** \brief Set the solver for use in calculating the chemical
        potential from the density 
    */
    void set_density_root(root_t &rp) {
      density_root=&rp;
      return;
    }

    /// The default solver for calc_density().
    root_t def_density_root;
    
    /// If true, verify the thermodynamic identity
    bool verify_ti;
    
    /// Return string denoting type ("fermion_nonrel")
    virtual const char *type() { return "fermion_nonrel"; }

    /** \brief Function to compute chemical potential from density

        Variable \c nog is the target baryon density divided by
        the spin degeneracy, and \c msT is the effective mass
        times the temperature.
    */
    fp_t solve_fun(fp_t x, fp_t nog, fp_t msT) {

      fp_t nden;

      //
    
      // If the argument to calc_1o2() is less than about -708 (in
      // double precision) then an underflow occurs. The value 0.4343 is
      // approximately log10(exp(1)). We set nden to zero in this case,
      // as this helps the solver find the right root.

      if (-x<std::numeric_limits<fp_t>::min_exponent10/0.4343 ||
          !isfinite(x)) {
        nden=0.0;
      } else {
        nden=this->fd_integ.calc_1o2(-x);
      
      }
  
      nden*=pow(2.0*msT,1.5)/4.0/o2scl_const::pi2;
      fp_t ret=nden/nog-1.0;
      return ret;
    }

  protected:

    /// Solver to compute chemical potential from density
    root_t *density_root;
    
  private:

    fermion_nonrel_tl(const fermion_nonrel_tl &);
    fermion_nonrel_tl& operator=(const fermion_nonrel_tl&);

  };

  /** \brief Double precision version of \ref o2scl::fermion_nonrel_tl
  */
  typedef fermion_nonrel_tl<> fermion_nonrel;

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
  
  /** \brief Long double version of 
      \ref o2scl::fermion_nonrel_tl 
  */
  typedef fermion_nonrel_tl
  <fermion_tl<long double>,
   fermi_dirac_integ_direct<long double,funct_cdf25,
                            cpp_dec_float_25>,
   bessel_K_exp_integ_boost<long double,
                            cpp_dec_float_25>,
   root_brent_gsl<funct_ld,long double>,
   funct_ld,long double> fermion_nonrel_ld;
  
  /** \brief 25-digit version of 
      \ref o2scl::fermion_nonrel_tl 
  */
  typedef fermion_nonrel_tl
  <fermion_tl<cpp_dec_float_25>,
   fermi_dirac_integ_direct<cpp_dec_float_25,funct_cdf25,
                            cpp_dec_float_25>,
   bessel_K_exp_integ_boost<cpp_dec_float_25,
                            cpp_dec_float_25>,
   root_brent_gsl<funct_cdf25,cpp_dec_float_25>,
   funct_cdf25,cpp_dec_float_25> fermion_nonrel_cdf25;

#endif

}

#endif
