/* -------------------------------------------------------------------
  
  Copyright (C) 2006-2019, Andrew W. Steiner
  
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
#ifndef O2SCL_REL_FERMION_H
#define O2SCL_REL_FERMION_H

/** \file fermion_rel.h
    \brief File defining \ref o2scl::fermion_rel
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/fermion.h>
#include <o2scl/root_brent_gsl.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Equation of state for a relativistic fermion

      This class computes the thermodynamics of a relativistic fermion
      either as a function of the density or the chemical potential.
      It employs direct integration, using two different integrators
      for the degenerate and non-degenerate regimes. The default
      integrators are inte_qag_gsl (for degenerate fermions) and
      inte_qagiu_gsl (for non-degenerate fermions). For the functions
      calc_mu() and calc_density(), if the temperature argument is
      less than or equal to zero, the functions \ref
      fermion_zerot::calc_mu_zerot() and \ref
      fermion_zerot::calc_density_zerot() will be used to compute the
      result.

      \hline 
      <b>Degeneracy parameter:</b>

      Define the degeneracy parameter 
      \f[
      \psi=(\nu-m^{*})/T 
      \f] 
      where \f$ \nu \f$ is the effective chemical potential (including
      the rest mass) and \f$
      m^{*} \f$ is the effective mass. For \f$ \psi \f$ smaller than
      \ref min_psi, the non-degenerate expansion in \ref
      fermion_thermo::calc_mu_ndeg() is attempted first. If that
      fails, then integration is used. For \f$ \psi \f$ greater than
      \ref deg_limit (degenerate regime), a finite interval integrator
      is used and for \f$ \psi \f$ less than \ref deg_limit
      (non-degenerate regime), an integrator over the interval from
      \f$ [0,\infty) \f$ is used. In the case where \ref
      part::inc_rest_mass is false, the degeneracy parameter is
      \f[
      \psi=(\nu+m-m^{*})/T 
      \f] 

      <b>Integration limits:</b>

      The upper limit on the degenerate integration is given by
      \f[
      \mathrm{upper~limit} = \sqrt{{\cal L}^2-m^{*,2}}
      \f]
      where \f$ {\cal L}\equiv u T+\nu \f$ and \f$ u \f$ is \ref
      fermion_rel::upper_limit_fac . In the case where \ref
      part::inc_rest_mass is false, the result is
      \f[
      \mathrm{upper~limit} = \sqrt{(m+{\cal L})^2-m^{*2}}
      \f]
      
      The entropy is only significant at the Fermi surface, thus
      in the degenerate case, the lower limit of the entropy
      integral can be given be determined by the value of \f$ k \f$ 
      which solves
      \f[
      - u = \frac{\sqrt{k^2+m^{* 2}}-\nu}{T}
      \f]
      The solution is 
      \f[
      \mathrm{lower~limit} = \sqrt{(-u T+{\nu})^2-m^{*,2}}
      \f]
      but this solution is only valid if \f$ (m^{*}-\nu)/T < -u \f$.
      In the case where part::inc_rest_mass is false, the result is
      \f[
      \mathrm{lower~limit} = \sqrt{(-u T + m +\nu)^2-m^{*,2}}
      \f]
      which is valid if \f$ (m^{*}-\nu - m)/T < -u \f$.

      <b>Entropy integrand:</b>

      In the degenerate regime, the entropy integrand
      \f[
      - k^2 \left[ f \log f + \left(1-f\right) \log 
      \left(1-f \right) \right]
      \f]
      where \f$ f \f$ is the fermionic distribution function can lose
      precision when \f$ (E^{*} - \nu)/T \f$ is negative and
      sufficiently large in absolute magnitude. Thus when \f$ (E^{*} -
      \nu)/T < S \f$ where \f$ S \f$ is stored in \ref deg_entropy_fac
      (default is -30), the integrand is written as
      \f[
      -k^2 \left( E/T-\nu/T \right) e^{E/T-\nu/T} \, .
      \f]
      If \f$ (E - \nu)/T < S \f$ is less than -1 times \ref exp_limit
      (e.g. less than -200), then the entropy integrand is assumed 
      to be zero.
      
      <b>Non-degenerate integrands:</b>
      
      \comment
      It's not at all clear that this dimensionless form is more
      accurate than other potential alternatives. On the other hand,
      it seems that the uncertainties in the integrations are larger
      than the errors made by the integrand at present.
      \endcomment
      The integrands in the non-degenerate regime are written
      in a dimensionless form, by defining \f$ u \f$ with
      the relation
      \f$ p = \sqrt{\left(T u + m^{*}\right)^2-m^{* 2}} \f$,
      \f$ y \equiv \nu/ T \f$, and 
      \f$ \mathrm{mx} \equiv m^{*}/T \f$. 
      The density integrand is 
      \f[
      \left(\mathrm{mx}+u\right) \sqrt{u^2+2 (\mathrm{mx}) u}
      \left(\frac{e^{y}}{e^{\mathrm{mx}+u}+e^{y}}\right) \, , 
      \f]
      the energy integrand is 
      \f[
      \left(\mathrm{mx}+u\right)^2 \sqrt{u^2+2 (\mathrm{mx}) u}
      \left(\frac{e^{y}}{e^{\mathrm{mx}+u}+e^{y}}\right) \, ,
      \f]
      and the entropy integrand is 
      \f[
      \left(\mathrm{mx}+u\right) \sqrt{u^2+2 (\mathrm{mx}) u} 
      \left(t_1+t_2\right) \, ,
      \f]
      where 
      \f{eqnarray*}
      t_1 &=& \log \left(1+e^{y-\mathrm{mx}-u}\right)/
      \left(1+e^{y-\mathrm{mx}-u}\right) \nonumber \\
      t_2 &=& \log \left(1+e^{\mathrm{mx}+u-y}\right)/
      \left(1+e^{\mathrm{mx}+u-y}\right) \, .
      \f}

      \hline 
      <b>Accuracy:</b>

      The default settings for for this class give an accuracy of at
      least 1 part in \f$ 10^6 \f$ (and frequently better than this).

      When the integrators provide numerical uncertainties, these
      uncertainties are stored in \ref unc. In the case of
      calc_density() and pair_density(), the uncertainty from the
      numerical accuracy of the solver is not included. (There is also
      a relatively small inaccuracy due to the mathematical evaluation
      of the integrands which is not included in \ref unc.)
     
      One can improve the accuracy to within 1 part in \f$ 10^{10} \f$ 
      using
      \code
      fermion_rel rf(1.0,2.0);
      rf.upper_limit_fac=40.0;
      rf.dit->tol_abs=1.0e-13;
      rf.dit->tol_rel=1.0e-13;
      rf.nit->tol_abs=1.0e-13;
      rf.nit->tol_rel=1.0e-13;
      rf.density_root->tol_rel=1.0e-10;
      \endcode
      which decreases the both the relative and absolute tolerances
      for both the degenerate and non-degenerate integrators and
      improves the accuracy of the solver which determines the
      chemical potential from the density. Of course if these
      tolerances are too small, the calculation may fail.

      \hline 
      <b>Todos:</b>

      \future The expressions which appear in in the integrand
      functions density_fun(), etc. could likely be improved,
      especially in the case where \ref o2scl::part::inc_rest_mass is
      <tt>false</tt>. There should not be a need to check if
      <tt>ret</tt> is finite.

      \future It appears this class doesn't compute the uncertainty in
      the chemical potential or density with calc_density(). This
      could be fixed.

      \future I'd like to change the lower limit on the entropy 
      integration, but the value in the code at the moment (stored
      in <tt>ll</tt>) makes bm_part2.cpp worse.

      \future The function pair_mu() should set the antiparticle
      integrators as done in fermion_deriv_rel.
  */
  class fermion_rel : public fermion_thermo {

  public:

    /// \name Numerical parameters
    //@{
    /** \brief If true, call the error handler when convergence 
	fails (default true)
    */
    bool err_nonconv;

    /** \brief The smallest value of \f$ (\mu-m)/T \f$ for which 
	integration is used
    */
    double min_psi;

    /** \brief The critical degeneracy at which to switch integration 
	techniques (default 2)
    */
    double deg_limit;
    
    /** \brief The limit for exponentials to ensure integrals are finite 
	(default 200)
    */
    double exp_limit;

    /// The factor for the degenerate upper limits (default 20)
    double upper_limit_fac;

    /// A factor for the degenerate entropy integration (default 30)
    double deg_entropy_fac;

    /// Verbosity parameter (default 0)
    int verbose;
    //@}

    /// Storage for the uncertainty
    fermion unc;

    /// If true, use expansions for extreme conditions (default true)
    bool use_expansions;

    /// Create a fermion with mass \c m and degeneracy \c g
    fermion_rel();

    virtual ~fermion_rel();

    /** \brief Calculate properties as function of chemical potential
     */
    virtual void calc_mu(fermion &f, double temper);

    /** \brief Calculate properties as function of density

        This function uses the current value of \c nu (or \c mu if the
	particle is non interacting) for an initial guess to solve for
	the chemical potential. If this guess is too small, then this
	function may fail.
    */
    virtual int calc_density(fermion &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	chemical potential
    */
    virtual void pair_mu(fermion &f, double temper);

    /** \brief Calculate properties with antiparticles as function of
	density
    */
    virtual int pair_density(fermion &f, double temper);

    /** \brief Calculate effective chemical potential from density
     */
    virtual int nu_from_n(fermion &f, double temper);
    
    /// The non-degenerate integrator
    std::shared_ptr<inte<> > nit;

    /// The degenerate integrator
    std::shared_ptr<inte<> > dit;

    /// The solver for calc_density()
    std::shared_ptr<root<> > density_root;

    /// Return string denoting type ("fermion_rel")
    virtual const char *type() { return "fermion_rel"; }

    /** \brief An integer indicating the last numerical method used

	In all functions
	- 0: no previous calculation or last calculation failed

        In \ref nu_from_n_tlate():
	- 1: default solver
	- 2: default solver with smaller tolerances
	- 3: bracketing solver

	In \ref calc_mu_tlate():
	- 4: non-degenerate expansion
	- 5: degenerate expansion
	- 6: exact integration, non-degenerate integrands
	- 7: exact integration, degenerate integrands, lower limit
	 on entropy integration
	- 8: exact integration, degenerate integrands, full
	entropy integration
	- 9: T=0 result

	In \ref calc_density_tlate(), the integer is a two-digit
	number. The first digit (1 to 3) is the method used by \ref
	nu_from_n_tlate() and the second digit is one of
	- 1: nondegenerate expansion
	- 2: degenerate expansion
	- 3: exact integration, non-degenerate integrands
	- 4: exact integration, degenerate integrands, lower limit
	 on entropy integration
	- 5: exact integration, degenerate integrands, full
	entropy integration
	If \ref calc_density_tlate() uses the T=0 code, then
	last_method is 40. 

	In \ref pair_mu_tlate(), the integer is a three-digit number.
	The third digit is always 0 (to ensure a value of last_method
	which is unique from the other values reported from other
	functions as described above). The first digit is the method
	used for particles from \ref calc_mu_tlate() above and the
	second digit is the method used for antiparticles. 

	In \ref pair_density_tlate(), the integer is a four-digit
	number. The first digit is from the list below and the
	remaining three digits, if nonzero, are from \ref
	pair_mu_tlate().
	- 1: T=0 result
	- 2: default solver
	- 3: bracketing solver
	- 4: default solver with smaller tolerances
	- 5: default solver with smaller tolerances in log units
	- 6: bracketing in log units
    */
    int last_method;
    
    /// \name Template versions of base functions
    //@{
    /** \brief Calculate the chemical potential from the density
	(template version)
     */
    template<class fermion_t>
      int nu_from_n_tlate(fermion_t &f, double temper) {

      last_method=0;
      
      double nex;

      // Try to ensure a good initial guess

      nex=f.nu/temper;
      double y=solve_fun(nex,f,temper);
      if (verbose>1) {
	std::cout << "nu_from_n(): initial guess " << nex << std::endl;
      }

      if (y>1.0-1.0e-6) {
	double scale=f.ms;
	if (temper>scale) scale=temper;
	for(size_t i=0;i<10;i++) {
	  if (nex<0.0) nex+=scale*1.0e5;
	  else nex*=10.0;
	  y=solve_fun(nex,f,temper);
	  if (y<1.0-1.0e-6) i=10;
	}
	if (verbose>1) {
	  std::cout << "nu_from_n(): adjusted guess to " << nex << std::endl;
	}
      }

      // If that didn't work, try a different guess
      if (y>1.0-1.0e-6) {
	if (f.inc_rest_mass) {
	  nex=f.ms/temper;
	} else {
	  nex=(f.ms-f.m)/temper;
	}
	y=solve_fun(nex,f,temper);
	if (verbose>1) {
	  std::cout << "nu_from_n(): adjusted guess (try 2) to "
		    << nex << std::endl;
	}
      }
  
      // If neither worked, call the error handler
      if (y==1.0 || !std::isfinite(y)) {
	O2SCL_CONV2_RET("Couldn't find reasonable initial guess in ",
			"fermion_rel::nu_from_n().",exc_einval,
			this->err_nonconv);
      }

      // Perform full solution
      funct mf=std::bind(std::mem_fn<double(double,fermion &,double)>
			 (&fermion_rel::solve_fun),
			 this,std::placeholders::_1,std::ref(f),temper);

      // The default o2scl::root object is of type root_cern,
      // and this solver has problems when the root is near 0.
      // Below, we switch to a root_brent_gsl object in the case
      // that the default solver fails.
  
      bool drec=density_root->err_nonconv;
      density_root->err_nonconv=false;
      int ret=density_root->solve(nex,mf);

      if (ret!=0) {
    
	if (verbose>1) {
	  std::cout << "nu_from_n(): density_root failed x="
		    << nex << " ." << std::endl;
	  std::cout << "\tTrying to make integrators more accurate."
		    << std::endl;
	}

	// If it fails, try to make the integrators more accurate
	double tol1=dit->tol_rel, tol2=dit->tol_abs;
	double tol3=nit->tol_rel, tol4=nit->tol_abs;
	dit->tol_rel/=1.0e2;
	dit->tol_abs/=1.0e2;
	nit->tol_rel/=1.0e2;
	nit->tol_abs/=1.0e2;
	ret=density_root->solve(nex,mf);

	if (ret!=0) {

	  if (verbose>1) {
	    std::cout << "nu_from_n(): density_root failed again x=" << nex
		      << " ." << std::endl;
	    std::cout << "Trying to bracket root." << std::endl;
	  }
      
	  double lg=std::max(fabs(f.nu),f.ms);
	  double bhigh=lg/temper, blow=-bhigh;
	  double yhigh=mf(bhigh), ylow=mf(blow);
	  for(size_t j=0;j<5 && yhigh>0.0;j++) {
	    bhigh*=1.0e2;
	    yhigh=mf(bhigh);
	  }
	  for(size_t j=0;j<5 && ylow<0.0;j++) {
	    blow*=1.0e2;
	    ylow=mf(blow);
	  }
	  if (yhigh<0.0 && ylow>0.0) {
	    o2scl::root_brent_gsl<> rbg;
	    rbg.err_nonconv=false;
	    ret=rbg.solve_bkt(blow,bhigh,mf);
	    if (ret==0) {
	      // Bracketing solver worked
	      last_method=3;
	      nex=blow;
	    } else {
	      if (verbose>1) {
		std::cout << "nu_from_n(): density_root failed fourth solver "
			  << blow << std::endl;
	      }
	    }
	  } else if (verbose>1) {
	    std::cout << "nu_from_n(): Failed to bracket." << std::endl;
	  }
	} else {
	  // Increasing tolerances worked
	  last_method=2;
	}

	// Return tolerances to their original values
	dit->tol_rel=tol1;
	dit->tol_abs=tol2;
	nit->tol_rel=tol3;
	nit->tol_abs=tol4;

      } else {
	// First solver worked
	last_method=1;
      }

      density_root->err_nonconv=drec;

      if (ret!=0) {
	O2SCL_CONV2_RET("Density solver failed in ",
			"fermion_rel::nu_from_n().",exc_efailed,
			this->err_nonconv);
      }

      f.nu=nex*temper;

      return success;
    }

    /** \brief Calculate properties as function of chemical potential
	(template version)
     */
    template<class fermion_t>
      void calc_mu_tlate(fermion_t &f, double temper) {

      last_method=0;
      
      // -----------------------------------------------------------------
      // Handle T<=0

      if (temper<0.0) {
	O2SCL_ERR2("Temperature less than zero in ",
		   "fermion_rel::calc_mu_tlate().",exc_einval);
      }
      if (temper==0.0) {
	calc_mu_zerot(f);
	last_method=9;
	return;
      }

      if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

      // Compute the degeneracy parameter
  
      bool deg=true;
      double psi;
      if (f.inc_rest_mass) {
	psi=(f.nu-f.ms)/temper;
      } else {
	psi=(f.nu+(f.m-f.ms))/temper;
      }
      if (psi<deg_limit) deg=false;
  
      // Try the non-degenerate expansion if psi is small enough
      if (use_expansions && psi<min_psi) {
	bool acc=calc_mu_ndeg(f,temper,1.0e-14);
	if (acc) {
	  unc.n=f.n*1.0e-14;
	  unc.ed=f.ed*1.0e-14;
	  unc.pr=f.pr*1.0e-14;
	  unc.en=f.en*1.0e-14;
	  last_method=4;
	  return;
	}
      }

      // Try the degenerate expansion if psi is large enough
      if (use_expansions && psi>20.0) {
	bool acc=calc_mu_deg(f,temper,1.0e-14);
	if (acc) {
	  unc.n=f.n*1.0e-14;
	  unc.ed=f.ed*1.0e-14;
	  unc.pr=f.pr*1.0e-14;
	  unc.en=f.en*1.0e-14;
	  last_method=5;
	  return;
	}
      }

      if (!deg) {
    
	// If the temperature is large enough, perform the full integral
    
	funct mfd=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::density_fun),
			    this,std::placeholders::_1,std::ref(f),temper);
	funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::energy_fun),
			    this,std::placeholders::_1,std::ref(f),temper);
	funct mfs=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::entropy_fun),
			    this,std::placeholders::_1,std::ref(f),temper);
      
	double prefac=f.g*pow(temper,3.0)/2.0/o2scl_const::pi2;

	// Compute the number density
    
	f.n=nit->integ(mfd,0.0,0.0);
	f.n*=prefac;
	unc.n=nit->get_error()*prefac;

	// Compute the energy density

	f.ed=nit->integ(mfe,0.0,0.0);
	f.ed*=prefac*temper;
	if (!f.inc_rest_mass) f.ed-=f.n*f.m;
	unc.ed=nit->get_error()*prefac*temper;
    
	// Compute the entropy

	f.en=nit->integ(mfs,0.0,0.0);
	f.en*=prefac;
	unc.en=nit->get_error()*prefac;

	last_method=6;

      } else {
    
	// Otherwise, apply a degenerate approximation, by making the
	// upper integration limit finite
    
	funct mfd=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::deg_density_fun),
			    this,std::placeholders::_1,std::ref(f),temper);
	funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::deg_energy_fun),
			    this,std::placeholders::_1,std::ref(f),temper);
	funct mfs=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::deg_entropy_fun),
			    this,std::placeholders::_1,std::ref(f),temper);

	double prefac=f.g/2.0/o2scl_const::pi2;
    
	// Compute the upper limit for degenerate integrals

	double arg;
	if (f.inc_rest_mass) {
	  arg=pow(upper_limit_fac*temper+f.nu,2.0)-f.ms*f.ms;
	} else {
	  arg=pow(upper_limit_fac*temper+f.nu+f.m,2.0)-f.ms*f.ms;
	}
	double ul;
	if (arg>0.0) {
	  ul=sqrt(arg);
	} else {
	  f.n=0.0;
	  f.ed=0.0;
	  f.pr=0.0;
	  f.en=0.0;
	  unc.n=0.0;
	  unc.ed=0.0;
	  unc.pr=0.0;
	  unc.en=0.0;
	  O2SCL_ERR2("Zero density in degenerate limit in fermion_rel::",
		     "calc_mu(). Variable deg_limit set improperly?",
		     exc_efailed);
	  return;
	}
    
	// Compute the number density

	f.n=dit->integ(mfd,0.0,ul);
	f.n*=prefac;
	unc.n=dit->get_error()*prefac;
    
	// Compute the energy density

	f.ed=dit->integ(mfe,0.0,ul);
	f.ed*=prefac;
	unc.ed=dit->get_error()*prefac;

	// Compute the lower limit for the entropy integration

	double ll;
	if (f.inc_rest_mass) {
	  arg=pow(-upper_limit_fac*temper+f.nu,2.0)-f.ms*f.ms;
	  if (arg>0.0 && (f.ms-f.nu)/temper<-upper_limit_fac) {
	    ll=sqrt(arg);
	  } else {
	    ll=-1.0;
	  }
	} else {
	  arg=pow(-upper_limit_fac*temper+f.nu+f.m,2.0)-f.ms*f.ms;
	  if (arg>0.0 && (f.ms-f.nu-f.m)/temper<-upper_limit_fac) {
	    ll=sqrt(arg);
	  } else {
	    ll=-1.0;
	  }
	}

	// Compute the entropy

	if (ll>0.0) {
	  f.en=dit->integ(mfs,ll,ul);
	  last_method=7;
	} else {
	  f.en=dit->integ(mfs,0.0,ul);
	  last_method=8;
	}
	f.en*=prefac;
	unc.en=dit->get_error()*prefac;
    
      }

      // Compute the pressure using the thermodynamic identity

      f.pr=-f.ed+temper*f.en+f.nu*f.n;
      unc.pr=sqrt(unc.ed*unc.ed+temper*unc.en*temper*unc.en+
		  f.nu*unc.n*f.nu*unc.n);

      return;
    }

    /** \brief Calculate properties as function of density
	(template version)

	\future There is still quite a bit of code duplication
	between this function and \ref calc_mu_tlate() .
     */
    template<class fermion_t>
      int calc_density_tlate(fermion_t &f, double temper) {

      last_method=0;
      
      // The code below may modify the density, which is confusing to
      // the user, so we store it here so we can restore it before
      // this function returns
      double density_temp=f.n;
  
      // -----------------------------------------------------------------
      // Handle T<=0

      if (temper<0.0) {
	O2SCL_ERR2("Temperature less than zero in ",
		   "fermion_rel::calc_density_tlate().",
		   exc_einval);
      }
      if (temper==0.0) {
	last_method=40;
	calc_density_zerot(f);
	return 0;
      }

#if !O2SCL_NO_RANGE_CHECK
      // This may not be strictly necessary, because it should be clear
      // that this function will produce gibberish if the density is not
      // finite, but I've found this extra checking of the inputs useful
      // for debugging.
      if (!std::isfinite(f.n)) {
	O2SCL_ERR2("Density not finite in ",
		   "fermion_rel::calc_density_tlate().",exc_einval);
      }
#endif

      // -----------------------------------------------------------------
      // First determine the chemical potential by solving for the density

      if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }
  
      int ret=nu_from_n(f,temper);
      last_method*=10;
      if (ret!=0) {
	O2SCL_CONV2_RET("Function calc_density() failed in fermion_rel::",
			"calc_density().",exc_efailed,this->err_nonconv);
      }

      if (f.non_interacting) { f.mu=f.nu; }

      // -----------------------------------------------------------------
      // Now use the chemical potential to compute the energy density,
      // pressure, and entropy

      bool deg=true;
      double psi;
      if (f.inc_rest_mass) {
	psi=(f.nu-f.ms)/temper;
      } else {
	psi=(f.nu+(f.m-f.ms))/temper;
      }
      if (psi<deg_limit) deg=false;

      // Try the non-degenerate expansion if psi is small enough
      if (use_expansions && psi<min_psi) {
	bool acc=calc_mu_ndeg(f,temper,1.0e-14);
	if (acc) {
	  unc.ed=f.ed*1.0e-14;
	  unc.pr=f.pr*1.0e-14;
	  unc.en=f.en*1.0e-14;
	  f.n=density_temp;
	  last_method+=1;
	  return 0;
	}
      }
  
      // Try the degenerate expansion if psi is large enough
      if (use_expansions && psi>20.0) {
	bool acc=calc_mu_deg(f,temper,1.0e-14);
	if (acc) {
	  unc.n=f.n*1.0e-14;
	  unc.ed=f.ed*1.0e-14;
	  unc.pr=f.pr*1.0e-14;
	  unc.en=f.en*1.0e-14;
	  f.n=density_temp;
	  last_method+=2;
	  return 0;
	}
      }

      if (!deg) {
    
	funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::energy_fun),
			    this,std::placeholders::_1,std::ref(f),temper);
	funct mfs=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::entropy_fun),
			    this,std::placeholders::_1,std::ref(f),temper);
    
	f.ed=nit->integ(mfe,0.0,0.0);
	f.ed*=f.g*pow(temper,4.0)/2.0/o2scl_const::pi2;
	if (!f.inc_rest_mass) f.ed-=f.n*f.m;
	unc.ed=nit->get_error()*f.g*pow(temper,4.0)/2.0/o2scl_const::pi2;
    
	f.en=nit->integ(mfs,0.0,0.0);
	f.en*=f.g*pow(temper,3.0)/2.0/o2scl_const::pi2;
	unc.en=nit->get_error()*f.g*pow(temper,3.0)/2.0/o2scl_const::pi2;
	last_method+=3;

      } else {

	funct mfe=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::deg_energy_fun),
			    this,std::placeholders::_1,std::ref(f),temper);
	funct mfs=std::bind(std::mem_fn<double(double,fermion &,double)>
			    (&fermion_rel::deg_entropy_fun),
			    this,std::placeholders::_1,std::ref(f),temper);
      
	double arg;
	if (f.inc_rest_mass) {
	  arg=pow(upper_limit_fac*temper+f.nu,2.0)-f.ms*f.ms;
	} else {
	  arg=pow(upper_limit_fac*temper+f.nu+f.m,2.0)-f.ms*f.ms;
	}
	double ul;
	if (arg>0.0) {
      
	  ul=sqrt(arg);
      
	  double ll;
	  if (f.inc_rest_mass) {
	    arg=pow(-upper_limit_fac*temper+f.nu,2.0)-f.ms*f.ms;
	    if (arg>0.0 && (f.ms-f.nu)/temper<-upper_limit_fac) {
	      ll=sqrt(arg);
	    } else {
	      ll=-1.0;
	    }
	  } else {
	    arg=pow(-upper_limit_fac*temper+f.nu+f.m,2.0)-f.ms*f.ms;
	    if (arg>0.0 && (f.ms-f.nu-f.m)/temper<-upper_limit_fac) {
	      ll=sqrt(arg);
	    } else {
	      ll=-1.0;
	    }
	  }
      
	  f.ed=dit->integ(mfe,0.0,ul);
	  f.ed*=f.g/2.0/o2scl_const::pi2;
	  unc.ed=dit->get_error()*f.g/2.0/o2scl_const::pi2;
      
	  if (ll>0.0) {
	    f.en=dit->integ(mfs,ll,ul);
	    last_method+=4;
	  } else {
	    f.en=dit->integ(mfs,0.0,ul);
	    last_method+=5;
	  }
	  f.en*=f.g/2.0/o2scl_const::pi2;
	  unc.en=dit->get_error()*f.g/2.0/o2scl_const::pi2;
      
	} else {

	  f.ed=0.0;
	  f.en=0.0;
	  unc.ed=0.0;
	  unc.en=0.0;
	  O2SCL_ERR2("Zero density in degenerate limit in fermion_rel::",
		     "calc_mu(). Variable deg_limit set improperly?",
		     exc_efailed);
      
	}
      }

      f.n=density_temp;
      f.pr=-f.ed+temper*f.en+f.nu*f.n;
      unc.pr=sqrt(unc.ed*unc.ed+temper*unc.en*temper*unc.en+
		  f.nu*unc.n*f.nu*unc.n);
  
      return 0;
    }

    /** \brief Calculate properties with antiparticles as function of
	chemical potential (template version)
     */
    template<class fermion_t>
      void pair_mu_tlate(fermion_t &f, double temper) {

      last_method=0;
      
      if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }
      
      // AWS: 2/12/19: I'm taking this out, similar to the removal
      // of the code in fermion_rel::pair_fun(). If I put it
      // back in, I need to find a new value for last_method
      // other than 9.
      if (false && use_expansions) {
	if (calc_mu_ndeg(f,temper,1.0e-14,true)) {
	  unc.n=1.0e-14*f.n;
	  unc.ed=1.0e-14*f.ed;
	  unc.en=1.0e-14*f.en;
	  unc.pr=1.0e-14*f.pr;
	  last_method=9;
	  return;
	}
      }

      fermion antip(f.m,f.g);
      f.anti(antip);

      // Particles
      calc_mu(f,temper);
      double unc_n=unc.n;
      double unc_pr=unc.pr;
      double unc_ed=unc.ed;
      double unc_en=unc.en;

      // Antiparticles
      int lm=last_method*10;
      calc_mu(antip,temper);
      last_method+=lm;
      last_method*=10;

      // Add up thermodynamic quantities
      if (f.inc_rest_mass) {
	f.ed+=antip.ed;
      } else {
	f.ed=f.ed+antip.ed+2.0*antip.n*f.m;
      }
      f.n-=antip.n;
      f.pr+=antip.pr;
      f.en+=antip.en;

      // Add up uncertainties
      unc.n=gsl_hypot(unc.n,unc_n);
      unc.ed=gsl_hypot(unc.ed,unc_ed);
      unc.pr=gsl_hypot(unc.pr,unc_pr);
      unc.en=gsl_hypot(unc.ed,unc_en);

      return;
    }

    /** \brief Calculate thermodynamic properties with antiparticles
	from the density (template version)
     */
    template<class fermion_t>
      int pair_density_tlate(fermion_t &f, double temper) {

      last_method=0;
      
      // -----------------------------------------------------------------
      // Handle T<=0

      if (temper<=0.0) {
	calc_density_zerot(f);
	last_method=1000;
	return success;
      }

      // Storage the input density
      double density_temp=f.n;
  
      if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

      double initial_guess=f.nu;
  
      double nex=f.nu/temper;
      
      funct mf=std::bind(std::mem_fn<double(double,fermion &,double,bool)>
			 (&fermion_rel::pair_fun),
			 this,std::placeholders::_1,std::ref(f),temper,false);

      // Begin by trying the user-specified guess
      bool drec=density_root->err_nonconv;
      density_root->err_nonconv=false;
      int ret=density_root->solve(nex,mf);

      // If that doesn't work, try bracketing the root
      if (ret==0) {
	last_method=2000;
      } else {
	double lg=std::max(fabs(f.nu),f.ms);
	double bhigh=lg/temper, blow=-bhigh;
	double yhigh=mf(bhigh), ylow=mf(blow);
	for(size_t j=0;j<5 && yhigh<0.0;j++) {
	  bhigh*=1.0e2;
	  yhigh=mf(bhigh);
	}
	for(size_t j=0;j<5 && ylow>0.0;j++) {
	  blow*=1.0e2;
	  ylow=mf(blow);
	}
	if (yhigh>0.0 && ylow<0.0) {
	  root_brent_gsl<> rbg;
	  rbg.err_nonconv=false;
	  ret=rbg.solve_bkt(blow,bhigh,mf);
	  if (ret==0) {
	    nex=blow;
	    last_method=3000;
	  }
	}
      }
  
      if (ret!=0) {

	// If that fails, try to make the integrators more accurate
	double tol1=dit->tol_rel, tol2=dit->tol_abs;
	double tol3=nit->tol_rel, tol4=nit->tol_abs;
	dit->tol_rel/=1.0e2;
	dit->tol_abs/=1.0e2;
	nit->tol_rel/=1.0e2;
	nit->tol_abs/=1.0e2;
	ret=density_root->solve(nex,mf);
	if (ret==0) last_method=4000;
    
	// AWS: 7/25/18: We work in log units below, so we ensure the
	// chemical potential is not negative
	if (nex<0.0) nex=1.0e-10;
    
	// If that failed, try working in log units

	// Function in log units
	funct lmf=std::bind(std::mem_fn<double(double,fermion &,
					       double,bool)>
			    (&fermion_rel::pair_fun),
			    this,std::placeholders::_1,std::ref(f),
			    temper,true);
    
	if (ret!=0) {
	  nex=log(nex);
	  ret=density_root->solve(nex,lmf);
	  nex=exp(nex);
	  if (ret==0) last_method=5000;
	}
	
	if (ret!=0) {
	  // If that failed, try a different solver
	  root_brent_gsl<> rbg;
	  rbg.err_nonconv=false;
	  nex=log(nex);
	  ret=rbg.solve(nex,lmf);
	  nex=exp(nex);
	  if (ret==0) last_method=6000;
	}

	// Return tolerances to their original values
	dit->tol_rel=tol1;
	dit->tol_abs=tol2;
	nit->tol_rel=tol3;
	nit->tol_abs=tol4;
      }

      // Restore value of err_nonconv
      density_root->err_nonconv=drec;

      if (ret!=0) {
	std::cout.precision(14);
	std::cout << "m,ms,n,T: " << f.m << " " << f.ms << " "
	     << f.n << " " << temper << std::endl;
	std::cout << "nu: " << initial_guess << std::endl;
	O2SCL_CONV2_RET("Density solver failed in fermion_rel::",
			"pair_density().",exc_efailed,this->err_nonconv);
      }

      f.nu=nex*temper;
  
      if (f.non_interacting==true) { f.mu=f.nu; }

      // Finally, now that we have the chemical potential, use pair_mu()
      // to evaluate the energy density, pressure, and entropy
      int lm=last_method;
      pair_mu(f,temper);
      last_method+=lm;

      // The function pair_mu() can modify the density, which would be
      // confusing to the user, so we return it to the user-specified
      // value.
      f.n=density_temp;

      return success;
    }
    //@}

  protected:
    
#ifndef DOXYGEN_INTERNAL

    /// The integrand for the density for non-degenerate fermions
    double density_fun(double u, fermion &f, double T);

    /// The integrand for the energy density for non-degenerate fermions
    double energy_fun(double u, fermion &f, double T);

    /// The integrand for the entropy density for non-degenerate fermions
    double entropy_fun(double u, fermion &f, double T);

    /// The integrand for the density for degenerate fermions
    double deg_density_fun(double u, fermion &f, double T);

    /// The integrand for the energy density for degenerate fermions
    double deg_energy_fun(double u, fermion &f, double T);

    /// The integrand for the entropy density for degenerate fermions
    double deg_entropy_fun(double u, fermion &f, double T);

    /// Solve for the chemical potential given the density
    double solve_fun(double x, fermion &f, double T);

    /** \brief Solve for the chemical potential given the density 
	with antiparticles
	
	\future Particles and antiparticles have different degeneracy
	factors, so we separately use the expansions one at a time. It
	is probably better to separately generate a new expansion
	function which automatically handles the sum of particles and
	antiparticles.
    */
    double pair_fun(double x, fermion &f, double T, bool log_mode);
    
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
