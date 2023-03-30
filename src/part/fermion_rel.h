/* -------------------------------------------------------------------
  
   Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifndef O2SCL_FERMION_REL_H
#define O2SCL_FERMION_REL_H

/** \file fermion_rel.h
    \brief File defining \ref o2scl::fermion_rel_tl
*/

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#ifdef O2SCL_LD_TYPES
#include <boost/multiprecision/cpp_dec_float.hpp>
#endif

#include <o2scl/constants.h>
#include <o2scl/mroot.h>
#include <o2scl/inte.h>
#include <o2scl/fermion.h>
#include <o2scl/root_brent_gsl.h>
#include <o2scl/inte_qagiu_gsl.h>
#include <o2scl/inte_qag_gsl.h>
#include <o2scl/polylog.h>
#include <o2scl/inte_kronrod_boost.h>
#include <o2scl/inte_adapt_cern.h>

namespace o2scl {

  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<25> > cpp_dec_float_25;
  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<35> > cpp_dec_float_35;
  
  /** \brief Integrands for \ref o2scl::fermion_rel_tl
   */
  class fermion_rel_integ_base {
    
  protected:
    
    /** \brief The limit for exponentials to ensure integrals are finite 
	(default 200.0)
    */
    double exp_limit;

    /** \brief A factor for the degenerate entropy integration
        (default 30.0)
    */
    double deg_entropy_fac;

    /** \brief If true, call the error handler if the integration
        does not succeed
    */
    bool err_nonconv;
    
  public:

    fermion_rel_integ_base() {
      exp_limit=200.0;
      deg_entropy_fac=30.0;
      err_nonconv=true;
    }      
    
    /// The integrand for the density for non-degenerate fermions
    template<class internal_fp_t>
    internal_fp_t density_fun(internal_fp_t u, internal_fp_t y,
                              internal_fp_t eta) {
      
      internal_fp_t ret;
        
      internal_fp_t arg1=u*u+2*eta*u;
      internal_fp_t arg2=eta+u-y;
      internal_fp_t arg3=eta+u;
      if (y-eta-u>exp_limit) {
	ret=arg3*sqrt(arg1);
      } else if (y>u+exp_limit && eta>u+exp_limit) {
	ret=arg3*sqrt(arg1)/(exp(arg2)+1);
      } else {
	ret=arg3*sqrt(arg1)*exp(y)/(exp(arg3)+exp(y));
      }
      
      if (!isfinite(ret)) {
	ret=0.0;
      }

      return ret;
    }

    /// The integrand for the pressure for non-degenerate fermions
    template<class internal_fp_t>
    internal_fp_t pressure_fun_old(internal_fp_t u, internal_fp_t y,
                               internal_fp_t eta) {

      internal_fp_t ret;
      
      internal_fp_t arg1=u*u+2*eta*u;
      internal_fp_t term1=sqrt(arg1);
      internal_fp_t arg3=eta+u;
      ret=term1*term1*term1*exp(y)/(exp(arg3)+exp(y))/3;
      
      if (!isfinite(ret)) {
	ret=0.0;
      }

      return ret;
    }

    /// The integrand for the pressure for non-degenerate fermions
    template<class internal_fp_t>
    internal_fp_t pressure_fun(internal_fp_t u, internal_fp_t y,
                               internal_fp_t eta) {

      internal_fp_t ret;
      
      internal_fp_t arg1=u*u+2*eta*u;
      internal_fp_t term1=sqrt(arg1);
      internal_fp_t arg3=eta+u;
      ret=term1*arg3*log1p(exp(y-arg3));
      
      if (!isfinite(ret)) {
	ret=0.0;
      }

      return ret;
    }

    /// The integrand for the energy density for non-degenerate fermions
    template<class internal_fp_t>
    internal_fp_t energy_fun(internal_fp_t u, internal_fp_t y,
                             internal_fp_t eta) {

      internal_fp_t ret;

      internal_fp_t arg1=u*u+2*eta*u;
      internal_fp_t arg2=eta+u-y;
      internal_fp_t arg3=eta+u;
      if (y>u+exp_limit && eta>u+exp_limit) {
	ret=arg3*arg3*sqrt(arg1)/(exp(arg2)+1);
      } else {
	ret=arg3*arg3*sqrt(arg1)*exp(y)/(exp(arg3)+exp(y));
      }
 
      if (!isfinite(ret)) {
	return 0;
      }

      return ret;
    }

    /// The integrand for the entropy density for non-degenerate fermions
    template<class internal_fp_t>
    internal_fp_t entropy_fun(internal_fp_t u, internal_fp_t y,
                              internal_fp_t eta) {

      internal_fp_t ret;

      internal_fp_t arg1=u*u+2*eta*u;
      internal_fp_t arg2=eta+u-y;
      internal_fp_t arg3=eta+u;
      internal_fp_t arg4=y-eta-u;
      internal_fp_t arg5=1+exp(arg4);
      internal_fp_t arg6=1+exp(arg2);
      internal_fp_t term1=log(arg5)/arg5;
      internal_fp_t term2=log(arg6)/arg6;
      ret=arg3*sqrt(arg1)*(term1+term2);
  
      if (!isfinite(ret)) {
	return 0.0;
      }

      return ret;
    }

    /// The integrand for the density for degenerate fermions
    template<class internal_fp_t>
    internal_fp_t deg_density_fun(internal_fp_t k, internal_fp_t T,
                                  internal_fp_t y, internal_fp_t eta,
                                  internal_fp_t mot, bool debug) {

      internal_fp_t ret;
      
      internal_fp_t E=hypot(k/T,eta)-mot;
      internal_fp_t arg1=E-y;
      ret=k*k/(1+exp(arg1));

      if (debug) {
        std::cout << k << " " << ret << std::endl;
      }
          
      if (!isfinite(ret)) {
	O2SCL_ERR2("Returned not finite result ",
		   "in fermion_rel::deg_density_fun().",exc_einval);
      }
      return ret;
    }

    /// The integrand for the energy density for degenerate fermions
    template<class internal_fp_t>
    internal_fp_t deg_energy_fun(internal_fp_t k, internal_fp_t T,
                                 internal_fp_t y, internal_fp_t eta,
                                 internal_fp_t mot) {

      internal_fp_t ret;
      internal_fp_t E=hypot(k/T,eta)-mot;
      internal_fp_t arg1=E-y;
      
      ret=k*k*E*T/(1+exp(arg1));
      
      if (!isfinite(ret)) {
	O2SCL_ERR2("Returned not finite result ",
		   "in fermion_rel::deg_energy_fun().",exc_einval);
      }
  
      return ret;
    }

    /// The integrand for the energy density for degenerate fermions
    template<class internal_fp_t>
    internal_fp_t deg_pressure_fun(internal_fp_t k, internal_fp_t T,
                                   internal_fp_t y, internal_fp_t eta,
                                   internal_fp_t mot, bool debug) {

      internal_fp_t ret;
      internal_fp_t E=hypot(k/T,eta)-mot;
      internal_fp_t arg1=E-y;
      
      ret=k*k*k*k/3/hypot(k/T,eta)/T/(1+exp(arg1));
      //ret=k*k*T*log1p(exp(-arg1));
      //ret=k*k*T*log(1+exp(-arg1));
      
      if (debug) {
        std::cout << "Z: " << k << " " << T << " " << y << " "
                  << eta << " "
                  << mot << " " << ret << std::endl;
      }

      if (!isfinite(ret)) {
        return 0.0;
	//O2SCL_ERR2("Returned not finite result ",
        //"in fermion_rel::deg_pressure_fun().",exc_einval);
      }
  
      return ret;
    }

    /// The integrand for the entropy density for degenerate fermions
    template<class internal_fp_t>
    internal_fp_t deg_entropy_fun(internal_fp_t k, internal_fp_t T,
                                  internal_fp_t y, internal_fp_t eta,
                                  internal_fp_t mot) {
  
      internal_fp_t ret;
      internal_fp_t E=hypot(k/T,eta)-mot;
      internal_fp_t arg1=E-y;

      // If the argument to the exponential is really small, then the
      // value of the integrand is just zero
      if (arg1<-exp_limit) {
	ret=0.0;
	// Otherwise, if the argument to the exponential is still small,
	// then addition of 1 makes us lose precision, so we use an
	// alternative:
      } else if (arg1<-deg_entropy_fac) {
	ret=-k*k*(-1+arg1)*exp(arg1);
      } else {
	internal_fp_t nx=(1+exp(arg1));
        nx=1/nx;
        internal_fp_t arg2=1-nx;
        internal_fp_t t1=nx*log(nx);
        internal_fp_t t2=arg2*log(arg2);
        internal_fp_t t3=t1+t2;
        ret=-k*k*t3;
	//ret=-k*k*(nx*log(nx)+arg2*log(arg2));
      }

      if (!isfinite(ret)) {
	O2SCL_ERR2("Returned not finite result ",
		   "in fermion_rel::deg_entropy_fun().",exc_einval);
      }

      return ret;
    }
    
  };
  
  /** \brief Default integrator for \ref o2scl::fermion_rel_tl

      \note This version uses the multiprecision integrator
      \ref inte_multip_double_exp_boost which automatically
      increases precision in order to achieve accuracy
   */
  template<class fp_t> class fermion_rel_integ_multip :
    public fermion_rel_integ_base {

  public:

    /** \brief Verbosity parameter
     */
    int verbose;
    
    /** \brief Default integrator
     */
    inte_multip_double_exp_boost it;
    
    /** \brief Secondary integrator
     */
    inte_adapt_cern it2;

    /** \brief Desc
     */
    double tol_rel;

    fermion_rel_integ_multip() {
      verbose=0;
      it.pow_tol_func=1.5;
      it2.set_nsub(10000);
      it.err_nonconv=false;
    }
    
    /// The integrand for the density for non-degenerate fermions
    template<class internal_fp_t>
    internal_fp_t density_fun(internal_fp_t u, fp_t y2, fp_t eta2) {
      
      internal_fp_t ret;

      internal_fp_t y=static_cast<internal_fp_t>(y2);
      internal_fp_t eta=static_cast<internal_fp_t>(eta2);
      
      internal_fp_t arg1=u*u+2*eta*u;
      internal_fp_t arg2=eta+u-y;
      internal_fp_t arg3=eta+u;
      if (y-eta-u>exp_limit) {
	ret=arg3*sqrt(arg1);
      } else if (y>u+exp_limit && eta>u+exp_limit) {
	ret=arg3*sqrt(arg1)/(exp(arg2)+1);
      } else {
	ret=arg3*sqrt(arg1)*exp(y)/(exp(arg3)+exp(y));
      }
      
      if (!isfinite(ret)) {
	ret=0;
      }

      return ret;
    }

    /// The integrand for the pressure for non-degenerate fermions
    template<class internal_fp_t>
    internal_fp_t pressure_fun(internal_fp_t u, fp_t y2,
                               fp_t eta2) {

      internal_fp_t ret;
      
      internal_fp_t y=static_cast<internal_fp_t>(y2);
      internal_fp_t eta=static_cast<internal_fp_t>(eta2);
      
      internal_fp_t arg1=u*u+2*eta*u;
      internal_fp_t term1=sqrt(arg1);
      internal_fp_t arg3=eta+u;
      ret=term1*arg3*log1p(exp(y-arg3));
      
      if (!isfinite(ret)) {
	ret=0;
      }

      return ret;
    }

    /// The integrand for the energy density for non-degenerate fermions
    template<class internal_fp_t>
    internal_fp_t energy_fun(internal_fp_t u, fp_t y2,
                             fp_t eta2) {

      internal_fp_t ret;

      internal_fp_t y=static_cast<internal_fp_t>(y2);
      internal_fp_t eta=static_cast<internal_fp_t>(eta2);
      
      internal_fp_t arg1=u*u+2*eta*u;
      internal_fp_t arg2=eta+u-y;
      internal_fp_t arg3=eta+u;
      if (y>u+exp_limit && eta>u+exp_limit) {
	ret=arg3*arg3*sqrt(arg1)/(exp(arg2)+1);
      } else {
	ret=arg3*arg3*sqrt(arg1)*exp(y)/(exp(arg3)+exp(y));
      }
 
      if (!isfinite(ret)) {
	return 0;
      }

      return ret;
    }

    /// The integrand for the entropy density for non-degenerate fermions
    template<class internal_fp_t>
    internal_fp_t entropy_fun(internal_fp_t u, fp_t y2,
                              fp_t eta2) {

      internal_fp_t ret;

      internal_fp_t y=static_cast<internal_fp_t>(y2);
      internal_fp_t eta=static_cast<internal_fp_t>(eta2);
      
      internal_fp_t arg1=u*u+2*eta*u;
      internal_fp_t arg2=eta+u-y;
      internal_fp_t arg3=eta+u;
      internal_fp_t arg4=y-eta-u;
      internal_fp_t arg5=1+exp(arg4);
      internal_fp_t arg6=1+exp(arg2);
      internal_fp_t term1=log(arg5)/arg5;
      internal_fp_t term2=log(arg6)/arg6;
      ret=arg3*sqrt(arg1)*(term1+term2);
  
      if (!isfinite(ret)) {
	return 0;
      }

      return ret;
    }

    /// The integrand for the density for degenerate fermions
    template<class internal_fp_t>
    internal_fp_t deg_density_fun(internal_fp_t k, fp_t T2,
                                  fp_t y2, fp_t eta2,
                                  fp_t mot2, bool debug) {

      internal_fp_t ret;
      
      internal_fp_t T=static_cast<internal_fp_t>(T2);
      internal_fp_t y=static_cast<internal_fp_t>(y2);
      internal_fp_t eta=static_cast<internal_fp_t>(eta2);
      internal_fp_t mot=static_cast<internal_fp_t>(mot2);
      
      internal_fp_t E=hypot(k/T,eta)-mot;
      internal_fp_t arg1=E-y;
      ret=k*k/(1+exp(arg1));

      if (debug) {
        std::cout << k << " " << ret << std::endl;
      }
          
      if (!isfinite(ret)) {
	O2SCL_ERR2("Returned not finite result ",
		   "in fermion_rel_integ_multip::deg_density_fun().",
                   exc_einval);
      }
      return ret;
    }

    /// The integrand for the energy density for degenerate fermions
    template<class internal_fp_t>
    internal_fp_t deg_energy_fun(internal_fp_t k, fp_t T2,
                                 fp_t y2, fp_t eta2,
                                 fp_t mot2) {

      internal_fp_t ret;
      
      internal_fp_t T=static_cast<internal_fp_t>(T2);
      internal_fp_t y=static_cast<internal_fp_t>(y2);
      internal_fp_t eta=static_cast<internal_fp_t>(eta2);
      internal_fp_t mot=static_cast<internal_fp_t>(mot2);
      
      internal_fp_t E=hypot(k/T,eta)-mot;
      internal_fp_t arg1=E-y;
      
      ret=k*k*E*T/(1+exp(arg1));
      
      if (!isfinite(ret)) {
	O2SCL_ERR2("Returned not finite result ",
		   "in fermion_rel_integ_multip::deg_energy_fun().",exc_einval);
      }
  
      return ret;
    }

    /// The integrand for the energy density for degenerate fermions
    template<class internal_fp_t>
    internal_fp_t deg_pressure_fun(internal_fp_t k, fp_t T2,
                                   fp_t y2, fp_t eta2,
                                   fp_t mot2, bool debug) {

      internal_fp_t ret;
      
      internal_fp_t T=static_cast<internal_fp_t>(T2);
      internal_fp_t y=static_cast<internal_fp_t>(y2);
      internal_fp_t eta=static_cast<internal_fp_t>(eta2);
      internal_fp_t mot=static_cast<internal_fp_t>(mot2);
      
      internal_fp_t E=hypot(k/T,eta)-mot;
      internal_fp_t arg1=E-y;
      
      ret=k*k*k*k/3/hypot(k/T,eta)/T/(1+exp(arg1));
      //ret=k*k*T*log1p(exp(-arg1));
      //ret=k*k*T*log(1+exp(-arg1));
      
      if (debug) {
        std::cout << "Z: " << k << " " << T << " " << y << " "
                  << eta << " "
                  << mot << " " << ret << std::endl;
      }

      if (!isfinite(ret)) {
        return 0;
	//O2SCL_ERR2("Returned not finite result ",
        //"in fermion_rel::deg_pressure_fun().",exc_einval);
      }
  
      return ret;
    }

    /// The integrand for the entropy density for degenerate fermions
    template<class internal_fp_t>
    internal_fp_t deg_entropy_fun(internal_fp_t k, fp_t T2,
                                  fp_t y2, fp_t eta2,
                                  fp_t mot2) {
  
      internal_fp_t ret;
      
      internal_fp_t T=static_cast<internal_fp_t>(T2);
      internal_fp_t y=static_cast<internal_fp_t>(y2);
      internal_fp_t eta=static_cast<internal_fp_t>(eta2);
      internal_fp_t mot=static_cast<internal_fp_t>(mot2);
      
      internal_fp_t E=hypot(k/T,eta)-mot;
      internal_fp_t arg1=E-y;

      // If the argument to the exponential is really small, then the
      // value of the integrand is just zero
      if (arg1<-exp_limit) {
	ret=0;
	// Otherwise, if the argument to the exponential is still small,
	// then addition of 1 makes us lose precision, so we use an
	// alternative:
      } else if (arg1<-deg_entropy_fac) {
	ret=-k*k*(-1+arg1)*exp(arg1);
      } else {
	internal_fp_t nx=(1+exp(arg1));
        nx=1/nx;
        internal_fp_t arg2=1-nx;
        internal_fp_t t1=nx*log(nx);
        internal_fp_t t2=arg2*log(arg2);
        internal_fp_t t3=t1+t2;
        ret=-k*k*t3;
	//ret=-k*k*(nx*log(nx)+arg2*log(arg2));
      }

      if (!isfinite(ret)) {
	O2SCL_ERR2("Returned not finite result ",
		   "in fermion_rel_integ_multip::deg_entropy_fun().",
                   exc_einval);
      }

      return ret;
    }
    
    /** \brief Evaluate the density in the nondegenerate limit
     */
    int eval_density(fp_t y, fp_t eta, fp_t &res, fp_t &err) {

      fp_t zero=0;

      if (verbose>1) {
        std::cout << "Calling non-degenerate integrator for density "
                  << "with tolerance: " << tol_rel << std::endl;
      }
      int iret=it.integ_iu_err_multip([this,y,eta](auto &&u) mutable {
        return this->density_fun(u,y,eta); },
        zero,res,err,tol_rel);
      if (iret!=0) {
        iret=it2.integ_iu_err_multip([this,y,eta](auto &&u) mutable {
          return this->density_fun(u,y,eta); },
          zero,res,err,tol_rel);
        if (iret!=0) {
          O2SCL_ERR2("Nondegenerate density failed ",
                     "in fermion_rel_integ_multip::eval_density().",
                     o2scl::exc_efailed);
        }
      }
      return 0;
    }

    /** \brief Evaluate the energy density in the nondegenerate limit
     */
    int eval_energy(fp_t y, fp_t eta, fp_t &res, fp_t &err) {

      fp_t zero=0;
      
      if (verbose>1) {
        std::cout << "Calling non-degenerate integrator for energy density "
                  << "with tolerance: " << tol_rel << std::endl;
      }
      int iret=it.integ_iu_err_multip([this,y,eta](auto &&u) mutable {
        return this->energy_fun(u,y,eta); },
        zero,res,err,tol_rel);
      if (iret!=0) {
        iret=it2.integ_iu_err_multip([this,y,eta](auto &&u) mutable {
          return this->energy_fun(u,y,eta); },
          zero,res,err,tol_rel);
        if (iret!=0) {
          O2SCL_ERR2("Nondegenerate energy density failed ",
                     "in fermion_rel_integ_multip::eval_density().",
                     o2scl::exc_efailed);
        }
      }
      
      return 0;
    }

    /** \brief Evaluate the entropy in the nondegenerate limit
     */
    int eval_entropy(fp_t y, fp_t eta, fp_t &res, fp_t &err) {

      fp_t zero=0;
      
      if (verbose>1) {
        std::cout << "Calling non-degenerate integrator for entropy density "
                  << "with tolerance: " << tol_rel << std::endl;
      }
      int iret=it.integ_iu_err_multip([this,y,eta](auto &&u) mutable {
        return this->entropy_fun(u,y,eta); },
        zero,res,err,tol_rel);
      if (iret!=0) {
        iret=it2.integ_iu_err_multip([this,y,eta](auto &&u) mutable {
          return this->entropy_fun(u,y,eta); },
          zero,res,err,tol_rel);
        if (iret!=0) {
          O2SCL_ERR2("Nondegenerate entropy failed ",
                     "in fermion_rel_integ_multip::eval_density().",
                     o2scl::exc_efailed);
        }
      }
      
      return 0;
    }

    /** \brief Evaluate the pressure in the nondegenerate limit
     */
    int eval_pressure(fp_t y, fp_t eta, fp_t &res, fp_t &err) {

      fp_t zero=0;
      
      if (verbose>1) {
        std::cout << "Calling non-degenerate integrator for pressure "
                  << "with tolerance: " << tol_rel << std::endl;
      }
      int iret=it.integ_iu_err_multip([this,y,eta](auto &&u) mutable {
        return this->pressure_fun(u,y,eta); },
        zero,res,err,tol_rel);
      if (iret!=0) {
        iret=it2.integ_iu_err_multip([this,y,eta](auto &&u) mutable {
          return this->pressure_fun(u,y,eta); },
          zero,res,err,tol_rel);
        if (iret!=0) {
          O2SCL_ERR2("Nondegenerate pressure failed ",
                     "in fermion_rel_integ_multip::eval_density().",
                     o2scl::exc_efailed);
        }
      }
      
      return 0;
    }

    /** \brief Evaluate the density in the nondegenerate limit
     */
    int eval_deg_density(fp_t T, fp_t y, fp_t eta, fp_t mot,
                         fp_t ul, fp_t &res, fp_t &err) {

      fp_t zero=0;
      
      if (verbose>1) {
        std::cout << "Calling degenerate integrator for density "
                  << "with tolerance: " << tol_rel << std::endl;
      }
      int iret=it.integ_err_multip([this,T,y,eta,mot](auto &&k) mutable {
        return this->deg_density_fun(k,T,y,eta,mot,false); },
        zero,ul,res,err,tol_rel);
      if (iret!=0) {
        O2SCL_ERR("Deg density failed.",o2scl::exc_efailed);
      }
      
      return 0;
    }
    
    /** \brief Desc
     */
    int eval_deg_energy(fp_t T, fp_t y, fp_t eta, fp_t mot,
                        fp_t ul, fp_t &res, fp_t &err) {

      fp_t zero=0;
      
      if (verbose>1) {
        std::cout << "Calling degenerate integrator for energy density "
                  << "with tolerance: " << tol_rel << std::endl;
      }
      int iret=it.integ_err_multip([this,T,y,eta,mot](auto &&k) mutable {
        return this->deg_energy_fun(k,T,y,eta,mot); },
        zero,ul,res,err,tol_rel);
      if (iret!=0) {
        O2SCL_ERR("Deg energy failed.",o2scl::exc_efailed);
      }
      
      return 0;
    }
    
    /** \brief Desc
     */
    int eval_deg_entropy(fp_t T, fp_t y, fp_t eta, fp_t mot,
                         fp_t ll, fp_t ul, fp_t &res, fp_t &err) {

      if (verbose>1) {
        std::cout << "Calling degenerate integrator for entropy density "
                  << "with tolerance: " << tol_rel << std::endl;
      }
      int iret=it.integ_err_multip([this,T,y,eta,mot](auto &&k) mutable {
        return this->deg_entropy_fun(k,T,y,eta,mot); },
        ll,ul,res,err,tol_rel);
      if (iret!=0) {
        O2SCL_ERR("Deg entropy failed.",o2scl::exc_efailed);
      }
      
      return 0;
    }
    
    /** \brief Desc
     */
    int eval_deg_pressure(fp_t T, fp_t y, fp_t eta, fp_t mot,
                         fp_t ul, fp_t &res, fp_t &err) {

      fp_t zero=0;
      
      if (verbose>1) {
        std::cout << "Calling degenerate integrator for pressure "
                  << "with tolerance: " << tol_rel << std::endl;
      }
      int iret=it.integ_err_multip([this,T,y,eta,mot](auto &&k) mutable {
        return this->deg_pressure_fun(k,T,y,eta,mot,false); },
        zero,ul,res,err,tol_rel);
      if (iret!=0) {
        O2SCL_ERR("Deg pressure failed.",o2scl::exc_efailed);
      }
      
      return 0;
    }
    
    
  };
    
  /** \brief Default integrator for \ref o2scl::fermion_rel_tl
   */
  template<class func_t, class fp_t> class fermion_rel_integ :
    public fermion_rel_integ_base {

  protected:
    
  public:

    /// Nondegenerate integrator
    inte_qagiu_gsl<> nit;

    /// Degenerate integrator
    inte_qag_gsl<> dit;

    /** \brief Evalulate the density integral in the nondegenerate limit
     */
    int eval_density(fp_t y, fp_t eta, fp_t &res, fp_t &err) {
      func_t mfd=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t)>
                           (&fermion_rel_integ<func_t,
                            fp_t>::density_fun<fp_t>),
                           this,std::placeholders::_1,y,eta);
      fp_t zero=0;
      int iret=nit.integ_iu_err(mfd,0,res,err);
      if (iret!=0 && err_nonconv) {
        O2SCL_ERR2("Density integration failed in ",
                  "fermion_rel_integ::eval_density().",
                  o2scl::exc_efailed);
      }
      return iret;
    }

    /** \brief Evalulate the density integral in the degenerate limit
     */
    int eval_deg_density(fp_t T, fp_t y, fp_t eta, fp_t mot,
                         fp_t ul, fp_t &res, fp_t &err) {
      func_t mfd=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t,fp_t,fp_t,bool)>
                           (&fermion_rel_integ<func_t,
                            fp_t>::deg_density_fun<fp_t>),
                           this,std::placeholders::_1,T,y,eta,mot,false);
      fp_t zero=0;
      int iret=dit.integ_err(mfd,zero,ul,res,err);
      if (iret!=0 && err_nonconv) {
        O2SCL_ERR2("Degenerate density integration failed in ",
                  "fermion_rel_integ::eval_density().",
                  o2scl::exc_efailed);
      }
      return iret;
    }

    /** \brief Evalulate the energy density integral in the 
        nondegenerate limit
     */
    int eval_energy(fp_t y, fp_t eta, fp_t &res, fp_t &err) {
      func_t mfd=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t)>
                           (&fermion_rel_integ<func_t,
                            fp_t>::energy_fun<fp_t>),
                           this,std::placeholders::_1,y,eta);
      fp_t zero=0;
      int iret=nit.integ_iu_err(mfd,zero,res,err);
      if (iret!=0 && err_nonconv) {
        O2SCL_ERR2("Energy integration failed in ",
                  "fermion_rel_integ::eval_density().",
                  o2scl::exc_efailed);
      }
      return iret;
    }

    /** \brief Evalulate the energy density integral in the degenerate limit
     */
    int eval_deg_energy(fp_t T, fp_t y, fp_t eta, fp_t mot,
                        fp_t ul, fp_t &res, fp_t &err) {
      func_t mfd=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t,fp_t,fp_t)>
                           (&fermion_rel_integ<func_t,
                            fp_t>::deg_energy_fun<fp_t>),
                           this,std::placeholders::_1,T,y,eta,mot);
      fp_t zero=0;
      int iret=dit.integ_err(mfd,0,ul,res,err);
      if (iret!=0 && err_nonconv) {
        O2SCL_ERR2("Degenerate energy integration failed in ",
                  "fermion_rel_integ::eval_density().",
                  o2scl::exc_efailed);
      }
      return iret;
    }

    /** \brief Evalulate the entropy integral in the nondegenerate limit
     */
    int eval_entropy(fp_t y, fp_t eta, fp_t &res, fp_t &err) {
      func_t mfd=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t)>
                           (&fermion_rel_integ<func_t,
                            fp_t>::entropy_fun<fp_t>),
                           this,std::placeholders::_1,y,eta);
      fp_t zero=0;
      int iret=nit.integ_iu_err(mfd,0,res,err);
      if (iret!=0 && err_nonconv) {
        O2SCL_ERR2("Entropy integration failed in ",
                  "fermion_rel_integ::eval_density().",
                  o2scl::exc_efailed);
      }
      return iret;
    }

    /** \brief Evalulate the entropy integral in the degenerate limit
     */
    int eval_deg_entropy(fp_t T, fp_t y, fp_t eta, fp_t mot,
                         fp_t ll, fp_t ul, fp_t &res, fp_t &err) {
      func_t mfd=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t,fp_t,fp_t)>
                           (&fermion_rel_integ<func_t,
                            fp_t>::deg_entropy_fun<fp_t>),
                           this,std::placeholders::_1,T,y,eta,mot);
      fp_t zero=0;
      int iret=dit.integ_err(mfd,0,ul,res,err);
      if (iret!=0 && err_nonconv) {
        O2SCL_ERR2("Degnerate entropy integration failed in ",
                  "fermion_rel_integ::eval_density().",
                  o2scl::exc_efailed);
      }
      return iret;
    }

    /** \brief Evalulate the entropy integral in the nondegenerate limit
     */
    int eval_pressure(fp_t y, fp_t eta, fp_t &res, fp_t &err) {
      func_t mfd=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t)>
                           (&fermion_rel_integ<func_t,
                            fp_t>::pressure_fun<fp_t>),
                           this,std::placeholders::_1,y,eta);
      fp_t zero=0;
      int iret=nit.integ_iu_err(mfd,0,res,err);
      if (iret!=0 && err_nonconv) {
        O2SCL_ERR2("Pressure integration failed in ",
                  "fermion_rel_integ::eval_density().",
                  o2scl::exc_efailed);
      }
      return iret;
    }

    /** \brief Evalulate the entropy integral in the degenerate limit
     */
    int eval_deg_pressure(fp_t T, fp_t y, fp_t eta, fp_t mot,
                          fp_t ul, fp_t &res, fp_t &err) {
      func_t mfd=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fp_t,fp_t,
                                            fp_t,bool)>
                           (&fermion_rel_integ<func_t,
                            fp_t>::deg_pressure_fun<fp_t>),
                           this,std::placeholders::_1,T,y,eta,mot,false);
      fp_t zero=0;
      int iret=dit.integ_err(mfd,0,ul,res,err);
      if (iret!=0 && err_nonconv) {
        O2SCL_ERR2("Degenerate pressure integration failed in ",
                  "fermion_rel_integ::eval_density().",
                  o2scl::exc_efailed);
      }
      return iret;
    }

  };
  
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
      \f$ \eta \equiv m^{*}/T \f$. Then, 
      \f$ p/T=\sqrt{u^2+2 \eta u} \f$
      \f$ E/T = \mathrm{mx+u} \f$ and 
      \f$ p/T^2 dp = 2(\eta+u) du \f$
      The density integrand is 
      \f[
      \left(\eta+u\right) \sqrt{u^2+2 (\eta) u}
      \left(\frac{e^{y}}{e^{\eta+u}+e^{y}}\right) \, , 
      \f]
      the energy integrand is 
      \f[
      \left(\eta+u\right)^2 \sqrt{u^2+2 (\eta) u}
      \left(\frac{e^{y}}{e^{\eta+u}+e^{y}}\right) \, ,
      \f]
      and the entropy integrand is 
      \f[
      \left(\eta+u\right) \sqrt{u^2+2 (\eta) u} 
      \left(t_1+t_2\right) \, ,
      \f]
      where 
      \f{eqnarray*}
      t_1 &=& \log \left(1+e^{y-\eta-u}\right)/
      \left(1+e^{y-\eta-u}\right) \nonumber \\
      t_2 &=& \log \left(1+e^{\eta+u-y}\right)/
      \left(1+e^{\eta+u-y}\right) \, .
      \f}

      \hline 
      <b>Accuracy:</b>

      The default settings for for this class give an accuracy of at
      least 1 part in \f$ 10^6 \f$ (and frequently better than this).

      When the integrators provide numerical uncertainties, these
      uncertainties are stored in \ref unc. In the case of
      calc_density() and pair_density(), the uncertainty from the
      numerical accuracy of the solver is not included. There is also
      a relatively small inaccuracy due to the mathematical evaluation
      of the integrands which is not included in \ref unc.)
      
      One can improve the accuracy to within 1 part in \f$ 10^{10} \f$
      using \code fermion_rel rf(1.0,2.0); rf.upper_limit_fac=40.0;
      rf.dit->tol_abs=1.0e-13; rf.dit->tol_rel=1.0e-13;
      rf.nit->tol_abs=1.0e-13; rf.nit->tol_rel=1.0e-13;
      rf.density_root->tol_rel=1.0e-10; \endcode which decreases the
      both the relative and absolute tolerances for both the
      degenerate and non-degenerate integrators and improves the
      accuracy of the solver which determines the chemical potential
      from the density. Of course if these tolerances are too small,
      the calculation may fail.
      
      \verbatim embed:rst
      
      .. todo::
      
         In class fermion_rel_tl:
      
         - Future: I had to remove the shared_ptr stuff because the
           default algorithm types don't support multiprecision, but it
           might be nice to restore the shared_ptr mechanism somehow.
           
         - Future: The expressions which appear in in the integrand
           functions density_fun(), etc. could likely be improved,
           especially in the case where \ref o2scl::part::inc_rest_mass
           is <tt>false</tt>. There should not be a need to check if
           <tt>ret</tt> is finite.
         
         - Future: It appears this class doesn't compute the
           uncertainty in the chemical potential or density with
           calc_density(). This could be fixed.
           
         - Future: I'd like to change the lower limit on the entropy
           integration, but the value in the code at the moment (stored
           in <tt>ll</tt>) makes bm_part2.cpp worse.
           
         - Future: The function pair_mu() should set the antiparticle
           integrators as done in fermion_deriv_rel.
         
      \endverbatim
  */
  template<class fermion_t=fermion_tl<double>,
	   class fd_inte_t=class o2scl::fermi_dirac_integ_gsl,
	   class be_inte_t=o2scl::bessel_K_exp_integ_gsl,
           class inte_t=fermion_rel_integ<funct,double>,
	   class density_root_t=root_cern<>,
	   class root_t=root_cern<>, class func_t=funct,
	   class fp_t=double>
  class fermion_rel_tl :
    public fermion_thermo_tl<fermion_t,fd_inte_t,be_inte_t,root_t,
			     func_t,fp_t> {
    
  public:

    /// The integrator 
    inte_t fri;
    
    /// \name Numerical parameters
    //@{
    /** \brief If true, call the error handler when convergence 
	fails (default true)
    */
    bool err_nonconv;
    
    /** \brief The smallest value of \f$ (\mu-m)/T \f$ for which 
	integration is used (default -4.0)
    */
    fp_t min_psi;

    /** \brief The critical degeneracy at which to switch integration 
	techniques (default 2.0)
    */
    fp_t deg_limit;
    
    /** \brief The limit for exponentials to ensure integrals are finite 
	(default 200.0)
    */
    fp_t exp_limit;

    /// The factor for the degenerate upper limits (default 20.0)
    fp_t upper_limit_fac;

    /// A factor for the degenerate entropy integration (default 30.0)
    fp_t deg_entropy_fac;

    /// Verbosity parameter (default 0)
    int verbose;

    /// If true, use expansions for extreme conditions (default true)
    bool use_expansions;

    /// Tolerance for expansions (default \f$ 10^{-14} \f$)
    fp_t tol_expan;

    /// If true, verify the thermodynamic identity (default false)
    bool verify_ti;
    
    /// Value for verifying the thermodynamic identity
    fp_t therm_ident;

    /// Alternate solver
    o2scl::root_brent_gsl<func_t,fp_t> alt_solver;
    //@}

    /// Storage for the uncertainty
    fermion_t unc;

    /// Create a fermion with mass \c m and degeneracy \c g
    fermion_rel_tl() {
      
      deg_limit=2.0;
      
      exp_limit=200.0;
      upper_limit_fac=20.0;
      deg_entropy_fac=30.0;
      min_psi=-4.0;
      err_nonconv=true;
      use_expansions=true;
      verbose=0;
      last_method=0;

      //dit=&def_dit;
      //nit=&def_nit;
      density_root=&def_density_root;
      
      density_root->tol_rel=4.0e-7;

      tol_expan=1.0e-14;
      verify_ti=false;
      therm_ident=0.0;

      alt_solver.err_nonconv=false;
      // AWS, 6/24/21: This appears to give better results than
      // either test_form=0 or test_form=1.
      alt_solver.test_form=2;
    }

    virtual ~fermion_rel_tl() {
    }

    /// The solver for calc_density()
    root<func_t,func_t,fp_t> *density_root;

    /// The default solver for the chemical potential given the density
    density_root_t def_density_root;
    
    /// Return string denoting type ("fermion_rel")
    virtual const char *type() { return "fermion_rel"; }

    /** \brief An integer indicating the last numerical method used

	In all functions
	- 0: no previous calculation or last calculation failed

	In \ref nu_from_n():
	- 1: default solver
	- 2: default solver with smaller tolerances
	- 3: bracketing solver

	In \ref calc_mu():
	- 4: non-degenerate expansion
	- 5: degenerate expansion
	- 6: exact integration, non-degenerate integrands
	- 7: exact integration, degenerate integrands, lower limit
	on entropy integration
	- 8: exact integration, degenerate integrands, full
	entropy integration
	- 9: T=0 result

	In \ref calc_density(), the integer is a two-digit
	number. The first digit (1 to 3) is the method used by \ref
	nu_from_n() and the second digit is one of
	- 1: nondegenerate expansion
	- 2: degenerate expansion
	- 3: exact integration, non-degenerate integrands
	- 4: exact integration, degenerate integrands, lower limit
	on entropy integration
	- 5: exact integration, degenerate integrands, full
	entropy integration
	If \ref calc_density() uses the T=0 code, then
	last_method is 40. 

	In \ref pair_mu(), the integer is a three-digit number.
	The third digit is always 0 (to ensure a value of last_method
	which is unique from the other values reported from other
	functions as described above). The first digit is the method
	used for particles from \ref calc_mu() above and the
	second digit is the method used for antiparticles. 

	In \ref pair_density(), the integer is a four-digit
	number. The first digit is from the list below and the
	remaining three digits, if nonzero, are from \ref
	pair_mu().
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
    int nu_from_n(fermion_t &f, fp_t temper) {

      last_method=0;
      
      fp_t nex;

      // Try to ensure a good initial guess

      nex=f.nu/temper;
      fp_t y=solve_fun(nex,f,temper);
      if (verbose>1) {
	std::cout << "fermion_rel::nu_from_n(): " 
		  << "initial guess " << nex << std::endl;
      }

      if (y>1.0-1.0e-6) {
	fp_t scale=f.ms;
	if (temper>scale) scale=temper;
	for(size_t i=0;i<10;i++) {
	  if (nex<0.0) nex+=scale*1.0e5;
	  else nex*=10.0;
	  y=solve_fun(nex,f,temper);
	  if (y<1.0-1.0e-6) i=10;
	}
	if (verbose>1) {
	  std::cout << "fermion_rel::nu_from_n(): "
		    << "adjusted guess to " << nex << std::endl;
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
	  std::cout << "fermion_rel::nu_from_n(): "
		    << "adjusted guess (try 2) to "
		    << nex << std::endl;
	}
      }
  
      // If neither worked, call the error handler
      if (y==1.0 || !isfinite(y)) {
	O2SCL_CONV2_RET("Couldn't find reasonable initial guess in ",
			"fermion_rel::nu_from_n().",exc_einval,
			this->err_nonconv);
      }

      // Perform full solution
      func_t mf=std::bind(std::mem_fn<fp_t(fp_t,fermion_t &,fp_t)>
			  (&fermion_rel_tl<fermion_t,fd_inte_t,be_inte_t,
			   inte_t,density_root_t,
			   root_t,func_t,fp_t>::solve_fun),
			  this,std::placeholders::_1,std::ref(f),temper);

      // The default o2scl::root object is of type root_cern,
      // and this solver has problems when the root is near 0.
      // Below, we switch to a root_brent_gsl object in the case
      // that the default solver fails.
  
      bool drec=density_root->err_nonconv;
      density_root->err_nonconv=false;
      int ret=density_root->solve(nex,mf);
      last_method=1;

      if (ret!=0) {
    
	if (verbose>1) {
	  std::cout << "nu_from_n(): density_root failed x="
		    << nex << " ." << std::endl;
	}
	O2SCL_CONV2_RET("Density solver failed in ",
			"fermion_rel::nu_from_n().",exc_efailed,
			this->err_nonconv);
      }

      density_root->err_nonconv=drec;

      f.nu=nex*temper;

      return success;
    }

    /** \brief Calculate properties as function of chemical potential
     */
    int calc_mu(fermion_t &f, fp_t temper) {

      if (verbose>1) {
	std::cout << "calc_mu(): start."
		  << std::endl;
      }
      
      last_method=0;
      
      fp_t zero=0;

      // -----------------------------------------------------------------
      // Handle T<=0

      if (temper<zero) {
	O2SCL_ERR2("Temperature less than zero in ",
		   "fermion_rel::calc_mu().",exc_einval);
      }
      if (temper==zero) {
	this->calc_mu_zerot(f);
	last_method=9;
	return 0;
      }

      if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

      // Compute the degeneracy parameter
  
      bool deg=true;
      fp_t psi;
      if (f.inc_rest_mass) {
	psi=(f.nu-f.ms)/temper;
      } else {
        if (f.non_interacting) {
          psi=f.nu/temper;
        } else {
          psi=(f.nu+(f.m-f.ms))/temper;
        }
      }
      if (psi<deg_limit) deg=false;
      
      if (verbose>1) {
	std::cout << "calc_mu(): psi,deg,deg_limit: " << psi << " "
		  << deg << " " << deg_limit << std::endl;
      }
      
      // Try the non-degenerate expansion if psi is small enough
      if (use_expansions && psi<min_psi) {
	bool acc=this->calc_mu_ndeg(f,temper,tol_expan);
        //std::cout << "cmn_cmu: " << f.nu << " " << temper << " "
        //<< tol_expan << " " << acc << std::endl;
	if (verbose>1) {
	  std::cout << "calc_mu(): non-deg expan " << acc
		    << std::endl;
	}
	if (acc) {
	  unc.n=f.n*tol_expan;
	  unc.ed=f.ed*tol_expan;
	  unc.pr=f.pr*tol_expan;
	  unc.en=f.en*tol_expan;
	  last_method=4;
	  return 0;
	}
      }

      // Try the degenerate expansion if psi is large enough
      if (use_expansions && psi>20.0) {
	bool acc=this->calc_mu_deg(f,temper,tol_expan);
	if (verbose>1) {
	  std::cout << "calc_mu(): deg expan " << acc
		    << std::endl;
	}
	if (acc) {
	  unc.n=f.n*tol_expan;
	  unc.ed=f.ed*tol_expan;
	  unc.pr=f.pr*tol_expan;
	  unc.en=f.en*tol_expan;
	  last_method=5;
	  return 0;
	}
      }

      if (!deg) {

	// If the temperature is large enough, perform the full integral

        fp_t y, eta;
        if (f.inc_rest_mass) {
          y=f.nu/temper;
        } else {
          y=(f.nu+f.m)/temper;
        }
        eta=f.ms/temper;
        
	fp_t prefac=f.g*pow(temper,3)/2/this->pi2;

	// Compute the number density
    
	if (verbose>1) {
	  std::cout << "calc_mu(): non-deg number density:"
		    << std::endl;
	}

        fri.eval_density(y,eta,f.n,unc.n);
        f.n*=prefac;
        unc.n*=prefac;

	// Compute the energy density

	if (verbose>1) {
	  std::cout << "calc_mu(): non-deg energy density:"
		    << std::endl;
	}
	
        fri.eval_energy(y,eta,f.ed,unc.ed);
        f.ed*=prefac*temper;
        unc.ed*=prefac*temper;
	if (!f.inc_rest_mass) f.ed-=f.n*f.m;
        
	// Compute the entropy

	if (verbose>1) {
	  std::cout << "calc_mu(): non-deg entropy:"
		    << std::endl;
	}
	
        fri.eval_entropy(y,eta,f.en,unc.en);
        f.en*=prefac;
        unc.en*=prefac;

        if (verify_ti) {
          // Compute the pressure

          if (verbose>1) {
            std::cout << "calc_mu(): non-deg pressure:"
                      << std::endl;
          }
          
          int iret=fri.eval_pressure(y,eta,f.pr,unc.pr);
          
          f.pr*=prefac*temper;
          unc.pr*=prefac*temper;
        }
        
	if (verbose>1) {
	  std::cout << "calc_mu(): non-deg integrals done."
		    << std::endl;
	}
	
	last_method=6;

      } else {

	// Otherwise, apply a degenerate approximation, by making the
	// upper integration limit finite

        fp_t y=f.nu/temper;
        fp_t eta=f.ms/temper;
        fp_t mot;
        if (f.inc_rest_mass) {
          mot=0;
        } else {
          mot=f.m/temper;
        }

	fp_t prefac=f.g/2/this->pi2;
    
	// Compute the upper limit for degenerate integrals

	fp_t arg;
	if (f.inc_rest_mass) {
	  arg=pow(upper_limit_fac*temper+f.nu,2)-f.ms*f.ms;
	} else {
	  arg=pow(upper_limit_fac*temper+f.nu+f.m,2)-f.ms*f.ms;
	}
	fp_t ul;
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
	  return 0;
	}
    
	// Compute the number density

	if (verbose>1) {
	  std::cout << "calc_mu(): deg number density, ul, ulf: "
                    << ul << " " << upper_limit_fac << std::endl;
	}

        fri.eval_deg_density(temper,y,eta,mot,ul,f.n,unc.n);
        f.n*=prefac;
        unc.n*=prefac;

	// Compute the energy density

	if (verbose>1) {
	  std::cout << "calc_mu(): deg energy density."
		    << std::endl;
	}
	
        fri.eval_deg_energy(temper,y,eta,mot,ul,f.ed,unc.ed);
        f.ed*=prefac;
        unc.ed*=prefac;

	// Compute the lower limit for the entropy integration

	fp_t ll;
	if (f.inc_rest_mass) {
	  arg=pow(-upper_limit_fac*temper+f.nu,2)-f.ms*f.ms;
	  if (arg>0.0 && (f.ms-f.nu)/temper<-upper_limit_fac) {
	    ll=sqrt(arg);
	  } else {
	    ll=-1.0;
	  }
	} else {
	  arg=pow(-upper_limit_fac*temper+f.nu+f.m,2)-f.ms*f.ms;
	  if (arg>0.0 && (f.ms-f.nu-f.m)/temper<-upper_limit_fac) {
	    ll=sqrt(arg);
	  } else {
	    ll=-1.0;
	  }
	}

	// Compute the entropy

	if (verbose>1) {
	  std::cout << "calc_mu(): deg entropy."
		    << std::endl;
	}
	
	if (ll>0.0) {
          fri.eval_deg_entropy(temper,y,eta,mot,ll,ul,f.en,unc.en);
          
	  //f.en=dit->integ(mfs,ll,ul);
	  last_method=7;
	} else {
          fri.eval_deg_entropy(temper,y,eta,mot,0,ul,f.en,unc.en);

	  //f.en=dit->integ(mfs,0.0,ul);
	  last_method=8;
	}
        f.en*=prefac;
        unc.en*=prefac;

        if (verify_ti) {
          // Compute the pressure
          
          if (verbose>1) {
            std::cout << "calc_mu(): deg pressure."
                      << std::endl;
          }

          fri.eval_deg_pressure(temper,y,eta,mot,ul,f.pr,unc.pr);
          f.pr*=prefac;
          unc.pr*=prefac;
        }

	if (verbose>1) {
	  std::cout << "calc_mu(): deg integrals done."
		    << std::endl;
	}
	
	
      }

      // Compute the pressure using the thermodynamic identity

      if (verify_ti==false) {
        f.pr=-f.ed+temper*f.en+f.nu*f.n;
        unc.pr=sqrt(unc.ed*unc.ed+temper*unc.en*temper*unc.en+
                    f.nu*unc.n*f.nu*unc.n);
      }

      return 0;
    }

    /** \brief Calculate properties as function of density

	This function uses the current value of \c nu (or \c mu if the
	particle is non interacting) for an initial guess to solve for
	the chemical potential. If this guess is too small, then this
	function may fail.

        \verbatim embed:rst

        .. todo::

        In function fermion_rel_tl::calc_density()

        - Future: There is still quite a bit of code duplication
        between this function and \ref calc_mu() .

        \endverbatim
    */
    int calc_density(fermion_t &f, fp_t temper) {

      last_method=0;
      
      // The code below may modify the density, which is confusing to
      // the user, so we store it here so we can restore it before
      // this function returns
      fp_t density_temp=f.n;
  
      // -----------------------------------------------------------------
      // Handle T<=0

      if (temper<0.0) {
	O2SCL_ERR2("Temperature less than zero in ",
		   "fermion_rel::calc_density().",
		   exc_einval);
      }
      if (temper==0.0) {
	last_method=40;
	this->calc_density_zerot(f);
	return 0;
      }

#if !O2SCL_NO_RANGE_CHECK
      // This may not be strictly necessary, because it should be clear
      // that this function will produce gibberish if the density is not
      // finite, but I've found this extra checking of the inputs useful
      // for debugging.
      if (!isfinite(f.n)) {
	O2SCL_ERR2("Density not finite in ",
		   "fermion_rel::calc_density().",exc_einval);
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
      fp_t psi;
      if (f.inc_rest_mass) {
	psi=(f.nu-f.ms)/temper;
      } else {
        if (f.non_interacting) {
          psi=f.nu/temper;
        } else {
          psi=(f.nu+(f.m-f.ms))/temper;
        }
      }
      if (psi<deg_limit) deg=false;

      // Try the non-degenerate expansion if psi is small enough
      if (use_expansions && psi<min_psi) {
	bool acc=this->calc_mu_ndeg(f,temper,tol_expan);
	if (acc) {
	  unc.ed=f.ed*tol_expan;
	  unc.pr=f.pr*tol_expan;
	  unc.en=f.en*tol_expan;
	  f.n=density_temp;
	  last_method+=1;
	  return 0;
	}
      }
  
      // Try the degenerate expansion if psi is large enough
      if (use_expansions && psi>20.0) {
	bool acc=this->calc_mu_deg(f,temper,tol_expan);
	if (acc) {
	  unc.n=f.n*tol_expan;
	  unc.ed=f.ed*tol_expan;
	  unc.pr=f.pr*tol_expan;
	  unc.en=f.en*tol_expan;
	  f.n=density_temp;
	  last_method+=2;
	  return 0;
	}
      }

      if (!deg) {
    
        fp_t y, eta;
        if (f.inc_rest_mass) {
          y=f.nu/temper;
        } else {
          y=(f.nu+f.m)/temper;
        }
        eta=f.ms/temper;

        fp_t prefac=f.g*pow(temper,4.0)/2.0/this->pi2;
        
        fri.eval_energy(y,eta,f.ed,unc.ed);
        f.ed*=prefac;
        unc.ed*=prefac;
	if (!f.inc_rest_mass) f.ed-=f.n*f.m;
        
        prefac=f.g*pow(temper,3.0)/2.0/this->pi2;

        fri.eval_entropy(y,eta,f.en,unc.en);
        f.en*=prefac;
        unc.en*=prefac;
        
	last_method+=3;

      } else {

        fp_t y=f.nu/temper;
        fp_t eta=f.ms/temper;
        fp_t mot;
        if (f.inc_rest_mass) {
          mot=0;
        } else {
          mot=f.m/temper;
        }

	fp_t arg;
	if (f.inc_rest_mass) {
	  arg=pow(upper_limit_fac*temper+f.nu,2)-f.ms*f.ms;
	} else {
	  arg=pow(upper_limit_fac*temper+f.nu+f.m,2)-f.ms*f.ms;
	}
	fp_t ul;
	if (arg>0.0) {
      
	  ul=sqrt(arg);
      
	  fp_t ll;
	  if (f.inc_rest_mass) {
	    arg=pow(-upper_limit_fac*temper+f.nu,2)-f.ms*f.ms;
	    if (arg>0.0 && (f.ms-f.nu)/temper<-upper_limit_fac) {
	      ll=sqrt(arg);
	    } else {
	      ll=-1.0;
	    }
	  } else {
	    arg=pow(-upper_limit_fac*temper+f.nu+f.m,2)-f.ms*f.ms;
	    if (arg>0.0 && (f.ms-f.nu-f.m)/temper<-upper_limit_fac) {
	      ll=sqrt(arg);
	    } else {
	      ll=-1.0;
	    }
	  }

          fp_t prefac=f.g/2.0/this->pi2;
          fri.eval_deg_energy(temper,y,eta,mot,ul,f.ed,unc.ed);
          f.ed*=prefac;
          unc.ed*=prefac;
          
	  if (ll>0.0) {
            fri.eval_deg_entropy(temper,y,eta,mot,ll,ul,f.en,unc.en);
	    last_method+=4;
	  } else {
            fri.eval_deg_entropy(temper,y,eta,mot,0.0,ul,f.en,unc.en);
	    last_method+=5;
	  }
          f.en*=prefac;
          unc.en*=prefac;

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
	chemical potential
    */
    void pair_mu(fermion_t &f, fp_t temper) {

      last_method=0;
      
      if (f.non_interacting) { f.nu=f.mu; f.ms=f.m; }

      // AWS: 6/26/21: Note that when the last argument is true, the
      // function calc_mu_ndeg() includes antiparticles. However, this
      // code caused problems for low densities. Additionally, I'm not
      // sure that the value of last_method=9 here is unambiguous.

      if (false && use_expansions) {
	if (this->calc_mu_ndeg(f,temper,tol_expan,true)) {
	  unc.n=tol_expan*f.n;
	  unc.ed=tol_expan*f.ed;
	  unc.en=tol_expan*f.en;
	  unc.pr=tol_expan*f.pr;
	  last_method=9;
	  return;
	}
      }

      
      fermion_t antip(f.m,f.g);
      f.anti(antip);

      // Particles
      calc_mu(f,temper);
      fp_t unc_n=unc.n;
      fp_t unc_pr=unc.pr;
      fp_t unc_ed=unc.ed;
      fp_t unc_en=unc.en;

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
      if (verbose>0) {
        std::cout << "pair_mu(), particles, antiparticles, total: "
                  << f.n << " " << antip.n << " "
                  << f.n-antip.n << std::endl;
      }
      f.n-=antip.n;
      f.pr+=antip.pr;
      f.en+=antip.en;

      // Add up uncertainties
      unc.n=hypot(unc.n,unc_n);
      unc.ed=hypot(unc.ed,unc_ed);
      unc.pr=hypot(unc.pr,unc_pr);
      unc.en=hypot(unc.ed,unc_en);

      return;
    }

    /** \brief Calculate thermodynamic properties with antiparticles
	from the density

        \verbatim embed:rst
        
        .. todo::

        In function pair_density():

        - This actually works for negative densities some of the
        time, but the solver probably doesn't work as well there and
        we need to document the density expectations for this 
        function.

        \endverbatim
    */
    int pair_density(fermion_t &f, fp_t temper) {
      
      if (verbose>0) {
        std::cout << "Value of verbose greater than zero in "
                  << "fermion_rel::pair_density()." << std::endl;
        std::cout << "Density: " << f.n << " temperature: "
                  << temper << std::endl;
        if (verbose>1) {
          std::cout << "Setting solver verbose parameters to 1."
                    << std::endl;
          density_root->verbose=1;
          alt_solver.verbose=1;
        }
      }
      last_method=0;
      
      // -----------------------------------------------------------------
      // Handle T<=0

      if (temper<=0.0) {
        if (verbose>0) {
          std::cout << "Value of T<=0, so using zero temperature code."
                    << std::endl;
        }
	this->calc_density_zerot(f);
	last_method=1000;
	return success;
      }

      // Store the input density
      fp_t density_match=f.n;
  
      if (f.non_interacting==true) { f.nu=f.mu; f.ms=f.m; }

      // Store the initial guess for the chemical potential
      fp_t initial_guess=f.nu;

      // We always work with mu/T instead of mu directly
      fp_t nex=f.nu/temper;
      
      func_t mf=std::bind(std::mem_fn<fp_t(fp_t,fp_t,fermion_t &,fp_t,bool)>
			  (&fermion_rel_tl<fermion_t,fd_inte_t,be_inte_t,
			   inte_t,density_root_t,
			   root_t,func_t,fp_t>::pair_fun),
			  this,std::placeholders::_1,density_match,
                          std::ref(f),temper,false);

      // Begin by trying the user-specified guess
      bool drec=density_root->err_nonconv;
      density_root->err_nonconv=false;
      int ret=density_root->solve(nex,mf);
      if (ret==0) {
        if (verbose>0) {
          std::cout << "Initial solver succeeded." << std::endl;
        }
        // If that worked, set last_method
	last_method=2000;
      }
        
      if (ret!=0) {
        
        // If that failed, try bracketing the root

	// (max doesn't work with boost::multiprecision?)
	fp_t lg;
	if (abs(f.nu)>f.ms) lg=abs(f.nu);
	lg=f.ms;

        if (verbose>0) {
          std::cout << "Initial solver returned ret=" << ret
                    << ". Trying to bracket." << std::endl;
        }
        
        // Construct an initial guess for the bracket
	fp_t b_high=lg/temper, b_low=-b_high;
	fp_t yhigh=mf(b_high), ylow=mf(b_low);

        // Increase the size of the interval to ensure a valid bracket
	for(size_t j=0;j<5 && yhigh<0.0;j++) {
	  b_high*=1.0e2;
	  yhigh=mf(b_high);
	}
	for(size_t j=0;j<5 && ylow>0.0;j++) {
	  b_low*=1.0e2;
	  ylow=mf(b_low);
	}

        // If we were successful in constructing a valid bracket,
        // then call the bracketing solver
	if (yhigh>0.0 && ylow<0.0) {

          if (verbose>0) {
            std::cout << "Bracket succeeded, trying solver." << std::endl;
          }
          
	  ret=alt_solver.solve_bkt(b_low,b_high,mf);
          // If it succeeded, then set nex to the new solution
          // and set last_method
	  if (ret==0) {
            if (verbose>0) {
              std::cout << "Alternate solver succeeded." << std::endl;
            }
	    nex=b_low;
	    last_method=3000;
	  }
	}
      }

      // Restore value of err_nonconv
      density_root->err_nonconv=drec;

      if (verbose>1) {
        density_root->verbose=0;
        alt_solver.verbose=0;
      }
        
      if (ret!=0) {
        
        // Make sure we don't print out anything unless we're going
        // to call the error handler anyway
        if (this->err_nonconv==true) {
          std::cout.precision(14);
          std::cout << "Function fermion_rel::pair_density() failed.\n  "
                    << "m,ms,n,T: " << f.m << " " << f.ms << " "
                    << f.n << " " << temper << std::endl;
          std::cout << "nu: " << initial_guess << std::endl;
        }
        
        // Return the density to the user-specified value
        f.n=density_match;

	O2SCL_CONV2_RET("Density solver failed in fermion_rel::",
			"pair_density().",exc_efailed,this->err_nonconv);
      }

      // If we succeeded (i.e. if ret==0), then continue
      
      f.nu=nex*temper;
  
      if (f.non_interacting==true) { f.mu=f.nu; }

      // Finally, now that we have the chemical potential, use pair_mu()
      // to evaluate the energy density, pressure, and entropy
      int lm=last_method;
      pair_mu(f,temper);
      last_method+=lm;

      if (false && fabs(f.n-density_match)/fabs(density_match)>1.0e-5) {
        std::cout << "last_method, ret: "
                  << last_method << " " << ret << std::endl;
        std::cout << "density_root tolerances: "
                  << density_root->tol_rel << " " << density_root->tol_abs
                  << std::endl;
        std::cout << "rbg tolerances: " << alt_solver.tol_rel << " "
                  << alt_solver.tol_abs
                  << std::endl;
        std::cout << "T,n,density_match: " << temper << " " << f.n << " "
                  << density_match << std::endl;
        std::cout << "Rel. dev.: "
                  << fabs(f.n-density_match)/fabs(density_match) << std::endl;
        nex=f.nu/temper;
        std::cout << "mf: " << mf(nex) << std::endl;
        O2SCL_ERR2("Secondary failure in ",
                   "fermion_rel::pair_density().",o2scl::exc_esanity);
      }
      
      // Return the density to the user-specified value
      f.n=density_match;

      // But now that the density has been modified, we need to
      // recompute the pressure so that the thermodynamic identity is
      // satisified
      if (verify_ti==false) {
        f.pr=-f.ed+f.n*f.nu+temper*f.en;
      }

      return success;
    }
    //@}

#ifndef DOXYGEN_INTERNAL

    /// Solve for the chemical potential given the density
    fp_t solve_fun(fp_t x, fermion_t &f, fp_t T) {

      // AWS, 2/28/22: I'm getting some uninitialized variable
      // warnings, so I'm setting nden to a large value to
      // make sure that they're not causing problems.
      fp_t nden=1.0e99, yy;
      
      f.nu=T*x;

      if (f.non_interacting) f.mu=f.nu;

      bool deg=true;
      fp_t psi;
      if (f.inc_rest_mass) {
        psi=(f.nu-f.ms)/T;
      } else {
        if (f.non_interacting) {
          psi=f.nu/T;
        } else {
          psi=(f.nu+(f.m-f.ms))/T;
        }
      }
      if (psi<deg_limit) deg=false;

      // Try the non-degenerate expansion if psi is small enough
      if (use_expansions && psi<min_psi) {
        fp_t ntemp=f.n;
        bool acc=this->calc_mu_ndeg(f,T,tol_expan);
        if (acc) {
          unc.n=f.n*tol_expan;
          yy=(ntemp-f.n)/ntemp;
          f.n=ntemp;
          return yy;
        }
        f.n=ntemp;
      }

      // Try the degenerate expansion if psi is large enough
      if (use_expansions && psi>20.0) {
        fp_t ntemp=f.n;
        bool acc=this->calc_mu_deg(f,T,tol_expan);
        if (acc) {
          unc.n=f.n*tol_expan;
          yy=(ntemp-f.n)/ntemp;
          f.n=ntemp;
          return yy;
        }
        f.n=ntemp;
      }

      // Otherwise, directly perform the integration
      if (!deg) {

        fp_t y, eta;
        if (f.inc_rest_mass) {
          y=f.nu/T;
        } else {
          y=(f.nu+f.m)/T;
        }
        eta=f.ms/T;

        fp_t prefac=f.g*pow(T,3.0)/2.0/this->pi2;
        
        fri.eval_density(y,eta,nden,unc.n);
        nden*=prefac;
        unc.n*=prefac;
        
        yy=(f.n-nden)/f.n;

      } else {
    
        fp_t y=f.nu/T;
        fp_t eta=f.ms/T;
        fp_t mot;
        if (f.inc_rest_mass) {
          mot=0;
        } else {
          mot=f.m/T;
        }

        fp_t arg;
        if (f.inc_rest_mass) {
          arg=pow(upper_limit_fac*T+f.nu,2)-f.ms*f.ms;
        } else {
          arg=pow(upper_limit_fac*T+f.nu+f.m,2)-f.ms*f.ms;
        }

        fp_t ul;

        if (arg>0.0) {

          ul=sqrt(arg);
      
          fri.eval_deg_density(T,y,eta,mot,ul,nden,unc.n);
          nden*=f.g/2.0/this->pi2;
          unc.n*=f.g/2.0/this->pi2;
          
        } else {

          nden=0.0;

        }

        yy=(f.n-nden)/f.n;
      }
  
      return yy;
    }

    /** \brief Solve for the chemical potential given the density 
	with antiparticles
	
        \verbatim embed:rst

        .. todo::

        In function fermion_rel_tl::calc_density()

        - Future: Particles and antiparticles have different
        degeneracy factors, so we separately use the expansions one
        at a time. It is probably better to separately generate a
        new expansion function which automatically handles the sum
        of particles and antiparticles.

        \endverbatim
    */
    fp_t pair_fun(fp_t x, fp_t density_match, fermion_t &f, fp_t T,
                  bool log_mode) {

      // Number density of particles and antiparticles
      // AWS, 2/28/22: I'm getting some uninitialized variable
      // warnings, so I'm setting these to a large value to
      // make sure that they're not causing problems.
      fp_t nden_p=1.0e99, nden_ap=1.0e99;

      // -----------------------------------------------------------------

      f.nu=T*x;
      if (log_mode) {
	f.nu=T*exp(x);
      }

      // Sometimes the exp() call above causes an overflow, so
      // we avoid extreme values
      if (!isfinite(f.nu)) return 3;

      if (f.non_interacting) f.mu=f.nu;

      // -----------------------------------------------------------------
      // First, try the non-degenerate expansion with both particles and
      // antiparticles together

      // AWS: 6/26/21: Note that calc_mu_ndeg() includes antiparticles
      // when the last argument is true. However, this section is
      // commented out because it caused problems for the n=0, T!=0
      // case and it causes the calibrate() test function to fail.

      if (false && use_expansions) {
	if (this->calc_mu_ndeg(f,T,1.0e-18,true) && isfinite(f.n)) {
          fp_t y1;
          if (density_match==0.0) {
            y1=f.n;
          } else {
            y1=(f.n-density_match)/fabs(density_match);
          }
	  if (!isfinite(y1)) {
	    O2SCL_ERR("Value 'y1' not finite (10) in fermion_rel::pair_fun().",
		      exc_einval);
	  }
	  // Make sure to restore the value of f.n to it's original value,
	  // nn_match
	  f.n=density_match;
	  return y1;
	}
      }

      // -----------------------------------------------------------------
      // Evaluate particles and antiparticles separately. This is the
      // contribution for particles

      bool deg=true;
      fp_t psi;
      if (f.inc_rest_mass) {
	psi=(f.nu-f.ms)/T;
      } else {
        if (f.non_interacting) {
          psi=f.nu/T;
        } else {
          psi=(f.nu+(f.m-f.ms))/T;
        }
      }
      if (psi<deg_limit) deg=false;

      bool particles_done=false;

      // Try the non-degenerate expansion if psi is small enough
      if (use_expansions && psi<min_psi) {
        bool acc=this->calc_mu_ndeg(f,T,1.0e-18);
	if (acc && isfinite(f.n)) {
	  particles_done=true;
	  nden_p=f.n;
	  if (!isfinite(nden_p)) {
	    O2SCL_ERR2("Value 'nden_p' not finite (1) in ",
		       "fermion_rel::pair_fun().",exc_einval);
	  }
	}
      }
  
      // Try the degenerate expansion if psi is large enough
      if (use_expansions && particles_done==false && psi>20.0) {
	if (this->calc_mu_deg(f,T,1.0e-8) && isfinite(f.n)) {
	  particles_done=true;
	  nden_p=f.n;
	  if (!isfinite(nden_p)) {
	    O2SCL_ERR2("Value 'nden_p' not finite (2) in",
		       "fermion_rel::pair_fun().",exc_einval);
	  }
	}
      }

      // If neither expansion worked, use direct integration
      if (particles_done==false) {
    
	if (!deg) {

          fp_t y, eta;
          if (f.inc_rest_mass) {
            y=f.nu/T;
          } else {
            y=(f.nu+f.m)/T;
          }
          eta=f.ms/T;
          
	  // Nondegenerate case
      
          fp_t prefac=f.g*pow(T,3.0)/2.0/this->pi2, unc2=0;

          //bool save=fri.nit.err_nonconv;
          //fri.nit.err_nonconv=false;
          int reti1=fri.eval_density(y,eta,nden_p,unc2);
          //fri.nit.err_nonconv=save;
          if (reti1!=0) return 1;
          nden_p*=prefac;
        
	  if (!isfinite(nden_p)) {
	    O2SCL_ERR2("Value 'nden_p' not finite (3) in",
		       "fermion_rel::pair_fun().",exc_einval);
	  }
      
	} else {
      
	  // Degenerate case
      
          fp_t y=f.nu/T;
          fp_t eta=f.ms/T;
          fp_t mot;
          if (f.inc_rest_mass) {
            mot=0;
          } else {
            mot=f.m/T;
          }
        
	  fp_t arg;
	  if (f.inc_rest_mass) {
	    arg=pow(upper_limit_fac*T+f.nu,2)-f.ms*f.ms;
	  } else {
	    arg=pow(upper_limit_fac*T+f.nu+f.m,2)-f.ms*f.ms;
	  }
      
	  fp_t ul, unc2=0;
	  if (arg>0.0) {
	    ul=sqrt(arg);

            //bool save=fri.dit.err_nonconv;
            //fri.dit.err_nonconv=false;
            int reti2=fri.eval_deg_density(T,y,eta,mot,ul,nden_p,unc2);
            //fri.dit.err_nonconv=save;
            if (reti2!=0) return 2;
            nden_p*=f.g/2.0/this->pi2;
            
	  } else {
	    nden_p=0.0;
	  }
      
	  if (!isfinite(nden_p)) {
	    O2SCL_ERR2("Value 'nden_p' not finite (4) in",
		       "fermion_rel::pair_fun().",exc_einval);
	  }

	}

	particles_done=true;

	// End of 'if (particles_done==false)'
      }

      // -----------------------------------------------------------------
      // Compute the contribution from the antiparticles

      if (f.inc_rest_mass) {
	f.nu=-T*x;
	if (log_mode) f.nu=-T*exp(x);
      } else {
	f.nu=-T*x-2.0*f.m;
	if (log_mode) f.nu=-T*exp(x)-2.0*f.m;
      }
      if (f.non_interacting) f.mu=f.nu;

      bool antiparticles_done=false;

      // Evaluate the degeneracy parameter
      deg=true;
      if (f.inc_rest_mass) {
	psi=(f.nu-f.ms)/T;
      } else {
        if (f.non_interacting) {
          psi=f.nu/T;
        } else {
          psi=(f.nu+f.m-f.ms)/T;
        }
      }
      if (psi<deg_limit) deg=false;

      // Try the non-degenerate expansion if psi is small enough
      if (use_expansions && psi<min_psi) {
        bool acc=this->calc_mu_ndeg(f,T,1.0e-18);
	if (acc) {
	  antiparticles_done=true;
	  nden_ap=f.n;
	  if (!isfinite(nden_ap)) {
	    O2SCL_ERR2("Value 'nden_ap' not finite (5) in",
		       "fermion_rel::pair_fun().",
		       exc_einval);
	  }
	}
      }

      // Try the degenerate expansion if psi is large enough
      if (use_expansions && antiparticles_done==false && psi>20.0) {
	if (this->calc_mu_deg(f,T,1.0e-8)) {
	  antiparticles_done=true;
	  nden_ap=f.n;
	  if (!isfinite(nden_ap)) {
	    O2SCL_ERR2("Value 'nden_ap' not finite (6) in",
		       "fermion_rel::pair_fun().",
		       exc_einval);
	  }
	}
      }

      // If neither expansion worked, use direct integration
      if (antiparticles_done==false) {
    
	if (!deg) {
      
          fp_t y, eta;
          if (f.inc_rest_mass) {
            y=f.nu/T;
          } else {
            y=(f.nu+f.m)/T;
          }
          eta=f.ms/T;
          
	  // Nondegenerate case
          
          fp_t prefac=f.g*pow(T,3.0)/2.0/this->pi2, unc2=0;
          
          fri.eval_density(y,eta,nden_ap,unc2);
          nden_ap*=prefac;
        
	  if (!isfinite(nden_ap)) {
	    O2SCL_ERR2("Value 'nden_ap' not finite (7) in",
		       "fermion_rel::pair_fun().",
		       exc_einval);
	  }
      
	} else {
      
	  // Degenerate case
      
          fp_t y=f.nu/T;
          fp_t eta=f.ms/T;
          fp_t mot;
          if (f.inc_rest_mass) {
            mot=0;
          } else {
            mot=f.m/T;
          }
        
	  fp_t arg;
	  if (f.inc_rest_mass) {
	    arg=pow(upper_limit_fac*T+f.nu,2)-f.ms*f.ms;
	  } else {
	    arg=pow(upper_limit_fac*T+f.nu+f.m,2)-f.ms*f.ms;
	  }
      
	  fp_t ul, unc2=0;
	  if (arg>0.0) {
	    ul=sqrt(arg);
            fri.eval_deg_density(T,y,eta,mot,ul,nden_ap,unc2);
            nden_ap*=f.g/2.0/this->pi2;

	  } else {
	    nden_ap=0.0;
	  }
	  if (!isfinite(nden_ap)) {
	    O2SCL_ERR2("Value 'nden_ap' not finite (8) in",
		       "fermion_rel::pair_fun().",
		       exc_einval);
	  }

	}

	antiparticles_done=true;
      }

      fp_t y2;
      // Finish computing the function value
      if (density_match==0.0) {
	y2=fabs(nden_p-nden_ap)/fabs(nden_p);
      } else {
        y2=(nden_p-nden_ap-density_match)/fabs(density_match);
	//y2=(nden_p-nden_ap)/density_match-1.0;
      }

      if (!isfinite(y2)) {
	O2SCL_ERR("Value 'y2' not finite (9) in fermion_rel::pair_fun().",
		  exc_einval);
      }
  
      // Make sure to restore the value of f.n to it's original value,
      // density_match
      f.n=density_match;
      return y2;
    }

    
#endif

  };

  /** \brief Double precision version of \ref o2scl::fermion_rel_tl
  */
  typedef fermion_rel_tl<> fermion_rel;

  /** \brief Long double precision version of \ref o2scl::fermion_rel_tl
   */
  class fermion_rel_ld : public
  fermion_rel_tl<
    // the fermion type
    fermion_tl<long double>,
    // the Fermi-Dirac integrator
    fermi_dirac_integ_direct<long double,funct_cdf25,25,
                             cpp_dec_float_25>,
    // the Bessel-exp integrator
    bessel_K_exp_integ_boost<long double,cpp_dec_float_25>,
    // The fermion integrator. Note that we can't use
    // fermion_rel_integ here because it is based on the GSL
    // integrators which only work with double types
    fermion_rel_integ_multip<long double>,
    // The density solver
    root_brent_gsl<funct_ld,long double>,
    // The parent solver for massless fermions
    root_brent_gsl<funct_ld,long double>,
    // The function type
    funct_ld,
    // The floating-point type
    long double> {

  public:
    
    fermion_rel_ld() {

      // The goal here is to get results to within 1 part in 10^18
      
      // Tolerance for the integrator for massless fermions
      this->fd_integ.set_tol(1.0e-18);

      // No tolerance needed for the boost version
      //this->be_integ.set_tol(1.0e-18);

      // Internal function tolerances

      // This could be as large as log(1.0e4932)=11400,
      // but only 200 is used for double, so we try this for now.
      this->exp_limit=4000.0;
      
      // log(1.0e18) is 41.4
      this->upper_limit_fac=42.0;
      this->deg_entropy_fac=42.0;
      this->tol_expan=1.0e-18;

      // Solver tolerances
      this->def_density_root.tol_abs=1.0e-18;
      this->def_massless_root.tol_abs=1.0e-18;

      // Integrator tolerances
      fri.tol_rel=1.0e-18;
      fri.it.tol_rel=1.0e-18;
      fri.it2.tol_rel=1.0e-18;
      
    }
    
  };

  /** \brief 25-digit precision version of \ref o2scl::fermion_rel_tl
   */
  class fermion_rel_cdf25 : public
  fermion_rel_tl<
    // the fermion type
    fermion_tl<cpp_dec_float_25>,
    // the Fermi-Dirac integrator
    fermi_dirac_integ_direct<
      cpp_dec_float_25,funct_cdf35,25,cpp_dec_float_35>,
    // the Bessel-exp integrator
    bessel_K_exp_integ_boost<cpp_dec_float_25,cpp_dec_float_35>,
    // The fermion integrator
    fermion_rel_integ_multip<cpp_dec_float_25>,
    // The density solver
    root_brent_gsl<funct_cdf25,cpp_dec_float_25>,
    // The parent solver for massless fermions
    root_brent_gsl<funct_cdf25,cpp_dec_float_25>,
    // The function type
    funct_cdf25,
    // The floating-point type
    cpp_dec_float_25> {

  public:
    
    fermion_rel_cdf25() {

      // See output of polylog_ts for numeric limit information
      
      // Tolerance for the integrator for massless fermions
      this->fd_integ.set_tol(1.0e-23);

      // No tolerance needed for the boost version
      //this->be_integ.set_tol(1.0e-23);

      // Internal function tolerances

      this->exp_limit=1000000.0;
      
      // log(1.0e25) is 57.5
      this->upper_limit_fac=75.0;
      this->deg_entropy_fac=80.0;
      this->tol_expan=1.0e-23;

      // Solver tolerances
      this->def_density_root.tol_abs=1.0e-23;
      this->def_massless_root.tol_abs=1.0e-23;

      // Integrator tolerances
      fri.tol_rel=1.0e-23;
      fri.it.tol_rel=1.0e-23;
      fri.it2.tol_rel=1.0e-23;
    }
    
  };

#ifdef O2SCL_NEVER_DEFINED  
  /** \brief Desc
   */
  class fermion_rel_ld_multip : public
  fermion_rel_tl<
    // the fermion type
    fermion_tl<long double>,
    // the Fermi-Dirac integrator
    fermi_dirac_multip,
    // the Bessel-exp integrator
    bessel_K_exp_integ_boost<long double,cpp_dec_float_25>,
    // The fermion integrator
    fermion_rel_integ_multip<long double>,
    // The density solver
    root_brent_gsl<funct_ld,long double>,
    // The parent solver for massless fermions
    root_brent_gsl<funct_ld,long double>,
    // The function type
    funct_ld,
    // The floating-point type
    long double> {

  public:
    
    fermion_rel_ld_multip() {

      // See output of polylog_ts for numeric limit information
      
      // Tolerance for the integrator for massless fermions
      this->fd_integ.set_tol(1.0e-21);

      // Tolerance for the integrator for the nondegenerate expansion
      this->be_integ.set_tol(1.0e-21);

      // Internal function tolerances

      // This could be as large as log(1.0e4932)=11400,
      // but only 200 is used for double, so we try this for now.
      this->exp_limit=4000.0;
      
      // log(1.0e18) is 41.4
      this->upper_limit_fac=42.0;
      this->deg_entropy_fac=42.0;
      this->tol_expan=1.0e-17;

      // Solver tolerances
      this->def_density_root.tol_abs=1.0e-18;
      this->def_massless_root.tol_abs=1.0e-18;

      // Integrator tolerances
      fri.tol_rel=1.0e-16;
    }
    
  };

#endif  
#ifdef O2SCL_NEVER_DEFINED
  
  /** \brief Desc
  */
  class fermion_rel_cdf252 : public
  fermion_rel_tl<
    // the fermion type
    fermion_tl<cpp_dec_float_25>,
    // the Fermi-Dirac integrator
    fermi_dirac_integ_direct<
      cpp_dec_float_25,funct_cdf35,25,cpp_dec_float_35>,
    // the Bessel-exp integrator
    bessel_K_exp_integ_direct<
      cpp_dec_float_25,funct_cdf35,25,cpp_dec_float_35>,
    // The fermion integrator
    fermion_rel_integ_multip<cpp_dec_float_25>,
    // The density solver
    root_brent_gsl<funct_cdf25,cpp_dec_float_25>,
    // The parent solver for massless fermions
    root_brent_gsl<funct_cdf25,cpp_dec_float_25>,
    // The function type
    funct_cdf25,
    // The floating-point type
    cpp_dec_float_25> {

  public:
    
    fermion_rel_cdf252() {

      // See output of polylog_ts for numeric limit information
      
      // Tolerance for the integrator for massless fermions
      this->fd_integ.set_tol(1.0e-25);

      // Tolerance for the integrator for the nondegenerate expansion
      this->be_integ.set_tol(1.0e-25);

      // Internal function tolerances

      this->exp_limit=1000000.0;
      
      // log(1.0e25) is 57.5
      this->upper_limit_fac=58.0;
      this->deg_entropy_fac=58.0;
      this->tol_expan=1.0e-24;

      // Solver tolerances
      this->def_density_root.tol_abs=1.0e-25;
      this->def_massless_root.tol_abs=1.0e-25;

      // Integrator tolerances
      fri.tol_rel=1.0e-23;
    }
    
  };

#endif
  
}

#endif
