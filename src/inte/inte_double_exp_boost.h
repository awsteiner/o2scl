/*
  -------------------------------------------------------------------
  
  Copyright (C) 2019-2023, Andrew W. Steiner
  
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
#ifndef O2SCL_INTE_DOUBLE_EXP_BOOST_H
#define O2SCL_INTE_DOUBLE_EXP_BOOST_H

/** \file inte_tanh_sinh_boost.h
    \brief File defining \ref o2scl::inte_tanh_sinh_boost
*/

#include <cmath>

#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/sinh_sinh.hpp>

#include <o2scl/inte.h>
#include <o2scl/funct_multip.h>

namespace o2scl {

  /** \brief Tanh-sinh integration class (Boost)
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .

      The native range of the integrator is -1 to 1, but supports
      infinite limits on either (or both) sides.
  */
  template<class func_t=funct, size_t max_refine=15, class fp_t=double>
  class inte_tanh_sinh_boost : public inte<func_t, fp_t> {
    
  protected:

    /// The boost integration object
    boost::math::quadrature::tanh_sinh<fp_t> it;
  
  public:

    inte_tanh_sinh_boost() : it(max_refine) {
    }
  
    virtual ~inte_tanh_sinh_boost() {
    }

    /// Return string denoting type ("inte_tanh_sinh_boost")
    virtual const char *type() { return "inte_tanh_sinh_boost"; }
    
    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			  fp_t &res, fp_t &err) {
      // Dropping the tolerance by a factor of 10 seems to help
      // the boost integrator succeed.
      res=it.integrate(func,a,b,this->tol_rel/10.0,&err,&L1norm,
                       &this->levels);
      if (err>this->tol_rel) {
        if (this->verbose>0) {
          std::cout << "Function inte_tanh_sinh_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,levels,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm << " " << this->levels << " " << max_refine
                    << std::endl;
        }
        O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_tanh_sinh_boost::integ_err().",
                        o2scl::exc_efailed,this->err_nonconv);
      }
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \f$ \infty \f$
	and place the result in \c res and the error in \c err
    */
    virtual int integ_iu_err(func_t &func, fp_t a, 
			     fp_t &res, fp_t &err) {
      return integ_err(func,a,std::numeric_limits<double>::infinity(),
		       res,err);
    }
    
    /** \brief Integrate function \c func from \f$ -\infty \f$ to \c b
	and place the result in \c res and the error in \c err
    */
    virtual int integ_il_err(func_t &func, fp_t b, 
			     fp_t &res, fp_t &err) {
      return integ_err(func,-std::numeric_limits<double>::infinity(),
		       b,res,err);
    }
    
    /** \brief Integrate function \c func from \f$ -\infty \f$ to \f$
	\infty \f$ and place the result in \c res and the error in \c
	err
    */
    virtual int integ_i_err(func_t &func, 
			    fp_t &res, fp_t &err) {
      return integ_err(func,std::numeric_limits<double>::infinity(),
		       -std::numeric_limits<double>::infinity(),res,err);
    }
  
    /** \brief Integrate function \c func from -1 to 1 and place
	the result in \c res and the error in \c err
    */
    virtual int integ_moo_err(func_t &func, fp_t &res, fp_t &err) {
      // Dropping the tolerance by a factor of 10 seems to help
      // the boost integrator succeed.
      res=it.integrate(func,this->tol_rel/10.0,&err,&L1norm,&levels);
      if (this->verbose>0) {
	std::cout << "Function inte_tanh_sinh_boost::integ_moo_err() failed."
                  << std::endl;
        std::cout << "Values err,tol_rel,L1norm,levels,max: "
		  << err << " " << this->tol_rel << " "
		  << L1norm << " " << levels << " " << max_refine
                  << std::endl;
	O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_tanh_sinh_boost::integ_moo_err().",
                        o2scl::exc_efailed,this->err_nonconv);
      }
      return 0;
    }
  
    /// L1 norm of the last integral computed
    fp_t L1norm;

    /// Number of refinement levels in last integral computed
    size_t levels;

  };
  
  /** \brief Exp-sinh integration class (Boost)
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .

      Native range is 0 to \f$ \infty \f$, but 
      any semi-infinite range is supported.
  */
  template<class func_t=funct, size_t max_refine=15, class fp_t=double>
  class inte_exp_sinh_boost : public inte<func_t, fp_t> {
    
  protected:

    /// The boost integration object
    boost::math::quadrature::exp_sinh<fp_t> it;
  
  public:

    inte_exp_sinh_boost() : it(max_refine) {
    }
  
    virtual ~inte_exp_sinh_boost() {
    }
    
    /// Return string denoting type ("inte_exp_sinh_boost")
    virtual const char *type() { return "inte_exp_sinh_boost"; }
    
    /** \brief Integrate function \c func from \c a to \c b and place
	the result in \c res and the error in \c err
    */
    virtual int integ_err(func_t &func, fp_t a, fp_t b, 
			  fp_t &res, fp_t &err) {
      // Dropping the tolerance by a factor of 10 seems to help
      // the boost integrator succeed.
      res=it.integrate(func,a,b,this->tol_rel/10.0,&err,&L1norm,&levels);
      if (err>this->tol_rel) {
        if (this->verbose>0) {
          std::cout << "Function inte_exp_sinh_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,levels,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm << " " << levels << " " << max_refine
                    << std::endl;
        }
	O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_exp_sinh_boost::integ_err().",
                        o2scl::exc_efailed,this->err_nonconv);
      }
      return 0;
    }
  
    /** \brief Integrate function \c func from \c a to \f$ \infty \f$
	and place the result in \c res and the error in \c err
    */
    virtual int integ_iu_err(func_t &func, fp_t a, 
			     fp_t &res, fp_t &err) {
      return integ_err(func,a,std::numeric_limits<double>::infinity(),
		       res,err);
    }
  
    /** \brief Integrate function \c func from \f$ -\infty \f$ to \c b
	and place the result in \c res and the error in \c err
    */
    virtual int integ_il_err(func_t &func, fp_t b, 
			     fp_t &res, fp_t &err) {
      return integ_err(func,-std::numeric_limits<double>::infinity(),
		       b,res,err);
    }
  
    /// L1 norm of the last integral computed
    fp_t L1norm;

    /// Number of refinement levels in last integral computed
    size_t levels;

  };
  
  /** \brief Sinh-sinh integration class (Boost)
      
      This class calls the error handler if the
      error returned by boost is larger than \ref inte::tol_rel .

      Only infinite limits (upper and lower) are supported.
  */
  template<class func_t=funct, size_t max_refine=15, class fp_t=double>
  class inte_sinh_sinh_boost : public inte<func_t,fp_t> {
    
  protected:

    /// The boost integration object
    boost::math::quadrature::sinh_sinh<fp_t> it;
  
  public:

    inte_sinh_sinh_boost() : it(max_refine) {
    }
  
    virtual ~inte_sinh_sinh_boost() {
    }
    
    /// Return string denoting type ("inte_sinh_sinh_boost")
    virtual const char *type() { return "inte_sinh_sinh_boost"; }
    
    /** \brief Integrate function \c func from \f$ -\infty \f$ to \f$
	\infty \f$ and place the result in \c res and the error in \c
	err
    */
    virtual int integ_i_err(func_t &func, 
			    fp_t &res, fp_t &err) {
      // Dropping the tolerance by a factor of 10 seems to help
      // the boost integrator succeed.
      res=it.integrate(func,this->tol_rel/10.0,&err,&L1norm,&levels);
      if (err>this->tol_rel) {
        if (this->verbose>0) {
          std::cout << "Function inte_sinh_sinh_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,levels,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm << " " << levels << " " << max_refine
                    << std::endl;
        }
	O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_sinh_sinh_boost::integ_err().",
                        o2scl::exc_efailed,this->err_nonconv);
      }
      return 0;
    }

    /// L1 norm of the last integral computed
    fp_t L1norm;

    /// Number of refinement levels in last integral computed
    size_t levels;

  };

  /** \brief Multiprecision integration class using Boost

      \note The uncertainties reported by this class depend on those
      returned by the boost integration object and are occasionally
      be underestimated. 
  */
  class inte_multip_double_exp_boost {
    
  protected:

    /// Desc
    size_t max_refine;
    
    template <typename func_t, class fp_t>
    int integ_err_funct(func_t &&func, fp_t a, fp_t b, 
                        fp_t &res, fp_t &err, fp_t &L1norm_loc,
                        double target_tol, double integ_tol) {
      
      boost::math::quadrature::tanh_sinh<fp_t> it_x(max_refine);
      res=it_x.integrate(func,a,b,target_tol,&err,&L1norm_loc,
                         &this->levels);

      if (verbose>1) {
        std::cout << "inte_multip_double_exp_boost::integ_err_funct() "
                  << "tols(target,integ,func),err:\n  "
                  << target_tol << " " << integ_tol << " "
                  << err << std::endl;
      }

      if (err/abs(res)>integ_tol) {
        return 1;
      }
      return 0;
    }
    
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err

        There are three tolerances:
        - \c target_tol is the target tolerance which is sent to
        the boost integration function. The error value returned 
        by the boost integration function is often larger than this
        - \c integ_tol is the desired final tolerance of the integration.
        This function regards the integration as a failure if the 
        error value is larger than \c integ_tol
        - \c func_tol is the tolerance for evaluations of the 
        integrand. This value is passed to \ref o2scl::funct_multip.

    */
    template <typename func_t, class fp_t>
    int integ_err_int(func_t &&func, fp_t a, fp_t b, 
                      fp_t &res, fp_t &err, fp_t &L1norm,
                      double target_tol, double integ_tol, double func_tol) {
      
      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;

      std::function<fp_t(fp_t)> fx=[fm2,func](fp_t x) mutable -> fp_t
      { return fm2(func,x); };

      boost::math::quadrature::tanh_sinh<fp_t> it_x(max_refine);
      res=it_x.integrate(fx,a,b,target_tol,&err,&L1norm,&this->levels);

      if (verbose>1) {
        std::cout << "inte_multip_double_exp_boost::integ_err() "
                  << "tols(target,integ,func),err:\n  "
                  << target_tol << " " << integ_tol << " "
                  << func_tol << " " << err << std::endl;
      }

      if (err/abs(res)>integ_tol) {
        return 1;
      }
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \f$ \infty \f$ 
        and place the result in \c res and the error in \c err

        There are three tolerances:
        - \c target_tol is the target tolerance which is sent to
        the boost integration function. The error value returned 
        by the boost integration function is often larger than this
        - \c integ_tol is the desired final tolerance of the integration.
        This function regards the integration as a failure if the 
        error value is larger than \c integ_tol
        - \c func_tol is the tolerance for evaluations of the 
        integrand. This value is passed to \ref o2scl::funct_multip.
     */
    template <typename func_t, class fp_t>
    int integ_iu_err_int(func_t &&func, fp_t a, 
                         fp_t &res, fp_t &err, fp_t &L1norm,
                         double target_tol, double integ_tol,
                         double func_tol) {
      
      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;
      
      std::function<fp_t(fp_t)> fx=[fm2,func](fp_t x) mutable -> fp_t
      { return fm2(func,x); };
      
      boost::math::quadrature::exp_sinh<fp_t> it_x(max_refine);
      res=it_x.integrate(fx,a,std::numeric_limits<double>::infinity(),
                       target_tol,&err,&L1norm,&this->levels);

      if (verbose>1) {
        std::cout << "inte_multip_double_exp_boost::integ_iu_err() "
                  << "tols(target,integ,func),err,L1norm:\n  "
                  << target_tol << " " << integ_tol << " "
                  << func_tol << " " << err << " "
                  << L1norm << std::endl;
      }

      if (err>integ_tol) {
        return 1;
      }
      return 0;
    }
    
    /** \brief Integrate function \c func from \c a to \f$ \infty \f$ 
        and place the result in \c res and the error in \c err

        There are three tolerances:
        - \c target_tol is the target tolerance which is sent to
        the boost integration function. The error value returned 
        by the boost integration function is often larger than this
        - \c integ_tol is the desired final tolerance of the integration.
        This function regards the integration as a failure if the 
        error value is larger than \c integ_tol
        - \c func_tol is the tolerance for evaluations of the 
        integrand. This value is passed to \ref o2scl::funct_multip.
     */
    template <typename func_t, class fp_t>
    int integ_il_err_int(func_t &&func, fp_t b, 
                         fp_t &res, fp_t &err, fp_t &L1norm,
                         double target_tol, double integ_tol,
                         double func_tol) {
      
      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;
      
      std::function<fp_t(fp_t)> fx=[fm2,func](fp_t x) mutable -> fp_t
      { return fm2(func,x); };
      
      boost::math::quadrature::exp_sinh<fp_t> it_x(max_refine);
      res=it_x.integrate(fx,-std::numeric_limits<double>::infinity(),b,
                       target_tol,&err,&L1norm,&this->levels);

      if (verbose>1) {
        std::cout << "inte_multip_double_exp_boost::integ_iu_err() "
                  << "tols(target,integ,func),err,L1norm:\n  "
                  << target_tol << " " << integ_tol << " "
                  << func_tol << " " << err << " "
                  << L1norm << std::endl;
      }

      if (err>integ_tol) {
        return 1;
      }
      return 0;
    }
    
    /** \brief Integrate function \c func from \c a to \f$ \infty \f$ 
        and place the result in \c res and the error in \c err

        There are three tolerances:
        - \c target_tol is the target tolerance which is sent to
        the boost integration function. The error value returned 
        by the boost integration function is often larger than this
        - \c integ_tol is the desired final tolerance of the integration.
        This function regards the integration as a failure if the 
        error value is larger than \c integ_tol
        - \c func_tol is the tolerance for evaluations of the 
        integrand. This value is passed to \ref o2scl::funct_multip.
     */
    template <typename func_t, class fp_t>
    int integ_i_err_int(func_t &&func, fp_t &res, fp_t &err, fp_t &L1norm,
                        double target_tol, double integ_tol,
                        double func_tol) {
      
      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;
      
      std::function<fp_t(fp_t)> fx=[fm2,func](fp_t x) mutable -> fp_t
      { return fm2(func,x); };
      
      if (verbose>1) {
        std::cout << "inte_multip_double_exp_boost::integ_i_err_int(): "
                  << std::endl;
      }
      boost::math::quadrature::sinh_sinh<fp_t> it_x(max_refine);
      res=it_x.integrate(fx,target_tol,&err,&L1norm,&this->levels);
      
      if (verbose>1) {
        std::cout << "inte_multip_double_exp_boost::integ_i_err_int() "
                  << "tols(target,integ,func),err,L1norm:\n  "
                  << target_tol << " " << integ_tol << " "
                  << func_tol << " " << err << " "
                  << L1norm << std::endl;
      }

      if (err>integ_tol) {
        return 1;
      }
      return 0;
    }
    
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<50>> cpp_dec_float_50;
    typedef boost::multiprecision::number<
      boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;

  public:

    /** \brief The maximum relative uncertainty for multipreicsion
	integrals (default \f$ -1 \f$)
    */
    double tol_rel_multip;

    /// Number of refinement levels in last integral computed
    size_t levels;
    
    /** \brief The maximum relative uncertainty 
	in the value of the integral (default \f$ 10^{-8} \f$)
    */
    double tol_rel;

    /** \brief The maximum absolute uncertainty 
	in the value of the integral (default \f$ 10^{-8} \f$)
    */
    double tol_abs;

    /** \brief Power for tolerance of function evaluations 
        (default 1.33)
    */
    double pow_tol_func;

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief If true, call the error handler if the integration
        does not succeed (default true)
    */
    bool err_nonconv;
    
    inte_multip_double_exp_boost() {
      tol_rel_multip=-1.0;
      verbose=0;
      pow_tol_func=1.33;
      err_nonconv=true;
      tol_rel=1.0e-8;
      tol_abs=1.0e-8;
      max_refine=15;
    }

    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_err(func_t &func, fp_t a, fp_t b, fp_t &res, fp_t &err) {
      
      fp_t L1norm_loc;
      int ret=integ_err_funct(func,a,b,res,err,L1norm_loc,
                              this->tol_rel/10.0,this->tol_rel);
      
      if (ret!=0) {
        if (this->verbose>0) {
          std::cout << "Function inte_kronrod_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm_loc
                    << std::endl;
        }
        O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_kronrod_boost::integ_err().",o2scl::exc_efailed,
                        this->err_nonconv);
      }
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_iu_err(func_t &func, fp_t a, fp_t &res, fp_t &err) {
      
      fp_t L1norm_loc;
      int ret=integ_iu_err_funct(func,a,res,err,L1norm_loc,
                              this->tol_rel/10.0,this->tol_rel);
      
      if (ret!=0) {
        if (this->verbose>0) {
          std::cout << "Function inte_kronrod_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm_loc
                    << std::endl;
        }
        O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_kronrod_boost::integ_err().",o2scl::exc_efailed,
                        this->err_nonconv);
      }
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_il_err(func_t &func, fp_t b, fp_t &res, fp_t &err) {
      
      fp_t L1norm_loc;
      int ret=integ_il_err_funct(func,b,res,err,L1norm_loc,
                              this->tol_rel/10.0,this->tol_rel);
      
      if (ret!=0) {
        if (this->verbose>0) {
          std::cout << "Function inte_kronrod_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm_loc
                    << std::endl;
        }
        O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_kronrod_boost::integ_err().",o2scl::exc_efailed,
                        this->err_nonconv);
      }
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_i_err(func_t &func, fp_t &res, fp_t &err) {
      
      fp_t L1norm_loc;
      int ret=integ_i_err_funct(func,res,err,L1norm_loc,
                              this->tol_rel/10.0,this->tol_rel);
      
      if (ret!=0) {
        if (this->verbose>0) {
          std::cout << "Function inte_kronrod_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,L1norm,max: "
                    << err << " " << this->tol_rel << " "
                    << L1norm_loc
                    << std::endl;
        }
        O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_kronrod_boost::integ_err().",o2scl::exc_efailed,
                        this->err_nonconv);
      }
      return 0;
    }

    /** \brief Integrate function \c func from \c a to \c b using
        multipreicsion, placing the result in \c res and the error in
        \c err

        \warning For sufficiently difficult integrands, this
        function may take a very long time to complete.
    */
    template <typename func_t, class fp_t>
    int integ_err_multip(func_t &&func, fp_t a, fp_t b, 
                         fp_t &res, fp_t &err, double integ_tol=-1.0) {
      
      if (b==std::numeric_limits<double>::infinity()) {
        if (a==-std::numeric_limits<double>::infinity()) {
          return integ_i_err_multip(func,res,err,integ_tol);
        } else {
          return integ_iu_err_multip(func,a,res,err,integ_tol);
        }
      } else if (a==-std::numeric_limits<double>::infinity()) {
        return integ_il_err_multip(func,b,res,err,integ_tol);
      }
        
      if (integ_tol<=0.0) {
        if (tol_rel_multip<=0.0) {
          integ_tol=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          integ_tol=tol_rel_multip;
        }
      } 

      if (verbose>0) {
        std::cout << "int_multip_double_exp_boost::integ_err_multip(): set "
                  << "tolerance to: " << integ_tol << std::endl;
      }
      
      // Demand that the function evaluations are higher precision
      double func_tol=pow(integ_tol,pow_tol_func);

      // We set the target tolerance an order of magnitude smaller
      // than the desired tolerance to make sure we achieve the
      // requested tolerance
      double target_tol=integ_tol/10.0;
      
      int ret;

      // We require that there are 3 more digits in the floating point
      // type than the required integration tolerance
      if (integ_tol>pow(10.0,-std::numeric_limits<double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<double>::digits10+3)
                    << "\n  for double integration." << std::endl;
        }
        double a_d=static_cast<double>(a);
        double b_d=static_cast<double>(b);
        double res_d, err_d, L1norm_d;
        
        ret=integ_err_int(func,a_d,b_d,res_d,err_d,L1norm_d,
                          target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_d/abs(res_d)<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<long double>::digits10+3)
                    << "\n  for long double integration." << std::endl;
        }
        long double a_ld=static_cast<long double>(a);
        long double b_ld=static_cast<long double>(b);
        long double res_ld, err_ld, L1norm_ld;
        
        ret=integ_err_int(func,a_ld,b_ld,res_ld,err_ld,L1norm_ld,
                          target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_ld/abs(res_ld)<integ_tol) {
          res=static_cast<fp_t>(res_ld);
          err=static_cast<fp_t>(err_ld);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_25>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_25>::digits10+3)
                    << "\n  for cpp_dec_float_25 integration." << std::endl;
        }
        cpp_dec_float_25 a_cdf25=static_cast<cpp_dec_float_25>(a);
        cpp_dec_float_25 b_cdf25=static_cast<cpp_dec_float_25>(b);
        cpp_dec_float_25 res_cdf25, err_cdf25, L1norm_cdf25;
        
        ret=integ_err_int(func,a_cdf25,b_cdf25,res_cdf25,
                          err_cdf25,L1norm_cdf25,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf25/abs(res_cdf25)<integ_tol) {
          res=static_cast<fp_t>(res_cdf25);
          err=static_cast<fp_t>(err_cdf25);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_35>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_35>::digits10+3)
                    << "\n  for cpp_dec_float_35 integration." << std::endl;
        }
        cpp_dec_float_35 a_cdf35=static_cast<cpp_dec_float_35>(a);
        cpp_dec_float_35 b_cdf35=static_cast<cpp_dec_float_35>(b);
        cpp_dec_float_35 res_cdf35, err_cdf35, L1norm_cdf35;
        
        ret=integ_err_int(func,a_cdf35,b_cdf35,res_cdf35,
                          err_cdf35,L1norm_cdf35,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf35/abs(res_cdf35)<integ_tol) {
          res=static_cast<fp_t>(res_cdf35);
          err=static_cast<fp_t>(err_cdf35);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_50>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_50>::digits10+3)
                    << "\n  for cpp_dec_float_50 integration." << std::endl;
        }
        cpp_dec_float_50 a_cdf50=static_cast<cpp_dec_float_50>(a);
        cpp_dec_float_50 b_cdf50=static_cast<cpp_dec_float_50>(b);
        cpp_dec_float_50 res_cdf50, err_cdf50, L1norm_cdf50;
        
        ret=integ_err_int(func,a_cdf50,b_cdf50,res_cdf50,
                          err_cdf50,L1norm_cdf50,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf50/abs(res_cdf50)<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_100>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_100>::digits10+3)
                    << "\n  for cpp_dec_float_100 integration." << std::endl;
        }
        cpp_dec_float_100 a_cdf100=static_cast<cpp_dec_float_100>(a);
        cpp_dec_float_100 b_cdf100=static_cast<cpp_dec_float_100>(b);
        cpp_dec_float_100 res_cdf100, err_cdf100, L1norm_cdf100;
        
        ret=integ_err_int(func,a_cdf100,b_cdf100,res_cdf100,
                          err_cdf100,L1norm_cdf100,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf100/abs(res_cdf100)<integ_tol) {
          res=static_cast<fp_t>(res_cdf100);
          err=static_cast<fp_t>(err_cdf100);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_multip_double_exp_boost::integ_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << integ_tol << std::endl;
      }
    
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in inte_multip_double_exp_boost::integ_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

    /** \brief Integrate function \c func from \c a to \f$ \infty \f$ using
        multipreicsion, placing the result in \c res and the error in
        \c err

        \warning For sufficiently difficult integrands, this
        function may take a very long time to complete.
    */
    template <typename func_t, class fp_t>
    int integ_iu_err_multip(func_t &&func, fp_t a, 
                            fp_t &res, fp_t &err, double integ_tol=-1.0) {
      
      if (integ_tol<=0.0) {
        if (tol_rel<=0.0) {
          integ_tol=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          integ_tol=tol_rel;
        }
      } 

      if (verbose>0) {
        std::cout << "int_multip_double_exp_boost::integ_iu_err(): set "
                  << "integ_tol to: " << integ_tol << std::endl;
      }
      
      // Demand that the function evaluations are higher precision
      double func_tol=pow(integ_tol,pow_tol_func);

      // We set the target tolerance an order of magnitude smaller
      // than the desired tolerance to make sure we achieve the
      // requested tolerance
      double target_tol=integ_tol/10.0;
      
      int ret;

      // We require that there are 3 more digits in the floating point
      // type than the required integration tolerance
      if (integ_tol>pow(10.0,-std::numeric_limits<double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_iu_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<double>::digits10+3)
                    << "\n  for double integration." << std::endl;
        }
        double a_d=static_cast<double>(a);
        double res_d, err_d, L1norm_d;
        
        ret=integ_iu_err_int(func,a_d,res_d,err_d,L1norm_d,
                             target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_d<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_iu_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<long double>::digits10+3)
                    << "\n  for long double integration." << std::endl;
        }
        long double a_ld=static_cast<long double>(a);
        long double res_ld, err_ld, L1norm_ld;
        
        ret=integ_iu_err_int(func,a_ld,res_ld,err_ld,L1norm_ld,
                             target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_ld<integ_tol) {
          res=static_cast<fp_t>(res_ld);
          err=static_cast<fp_t>(err_ld);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_25>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_iu_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_25>::digits10+3)
                    << "\n  for cpp_dec_float_25 integration." << std::endl;
        }
        cpp_dec_float_25 a_cdf25=static_cast<cpp_dec_float_25>(a);
        cpp_dec_float_25 res_cdf25, err_cdf25, L1norm_cdf25;
        
        ret=integ_iu_err_int(func,a_cdf25,res_cdf25,
                             err_cdf25,L1norm_cdf25,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf25<integ_tol) {
          res=static_cast<fp_t>(res_cdf25);
          err=static_cast<fp_t>(err_cdf25);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_35>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_iu_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_35>::digits10+3)
                    << "\n  for cpp_dec_float_35 integration." << std::endl;
        }
        cpp_dec_float_35 a_cdf35=static_cast<cpp_dec_float_35>(a);
        cpp_dec_float_35 res_cdf35, err_cdf35, L1norm_cdf35;
        
        ret=integ_iu_err_int(func,a_cdf35,res_cdf35,
                             err_cdf35,L1norm_cdf35,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf35<integ_tol) {
          res=static_cast<fp_t>(res_cdf35);
          err=static_cast<fp_t>(err_cdf35);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_50>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_iu_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_50>::digits10+3)
                    << "\n  for cpp_dec_float_50 integration." << std::endl;
        }
        cpp_dec_float_50 a_cdf50=static_cast<cpp_dec_float_50>(a);
        cpp_dec_float_50 res_cdf50, err_cdf50, L1norm_cdf50;
        
        ret=integ_iu_err_int(func,a_cdf50,res_cdf50,
                             err_cdf50,L1norm_cdf50,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf50<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_100>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_iu_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_100>::digits10+3)
                    << "\n  for cpp_dec_float_100 integration." << std::endl;
        }
        cpp_dec_float_100 a_cdf100=static_cast<cpp_dec_float_100>(a);
        cpp_dec_float_100 res_cdf100, err_cdf100, L1norm_cdf100;
        
        ret=integ_iu_err_int(func,a_cdf100,res_cdf100,
                             err_cdf100,L1norm_cdf100,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf100<integ_tol) {
          res=static_cast<fp_t>(res_cdf100);
          err=static_cast<fp_t>(err_cdf100);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_multip_double_exp_boost::integ_iu_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << integ_tol << std::endl;
      }
    
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in inte_multip_double_exp_boost::integ_iu_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

    /** \brief Integrate function \c func from \f$ -\infty \f$ to \c b using
        multipreicsion, placing the result in \c res and the error in
        \c err

        \warning For sufficiently difficult integrands, this
        function may take a very long time to complete.
    */
    template <typename func_t, class fp_t>
    int integ_il_err_multip(func_t &&func, fp_t b, 
                            fp_t &res, fp_t &err, double integ_tol=-1.0) {
      
      if (integ_tol<=0.0) {
        if (tol_rel<=0.0) {
          integ_tol=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          integ_tol=tol_rel;
        }
      } 

      if (verbose>0) {
        std::cout << "int_multip_double_exp_boost::integ_il_err(): set "
                  << "integ_tol to: " << integ_tol << std::endl;
      }
      
      // Demand that the function evaluations are higher precision
      double func_tol=pow(integ_tol,pow_tol_func);

      // We set the target tolerance an order of magnitude smaller
      // than the desired tolerance to make sure we achieve the
      // requested tolerance
      double target_tol=integ_tol/10.0;
      
      int ret;

      // We require that there are 3 more digits in the floating point
      // type than the required integration tolerance
      if (integ_tol>pow(10.0,-std::numeric_limits<double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_il_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<double>::digits10+3)
                    << "\n  for double integration." << std::endl;
        }
        double b_d=static_cast<double>(b);
        double res_d, err_d, L1norm_d;
        
        ret=integ_il_err_int(func,b_d,res_d,err_d,L1norm_d,
                             target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_d<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_il_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<long double>::digits10+3)
                    << "\n  for long double integration." << std::endl;
        }
        long double b_ld=static_cast<long double>(b);
        long double res_ld, err_ld, L1norm_ld;
        
        ret=integ_il_err_int(func,b_ld,res_ld,err_ld,L1norm_ld,
                             target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_ld<integ_tol) {
          res=static_cast<fp_t>(res_ld);
          err=static_cast<fp_t>(err_ld);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_25>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_il_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_25>::digits10+3)
                    << "\n  for cpp_dec_float_25 integration." << std::endl;
        }
        cpp_dec_float_25 b_cdf25=static_cast<cpp_dec_float_25>(b);
        cpp_dec_float_25 res_cdf25, err_cdf25, L1norm_cdf25;
        
        ret=integ_il_err_int(func,b_cdf25,res_cdf25,
                             err_cdf25,L1norm_cdf25,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf25<integ_tol) {
          res=static_cast<fp_t>(res_cdf25);
          err=static_cast<fp_t>(err_cdf25);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_35>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_il_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_35>::digits10+3)
                    << "\n  for cpp_dec_float_35 integration." << std::endl;
        }
        cpp_dec_float_35 b_cdf35=static_cast<cpp_dec_float_35>(b);
        cpp_dec_float_35 res_cdf35, err_cdf35, L1norm_cdf35;
        
        ret=integ_il_err_int(func,b_cdf35,res_cdf35,
                             err_cdf35,L1norm_cdf35,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf35<integ_tol) {
          res=static_cast<fp_t>(res_cdf35);
          err=static_cast<fp_t>(err_cdf35);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_50>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_il_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_50>::digits10+3)
                    << "\n  for cpp_dec_float_50 integration." << std::endl;
        }
        cpp_dec_float_50 b_cdf50=static_cast<cpp_dec_float_50>(b);
        cpp_dec_float_50 res_cdf50, err_cdf50, L1norm_cdf50;
        
        ret=integ_il_err_int(func,b_cdf50,res_cdf50,
                             err_cdf50,L1norm_cdf50,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf50<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_100>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_il_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_100>::digits10+3)
                    << "\n  for cpp_dec_float_100 integration." << std::endl;
        }
        cpp_dec_float_100 b_cdf100=static_cast<cpp_dec_float_100>(b);
        cpp_dec_float_100 res_cdf100, err_cdf100, L1norm_cdf100;
        
        ret=integ_il_err_int(func,b_cdf100,res_cdf100,
                             err_cdf100,L1norm_cdf100,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf100<integ_tol) {
          res=static_cast<fp_t>(res_cdf100);
          err=static_cast<fp_t>(err_cdf100);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_multip_double_exp_boost::integ_il_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << integ_tol << std::endl;
      }
    
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in inte_multip_double_exp_boost::integ_il_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

    /** \brief Integrate function \c func from \f$ -\infty \f$ to \f$
        \infty \f$ using multipreicsion, placing the result in \c res
        and the error in \c err

        \warning For sufficiently difficult integrands, this
        function may take a very long time to complete.
    */
    template <typename func_t, class fp_t>
    int integ_i_err_multip(func_t &&func, 
                            fp_t &res, fp_t &err, double integ_tol=-1.0) {
      
      if (integ_tol<=0.0) {
        if (tol_rel<=0.0) {
          integ_tol=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          integ_tol=tol_rel;
        }
      } 

      if (verbose>0) {
        std::cout << "int_multip_double_exp_boost::integ_i_err_multip(): set "
                  << "integ_tol to: " << integ_tol << std::endl;
      }
      
      // Demand that the function evaluations are higher precision
      double func_tol=pow(integ_tol,pow_tol_func);

      // We set the target tolerance an order of magnitude smaller
      // than the desired tolerance to make sure we achieve the
      // requested tolerance
      double target_tol=integ_tol/10.0;
      
      int ret;

      // We require that there are 3 more digits in the floating point
      // type than the required integration tolerance
      if (integ_tol>pow(10.0,-std::numeric_limits<double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_i_err_multip(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<double>::digits10+3)
                    << "\n  for double integration." << std::endl;
        }
        double res_d, err_d, L1norm_d;
        
        ret=integ_i_err_int(func,res_d,err_d,L1norm_d,
                            target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_d<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_i_err_multip(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<long double>::digits10+3)
                    << "\n  for long double integration." << std::endl;
        }
        long double res_ld, err_ld, L1norm_ld;
        
        ret=integ_i_err_int(func,res_ld,err_ld,L1norm_ld,
                             target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_ld<integ_tol) {
          res=static_cast<fp_t>(res_ld);
          err=static_cast<fp_t>(err_ld);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_25>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_i_err_multip(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_25>::digits10+3)
                    << "\n  for cpp_dec_float_25 integration." << std::endl;
        }
        cpp_dec_float_25 res_cdf25, err_cdf25, L1norm_cdf25;
        
        ret=integ_i_err_int(func,res_cdf25,
                             err_cdf25,L1norm_cdf25,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf25<integ_tol) {
          res=static_cast<fp_t>(res_cdf25);
          err=static_cast<fp_t>(err_cdf25);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_35>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_i_err_multip(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_35>::digits10+3)
                    << "\n  for cpp_dec_float_35 integration." << std::endl;
        }
        cpp_dec_float_35 res_cdf35, err_cdf35, L1norm_cdf35;
        
        ret=integ_i_err_int(func,res_cdf35,
                             err_cdf35,L1norm_cdf35,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf35<integ_tol) {
          res=static_cast<fp_t>(res_cdf35);
          err=static_cast<fp_t>(err_cdf35);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_50>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_i_err_multip(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_50>::digits10+3)
                    << "\n  for cpp_dec_float_50 integration." << std::endl;
        }
        cpp_dec_float_50 res_cdf50, err_cdf50, L1norm_cdf50;
        
        ret=integ_i_err_int(func,res_cdf50,
                             err_cdf50,L1norm_cdf50,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf50<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_100>::digits10+3)) {
        if (verbose>0) {
          std::cout << "int_multip_double_exp_boost::integ_i_err_multip(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_100>::digits10+3)
                    << "\n  for cpp_dec_float_100 integration." << std::endl;
        }
        cpp_dec_float_100 res_cdf100, err_cdf100, L1norm_cdf100;
        
        ret=integ_i_err_int(func,res_cdf100,
                             err_cdf100,L1norm_cdf100,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf100<integ_tol) {
          res=static_cast<fp_t>(res_cdf100);
          err=static_cast<fp_t>(err_cdf100);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_multip_double_exp_boost::integ_i_err_multip() "
                  << "failed after cpp_dec_float_100:\n  "
                  << integ_tol << std::endl;
      }
    
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in inte_multip_double_exp_boost::integ_i_err_multip().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

  };
  
}

#endif
