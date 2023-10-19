/*
  ───────────────────────────────────────────────────────────────────
  
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

  ───────────────────────────────────────────────────────────────────
*/
/** \file inte_adapt_cern.h
    \brief File defining \ref o2scl::inte_adapt_cern
*/
#ifndef O2SCL_INTE_ADAPT_CERN_H
#define O2SCL_INTE_ADAPT_CERN_H

#include <o2scl/misc.h>
#include <o2scl/funct_multip.h>
#include <o2scl/inte.h>
#include <o2scl/inte_gauss56_cern.h>
#include <o2scl/string_conv.h>
#include <o2scl/vector_derint.h>
 
namespace o2scl {

  /** \brief Integration subdivision object for \ref
      o2scl::inte_adapt_cern
  */
  template<class fp_t> class inte_subdiv {
  public:
    
    /** \brief Constructor
     */
    inte_subdiv(int n) {
      resize(n);
    }

    /// Resize the object to allow at least \c n subdivisions
    void resize(int n) {
      prev_subdiv=0;
      nsub=n;
      xlo.resize(n);
      xhi.resize(n);
      tval.resize(n);
      ters.resize(n);
      return;
    }

    /// The number of subdivisions (set in resize())
    int nsub;
    
    /// The previous number of subdivisions
    int prev_subdiv;
    
    /// Lower end of subdivision
    std::vector<fp_t> xlo;
    
    /// High end of subdivision
    std::vector<fp_t> xhi;
    
    /// Value of integral for subdivision
    std::vector<fp_t> tval;
    
    /// Squared error for subdivision
    std::vector<fp_t> ters;
    
    /// Return the number of subdivisions used in the last integration
    size_t get_nsubdivisions() {
      return prev_subdiv;
    }

    /// Return the ith subdivision
    int get_ith_subdivision(size_t i, fp_t &xlow, fp_t &xhigh, 
                            fp_t &value, fp_t &errsq) {
      if (i<prev_subdiv) {
        xlow=xlo[i];
        xhigh=xhi[i];
        value=tval[i];
        errsq=ters[i];
      }
      return 0;
    }
      
    /// Return all of the subdivisions
    template<class vec_t> 
    int get_subdivisions(vec_t &xlow, vec_t &xhigh, vec_t &value, 
                         vec_t &errsq) {

      for(int i=0;i<prev_subdiv;i++) {
        xlow[i]=xlo[i];
        xhigh[i]=xhi[i];
        value[i]=tval[i];
        errsq[i]=ters[i];
      }
      return 0;
    }
    
  };
  
  /** \brief Adaptive integration (CERNLIB)
    
      Uses a base integration object (default is \ref
      inte_gauss56_cern) to perform adaptive integration by
      automatically subdividing the integration interval. At each
      step, the interval with the largest absolute uncertainty is
      divided in half. The routine succeeds if the absolute tolerance
      is less than \ref tol_abs or if the relative tolerance is less
      than \ref tol_rel, i.e.
      \f[
      \mathrm{err}\leq\mathrm{tol\_abs}~\mathrm{or}~
      \mathrm{err}\leq\mathrm{tol\_rel}\cdot|I|
      \f]
      where \f$ I \f$ is the current estimate for the integral and \f$
      \mathrm{err} \f$ is the current estimate for the uncertainty. If
      the number of subdivisions exceeds the template parameter \c
      nsub, the error handler is called, since the integration may not
      have been successful. The number of subdivisions used in the
      last integration can be obtained from get_nsubdivisions().

      The template parameter \c nsub, is the maximum number of
      subdivisions. It is automatically set to 100 in the original
      CERNLIB routine, and defaults to 100 here. The default base
      integration object is of type \ref inte_gauss56_cern. This is the
      CERNLIB default, but can be modified by calling set_inte().

      This class is based on the CERNLIB routines RADAPT and
      DADAPT which are documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d102/top.html
      
      \future 
      - Allow user to set the initial subdivisions?
      - It might be interesting to directly compare the performance
      of this class to \ref o2scl::inte_qag_gsl .
      - There is a fixme entry in the code which could be resolved.
      - Output the point where most subdividing was required?
  */
  template<class fp_25_t=o2fp_25, class fp_35_t=o2fp_35,
           class fp_50_t=o2fp_50, class fp_100_t=o2fp_100>
  class inte_adapt_cern_tl {

  protected:

    /// \name Internal integration functions [protected]
    //@{
    /** \brief Integrate function \c func from \c a to \c b
        giving result \c res and error \c err
    */
    template<typename func_t, class fp_t>
    int integ_err_funct(func_t &func, fp_t a, fp_t b,
                        fp_t &res, fp_t &err, double target_tol,
                        double integ_tol, inte_subdiv<fp_t> &is) {

      inte_gauss56_cern<func_t,fp_t> it;
      
      fp_t tvals=0.0, terss, xlob, xhib, yhib=0.0, te, root=0.0;
      int i, nsubdivd;
      fp_t two=2.0;
      
      if (nsubdiv==0) {
        if (is.prev_subdiv==0) {
          // If the previous binning was requested, but
          // there is no previous binning stored, then
          // just shift to automatic binning
          nsubdivd=1;
        } else {
          tvals=0.0;
          terss=0.0;
          for(i=0;i<is.prev_subdiv;i++) {
            it.integ_err(func,is.xlo[i],is.xhi[i],is.tval[i],te);
            is.ters[i]=te*te;
            tvals+=is.tval[i];
            terss+=is.ters[i];
          }
          err=sqrt(two*terss);
          res=tvals;
          return 0;
        }
      }
  
      // Ensure we're not asking for too many divisions
      if (nsub<((int)nsubdiv)) {
        nsubdivd=nsub;
      } else {
        nsubdivd=nsubdiv;
      }

      // Compute the initial set of intervals and integral values
      xhib=a;
      fp_t bin=(b-a)/((fp_t)nsubdivd);
      for(i=0;i<nsubdivd;i++) {
        is.xlo[i]=xhib;
        xlob=is.xlo[i];
        is.xhi[i]=xhib+bin;
        if (i==nsubdivd-1) is.xhi[i]=b;
        xhib=is.xhi[i];
        it.integ_err(func,xlob,xhib,is.tval[i],te);
        is.ters[i]=te*te;
      }
      is.prev_subdiv=nsubdivd;

      for(size_t iter=1;iter<=((size_t)nsub);iter++) {

        // Compute the total value of the integrand
        // and the squared uncertainty
        tvals=is.tval[0];
        terss=is.ters[0];
        for(i=1;i<is.prev_subdiv;i++) {
          tvals+=is.tval[i];
          terss+=is.ters[i];
        }
          
        // Output iteration information
        if (this->verbose>0) {
          std::cout << "inte_adapt_cern_tl iter: " << iter;
          std::cout << " nsub: " << nsub;
          std::cout.setf(std::ios::showpos);
          std::cout << " Res: " << tvals;
          std::cout.unsetf(std::ios::showpos);
          std::cout << " Err: " << sqrt(two*terss);
          if (this->tol_abs>integ_tol*abs(tvals)) {
            std::cout << " Tol: " << this->tol_abs << std::endl;
          } else {
            std::cout << " Tol: " << integ_tol*abs(tvals)
                      << std::endl;
          }
          if (this->verbose>1) {
            char ch;
            std::cout << "Press a key and type enter to continue. " ;
            std::cin >> ch;
          }
        }

        // See if we're finished
        root=sqrt(two*terss);
        if (root<=this->tol_abs || root<=integ_tol*abs(tvals)) {
          res=tvals;
          err=root;
          this->last_iter=iter;
          return 0;
        }

        // Test if we've run out of intervals
        if (is.prev_subdiv==nsub) {
          res=tvals;
          err=root;
          this->last_iter=iter;
          std::string s="Reached maximum number ("+itos(nsub)+
            ") of subdivisions in inte_adapt_cern_tl::integ_err().";
          O2SCL_CONV_RET(s.c_str(),exc_etable,this->err_nonconv);
        }

        // Find the subdivision with the largest error
        fp_t bige=is.ters[0];
        int ibig=0;
        for(i=1;i<is.prev_subdiv;i++) {
          if (is.ters[i]>bige) {
            bige=is.ters[i];
            ibig=i;
          }
        }

        // Subdivide that subdivision further
        is.xhi[is.prev_subdiv]=is.xhi[ibig];
        fp_t xnew=(is.xlo[ibig]+is.xhi[ibig])/two;
        is.xhi[ibig]=xnew;
        is.xlo[is.prev_subdiv]=xnew;
        it.integ_err(func,is.xlo[ibig],is.xhi[ibig],is.tval[ibig],te);
        is.ters[ibig]=te*te;
        it.integ_err(func,is.xlo[is.prev_subdiv],
                     is.xhi[is.prev_subdiv],is.tval[is.prev_subdiv],te);
        is.ters[is.prev_subdiv]=te*te;
        is.prev_subdiv++;
        
      }

      // FIXME: Should we set an error here, or does this
      // only happen if we happen to need exactly nsub
      // intervals?
      res=tvals;
      err=root;
      return 0;
    }
    //@}

    /** \brief Desc
     */
    template <typename func_t, class fp_t>
    int integ_iu_err_int(func_t &&func, fp_t a, fp_t &res,
                         fp_t &err, double target_tol, 
                         double integ_tol, double func_tol) {
      
      inte_subdiv<fp_t> is(nsub);
      
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      funct_multip_transform<fp_t> fm2;
      fm2.lower_lim=a;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;

      std::function<fp_t(fp_t)> ft=[fm2,func](auto &&t) mutable 
      {
        return fm2.eval_iu(func,t);
      };

      fp_t one=1;
      fp_t zero=0;
      integ_err_funct(ft,zero,one,res,err,target_tol,integ_tol,is);
      
      if (verbose>1) {
        std::cout << "inte_kronrod_boost::integ_err() "
                  << "tols(target,integ,func),err:\n  "
                  << target_tol << " " << integ_tol << " "
                  << func_tol << " " << err << std::endl;
      }

#else
      err=std::numeric_limits<fp_t>::infinity();
#endif
      
      if (err/abs(res)>integ_tol) {
        std::cout << "Ret 1" << std::endl;
        return 1;
      }
      return 0;
    }
    
    /** \brief Desc
     */
    template <typename func_t, class fp_t>
    int integ_il_err_int(func_t &&func, fp_t b, fp_t &res, fp_t &err, 
                         double target_tol, double integ_tol,
                         double func_tol) {
      
      inte_subdiv<fp_t> is(nsub);
      
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      funct_multip_transform<fp_t> fm2;
      fm2.upper_lim=b;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;

      std::function<fp_t(fp_t)> ft=[fm2,func](auto &&t) mutable 
      {
        return fm2.eval_il(func,t);
      };
      
      fp_t one=1;
      fp_t zero=0;
      integ_err_funct(ft,zero,one,res,err,target_tol,integ_tol,is);
      
      if (verbose>1) {
        std::cout << "inte_kronrod_boost::integ_err() "
                  << "tols(target,integ,func),err:\n  "
                  << target_tol << " " << integ_tol << " "
                  << func_tol << " " << err << std::endl;
      }

#else
      err=std::numeric_limits<fp_t>::infinity();
#endif

      if (err/abs(res)>integ_tol) {
        return 1;
      }
      return 0;
    }

    /** \brief Desc
     */
    template <typename func_t, class fp_t>
    int integ_i_err_int(func_t &&func, fp_t &res, fp_t &err, 
                        double target_tol, double integ_tol,
                        double func_tol) {
      
      inte_subdiv<fp_t> is(nsub);
      
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      funct_multip_transform<fp_t> fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;

      std::function<fp_t(fp_t)> ft=[fm2,func](auto &&t) mutable 
      {
        return fm2.eval_i(func,t);
      };
      
      fp_t one=1;
      fp_t zero=0;
      integ_err_funct(ft,zero,one,res,err,target_tol,integ_tol,is);
      
      if (verbose>1) {
        std::cout << "inte_kronrod_boost::integ_err() "
                  << "tols(target,integ,func),err:\n  "
                  << target_tol << " " << integ_tol << " "
                  << func_tol << " " << err << std::endl;
      }

#else
      err=std::numeric_limits<fp_t>::infinity();
#endif

      if (err/abs(res)>integ_tol) {
        return 1;
      }
      return 0;
    }

    /** \brief Internal version of integration function
     */
    template <typename func_t, class fp_t>
    int integ_err_int(func_t &&func, fp_t a, fp_t b, 
                      fp_t &res, fp_t &err, double target_tol,
                      double integ_tol, double func_tol) {

      if (b==std::numeric_limits<double>::infinity()) {
        if (a==-std::numeric_limits<double>::infinity()) {
          return integ_i_err_int(func,res,err,target_tol,integ_tol,
                                 func_tol);
        } else {
          return integ_iu_err_int(func,a,res,err,target_tol,integ_tol,
                                  func_tol);
        }
      } else if (a==-std::numeric_limits<double>::infinity()) {
        return integ_il_err_int(func,b,res,err,target_tol,integ_tol,
                                func_tol);
      }
      
      inte_subdiv<fp_t> is(nsub);
      
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;

      std::function<fp_t(fp_t)> fx=[fm2,func](fp_t x) mutable -> fp_t
      { return fm2(func,x); };
      
      integ_err_funct(fx,a,b,res,err,target_tol,integ_tol,is);

      if (verbose>1) {
        std::cout << "inte_kronrod_boost::integ_err() "
                  << "tols(target,integ,func),err:\n  "
                  << target_tol << " " << integ_tol << " "
                  << func_tol << " " << err << std::endl;
      }

#else
      err=std::numeric_limits<fp_t>::infinity();
#endif
      
      if (err/abs(res)>integ_tol) {
        return 1;
      }
      return 0;
    }
    //@}

  public:

    /// \name Constructor
    //@{
    inte_adapt_cern_tl() {
      tol_rel_multip=-1.0;
      verbose=0;
      pow_tol_func=1.33;
      err_nonconv=true;
      tol_rel=1.0e-8;
      tol_abs=0.0;
      nsubdiv=1;
      nsub=100;
    }
    //@}

    /// \name Integration settings
    //@{
    /** \brief Set the number of subdivisions for the next integration
     */
    void set_nsub(int n) {
      nsub=n;
      return;
    }

    /// The number of subdivisions for the next integration
    int nsub;
    
    /** \brief The maximum relative uncertainty for multipreicsion
	integrals (default \f$ -1 \f$)
    */
    double tol_rel_multip;

    /** \brief Power for tolerance of function evaluations in
        multiprecision integrations (default 1.33)
    */
    double pow_tol_func;

    /** \brief The maximum relative uncertainty 
	in the value of the integral (default \f$ 10^{-8} \f$)
    */
    double tol_rel;

    /** \brief The maximum absolute uncertainty 
	in the value of the integral (default \f$ 10^{-8} \f$)
    */
    double tol_abs;

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief If true, call the error handler if the integration
        does not succeed (default true)
    */
    bool err_nonconv;

    /// The last number of required iterations
    int last_iter;
    //@}

    /// \name Subdivisions
    //@{
    /** \brief Number of subdivisions 

        The options are
        - 0: Use previous binning and do not subdivide further \n
        - 1: Automatic - adapt until tolerance is attained (default) \n
        - n: (n>1) split first in n equal subdivisions, then adapt
        until tolerance is obtained.
    */
    size_t nsubdiv;
    //@}

    /// \name Integration functions
    //@{
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_err(func_t &func, fp_t a, fp_t b, fp_t &res, fp_t &err) {
      
      if (b==std::numeric_limits<double>::infinity()) {
        if (a==-std::numeric_limits<double>::infinity()) {
          return integ_i_err(func,res,err);
        } else {
          return integ_iu_err(func,a,res,err);
        }
      } else if (a==-std::numeric_limits<double>::infinity()) {
        return integ_il_err(func,b,res,err);
      }
      
      inte_subdiv<fp_t> is(nsub);
      
      return integ_err_is(func,a,b,res,err,is);
    }
    
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_il_err(func_t &func, fp_t b, fp_t &res, fp_t &err) {
      
      inte_subdiv<fp_t> is(nsub);

      return integ_il_err_is(func,b,res,err,is);
    }
    
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_iu_err(func_t &func, fp_t a, fp_t &res, fp_t &err) {
      
      inte_subdiv<fp_t> is(nsub);

      return integ_il_err_is(func,a,res,err,is);
    }
    
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_i_err(func_t &func, fp_t &res, fp_t &err) {
      
      inte_subdiv<fp_t> is(nsub);

      return integ_i_err_is(func,res,err,is);
    }
    
    /** \brief Integrate function \c func from \c a to \c b.
     */
    template<typename func_t, class fp_t>
    fp_t integ(func_t &func, fp_t a, fp_t b) {
      fp_t res, interror;
      int ret=integ_err(func,a,b,res,interror);
      if (ret!=0) {
	O2SCL_ERR2("Integration failed in inte::integ(), ",
		   "but cannot return int.",o2scl::exc_efailed);
      }
      return res;
    }
    //@}

    /// \name Integration specifying subdivision
    //@{
    /** \brief Integrate function \c func from \c a to \c b and place
        the result in \c res and the error in \c err
    */
    template<typename func_t, class fp_t>
    int integ_err_is(func_t &func, fp_t a, fp_t b, fp_t &res, fp_t &err,
                     inte_subdiv<fp_t> &is) {
      
      if (b==std::numeric_limits<double>::infinity()) {
        if (a==-std::numeric_limits<double>::infinity()) {
          return integ_i_err_is(func,res,err,is);
        } else {
          return integ_iu_err_is(func,a,res,err,is);
        }
      } else if (a==-std::numeric_limits<double>::infinity()) {
        return integ_il_err_is(func,b,res,err,is);
      }
      
      if (is.nsub!=nsub) is.resize(nsub);
      
      int ret=integ_err_funct(func,a,b,res,err,
                              this->tol_rel,this->tol_rel,is);
      
      if (ret!=0) {
        if (this->verbose>0) {
          std::cout << "Function inte_kronrod_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,max: "
                    << err << " " << this->tol_rel << " "
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
    int integ_il_err_is(func_t &func, fp_t b, fp_t &res, fp_t &err,
                        inte_subdiv<fp_t> &is) {
      
      if (is.nsub!=nsub) is.resize(nsub);
      
      std::function<fp_t(fp_t)> ft=[func,b](fp_t t) mutable -> fp_t
      {
        fp_t x=b-(1-t)/t;
        fp_t y=func(x);
        return y/t/t;
      };
      
      int ret=integ_err_funct(ft,((fp_t)0),((fp_t)1),res,err,
                              this->tol_rel,this->tol_rel,is);
      
      if (ret!=0) {
        if (this->verbose>0) {
          std::cout << "Function inte_kronrod_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,max: "
                    << err << " " << this->tol_rel << " "
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
    int integ_iu_err_is(func_t &func, fp_t a, fp_t &res, fp_t &err,
                        inte_subdiv<fp_t> &is) {
      
      if (is.nsub!=nsub) is.resize(nsub);
      
      std::function<fp_t(fp_t)> ft=[func,a](fp_t t) mutable -> fp_t
      {
        fp_t x=a+(1-t)/t;
        fp_t y=func(x);
        return y/t/t;
      };
      
      int ret=integ_err_funct(ft,((fp_t)0),((fp_t)1),res,err,
                              this->tol_rel,this->tol_rel,is);
      
      if (ret!=0) {
        if (this->verbose>0) {
          std::cout << "Function inte_kronrod_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,max: "
                    << err << " " << this->tol_rel << " "
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
    int integ_i_err_is(func_t &func, fp_t &res, fp_t &err,
                       inte_subdiv<fp_t> &is) {
      
      if (is.nsub!=nsub) is.resize(nsub);
      
      std::function<fp_t(fp_t)> ft=[func](fp_t t) mutable -> fp_t
      {
        fp_t x=(1-t)/t;
        fp_t x2=-(1-t)/t;
        fp_t y=func(x)+func(x2);
        return y/t/t;
      };
      
      int ret=integ_err_funct(ft,((fp_t)0),((fp_t)1),res,err,
                              this->tol_rel,this->tol_rel,is);
      
      if (ret!=0) {
        if (this->verbose>0) {
          std::cout << "Function inte_kronrod_boost::integ_err() failed."
                    << std::endl;
          std::cout << "Values err,tol_rel,max: "
                    << err << " " << this->tol_rel << " "
                    << std::endl;
        }
        O2SCL_CONV2_RET("Failed to achieve tolerance in ",
                        "inte_kronrod_boost::integ_err().",o2scl::exc_efailed,
                        this->err_nonconv);
      }
      return 0;
    }
    //@}
    
    /// \name Multipreicison integration functions
    //@{
    /** \brief Integrate function \c func from \c a to \c b using
        multiprecision, placing the result in \c res and the error in
        \c err
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
        std::cout << "inte_adapt_cern_tl::integ_err(): set "
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<double>::digits10+3)
                    << "\n  for double integration." << std::endl;
        }
        double a_d=static_cast<double>(a);
        double b_d=static_cast<double>(b);
        double res_d, err_d;
        
        set_nsub(1000);
        ret=integ_err_int(func,a_d,b_d,res_d,err_d,
                          target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_d/abs(res_d)<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          target_tol/=10;
        }
      }
      
      if (integ_tol>pow(10.0,
                        -std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,
                           -std::numeric_limits<long double>::digits10+3)
                    << "\n  for long double integration." << std::endl;
        }
        long double a_ld=static_cast<long double>(a);
        long double b_ld=static_cast<long double>(b);
        long double res_ld, err_ld;
        
        set_nsub(1000);
        ret=integ_err_int(func,a_ld,b_ld,res_ld,err_ld,
                          target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_ld/abs(res_ld)<integ_tol) {
          res=static_cast<fp_t>(res_ld);
          err=static_cast<fp_t>(err_ld);
          return 0;
        } else {
          target_tol/=10;
        }
      }

#ifndef O2SCL_NO_BOOST_MULTIPRECISION

      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_25>::digits10+3)) {
        if (verbose>0) {
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_25>::digits10+3)
                    << "\n  for cpp_dec_float_25 integration." << std::endl;
        }
        cpp_dec_float_25 a_cdf25=static_cast<cpp_dec_float_25>(a);
        cpp_dec_float_25 b_cdf25=static_cast<cpp_dec_float_25>(b);
        cpp_dec_float_25 res_cdf25, err_cdf25;

        set_nsub(10000);
        ret=integ_err_int(func,a_cdf25,b_cdf25,res_cdf25,
                          err_cdf25,target_tol,
                          integ_tol,func_tol);

        if (verbose>1) {
          std::cout << "ret,res,err,tol: " << ret << " "
                    << res_cdf25 << " " << err_cdf25 << " "
                    << integ_tol << std::endl;
        }
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_35>::digits10+3)
                    << "\n  for cpp_dec_float_35 integration." << std::endl;
        }
        cpp_dec_float_35 a_cdf35=static_cast<cpp_dec_float_35>(a);
        cpp_dec_float_35 b_cdf35=static_cast<cpp_dec_float_35>(b);
        cpp_dec_float_35 res_cdf35, err_cdf35;
        
        set_nsub(10000);
        ret=integ_err_int(func,a_cdf35,b_cdf35,res_cdf35,
                          err_cdf35,target_tol,
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_50>::digits10+3)
                    << "\n  for cpp_dec_float_50 integration." << std::endl;
        }
        cpp_dec_float_50 a_cdf50=static_cast<cpp_dec_float_50>(a);
        cpp_dec_float_50 b_cdf50=static_cast<cpp_dec_float_50>(b);
        cpp_dec_float_50 res_cdf50, err_cdf50;
        
        set_nsub(100000);
        ret=integ_err_int(func,a_cdf50,b_cdf50,res_cdf50,
                          err_cdf50,target_tol,
                          integ_tol,func_tol);
        
        if (ret==0 && err_cdf50/abs(res_cdf50)<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_kronrod_boost::integ_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << integ_tol << std::endl;
      }

#endif      
      
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in inte_kronrod_boost::integ_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

    /** \brief Integrate function \c func from \c a to \c b using
        multiprecision, placing the result in \c res and the error in
        \c err
    */
    template <typename func_t, class fp_t>
    int integ_iu_err_multip(func_t &&func, fp_t a, 
                            fp_t &res, fp_t &err, double integ_tol=-1.0) {
      
      if (integ_tol<=0.0) {
        if (tol_rel_multip<=0.0) {
          integ_tol=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          integ_tol=tol_rel_multip;
        }
      } 
      
      if (verbose>0) {
        std::cout << "inte_adapt_cern_tl::integ_err(): set "
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<double>::digits10+3)
                    << "\n  for double integration." << std::endl;
        }
        double a_d=static_cast<double>(a);
        double res_d, err_d;

        set_nsub(1000);
        ret=integ_iu_err_int(func,a_d,res_d,err_d,
                             target_tol,integ_tol,func_tol);
        if (ret==0 && err_d/abs(res_d)<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          target_tol/=10;
        }
      }
      
      if (integ_tol>pow(10.0,
                        -std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,
                           -std::numeric_limits<long double>::digits10+3)
                    << "\n  for long double integration." << std::endl;
        }
        long double a_ld=static_cast<long double>(a);
        long double res_ld, err_ld;
        
        set_nsub(1000);
        ret=integ_iu_err_int(func,a_ld,res_ld,err_ld,
                             target_tol,integ_tol,func_tol);
        if (ret==0 && err_ld/abs(res_ld)<integ_tol) {
          res=static_cast<fp_t>(res_ld);
          err=static_cast<fp_t>(err_ld);
          return 0;
        } else {
          target_tol/=10;
        }
      }

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      
      if (integ_tol>pow(10.0,-std::numeric_limits
                        <cpp_dec_float_25>::digits10+3)) {
        if (verbose>0) {
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_25>::digits10+3)
                    << "\n  for cpp_dec_float_25 integration." << std::endl;
        }
        cpp_dec_float_25 a_cdf25=static_cast<cpp_dec_float_25>(a);
        cpp_dec_float_25 res_cdf25, err_cdf25;

        set_nsub(10000);
        ret=integ_iu_err_int(func,a_cdf25,res_cdf25,
                             err_cdf25,target_tol,
                             integ_tol,func_tol);

        if (verbose>1) {
          std::cout << "ret,res,err,tol: " << ret << " "
                    << res_cdf25 << " " << err_cdf25 << " "
                    << integ_tol << std::endl;
        }
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_35>::digits10+3)
                    << "\n  for cpp_dec_float_35 integration." << std::endl;
        }
        cpp_dec_float_35 a_cdf35=static_cast<cpp_dec_float_35>(a);
        cpp_dec_float_35 res_cdf35, err_cdf35;
        
        set_nsub(10000);
        ret=integ_iu_err_int(func,a_cdf35,res_cdf35,
                             err_cdf35,target_tol,
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_50>::digits10+3)
                    << "\n  for cpp_dec_float_50 integration." << std::endl;
        }
        cpp_dec_float_50 a_cdf50=static_cast<cpp_dec_float_50>(a);
        cpp_dec_float_50 res_cdf50, err_cdf50;
        
        set_nsub(100000);
        ret=integ_iu_err_int(func,a_cdf50,res_cdf50,
                             err_cdf50,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf50/abs(res_cdf50)<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_kronrod_boost::integ_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << integ_tol << std::endl;
      }

#endif
      
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in inte_kronrod_boost::integ_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

    /** \brief Integrate function \c func from \c a to \c b using
        multiprecision, placing the result in \c res and the error in
        \c err
    */
    template <typename func_t, class fp_t>
    int integ_il_err_multip(func_t &&func, fp_t b, 
                            fp_t &res, fp_t &err, double integ_tol=-1.0) {
      
      if (integ_tol<=0.0) {
        if (tol_rel_multip<=0.0) {
          integ_tol=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          integ_tol=tol_rel_multip;
        }
      } 
      
      if (verbose>0) {
        std::cout << "inte_adapt_cern_tl::integ_err(): set "
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<double>::digits10+3)
                    << "\n  for double integration." << std::endl;
        }
        double b_d=static_cast<double>(b);
        double res_d, err_d;
        
        ret=integ_il_err_int(func,b_d,res_d,err_d,
                             target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_d/abs(res_d)<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          target_tol/=10;
        }
      }
      
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      
      if (integ_tol>pow(10.0,
                        -std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,
                           -std::numeric_limits<long double>::digits10+3)
                    << "\n  for long double integration." << std::endl;
        }
        long double b_ld=static_cast<long double>(b);
        long double res_ld, err_ld;
        
        ret=integ_il_err_int(func,b_ld,res_ld,err_ld,
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_25>::digits10+3)
                    << "\n  for cpp_dec_float_25 integration." << std::endl;
        }
        cpp_dec_float_25 b_cdf25=static_cast<cpp_dec_float_25>(b);
        cpp_dec_float_25 res_cdf25, err_cdf25;

        ret=integ_il_err_int(func,b_cdf25,res_cdf25,
                             err_cdf25,target_tol,
                             integ_tol,func_tol);

        if (verbose>1) {
          std::cout << "ret,res,err,tol: " << ret << " "
                    << res_cdf25 << " " << err_cdf25 << " "
                    << integ_tol << std::endl;
        }
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_35>::digits10+3)
                    << "\n  for cpp_dec_float_35 integration." << std::endl;
        }
        cpp_dec_float_35 b_cdf35=static_cast<cpp_dec_float_35>(b);
        cpp_dec_float_35 res_cdf35, err_cdf35;
        
        ret=integ_il_err_int(func,b_cdf35,res_cdf35,
                             err_cdf35,target_tol,
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_50>::digits10+3)
                    << "\n  for cpp_dec_float_50 integration." << std::endl;
        }
        cpp_dec_float_50 b_cdf50=static_cast<cpp_dec_float_50>(b);
        cpp_dec_float_50 res_cdf50, err_cdf50;
        
        ret=integ_il_err_int(func,b_cdf50,res_cdf50,
                             err_cdf50,target_tol,
                             integ_tol,func_tol);
        
        if (ret==0 && err_cdf50/abs(res_cdf50)<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_kronrod_boost::integ_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << integ_tol << std::endl;
      }

#endif
      
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in inte_kronrod_boost::integ_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

    /** \brief Integrate function \c func from \c a to \c b using
        multiprecision, placing the result in \c res and the error in
        \c err
    */
    template <typename func_t, class fp_t>
    int integ_i_err_multip(func_t &&func, 
                           fp_t &res, fp_t &err, double integ_tol=-1.0) {
      
      if (integ_tol<=0.0) {
        if (tol_rel_multip<=0.0) {
          integ_tol=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        } else {
          integ_tol=tol_rel_multip;
        }
      } 
      
      if (verbose>0) {
        std::cout << "inte_adapt_cern_tl::integ_err(): set "
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits<double>::digits10+3)
                    << "\n  for double integration." << std::endl;
        }
        double res_d, err_d;
        
        ret=integ_i_err_int(func,res_d,err_d,
                            target_tol,integ_tol,func_tol);
        
        if (ret==0 && err_d/abs(res_d)<integ_tol) {
          res=static_cast<fp_t>(res_d);
          err=static_cast<fp_t>(err_d);
          return 0;
        } else {
          target_tol/=10;
        }
      }
      
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      
      if (integ_tol>pow(10.0,
                        -std::numeric_limits<long double>::digits10+3)) {
        if (verbose>0) {
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,
                           -std::numeric_limits<long double>::digits10+3)
                    << "\n  for long double integration." << std::endl;
        }
        long double res_ld, err_ld;
        
        ret=integ_i_err_int(func,res_ld,err_ld,
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_25>::digits10+3)
                    << "\n  for cpp_dec_float_25 integration." << std::endl;
        }
        cpp_dec_float_25 res_cdf25, err_cdf25;

        ret=integ_i_err_int(func,res_cdf25,
                            err_cdf25,target_tol,
                            integ_tol,func_tol);

        if (verbose>1) {
          std::cout << "ret,res,err,tol: " << ret << " "
                    << res_cdf25 << " " << err_cdf25 << " "
                    << integ_tol << std::endl;
        }
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_35>::digits10+3)
                    << "\n  for cpp_dec_float_35 integration." << std::endl;
        }
        cpp_dec_float_35 res_cdf35, err_cdf35;
        
        ret=integ_i_err_int(func,res_cdf35,
                            err_cdf35,target_tol,
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
          std::cout << "inte_adapt_cern_tl::integ_err(): "
                    << integ_tol << " > "
                    << pow(10.0,-std::numeric_limits
                           <cpp_dec_float_50>::digits10+3)
                    << "\n  for cpp_dec_float_50 integration." << std::endl;
        }
        cpp_dec_float_50 res_cdf50, err_cdf50;
        
        ret=integ_i_err_int(func,res_cdf50,
                            err_cdf50,target_tol,
                            integ_tol,func_tol);
        
        if (ret==0 && err_cdf50/abs(res_cdf50)<integ_tol) {
          res=static_cast<fp_t>(res_cdf50);
          err=static_cast<fp_t>(err_cdf50);
          return 0;
        } else {
          target_tol/=10;
        }
      }

      if (verbose>0) {
        std::cout << "inte_kronrod_boost::integ_err() "
                  << "failed after cpp_dec_float_100:\n  "
                  << integ_tol << std::endl;
      }

#endif
      
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in inte_kronrod_boost::integ_err().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }

  };

  typedef inte_adapt_cern_tl<> inte_adapt_cern;
  
}

#endif
