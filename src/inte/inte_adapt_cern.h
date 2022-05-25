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
/** \file inte_adapt_cern.h
    \brief File defining \ref o2scl::inte_adapt_cern
*/
#ifndef O2SCL_CERN_ADAPT_H
#define O2SCL_CERN_ADAPT_H

#include <o2scl/misc.h>
#include <o2scl/inte.h>
#include <o2scl/inte_gauss56_cern.h>
#include <o2scl/string_conv.h>
#include <o2scl/vector_derint.h>
 
namespace o2scl {

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
  template<class func_t=funct,
           class def_inte_t=inte_gauss56_cern<funct,double,
                                              inte_gauss56_coeffs_double>,
           size_t nsub=100, class fp_t=double>
  class inte_adapt_cern : public inte<func_t,fp_t> {

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Lower end of subdivision
    fp_t xlo[nsub];

    /// High end of subdivision
    fp_t xhi[nsub];

    /// Value of integral for subdivision
    fp_t tval[nsub];

    /// Squared error for subdivision
    fp_t ters[nsub];
      
    /// Previous number of subdivisions
    int prev_subdiv;

    /// The base integration object
    inte<func_t,fp_t> *it;
      
#endif
      
  public:
  
    inte_adapt_cern() {
      nsubdiv=1;
      prev_subdiv=0;

      it=&def_inte;
    }

    /// \name Basic usage
    //@{
    /** \brief Integrate function \c func from \c a to \c b
        giving result \c res and error \c err
    */
    virtual int integ_err(func_t &func, fp_t a, fp_t b,
                          fp_t &res, fp_t &err) {
        
      fp_t tvals=0.0, terss, xlob, xhib, yhib=0.0, te, root=0.0;
      int i, nsubdivd;
      fp_t two=2.0;

      if (nsubdiv==0) {
        if (prev_subdiv==0) {
          // If the previous binning was requested, but
          // there is no previous binning stored, then
          // just shift to automatic binning
          nsubdivd=1;
        } else {
          tvals=0.0;
          terss=0.0;
          for(i=0;i<prev_subdiv;i++) {
            it->integ_err(func,xlo[i],xhi[i],tval[i],te);
            ters[i]=te*te;
            tvals+=tval[i];
            terss+=ters[i];
          }
          err=sqrt(two*terss);
          res=tvals;
          return 0;
        }
      }
  
      // Ensure we're not asking for too many divisions
      if (nsub<nsubdiv) {
        nsubdivd=nsub;
      } else {
        nsubdivd=nsubdiv;
      }

      // Compute the initial set of intervals and integral values
      xhib=a;
      fp_t bin=(b-a)/((fp_t)nsubdivd);
      for(i=0;i<nsubdivd;i++) {
        xlo[i]=xhib;
        xlob=xlo[i];
        xhi[i]=xhib+bin;
        if (i==nsubdivd-1) xhi[i]=b;
        xhib=xhi[i];
        it->integ_err(func,xlob,xhib,tval[i],te);
        ters[i]=te*te;
      }
      prev_subdiv=nsubdivd;

      for(size_t iter=1;iter<=nsub;iter++) {

        // Compute the total value of the integrand
        // and the squared uncertainty
        tvals=tval[0];
        terss=ters[0];
        for(i=1;i<prev_subdiv;i++) {
          tvals+=tval[i];
          terss+=ters[i];
        }
          
        // Output iteration information
        if (this->verbose>0) {
          std::cout << "inte_adapt_cern Iter: " << iter;
          std::cout.setf(std::ios::showpos);
          std::cout << " Res: " << tvals;
          std::cout.unsetf(std::ios::showpos);
          std::cout << " Err: " << sqrt(two*terss);
          if (this->tol_abs>this->tol_rel*abs(tvals)) {
            std::cout << " Tol: " << this->tol_abs << std::endl;
          } else {
            std::cout << " Tol: " << this->tol_rel*abs(tvals)
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
        if (root<=this->tol_abs || root<=this->tol_rel*abs(tvals)) {
          res=tvals;
          err=root;
          this->last_iter=iter;
          return 0;
        }

        // Test if we've run out of intervals
        if (prev_subdiv==nsub) {
          res=tvals;
          err=root;
          this->last_iter=iter;
          std::string s="Reached maximum number ("+itos(nsub)+
            ") of subdivisions in inte_adapt_cern::integ_err().";
          O2SCL_CONV_RET(s.c_str(),exc_etable,this->err_nonconv);
        }

        // Find the subdivision with the largest error
        fp_t bige=ters[0];
        int ibig=0;
        for(i=1;i<prev_subdiv;i++) {
          if (ters[i]>bige) {
            bige=ters[i];
            ibig=i;
          }
        }

        // Subdivide that subdivision further
        xhi[prev_subdiv]=xhi[ibig];
        fp_t xnew=(xlo[ibig]+xhi[ibig])/two;
        xhi[ibig]=xnew;
        xlo[prev_subdiv]=xnew;
        it->integ_err(func,xlo[ibig],xhi[ibig],tval[ibig],te);
        ters[ibig]=te*te;
        it->integ_err(func,xlo[prev_subdiv],
                      xhi[prev_subdiv],tval[prev_subdiv],te);
        ters[prev_subdiv]=te*te;
        prev_subdiv++;

      }

      // FIXME: Should we set an error here, or does this
      // only happen if we happen to need exactly nsub
      // intervals?
      res=tvals;
      err=root;
      return 0;
    }
    //@}

    /// \name Integration object
    //@{
    /// Set the base integration object to use
    int set_inte(inte<func_t,fp_t> &i) {
      it=&i;
      return 0;
    }
      
    /// Default integration object
    def_inte_t def_inte;
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
    //@}
    
  };
  
  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<25>> cpp_dec_float_25;
  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<35>> cpp_dec_float_35;
  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<50>> cpp_dec_float_50;
  typedef boost::multiprecision::number<
    boost::multiprecision::cpp_dec_float<100>> cpp_dec_float_100;

  typedef
  inte_adapt_cern<funct_ld,inte_gauss56_cern
                  <funct_ld,long double,
                   inte_gauss56_coeffs_long_double>,100,
                  long double> inte_adapt_cern_ld;

  typedef
  inte_adapt_cern<funct_cdf25,inte_gauss56_cern
                  <funct_cdf25,cpp_dec_float_25,
                   inte_gauss56_coeffs_float_50<cpp_dec_float_25>>,1000,
                   cpp_dec_float_25> inte_adapt_cern_cdf25;
  
  typedef
  inte_adapt_cern<funct_cdf35,inte_gauss56_cern
                  <funct_cdf35,cpp_dec_float_35,
                   inte_gauss56_coeffs_float_50<cpp_dec_float_35>>,1000,
                   cpp_dec_float_35> inte_adapt_cern_cdf35;
  
  typedef
  inte_adapt_cern<funct_cdf50,inte_gauss56_cern
                  <funct_cdf50,cpp_dec_float_50,
                   inte_gauss56_coeffs_float_50<cpp_dec_float_50>>,1000,
                   cpp_dec_float_50> inte_adapt_cern_cdf50;
  
  typedef std::function<double(const double &)> funct_cr;
  typedef std::function<long double(const long double &)> funct_cr_ld;
  typedef std::function<cpp_dec_float_25(const cpp_dec_float_25 &)>
  funct_cr_cdf25;
  typedef std::function<cpp_dec_float_35(const cpp_dec_float_35 &)>
  funct_cr_cdf35;
  typedef std::function<cpp_dec_float_50(const cpp_dec_float_50 &)>
  funct_cr_cdf50;
  typedef std::function<cpp_dec_float_100(const cpp_dec_float_100 &)>
  funct_cr_cdf100;
  
  typedef
  inte_adapt_cern<funct_cr,inte_gauss56_cern
                  <funct_cr,double,
                   inte_gauss56_coeffs_double>,100,
                  double> inte_adapt_cern_cr;
  
  typedef
  inte_adapt_cern<funct_cr_ld,inte_gauss56_cern
                  <funct_cr_ld,long double,
                   inte_gauss56_coeffs_long_double>,100,
                  long double> inte_adapt_cern_cr_ld;

  typedef
  inte_adapt_cern<funct_cr_cdf25,inte_gauss56_cern
                  <funct_cr_cdf25,cpp_dec_float_25,
                   inte_gauss56_coeffs_float_50<cpp_dec_float_25>>,1000,
                   cpp_dec_float_25> inte_adapt_cern_cr_cdf25;
  
  typedef
  inte_adapt_cern<funct_cr_cdf35,inte_gauss56_cern
                  <funct_cr_cdf35,cpp_dec_float_35,
                   inte_gauss56_coeffs_float_50<cpp_dec_float_35>>,1000,
                   cpp_dec_float_35> inte_adapt_cern_cr_cdf35;
  
  typedef
  inte_adapt_cern<funct_cr_cdf50,inte_gauss56_cern
                  <funct_cr_cdf50,cpp_dec_float_50,
                   inte_gauss56_coeffs_float_50<cpp_dec_float_50>>,1000,
                   cpp_dec_float_50> inte_adapt_cern_cr_cdf50;
  
  
  template<class func_t=funct_multip<>>
  class inte_multip_adapt_cern {
    
  protected:
    
    /// \name The derivative objects for varying levels of precision
    //@{
    inte_adapt_cern_cr iac_d;
    inte_adapt_cern_cr_ld iac_ld;
    inte_adapt_cern_cr_cdf25 iac_cdf25;
    inte_adapt_cern_cr_cdf35 iac_cdf35;
    inte_adapt_cern_cr_cdf50 iac_cdf50;
    //@}
    
  public:

    /** \brief Relative tolerance
     */
    double tol_rel;

    /** \brief Power for tolerance of function evaluations 
        (default 1.33)
     */
    double pow_tol_func;

    /** \brief Verbosity parameter
     */
    int verbose;

    inte_multip_adapt_cern() {
      tol_rel=-1.0;
      verbose=0;
      pow_tol_func=1.33;
      iac_d.err_nonconv=false;
      iac_ld.err_nonconv=false;
      iac_cdf25.err_nonconv=false;
      iac_cdf35.err_nonconv=false;
      iac_cdf50.err_nonconv=false;
    }
    
  /** \brief Calculate the first derivative of \c func  w.r.t. x and 
      uncertainty
  */
  template<class fp_t>
  int integ_err(func_t &func, fp_t a, fp_t b, 
                fp_t &res, fp_t &err, double tol_loc=-1.0) {
    
    if (tol_loc<=0.0) {
      if (tol_rel<=0.0) {
        tol_loc=pow(10.0,-std::numeric_limits<fp_t>::digits10);
      } else {
        tol_loc=tol_rel;
      }
    } 
    
    if (verbose>0) {
      std::cout << "Function deriv_multi_gsl::deriv_err(): set "
                << "tolerance to: " << tol_loc << std::endl;
    }
    
    // Demand that the function evaluations are higher precision
    func.tol_rel=pow(tol_loc,pow_tol_func);
    
    int ret;
    
    if (tol_loc>pow(10.0,-std::numeric_limits<double>::digits10+3)) {
      double a_d=static_cast<double>(a);
      double b_d=static_cast<double>(b);
      double res_d, err_d;
      
      iac_d.tol_rel=tol_loc;
      ret=iac_d.integ_err(func,a_d,b_d,res_d,err_d);
      
      if (ret==0 && err_d<tol_loc) {
        res=static_cast<fp_t>(res_d);
        err=static_cast<fp_t>(err_d);
        return 0;
      }
    }
    
    if (tol_loc>pow(10.0,-std::numeric_limits<long double>::digits10+3)) {
      long double a_ld=static_cast<long double>(a);
      long double b_ld=static_cast<long double>(b);
      long double res_ld, err_ld;
      
      iac_ld.tol_rel=tol_loc;
      ret=iac_ld.integ_err(func,a_ld,b_ld,res_ld,err_ld);
        
      if (ret==0 && err_ld<tol_loc) {
        res=static_cast<fp_t>(res_ld);
        err=static_cast<fp_t>(err_ld);
        return 0;
      }
    }

    if (tol_loc>pow(10.0,-std::numeric_limits
                    <cpp_dec_float_25>::digits10+3)) {
      cpp_dec_float_25 a_cdf25=static_cast<cpp_dec_float_25>(a);
      cpp_dec_float_25 b_cdf25=static_cast<cpp_dec_float_25>(b);
      cpp_dec_float_25 res_cdf25, err_cdf25;
        
      iac_cdf25.tol_rel=tol_loc;
      ret=iac_cdf25.integ_err(func,a_cdf25,b_cdf25,res_cdf25,err_cdf25);
        
      if (ret==0 && err_cdf25<tol_loc) {
        res=static_cast<fp_t>(res_cdf25);
        err=static_cast<fp_t>(err_cdf25);
        return 0;
      }
    }

    if (tol_loc>pow(10.0,-std::numeric_limits
                    <cpp_dec_float_35>::digits10+3)) {
      cpp_dec_float_35 a_cdf35=static_cast<cpp_dec_float_35>(a);
      cpp_dec_float_35 b_cdf35=static_cast<cpp_dec_float_35>(b);
      cpp_dec_float_35 res_cdf35, err_cdf35;
        
      iac_cdf35.tol_rel=tol_loc;
      ret=iac_cdf35.integ_err(func,a_cdf35,b_cdf35,res_cdf35,err_cdf35);
        
      if (ret==0 && err_cdf35<tol_loc) {
        res=static_cast<fp_t>(res_cdf35);
        err=static_cast<fp_t>(err_cdf35);
        return 0;
      }
    }

    if (tol_loc>pow(10.0,-std::numeric_limits
                    <cpp_dec_float_50>::digits10+3)) {
      cpp_dec_float_50 a_cdf50=static_cast<cpp_dec_float_50>(a);
      cpp_dec_float_50 b_cdf50=static_cast<cpp_dec_float_50>(b);
      cpp_dec_float_50 res_cdf50, err_cdf50;
        
      iac_cdf50.tol_rel=tol_loc;
      ret=iac_cdf50.integ_err(func,a_cdf50,b_cdf50,res_cdf50,err_cdf50);
        
      if (ret==0 && err_cdf50<tol_loc) {
        res=static_cast<fp_t>(res_cdf50);
        err=static_cast<fp_t>(err_cdf50);
        return 0;
      }
    }

    if (verbose>0) {
      std::cout << "Function inte_multip_adapt_cern::deriv_err() "
                << "failed after cpp_dec_float_50:\n  "
                << tol_loc << std::endl;
    }
    
    O2SCL_ERR2("Failed to compute with requested accuracy ",
               "in inte_multip_adapt_cern_gsl::deriv_err().",
               o2scl::exc_efailed);
    return o2scl::exc_efailed;
  }

};
  
}

#endif
