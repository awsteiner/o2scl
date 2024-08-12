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
/* roots/brent.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, 
 * Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */

#ifndef O2SCL_ROOT_BRENT_GSL_H
#define O2SCL_ROOT_BRENT_GSL_H

/** \file root_brent_gsl.h
    \brief File defining \ref o2scl::root_brent_gsl 
*/

#include <limits>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <o2scl/funct.h>
#include <o2scl/funct_multip.h>
#include <o2scl/root.h>

namespace o2scl {
  
  /** \brief One-dimensional root-finding (GSL)

      This class finds the root of a user-specified function. If \ref
      test_form is 0 (the default), then solve_bkt() stops when the
      size of the bracket is smaller than \ref root::tol_abs. If \ref
      test_form is 1, then the function stops when the residual is
      less than \ref root::tol_rel. If test_form is 2, then both tests
      are applied.

      \verbatim embed:rst
      See the :ref:`One-dimensional solvers` section of the User's
      guide for general information about O2scl solvers. An example
      demonstrating the usage of this class is given in
      ``examples/ex_fptr.cpp`` and the :ref:`First function object
      example`.

      .. todo::

         class root_brent_gsl

         Future:

         - There is some duplication in the variables \c x_lower, 
         \c x_upper, \c a, and \c b, which could be removed. Some
         better variable names would also be helpful.
         - Create a meaningful enum list for \ref
         o2scl::root_brent_gsl::test_form.
         - There is code duplication between the test_interval here
         and in root_toms748.

      \endverbatim

      \comment
      Note that I used \c instead of \ref to refer to variables above
      since the variables a and b are protected, and not available if
      compiling the documentation without the internal portion.

      Also, at first I got confused that GSL didn't require
      lower<upper, but it turns out that this is indeed a requirement
      in GSL, but I didn't see it because it was in roots/fsolver.c
      rather than in roots/brent.c . Thus, everything looks fine now.
      \endcomment
  */
  template<class func_t=funct, class fp_t=double> class root_brent_gsl : 
    public root_bkt<func_t,func_t,fp_t> {

  protected:

    /** \brief Floating point-type agnostic version of
	\c gsl_root_test_interval() .
    */
    template<class fp2_t>
    int test_interval(fp2_t xx_lower, fp2_t xx_upper, fp2_t epsabs,
		      fp2_t epsrel, fp2_t &tolerance, fp2_t &interval) {
      
      fp2_t abs_lower, abs_upper;

      if (xx_lower<0.0) abs_lower=-xx_lower;
      else abs_lower=xx_lower;
      if (xx_upper<0.0) abs_upper=-xx_upper;
      else abs_upper=xx_upper;
      
      fp2_t min_abs;
      if (epsrel<0.0) {
	O2SCL_ERR2("Relative tolerance is negative in ",
		   "root_brent_gsl::test_interval().",o2scl::exc_ebadtol);
      }
      if (epsabs<0.0) {
	O2SCL_ERR2("Absolute tolerance is negative in ",
		   "root_brent_gsl::test_interval().",o2scl::exc_ebadtol);
      }
      if (xx_lower>xx_upper) {
	O2SCL_ERR2("Lower bound larger than upper bound in ",
		   "root_brent_gsl::test_interval().",o2scl::exc_einval);
      }

      if ((xx_lower>0.0 && xx_upper>0.0) ||
	  (xx_lower<0.0 && xx_upper<0.0)) {
	if (abs_lower<abs_upper) min_abs=abs_lower;
	else min_abs=abs_upper;
      } else {
	min_abs=0;
      }

      tolerance=epsabs+epsrel*min_abs;
      
      if (xx_lower<xx_upper) {
        interval=xx_upper-xx_lower;
      } else {
        interval=xx_lower-xx_upper;
      }
      if (interval<tolerance) {
        return o2scl::success;
      }
      
      return o2scl::gsl_continue;
    }

    /** \brief Solve \c func in region \f$ x_1<x<x_2 \f$
        returning \f$ x_1 \f$ (internal multiprecision version)
     */
    template<typename func2_t, class fp2_t>
    int solve_bkt_int_multip(fp2_t &x1, fp2_t x2, func2_t &f,
                             double root_tol, double func_tol) {
                             
      fp2_t storage[11];
      
      fp2_t &lroot=storage[0];
      fp2_t &lx_lower=storage[1];
      fp2_t &lx_upper=storage[2];
      fp2_t &la=storage[3];
      fp2_t &lb=storage[4];
      fp2_t &lc=storage[5];
      fp2_t &ld=storage[6];
      fp2_t &le=storage[7];
      fp2_t &lfa=storage[8];
      fp2_t &lfb=storage[9];
      fp2_t &lfc=storage[10];
      
      int status, iter=0;

      int set_ret=set_multip(f,x1,x2,func_tol,storage);
      if (set_ret!=0) return set_ret;
    
      if (test_form==0) {

        // Test the bracket size

        status=gsl_continue;
        while (status==gsl_continue && iter<this->ntrial) {
      
          iter++;
          iterate_multip(f,func_tol,storage);
          fp2_t tol, interval;
          fp2_t ta=static_cast<fp2_t>(0);
          fp2_t tr=static_cast<fp2_t>(root_tol);
          status=test_interval(lx_lower,lx_upper,ta,tr,tol,interval);
      
          if (this->verbose>0) {
            fp2_t y;

            y=f(lroot);

            // This additional temporary seems to be required
            // for boost::multiprecision types
            fp2_t x_diff=lx_upper-lx_lower;
            //this->print_iter(lroot,y,iter,interval,tol,
            //"root_brent_gsl (test_interval)");
          }
        }
	
      } else if (test_form==1) {

        // Test the residual

        status=gsl_continue;
        while (status==gsl_continue && iter<this->ntrial) {
      
          iter++;
          iterate_multip(f,func_tol,storage);

          fp2_t y=f(lroot);

          if (abs(y)<this->tol_rel) status=o2scl::success;
      
          if (this->verbose>0) {
            //this->print_iter(lroot,y,iter,abs(y),this->tol_rel,
            //"root_brent_gsl (relative deviation)");
          }
        }


      } else {

        // Test the bracket size and the residual

        status=gsl_continue;
        while (status==gsl_continue && iter<this->ntrial) {
      
          iter++;
          iterate_multip(f,func_tol,storage);
        
          fp2_t tol, interval;
          fp2_t ta=static_cast<fp2_t>(0);
          fp2_t tr=static_cast<fp2_t>(root_tol);
          status=test_interval(lx_lower,lx_upper,ta,tr,tol,interval);
        
          if (status==o2scl::success) {
            fp2_t y=f(lroot);
            if (abs(y)>=this->tol_rel) status=gsl_continue;
            if (this->verbose>0) {
              //this->print_iter(root,y,iter,abs(y),this->tol_rel,
              //"root_brent_gsl (relative deviation 2)");
            }
          } else {
            if (this->verbose>0) {
              fp2_t y=f(lroot);
              // This additional temporary seems to be required
              // for boost::multiprecision types
              fp2_t x_diff=lx_upper-lx_lower;
              std::cout << "lower,root,upper: "
                        << lx_lower << " " << lroot << " "
                        << lx_upper << std::endl;
              //this->print_iter(root,y,iter,interval,tol,
              //"root_brent_gsl (test_interval 2)");
            }
          }
        }

      }

      x1=lroot;
  
      if (iter>=this->ntrial) {
        O2SCL_CONV2_RET("Function root_brent_gsl::solve_bkt() exceeded ",
                        "maximum number of iterations.",o2scl::exc_emaxiter,
                        this->err_nonconv);
      }
  
      if (status!=o2scl::success) {
        return status;
      }

      return o2scl::success;
    }

  public:

    root_brent_gsl() {
      test_form=0;
      pow_tol_func=1.33;
    }
    
    /// Return the type, \c "root_brent_gsl".
    virtual const char *type() { return "root_brent_gsl"; }

    /** \brief Power for tolerance of function evaluations
        with adaptive multiprecision (default 1.33)
    */
    double pow_tol_func;
    
    /** \brief Perform an iteration

        This function currently always returns \ref success.
    */
    int iterate(func_t &f) {
      
      fp_t tol, m, two=2;
	
      int ac_equal=0;
	
      if ((gfb < 0 && gfc < 0) || (gfb > 0 && gfc > 0)) {
        ac_equal=1;
        gc=ga;
        gfc=gfa;
        gd=gb-ga;
        ge=gb-ga;
      }
  
      if (abs(gfc) < abs(gfb)) {
        ac_equal=1;
        ga=gb;
        gb=gc;
        gc=ga;
        gfa=gfb;
        gfb=gfc;
        gfc=gfa;
      }
    
      tol=abs(gb)*std::numeric_limits<fp_t>::epsilon()/two;
      m=(gc-gb)/two;
  
      if (gfb == 0) {
        groot=gb;
        gx_lower=gb;
        gx_upper=gb;
    
        return o2scl::success;
      }
      if (abs(m) <= tol) {
        groot=gb;
    
        if (gb < gc) {
          gx_lower=gb;
          gx_upper=gc;
        } else {
          gx_lower=gc;
          gx_upper=gb;
        }
    
        return o2scl::success;
      }
  
      if (abs(ge) < tol || abs(gfa) <= abs(gfb)) {
        // [GSL] Use bisection 
        gd=m;            
        ge=m;
      } else {

        // [GSL] Use inverse cubic interpolation 
        fp_t p, q, r;   
        fp_t s=gfb/gfa;
    
        if (ac_equal) {
          p=2*m*s;
          q=1-s;
        } else {
          q=gfa/gfc;
          r=gfb/gfc;
          p=s*(2*m*q*(q-r)-(gb-ga)*(r-1));
          q=(q-1)*(r-1)*(s-1);
        }
      
        if (p > 0) {
          q=-q;
        } else {
          p=-p;
        }
        fp_t dtmp;
        fp_t ptmp=ge*q;
        fp_t ptmp2=tol*q;
        if (3*m*q-abs(ptmp2)<abs(ptmp)) {
          dtmp=3*m*q-abs(ptmp2);
        } else {
          dtmp=abs(ptmp);
        }
        if (2*p<dtmp) {
          ge=gd;
          gd=p/q;
        } else {
          // [GSL] Interpolation failed, fall back to bisection.
          gd=m;
          ge=m;
        }
      }
  
      ga=gb;
      gfa=gfb;
  
      if (abs(gd) > tol) {
        gb+=gd;
      } else {
        gb+=(m > 0 ? +tol : -tol);
      }
	
      gfb=f(gb);
  
      // Update the best estimate of the root and bounds on each
      // iteration
	
      groot=gb;
  
      if ((gfb < 0 && gfc < 0) || (gfb > 0 && gfc > 0)) {
        gc=ga;
      }
      if (gb < gc) {
        gx_lower=gb;
        gx_upper=gc;
      } else {
        gx_lower=gc;
        gx_upper=gb;
      }
  
      return o2scl::success;
    }

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
    
    /** \brief Perform an iteration (adaptive multiprecision version)
    */
    template<typename func2_t, class fp2_t>
    int iterate_multip(func2_t &f, double func_tol, fp2_t storage[11]) {
      
      fp2_t &lroot=storage[0];
      fp2_t &lx_lower=storage[1];
      fp2_t &lx_upper=storage[2];
      fp2_t &la=storage[3];
      fp2_t &lb=storage[4];
      fp2_t &lc=storage[5];
      fp2_t &ld=storage[6];
      fp2_t &le=storage[7];
      fp2_t &lfa=storage[8];
      fp2_t &lfb=storage[9];
      fp2_t &lfc=storage[10];

      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;
      
      fp2_t tol, m, two=2;
	
      int ac_equal=0;
	
      if ((lfb < 0 && lfc < 0) || (lfb > 0 && lfc > 0)) {
        ac_equal=1;
        lc=la;
        lfc=lfa;
        ld=lb-la;
        le=lb-la;
      }
  
      if (abs(lfc) < abs(lfb)) {
        ac_equal=1;
        la=lb;
        lb=lc;
        lc=la;
        lfa=lfb;
        lfb=lfc;
        lfc=lfa;
      }
    
      tol=abs(lb)*std::numeric_limits<fp2_t>::epsilon()/two;
      m=(lc-lb)/two;
  
      if (lfb == 0) {
        lroot=lb;
        lx_lower=lb;
        lx_upper=lb;
    
        return o2scl::success;
      }
      if (abs(m) <= tol) {
        lroot=lb;
    
        if (lb < lc) {
          lx_lower=lb;
          lx_upper=lc;
        } else {
          lx_lower=lc;
          lx_upper=lb;
        }
    
        return o2scl::success;
      }
  
      if (abs(le) < tol || abs(lfa) <= abs(lfb)) {
        // [GSL] Use bisection 
        ld=m;            
        le=m;
      } else {

        // [GSL] Use inverse cubic interpolation 
        fp2_t p, q, r;   
        fp2_t s=lfb/lfa;
    
        if (ac_equal) {
          p=2*m*s;
          q=1-s;
        } else {
          q=lfa/lfc;
          r=lfb/lfc;
          p=s*(2*m*q*(q-r)-(lb-la)*(r-1));
          q=(q-1)*(r-1)*(s-1);
        }
      
        if (p > 0) {
          q=-q;
        } else {
          p=-p;
        }
        fp2_t dtmp;
        fp2_t ptmp=le*q;
        fp2_t ptmp2=tol*q;
        if (3*m*q-abs(ptmp2)<abs(ptmp)) {
          dtmp=3*m*q-abs(ptmp2);
        } else {
          dtmp=abs(ptmp);
        }
        if (2*p<dtmp) {
          le=ld;
          ld=p/q;
        } else {
          // [GSL] Interpolation failed, fall back to bisection.
          ld=m;
          le=m;
        }
      }
  
      la=lb;
      lfa=lfb;
  
      if (abs(ld) > tol) {
        lb+=ld;
      } else {
        lb+=(m > 0 ? +tol : -tol);
      }

      fp2_t err;
      int fm_ret1=fm2.eval_tol_err(f,lb,lfb,err);
      if (fm_ret1!=0) return 1;
  
      // Update the best estimate of the root and bounds on each
      // iteration
	
      lroot=lb;
  
      if ((lfb < 0 && lfc < 0) || (lfb > 0 && lfc > 0)) {
        lc=la;
      }
      if (lb < lc) {
        lx_lower=lb;
        lx_upper=lc;
      } else {
        lx_lower=lc;
        lx_upper=lb;
      }
  
      return o2scl::success;
    }

#endif
      
    /** \brief Solve \c func in region \f$ x_1<x<x_2 \f$ returning
        \f$ x_1 \f$.
     */
    virtual int solve_bkt(fp_t &x1, fp_t x2, func_t &f) {
	
      int status, iter=0;

      int set_ret=set(f,x1,x2);
      if (set_ret!=0) return set_ret;
    
      if (test_form==0) {

        // Test the bracket size

        status=gsl_continue;
        while (status==gsl_continue && iter<this->ntrial) {
      
          iter++;
          iterate(f);
          fp_t tol, interval;
          status=test_interval(gx_lower,gx_upper,
                               this->tol_abs,this->tol_rel,tol,interval);
      
          if (this->verbose>0) {
            fp_t y;

            y=f(groot);

            // This additional temporary seems to be required
            // for boost::multiprecision types
            fp_t x_diff=gx_upper-gx_lower;
            this->print_iter(groot,y,iter,interval,tol,
                             "root_brent_gsl (test_interval)");
          }
        }
	
      } else if (test_form==1) {

        // Test the residual

        status=gsl_continue;
        while (status==gsl_continue && iter<this->ntrial) {
      
          iter++;
          iterate(f);

          fp_t y=f(groot);

          if (abs(y)<this->tol_rel) status=o2scl::success;
      
          if (this->verbose>0) {
            this->print_iter(groot,y,iter,abs(y),this->tol_rel,
                             "root_brent_gsl (relative deviation)");
          }
        }


      } else {

        // Test the bracket size and the residual

        status=gsl_continue;
        while (status==gsl_continue && iter<this->ntrial) {
      
          iter++;
          iterate(f);
        
          fp_t tol, interval;
          status=test_interval(gx_lower,gx_upper,
                               this->tol_abs,this->tol_rel,tol,interval);
        
          if (status==o2scl::success) {
            fp_t y=f(groot);
            if (abs(y)>=this->tol_rel) status=gsl_continue;
            if (this->verbose>0) {
              this->print_iter(groot,y,iter,abs(y),this->tol_rel,
                               "root_brent_gsl (relative deviation 2)");
            }
          } else {
            if (this->verbose>0) {
              fp_t y=f(groot);
              // This additional temporary seems to be required
              // for boost::multiprecision types
              fp_t x_diff=gx_upper-gx_lower;
              std::cout << "lower,root,upper: "
                        << gx_lower << " " << groot << " "
                        << gx_upper << std::endl;
              this->print_iter(groot,y,iter,interval,tol,
                               "root_brent_gsl (test_interval 2)");
            }
          }
        }

      }

      x1=groot;
  
      if (iter>=this->ntrial) {
        O2SCL_CONV2_RET("Function root_brent_gsl::solve_bkt() exceeded ",
                        "maximum number of iterations.",o2scl::exc_emaxiter,
                        this->err_nonconv);
      }
  
      if (status!=o2scl::success) {
        return status;
      }

      return o2scl::success;
    }

    /** \brief Solve \c func in region \f$ x_1<x<x_2 \f$ returning
        \f$ x_1 \f$ (adaptive multiprecision version)
     */
    template<typename func2_t, class fp2_t>
    int solve_bkt_multip(fp2_t &x1, fp2_t x2, func2_t &&f, fp2_t &err,
                         double tol_loc=-1.0) {
      
      if (tol_loc<=0.0) {
        tol_loc=pow(10.0,-std::numeric_limits<fp2_t>::digits10);
      } 
      
      if (this->verbose>0) {
        std::cout << "solve_bkt_multip(): set "
                  << "tolerance to: " << tol_loc << std::endl;
      }
      
      // Demand that the function evaluations are higher precision
      double func_tol=pow(tol_loc,pow_tol_func);

      // ─────────────────────────────────────────────────────────────────
      // Double precision derivative evaluation
      
      bool called_d=false;
      double x1_d;

      // Begin with a double precision derivative, but only
      // if double precision is enough for the function evaluations
      if (tol_loc>pow(10.0,-std::numeric_limits<double>::digits10)) {
        
        x1_d=static_cast<double>(x1);
        solve_bkt_int_multip(x1_d,static_cast<double>(x2),f,
                             tol_loc,func_tol);
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "double." << std::endl;
        }
        
        called_d=true;

        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "failed after double:\n  "
                    << dtos(x1_d,0) << " "
                    << tol_loc << std::endl;
        }
        
      } else if (this->verbose>0) {
        std::cout << "double: " << tol_loc << " "
                  << pow(10.0,-std::numeric_limits<double>::digits10)
                  << std::endl;
      }

      // ─────────────────────────────────────────────────────────────────
      // Long double precision derivative evaluation
      
      bool called_ld=false;
      long double x1_ld, err_ld=0;
      
      // Attempt to evaluate at long double precision, but only if 
      // long double precision is enough for the function evaluations
      if (tol_loc>pow(10.0,-std::numeric_limits<long double>::digits10)) {
        
        x1_ld=static_cast<double>(x1);
        solve_bkt_int_multip(x1_ld,static_cast<long double>(x2),f,
                             tol_loc,func_tol);
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "long double." << std::endl;
        }
        
        called_ld=true;
        
        // If the comparison between the double and long double
        // precision results shows an accurate result, then return
        if (called_d) {
          err=static_cast<fp2_t>(abs(x1_ld-x1_d)/abs(x1_ld));
          if (err<tol_loc) {
            x1=static_cast<fp2_t>(x1_ld);
            return 0;
          }
        }
        
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "after long double:\n  ";
          if (called_d) {
            std::cout << dtos(x1_d,0) << " ";
          }
          std::cout << dtos(x1_ld,0) << " "
                    << tol_loc << std::endl;
        }

      } else if (this->verbose>0) {
        std::cout << "long double: " << tol_loc << " "
                  << pow(10.0,-std::numeric_limits<double>::digits10)
                  << std::endl;
      }
      
#ifndef O2SCL_NO_BOOST_MULTIPRECISION
      
      // ─────────────────────────────────────────────────────────────────
      // 25-digit precision derivative evaluation
      
      bool called_cdf25=false;
      cpp_dec_float_25 x1_cdf25;
      
      // Attempt to evaluate at 25-digit precision, but only if
      // 25-digit precision is enough for the function evaluations
      if (tol_loc>pow(10.0,-std::numeric_limits
                       <cpp_dec_float_25>::digits10)) {
        
        x1_cdf25=static_cast<double>(x1);
        solve_bkt_int_multip(x1_cdf25,
                             static_cast<cpp_dec_float_25>(x2),f,
                             tol_loc,func_tol);
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "[cpp_dec_float_25]\n "
                    << "x1: "
                    << o2scl::dtos(x1_cdf25,0) << std::endl;
        }
        
        called_cdf25=true;
        
        // If the comparison between the long double and 25-digit
        // precision results shows an accurate result, then return
        if (called_ld) {
          err=static_cast<fp2_t>(abs(x1_cdf25-x1_ld)/abs(x1_cdf25));
          if (this->verbose>0) {
            std::cout << "root_brent_gsl::solve_bkt_multip() "
                      << "[cpp_dec_float_25]\n "
                      << "err,tol_loc: " << err << " " << tol_loc
                      << std::endl;
          }
          if (err<tol_loc) {
            x1=static_cast<fp2_t>(x1_cdf25);
            return 0;
          }
        }
        
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "after cpp_dec_25:\n  "
                    << dtos(x1_ld,0) << " "
                    << dtos(x1_cdf25,0) << " "
                    << dtos(err,0) << " " 
                    << tol_loc << std::endl;
        }
        
      }
    
      // ─────────────────────────────────────────────────────────────────
      // 35-digit precision derivative evaluation

      bool called_cdf35=false;
      cpp_dec_float_35 x1_cdf35;

      // Attempt to evaluate at 35-digit precision, but only if
      // 35-digit precision is enough for the function evaluations
      if (tol_loc>pow(10.0,-std::numeric_limits
                       <cpp_dec_float_35>::digits10)) {
        
        x1_cdf35=static_cast<double>(x1);
        solve_bkt_int_multip(x1_cdf35,
                             static_cast<cpp_dec_float_35>(x2),f,
                             tol_loc,func_tol);
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "[cpp_dec_float_35]\n "
                    << "x1,err: "
                    << x1_cdf35 << std::endl;
        }
        
        called_cdf35=true;
        
        // If the comparison between the 25-digit and 35-digit
        // precision results shows an accurate result, then return
        if (called_cdf25 && x1_cdf35!=0) {
          err=static_cast<fp2_t>(abs(x1_cdf35-x1_cdf25)/abs(x1_cdf35));
          if (err<tol_loc) {
            x1=static_cast<fp2_t>(x1_cdf35);
            return 0;
          }
        }
        
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "failed after cpp_dec_float_35:\n  "
                    << dtos(x1_cdf25,0) << " "
                    << dtos(x1_cdf35,0) << " "
                    << dtos(err,0) << " " 
                    << tol_loc << std::endl;
        }

      }
    
      // ─────────────────────────────────────────────────────────────────
      // 50-digit precision derivative evaluation
      
      bool called_cdf50=false;
      cpp_dec_float_50 x1_cdf50;
      
      // Attempt to evaluate at 50-digit precision, but only if
      // 50-digit precision is enough for the function evaluations
      if (tol_loc>pow(10.0,-std::numeric_limits
                       <cpp_dec_float_50>::digits10)) {
        
        x1_cdf50=static_cast<double>(x1);
        solve_bkt_int_multip(x1_cdf50,
                             static_cast<cpp_dec_float_50>(x2),f,
                             tol_loc,func_tol);
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "[cpp_dec_float_50]\n "
                    << "x1,err: "
                    << x1_cdf50 << std::endl;
        }
        
        called_cdf50=true;
        
        // If the comparison between the 35-digit and 50-digit
        // precision results shows an accurate result, then return
        if (called_cdf35 && x1_cdf50!=0) {
          err=static_cast<fp2_t>(abs(x1_cdf50-x1_cdf35)/abs(x1_cdf50));
          if (err<tol_loc) {
            x1=static_cast<fp2_t>(x1_cdf50);
            return 0;
          }
        }
      
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "failed after cpp_dec_float_50:\n  "
                    << dtos(x1_cdf35,0) << " "
                    << dtos(x1_cdf50,0) << " "
                    << dtos(err,0) << " " 
                    << tol_loc << std::endl;
        }

      }
      
      // ─────────────────────────────────────────────────────────────────
      // 100-digit precision derivative evaluation
      
      bool called_cdf100=false;
      cpp_dec_float_100 x1_cdf100;
      
      // Attempt to evaluate at 100-digit precision, but only if
      // 100-digit precision is enough for the function evaluations
      if (tol_loc>pow(10.0,-std::numeric_limits
                       <cpp_dec_float_100>::digits10)) {
        
        x1_cdf100=static_cast<double>(x1);
        solve_bkt_int_multip(x1_cdf100,
                             static_cast<cpp_dec_float_100>(x2),f,
                             tol_loc,func_tol);
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "[cpp_dec_float_100]\n "
                    << "x1,err: "
                    << x1_cdf100 << std::endl;
        }
        
        called_cdf100=true;
        
        // If the comparison between the 50-digit and 100-digit
        // precision results shows an accurate result, then return
        if (called_cdf50 && x1_cdf100!=0) {
          err=static_cast<fp2_t>(abs(x1_cdf100-x1_cdf50)/
                                 abs(x1_cdf100));
          if (err<tol_loc) {
            x1=static_cast<fp2_t>(x1_cdf100);
            return 0;
          }
        }
      
        if (this->verbose>0) {
          std::cout << "root_brent_gsl::solve_bkt_multip() "
                    << "failed after cpp_dec_float_100:\n  "
                    << dtos(x1_cdf50,0) << " "
                    << dtos(x1_cdf100,0) << " "
                    << dtos(err,0) << " " 
                    << tol_loc << std::endl;
        }

      }

#endif
      
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in root_brent_gsl::solve_bkt_multip().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }
    
    /// The type of convergence test applied: 0, 1, or 2 (default 0)
    int test_form;
     
    /// Get the most recent value of the root
    fp_t get_root() { return groot; }
      
    /// Get the lower limit
    fp_t get_lower() { return gx_lower; }
      
    /// Get the upper limit
    fp_t get_upper() { return gx_upper; }
    
    /** \brief Set the information for the solver

        This function currently always returns \ref success.
    */
    int set(func_t &ff, fp_t lower, fp_t upper) {
      
      if (lower > upper) {
        fp_t tmp=lower;
        lower=upper;
        upper=tmp;
      }
	
      gx_lower=lower;
      gx_upper=upper;
      fp_t two=2;
      groot=(lower+upper)/two;
  
      fp_t f_lower, f_upper;
    
      f_lower=ff(gx_lower);
      f_upper=ff(gx_upper);
	
      ga=gx_lower;
      gfa=f_lower;
	
      gb=gx_upper;
      gfb=f_upper;
	
      gc=gx_upper;
      gfc=f_upper;
	
      gd=gx_upper-gx_lower;
      ge=gx_upper-gx_lower;
	
      if ((f_lower<0.0 && f_upper<0.0) || 
          (f_lower>0.0 && f_upper>0.0)) {
        O2SCL_CONV2_RET("Endpoints don't straddle y=0 in ",
                        "root_brent_gsl::set().",o2scl::exc_einval,
                        this->err_nonconv);
      }
	
      return o2scl::success;
	
    }

#ifndef O2SCL_NO_BOOST_MULTIPRECISION
    
    /** \brief Set the information for the solver
        (adaptive multiprecision version)
    */
    template<typename func2_t, class fp2_t>
    int set_multip(func2_t &ff, fp2_t lower, fp2_t upper,
                   double func_tol, fp2_t storage[11]) {

      fp2_t &lroot=storage[0];
      fp2_t &lx_lower=storage[1];
      fp2_t &lx_upper=storage[2];
      fp2_t &la=storage[3];
      fp2_t &lb=storage[4];
      fp2_t &lc=storage[5];
      fp2_t &ld=storage[6];
      fp2_t &le=storage[7];
      fp2_t &lfa=storage[8];
      fp2_t &lfb=storage[9];
      fp2_t &lfc=storage[10];

      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;
      
      if (lower > upper) {
        fp2_t tmp=lower;
        lower=upper;
        upper=tmp;
      }
	
      lx_lower=lower;
      lx_upper=upper;
      fp2_t two=2;
      lroot=(lower+upper)/two;
  
      fp2_t f_lower, f_upper;
    
      fp2_t err;
      int fm_ret1=fm2.eval_tol_err(ff,lx_lower,f_lower,err);
      int fm_ret2=fm2.eval_tol_err(ff,lx_upper,f_upper,err);
      
      la=lx_lower;
      lfa=f_lower;
	
      lb=lx_upper;
      lfb=f_upper;
	
      lc=lx_upper;
      lfc=f_upper;
	
      ld=lx_upper-lx_lower;
      le=lx_upper-lx_lower;
	
      if ((f_lower<0.0 && f_upper<0.0) || 
          (f_lower>0.0 && f_upper>0.0)) {
        O2SCL_CONV2_RET("Endpoints don't straddle y=0 in ",
                        "root_brent_gsl::set().",o2scl::exc_einval,
                        this->err_nonconv);
      }
	
      return o2scl::success;
	
    }

#endif

  protected:
      
    /// The present solution estimate
    fp_t groot;
    /// The present lower limit
    fp_t gx_lower;
    /// The present upper limit
    fp_t gx_upper;

    /// \name Storage for solver state
    //@{
    fp_t ga, gb, gc, gd, ge;
    fp_t gfa, gfb, gfc;
    //@}
      
  };

}

#endif
