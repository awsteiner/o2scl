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
/** \file root_cern.h
    \brief File defining \ref o2scl::root_cern
*/
#ifndef O2SCL_ROOT_CERN_H
#define O2SCL_ROOT_CERN_H

#include <string>

#include <o2scl/string_conv.h>
#include <o2scl/root.h>
#include <o2scl/table.h>

namespace o2scl {

  /** \brief One-dimensional version of cern_mroot
      
      This one-dimensional root-finding routine, based on \ref
      o2scl::mroot_cern. Previous testing has suggested that it is
      probably slower than the more typical 1-D solvers, but also
      tends to converge for a larger class of functions than \ref
      o2scl::root_bkt_cern, \ref o2scl::root_cern, or \ref
      o2scl::root_stef. It has been modified from \ref
      o2scl::mroot_cern and slightly optimized, but has the same basic
      behavior.

      If \f$ x_i \f$ denotes the current iteration, and \f$
      x^{\prime}_i \f$ denotes the previous iteration, then the
      calculation is terminated if either (or both) of the following
      tests is successful
      \f[
      1:\quad \mathrm{max} | f_i(x) | \leq \mathrm{tol\_rel}
      \f]
      \f[
      2:\quad \mathrm{max} |x_i-x^{\prime}_i| \leq
      \mathrm{tol\_abs} \times \mathrm{max} | x_i |
      \f]

      \note This class has difficulty finding the root when
      the desired root is near 0 (AWS 1/22/19)
      
      \note This code has not been checked to ensure that it cannot
      fail to solve the equations without calling the error handler
      and returning a non-zero value. Until then, the solution may
      need to be checked explicitly by the caller.

      \verbatim embed:rst
      See the :ref:`One-dimensional solvers` section of the User's
      guide for general information about O2scl solvers.

      .. todo::

         class root_cern

         Future:

         - Double-check this class to make sure it cannot fail
           while returning 0 for success.

      \endverbatim
  */
  template<class func_t=funct, class fp_t=double> class root_cern :
    public root<func_t,func_t,fp_t> {
    
    public:
    
    root_cern() {
      info=0;
      // The original value in the CERNLIB code was:
      // eps=0.1490116119384766e-07;
      // This is an approximately equivalent replacement.
      eps=sqrt(std::numeric_limits<fp_t>::epsilon());
      scale=10.0;
      maxf=0;
      store_funcs=false;
      pow_tol_func=1.33;
    }

    virtual ~root_cern() {}
    
    /// For storing function evaluations
    o2scl::table<std::vector<fp_t>,fp_t> tab;

    /** \brief Power for tolerance of function evaluations
        with adaptive multiprecision (default 1.33)
    */
    double pow_tol_func;

    /** \brief Get the value of \c INFO from the last call to solve() 
	(default 0)
	
	The value of info is assigned according to the following list.
	The values 1-8 are the standard behavior from CERNLIB.
	0 - The function solve() has not been called.
	1 - Test 1 was successful. \n
	2 - Test 2 was successful. \n
	3 - Both tests were successful. \n
	4 - Number of iterations is greater than root_cern::maxf. \n
	5 - Approximate (finite difference) Jacobian matrix is
	singular. \n
	6 - Iterations are not making good progress. \n
	7 - Iterations are diverging. \n
	8 - Iterations are converging, but either root_cern::tol_abs
	is too small or the Jacobian is nearly singular
	or the variables are badly scaled. \n
	9 - Either root::tol_rel or root::tol_abs is not greater than zero.

    */
    int get_info() { return info; }

    /** \brief Maximum number of function evaluations
	
	If \f$ \mathrm{maxf}\leq 0 \f$, then 200 (which is the CERNLIB
	default) is used.  The default value of \c maxf is zero which
	then implies the default from CERNLIB.
    */
    int maxf;

    /// Return the type, \c "root_cern".
    virtual const char *type() { return "root_cern"; }

    /** \brief The original scale parameter from CERNLIB (default 10.0)
     */
    fp_t scale;
    
    /** \brief The square root of epsilon

	The constructor sets this value to
        <tt>sqrt(std::numeric_limits<fp_t>::epsilon())</tt>.
        The original prescription from CERNLIB for \c eps is
	given below:
	\verbatim
	#if !defined(CERNLIB_DOUBLE)
	PARAMETER (EPS =  0.84293 69702 17878 97282 52636 392E-07)
	#endif
	#if defined(CERNLIB_IBM)
	PARAMETER (EPS =  0.14901 16119 38476 562D-07)
	#endif
	#if defined(CERNLIB_VAX)
	PARAMETER (EPS =  0.37252 90298 46191 40625D-08)
	#endif
	#if (defined(CERNLIB_UNIX))&&(defined(CERNLIB_DOUBLE))
	PARAMETER (EPS =  0.14901 16119 38476 600D-07)
	#endif
        \endverbatim

    */
    fp_t eps;
    
    /// If true, store function evaluations
    bool store_funcs;
    
    /// Internal multiprecision version of the solve function
    template<typename func2_t, class fp2_t>
    int solve_int_multip(fp2_t &ux, func2_t &func,
                         double root_tol, double func_tol) {
      
      funct_multip fm2;
      fm2.err_nonconv=false;
      fm2.tol_rel=func_tol;
      
      if (store_funcs) {
        tab.clear();
        tab.line_of_names("x y");
      }
      
      int it=0;
      fp2_t fky;

      int lmaxf;
      if (maxf<=0) lmaxf=200;
      else lmaxf=maxf;
      
      info=0;
  
      if (this->tol_rel<=0 || this->tol_abs<=0) {
	info=9;
	std::string str="Invalid value of tol_rel ("+dtos(this->tol_rel)+
	  ") or tol_abs ("+dtos(this->tol_abs)+
          ") in root_cern::solve_int_multip().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }
  
      int iflag=0, numf=0, nfcall=0, nier6=-1, nier7=-1, nier8=0;
      fp2_t fnorm=0, difit=0, xnorm=0;
      bool set=false;
	
      if (xnorm<abs(ux)) {
	xnorm=abs(ux);
	set=true;
      }

      fp2_t delta=scale*xnorm;
      if (set==false) delta=scale;
	
      fp2_t wmat, farr, w0arr, w1arr, w2arr;
	
      bool solve_done=false;
      while (solve_done==false) {
    
	int nsing=1;
	fp2_t fnorm1=fnorm;
	fp2_t difit1=difit;
	fnorm=0;
    
	// Compute step H for the divided difference which approximates
	// the K-th row of the Jacobian matrix
	
	fp2_t h=eps*xnorm;
	if (h==0) h=eps;

	wmat=h;
	w1arr=ux;

	// Enter a subiteration
    
	iflag=0;

        fp2_t err;
        int fm_ret1=fm2.eval_tol_err(func,w1arr,farr,err);
        if (fm_ret1!=0) return 1;

        if (store_funcs) {
          std::vector<fp_t> line={static_cast<fp_t>(w1arr),
                                  static_cast<fp_t>(farr)};
          tab.line_of_data(line);
        }

	fky=farr;
	nfcall++;
	numf=nfcall;
	if (fnorm<abs(fky)) fnorm=abs(fky);
      
	// Compute the K-th row of the Jacobian matrix

	w2arr=w1arr+wmat;
	    
        int fm_ret2=fm2.eval_tol_err(func,w2arr,farr,err);
        if (fm_ret2!=0) return 2;
	    
        if (store_funcs) {
          std::vector<fp_t> line={static_cast<fp_t>(w2arr),
                                  static_cast<fp_t>(farr)};
          tab.line_of_data(line);
        }
        
	fp2_t fkz=farr;
	nfcall++;
	numf=nfcall;
	w0arr=fkz-fky;
	    
	farr=fky;
      
	// Compute the Householder transformation to reduce the K-th row
	// of the Jacobian matrix to a multiple of the K-th unit vector

	fp2_t eta=0;
	if (eta<abs(w0arr)) eta=abs(w0arr);
	
	if (eta!=0) {
	  nsing--;
	  fp2_t sknorm=0;
	      
	  w0arr/=eta;
	  sknorm+=w0arr*w0arr;

	  sknorm=sqrt(sknorm);
	  if (w0arr<0) sknorm=-sknorm;
	  w0arr+=sknorm;
	  
	  // Apply the transformation

	  w2arr=0;
	  w2arr+=w0arr*wmat;
	  fp2_t temp=w0arr/(sknorm*w0arr);
	  wmat-=temp*w2arr;

	  // Compute the subiterate

	  w0arr=sknorm*eta;
	  fp2_t temp2=fky/w0arr;
	  if (h*abs(temp2)>delta) {
            fp2_t arg1=abs(delta/h);
            fp2_t arg2=-abs(delta/h);
	    temp2=(temp2>=0) ? arg1 : arg2;
          }
	  w1arr+=temp2*wmat;
	}

	// Compute the norms of the iterate and correction vector

	xnorm=0;
	difit=0;
	  
	if (xnorm<abs(w1arr)) xnorm=abs(w1arr);
	if (difit<abs(ux-w1arr)) difit=abs(ux-w1arr);
	ux=w1arr;
	  
	// Update the bound on the correction vector

	if(delta<scale*xnorm) delta=scale*xnorm;
    
	// Determine the progress of the iteration

	bool lcv=(fnorm<fnorm1 && difit<difit1 && nsing==0);
	nier6++;
	nier7++;
	nier8++;
	if (lcv) nier6=0;
	if (fnorm<fnorm1 || difit<difit1) nier7=0;
	if (difit>eps*xnorm) nier8=0;

	// Print iteration information
	  
	if (this->verbose>0) {
	  this->print_iter(ux,farr,++it,fnorm,
                           static_cast<fp2_t>(this->tol_rel),
                           "root_cern");
	}
	
	// Tests for convergence

	if (fnorm<=this->tol_rel) info=1;
	if (difit<=this->tol_abs*xnorm && lcv) info=2;
	if (fnorm<=this->tol_rel && info==2) info=3;
	if (info!=0) {
	  if (!boost::math::isfinite(ux)) {
	    O2SCL_CONV2_RET("Solver converged to non-finite value ",
			    "in root_cern::solve_int_multip().",
                            exc_erange,this->err_nonconv);
	  }
	  return 0;
	}

	// Tests for termination

	if (numf>=lmaxf) {
	  info=4;
	  O2SCL_CONV2_RET("Too many iterations in ",
                          "root_cern::solve_int_multip().",
                          exc_emaxiter,this->err_nonconv);
	}
	if (nsing==1) {
	  info=5;
	  O2SCL_CONV2_RET("Jacobian matrix singular in ",
			  "root_cern::solve_int_multip().",
			  exc_esing,this->err_nonconv);
	}
	if (nier6==5) {
	  info=6;
	  O2SCL_CONV_RET("No progress in root_cern::solve_int_multip().",
			 exc_enoprog,this->err_nonconv);
	}
	if (nier7==3) {
	  info=7;
	  O2SCL_CONV_RET("Iterations diverging in root_cern::solve_int_multip().",
			 exc_erunaway,this->err_nonconv);
	}
	if (nier8==4) {
	  info=8;
	  std::string s2="Variable tol_abs too small, J singular, ";
	  s2+="or bad scaling in root_cern::solve_int_multip().";
	  O2SCL_CONV(s2.c_str(),exc_efailed,this->err_nonconv);
	}
	
	// Exit if necessary

	if (info!=0) return exc_efailed;
	  
      }
      
      if (!boost::math::isfinite(ux)) {
	O2SCL_CONV2_RET("Solver converged to non-finite value ",
			"in root_cern::solve_int_multip() (2).",exc_erange,
			this->err_nonconv);
      }
      return 0;
    }
    
    /** \brief Solve \c f in returning \c x1 (adaptive multiprecision
        version)
     */
    template<typename func2_t, class fp2_t>
    int solve_multip(fp2_t &x1, func2_t &&f, fp2_t &err,
                         double tol_loc=-1.0) {
      
      if (tol_loc<=0.0) {
        tol_loc=pow(10.0,-std::numeric_limits<fp2_t>::digits10);
      } 
      
      if (this->verbose>0) {
        std::cout << "solve_multip(): set "
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
        solve_int_multip(x1_d,f,tol_loc,func_tol);
                         
        if (this->verbose>0) {
          std::cout << "root_cern::solve_multip() "
                    << "double." << std::endl;
        }
        
        called_d=true;

        if (this->verbose>0) {
          std::cout << "root_cern::solve_multip() "
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
        solve_int_multip(x1_ld,f,tol_loc,func_tol);
                         
        if (this->verbose>0) {
          std::cout << "root_cern::solve_multip() "
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
          std::cout << "root_cern::solve_multip() "
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
        solve_int_multip(x1_cdf25,f,tol_loc,func_tol);
                         
        if (this->verbose>0) {
          std::cout << "root_cern::solve_multip() "
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
            std::cout << "root_cern::solve_multip() "
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
          std::cout << "root_cern::solve_multip() "
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
        solve_int_multip(x1_cdf35,f,tol_loc,func_tol);
                         
        if (this->verbose>0) {
          std::cout << "root_cern::solve_multip() "
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
          std::cout << "root_cern::solve_multip() "
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
        solve_int_multip(x1_cdf50,f,tol_loc,func_tol);
                         
        if (this->verbose>0) {
          std::cout << "root_cern::solve_multip() "
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
          std::cout << "root_cern::solve_multip() "
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
        solve_int_multip(x1_cdf100,f,tol_loc,func_tol);
                         
        if (this->verbose>0) {
          std::cout << "root_cern::solve_multip() "
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
          std::cout << "root_cern::solve_multip() "
                    << "failed after cpp_dec_float_100:\n  "
                    << dtos(x1_cdf50,0) << " "
                    << dtos(x1_cdf100,0) << " "
                    << dtos(err,0) << " " 
                    << tol_loc << std::endl;
        }

      }

#endif
      
      O2SCL_ERR2("Failed to compute with requested accuracy ",
                 "in root_cern::solve_multip().",
                 o2scl::exc_efailed);
      return o2scl::exc_efailed;
    }
    
    /// Solve \c func using \c x as an initial guess, returning \c x.
    virtual int solve(fp_t &ux, func_t &func) {
      
      if (store_funcs) {
        tab.clear();
        tab.line_of_names("x y");
      }
      
      int it=0;
      fp_t fky;

      int lmaxf;
      if (maxf<=0) lmaxf=200;
      else lmaxf=maxf;
      
      info=0;
  
      if (this->tol_rel<=0 || this->tol_abs<=0) {
	info=9;
	std::string str="Invalid value of tol_rel ("+dtos(this->tol_rel)+
	  ") or tol_abs ("+dtos(this->tol_abs)+") in root_cern::solve().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }
  
      int iflag=0, numf=0, nfcall=0, nier6=-1, nier7=-1, nier8=0;
      fp_t fnorm=0, difit=0, xnorm=0;
      bool set=false;
	
      if (xnorm<abs(ux)) {
	xnorm=abs(ux);
	set=true;
      }

      fp_t delta=scale*xnorm;
      if (set==false) delta=scale;
	
      fp_t wmat, farr, w0arr, w1arr, w2arr;
	
      bool solve_done=false;
      while (solve_done==false) {
    
	int nsing=1;
	fp_t fnorm1=fnorm;
	fp_t difit1=difit;
	fnorm=0;
    
	// Compute step H for the divided difference which approximates
	// the K-th row of the Jacobian matrix
	
	fp_t h=eps*xnorm;
	if (h==0) h=eps;

	wmat=h;
	w1arr=ux;

	// Enter a subiteration
    
	iflag=0;

	farr=func(w1arr);

        if (store_funcs) {
          std::vector<fp_t> line={w1arr,farr};
          tab.line_of_data(line);
        }

	fky=farr;
	nfcall++;
	numf=nfcall;
	if (fnorm<abs(fky)) fnorm=abs(fky);
      
	// Compute the K-th row of the Jacobian matrix

	w2arr=w1arr+wmat;
	    
	farr=func(w2arr);
	    
        if (store_funcs) {
          std::vector<fp_t> line={w2arr,farr};
          tab.line_of_data(line);
        }
        
	fp_t fkz=farr;
	nfcall++;
	numf=nfcall;
	w0arr=fkz-fky;
	    
	farr=fky;
      
	// Compute the Householder transformation to reduce the K-th row
	// of the Jacobian matrix to a multiple of the K-th unit vector

	fp_t eta=0;
	if (eta<abs(w0arr)) eta=abs(w0arr);
	
	if (eta!=0) {
	  nsing--;
	  fp_t sknorm=0;
	      
	  w0arr/=eta;
	  sknorm+=w0arr*w0arr;

	  sknorm=sqrt(sknorm);
	  if (w0arr<0) sknorm=-sknorm;
	  w0arr+=sknorm;
	  
	  // Apply the transformation

	  w2arr=0;
	  w2arr+=w0arr*wmat;
	  fp_t temp=w0arr/(sknorm*w0arr);
	  wmat-=temp*w2arr;

	  // Compute the subiterate

	  w0arr=sknorm*eta;
	  fp_t temp2=fky/w0arr;
	  if (h*abs(temp2)>delta) {
            fp_t arg1=abs(delta/h);
            fp_t arg2=-abs(delta/h);
	    temp2=(temp2>=0) ? arg1 : arg2;
          }
	  w1arr+=temp2*wmat;
	}

	// Compute the norms of the iterate and correction vector

	xnorm=0;
	difit=0;
	  
	if (xnorm<abs(w1arr)) xnorm=abs(w1arr);
	if (difit<abs(ux-w1arr)) difit=abs(ux-w1arr);
	ux=w1arr;
	  
	// Update the bound on the correction vector

	if(delta<scale*xnorm) delta=scale*xnorm;
    
	// Determine the progress of the iteration

	bool lcv=(fnorm<fnorm1 && difit<difit1 && nsing==0);
	nier6++;
	nier7++;
	nier8++;
	if (lcv) nier6=0;
	if (fnorm<fnorm1 || difit<difit1) nier7=0;
	if (difit>eps*xnorm) nier8=0;

	// Print iteration information
	  
	if (this->verbose>0) {
	  this->print_iter(ux,farr,++it,fnorm,this->tol_rel,
		     "root_cern");
	}
	
	// Tests for convergence

	if (fnorm<=this->tol_rel) info=1;
	if (difit<=this->tol_abs*xnorm && lcv) info=2;
	if (fnorm<=this->tol_rel && info==2) info=3;
	if (info!=0) {
	  if (!boost::math::isfinite(ux)) {
	    O2SCL_CONV2_RET("Solver converged to non-finite value ",
			    "in root_cern::solve().",exc_erange,
			    this->err_nonconv);
	  }
	  return 0;
	}

	// Tests for termination

	if (numf>=lmaxf) {
	  info=4;
	  O2SCL_CONV_RET("Too many iterations in root_cern::solve().",
			 exc_emaxiter,this->err_nonconv);
	}
	if (nsing==1) {
	  info=5;
	  O2SCL_CONV2_RET("Jacobian matrix singular in ",
			  "root_cern::solve().",
			  exc_esing,this->err_nonconv);
	}
	if (nier6==5) {
	  info=6;
	  O2SCL_CONV_RET("No progress in root_cern::solve().",
			 exc_enoprog,this->err_nonconv);
	}
	if (nier7==3) {
	  info=7;
	  O2SCL_CONV_RET("Iterations diverging in root_cern::solve().",
			 exc_erunaway,this->err_nonconv);
	}
	if (nier8==4) {
	  info=8;
	  std::string s2="Variable tol_abs too small, J singular, ";
	  s2+="or bad scaling in root_cern::solve().";
	  O2SCL_CONV(s2.c_str(),exc_efailed,this->err_nonconv);
	}
	
	// Exit if necessary

	if (info!=0) return exc_efailed;
	  
      }
      
      if (!boost::math::isfinite(ux)) {
	O2SCL_CONV2_RET("Solver converged to non-finite value ",
			"in root_cern::solve() (2).",exc_erange,
			this->err_nonconv);
      }
      return 0;
    }
      
    protected:
      
    /// Internal storage for the value of \c info
    int info;
      
  };
  
}
  
#endif
