/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
/** \file root_cern.h
    \brief File defining \ref o2scl::root_cern
*/
#ifndef O2SCL_ROOT_CERN_H
#define O2SCL_ROOT_CERN_H

#include <string>

#include <o2scl/string_conv.h>
#include <o2scl/root.h>

namespace o2scl {

  /** \brief One-dimensional version of cern_mroot
      
      This one-dimensional root-finding routine, based on \ref
      o2scl::mroot_cern, is probably slower than the more typical 1-d
      routines, but also tends to converge for a larger class of
      functions than \ref o2scl::root_bkt_cern, \ref
      o2scl::root_brent_gsl, or \ref o2scl::root_stef. It has been
      modified from \ref o2scl::mroot_cern and slightly optimized, but
      has the same basic behavior.

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
      
      \note This code has not been checked to ensure that it cannot
      fail to solve the equations without calling the error handler
      and returning a non-zero value. Until then, the solution may
      need to be checked explicitly by the caller.

      See the \ref onedsolve_subsect section of the User's guide for
      general information about \o2 solvers. 

      \future Double-check this class to make sure it cannot fail
      while returning 0 for success.
  */
#ifdef DOXYGEN_NO_O2NS
  template<class func_t=funct11> class root_cern : public root
#else
    template<class func_t=funct11> class root_cern : 
  public root<func_t> 
#endif
  {
    
    public:
    
    root_cern() {
      info=0;
      eps=0.1490116119384766e-07;
      scale=10.0;
      maxf=0;
    }

    virtual ~root_cern() {}
    
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
    double scale;
    
    /** \brief The smallest floating point number
	(default \f$ \sim 1.49012 \times 10^{-8} \f$)

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

	\future This number should probably default to one of the
	GSL tolerances.
    */
    double eps;
    
    /// Solve \c func using \c x as an initial guess, returning \c x.
    virtual int solve(double &ux, func_t &func) {
      
      int it=0;
      double fky;

      int lmaxf;
      if (maxf<=0) lmaxf=200;
      else lmaxf=maxf;
      
      info=0;
  
      if (this->tol_rel<=0.0 || this->tol_abs<=0.0) {
	info=9;
	std::string str="Invalid value of tol_rel ("+dtos(this->tol_rel)+
	  ") or tol_abs ("+dtos(this->tol_abs)+") in root_cern::solve().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }
  
      int iflag=0, numf=0, nfcall=0, nier6=-1, nier7=-1, nier8=0;
      double fnorm=0.0, difit=0.0, xnorm=0.0;
      bool set=false;
	
      if (xnorm<fabs(ux)) {
	xnorm=fabs(ux);
	set=true;
      }

      double delta=scale*xnorm;
      if (set==false) delta=scale;
	
      double wmat, farr, w0arr, w1arr, w2arr;
	
      bool solve_done=false;
      while (solve_done==false) {
    
	int nsing=1;
	double fnorm1=fnorm;
	double difit1=difit;
	fnorm=0.0;
    
	// Compute step H for the divided difference which approximates
	// the K-th row of the Jacobian matrix
	
	double h=eps*xnorm;
	if (h==0.0) h=eps;

	wmat=h;
	w1arr=ux;

	// Enter a subiteration
    
	iflag=0;

	farr=func(w1arr);

	fky=farr;
	nfcall++;
	numf=nfcall;
	if (fnorm<fabs(fky)) fnorm=fabs(fky);
      
	// Compute the K-th row of the Jacobian matrix

	w2arr=w1arr+wmat;
	    
	farr=func(w2arr);
	    
	double fkz=farr;
	nfcall++;
	numf=nfcall;
	w0arr=fkz-fky;
	    
	farr=fky;
      
	// Compute the Householder transformation to reduce the K-th row
	// of the Jacobian matrix to a multiple of the K-th unit vector

	double eta=0.0;
	if (eta<fabs(w0arr)) eta=fabs(w0arr);
	
	if (eta!=0.0) {
	  nsing--;
	  double sknorm=0.0;
	      
	  w0arr/=eta;
	  sknorm+=w0arr*w0arr;

	  sknorm=sqrt(sknorm);
	  if (w0arr<0.0) sknorm=-sknorm;
	  w0arr+=sknorm;
	  
	  // Apply the transformation

	  w2arr=0.0;
	  w2arr+=w0arr*wmat;
	  double temp=w0arr/(sknorm*w0arr);
	  wmat-=temp*w2arr;

	  // Compute the subiterate

	  w0arr=sknorm*eta;
	  double temp2=fky/w0arr;
	  if (h*fabs(temp2)>delta) 
	    temp2=(temp2>=0.0) ? fabs(delta/h) : -fabs(delta/h);
	  w1arr+=temp2*wmat;
	}

	// Compute the norms of the iterate and correction vector

	xnorm=0.0;
	difit=0.0;
	  
	if (xnorm<fabs(w1arr)) xnorm=fabs(w1arr);
	if (difit<fabs(ux-w1arr)) difit=fabs(ux-w1arr);
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
	  s2+="or bad scaling in cerm_mroot_root::solve().";
	  O2SCL_CONV(s2.c_str(),exc_efailed,this->err_nonconv);
	}
	
	// Exit if necessary

	if (info!=0) return exc_efailed;
	  
      }
	
      return 0;
    }
      
#ifndef DOXYGEN_INTERNAL

    protected:
      
    /// Internal storage for the value of \c info
    int info;
      
#endif
      
  };
  
}
  
#endif
