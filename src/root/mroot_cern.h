/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
/** \file mroot_cern.h
    \brief File defining \ref o2scl::mroot_cern
*/
#ifndef O2SCL_MROOT_CERN_H
#define O2SCL_MROOT_CERN_H

#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/vector.h>
#include <o2scl/mroot.h>

namespace o2scl {

  /** \brief Multi-dimensional mroot-finding routine (CERNLIB)
	
      If \f$ x_i \f$ denotes the current iteration, and \f$
      x^{\prime}_i \f$ denotes the previous iteration, then the
      calculation is terminated if either of the following tests is
      successful
      \f[
      1:\quad \mathrm{max} | f_i(x) | \leq \mathrm{tol\_rel}
      \f]
      \f[
      2:\quad \mathrm{max} |x_i-x^{\prime}_i| \leq
      \mathrm{tol\_abs} \times \mathrm{max} | x_i |
      \f]

      This routine treats the functions specified as a \ref mm_funct11
      object slightly differently than \ref o2scl::mroot_hybrids. First
      the equations should be numbered (as much as is possible) in
      order of increasing nonlinearity. Also, instead of calculating
      all of the equations on each function call, only the equation
      specified by the \c size_t parameter needs to be calculated. If
      the equations are specified as
      \f{eqnarray*}
      &0=f_0(x_0,x_1,...,x_{n-1})& \\
      &0=f_1(x_0,x_1,...,x_{n-1})& \\
      &...& \\
      &0=f_{n-1}(x_0,x_1,...,x_{n-1})& \\
      \f}
      then when the \c size_t argument is given as \c i, then
      only the function \f$ f_i \f$ needs to be calculated.

      \warning This code has not been checked to ensure that it cannot
      fail to solve the equations without calling the error handler
      and returning a non-zero value. Until then, the solution may
      need to be checked explicitly by the caller.

      See the \ref multisolve_subsect section of the User's guide for
      general information about \o2 solvers. 
      There is an example for the usage of the multidimensional solver
      classes given in <tt>examples/ex_mroot.cpp</tt>, see \ref
      ex_mroot_sect .

      \future Modify this so it handles functions which return
      non-zero values.
      \future Move some of the memory allocation out of msolve()
      \future Give the user access to the number of function
      calls
      \future Rename nier6, nier7, and nier8 to something sensible.
      \future It may be that the \o2 native Householder transformations
      should be used here instead of the inline version given here.

      Based on the CERNLIB routines RSNLEQ and DSNLEQ, which was 
      based on \ref More79 and \ref More80 and is documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/c201/top.html
  */
#ifdef O2SCL_NO_CPP11 
  template<class func_t=mm_funct<boost::numeric::ublas::vector<double> >, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class jfunc_t=jac_funct<vec_t,boost::numeric::ublas::matrix<double> > > 
    class mroot_cern : public mroot<func_t,vec_t,jfunc_t>
#else
    template<class func_t=mm_funct11,
    class vec_t=boost::numeric::ublas::vector<double>, 
    class jfunc_t=jac_funct11> class mroot_cern : 
    public mroot<func_t,vec_t,jfunc_t>
#endif
    {
      
#ifndef DOXYGEN_INTERNAL
      
    protected:
    
    /// Desc
    boost::numeric::ublas::matrix<double> w;
    
#endif
    
    public:
    
    mroot_cern() {
      info=0;
      eps=0.1490116119384766e-07;
      scale=10.0;
      maxf=0;
	
      int tmp_mpt[289]=
	{0,1,2,3,3,3,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,9,10,
	 10,10,10,10,11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,14,14,
	 14,14,14,15,15,15,15,15,15,16,16,16,16,16,16,17,17,17,17,17,18,
	 18,18,18,18,18,19,19,19,19,19,19,20,20,20,20,20,20,21,21,21,21,
	 21,21,21,22,22,22,22,22,22,23,23,23,23,23,23,24,24,24,24,24,24,
	 24,25,25,25,25,25,25,26,26,26,26,26,26,26,27,27,27,27,27,27,28,
	 28,28,28,28,28,28,29,29,29,29,29,29,29,30,30,30,30,30,30,30,31,
	 31,31,31,31,31,31,32,32,32,32,32,32,32,33,33,33,33,33,33,33,34,
	 34,34,34,34,34,34,35,35,35,35,35,35,35,36,36,36,36,36,36,36,37,
	 37,37,37,37,37,37,37,38,38,38,38,38,38,38,39,39,39,39,39,39,39,
	 40,40,40,40,40,40,40,40,41,41,41,41,41,41,41,42,42,42,42,42,42,
	 42,42,43,43,43,43,43,43,43,44,44,44,44,44,44,44,44,45,45,45,45,
	 45,45,45,45,46,46,46,46,46,46,46,47,47,47,47,47,47,47,47,48,48,
	 48,48,48,48,48,48};
      // The for loop is just a convenient way of using
      // aggregate initialization
      for(size_t i=0;i<289;i++) mpt[i]=tmp_mpt[i];
    }
      
    /** \brief Get the value of \c INFO from the last call to msolve()
	  
	The value of info is assigned according to the following list.
	The values 1-8 are the standard behavior from CERNLIB.
	0 - The function solve() has not been called.
	1 - Test 1 was successful. \n
	2 - Test 2 was successful. \n
	3 - Both tests were successful. \n
	4 - Number of iterations is greater than mroot_cern_root::maxf. \n
	5 - Approximate (finite difference) Jacobian matrix is
	singular. \n
	6 - Iterations are not making good progress. \n
	7 - Iterations are diverging. \n
	8 - Iterations are converging, but either mroot_cern_root::tol_abs
	is too small or the Jacobian is nearly singular
	or the variables are badly scaled. \n
	9 - Either root::tol_rel or root::tol_abs is not greater than zero
	or the specified number of variables is \f$ \leq 0\f$.

	The return values returned by msolve() corresponding
	to the values of \c INFO above are
	1 - \ref success
	2 - \ref success
	3 - \ref success
	4 - \ref exc_emaxiter
	5 - \ref exc_esing
	6 - \ref exc_enoprog
	7 - \ref exc_erunaway
	8 - \ref exc_efailed
	9 - \ref exc_einval
    */
    int get_info() { return info; }

    /** \brief Get the a string corresponding to the integer returned
	by \ref mroot_cern::get_info().
    */
    std::string get_info_string() {
      if (info==0) {
	return "The function solve() has not been called.";
      } else if (info==1) {
	return "Test 1 was successful.";
      } else if (info==2) {
	return "Test 2 was successful.";
      } else if (info==3) {
	return "Both tests were successful.";
      } else if (info==4) {
	return "Number of iterations is greater than maxf.";
      } else if (info==5) {
	return "Approximate Jacobian matrix is singular.";
      } else if (info==6) {
	return "Iterations are not making good progress.";
      } else if (info==7) {
	return "Iterations are diverging.";
      } else if (info==8) {
	return "Either tol_abs is too small or Jacobian is nearly singular.";
      } else if (info==9) {
	return "Either tol_rel, tol_abs, or the number of vars is not positive.";
      }
    }

    /** \brief Maximum number of function evaluations

	If \f$ \mathrm{maxf}\leq 0 \f$ , then \f$ 50(\mathrm{nv}+3) \f$ 
	(which is the CERNLIB default) is used.  The default value of
	\c maxf is zero which then implies the default from CERNLIB.
    */
    int maxf;

    /// Return the type, \c "mroot_cern".
    virtual const char *type() { return "mroot_cern"; }

    /** \brief The original scale parameter from CERNLIB (default 10.0)
     */
    double scale;
    
    /** \brief The smallest floating point number
	(default \f$ \sim 1.49012 \times 10^{-8} \f$ )

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
    double eps;
    
    /// Solve \c func using \c x as an initial guess, returning \c x.
    virtual int msolve(size_t nvar, vec_t &x, func_t &func) {

      int mopt=0, i, j, k, it=0;
      double fky;

      int lmaxf;
      if (maxf<=0) lmaxf=50*(nvar+3);
      else lmaxf=maxf;
  
      info=0;
  
      if (nvar<=0 || this->tol_rel<=0.0 || this->tol_abs<=0.0) {
	info=9;
	std::string str="Invalid value of tol_rel ("+dtos(this->tol_rel)+
	  "), tol_abs ("+dtos(this->tol_abs)+"), or nvar ("+itos(nvar)+
	  " in mroot_cern::msolve().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }
  
      // Find optimal value of mopt for iterative refinement

      if (nvar<=288) mopt=mpt[nvar-1];
      else {
	bool done=false;
	double h=0.0;
	for(i=49;i<=((int)nvar) && done==false;i++) {
	  double temp=log(((double)i)+1.0)/((double)(nvar+2*i+1));
	  if (temp<h) {
	    mopt=i-1;
	    done=true;
	  }
	  if (!done) h=temp;
	}
      }
  
      int iflag=0, numf=0, nfcall=0, nier6=-1, nier7=-1, nier8=0;
      double fnorm=0.0, difit=0.0, xnorm=0.0;
      bool set=false;

      for(i=0;i<((int)nvar);i++) {
	if (xnorm<fabs(x[i])) {
	  xnorm=fabs(x[i]);
	  set=true;
	}
      }
      double delta=scale*xnorm;
      if (set==false) delta=scale;

      w.resize(nvar,nvar);

      vec_t f(nvar), w0(nvar), w1(nvar), w2(nvar);

      bool solve_done=false;
      while (solve_done==false) {
	bool bskip=false;
    
	int nsing=nvar;
	double fnorm1=fnorm;
	double difit1=difit;
	fnorm=0.0;
    
	// Compute step H for the divided difference which approximates
	// the K-th row of the Jacobian matrix
    
	double h=eps*xnorm;
	if (h==0.0) h=eps;
	for(j=0;j<((int)nvar);j++) {
	  for(i=0;i<((int)nvar);i++) {
	    w(j,i)=0.0;
	  }
	  w(j,j)=h;
	  w1[j]=x[j];
	}

	// Enter a subiteration
    
	for(k=0;k<((int)nvar);k++) {
	  iflag=k;
	    
	  func(iflag,w1,f);
	    
	  fky=f[k];
	  nfcall++;
	  numf=(int)(((double)nfcall)/((double)nvar));
	  if (fnorm<fabs(fky)) fnorm=fabs(fky);
      
	  // Compute the K-th row of the Jacobian matrix

	  for(j=k;j<((int)nvar);j++) {
	    for(i=0;i<((int)nvar);i++) {
	      w2[i]=w1[i]+w(j,i);
	    }

	    func(iflag,w2,f);

	    double fkz=f[k];
	    nfcall++;
	    numf=(int)(((double)nfcall)/((double)nvar));
	    w0[j]=fkz-fky;
	  }

	  f[k]=fky;
      
	  // Compute the Householder transformation to reduce the K-th row
	  // of the Jacobian matrix to a multiple of the K-th unit vector

	  double eta=0.0;
	  for(i=k;i<((int)nvar);i++) if (eta<fabs(w0[i])) eta=fabs(w0[i]);
	
	  if (eta!=0.0) {
	    nsing--;
	    double sknorm=0.0;
	    for(i=k;i<((int)nvar);i++) {
	      w0[i]/=eta;
	      sknorm+=w0[i]*w0[i];
	    }
	    sknorm=sqrt(sknorm);
	    if (w0[k]<0.0) sknorm=-sknorm;
	    w0[k]+=sknorm;
	  
	    // Apply the transformation

	    for(i=0;i<((int)nvar);i++) {
	      w2[i]=0.0;
	    }
	    for(j=k;j<((int)nvar);j++) {
	      for(i=0;i<((int)nvar);i++) {
		w2[i]+=w0[j]*w(j,i);
	      }
	    }
	    for(j=k;j<((int)nvar);j++) {
	      double temp=w0[j]/(sknorm*w0[k]);
	      for(i=0;i<((int)nvar);i++) {
		w(j,i)-=temp*w2[i];
	      }
	    }

	    // Compute the subiterate

	    w0[k]=sknorm*eta;
	    double temp2=fky/w0[k];
	    if (h*fabs(temp2)>delta) 
	      temp2=(temp2>=0.0) ? fabs(delta/h) : -fabs(delta/h);
	    for(i=0;i<((int)nvar);i++) {
	      w1[i]+=temp2*w(k,i);
	    }
	  }
	}

	// Compute the norms of the iterate and correction vector

	xnorm=0.0;
	difit=0.0;
	for(i=0;i<((int)nvar);i++) {
	  if (xnorm<fabs(w1[i])) xnorm=fabs(w1[i]);
	  if (difit<fabs(x[i]-w1[i])) difit=fabs(x[i]-w1[i]);
	  x[i]=w1[i];
	}

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
	  this->print_iter(nvar,x,f,++it,fnorm,this->tol_rel,"mroot_cern");
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
	  O2SCL_CONV_RET("Too many iterations in mroot_cern::msolve().",
			 exc_emaxiter,this->err_nonconv);
	}
	if (nsing==((int)nvar)) {
	  info=5;
	  O2SCL_CONV_RET("Jacobian matrix singular in mroot_cern::msolve().",
			 exc_esing,this->err_nonconv);
	}
	if (nier6==5) {
	  info=6;
	  O2SCL_CONV_RET("No progress in mroot_cern::msolve().",
			 exc_enoprog,this->err_nonconv);
	}
	if (nier7==3) {
	  info=7;
	  O2SCL_CONV_RET("Iterations diverging in mroot_cern::msolve().",
			 exc_erunaway,this->err_nonconv);
	}
	if (nier8==4) {
	  info=8;
	  std::string s="Variable tol_abs too small, J singular, or bad ";
	  s+="scaling in mroot_cern::msolve().";
	  O2SCL_CONV_RET(s.c_str(),exc_efailed,this->err_nonconv);
	}

	// Exit if necessary

	if (info!=0) {
	  O2SCL_ERR("Unspecified error in mroot_cern::msolve().",
			exc_efailed);
	}

	if (!((!lcv) || difit>0.05*xnorm)) {
	  // 8/20/08: Could this just be rewritten? 
	  // if (lcv && difit<=0.05*xnorm)
      
	  // Iterative refinement (if the iteration is converging)

	  for(int m=2;m<=mopt && bskip==false;m++) {
	    fnorm1=fnorm;
	    fnorm=0.0;
	    for(k=0;k<((int)nvar) && bskip==false;k++) {
	      iflag=k;

	      func(iflag,w1,f);
	  
	      fky=f[k];
	      nfcall++;
	      numf=(int)(((double)nfcall)/((double)nvar));

	      if (fnorm<fabs(fky)) fnorm=fabs(fky);
	
	      // Iterative refinement is terminated if it does not give a
	      // reduction on residuals
	  
	      if (fnorm>=fnorm1) {
		fnorm=fnorm1;
		bskip=true;
	      } 

	      if (!bskip) {
		double temp3=fky/w0[k];
	    
		for(i=0;i<((int)nvar);i++) {
		  w1[i]+=temp3*w(k,i);
		}
	      }
	    }
	
	    if (!bskip) {

	      // Compute the norms of the iterate and correction vector

	      xnorm=0.0;
	      difit=0.0;

	      for(i=0;i<((int)nvar);i++) {
		if (xnorm<fabs(w1[i])) xnorm=fabs(w1[i]);
		if (difit<fabs(x[i]-w1[i])) difit=fabs(x[i]-w1[i]);
		x[i]=w1[i];
	      }

	      // Stopping criteria for iterative refinement

	      if (fnorm<=this->tol_rel) info=1;
	      if (difit<=xnorm*this->tol_abs) info=2;
	      if (fnorm<=this->tol_rel && info==2) info=3;
	      if (numf>=lmaxf && info==0) {
		info=4;
		O2SCL_CONV_RET("Too many iterations in mroot_cern::msolve().",
			       exc_emaxiter,this->err_nonconv);
	      }

	      if (info!=0) {
		return 0;
	      }
	    }
	  }
	}
      }
      
      return 0;
    }

#ifndef DOXYGEN_INTERNAL

    protected:
      
    /// Internal storage for the value of \c info
    int info;
      
    /// Store the number of function evaluations
    int mpt[289];
      
#endif
      
    };
  

}

#endif

