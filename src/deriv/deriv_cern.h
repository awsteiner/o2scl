/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
/** \file deriv_cern.h
    \brief File defining \ref o2scl::deriv_cern
*/
#ifndef O2SCL_DERIV_CERN_H
#define O2SCL_DERIV_CERN_H

#include <o2scl/deriv.h>
#include <o2scl/funct.h>
#include <o2scl/string_conv.h>
#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Numerical differentiation routine (CERNLIB)

      This uses Romberg extrapolation to compute the 
      derivative with the finite-differencing formula
      \f[
      f^{\prime}(x) = [f(x+h)-f(x-h)]/(2h)
      \f]

      If \ref deriv_base::verbose is greater than zero, then each iteration
      prints out the extrapolation table, and if \ref deriv_base::verbose
      is greater than 1, then a keypress is required at the end of
      each iteration.

      For sufficiently difficult functions, the derivative computation 
      can fail, and will call the error handler and return zero with
      zero error.

      Based on the CERNLIB routine DERIV, which was 
      based on \ref Rutishauser63 and is documented at
      http://wwwasdoc.web.cern.ch/wwwasdoc/shortwrupsdir/d401/top.html
      
      An example demonstrating the usage of this class is 
      given in <tt>examples/ex_deriv.cpp</tt> and the \ref ex_deriv_sect .

      If \ref deriv_base::verbose is greater than zero, at each iteration
      this class prints something similar to
      \verbatim
      deriv_cern, iteration: 1
      (hh, a[k], derivative) list: 
      -4.193459e-05 4.387643e-14 8.775286e-01
      -2.995402e-05 4.387792e-14 8.775585e-01
      -1.048405e-05 4.387845e-14 8.775690e-01
      -7.488654e-06 4.387882e-14 8.775765e-01
      -2.621038e-06 4.387895e-14 8.775791e-01
      -1.872173e-06 4.387905e-14 8.775810e-01
      -6.552611e-07 4.387908e-14 8.775817e-01
      -4.680438e-07 4.387910e-14 8.775821e-01
      -1.638153e-07 4.387911e-14 8.775823e-01
      \endverbatim
      If \ref deriv_base::verbose is greater than 1, a keypress is required
      after each iteration.

      \note Second and third derivatives are computed by naive nested
      applications of the formula for the first derivative.
      No uncertainty for these derivatives is provided. 

      \future All of the coefficients appear to be fractions which
      could be replaced with exact representation?
      \future Record the number of function calls?
      
      \comment
      - Maybe we should consider moving the table size to a template 
      parameter? (1/29/07 - Probably not, as we'd have to re-derive
      the coefficients for sizes other than 10)
      \endcomment
  */
  template<class func_t=funct> class deriv_cern : 
  public deriv_base<func_t> {

    public:
  
    /// A scaling factor (default 1.0)
    double delta;
  
    /// Extrapolation tolerance (default is \f$ 5 \times 10^{14} \f$)
    double eps;

    deriv_cern() {
      dx[0]=0.0256;
      dx[1]=0.0192;
      dx[2]=0.0128;
      dx[3]=0.0096;
      dx[4]=0.0064;
      dx[5]=0.0048;
      dx[6]=0.0032;
      dx[7]=0.0024;
      dx[8]=0.0016;
      dx[9]=0.0012;
  
      w[1][1]=1.3333333333333333;
      w[3][1]=1.0666666666666667;
      w[5][1]=1.0158730158730159;
      w[7][1]=1.0039215686274510;
  
      w[2][1]=3.3333333333333333e-1;
      w[4][1]=6.6666666666666667e-2;
      w[6][1]=1.5873015873015873e-2;
      w[8][1]=3.9215686274509804e-3;

      w[0][2]=2.2857142857142857;
      w[2][2]=1.1636363636363636;
      w[4][2]=1.0364372469635628;
      w[6][2]=1.0088669950738916;
      w[8][2]=1.0022021042329337;

      w[1][2]=1.2857142857142857;
      w[3][2]=1.6363636363636364e-1;
      w[5][2]=3.6437246963562753e-2;
      w[7][2]=8.8669950738916256e-3;
      w[9][2]=2.2021042329336922e-3;
  
      w[0][3]=1.8000000000000000;
      w[2][3]=1.1250000000000000;
      w[4][3]=1.0285714285714286;
      w[6][3]=1.0069930069930070;
      w[8][3]=1.0017391304347826;
  
      w[1][3]=8.0000000000000000e-1;
      w[3][3]=1.2500000000000000e-1;
      w[5][3]=2.8571428571428571e-2;
      w[7][3]=6.9930069930069930e-3;
      w[9][3]=1.7391304347826087e-3;
  
      delta=1.0;
      eps=5.0e-14;
    }
    
    /** \brief Calculate the first derivative of \c func w.r.t. x and the
	uncertainty
    */
    virtual int deriv_err(double x, func_t &func,
			 double &dfdx, double &err) {
      return deriv_tlate<func_t>(x,func,dfdx,err);
    }

    /// Return string denoting type ("deriv_cern")
    virtual const char *type() { return "deriv_cern"; }

  protected:

#ifndef DOXYGEN_INTERNAL
 
    /** \brief Internal template version of the derivative function
    */
    template<class func2_t> int deriv_tlate(double x, func2_t &func,
					  double &dfdx, double &err) {
      
      double t[10][10], a[10], del, hh;
      bool lev[10]={1,0,1,0,1,0,1,0,1,0}, lmt;
      int is, k, m;

      del=10.0*fabs(delta);
      is=10;
  
      do {
	is--;
	del=del/10.0;

	if (is==0 || x+del*dx[9]==x) {
	  delta=del;
	  dfdx=0.0;
	  err=0.0;
	  std::string str="Calculation of derivative failed (is="+
	    itos(is)+" and del*dx[9]="+dtos(del*dx[9])+
	    ") in deriv_cern::deriv_tlate().";
	  O2SCL_ERR_RET(str.c_str(),exc_efailed);
	}
    
	for(k=0;k<=9;k++) {
	  hh=del*dx[k];
	  t[k][0]=(func(x+hh)-func(x-hh))/(hh+hh);
	  a[k]=t[k][0];
	}
    
	if (a[0]>=a[9]) {
	  for(k=0;k<=9;k++) a[k]=-a[k];
	}
    
	lmt=true;
	for(k=1;k<=9;k++) {
	  hh=a[k-1]-a[k];
	  lmt=(lmt && (hh<=0.0 || fabs(hh)<=eps*fabs(a[k])));
	}
	
	if (this->verbose>0) {
	  std::cout << "deriv_cern, iteration: " << 10-is << std::endl;
	  std::cout << "(hh, a[k], derivative) list: " 
		    << std::endl;
	  for(k=1;k<=9;k++) {
	    std::cout << a[k-1]-a[k] << " " << eps*fabs(a[k]) << " "
		      << t[k][0] << std::endl;
	  }
	  std::cout << "Converged: " << lmt << std::endl;
	  if (this->verbose>1) {
	    char ch;
	    std::cin >> ch;
	  }
	}

      } while (lmt==false);
  
      for(m=1;m<=9;m++) {
	for(k=0;k<=9-m;k++) {
	  if (lev[m]) {
	    t[k][m]=w[m-1][1]*t[k+1][m-1]-w[m][1]*t[k][m-1];
	  } else if (lev[k]) {
	    t[k][m]=w[m-1][2]*t[k+1][m-1]-w[m][2]*t[k][m-1];
	  } else {
	    t[k][m]=w[m-1][3]*t[k+1][m-1]-w[m][3]*t[k][m-1];
	  }
	}
      }
      dfdx=t[0][9];
      if (dfdx!=0.0) err=(dfdx-t[0][8])/dfdx;
      else err=0.0;
      delta=del;
  
      return 0;
    }

    /** \brief Calculate the first derivative of \c func w.r.t. x

	This is an internal version of deriv() which is used in
	computing second and third derivatives
     */
    virtual int deriv_err_int(double x, funct &func, double &dfdx, 
			     double &err) {
      return deriv_tlate<funct>(x,func,dfdx,err);
    }
 
    /// \name Storage for the fixed coefficients
    //@{
    double dx[10];
    double w[10][4];
    //@}

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif


