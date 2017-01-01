/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2017, Andrew W. Steiner
  
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
#ifndef O2SCL_DERIV_EQI_H
#define O2SCL_DERIV_EQI_H

/** \file deriv_eqi.h
    \brief File defining \ref o2scl::deriv_eqi
*/

#include <cmath>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/deriv.h>
#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Derivatives for equally-spaced abscissas

      This is an implementation of the formulas for equally-spaced
      abscissas as indicated below. The level of approximation is
      specified in set_npoints(). The value of \f$ p \times h \f$ 
      can be specified in \c xoff (default is zero).

      \note Uncertainties are not computed and the code
      for second and third derivatives is unfinished. 

      \note The derivatives given, for example, from the
      five-point formula can sometimes be more accurate
      than computing the derivative from the interpolation class.
      This is especially true near the boundaries of the interpolated
      region.
      
      \future Finish the second and third derivative formulas.

      Two-point formula (note that this is independent of p).
      \f[
      f^{\prime}(x_0+p h)=\frac{1}{h}\left[
      f_{1}-f_{0} \right]
      \f]
      Three-point formula from Abramowitz and Stegun
      \f[
      f^{\prime}(x_0+p h)=\frac{1}{h}\left[
      \frac{2p-1}{2}f_{-1}-2 p f_{0}+\frac{2p+1}{2}f_{1}\right]
      \f]
      Four-point formula from Abramowitz and Stegun
      \f[
      f^{\prime}(x_0+p h)=\frac{1}{h}\left[
      -\frac{3 p^2-6 p+2}{6}f_{-1}
      +\frac{3 p^2-4 p -1}{2}f_{0}
      -\frac{3 p^2-2 p-2}{2}f_{1}
      +\frac{3 p^2-1}{6}f_{2}
      \right]
      \f]
      Five-point formula from Abramowitz and Stegun
      \f{eqnarray*}
      f^{\prime}(x_0+p h)&=&\frac{1}{h}\left[
      \frac{2 p^3-3 p^2-p+1}{12}f_{-2}
      -\frac{4 p^3-3p^2-8p+4}{6}f_{-1}
      \right. \\ && \left. 
      +\frac{2p^3-5p}{2}f_{0}
      -\frac{4p^3+3p^2-8p-4}{6}f_{1}
      \right. \\ && \left. 
      +\frac{2p^3+3p^2-p-1}{12}f_{2}
      \right]
      \f}

      The relations above can be confined to give formulas
      for second derivative formulas:
      Three-point formula 
      \f[
      f^{\prime}(x_0+p h)=\frac{1}{h^2}
      \left[f_{-1}-2 f_0+f_1\right]
      \f]
      Four-point formula:
      \f[
      f^{\prime}(x_0+p h)=\frac{1}{2 h^2}
      \left[\left(1-2p\right)f_{-1}-\left(1-6p\right)f_0
      -\left(1+6p\right)f_1+\left(1+2p\right)f_2\right]
      \f]
      Five-point formula:
      \f[
      f^{\prime}(x_0+p h)=\frac{1}{4 h^2}
      \left[\left(1-2p\right)^2f_{-2}
      +\left(8p-16 p^2\right)f_{-1}
      -\left(2-24 p^2\right)f_0
      -\left(8p+16p^2\right)f_1
      +\left(1+2p\right)^2 f_2\right]
      \f]
      Six-point formula:
      \f{eqnarray*}
      f^{\prime}(x_0+p h)&=&\frac{1}{12 h^2}\left[
      \left(2-10p+15 p^2-6p^3\right)f_{-2}
      +\left(3+14p-57p^2+30p^3\right)f_{-1}
      \right. \\ && \left. 
      +\left(-8+20p+78 p^2-60p^3\right)f_0
      +\left(-2-44p-42p^2+60p^3\right)f_1
      \right. \\ && \left. 
      +\left(6+22p+3p^2-30p^3\right)f_2
      +\left(-1-2p+3p^2+6p^3\right)f_3
      \right]
      \f}
      Seven-point formula:
      \f{eqnarray*}
      f^{\prime}(x_0+p h)&=&\frac{1}{36 h^2}\left[
      \left(4-24p+48p^2-36p^3+9p^4\right)f_{-3}
      +\left(12+12p-162p^2+180p^3-54p^4\right)f_{-2}
      \right. \\ && \left. 
      +\left(-15+120p+162p^2-360p^3+135p^4\right)f_{-1} 
      -4\left(8+48p-3p^2-90p^3+45p^4\right)f_0
      \right. \\ && \left. 
      +3\left(14+32p-36p^2-60p^3+45p^4\right)f_1
      +\left(-12-12p+54p^2+36p^3-54p^4\right)f_2 
      \right. \\ && \left. 
      +\left(1-6p^2+9p^4\right)f_3 
      \right]
      \f}
  */
  template<class func_t=funct11, 
    class vec_t=boost::numeric::ublas::vector<double> > 
    class deriv_eqi : public deriv_base<func_t> {
    public:

    deriv_eqi() {
      h=1.0e-4;
      xoff=0.0;
      cap=&deriv_eqi::deriv_vector5;
      cp=&deriv_eqi::derivp5;
      c2p=&deriv_eqi::deriv2p5;
    }

    /// Stepsize (Default  \f$ 10^{-4} \f$ ).
    double h;

    /// Offset (default 0.0)
    double xoff;

    /** \brief Set the number of points to use for first derivatives 
	(default 5)

	Acceptable values are 2-5 (see above).
    */
    int set_npoints(int npoints) {
      if (npoints==2) {
	cap=&deriv_eqi::deriv_vector3;
	cp=&deriv_eqi::derivp2;
      } else if (npoints==3) {
	cap=&deriv_eqi::deriv_vector3;
	cp=&deriv_eqi::derivp3;
      } else if (npoints==4) {
	cap=&deriv_eqi::deriv_vector4;
	cp=&deriv_eqi::derivp4;
      } else {
	cap=&deriv_eqi::deriv_vector5;
	cp=&deriv_eqi::derivp5;
      }
      if (npoints<=1 || npoints>5) {
	O2SCL_ERR("Invalid # of points in set_npoints(). Using default",
		      exc_einval);
      }
      return 0;
    }

    /** \brief Set the number of points to use for second derivatives
	(default 5)

	Acceptable values are 3-5 (see above).
    */
    int set_npoints2(int npoints) {
      if (npoints==3) {
	c2p=&deriv_eqi::deriv2p3;
      } else if (npoints==4) {
	c2p=&deriv_eqi::deriv2p4;
      } else {
	c2p=&deriv_eqi::deriv2p5;
      }
      if (npoints<=2 || npoints>5) {
	O2SCL_ERR("Invalid # of points in set_npoints2(). Using default",
		      exc_einval);
      }
      return 0;
    }

    /** \brief Calculate the first derivative of \c func w.r.t. x
    */
    virtual int deriv_err(double x, func_t &func,
			 double &dfdx, double &err) {
      double p=xoff/h;
      dfdx=(this->*cp)(x,p,func)/h;
      err=0.0;
      return success;
    }

    /** \brief Calculate the second derivative of \c func w.r.t. x
     */
    virtual int deriv2_err(double x, func_t &func,
			  double &dfdx, double &err) {
      double p=xoff/h;
      dfdx=(this->*c2p)(x,p,func)/h/h;
      err=0.0;
      return success;;
    }
    
    /** \brief Calculate the third derivative of \c func w.r.t. x
     */
    virtual int deriv3_err(double x, func_t &func,
			  double &dfdx, double &err) {
      double p=xoff/h;
      dfdx=(this->*c3p)(x,h,p,func)/h;
      err=0.0;
      return success;;
    }
    
    /** \brief Calculate the derivative at \c x given an array

	This calculates the derivative at \c x given a function
	specified in an array \c y of size \c nx with equally spaced
	abscissas. The first abscissa should be given as \c x0
	and the distance between adjacent abscissas should be
	given as \c dx. The value \c x need not be one of the
	abscissas (i.e. it can lie in between an interval). The 
	appropriate offset is calculated automatically.
    */
    double deriv_vector(double x, double x0, double dx, 
		      size_t nx, const vec_t &y) {
      size_t ix=(size_t)((x-x0)/dx);
      return (this->*cap)(x,x0,dx,nx,y,ix)/dx;
    }

    /** \brief Calculate the second derivative at \c x given an array

	This calculates the second derivative at \c x given a function
	specified in an array \c y of size \c nx with equally spaced
	abscissas. The first abscissa should be given as \c x0
	and the distance between adjacent abscissas should be
	given as \c dx. The value \c x need not be one of the
	abscissas (i.e. it can lie in between an interval). The 
	appropriate offset is calculated automatically.
    */
    double deriv2_vector(double x, double x0, double dx, 
		       size_t nx, const vec_t &y) 
    {
      size_t ix=(size_t)((x-x0)/dx);
      return (this->*c2ap)(x,x0,dx,nx,y,ix)/dx;
    }

    /** \brief Calculate the third derivative at \c x given an array

	This calculates the third derivative at \c x given a function
	specified in an array \c y of size \c nx with equally spaced
	abscissas. The first abscissa should be given as \c x0 and the
	distance between adjacent abscissas should be given as \c
	dx. The value \c x need not be one of the abscissas (i.e. it
	can lie in between an interval). The appropriate offset is
	calculated automatically.
    */
    double deriv3_vector(double x, double x0, double dx, 
		       size_t nx, const vec_t &y) 
    {
      size_t ix=(size_t)((x-x0)/dx);
      return (this->*c3ap)(x,x0,dx,nx,y,ix)/dx;
    }

    /** \brief Calculate the derivative of an entire array

	Right now this uses np=5.

	\todo generalize to other values of npoints.
    */
    int deriv_vector(size_t nv, double dx, const vec_t &y, 
		     vec_t &dydx) 
    {
      dydx[0]=(-25.0/12.0*y[0]+4.0*y[1]-3.0*y[2]+4.0/3.0*y[3]-0.25*y[4])/dx;
      dydx[1]=(-0.25*y[0]-5.0/6.0*y[1]+1.5*y[2]-0.5*y[3]+1.0/12.0*y[4])/dx;
      for(size_t i=2;i<nv-2;i++) {
	dydx[i]=(1.0/12.0*y[i-2]-2.0/3.0*y[i-1]+2.0/3.0*y[i+1]-
		 1.0/12.0*y[i+2])/dx;
      }
      dydx[nv-2]=(-1.0/12.0*y[nv-5]+0.5*y[nv-4]-1.5*y[nv-3]+
		  5.0/6.0*y[nv-2]+0.25*y[nv-1])/dx;
      dydx[nv-1]=(0.25*y[nv-5]-4.0/3.0*y[nv-4]+3.0*y[nv-3]-
		  4.0*y[nv-2]+25.0/12.0*y[nv-1])/dx;
      return 0;
    }

    /// Return string denoting type ("deriv_eqi")
    virtual const char *type() { return "deriv_eqi"; }

#ifndef DOXYGEN_NO_O2NS  

    protected:

    /** \brief Calculate the first derivative of \c func w.r.t. x and the
	uncertainty

	This function doesn't do anything, and isn't required for 
	this class since it computes higher-order derivatives directly.
    */
    virtual int deriv_err_int
    (double x, funct11 &func, double &dfdx, double &err) {
      return success;
    }
    
    /// Two-point first derivative
    double derivp2(double x, double p, func_t &func) {
      return (func(x+h)-func(x));
    }
      
    /// Three-point first derivative
    double derivp3(double x, double p, func_t &func)    {
      if (p==0.0) {
	return ((-0.5)*func(x-h)+(0.5)*func(x+h));
      }
      return ((p-0.5)*func(x-h)-2.0*p*func(x)+(p+0.5)*func(x+h));
    }

    /// Four-point first derivative
    double derivp4(double x, double p, 
		  func_t &func) {
      double p2=p*p;
      return (-(3.0*p2-6.0*p+2.0)/6.0*func(x-h)+
	      (1.5*p2-2.0*p-0.5)*func(x)-
	      (1.5*p2-p-1.0)*func(x+h)+
	      (3.0*p2-1.0)/6.0*func(x+2.0*h));
    }

    /// Five-point first derivative
    double derivp5(double x, double p, 
		  func_t &func) {
      double p2=p*p, p3=p*p*p;
      if (p==0.0) {
	return ((1.0)/12.0*func(x-2.0*h)-
		(4.0)/6.0*func(x-h)-
		(-4.0)/6.0*func(x+h)+
		(-1.0)/12.0*func(x+2.0*h));
      }
      return ((2.0*p3-3.0*p2-p+1.0)/12.0*func(x-2.0*h)-
	      (4.0*p3-3.0*p2-8.0*p+4.0)/6.0*func(x-h)+
	      (p3-2.5*p)*func(x)-
	      (4.0*p3+3.0*p2-8.0*p-4.0)/6.0*func(x+h)+
	      (2.0*p3+3.0*p2-p-1.0)/12.0*func(x+2.0*h));
    }

    
    /// Three-point first derivative for arrays
    double deriv_vector3(double x, double x0, double dx, size_t nx,
		       const vec_t &y, size_t ix) {
      double p;
      if (ix>0 && ix<nx-1) {
	p=x-(x0+ix*dx);
	return ((p-0.5)*y[ix-1]-2.0*p*y[ix]+(p+0.5)*y[ix+1]);
      } else if (ix==0) {
	p=x-(x0+dx);
	return ((p-0.5)*y[0]-2.0*p*y[1]+(p+0.5)*y[2]);
      } 
      p=x-(x0+(nx-2)*dx);
      return ((p-0.5)*y[nx-3]-2.0*p*y[nx-2]+(p+0.5)*y[nx-1]);
    }

    /// Four-point first derivative for arrays
    double deriv_vector4(double x, double x0, double dx, size_t nx,
		       const vec_t &y, size_t ix) {
      double p, p2;
      if (ix>0 && ix<nx-2) {
	p=x-(x0+ix*dx);
	p2=p*p;
	return (-(3.0*p2-6.0*p+2.0)/6.0*y[ix-1]+
		(1.5*p2-2.0*p-0.5)*y[ix]-
		(1.5*p2-p-1.0)*y[ix+1]+
		(3.0*p2-1.0)/6.0*y[ix+2]);
      } else if (ix==0) {
	p=x-(x0+dx);
	p2=p*p;
	return (-(3.0*p2-6.0*p+2.0)/6.0*y[0]+
		(1.5*p2-2.0*p-0.5)*y[1]-
		(1.5*p2-p-1.0)*y[2]+
		(3.0*p2-1.0)/6.0*y[3]);
      }
      p=x-(x0+(nx-3)*dx);
      p2=p*p;
      return (-(3.0*p2-6.0*p+2.0)/6.0*y[nx-4]+
	      (1.5*p2-2.0*p-0.5)*y[nx-3]-
	      (1.5*p2-p-1.0)*y[nx-2]+
	      (3.0*p2-1.0)/6.0*y[nx-1]);
    }

    /// Five-point first derivative for arrays
    double deriv_vector5(double x, double x0, 
		       double dx, size_t nx,
		       const vec_t &y, size_t ix) {
      double p, p2, p3;
      if (ix>1 && ix<nx-2) {
	p=x-(x0+ix*dx);
	p2=p*p, p3=p*p*p;
	return ((2.0*p3-3.0*p2-p+1.0)/12.0*y[ix-2]-
		(4.0*p3-3.0*p2-8.0*p+4.0)/6.0*y[ix-1]+
		(p3-2.5*p)*y[ix]-
		(4.0*p3+3.0*p2-8.0*p-4.0)/6.0*y[ix+1]+
		(2.0*p3+3.0*p2-p-1.0)/12.0*y[ix+2]);
      } else if (ix<=1) {
	p=x-(x0+2*dx);
	p2=p*p, p3=p*p*p;
	return ((2.0*p3-3.0*p2-p+1.0)/12.0*y[0]-
		(4.0*p3-3.0*p2-8.0*p+4.0)/6.0*y[1]+
		(p3-2.5*p)*y[2]-
		(4.0*p3+3.0*p2-8.0*p-4.0)/6.0*y[3]+
		(2.0*p3+3.0*p2-p-1.0)/12.0*y[4]);
      }
      p=x-(x0+(nx-3)*dx);
      p2=p*p, p3=p*p*p;
      return ((2.0*p3-3.0*p2-p+1.0)/12.0*y[nx-5]-
	      (4.0*p3-3.0*p2-8.0*p+4.0)/6.0*y[nx-4]+
	      (p3-2.5*p)*y[nx-3]-
	      (4.0*p3+3.0*p2-8.0*p-4.0)/6.0*y[nx-2]+
	      (2.0*p3+3.0*p2-p-1.0)/12.0*y[nx-1]);
    }

    /// Three-point second derivative
    double deriv2p3(double x, double p, func_t &func) {
      return (func(x-h)-2*func(x)+func(x+h));
    }

    /// Four-point second derivative
    double deriv2p4(double x, double p, func_t &func) {
      return ((1.0-2.0*p)*func(x-h)-(1.0-6.0*p)*func(x)
	      -(1.0-6.0*p)*func(x+h)+(1.0+2.0*p)*func(x+2.0*h))/2.0;
    }

    /// Five-point second derivative
    double deriv2p5(double x, double p, func_t &func) {
      return ((1.0-2.0*p)*(1.0-2.0*p)*func(x-2.0*h)
	      +(8.0*p-16.0*p*p)*func(x-h)
	      -(2.0-24.0*p*p)*func(x)
	      -(8.0*p+16.0*p*p)*func(x+h)
	      +(1.0+2.0*p)*(1.0+2.0*p)*func(x+2.0*h))/4.0;
    }

    /// Pointer to the first derivative function
    double (deriv_eqi::*cp)(double x, double p, 
			    func_t &func);

    /// Pointer to the first derivative for arrays function
    double (deriv_eqi::*cap)(double x, double x0, 
			     double dx, size_t nx,
			     const vec_t &y, size_t ix);

    /// Pointer to the second derivative function
    double (deriv_eqi::*c2p)(double x, double p,
			     func_t &func);

    /// Pointer to the second derivative for arrays function
    double (deriv_eqi::*c2ap)(double x, double x0, 
			      double dx, size_t nx,
			      const vec_t &y, size_t ix);

    /// Pointer to the third derivative function
    double (deriv_eqi::*c3p)(double x, double h, double p, 
			     func_t &func);

    /// Pointer to the third derivative for arrays function
    double (deriv_eqi::*c3ap)(double x, double x0, 
			      double dx, size_t nx,
			      const vec_t &y, size_t ix);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
