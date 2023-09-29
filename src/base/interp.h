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
/* interpolation/linear.c
 * interpolation/cspline.c
 * interpolation/akima.c
 * interpolation/steffen.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman
 * Copyright (C) 2014 Jean-François Caron
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
#ifndef O2SCL_INTERP_H
#define O2SCL_INTERP_H

/** \file interp.h
    \brief One-dimensional interpolation classes and interpolation types
*/

#include <iostream>
#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <o2scl/search_vec.h>
#include <o2scl/tridiag.h>
#include <o2scl/vector.h>
#include <o2scl/invert.h>

namespace o2scl {
  
  /** \brief Interpolation types in src/base/interp.h
   */
  enum {
    /// Linear
    itp_linear=1,
    /// Cubic spline for natural boundary conditions
    itp_cspline=2,
    /// Cubic spline for periodic boundary conditions
    itp_cspline_peri=3,
    /// Akima spline for natural boundary conditions
    itp_akima=4,
    /// Akima spline for periodic boundary conditions
    itp_akima_peri=5,
    /// Monotonicity-preserving interpolation (unfinished)
    itp_monotonic=6,
    /// Steffen's monotonic method
    itp_steffen=7,
    /// Nearest-neighbor lookup
    itp_nearest_neigh=8,
    /// Gaussian process with leave-one out cross validation (experimental)
    itp_gp_rbf_noise_loo_cv=9,
    /// Gaussian process with maximum log marginal likelihood (experimental)
    itp_gp_rbf_noise_max_lml=10
  };

  /** \brief Base low-level interpolation class [abstract base]

      \verbatim embed:rst
      See also the :ref:`Interpolation` section of the 
      O\ :sub:`2`\ scl User's guide. 
      \endverbatim
      
      To interpolate a set vectors \c x and \c y, call set() and then
      the interpolation functions eval(), deriv(), deriv2() and
      integ(). If the \c x and \c y vectors do not change, then you
      may call the interpolation functions multiple times in
      succession. These classes do not copy the user-specified vectors
      but store pointers to them. Thus, if the vector is changed
      without a new call to \ref interp_base::set(), the behavior of
      the interpolation functions is undefined.

      \comment
      AWS, 12/27/13: Copy constructors might be ill-advised for
      this class since we store pointers. For now, we don't 
      allow the user to use them.
      \endcomment
  */
  template<class vec_t, class vec2_t=vec_t> class interp_base {

#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /** \brief To perform binary searches
        
        This pointer is set to zero in the constructor and should be
        non-zero only if it has been allocated with \c new.
    */
    search_vec<const vec_t> svx;
    
    /// Independent vector
    const vec_t *px;
    
    /// Dependent vector
    const vec2_t *py;
    
    /// Vector size
    size_t sz;
    
    /** \brief An internal function to assist in computing the 
        integral for both the cspline and Akima types
    */
    double integ_eval(double ai, double bi, double ci, double di, double xi, 
                      double a, double b) const {
      
      double r1=a-xi;
      double r2=b-xi;
      double r12=r1+r2;
      double bterm=0.5*bi*r12;
      double cterm=ci*(r1*r1+r2*r2+r1*r2)/3.0;
      double dterm=0.25*di*r12*(r1*r1+r2*r2);
      
      return (b-a)*(ai+bterm+cterm+dterm);
    }
    
#endif
    
  public:
    
    interp_base() {
      sz=0;
    }
    
    virtual ~interp_base() {
    }

    /** \brief The minimum size of the vectors to interpolate between

        This variable must be set in the constructor of the children
        for access by the class user.
    */
    size_t min_size;

    /// Initialize interpolation routine
    virtual void set(size_t size, const vec_t &x, const vec2_t &y)=0;
    
    /// Give the value of the function \f$ y(x=x_0) \f$  
    virtual double eval(double x0) const=0;
    
    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double operator()(double x0) const {
      return eval(x0);
    }

    /// Give the value of the derivative \f$ y^{\prime}(x=x_0) \f$ .
    virtual double deriv(double x0) const=0;
    
    /** \brief Give the value of the second derivative 
        \f$ y^{\prime \prime}(x=x_0) \f$ .
    */
    virtual double deriv2(double x0) const=0;
    
    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(double a, double b) const=0;

    /// Return the type
    virtual const char *type() const=0;
 
#ifndef DOXYGEN_INTERNAL

  private:

    interp_base<vec_t,vec2_t>(const interp_base<vec_t,vec2_t> &);
    interp_base<vec_t,vec2_t>& operator=(const interp_base<vec_t,vec2_t>&);

#endif

  };
  
  /** \brief Linear interpolation (GSL) 

      \verbatim embed:rst
      See also the :ref:`Interpolation` section of the 
      O\ :sub:`2`\ scl User's guide. 
      \endverbatim

      Linear interpolation requires no calls to allocate() or free()
      as there is no internal storage required.
  */
  template<class vec_t, class vec2_t=vec_t> class interp_linear : 
  public interp_base<vec_t,vec2_t> {
    
  public:
    
    interp_linear() {
      this->min_size=2;
    }
    
    virtual ~interp_linear() {}
    
    /// Initialize interpolation routine
    virtual void set(size_t size, const vec_t &x, const vec2_t &y) {
      if (size<this->min_size) {
        O2SCL_ERR((((std::string)"Vector size, ")+szttos(size)+", is less"+
                   " than "+szttos(this->min_size)+" in interp_linear::"+
                   "set().").c_str(),exc_einval);
      }
      this->svx.set_vec(size,x);
      this->px=&x;
      this->py=&y;
      this->sz=size;
      return;
    }
    
    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
      
      double x_lo=(*this->px)[index];
      double x_hi=(*this->px)[index+1];
      double y_lo=(*this->py)[index];
      double y_hi=(*this->py)[index+1];
      double dx=x_hi-x_lo;
      
      return y_lo+(x0-x_lo)/dx*(y_hi-y_lo);
    }
    
    /// Give the value of the derivative \f$ y^{\prime}(x=x_0) \f$ .
    virtual double deriv(double x0) const {
      
      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
      
      double x_lo=(*this->px)[index];
      double x_hi=(*this->px)[index+1];
      double y_lo=(*this->py)[index];
      double y_hi=(*this->py)[index+1];
      double dx=x_hi-x_lo;
      double dy=y_hi-y_lo;

      return dy/dx;
    }

    /** \brief Give the value of the second derivative  
        \f$ y^{\prime \prime}(x=x_0) \f$ (always zero)
    */
    virtual double deriv2(double x0) const {
      return 0.0;
    }

    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(double a, double b) const {

      size_t cache=0;
      size_t i, index_a, index_b;
      
      bool flip=false;
      if (((*this->px)[0]<(*this->px)[this->sz-1] && a>b) ||
          ((*this->px)[0]>(*this->px)[this->sz-1] && a<b)) {
        double tmp=a;
        a=b;
        b=tmp;
        flip=true;
      }

      index_a=this->svx.find_const(a,cache);
      index_b=this->svx.find_const(b,cache);
      
      double result=0.0;
      for(i=index_a; i<=index_b; i++) {

        double x_lo=(*this->px)[i];
        double x_hi=(*this->px)[i+1];
        double y_lo=(*this->py)[i];
        double y_hi=(*this->py)[i+1];
        double dx=x_hi-x_lo;
        
        if(dx != 0.0) {

          if (i == index_a || i == index_b) {
            double x1=(i == index_a) ? a : x_lo;
            double x2=(i == index_b) ? b : x_hi;
            double D=(y_hi-y_lo)/dx;
            result += (x2-x1)*(y_lo+0.5*D*((x2-x_lo)+(x1-x_lo)));
          } else  {
            result += 0.5*dx*(y_lo+y_hi);
          }

        } else {
          std::string str=((std::string)"Interval of length zero ")+
            "between ("+o2scl::dtos(x_lo)+","+o2scl::dtos(y_lo)+
            ") and ("+o2scl::dtos(x_hi)+","+o2scl::dtos(y_hi)+
            " in interp_linear::integ().";
          O2SCL_ERR(str.c_str(),exc_einval);
        }
      }
      
      if (flip) result=-result;
      return result;
    }

    /// Return the type, \c "interp_linear".
    virtual const char *type() const { return "interp_linear"; }

#ifndef DOXYGEN_INTERNAL

  private:

    interp_linear<vec_t,vec2_t>(const interp_linear<vec_t,vec2_t> &);
    interp_linear<vec_t,vec2_t>& operator=(const interp_linear<vec_t,vec2_t>&);

#endif

  };
  
  /** \brief Nearest-neighbor interpolation

      Nearest interpolation requires no calls to allocate() or free()
      as there is no internal storage required.
  */
  template<class vec_t, class vec2_t=vec_t> class interp_nearest_neigh : 
  public interp_base<vec_t,vec2_t> {
    
  public:
    
    interp_nearest_neigh() {
      this->min_size=1;
    }
    
    virtual ~interp_nearest_neigh() {}
    
    /// Initialize interpolation routine
    virtual void set(size_t size, const vec_t &x, const vec2_t &y) {
      if (size<this->min_size) {
        O2SCL_ERR((((std::string)"Vector size, ")+szttos(size)+", is less"+
                   " than "+szttos(this->min_size)+" in interp_nearest_neigh::"+
                   "set().").c_str(),exc_einval);
      }
      this->svx.set_vec(size,x);
      this->px=&x;
      this->py=&y;
      this->sz=size;
      return;
    }
    
    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(double x0) const {

      size_t index=this->svx.ordered_lookup(x0);
      return (*this->py)[index];
    }
    
    /// Give the value of the derivative \f$ y^{\prime}(x=x_0) \f$ .
    virtual double deriv(double x0) const {
      return 0.0;
    }

    /** \brief Give the value of the second derivative  
        \f$ y^{\prime \prime}(x=x_0) \f$ (always zero)
    */
    virtual double deriv2(double x0) const {
      return 0.0;
    }

    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(double a, double b) const {
      return 0.0;
    }

    /// Return the type, \c "interp_nearest_neigh".
    virtual const char *type() const { return "interp_nearest_neigh"; }

#ifndef DOXYGEN_INTERNAL

  private:

    interp_nearest_neigh<vec_t,vec2_t>
      (const interp_nearest_neigh<vec_t,vec2_t> &);
    interp_nearest_neigh<vec_t,vec2_t>& operator=
      (const interp_nearest_neigh<vec_t,vec2_t>&);

#endif

  };
  
  /** \brief Cubic spline interpolation (GSL)
      
      \verbatim embed:rst
      See also the :ref:`Interpolation` section of the 
      O\ :sub:`2`\ scl User's guide. 
      \endverbatim

      By default, this class uses natural boundary conditions, where
      the second derivative vanishes at each end point. Extrapolation
      effectively assumes that the second derivative is linear outside
      of the endpoints. 
  */
  template<class vec_t, class vec2_t=vec_t> class interp_cspline : 
  public interp_base<vec_t,vec2_t> {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector_slice<ubvector> ubvector_slice;
    typedef boost::numeric::ublas::vector_range<ubvector> ubvector_range;
    typedef boost::numeric::ublas::slice slice;
    typedef boost::numeric::ublas::range range;
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /// \name Storage for cubic spline interpolation
    //@{
    ubvector c;
    ubvector g;
    ubvector diag;
    ubvector offdiag;
    //@}
    
    /// Memory for the tridiagonalization
    o2scl_linalg::ubvector_4_mem<double> p4m;
    
    /// Compute coefficients for cubic spline interpolation
    void coeff_calc(const ubvector &c_array, double dy, double dx, 
                    size_t index, double &b, double &c2, double &d) const {
      
      double c_i=c_array[index];
      double c_ip1=c_array[index+1];
      b=(dy/dx)-dx*(c_ip1+2.0*c_i)/3.0;
      c2=c_i;
      d=(c_ip1-c_i)/(3.0*dx);

      return;
    }

#endif
    
  public:

    /** \brief Create a base interpolation object with natural or
        periodic boundary conditions
    */
    interp_cspline() {
      this->min_size=3;
    }

    virtual ~interp_cspline() {
    }

    /** \brief Initialize interpolation routine
     */
    virtual void set(size_t size, const vec_t &xa, const vec2_t &ya) {

      if (size<this->min_size) {
        O2SCL_ERR((((std::string)"Vector size, ")+szttos(size)+", is less"+
                   " than "+szttos(this->min_size)+" in interp_cspline::"+
                   "set().").c_str(),exc_einval);
      }

      if (size!=this->sz) {
        c.resize(size);
        g.resize(size);
        diag.resize(size);
        offdiag.resize(size);
        p4m.resize(size);
      }
      
      this->px=&xa;
      this->py=&ya;
      this->sz=size;

      this->svx.set_vec(size,xa);

      /// Natural boundary conditions

      size_t i;
      size_t num_points=size;
      size_t max_index=num_points-1;
      size_t sys_size=max_index-1;

      c[0]=0.0;
      c[max_index]=0.0;

      for (i=0; i < sys_size; i++) {
        double h_i=xa[i+1]-xa[i];
        double h_ip1=xa[i+2]-xa[i+1];
        double ydiff_i=ya[i+1]-ya[i];
        double ydiff_ip1=ya[i+2]-ya[i+1];
        double g_i=(h_i != 0.0) ? 1.0/h_i : 0.0;
        double g_ip1=(h_ip1 != 0.0) ? 1.0/h_ip1 : 0.0;
        offdiag[i]=h_ip1;
        diag[i]=2.0*(h_ip1+h_i);
        g[i]=3.0*(ydiff_ip1*g_ip1-ydiff_i*g_i);
      }

      if (sys_size == 1) {

        c[1]=g[0]/diag[0];

        return;
      }

      ubvector_range cp1(c,range(1,c.size()));
      o2scl_linalg::solve_tridiag_sym<ubvector,ubvector,ubvector,
				      ubvector_range,
				      o2scl_linalg::ubvector_4_mem<double>,
				      ubvector,double>
        (diag,offdiag,g,cp1,sys_size,p4m);
      
      return;
    }

    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);

      double x_lo=(*this->px)[index];
      double x_hi=(*this->px)[index+1];
      double dx=x_hi-x_lo;

      double y_lo=(*this->py)[index];
      double y_hi=(*this->py)[index+1];
      double dy=y_hi-y_lo;
      double delx=x0-x_lo;
      double b_i, c_i, d_i; 

      coeff_calc(c,dy,dx,index,b_i,c_i,d_i);
      
      return y_lo+delx*(b_i+delx*(c_i+delx*d_i));
    }

    /// Give the value of the derivative \f$ y^{\prime}(x=x_0) \f$ .
    virtual double deriv(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
  
      double x_lo=(*this->px)[index];
      double x_hi=(*this->px)[index+1];
      double dx=x_hi-x_lo;

      double y_lo=(*this->py)[index];
      double y_hi=(*this->py)[index+1];
      double dy=y_hi-y_lo;
      double delx=x0-x_lo;
      double b_i, c_i, d_i; 

      coeff_calc(c,dy,dx,index,b_i,c_i,d_i);

      return b_i+delx*(2.0*c_i+3.0*d_i*delx);
    }

    /** \brief Give the value of the second derivative  
        \f$ y^{\prime \prime}(x=x_0) \f$ .
    */
    virtual double deriv2(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
  
      double x_lo=(*this->px)[index];
      double x_hi=(*this->px)[index+1];
      double dx=x_hi-x_lo;

      double y_lo=(*this->py)[index];
      double y_hi=(*this->py)[index+1];
      double dy=y_hi-y_lo;
      double delx=x0-x_lo;
      double b_i, c_i, d_i;

      coeff_calc(c,dy,dx,index,b_i,c_i,d_i);

      return 2.0*c_i+6.0*d_i*delx;
    }

    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(double a, double b) const {

      size_t i, index_a, index_b;
  
      bool flip=false;
      if (((*this->px)[0]<(*this->px)[this->sz-1] && a>b) ||
          ((*this->px)[0]>(*this->px)[this->sz-1] && a<b)) {
        double tmp=a;
        a=b;
        b=tmp;
        flip=true;
      }

      size_t cache=0;
      index_a=this->svx.find_const(a,cache);
      index_b=this->svx.find_const(b,cache);

      double result=0.0;
  
      for(i=index_a; i<=index_b; i++) {

        double x_lo=(*this->px)[i];
        double x_hi=(*this->px)[i+1];
        double y_lo=(*this->py)[i];
        double y_hi=(*this->py)[i+1];
        double dx=x_hi-x_lo;
        double dy=y_hi-y_lo;

        if(dx != 0.0) {
          double b_i, c_i, d_i; 
          coeff_calc(c,dy,dx,i,b_i,c_i,d_i);
          if (i == index_a || i == index_b) {
            double x1=(i == index_a) ? a : x_lo;
            double x2=(i == index_b) ? b : x_hi;
            result += this->integ_eval(y_lo,b_i,c_i,d_i,x_lo,x1,x2);
          } else {
            result += dx*(y_lo+dx*(0.5*b_i+
                                   dx*(c_i/3.0+0.25*d_i*dx)));
          }
        } else {
          std::string str=((std::string)"Interval of length zero ")+
            "between ("+o2scl::dtos(x_lo)+","+o2scl::dtos(y_lo)+
            ") and ("+o2scl::dtos(x_hi)+","+o2scl::dtos(y_hi)+
            " in interp_cspline::integ().";
          O2SCL_ERR(str.c_str(),exc_einval);
        }

      }
  
      if (flip) result*=-1.0;

      return result;
    }

    /// Return the type, \c "interp_cspline".
    virtual const char *type() const { return "interp_cspline"; }

#ifndef DOXYGEN_INTERNAL

  private:
  
  interp_cspline<vec_t,vec2_t>(const interp_cspline<vec_t,vec2_t> &);
  interp_cspline<vec_t,vec2_t>& operator=
  (const interp_cspline<vec_t,vec2_t>&);
  
#endif

  };
  
  /** \brief Cubic spline interpolation with periodic 
      boundary conditions (GSL)

      \verbatim embed:rst
      See also the :ref:`Interpolation` section of the 
      O\ :sub:`2`\ scl User's guide. 
      \endverbatim
  */
  template<class vec_t, class vec2_t=vec_t> class interp_cspline_peri : 
  public interp_cspline<vec_t,vec2_t> {

  protected:
    
    /// Memory for the tridiagonalization
    o2scl_linalg::ubvector_5_mem<double> p5m;
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector_slice<ubvector> ubvector_slice;
    typedef boost::numeric::ublas::vector_range<ubvector> ubvector_range;
    typedef boost::numeric::ublas::slice slice;
    typedef boost::numeric::ublas::range range;

  interp_cspline_peri() : interp_cspline<vec_t,vec2_t>() {
      this->min_size=2;
    }
      
    virtual ~interp_cspline_peri() {
    }

    /// Return the type, \c "interp_cspline_peri".
    virtual const char *type() const { return "interp_cspline_peri"; }

    /** \brief Initialize interpolation routine
     */
    virtual void set(size_t size, const vec_t &xa, const vec2_t &ya) {

      if (size<this->min_size) {
        O2SCL_ERR((((std::string)"Vector size, ")+szttos(size)+", is less"+
                   " than "+szttos(this->min_size)+" in interp_cspline"+
                   "_peri::set().").c_str(),exc_einval);
      }

      if (size!=this->sz) {
        this->c.resize(size);
        this->g.resize(size);
        this->diag.resize(size);
        this->offdiag.resize(size);
        p5m.resize(size);
      }

      this->px=&xa;
      this->py=&ya;
      this->sz=size;

      this->svx.set_vec(size,xa);

      /// Periodic boundary conditions
         
      size_t i;
      size_t num_points=size;
      // Engeln-Mullges+Uhlig "n" 
      size_t max_index=num_points-1;  
      // linear system is sys_size x sys_size 
      size_t sys_size=max_index;    
        
      if (sys_size == 2) {

        // solve 2x2 system 
          
        double h0=xa[1]-xa[0];
        double h1=xa[2]-xa[1];
        double h2=xa[3]-xa[2];
        double A=2.0*(h0+h1);
        double B=h0+h1;
        double gx[2];
        double det;
          
        gx[0]=3.0*((ya[2]-ya[1])/h1-(ya[1]-ya[0])/h0);
        gx[1]=3.0*((ya[1]-ya[2])/h2-(ya[2]-ya[1])/h1);
          
        det=3.0*(h0+h1)*(h0+h1);
        this->c[1]=( A*gx[0]-B*gx[1])/det;
        this->c[2]=(-B*gx[0]+A*gx[1])/det;
        this->c[0]=this->c[2];
          
        return;

      } else {
          
        for (i=0; i < sys_size-1; i++) {
          double h_i=xa[i+1]-xa[i];
          double h_ip1=xa[i+2]-xa[i+1];
          double ydiff_i=ya[i+1]-ya[i];
          double ydiff_ip1=ya[i+2]-ya[i+1];
          double g_i=(h_i != 0.0) ? 1.0/h_i : 0.0;
          double g_ip1=(h_ip1 != 0.0) ? 1.0/h_ip1 : 0.0;
          this->offdiag[i]=h_ip1;
          this->diag[i]=2.0*(h_ip1+h_i);
          this->g[i]=3.0*(ydiff_ip1*g_ip1-ydiff_i*g_i);
        }
          
        i=sys_size-1;
        {
          double h_i=xa[i+1]-xa[i];
          double h_ip1=xa[1]-xa[0];
          double ydiff_i=ya[i+1]-ya[i];
          double ydiff_ip1=ya[1]-ya[0];
          double g_i=(h_i != 0.0) ? 1.0/h_i : 0.0;
          double g_ip1=(h_ip1 != 0.0) ? 1.0/h_ip1 : 0.0;
          this->offdiag[i]=h_ip1;
          this->diag[i]=2.0*(h_ip1+h_i);
          this->g[i]=3.0*(ydiff_ip1*g_ip1-ydiff_i*g_i);
        }
        
        ubvector_range cp1(this->c,range(1,this->c.size()));
        o2scl_linalg::solve_cyc_tridiag_sym
	  <ubvector,ubvector,ubvector,ubvector_range,
	   o2scl_linalg::ubvector_5_mem<double>,ubvector,
	   double>
          (this->diag,this->offdiag,this->g,cp1,sys_size,p5m);
        this->c[0]=this->c[max_index];

        return;
      }

    }

#ifndef DOXYGEN_INTERNAL

  private:

    interp_cspline_peri<vec_t,vec2_t>
  (const interp_cspline_peri<vec_t,vec2_t> &);
    interp_cspline_peri<vec_t,vec2_t>& operator=
      (const interp_cspline_peri<vec_t,vec2_t>&);

#endif

  };

  /** \brief Akima spline interpolation (GSL)
      
      \verbatim embed:rst
      See also the :ref:`Interpolation` section of the 
      O\ :sub:`2`\ scl User's guide. 
      \endverbatim

      This class uses natural boundary conditions, where the second
      derivative vanishes at each end point. Extrapolation effectively
      assumes that the second derivative is linear outside of the
      endpoints. Use \ref interp_akima_peri for perodic boundary
      conditions.
  */
  template<class vec_t, class vec2_t=vec_t> class interp_akima : 
  public interp_base<vec_t,vec2_t> {

  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector_slice<ubvector> ubvector_slice;
    typedef boost::numeric::ublas::vector_range<ubvector> ubvector_range;
    typedef boost::numeric::ublas::slice slice;
    typedef boost::numeric::ublas::range range;

#ifndef DOXYGEN_INTERNAL
    
  protected:

    /// \name Storage for Akima spline interpolation
    //@{
    ubvector b;
    ubvector c;
    ubvector d;
    ubvector um;
    //@}

    /// For initializing the interpolation
    void akima_calc(const vec_t &x_array, size_t size, 
                    ubvector &umx) {
      
      for(size_t i=0;i<this->sz-1;i++) {
        
        double NE=fabs(umx[3+i]-umx[2+i])+fabs(umx[1+i]-umx[i]);
        
        if (NE == 0.0) {
          b[i]=umx[2+i];
          c[i]=0.0;
          d[i]=0.0;
        } else {
          double h_i=(*this->px)[i+1]-(*this->px)[i];
          double NE_next=fabs(umx[4+i]-umx[3+i])+
            fabs(umx[2+i]-umx[1+i]);
          double alpha_i=fabs(umx[1+i]-umx[i])/NE;
          double alpha_ip1;
          double tL_ip1;
          if (NE_next == 0.0) {
            tL_ip1=umx[2+i];
          } else {
            alpha_ip1=fabs(umx[2+i]-umx[1+i])/NE_next;
            tL_ip1=(1.0-alpha_ip1)*umx[2+i]+alpha_ip1*umx[3+i];
          }
          b[i]=(1.0-alpha_i)*umx[1+i]+alpha_i*umx[2+i];
          c[i]=(3.0*umx[2+i]-2.0*b[i]-tL_ip1)/h_i;
          d[i]=(b[i]+tL_ip1-2.0*umx[2+i])/(h_i*h_i);
        }
      }
    }
    
#endif
    
  public:

    /** \brief Create a base interpolation object with or without
        periodic boundary conditions
    */
    interp_akima() {
      this->min_size=5;
    }

    virtual ~interp_akima() {
    }
    
    /** \brief Initialize interpolation routine
     */
    virtual void set(size_t size, const vec_t &xa, const vec2_t &ya) {
      
      if (size<this->min_size) {
        O2SCL_ERR((((std::string)"Vector size, ")+szttos(size)+", is less"+
                   " than "+szttos(this->min_size)+" in interp_akima::"+
                   "set().").c_str(),exc_einval);
      }

      if (size!=this->sz) {
        b.resize(size);
        c.resize(size);
        d.resize(size);
        um.resize(size+4);
      }

      this->px=&xa;
      this->py=&ya;
      this->sz=size;

      this->svx.set_vec(size,xa);

      // Non-periodic boundary conditions

      ubvector_range m(um,range(2,um.size()));
      size_t i;
      for (i=0;i<=size-2;i++) {
        m[i]=(ya[i+1]-ya[i])/(xa[i+1]-xa[i]);
      }
        
      um[0]=3.0*m[0]-2.0*m[1];
      um[1]=2.0*m[0]-m[1];
      m[this->sz-1]=2.0*m[size-2]-m[size-3];
      m[size]=3.0*m[size-2]-2.0*m[size-3];
        
      akima_calc(xa,size,um);

      return;
    }
          
    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(double x0) const {
      
      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
  
      double x_lo=(*this->px)[index];
      double delx=x0-x_lo;
      double bb=b[index];
      double cc=c[index];
      double dd=d[index];

      return (*this->py)[index]+delx*(bb+delx*(cc+dd*delx));
    }

    /// Give the value of the derivative \f$ y^{\prime}(x=x_0) \f$ .
    virtual double deriv(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);

      double x_lo=(*this->px)[index];
      double delx=x0-x_lo;
      double bb=b[index];
      double cc=c[index];
      double dd=d[index];

      return bb+delx*(2.0*cc+3.0*dd*delx);
    }

    /** \brief Give the value of the second derivative  
        \f$ y^{\prime \prime}(x=x_0) \f$ .
    */
    virtual double deriv2(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
  
      double x_lo=(*this->px)[index];
      double delx=x0-x_lo;
      double cc=c[index];
      double dd=d[index];

      return 2.0*cc+6.0*dd*delx;
    }

    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(double aa, double bb) const {

      size_t i, index_a, index_b;
  
      bool flip=false;
      if (((*this->px)[0]<(*this->px)[this->sz-1] && aa>bb) ||
          ((*this->px)[0]>(*this->px)[this->sz-1] && aa<bb)) {
        double tmp=aa;
        aa=bb;
        bb=tmp;
        flip=true;
      }

      size_t cache=0;
      index_a=this->svx.find_const(aa,cache);
      index_b=this->svx.find_const(bb,cache);

      double result=0.0;
  
      for(i=index_a; i<=index_b; i++) {

        double x_lo=(*this->px)[i];
        double x_hi=(*this->px)[i+1];
        double y_lo=(*this->py)[i];
        double dx=x_hi-x_lo;

        if (dx != 0.0) {
          
          if (i==index_a || i==index_b) {
            double x1=(i==index_a) ? aa : x_lo;
            double x2=(i==index_b) ? bb : x_hi;
            result += this->integ_eval(y_lo,b[i],c[i],d[i],x_lo,x1,x2);
          } else {
            result+=dx*(y_lo+dx*(0.5*b[i]+dx*(c[i]/3.0+0.25*d[i]*dx)));
          }

        } else {
          double y_hi=(*this->py)[i+1];
          std::string str=((std::string)"Interval of length zero ")+
            "between ("+o2scl::dtos(x_lo)+","+o2scl::dtos(y_lo)+
            ") and ("+o2scl::dtos(x_hi)+","+o2scl::dtos(y_hi)+
            " in interp_akima::integ().";
          O2SCL_ERR(str.c_str(),exc_einval);
        }
      }
  
      if (flip) result*=-1.0;

      return result;
    }

    /// Return the type, \c "interp_akima".
    virtual const char *type() const { return "interp_akima"; }

#ifndef DOXYGEN_INTERNAL

  private:

    interp_akima<vec_t,vec2_t>(const interp_akima<vec_t,vec2_t> &);
    interp_akima<vec_t,vec2_t>& operator=(const interp_akima<vec_t,vec2_t>&);

#endif

  };
  
  /** \brief Akima spline interpolation with periodic  
      boundary conditions (GSL)

      \verbatim embed:rst
      See also the :ref:`Interpolation` section of the 
      O\ :sub:`2`\ scl User's guide. 
      \endverbatim
  */
  template<class vec_t, class vec2_t=vec_t> class interp_akima_peri : 
  public interp_akima<vec_t,vec2_t> {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector_slice<ubvector> ubvector_slice;
    typedef boost::numeric::ublas::vector_range<ubvector> ubvector_range;
    typedef boost::numeric::ublas::slice slice;
    typedef boost::numeric::ublas::range range;

  public:
    
  interp_akima_peri() : interp_akima<vec_t,vec2_t>() {
    }
    
    virtual ~interp_akima_peri() {
    }
    
    /// Return the type, \c "interp_akima_peri".
    virtual const char *type() const { return "interp_akima_peri"; }

    /// Initialize interpolation routine
    virtual void set(size_t size, const vec_t &xa, const vec2_t &ya) {
      
      if (size<this->min_size) {
        O2SCL_ERR((((std::string)"Vector size, ")+szttos(size)+", is less"+
                   " than "+szttos(this->min_size)+" in interp_akima"+
                   "_peri::set().").c_str(),exc_einval);
      }
      
      if (size!=this->sz) {
        this->b.resize(size);
        this->c.resize(size);
        this->d.resize(size);
        this->um.resize(size+4);
      }

      this->px=&xa;
      this->py=&ya;
      this->sz=size;

      this->svx.set_vec(size,xa);

      // Periodic boundary conditions
      
      ubvector_range m(this->um,range(2,this->um.size()));

      // Form the required set of divided differences
      for (size_t i=0;i<=size-2;i++) {
        m[i]=(ya[i+1]-ya[i])/(xa[i+1]-xa[i]);
      }
      
      this->um[0]=m[this->sz-3];
      this->um[1]=m[this->sz-2];
      m[this->sz-1]=m[0];
      m[this->sz]=m[1];
      
      this->akima_calc(xa,size,this->um);
      
      return;
    }
      
#ifndef DOXYGEN_INTERNAL

  private:

    interp_akima_peri<vec_t,vec2_t>(const interp_akima_peri<vec_t,vec2_t> &);
    interp_akima_peri<vec_t,vec2_t>& operator=
      (const interp_akima_peri<vec_t,vec2_t>&);

#endif

  };

  /** \brief Steffen's monotonicity-preserving interpolation

      \verbatim embed:rst
      Adapted from the GSL version by J.-F. Caron which
      was based on [Steffen90]_.
      \endverbatim
   */
  template<class vec_t, class vec2_t=vec_t> class interp_steffen : 
  public interp_base<vec_t,vec2_t> {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector_slice<ubvector> ubvector_slice;
    typedef boost::numeric::ublas::vector_range<ubvector> ubvector_range;
    typedef boost::numeric::ublas::slice slice;
    typedef boost::numeric::ublas::range range;
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /// \name Storage for cubic spline interpolation
    //@{
    ubvector a;
    ubvector b;
    ubvector c;
    ubvector d;
    ubvector y_prime;
    //@}
    
    /** \brief Flip the sign of x if x and y have different signs
     */
    double copysign(const double x, const double y) {
      if ((x < 0 && y > 0) || (x > 0 && y < 0)) {
        return -x;
      }
      return x;
    }

#endif
    
  public:

    /** \brief Create a base interpolation object */
    interp_steffen() {
      this->min_size=3;
    }

    virtual ~interp_steffen() {
    }

    /** \brief Initialize interpolation routine
     */
    virtual void set(size_t size, const vec_t &xa, const vec2_t &ya) {
      
      if (size<this->min_size) {
        O2SCL_ERR((((std::string)"Vector size, ")+szttos(size)+", is less"+
                   " than "+szttos(this->min_size)+" in interp_steffen::"+
                   "set().").c_str(),exc_einval);
      }
      
      if (size!=this->sz) {
        a.resize(size);
        b.resize(size);
        c.resize(size);
        d.resize(size);
        y_prime.resize(size);
      }
      
      this->px=&xa;
      this->py=&ya;
      this->sz=size;
      
      this->svx.set_vec(size,xa);
      
      /*
       * first assign the interval and slopes for the left boundary.
       * We use the "simplest possibility" method described in the paper
       * in section 2.2
       */
      double h0=(xa[1]-xa[0]);
      double s0=(ya[1]-ya[0]) / h0;
      
      y_prime[0]=s0;
      
      /* Now we calculate all the necessary s, h, p, and y' variables 
         from 1 to N-2 (0 to size-2 inclusive) */
      for (size_t i=1; i < (size-1); i++) {
        
        double pi;
        
        /* equation 6 in the paper */
        double hi=(xa[i+1]-xa[i]);
        double him1=(xa[i]-xa[i-1]);
        
        /* equation 7 in the paper */
        double si=(ya[i+1]-ya[i]) / hi;
        double sim1=(ya[i]-ya[i-1]) / him1;
        
        /* equation 8 in the paper */
        pi=(sim1*hi + si*him1) / (him1 + hi);
        
        /* This is a C equivalent of the FORTRAN statement below eqn 11 */
        y_prime[i]=(copysign(1.0,sim1)+copysign(1.0,si))*
          std::min(fabs(sim1),std::min(fabs(si),0.5*fabs(pi))); 

      }

      /*
       * we also need y' for the rightmost boundary; we use the
       * "simplest possibility" method described in the paper in
       * section 2.2
       *
       * y'=s_{n-1}
       */
      y_prime[size-1]=(ya[size-1]-ya[size-2])/
        (xa[size-1]-xa[size-2]);
      
      /* Now we can calculate all the coefficients for the whole range. */
      for (size_t i=0; i < (size-1); i++) {
        double hi=(xa[i+1]-xa[i]);
        double si=(ya[i+1]-ya[i]) / hi;
        
        /* These are from equations 2-5 in the paper. */
        a[i]=(y_prime[i] + y_prime[i+1]-2*si) / hi / hi;
        b[i]=(3*si-2*y_prime[i]-y_prime[i+1]) / hi;
        c[i]=y_prime[i];
        d[i]=ya[i];
      }
      
      return;
    }
    
    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
      double x_lo=(*this->px)[index];
      double delx=x0-x_lo;
      
      /* Use Horner's scheme for efficient evaluation of polynomials. */
      double y = d[index]+delx*(c[index]+delx*(b[index]+delx*a[index]));
      
      return y;
    }

    /// Give the value of the derivative \f$ y^{\prime}(x=x_0) \f$ .
    virtual double deriv(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
      double x_lo=(*this->px)[index];
      double delx=x0-x_lo;

      return c[index]+delx*(2.0*b[index]+delx*3.0*a[index]);
    }

    /** \brief Give the value of the second derivative  
        \f$ y^{\prime \prime}(x=x_0) \f$ .
    */
    virtual double deriv2(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
      double x_lo=(*this->px)[index];
      double delx=x0-x_lo;

      return 2.0*b[index]+delx*6.0*a[index];
    }

    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(double al, double bl) const {

      size_t i, index_a, index_b;
  
      bool flip=false;
      if (((*this->px)[0]<(*this->px)[this->sz-1] && al>bl) ||
          ((*this->px)[0]>(*this->px)[this->sz-1] && al<bl)) {
        double tmp=al;
        al=bl;
        bl=tmp;
        flip=true;
      }

      size_t cache=0;
      index_a=this->svx.find_const(al,cache);
      index_b=this->svx.find_const(bl,cache);

      double result=0.0;
  
      for(i=index_a; i<=index_b; i++) {

        double x_lo=(*this->px)[i];
        double x_hi=(*this->px)[i+1];
        double y_lo=(*this->py)[i];
        double y_hi=(*this->py)[i+1];
        double dx=x_hi-x_lo;

        if(dx != 0.0) {

          double x1=(i == index_a) ? al-x_lo : 0.0;
          double x2=(i == index_b) ? bl-x_lo : x_hi-x_lo;
          result += (1.0/4.0)*a[i]*(x2*x2*x2*x2-x1*x1*x1*x1)+
            (1.0/3.0)*b[i]*(x2*x2*x2-x1*x1*x1)+
            (1.0/2.0)*c[i]*(x2*x2-x1*x1)+d[i]*(x2-x1);

        } else {
          std::string str=((std::string)"Interval of length zero ")+
            "between ("+o2scl::dtos(x_lo)+","+o2scl::dtos(y_lo)+
            ") and ("+o2scl::dtos(x_hi)+","+o2scl::dtos(y_hi)+
            " in interp_steffen::integ().";
          O2SCL_ERR(str.c_str(),exc_einval);
        }

      }
  
      if (flip) result*=-1.0;

      return result;
    }

    /// Return the type, \c "interp_steffen".
    virtual const char *type() const { return "interp_steffen"; }

#ifndef DOXYGEN_INTERNAL

  private:
  
  interp_steffen<vec_t,vec2_t>(const interp_steffen<vec_t,vec2_t> &);
  interp_steffen<vec_t,vec2_t>& operator=
  (const interp_steffen<vec_t,vec2_t>&);
  
#endif

  };
  
  /** \brief Monotonicity-preserving interpolation

      \warning This class is experimental. Integrals don't work yet.

      \verbatim embed:rst
      This class uses a method based on cubic Hermite interpolation,
      modifying the slopes to guarantee monotonicity. In the
      notation of [Fritsch80]_, if 
      \endverbatim
      \f[
      \alpha_i^2+\beta_i^2 \geq 9 \, ,
      \f]
      then \f$ \alpha \f$ and \f$ \beta \f$ are decreased by 
      the factor \f$ \tau \f$ as described at the end of 
      section 4. 

      \note The results of the interpolation will only be monotonic in
      the regions where the original data set is also monotonic. Also,
      this interpolation routine does not enforce strict monotonicity,
      and the results of the interpolation will be flat where the data
      is also flat.

      \verbatim embed:rst
      Based on [Fritsch80]_.
      \endverbatim

      \future Convert into fewer loops over the data
  */
  template<class vec_t, class vec2_t=vec_t> class interp_monotonic :
  public interp_base<vec_t,vec2_t> {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    
  protected:
    
    /// Slopes
    ubvector m;
    /// Finite differences
    ubvector Delta;
    /// Ratio 
    ubvector alpha;
    /// Staggered ratio
    ubvector beta;
    
  public:

    interp_monotonic() {
      this->min_size=2;
    }
  
    virtual ~interp_monotonic() {
    }
    
    /// Initialize interpolation routine
    virtual void set(size_t size, const vec_t &x, const vec2_t &y) {

      // Verify size
      if (size<this->min_size) {
        O2SCL_ERR((((std::string)"Vector size, ")+szttos(size)+", is less"+
                   " than "+szttos(this->min_size)+" in interp_monotonic::"+
                   "set().").c_str(),exc_einval);
      }
      
      // Setup search_vec object
      this->svx.set_vec(size,x);

      // Resize internal vectors
      if (this->sz!=size) {
        m.resize(size);
        Delta.resize(size-1);
        alpha.resize(size-1);
        beta.resize(size-1);
      }
      
      // Copy pointers
      this->px=&x;
      this->py=&y;
      this->sz=size;

      // Create the interpolation arrays in this sub-interval
      
      // Compute Delta and m
      for(size_t i=0;i<size-1;i++) {
        Delta[i]=(y[i+1]-y[i])/(x[i+1]-x[i]);
        if (i>0) {
          m[i]=(Delta[i]+Delta[i-1])/2.0;
        }
      }
      m[0]=Delta[0];
      m[size-1]=Delta[size-2];
      
      // Check to see if the data is flat anywhere
      for(size_t i=0;i<size-1;i++) {
        if (y[i]==y[i+1]) {
          m[i]=0.0;
          m[i+1]=0.0;
        }
      }
      
      // Compute alpha and beta
      for(size_t i=0;i<size-1;i++) {
        alpha[i]=m[i]/Delta[i];
        beta[i]=m[i+1]/Delta[i];
      }

      // Constrain m to ensure monotonicity
      for(size_t i=0;i<size-1;i++) {
        double norm2=alpha[i]*alpha[i]+beta[i]*beta[i];
        if (norm2>9.0) {
          double tau=3.0/sqrt(norm2);
          m[i]=tau*alpha[i]*Delta[i];
          m[i+1]=tau*beta[i]*Delta[i];
        }
      }
            
      return;
    }
    
    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
      
      double x_lo=(*this->px)[index];
      double x_hi=(*this->px)[index+1];
      double y_lo=(*this->py)[index];
      double y_hi=(*this->py)[index+1];
      double h=x_hi-x_lo;
      double t=(x0-x_lo)/h;
      double t2=t*t, t3=t2*t;

      double h00=2.0*t3-3.0*t2+1.0;
      double h10=t3-2.0*t2+t;
      double h01=-2.0*t3+3.0*t2;
      double h11=t3-t2;
      double interp=y_lo*h00+h*m[index]*h10+y_hi*h01+h*m[index+1]*h11;
      
      return interp;
    }
    
    /// Give the value of the derivative \f$ y^{\prime}(x=x_0) \f$ .
    virtual double deriv(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
      
      double x_lo=(*this->px)[index];
      double x_hi=(*this->px)[index+1];
      double y_lo=(*this->py)[index];
      double y_hi=(*this->py)[index+1];
      double h=x_hi-x_lo;
      double t=(x0-x_lo)/h;
      double t2=t*t;

      double dh00=6.0*t2-6.0*t;
      double dh10=3.0*t2-4.0*t+1.0;
      double dh01=-6.0*t2+6.0*t;
      double dh11=3.0*t2-2.0*t;
      double deriv=(y_lo*dh00+h*m[index]*dh10+y_hi*dh01+
                    h*m[index+1]*dh11)/h;

      return deriv;
    }

    /** \brief Give the value of the second derivative  
        \f$ y^{\prime \prime}(x=x_0) \f$
    */
    virtual double deriv2(double x0) const {

      size_t cache=0;
      size_t index=this->svx.find_const(x0,cache);
      
      double x_lo=(*this->px)[index];
      double x_hi=(*this->px)[index+1];
      double y_lo=(*this->py)[index];
      double y_hi=(*this->py)[index+1];
      double h=x_hi-x_lo;
      double t=(x0-x_lo)/h;

      double ddh00=12.0*t-6.0;
      double ddh10=6.0*t-4.0;
      double ddh01=-12.0*t+6.0;
      double ddh11=6.0*t-2.0;
      double deriv2=(y_lo*ddh00+h*m[index]*ddh10+y_hi*ddh01+
                     h*m[index+1]*ddh11)/h/h;

      return deriv2;
    }

    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(double a, double b) const {
      
      size_t i, index_a, index_b;
      
      bool flip=false;
      if (((*this->px)[0]<(*this->px)[this->sz-1] && a>b) ||
          ((*this->px)[0]>(*this->px)[this->sz-1] && a<b)) {
        double tmp=a;
        a=b;
        b=tmp;
        flip=true;
      }

      size_t cache=0;
      index_a=this->svx.find_const(a,cache);
      index_b=this->svx.find_const(b,cache);

      double result=0.0;
  
      for(i=index_a; i<=index_b; i++) {

        double x_hi=(*this->px)[i+1];
        double x_lo=(*this->px)[i];
        double y_lo=(*this->py)[i];
        double y_hi=(*this->py)[i+1];
        double h=x_hi-x_lo;
        
        if (h != 0.0) {
          
          if (i == index_a) {
            x_lo=a;
          }
          if (i == index_b) {
            x_hi=b;
          }

          double t=(x_hi-x_lo)/h;
          double t2=t*t, t3=t2*t, t4=t3*t;
          
          double ih00=t4/2.0-t3+t;
          double ih10=t4/4.0-2.0*t3/3.0+t2/2.0;
          double ih01=-t4/2.0+t3;
          double ih11=t4/4.0-t3/3.0;
          double intres=h*(y_lo*ih00+h*m[i]*ih10+y_hi*ih01+
                           h*m[i+1]*ih11);
          result+=intres;

        } else {
          std::string str=((std::string)"Interval of length zero ")+
            "between ("+o2scl::dtos(x_lo)+","+o2scl::dtos(y_lo)+
            ") and ("+o2scl::dtos(x_hi)+","+o2scl::dtos(y_hi)+
            " in interp_monotonic::integ().";
          O2SCL_ERR(str.c_str(),exc_einval);
        }

      }
  
      if (flip) result*=-1.0;

      return result;
    }

    /// Return the type, \c "interp_monotonic".
    virtual const char *type() const { return "interp_monotonic"; }

#ifndef DOXYGEN_INTERNAL

  private:

    interp_monotonic<vec_t,vec2_t>(const interp_monotonic<vec_t,vec2_t> &);
    interp_monotonic<vec_t,vec2_t>& operator=
      (const interp_monotonic<vec_t,vec2_t>&);

#endif

  };

}

#endif
