/*
  -------------------------------------------------------------------
  
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

  -------------------------------------------------------------------
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
#ifndef O2SCL_INTERP_VEC_H
#define O2SCL_INTERP_VEC_H

/** \file interp.h
    \brief One-dimensional interpolation classes and interpolation types
*/

#include <o2scl/interp.h>
#include <o2scl/interp_krige.h>

#ifdef O2SCL_EIGEN
#include <eigen3/Eigen/Dense>
#endif

namespace o2scl {

  /** \brief Interpolation class for pre-specified vector
      
      \verbatim embed:rst
      See also the :ref:`Interpolation` section of the 
      O\ :sub:`2`\ scl User's guide. 
      \endverbatim

      This interpolation class is intended to be sufficiently general
      to handle any vector type. Interpolation of ublas vector like
      objects is performed with the default template parameters.

      This class does not double check the vector to ensure that 
      all of the intervals for the abcissa are all increasing or
      all decreasing. 

      The type of interpolation to be performed can be specified 
      using the set_type() function. The default is cubic splines
      with natural boundary conditions. 

      \future Make version which copies vectors rather than storing
      pointers might be better and then has copy constructors.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class vec2_t=vec_t> class interp_vec : 
    public interp_base<vec_t,vec2_t> {
      
#ifndef DOXYGEN_INTERNAL
      
  protected:
    
    /// Base interpolation object
    interp_base<vec_t,vec2_t> *itp;
    
    /// Interpolation type
    size_t itype;

    /// Covariance function for Gaussian process interpolation
    covar_funct *cf;

    /// Parameter list for Gaussian process interpolation
    std::vector<std::vector<double>> param_list;
    
#endif
      
  public:

    /** \brief Create an interpolation object with interpolation type
        \c interp_type
    */
    interp_vec(size_t interp_type=itp_cspline) {
      itp=0;
      cf=0;
      itype=interp_type;
    }

    /** \brief Create an interpolation object with interpolation type
        \c interp_type based on the first \c n entries of vectors
        \c x and \c y
    */
    interp_vec(size_t n, const vec_t &x, 
               const vec2_t &y, size_t interp_type=itp_cspline) {
      cf=0;
      itp=0;
      set(n,x,y,interp_type);
    }

    /// Free the memory associated with the interpolation object
    virtual ~interp_vec() {
      clear();
    }

    /// \name Set methods
    //@{
    /** \brief Set the interpolation type
     */
    void set_type(size_t interp_type=itp_cspline) {
      itype=interp_type;
      return;
    }

    /** \brief Modify the interpolation object to operate on the first
        \c n entries of vectors \c x and \c y
    */
    void set(size_t n, const vec_t &x, const vec2_t &y) {
      set(n,x,y,itype);
      return;
    }

    /** \brief Set a new vector to interpolate
     */
    void set(size_t n, const vec_t &x, 
             const vec2_t &y, size_t interp_type) {

      if (n==0) {
        O2SCL_ERR("Cannot give vector of length 0 to interp_vec::set().",
                  o2scl::exc_einval);
      }
    
      if (x[0]==x[n-1]) {
        O2SCL_ERR((((std::string)"Vector endpoints equal (value=")+
                   o2scl::dtos(x[0])+") in interp_vec()::"+
                   "interp_vec().").c_str(),exc_einval);
      }
      
      clear();
    
      if (interp_type==itp_linear) {
        itp=new interp_linear<vec_t,vec2_t>;
      } else if (interp_type==itp_cspline) {
        itp=new interp_cspline<vec_t,vec2_t>;
      } else if (interp_type==itp_cspline_peri) {
        itp=new interp_cspline_peri<vec_t,vec2_t>;
      } else if (interp_type==itp_akima) {
        itp=new interp_akima<vec_t,vec2_t>;
      } else if (interp_type==itp_akima_peri) {
        itp=new interp_akima_peri<vec_t,vec2_t>;
      } else if (interp_type==itp_monotonic) {
        itp=new interp_monotonic<vec_t,vec2_t>;
      } else if (interp_type==itp_steffen) {
        itp=new interp_steffen<vec_t,vec2_t>;
      } else if (interp_type==itp_nearest_neigh) {
        itp=new interp_nearest_neigh<vec_t,vec2_t>;
        
      } else if (interp_type==itp_gp_rbf_noise_loo_cv ||
                 interp_type==itp_gp_rbf_noise_max_lml) {

#ifdef O2SCL_EIGEN
        interp_krige_optim<vec_t,vec2_t,covar_funct_rbf_noise,
                           Eigen::MatrixXd,
                           o2scl_linalg::matrix_invert_det_eigen<>> *ikon=
          new interp_krige_optim<vec_t,vec2_t,covar_funct_rbf_noise,
                                 Eigen::MatrixXd,
                                 o2scl_linalg::matrix_invert_det_eigen<>>;
#else        
        interp_krige_optim<vec_t,vec2_t,covar_funct_rbf_noise> *ikon=
          new interp_krige_optim<vec_t,vec2_t,covar_funct_rbf_noise>;
#endif
        
        itp=ikon;
        covar_funct_rbf_noise *cfrn=new covar_funct_rbf_noise;
        cf=cfrn;

        if (interp_type==itp_gp_rbf_noise_loo_cv) {
          ikon->mode=ikon->mode_loo_cv;
        } else {
          ikon->mode=ikon->mode_max_lml;
        }

        std::vector<double> diff;
        o2scl::vector_diffs(n,x,diff);
        double min_len=vector_min_value<std::vector<double>,double>
          (diff.size(),diff)/10.0;
        double max_len=vector_max_value<std::vector<double>,double>
          (diff.size(),diff)*10.0;
        std::vector<double> len_list;
        for(size_t i=0;i<10;i++) {
          len_list.push_back(min_len*pow(max_len/min_len,
                                         ((double)i)/9.0));
        }
        param_list.clear();
        param_list.push_back(len_list);
        double min_y=vector_min_value<vec2_t,double>(n,y);
        if (min_y<=0.0) {
          param_list.push_back({-15,-13,-11,-9});
        } else {
          param_list.push_back({log10(min_y)-15.0,log10(min_y)-13.0,
              log10(min_y)-11.0,log10(min_y)-9.0});
        }

        ikon->verbose=2;
        ikon->set_covar(*cfrn,param_list,true);
        
      } else {
        O2SCL_ERR((((std::string)"Invalid interpolation type, ")+
                   o2scl::szttos(interp_type)+", in "+
                   "interp_vec::set().").c_str(),exc_einval);
      }
      itype=interp_type;
    
      itp->set(n,x,y);

      return;
    }
    //@}

    /// \name Interpolation methods
    //@{
    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(const double x0) const {
      if (itp==0) {
        O2SCL_ERR("No vector set in interp_vec::eval().",
                  exc_einval);
      }
      return itp->eval(x0);
    }                   
    
    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double operator()(double x0) const {
      if (itp==0) {
        O2SCL_ERR("No vector set in interp_vec::operator().",
                  exc_einval);
      }
      return itp->eval(x0);
    }
    
    /// Give the value of the derivative \f$ y^{\prime}(x=x_0) \f$ .
    virtual double deriv(const double x0) const {
      if (itp==0) {
        O2SCL_ERR("No vector set in interp_vec::deriv().",
                  exc_einval);
      }
      return itp->deriv(x0);
    }                   
    
    /** \brief Give the value of the second derivative  
        \f$ y^{\prime \prime}(x=x_0) \f$ .
    */
    virtual double deriv2(const double x0) const {
      if (itp==0) {
        O2SCL_ERR("No vector set in interp_vec::deriv2().",
                  exc_einval);
      }
      return itp->deriv2(x0);
    }                   
    
    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(const double x1, const double x2) const {
      if (itp==0) {
        O2SCL_ERR("No vector set in interp_vec::integ().",
                  exc_einval);
      }
      return itp->integ(x1,x2);
    }
    //@}

    /// \name Other methods
    //@{
    /** \brief Clear the base interpolation object and covariance function
        (if necessary)
     */
    void clear() {
      if (itp!=0) {
        if (cf!=0) {
          delete cf;
          cf=0;
        }
        delete itp;
        itp=0;
      }
      return;
    }

    /// Return the type, "interp_vec"
    virtual const char *type() const {
      return "interp_vec";
    }
    //@}

#ifndef DOXYGEN_INTERNAL

  private:
  
    interp_vec<vec_t,vec2_t>(const interp_vec<vec_t,vec2_t> &it);
    interp_vec<vec_t,vec2_t> &operator=
    (const interp_vec<vec_t,vec2_t> &it);

#endif
  
  };

  /// \name A function for inverse interpolation in src/base/interp.h
  //@{
  /** \brief Count level crossings

      This function counts the number of times the function \f$ y(x) =
      \mathrm{level} \f$ where the function is defined by the vectors
      \c x and \c y.

      The value returned is exactly the same as the size of the
      \c locs vector computed by \ref vector_find_level().

      This function will call the error handler if \c n is
      less than two. 

      If one of the entries in \c y is either larger or smaller than
      its neighbors (i.e. if the function \f$ y(x) \f$ has an
      extremum), and if the value of \c level is exactly equal to the
      extremum, then this is counted as 1 level crossing. The same
      applies if either the first or last entry in \c y is exactly
      equal to \c level .  Finite-precision errors may affect
      the result of this function when \c level is nearly
      equal to one of the value in the vector \c y.
  */
  template<class vec_t, class vec2_t> size_t vector_level_count
  (double level, size_t n, vec_t &x, vec2_t &y) {

    if (n<=1) {
      O2SCL_ERR2("Need at least two data points in ",
                 "vector_find_count().",exc_einval);
    }
    
    size_t count=0;

    // Handle left boundary 
    if (y[0]==level) count++;

    // Find intersections
    for(size_t i=0;i<n-1;i++) {

      if ((y[i]<level && y[i+1]>=level) ||
          (y[i]>level && y[i+1]<=level)) {
        count++;
      }
    }

    return count;
  }
  //@}

  /// \name Derivatives and integrals of vectors in src/base/interp.h
  //@{
  /** \brief Compute derivative at all points from an 
      interpolation object
  */
  template<class ovec_t, class vec2_t>
  void vector_deriv_interp(size_t n, ovec_t &v, vec2_t &dv, 
                           size_t interp_type=itp_linear) {
    ovec_t grid(n);
    for(size_t i=0;i<n;i++) grid[i]=((double)i);
    interp_vec<ovec_t,ovec_t> oi(n,grid,v,interp_type);
    for(size_t i=0;i<n;i++) dv[i]=oi.deriv(((double)i));
    return;
  }

  /** \brief Compute second derivative at all points from an 
      interpolation object
  */
  template<class ovec_t, class vec2_t>
  void vector_deriv2_interp(size_t n, ovec_t &v, vec2_t &dv, 
                            size_t interp_type=itp_linear) {
    ovec_t grid(n);
    for(size_t i=0;i<n;i++) grid[i]=((double)i);
    interp_vec<ovec_t,ovec_t> oi(n,grid,v,interp_type);
    for(size_t i=0;i<n;i++) dv[i]=oi.deriv2(((double)i));
    return;
  }

  /** \brief Compute derivative at all points from an 
      interpolation object
  */
  template<class vec_t, class vec2_t, class vec3_t>
  void vector_deriv_xy_interp(size_t n, vec_t &vx, vec2_t &vy, vec3_t &dv, 
                              size_t interp_type=itp_linear) {
    interp_vec<vec_t,vec2_t> oi(n,vx,vy,interp_type);
    for(size_t i=0;i<n;i++) dv[i]=oi.deriv(vx[i]);
    return;
  }

  /** \brief Compute second derivative at all points from an 
      interpolation object
  */
  template<class vec_t, class vec2_t, class vec3_t>
  void vector_deriv2_xy_interp(size_t n, vec_t &vx, vec2_t &vy, vec3_t &dv, 
                               size_t interp_type=itp_linear) {
    interp_vec<vec_t,vec2_t> oi(n,vx,vy,interp_type);
    for(size_t i=0;i<n;i++) dv[i]=oi.deriv(vx[i]);
    return;
  }

  /** \brief Integral of a vector from interpolation object
   */
  template<class ovec_t>
  double vector_integ_interp(size_t n, ovec_t &v, size_t interp_type) {
    ovec_t grid(n);
    for(size_t i=0;i<n;i++) grid[i]=((double)i);
    interp_vec<ovec_t> oi(n,grid,v,interp_type);
    return oi.integ(0.0,((double)(n-1)));
  }

  /** \brief Compute the integral over <tt>y(x)</tt> using 
      interpolation

      Assuming vectors \c y and \c x define a function \f$ y(x) \f$
      then this computes 
      \f[
      \int_{x_0}^{x^{n-1}} y(x)~dx
      \f]
      using interpolation to approximate the result.

      See also \ref vector_deriv_interp() and 
      \ref vector_deriv2_interp() in \ref vector_derint.h .
  */
  template<class vec_t, class vec2_t> 
  double vector_integ_xy_interp(size_t n, const vec_t &x, const vec2_t &y,
                                size_t interp_type=itp_linear) {
    
    // Interpolation object
    interp_vec<vec_t,vec2_t> si(interp_type);

    // Compute full integral
    si.set(n,x,y);
    double total=si.integ(x[0],x[n-1]);//,n,x,y);

    return total;
  }

  /** \brief Compute integral over <tt>y(x)</tt> and store result
      in a vector using interpolation
  */
  template<class vec_t, class vec2_t, class vec3_t> 
  void vector_integ_xy_interp(size_t n, const vec_t &x, const vec2_t &y,
                              vec3_t &iy, size_t interp_type=itp_linear) {
    
    // Interpolation object
    interp_vec<vec_t,vec2_t> si(interp_type);
    
    // Compute full integral
    iy[0]=0.0;
    si.set(n,x,y);
    for(size_t i=1;i<n;i++) {
      iy[i]=si.integ(x[0],x[i]);//,n,x,y);
    }

    return;
  }

  /** \brief Compute the integral of a vector using
      interpolation up to a specified upper limit
  */
  template<class ovec_t>
  double vector_integ_ul_interp(size_t n, double x2,
                                ovec_t &v, size_t interp_type) {
    ovec_t grid(n);
    for(size_t i=0;i<n;i++) grid[i]=((double)i);
    interp_vec<ovec_t> oi(n,grid,v,interp_type);
    return oi.integ(0.0,x2);
  }

  /** \brief Compute the integral over <tt>y(x)</tt> using 
      interpolation up to a specified upper limit
  */
  template<class vec_t, class vec2_t> 
  double vector_integ_ul_xy_interp(size_t n, double x2,
                                   const vec_t &x, const vec2_t &y, 
                                   size_t interp_type=itp_linear) {
    
    // Interpolation object
    interp_vec<vec_t,vec2_t> si(interp_type);

    // Compute full integral
    si.set(n,x,y);
    double total=si.integ(x[0],x2);//,n,x,y);

    return total;
  }
  //@}

  /// \name Inverse interpolation and related in src/base/interp.h
  //@{
  /** \brief Perform inverse linear interpolation

      This function performs inverse linear interpolation of the data
      defined by \c x and \c y, finding all points in \c x which have
      the property \f$ \mathrm{level} = y(x) \f$. All points for which
      this relation holds are put into the vector \c locs. The
      previous information contained in vector \c locs before the
      function call is destroyed.

      This is the 1-dimensional analog of finding contour levels as
      the \ref contour class does for 2 dimensions.

      This function will call the error handler if \c n is
      less than two. 

      This function is inclusive of the endpoints, so that if either
      <tt>y[0]</tt> and/or <tt>y[n-1]</tt> is equal to level, then
      <tt>x[0]</tt> and/or <tt>x[n-1]</tt> are added to \c locs,
      respectively.
  */
  template<class vec_t, class vec2_t> void vector_find_level
  (double level, size_t n, vec_t &x, vec2_t &y, std::vector<double> &locs) {
    
    if (n<=1) {
      O2SCL_ERR2("Need at least two data points in ",
                 "vector_find_level().",exc_einval);
    }
    
    // Ensure that the location vector is empty
    locs.clear();

    // Handle left boundary 
    if (y[0]==level) {
      locs.push_back(x[0]);
    }
    
    // Find intersections
    for(size_t i=0;i<n-1;i++) {
      
      if ((y[i]<level && y[i+1]>level) ||
          (y[i]>level && y[i+1]<level)) {
        // For each intersection, add the location using linear 
        // interpolation
        double x0=x[i]+(x[i+1]-x[i])*(level-y[i])/(y[i+1]-y[i]);
        locs.push_back(x0);
      } else if (y[i+1]==level) {
        locs.push_back(x[i+1]);
      }
    }

    return;
  }

  /** \brief Compute the endpoints which enclose the regions whose
      integral is equal to \c sum

      Defining a new function, \f$ g(y_0) \f$ which takes as input any
      y-value, \f$ y_0 \f$ from the function \f$ y(x) \f$ (specified
      with the parameters \c x and \c y) and outputs the integral of
      the function \f$ y(x) \f$ over all regions where \f$ y(x) > y_0
      \f$. This function inverts \f$ g(y) \f$, taking the value
      of an integral as input, and returns the corresponding y-value
      in the variable \c lev. 

      This function is particularly useful, for example, in computing 
      the region which defines 68\% around a peak of data, thus 
      providing \f$ 1~\sigma \f$ confidence limits. 

      By default, this function does not allow any enclosed regions to
      go beyond the x region specified by the data. In some cases, it
      is useful to fix the boundaries to zero to ensure the integral
      is well-defined. If \c boundaries is set to 1, then the LHS
      boundary is set to zero, if \c boundaries is set to 2, then the
      RHS boundary is set to zero, and if \c boundaries is set to 3,
      then both boundaries are set to zero.

      Even if the boundaries are set to zero, the region enclosing a
      particular integral may not be well-defined, and this function
      can fail to find a region given a specified value of \c sum.
      Linear interpolation is used to describe the function \f$ g \f$,
      and the precision of this function is limited by this
      assumption. This function may also sometimes fail if \c sum is
      very close to the minimum or maximum value of the function \f$ g
      \f$.

      \comment
      Note that the two vector types for x and y must be the
      same in order to use o2scl_interp.
      \endcomment
  */
  template<class vec_t, class vec2_t> int vector_invert_enclosed_sum
  (double sum, size_t n, vec_t &x, vec2_t &y, double &lev,
   int boundaries=0, int verbose=0, bool err_on_fail=true) {
    
    if (n<=1) {
      O2SCL_ERR2("Need at least two data points in ",
                 "vector_invert_enclosed_sum().",exc_einval);
    }

    typedef boost::numeric::ublas::vector<double> ubvector;

    // Handle boundaries
    std::vector<double> x2, y2;
    size_t n2;
    if (boundaries==1) {
      if (verbose>0) {
        std::cout << "Fix left boundary to zero." << std::endl;
      }
      x2.resize(n+1);
      y2.resize(n+1);
      x2[0]=x[0]-(x[1]-x[0])/1.0e6;
      y2[0]=0.0;
      for(size_t i=0;i<n;i++) {
        x2[i+1]=x[i];
        y2[i+1]=y[i];
      }
      n2=n+1;
    } else if (boundaries==2) {
      if (verbose>0) {
        std::cout << "Fix right boundary to zero." << std::endl;
      }
      x2.resize(n+1);
      y2.resize(n+1);
      for(size_t i=0;i<n;i++) {
        x2[i]=x[i];
        y2[i]=y[i];
      }
      x2[n]=x[n-1]+(x[n-1]-x[n-2])/1.0e6;
      y2[n]=0.0;
      n2=n+1;
    } else if (boundaries==3) {
      if (verbose>0) {
        std::cout << "Fix both boundaries to zero." << std::endl;
      }
      x2.resize(n+2);
      y2.resize(n+2);
      x2[0]=x[0]-(x[1]-x[0])/1.0e6;
      y2[0]=0.0;
      for(size_t i=0;i<n;i++) {
        x2[i+1]=x[i];
        y2[i+1]=y[i];
      }
      x2[n+1]=x[n-1]+(x[n-1]-x[n-2])/1.0e6;
      y2[n+1]=0.0;
      n2=n+2;
    } else {
      if (verbose>0) {
        std::cout << "No boundary extrapolation." << std::endl;
      }
      x2.resize(n);
      y2.resize(n);
      for(size_t i=0;i<n;i++) {
        x2[i]=x[i];
        y2[i]=y[i];
      }
      n2=n;
    }

    // Construct a sorted list of function values 
    ubvector ysort(n2);
    vector_copy(n2,y2,ysort);
    vector_sort<ubvector,double>(ysort.size(),ysort);

    // Create list of y-values to perform y-value and integral
    // interpolation. We choose values in between the grid points to
    // ensure that we don't accidentally land precisely on a peak or
    // valley.
    std::vector<double> ylist;
    for(size_t i=0;i<ysort.size()-1;i++) {
      if (ysort[i]!=ysort[i+1]) {
        ylist.push_back((ysort[i+1]+3.0*ysort[i])/4.0);
        ylist.push_back((ysort[i+1]*3.0+ysort[i])/4.0);
      }
    }
    
    // Interpolation objects. We need two, one for the
    // user-specified vector type, and one for std::vector<double>
    interp_vec<vec_t,vec2_t> itp(itp_linear);
    interp_vec<std::vector<double>,std::vector<double> > itp2(itp_linear);
    
    // Construct vectors for interpolation
    std::vector<double> xi, yi;
    
    size_t nfail=0;

    for(size_t k=0;k<ylist.size();k++) {
      double lev_tmp=ylist[k];
      std::vector<double> locs;
      vector_find_level(lev_tmp,n2,x2,y2,locs);
      if ((locs.size()%2)!=0) {
        nfail++;
        if (verbose>0) {
          std::cout << k << " " << lev_tmp << " " << 0.0 << " "
                    << locs.size() << " (fail)" << std::endl;
        }
      } else {
        double sum_temp=0.0;
        itp2.set(n2,x2,y2);
        for(size_t i=0;i<locs.size()/2;i++) {
          double x0=locs[2*i];
          double x1=locs[2*i+1];
          sum_temp+=itp2.integ(x0,x1);//,n2,x2,y2);
        }
        xi.push_back(sum_temp);
        yi.push_back(lev_tmp);
        if (verbose>0) {
          std::cout << k << " " << lev_tmp << " " << sum_temp << " "
                    << locs.size() << " ";
          for(size_t i=0;i<locs.size();i++) {
            std::cout << locs[i] << " ";
          }
          std::cout << std::endl;
        }
      }
    }
    if (nfail>10) {
      if (err_on_fail) {
        O2SCL_ERR2("More than 10 failures in ",
                   "vector_invert_enclosed_sum().",o2scl::exc_efailed);
      }
      return o2scl::exc_efailed;
    }

    itp2.set(xi.size(),xi,yi);
    lev=itp2.eval(sum);
    //lev=itp2.eval(sum,xi.size(),xi,yi);
    
    return 0;
  }
  
  /** \brief Find the region enclosing an integral
   */
  template<class vec_t, class vec2_t> int vector_region_int
  (size_t n, vec_t &x, vec2_t &y, double intl, std::vector<double> &locs,
   int boundaries=0, int verbose=0, bool err_on_fail=true) {

    // Total integral
    double total=vector_integ_xy_interp(n,x,y,itp_linear);
    if (verbose>0) {
      std::cout << "Total integral: " << total << std::endl;
    }
    // Specified fractional integral
    if (verbose>0) {
      std::cout << "Target integral: " << intl << std::endl;
    }
    // Find correct level
    double lev;
    int ret=vector_invert_enclosed_sum(intl,n,x,y,lev,
                                       boundaries,verbose,err_on_fail);
    if (ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Failed to find a level which enclosed the ",
                   "specified integral in vector_region_int().",
                   o2scl::exc_efailed);
      }
      return o2scl::exc_efailed;
    }
    if (verbose>0) {
      std::cout << "Level from vector_invert: " << lev << std::endl;
    }

    // Inverse interpolate to find locations corresponding to level
    vector_find_level(lev,n,x,y,locs);
    if (verbose>0) {
      std::cout << "Locations from vector_find_level: ";
      for(size_t i=0;i<locs.size();i++) {
        std::cout << locs[i];
        if (i!=locs.size()-1) {
          std::cout << " ";
        }
      }
      std::cout << std::endl;
    }
    return 0;
  }

  /** \brief Find the region enclosing a partial integral
   */
  template<class vec_t, class vec2_t> int vector_region_fracint
  (size_t n, vec_t &x, vec2_t &y, double frac, std::vector<double> &locs,
   int boundaries=0, int verbose=0, bool err_on_fail=true) {
    double total=vector_integ_xy_interp(n,x,y,itp_linear);
    return vector_region_int(n,x,y,frac*total,locs,boundaries,
                             verbose,err_on_fail);
  }

  /** \brief Find the boundaries of the region enclosing a integral
      
      This function finds the boundaries of the region which
      has integral equal to <tt>frac</tt> times the full
      integral from the lower x limit to the upper x limit.
  */
  template<class vec_t, class vec2_t> int vector_bound_fracint
  (size_t n, vec_t &x, vec2_t &y, double frac, double &low, double &high,
   int boundaries=0, int verbose=0, bool err_on_fail=true) {
    
    std::vector<double> locs;
    int ret=vector_region_fracint(n,x,y,frac,locs,boundaries,
                                  verbose,err_on_fail);
    if (locs.size()==0 || ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Zero level crossings or vector_region_fracint() ",
                   "failed in vector_bound_sigma().",exc_efailed);
      }
      return o2scl::exc_efailed;
    }
    // Return min and max location
    size_t ix;
    vector_min(locs.size(),locs,ix,low);
    vector_max(locs.size(),locs,ix,high);
    return 0;
  }
  
  /** \brief Find the boundaries of the region enclosing a integral
      
      This function finds the boundaries of the region which
      has integral equal to <tt>frac</tt> times the full
      integral from the lower x limit to the upper x limit.
  */
  template<class vec_t, class vec2_t> int vector_bound_int
  (size_t n, vec_t &x, vec2_t &y, double frac, double &low, double &high,
   int boundaries=0, int verbose=0, bool err_on_fail=true) {
    
    std::vector<double> locs;
    int ret=vector_region_int(n,x,y,frac,locs,boundaries,
                              verbose,err_on_fail);
    if (locs.size()==0 || ret!=0) {
      if (err_on_fail) {
        O2SCL_ERR2("Zero level crossings or vector_region_int() ",
                   "failed in vector_bound_sigma().",exc_efailed);
      }
      return o2scl::exc_efailed;
    }
    // Return min and max location
    size_t ix;
    vector_min(locs.size(),locs,ix,low);
    vector_max(locs.size(),locs,ix,high);
    return 0;
  }
  //@}

  /// \name Binning and log vs. linear in src/base/interp.h
  //@{
  /** \brief From a pair of vectors, \f$ (x,y) \f$ create a new pair
      of vectors, \f$ (x_{\mathrm{out}},y_{\mathrm{out}}) \f$ where
      \f$ x_{\mathrm{out}} \f$ is uniformly-spaced.
      
      This function uses interpolation (method determined by
      \c interp_type) to create new vectors \c x_out and 
      \c y_out which will be resized using the <tt>resize()</tt>
      method to ensure they are of size \c n_pts. This 
      function creates an object of type \ref o2scl::interp_vec
      to do the interpolation.

      This function is used in \ref linear_or_log_chi2() .

      \note The two vectors, x and y, must have the same size, 
      as reported by their <tt>size()</tt> method.

      \note In order to use interpolation, the input vector 
      \c x must be either increasing or decreasing.
  */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t>
  void rebin_xy(const vec_t &x, const vec2_t &y,
                vec3_t &x_out, vec4_t &y_out, size_t n_pts,
                size_t interp_type) {
    
    if (x.size()!=y.size()) {
      O2SCL_ERR2("The x and y vectors must have the same size ",
                 "in rebin_xy().",o2scl::exc_einval);
    }
    if (n_pts<2) {
      O2SCL_ERR2("Number of points must be at least 2 ",
                 "in rebin_xy().",o2scl::exc_einval);
    }

    // Vector sizes
    size_t n=x.size();
    
    // Create the grid and setup x_out
    uniform_grid_end<double> ugx(x[0],x[n-1],n_pts-1);
    ugx.vector(x_out);
    
    // Allocate space and interpolate into y_out
    y_out.resize(n_pts);
    interp_vec<vec_t,vec2_t> itp(n,x,y,interp_type);
    for(size_t i=0;i<n_pts;i++) {
      y_out[i]=itp.eval(x_out[i]);
    }
    
    return;
  }

  /** \brief From a pair of vectors, \f$ (x,y) \f$ create a new pair
      of vectors, \f$ (x_{\mathrm{out}},y_{\mathrm{out}}) \f$ where
      \f$ x_{\mathrm{out}} \f$ is uniformly-spaced 

      \note Experimental

      This function rebins the vectors \c x and \c y twice,
      with two different interpolation types, and 
      then compares the accuracy of the result to the 
      user-specified value \c acc. 

      \verbatim embed:rst
      .. todo:: 
      
      In function rebin_xy(): I'm not sure what the purpose of this
      function was originally.

      \endverbatim
  */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t>
  int rebin_xy(const vec_t &x, const vec2_t &y,
               vec3_t &x_out, vec4_t &y_out, size_t n_pts,
               size_t interp_type1, size_t interp_type2,
               double acc=1.0e-4) {
    
    if (x.size()!=y.size()) {
      O2SCL_ERR2("The x and y vectors must have the same size ",
                 "in rebin_xy().",o2scl::exc_einval);
    }
    if (n_pts<2) {
      O2SCL_ERR2("Number of points must be at least 2 ",
                 "in rebin_xy().",o2scl::exc_einval);
    }

    // Vector sizes
    size_t n=x.size();
    
    // Create the grid and setup x_out
    uniform_grid_end<double> ugx(x[0],x[n-1],n_pts-1);
    ugx.vector(x_out);
    
    // Allocate space and interpolate into y_out
    std::vector<double> y_out2(n_pts);
    y_out.resize(n_pts);
    interp_vec<vec3_t,vec4_t> itp1(n,x,y,interp_type1);
    interp_vec<vec3_t,vec4_t> itp2(n,x,y,interp_type2);
    for(size_t i=0;i<n_pts;i++) {
      y_out[i]=itp1.eval(x_out[i]);
      y_out2[i]=itp2.eval(x_out[i]);
    }

    double max_y, min_y;
    vector_minmax_value(n_pts,y_out,min_y,max_y);
    for(size_t i=0;i<n_pts;i++) {
      if (fabs(y_out2[i]-y_out[i])/(max_y-min_y)>acc) {
        return 1;
      }
    }
    
    return 0;
  }

  /** \brief Rebin, rescale, sort, and match to \f$ y=x \f$
      
      This function rebins, rescales, and sorts two vectors \c x and
      \c y in order to compare them to the line \f$ y=x \f$, returning
      the value of \f$ \chi^2 \f$.

      Principally used by \ref linear_or_log() .

      \future Rename this function. It's a bit confusing because
      there are no logarithms here. Maybe just rescale_fit_identity()?
  */
  template<class vec_t, class vec2_t>
  double linear_or_log_chi2(const vec_t &x, const vec2_t &y) {

    // Just choose a fiducial number of points so that we
    // can easily rescale to [0,1]
    static const size_t n2=11;
    
    // Rebin into vectors of length 11
    std::vector<double> xn, yn;
    rebin_xy(x,y,xn,yn,n2,itp_linear);
    
    // Rescale to [0,1]
    double min_y, max_y;
    vector_minmax_value(n2,yn,min_y,max_y);
    for(size_t i=0;i<n2;i++) {
      yn[i]=(yn[i]-min_y)/(max_y-min_y);
    }
    
    // Sort and match to the line y=x
    vector_sort_double(n2,yn);
    double chi2=0.0;
    for(size_t i=0;i<n2;i++) {
      double ty=((double)i)/10.0;
      chi2+=pow(yn[i]-ty,2.0);
    }

    return chi2;
  }

  /** \brief Attempt to determine if data represented by two arrays,
      (x,y), would be better plotted on a linear, semi-log, or log-log
      plot

      \note Experimental.

      This function attempts to guess whether the data stored in \c x
      and \c y might be best plotted on a log scale. This algorithm
      will fail for poorly sampled or highly oscillatory data.

      \note The two vectors, x and y, must have the same size, 
      as reported by their <tt>size()</tt> method.

      This uses the \ref linear_or_log_chi2() function to determine
      whether a log of x ore y improves the fit to the line \f$ y=x
      \f$.
  */
  template<class vec_t, class vec2_t>
  void linear_or_log(vec_t &x, vec2_t &y, bool &log_x, bool &log_y) {
    if (x.size()!=y.size()) {
      O2SCL_ERR2("The x and y vectors must have the same size ",
                 "in linear_or_log().",o2scl::exc_einval);
    }
    
    // Vector sizes
    size_t n=x.size();
    if (n<2) {
      O2SCL_ERR2("Vector size must be at least 2 ",
                 "in linear_or_log().",o2scl::exc_einval);
    }
    
    // Initial values of log_x and log_y
    log_x=false;
    log_y=false;

    // Check if vectors are positive
    bool x_positive=true;
    bool y_positive=true;
    for(size_t i=0;i<n;i++) {
      if (x[i]<=0.0) x_positive=false;
      if (y[i]<=0.0) y_positive=false;
    }

    // If not, presume linear
    if (x_positive==false && y_positive==false) return;

    double chi2=linear_or_log_chi2(x,y);
    
    // Take the log
    std::vector<double> lx(n), ly(n);
    if (x_positive) {
      for(size_t i=0;i<n;i++) {
        lx[i]=log(x[i]);
      }
    }
    if (y_positive) {
      for(size_t i=0;i<n;i++) {
        ly[i]=log(y[i]);
      }
    }

    // Depending on whether or not they are positive, find the best match
    if (x_positive) {
      if (y_positive) {
        double chi2_x=linear_or_log_chi2(lx,y);
        double chi2_y=linear_or_log_chi2(x,ly);
        double chi2_xy=linear_or_log_chi2(lx,ly);
        if (chi2_xy<chi2_x && chi2_xy<chi2_y && chi2_xy<chi2) {
          log_x=true;
          log_y=true;
        } else if (chi2_x<chi2_xy && chi2_x<chi2_y && chi2_x<chi2) {
          log_x=true;
        } else if (chi2_y<chi2_xy && chi2_y<chi2_x && chi2_y<chi2) {
          log_y=true;
        }
      } else {
        double chi2_x=linear_or_log_chi2(lx,y);
        if (chi2_x<chi2) log_x=true;
      }
    } else {
      double chi2_y=linear_or_log_chi2(x,ly);
      if (chi2_y<chi2) log_y=true;
    }

    return;
  }

  /** \brief Fill in between vector entries with linear interpolation
   */
  template<class vec_t>
  void vector_refine_inplace(vec_t &v, size_t factor, bool log=false) {
    if (factor<2) {
      O2SCL_ERR("Factor less than 2 in vector_refine_inplace().",
                o2scl::exc_einval);
    }
    vec_t v2((v.size()-1)*factor+1);
    v2[v2.size()-1]=v[v.size()-1];
    if (log) {
      size_t ix=0;
      for(size_t i=0;i<v.size()-1;i++) {
        for(size_t j=0;j<factor;j++) {
          v2[ix]=v[i]*pow(v[i+1]/v[i],((double)j)/((double)factor));
          ix++;
        }
      }
    } else {
      size_t ix=0;
      for(size_t i=0;i<v.size()-1;i++) {
        for(size_t j=0;j<factor;j++) {
          v2[ix]=((double)j)/((double)factor)*(v[i+1]-v[i])+v[i];
          ix++;
        }
      }
    }
    v=v2;
    return;
  }
  
  /** \brief Refine a vector by interpolating with a second
      index vector

      \warning Untested.
  */
  template<class vec_t, class vec2_t, class data_t>
  void vector_refine(size_t n, const vec_t &index, vec2_t &data,
                     size_t factor, size_t interp_type=itp_linear) {
    interp_vec<vec_t,vec2_t> iv(n,index,data,interp_type);
    vec2_t copy=data;
    data.resize((n-1)*factor+1);
    for (size_t j=0;j<n-1;j++) {
      for(size_t k=0;k<factor;k++) {
        data[j*factor+k]=iv.eval(index[j]+
                                 ((data_t)k)/((data_t)factor)*
                                 (index[j+1]-index[j]));
      }
    }
    data[data.size()-1]=copy[copy.size()-1];
    return;
  }
  
  /** \brief Attempt to determine if data stored in \c y
      is more nearly linear or logarithmic
      
      \note Experimental.

      \future There is a bit of extra calculation because the rebin
      function creates a new x vector with a uniform grid, so this
      could be streamlined.
  */
  template<class vec_t>
  void linear_or_log(vec_t &y, bool &log_y) {
    std::vector<double> x(y.size());
    for(size_t i=0;i<y.size();i++) {
      x[i]=((double)i);
    }
    bool log_x;
    linear_or_log(x,y,log_x,log_y);
    return;
  }
  //@}
  
}

#endif
