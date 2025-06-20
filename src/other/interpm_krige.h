/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2017-2025, Andrew W. Steiner and Josue Bautista
  
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
#ifndef O2SCL_INTERPM_KRIGE_H
#define O2SCL_INTERPM_KRIGE_H

/** \file interpm_krige.h
    \brief File defining \ref o2scl::interpm_krige_optim
*/

#include <iostream>
#include <string>
#include <cmath>
#include <ctime>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <gsl/gsl_combination.h>

#include <o2scl/interpm_base.h>
#include <o2scl/vector.h>
#include <o2scl/vec_stats.h>
#include <o2scl/linear_solver.h>
#include <o2scl/columnify.h>
#include <o2scl/cholesky.h>
#include <o2scl/min_brent_gsl.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/mmin_conf.h>
#include <o2scl/diff_evo_adapt.h>
#include <o2scl/mmin_bfgs2.h>
#include <o2scl/cblas.h>
#include <o2scl/invert.h>

namespace o2scl {

  typedef boost::numeric::ublas::vector<double> ubvector;

  /** \brief Multidimensional covariance function base

      It is expected that the covariance function is symmetric
      with respect to its arguments. The derivative is then
      antisymmetric. 
   */
  template<class vec_t, class vec2_t> class mcovar_base {
    
  public:
    
    /// Get the number of parameters
    virtual size_t get_n_params()=0;
    
    /// Set the parameters
    virtual void set_params(const ubvector &p)=0;

    /// The covariance function
    virtual double operator()(const vec_t &x1, const vec2_t &x2)=0;
    
    /// The covariance function
    virtual double covar(const vec_t &x1, const vec_t &x2)=0;
    
    /// The covariance function
    virtual double covar2(const vec2_t &x1, const vec2_t &x2)=0;
    
    /** \brief The derivative of the covariance function with
        respect to the first argument
    */
    virtual double deriv(const vec_t &x1, const vec2_t &x2,
                         size_t ix)=0;
    
    /** \brief The second derivative of the covariance function with
        respect to the first argument
    */
    virtual double deriv2(const vec_t &x1, const vec2_t &x2,
                          size_t ix, size_t iy)=0;
    
  };
  
  /** \brief Multidimensional covariance function using radial
      basis functions and a noise term
  */
  template<class vec_t, class vec2_t>
  class mcovar_funct_rbf_noise : public mcovar_base<vec_t,vec2_t> {
    
  public:

    /// Length parameters
    std::vector<double> len;
    
    /// Noise parameter
    double log10_noise;
    
    /// Get the number of parameters
    size_t get_n_params() {
      return len.size()+1;
    }
    
    /// Set the parameters
    void set_params(const vec_t &p) {
      for(size_t j=0;j<len.size();j++) {
        len[j]=p[j];
      }
      log10_noise=p[len.size()];
      return;
    }
    
    /// The covariance function template
    template<class vec3_t, class vec4_t>
    double covar_tl(const vec3_t &x1, const vec4_t &x2){
      double sum=0.0;
      bool equal=true;
      for(size_t j=0;j<len.size();j++) {
        if (x1[j]!=x2[j]) equal=false;
        sum+=-(x1[j]-x2[j])*(x1[j]-x2[j])/len[j]/len[j]/2.0;
      }
      if (equal) return exp(sum)+pow(10.0,log10_noise);
      return exp(sum);
    }

    /// The covariance function
    virtual double operator()(const vec_t &x1, const vec2_t &x2){
      return covar_tl(x1,x2);
    }
    
    /// The covariance function
    virtual double covar(const vec_t &x1, const vec_t &x2){
      return covar_tl(x1,x2);
    }
    
    /// The covariance function
    virtual double covar2(const vec2_t &x1, const vec2_t &x2){
      return covar_tl(x1,x2);
    }
    
    /** \brief The derivative of the covariance function with
        respect to the first argument
    */
    virtual double deriv(const vec_t &x1, const vec2_t &x2, size_t ix){
      double sum=0.0;
      for(size_t j=0;j<len.size();j++) {
        sum+=-(x1[j]-x2[j])*(x1[j]-x2[j])/len[j]/len[j]/2.0;
      }
      return -exp(sum)/len[ix]/len[ix]*(x1[ix]-x2[ix]);
    }

    /** \brief The second derivative of the covariance function with
        respect to the first argument
    */
    virtual double deriv2(const vec_t &x1, const vec2_t &x2, size_t ix,
                          size_t iy){
      double sum=0.0;
      for(size_t j=0;j<len.size();j++) {
        sum+=-(x1[j]-x2[j])*(x1[j]-x2[j])/len[j]/len[j]/2.0;
      }
      if (ix==iy) {
        return exp(sum)/len[ix]/len[ix]/len[ix]/len[ix]*
          ((x1[ix]-x2[ix])*(x1[ix]-x2[ix])-len[ix]*len[ix]);
      }
      return exp(sum)/len[ix]/len[ix]*(x1[ix]-x2[ix])/
        len[iy]/len[iy]*(x1[iy]-x2[iy]);
    }
    
  };

  /** \brief Multidimensional covariance function using radial basis
      functions with quadratic correlation length and a noise term
      
      The covariance function is
      \f[
      K_{ij} = \prod_k \exp \left[
      \frac{ - \left(x_{ik} - x_{jk}\right)^2}
      {2 d_k^2} \right] 
      \f]
      where
      \f[
      d_k^2 \equiv \left[\ell_k^2 + m_k^2 \left( x_{ik} + x_{jk} -
      n_k\right)^2\right] \, .
      \f]
      The values of \f$ \ell_k \f$ are stored in \ref len, the
      values of \f$ m_k \f$ are stored in \ref slope, and the values
      of \f$ n_k \f$ are stored in \ref pos.
  */
  template<class vec_t, class vec2_t>
  class mcovar_funct_quad_correl : public mcovar_base<vec_t,vec2_t> {
    
  public:
    
    /// Length parameters
    std::vector<double> len;

    /// Slope parameters
    std::vector<double> slope;
    
    /// Position parameters
    std::vector<double> pos;
    
    /// Noise parameter
    double log10_noise;
    
    /// Get the number of parameters
    size_t get_n_params() {
      return len.size()*3+1;
    }
    
    /// Set the parameters
    void set_params(const vec_t &p) {
      for(size_t j=0;j<len.size();j++) {
        len[j]=p[j];
      }
      for(size_t j=len.size();j<len.size()*2;j++) {
        slope[j-len.size()]=p[j];
      }
      for(size_t j=len.size()*2;j<len.size()*3;j++) {
        pos[j-len.size()*2]=p[j];
      }
      log10_noise=p[len.size()*3];
      return;
    }

    /// The covariance function template
    template<class vec3_t, class vec4_t>
    double covar_tl(const vec3_t &x1, const vec4_t &x2){
      double sum=0.0;
      bool equal=true;
      for(size_t j=0;j<len.size();j++) {
        double ell=len[j];
        double m=slope[j];
        double n=pos[j];
        if (x1[j]!=x2[j]) equal=false;
        double den2=ell*ell+m*m*pow((x1[j]+x2[j])-n,2);
        sum+=-(x1[j]-x2[j])*(x1[j]-x2[j])/den2/2.0;
      }
      if (equal) return exp(sum)+pow(10.0,log10_noise);
      return exp(sum);
    }

    /** \brief The covariance function
     */
    virtual double operator()(const vec_t &x1, const vec2_t &x2){
      return covar_tl(x1,x2);
    }
    
    /** \brief The covariance function
     */
    virtual double covar(const vec_t &x1, const vec_t &x2){
      return covar_tl(x1,x2);
    }
    
    /** \brief The covariance function
     */
    virtual double covar2(const vec2_t &x1, const vec2_t &x2){
      return covar_tl(x1,x2);
    }

    /** \brief The template for the derivative of the
        covariance function

        This function computes
        \f[
        \frac{\partial K_{ij}}{\partial x_{ik}} =
        \frac{K_{ij} z_k}{d_k^4} \, .
        \f]
        where
        \f[
        z_k \equiv \left[m_k^2 \left( n_k - 2 x_{jk}
        \right) \left( x_{ik} + x_{jk} - n_k\right) - \ell_k^2\right]
        \left( x_{ik} - x_{jk} \right)
        \f]
    */
    template<class vec3_t, class vec4_t>
    double deriv_tl(const vec3_t &x1, const vec4_t &x2, size_t ix){
      double ell=len[ix];
      double m=slope[ix];
      double n=pos[ix];
      double K=covar_tl(x1,x2), deriv;
      double d2ix=ell*ell+m*m*pow((x1[ix]+x2[ix])-n,2);
      double zix=(m*m*(n-2.0*x2[ix])*(x1[ix]+x2[ix]-n)-ell*ell)*
        (x1[ix]-x2[ix]);
      deriv=K*zix/d2ix/d2ix;
        
      return deriv;
    }
    
    /** \brief The derivative of the covariance function with
        respect to one element of the first argument

        See \ref deriv_tl().
    */
    virtual double deriv(const vec_t &x1, const vec2_t &x2, size_t ix){
      return deriv_tl(x1,x2,ix);
    }
    
    /** \brief Template for the second derivative of the covariance
        function with respect to the first argument

        This function computes
        \f[
        \frac{\partial^2 K_{ij}}{\partial x_{ik} \partial x_{im}} \, .
        \f]
        If \f$ k \neq m \f$, then
        \f[
        \frac{\partial^2 K_{ij}}{\partial x_{ik} \partial x_{im}}
        = \frac{K_{ij} z_k z_m}{d_k^4 d_m^4} \, .
        \f]
        Otherwise
        \f[
        \frac{\partial^2 K_{ij}}{\partial x_{ik}^2} =
        \frac{K_{ij}}{d_k^4}
        \left[
        \frac{z_k^2}{d_k^4} + \frac{\partial z_k}{\partial x_{ik}}
        - \frac{4 z_k m_k^2}{d_k^2} \left( x_{ik} + x_{jk} - n_k \right)
        \right]
        \f]
        where
        \f[
        \frac{\partial z_k}{\partial x_{ik}} =
        m_k^2 \left( n_k - 2 x_{jk} \right) \left( x_{ik} - x_{jk} \right)
        - \left[ \ell_k^2 - m_k^2 \left( n_k - 2 x_{jk} \right)
        \left( x_{ik} + x_{jk} - n_k\right) \right]
        \f]
    */
    template<class vec3_t, class vec4_t>
    double deriv2_tl(const vec3_t &x1, const vec4_t &x2, size_t ix,
                     size_t iy) {
      double ellx=len[ix];
      double mx=slope[ix];
      double nx=pos[ix];
      double elly=len[iy];
      double my=slope[iy];
      double ny=pos[iy];
      double sum=0.0;
      double K=covar_tl(x1,x2), deriv2;
      double d2ix=ellx*ellx+mx*mx*pow((x1[ix]+x2[ix])-nx,2);
      double zix=(mx*mx*(nx-2.0*x2[ix])*(x1[ix]+x2[ix]-nx)-ellx*ellx)*
        (x1[ix]-x2[ix]);
      if (ix!=iy) {
        double d2iy=elly*elly+my*my*pow((x1[iy]+x2[iy])-ny,2);
        double ziy=(my*my*(ny-2.0*x2[iy])*(x1[iy]+x2[iy]-ny)-elly*elly)*
          (x1[iy]-x2[iy]);
        deriv2=K*zix*ziy/d2ix/d2ix/d2iy/d2iy;
      } else {
        double dzkdxik=mx*mx*(nx-2.0*x2[ix])*(x1[ix]-x2[ix])-
          (ellx*ellx+mx*mx*(nx-2.0*x2[ix])*(nx-x1[ix]-x2[ix]));
        deriv2=K/d2ix/d2ix*(zix*zix/d2ix/d2ix+dzkdxik-
                            2.0*zix/d2ix*mx*mx*2.0*(x1[ix]+x2[ix]-nx));
      }
      return deriv2;
    }
    
    /** \brief The second derivative of the covariance
        function with respect to the first argument
        
        See \ref deriv2_tl().
    */
    virtual double deriv2(const vec_t &x1, const vec2_t &x2, size_t ix,
                          size_t iy) {
      return deriv2_tl(x1,x2,ix,iy);
    }
    
  };

  /** \brief Multi-dimensional interpolation using Kriging (a.k.a.
      Gaussian process interpolation)

      \note The type mat_y_t should be const_matrix_view_table, but the
      code would need to be modified because when the data is rescaled
      the user copy of the data is modified. This should be fixed in
      the future, as there's no need to change the user-specified data,
      it would be better to just modify local storage.
      
      \verbatim embed:rst
      See also the :ref:`Higher-dimensional Interpolation` 
      section of the User's guide. 
      \endverbatim
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_x_t=o2scl::const_matrix_view_table<>,
           class mat_x_row_t=const const_matrix_row_gen
           <o2scl::const_matrix_view_table<>>, 
           class mat_y_t=o2scl::matrix_view_table<>,
           class mat_y_col_t=const matrix_column_gen<
             o2scl::matrix_view_table<>>,
           class mat_inv_kxx_t=boost::numeric::ublas::matrix<double>,
           class mat_inv_t=
           o2scl_linalg::matrix_invert_det_cholesky<
             boost::numeric::ublas::matrix<double>>,
           class vec3_t=std::vector<std::vector<std::vector<double>>> >
  class interpm_krige_optim :
    public interpm_base<vec_t,mat_x_t,mat_y_t> {

  public:

    /// Internal matrix type for the Nelder-Mead simplex
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    
    /// Typedef for this type
    typedef interpm_krige_optim<vec_t,mat_x_t,mat_x_row_t, 
                                mat_y_t,mat_y_col_t,mat_inv_kxx_t,
                                mat_inv_t,vec3_t> class_t;
    
  protected:

    /** \brief Inverse covariance matrix times function vector
     */
    std::vector<ubvector> Kinvf;
    
    /** \brief The inverse of the covariance matrix for each output
        quantity
    */
    std::vector<mat_inv_kxx_t> inv_KXX;

    /// The matrix inversion object
    mat_inv_t mi;

    /** \brief The list of parameter values

        The first index runs over output quantities, the second over
        different sets of values, and the third over parameters.
     */
    vec3_t plists;
  
     /// The quality factor of the optimization for each output function
    std::vector<double> qual;

    /// Pointer to the user-specified minimizer
    mmin_base<multi_funct,multi_funct,ubvector> *mp;

    /// The input data
    mat_x_t x;
    /// The output data
    mat_y_t y;
    /// True if the data has been specified
    bool data_set;
    /// The output means for rescaling
    ubvector mean_y;
    /// The output standard deviations for rescaling
    ubvector std_y;
  
  public:

    /// \name Settings
    //@{
    /** \brief If true, output timing results
     */
    bool timing;
    
    /// If true, use the full minimizer
    bool full_min;

    /// If true, skip optimization
    bool skip_optim;
    
    /** \brief If true, then the data will be automatically rescaled
        (default true)
    */
    bool rescale;
    
    /// If true, keep \f$ K^{-1} \f$ (default true)
    bool keep_matrix;
    //@}

    /// \name Minimizers and settings
    //@{
    /// Default minimizer
    mmin_simp2<multi_funct,ubvector> def_mmin;

    /// Alternate minimizer
    diff_evo_adapt<> alt_mmin;

    /// If true, use the alternate minimizer
    bool use_alt_mmin;
    
    /// Set the minimizer to use
    void set_mmin(mmin_base<multi_funct,multi_funct,ubvector> &mb) {
      mp=&mb;
      return;
    }
    //@}

    /// \name Other functions
    //@{
    /** \brief Additional constraints to add to the fit
     */
    virtual int addl_const(size_t iout, double &ret) {
      return 0;
    }
    
    /** \brief Function to optimize the covariance parameters
     */
    virtual double qual_fun(size_t iout, int &success) {

      // Select the row of the data matrix
      mat_y_col_t yiout2(this->y,iout);
      
      double ret=0.0;
    
      success=0;

      size_t size=this->x.size1();

      time_t t1=0, t2=0, t3=0, t4=0;

      if (mode==mode_loo_cv_bf) {

        for(size_t k=0;k<size;k++) {
          // Leave one observation out

          mat_x_row_t xk(this->x,k);
          
          ubvector y2(size-1);
          o2scl::vector_copy_jackknife(size,yiout2,k,y2);

          // Construct the inverse of the KXX matrix. Note that we
          // don't use the inv_KXX data member because of the size
          // mismatch
          mat_inv_kxx_t inv_KXX2(size-1,size-1);
          for(size_t irow=0;irow<size-1;irow++) {
            size_t irow2=irow;
            if (irow>=k) irow2++;        
            mat_x_row_t xrow(this->x,irow2);
            for(size_t icol=0;icol<size-1;icol++) {
              size_t icol2=icol;
              if (icol>=k) icol2++;        
              mat_x_row_t xcol(this->x,icol2);
              if (irow2>icol2) {
                inv_KXX2(irow,icol)=inv_KXX2(icol,irow);
              } else {
                inv_KXX2(irow,icol)=cf[iout]->covar2(xrow,xcol);
              }
            }
          }
          
          // Construct the inverse of KXX
          this->mi.invert_inplace(size-1,inv_KXX2);
          
          // Inverse covariance matrix times function vector
          ubvector Kinvf2(size-1);
          o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                             o2scl_cblas::o2cblas_NoTrans,
                             size-1,size-1,1.0,inv_KXX2,y2,0.0,
                             Kinvf2);
          

          // The actual value
          double yact=yiout2[k];

          // Compute the predicted value
          double ypred=0.0;
          ubvector kxx0(size-1);
          for(size_t i=0;i<size-1;i++) {
            size_t i2=i;
            if (i>=k) i2++;        
            mat_x_row_t xi2(this->x,i2);
            kxx0[i]=cf[iout]->covar2(xk,xi2);
            ypred+=kxx0[i]*Kinvf2[i];
          }

          // AWS 12/31/22: This uses absolute, rather than relative
          // differences to evaluate the quality, but this seems
          // sensible to me because we are presuming the data has zero
          // mean and unit standard deviation.
          ret+=pow(yact-ypred,2.0);
          
        }

      } else if (mode==mode_loo_cv) {

        // Construct the KXX matrix
        this->inv_KXX[iout].resize(size,size);
        for(size_t irow=0;irow<size;irow++) {
          mat_x_row_t xrow(this->x,irow);
          for(size_t icol=0;icol<size;icol++) {
            mat_x_row_t xcol(this->x,icol);
            if (irow<icol) {
              this->inv_KXX[iout](irow,icol)=this->inv_KXX[iout](icol,irow);
            } else {
              this->inv_KXX[iout](irow,icol)=cf[iout]->covar2(xrow,xcol);
            }
          }
        }

        if (this->verbose>2) {
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Done creating covariance matrix with size "
                    << size << std::endl;
        }

        if (timing) {
          t1=time(0);          
        }
        
        // Construct the inverse of KXX
        if (this->verbose>2) {
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Performing matrix inversion with size "
                    << size << std::endl;
        }
        int cret=this->mi.invert_inplace(size,this->inv_KXX[iout]);
        if (cret!=0) {
          success=1;
          return 1.0e99;
        }
	
        if (timing) {
          t2=time(0);          
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Matrix inversion took "
                    << t2-t1 << " seconds." << std::endl;
        }
        
        if (this->verbose>2) {
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Done performing matrix inversion with size "
                    << size << std::endl;
        }
        
        // Inverse covariance matrix times function vector
        this->Kinvf[iout].resize(size);
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           size,size,1.0,this->inv_KXX[iout],
                           yiout2,0.0,this->Kinvf[iout]);
        
        if (timing) {
          t3=time(0);          
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Matrix vector multiply took "
                    << t3-t2 << " seconds." << std::endl;
        }
        
        ret=0.0;
        
        // Select the row of the data matrix
        mat_y_col_t yiout(this->y,iout);
        
        for(size_t ii=0;ii<size;ii++) {
          
          double yact=yiout[ii];
          
          // Compute sigma and ypred from Eq. 5.12
          double sigma2=1.0/this->inv_KXX[iout](ii,ii);
          double ypred=yact-this->Kinvf[iout][ii]*sigma2;
          
          // Then use Eq. 5.10
          ret+=pow(yact-ypred,2.0)/sigma2/2.0;
          ret+=0.5*log(sigma2);

        }

        if (timing) {
          t4=time(0);          
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Final evaluation took "
                    << t4-t3 << " seconds." << std::endl;
        }
        
        if (this->verbose>2) {
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "ret: " << ret << std::endl;
        }
      
      } else if (mode==mode_max_lml || mode==mode_final) {

        if (this->verbose>2) {
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Creating covariance matrix with size "
                    << size << std::endl;
        }

        // Construct the KXX matrix
        mat_inv_kxx_t KXX(size,size);
        for(size_t irow=0;irow<size;irow++) {
          mat_x_row_t xrow(this->x,irow);
          for(size_t icol=0;icol<size;icol++) {
            mat_x_row_t xrow2(this->x,icol);
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else {
              KXX(irow,icol)=cf[iout]->covar2(xrow,xrow2);
            }
          }
        }

        if (this->verbose>2) {
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Done creating covariance matrix with size "
                    << size << std::endl;
        }

        // Perform the matrix inversion and compute the determinant
        
        double lndet;
	
        if (timing) {
          t1=time(0);          
        }
        
        // Construct the inverse of KXX
        if (this->verbose>2) {
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Performing matrix inversion with size "
                    << size << std::endl;
        }
        this->inv_KXX[iout].resize(size,size);
        int cret=this->mi.invert_det(size,KXX,this->inv_KXX[iout],lndet);
        if (cret!=0) {
          success=2;
          return 1.0e99;
        }
	
        lndet=log(lndet);

        // The determinant is only required if we're marginalizing the
        // maximum likelihood.
        if (mode!=mode_final && !std::isfinite(lndet)) {
          success=5;
          return 1.0e99;
        }
	
        if (this->verbose>2) {
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Done performing matrix inversion with size "
                    << size << std::endl;
        }

        if (timing) {
          t2=time(0);          
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Matrix inversion took "
                    << t2-t1 << " seconds." << std::endl;
        }
        
        // Inverse covariance matrix times function vector
        this->Kinvf[iout].resize(size);
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           size,size,1.0,this->inv_KXX[iout],
                           yiout2,0.0,this->Kinvf[iout]);
	
        if (timing) {
          t3=time(0);          
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Matrix vector multiply took "
                    << t3-t2 << " seconds." << std::endl;
        }
        
        if (mode==mode_max_lml) {
          // Compute the log of the marginal likelihood, without
          // the constant term
          for(size_t i=0;i<size;i++) {
            ret+=0.5*yiout2[i]*this->Kinvf[iout][i];
          }
          ret+=0.5*lndet;
        }

        if (timing) {
          t4=time(0);          
          std::cout << "interpm_krige_optim::qual_fun(): ";
          std::cout << "Final evaluation took "
                    << t4-t3 << " seconds." << std::endl;
        }
        
      }

      double qual_ac=0.0;
      int ret_ac=addl_const(iout,qual_ac);
      if (ret_ac!=0) {
        success=4;
      }
      ret+=qual_ac;
      
      if (!isfinite(ret)) success=3;
      
      return ret;
    }

    /** \brief Minimization function for the covariance parameters
     */
    double min_fun(size_t iout, size_t n, const ubvector &v,
                   double max_val) {
      cf[iout]->set_params(v);
      int success;
      double ret=qual_fun(iout,success);
      if (success!=0) {
        ret=max_val;
      }
      return ret;
    }
    //@}

    /// \name Constructor and destructor
    //@{
    interpm_krige_optim() {
      full_min=false;
      def_mmin.ntrial*=10;
      mp=&def_mmin;
      mode=mode_loo_cv;
      loo_npts=100;
      timing=false;
      keep_matrix=true;
      use_alt_mmin=false;
      skip_optim=false;
      rescale=true;
    }

    virtual ~interpm_krige_optim() {
    }
    //@}
    
    /// \name Function to minimize and related options
    //@{
    /// Leave-one-out cross validation (brute force version)
    static const size_t mode_loo_cv_bf=1;
    /// Minus Log-marginal-likelihood
    static const size_t mode_max_lml=2;
    /// Leave-one-out cross validation (default)
    static const size_t mode_loo_cv=3;
    /// No optimization (for internal use)
    static const size_t mode_final=10;
    /// Function to minimize (default \ref mode_loo_cv)
    size_t mode;
    /// Number of points to test for cross validation (default 100)
    size_t loo_npts;
    ///@}

    /// \name Covariance function
    //@{
    /// Pointer to the covariance function
    std::vector<std::shared_ptr<mcovar_base<vec_t,mat_x_row_t>>> cf;

    /** \brief Set the covariance function and parameter lists
     */
    int set_covar(std::vector<std::shared_ptr<mcovar_base<vec_t,
                  mat_x_row_t>>> covar,
                  vec3_t &param_lists) {
      cf=covar;
      plists=param_lists;
      return 0;
    }
    //@}

    /// \name Basic interpolation functions
    //@{
    /** \brief Initialize the data for the interpolation
     */
    virtual int set_data(size_t n_in, size_t n_out, size_t n_pts,
                         mat_x_t &user_x, mat_y_t &user_y) {
                       
      if (n_pts<2) {
        O2SCL_ERR2("Must provide at least two points in ",
                   "interpm_krige_optim::set_data()",
                   exc_efailed);
      }
      if (n_in<1) {
        O2SCL_ERR2("Must provide at least one input column in ",
                   "interpm_krige_optim::set_data()",
                   exc_efailed);
      }
      if (n_out<1) {
        O2SCL_ERR2("Must provide at least one output column in ",
                   "interpm_krige_optim::set_data()",
                   exc_efailed);
      }
   
      // Check that the data is properly sized
      if (user_x.size1()!=n_pts || user_x.size2()!=n_in) {
        std::cout << "Object user_x, function size1() and size2(): "
                  << user_x.size1() << " " << user_x.size2() << std::endl;
        std::cout << "n_points: " << this->n_points << " n_in: "
                  << n_in << std::endl;
        O2SCL_ERR2("Size of x not correct in ",
                   "interpm_krige_new::set_data().",
                   o2scl::exc_efailed);
      }
    
      if (user_y.size1()!=n_pts || user_y.size2()!=n_out) {
        std::cout << "Object user_y, function size1() and size2(): "
                  << user_y.size1() << " " << user_y.size2() << std::endl;
        O2SCL_ERR2("Size of y not correct in ",
                   "interpm_krige_new::set_data().",
                   o2scl::exc_efailed);
      }

      // Set parent data members
      this->n_points=n_pts;
      this->n_params=n_in;
      this->n_outputs=n_out;
      this->rescale=rescale;
      this->data_set=true;
      
      time_t t1=0, t2=0, t3=0, t4=0, t5=0;
      if (timing) {
        t1=time(0);
      }
      
      std::swap(this->x,user_x);
      std::swap(this->y,user_y);

      if (timing) {
        t2=time(0);
        std::cout << "interpm_krige_optim::set_data(): ";
        std::cout << "Swap took " << t2-t1 << " seconds." << std::endl;
      }
      
      if (this->verbose>0) {
        std::cout << "interpm_krige_optim::set_data(): "
                  << "Using " << this->n_points
                  << " points with\n " << n_in << " input variables and "
                  << this->n_outputs << " output variables." << std::endl;
      }

      if (rescale==true) {
        if (this->verbose>1) {
          std::cout << "interpm_krige_optim::set_data(): "
                    << "Rescaling data:" << std::endl;
        }
        this->mean_y.resize(n_out);
        this->std_y.resize(n_out);
        for(size_t j=0;j<n_out;j++) {
          mat_y_col_t vec(this->y,j);
          this->mean_y[j]=vector_mean<mat_y_col_t,double>
            (this->n_points,vec);
          this->std_y[j]=vector_stddev(this->n_points,vec);
          if (this->verbose>1) {
            std::cout << "interpm_krige_optim::set_data(): ";
            std::cout << "Mean, std. dev. of output " << j << " of "
                      << n_out << " is\n  " << this->mean_y[j] << " and "
                      << this->std_y[j] << "." << std::endl;
          }
          for(size_t i=0;i<this->n_points;i++) {
            this->y(i,j)=(this->y(i,j)-this->mean_y[j])/this->std_y[j];
          }
        }
      }

      if (timing) {
        t3=time(0);
        std::cout << "interpm_krige_optim::set_data(): ";
        std::cout << "Rescale took " << t3-t2 << " seconds." << std::endl;
      }
      
      int success=0;

      this->Kinvf.resize(n_out);
      this->inv_KXX.resize(n_out);

      if (skip_optim) return 0;
      
      qual.resize(n_out);

      // Loop over all output functions
      for(size_t iout=0;iout<n_out;iout++) {
        
        if (this->verbose>0) {
          std::cout << "interpm_krige_optim::set_data(): "
                    << "Output " << iout+1 << " of " << n_out
                    << std::endl;
        }
        
        if (timing) {
          t4=time(0);
        }
        
        // Initialize to zero to prevent uninit'ed var. warnings
        double min_qual=0.0;

        double max_val=0.0;
        size_t np_covar=cf[iout]->get_n_params();
        vec_t params(np_covar), min_params(np_covar);
        
        if (full_min) {

          if (this->verbose>0) {
            std::cout << "interpm_krige_optim::set_data(): "
                      << "Full minimization." << std::endl;
          }
          
          // Create the simplex
          ubmatrix sx(np_covar+1,np_covar);
          ubvector sv(np_covar);
          
          if (use_alt_mmin==false) {
            
            if (this->verbose>0) {
              std::cout << "interpm_krige_optim::set_data(): "
                        << "Creating simplex." << std::endl;
            }
            
            for(size_t j=0;j<np_covar;j++) {
              sx(0,j)=plists[iout][j][plists[iout][j].size()/2];
            }
            for(size_t i=0;i<np_covar;i++) {
              for(size_t j=0;j<np_covar;j++) {
                if (i==j) {
                  sx(i+1,j)=plists[iout][j][plists[iout][j].size()-1];
                } else {
                  sx(i+1,j)=plists[iout][j][0];
                } 
              }
            }
            
            // Construct a maximum value to use if qual_fun() fails
            bool max_val_set=false;
            for(size_t i=0;i<np_covar+1;i++) {
              for(size_t j=0;j<np_covar;j++) {
                params[j]=sx(i,j);
              }
              cf[iout]->set_params(params);
              double qtmp=qual_fun(iout,success);
              if (success==0) {
                if (max_val_set==false || qtmp>max_val) {
                  max_val=qtmp;
                  max_val_set=true;
                }
              }
            }
            if (max_val_set==false) {
              O2SCL_ERR("Max val failed.",o2scl::exc_efailed);
            }
            
          } else {
            
            for(size_t j=0;j<np_covar;j++) {
              sv(j)=plists[iout][j][plists[iout][j].size()/2];
            }
            cf[iout]->set_params(sv);
            max_val=qual_fun(iout,success);
            if (success!=0) {
              O2SCL_ERR("Max val failed 2.",o2scl::exc_efailed);
            }
            
          }

          if (this->verbose>0) {
            std::cout << "interpm_krige_optim::set_data(): "
                      << "full minimization with\n  " << np_covar
                      << " parameters." << std::endl;
            if (use_alt_mmin==false) {
              std::cout << "  Simplex:" << std::endl;
              matrix_out(std::cout,np_covar+1,np_covar,sx,"  ");
            } else {
              std::cout << "  Initial point:" << std::endl;
              std::cout << "  ";
              vector_out(std::cout,sv,true);
            }
          }
          
          multi_funct mf=std::bind
            (std::mem_fn<double(size_t,size_t,const ubvector &,double)>
             (&class_t::min_fun),this,iout,
             std::placeholders::_1,std::placeholders::_2,max_val);

          if (use_alt_mmin==false) {
            int mret=def_mmin.mmin_simplex(np_covar,sx,min_qual,mf);
            if (mret!=0) {
              mret=def_mmin.mmin_simplex(np_covar,sx,min_qual,mf);
            }
            for(size_t j=0;j<np_covar;j++) {
              min_params[j]=sx(0,j);
            }
            if (false && mret!=0) {
              O2SCL_CONV_RET("Default minimizer failed in optim.",
                             o2scl::exc_einval,this->err_nonconv);
            }
          } else {
            int mret=alt_mmin.mmin(np_covar,sv,min_qual,mf);
            for(size_t j=0;j<np_covar;j++) {
              min_params[j]=sv[j];
            }
            if (mret!=0) {
              O2SCL_CONV_RET("Alternate minimizer failed in optim.",
                             o2scl::exc_einval,this->err_nonconv);
            }
          }
        
        } else {
          
          bool min_set=false, done=false;
          
          if (this->verbose>1) {
            std::cout << "interpm_krige_optim::set_data(): "
                      << "simple minimization with " << np_covar
                      << " parameters." << std::endl;
            for(size_t jk=0;jk<plists[iout].size();jk++) {
              std::cout << "  " << jk << " ";
              o2scl::vector_out(std::cout,plists[iout][jk],true);
            }
          }

          if (this->verbose>1) {
            std::cout << "interpm_krige_optim::set_data(): "
                      << "indexes params qual fail min_qual" << std::endl;
          }
          
          std::vector<size_t> index_list(np_covar);
          vector_set_all(np_covar,index_list,0);
          
          while (done==false) {
            
            for(size_t i=0;i<np_covar;i++) {
              params[i]=plists[iout][i][index_list[i]];
            }
            cf[iout]->set_params(params);
            
            qual[iout]=qual_fun(iout,success);
            
            if (success==0 && (min_set==false || qual[iout]<min_qual)) {
              min_params=params;
              min_qual=qual[iout];
              min_set=true;
            }
            
            if (this->verbose>1) {
              std::cout << "interpm_krige_optim::set_data(): ";
              o2scl::vector_out(std::cout,index_list);
              std::cout << " ";
              o2scl::vector_out(std::cout,params);
              std::cout << " " << qual[iout] << " "
                        << success << " " << min_qual << std::endl;
              if (this->verbose>2) {
                char ch;
                std::cin >> ch;
              }
            }
            
            index_list[0]++;
            for(size_t k=0;k<np_covar;k++) {
              if (index_list[k]==plists[iout][k].size()) {
                if (k==np_covar-1) {
                  done=true;
                } else {
                  index_list[k]=0;
                  index_list[k+1]++;
                }
              }
            }
            
          }
          

        }
        
        if (this->verbose>1) {
          std::cout << "interpm_krige_optim::set_data(): ";
          std::cout << "Minimum: " << min_qual << std::endl;
        }
    
        std::cout << "interpm_krige_optim::set_data(): "
                  << "Mode final: " << std::endl;
        cf[iout]->set_params(min_params);
        size_t mode_temp=mode;
        mode=mode_final;
        qual[iout]=qual_fun(iout,success);
        if (success!=0) {
          std::cout << "interpm_krige_optim::set_data(): "
                    << "Failed. Integer 'success' is "
                    << success << "." << std::endl;
          O2SCL_ERR2("Final calculation of Kinvf failed in ",
                     "interpm_krige_optim::set_data().",
                     o2scl::exc_efailed);
        }
        mode=mode_temp;
        std::cout << "interpm_krige_optim::set_data(): ";
        std::cout << "Mode final done: " << std::endl;
	
        if (this->verbose>0) {
          std::cout << "interpm_krige_optim::set_data():\n  "
                    << "Optimal parameters: ";
          o2scl::vector_out(std::cout,min_params,true);
        }
        
        if (timing) {
          t5=time(0);
          std::cout << "interpm_krige_optim::set_data(): ";
          std::cout << "Optimization of output " << iout << " took "
                    << t5-t4 << " seconds." << std::endl;
        }
        
        // End of loop over iout
        //std::cout << "Going to next iout." << std::endl;
      }

      return 0;
    }
    
    /** \brief Evaluate the interpolation at point \c x,
        returning \c y
    */
    virtual int eval(const vec_t &xp, vec_t &yp) const {
      return eval_tl(xp,yp);
    }

    /** \brief Evaluate the interpolation at point \c x,
        returning \c y and the uncertainties in \c y_unc
    */
    virtual int eval_unc(const vec_t &xp, vec_t &yp, vec_t &y_unc) const {
      int ret=eval(xp,yp);
      sigma(xp,y_unc);
      return ret;
    }
    
    /** \brief Given input vector \c x
        store the result of the interpolation in \c y
    */
    template<class vec4_t>
    int eval_tl(const vec_t &x0, vec4_t &y0) const {

      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_krige::eval_covar().",
                  exc_einval);
      }

      // Evaluate the interpolated result
      for(size_t iout=0;iout<this->n_outputs;iout++) {
        y0[iout]=0.0;
        for(size_t ipoints=0;ipoints<this->n_points;ipoints++) {
          mat_x_row_t xrow(x,ipoints);
          double covar_val=(*cf[iout])(x0,xrow);
          y0[iout]+=covar_val*Kinvf[iout][ipoints];
        }
        if (rescale) {
          y0[iout]*=std_y[iout];
          y0[iout]+=mean_y[iout];
        }
      }
      
      return 0;
    }
    //@}

    /// \name Uncertainties and derivatives
    //@{
    /** \brief Return the interpolation uncertainty from the 
        Gaussian process
    */
    template<class vec4_t>
    void sigma(const vec_t &x0, vec4_t &dy0) const {

      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_krige::sigma_covar().",
                  exc_einval);
      }
      if (!keep_matrix) {
        O2SCL_ERR2("Matrix information missing (keep_matrix==false) in ",
                   "interpm_krige::sigma_covar().",o2scl::exc_einval);
      }
      
      // Evaluate the interpolated result
      for(size_t iout=0;iout<this->n_outputs;iout++) {
        
        double kx0x0=cf[iout]->covar(x0,x0);
        
        vec_t kxx0(this->n_points), prod(this->n_points);
        
        for(size_t ipoints=0;ipoints<this->n_points;ipoints++) {
          mat_x_row_t xrow(x,ipoints);
          kxx0[ipoints]=(*cf[iout])(x0,xrow);
        }
        
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           this->n_points,this->n_points,1.0,
                           inv_KXX[iout],kxx0,0.0,prod);
        dy0[iout]=kx0x0-o2scl_cblas::ddot(this->n_points,kxx0,prod);
        
        if (rescale) {
          dy0[iout]*=std_y[iout];
        }
      }

      return;
    }
    
    /** \brief Given input vector \c x
        store the result of the interpolation in \c y
    */
    template<class vec2_t, class vec4_t>
    void deriv(const vec2_t &x0, vec4_t &y0, size_t ix) {

      // Evaluate the interpolated result
      for(size_t iout=0;iout<this->n_outputs;iout++) {
        y0[iout]=0.0;
        for(size_t ipoints=0;ipoints<this->n_points;ipoints++) {
          mat_x_row_t xrow(this->x,ipoints);
          double covar_val=(*cf[iout]).deriv(x0,xrow,ix);
          y0[iout]+=covar_val*this->Kinvf[iout][ipoints];
        }
        if (this->rescale) {
          y0[iout]*=this->std_y[iout];
        }
      }

      return;
    }

    /** \brief Given input vector \c x
        store the result of the interpolation in \c y
    */
    template<class vec2_t, class vec4_t>
    void deriv2(const vec2_t &x0, vec4_t &y0, size_t ix, size_t iy) {

      // Evaluate the interpolated result
      for(size_t iout=0;iout<this->n_outputs;iout++) {
        y0[iout]=0.0;
        for(size_t ipoints=0;ipoints<this->n_points;ipoints++) {
          mat_x_row_t xrow(this->x,ipoints);
          double covar_val=(*cf[iout]).deriv2(x0,xrow,ix,iy);
          y0[iout]+=covar_val*this->Kinvf[iout][ipoints];
        }
        if (this->rescale) {
          y0[iout]*=this->std_y[iout];
        }
      }

      return;
    }
    //@}

  };

#ifdef O2SCL_NEVER_DEFINED

  // AWS, 10/28/24, These don't work yet I'm not sure why...
  
  /// An Eigen specialization for \ref interpm_krige_optim
  typedef interpm_krige_optim
  <class vec_t=boost::numeric::ublas::vector<double>,
   class mat_x_t=o2scl::const_matrix_view_table<>,
   class mat_x_row_t=const const_matrix_row_gen
   <o2scl::const_matrix_view_table<>>, 
   class mat_y_t=o2scl::matrix_view_table<>,
   class mat_y_col_t=const matrix_column_gen<
     o2scl::matrix_view_table<>>,Eigen::MatrixXd,
   o2scl_linalg::matrix_invert_det_eigen<> > interpm_krige_optim_eigen;
  
  /// An Armadillo specialization for \ref interpm_krige_optim
  typedef interpm_krige_optim
  <class vec_t=boost::numeric::ublas::vector<double>,
   class mat_x_t=o2scl::const_matrix_view_table<>,
   class mat_x_row_t=const const_matrix_row_gen
   <o2scl::const_matrix_view_table<>>, 
   class mat_y_t=o2scl::matrix_view_table<>,
   class mat_y_col_t=const matrix_column_gen<
     o2scl::matrix_view_table<>>,arma::mat,
   o2scl_linalg::matrix_invert_det_sympd_arma<> > interpm_krige_optim_arma;

#endif
  
}
    
#endif



