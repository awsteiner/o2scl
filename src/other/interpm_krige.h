/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_INTERPM_KRIGE_H
#define O2SCL_INTERPM_KRIGE_H

/** \file interpm_krige.h
    \brief File defining \ref o2scl::interpm_krige,
    \ref o2scl::interpm_krige_optim and \ref o2scl::interpm_krige_optim
*/

#include <iostream>
#include <string>
#include <cmath>
#include <ctime>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <gsl/gsl_combination.h>

#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>
#include <o2scl/vec_stats.h>
#include <o2scl/linear_solver.h>
#include <o2scl/columnify.h>
#include <o2scl/cholesky.h>
#include <o2scl/min_brent_gsl.h>
#include <o2scl/mmin_simp2.h>
#include <o2scl/cblas.h>
#include <o2scl/invert.h>

namespace o2scl {

  /** \brief Covariance function: 1D RBF with a noise term

      \note There's no point making a base class, since there
      aren't really any virtual functions. The covariance functions
      have to be templates, to handle multiple vector types, so 
      no virtual functions are allowed.
   */
  class mcovar_funct_rbf {
    
  public:

    /// Length parameters
    std::vector<double> len;
    
    /// Noise parameter
    double noise;

    mcovar_funct_rbf() {
      noise=0.0;
    }
    
    /// Get the number of parameters
    size_t get_n_params(size_t ic) {
      return len.size();
    }
    
    /// Set the parameters
    template<class vec_t>
    void set_params(size_t ic, vec_t &p) {
      for(size_t j=0;j<len.size();j++) {
        len[j]=p[j];
      }
      return;
    }
    
    /// The covariance function
    template<class vec_t, class vec2_t>
    double operator()(size_t ic, const vec_t &x1, const vec2_t &x2) {
      double sum=0.0;
      bool equal=true;
      for(size_t j=0;j<len.size();j++) {
        if (x1[j]!=x2[j]) equal=false;
        sum+=-(x1[j]-x2[j])*(x1[j]-x2[j])/len[j]/len[j]/2.0;
      }
      if (equal) return exp(sum)+noise;
      return exp(sum);
    }

    /** \brief The derivative of the covariance function with
        respect to the first argument
    */
    template<class vec_t, class vec2_t>
    double deriv(size_t ic, vec_t &x1, vec2_t &x2, size_t ix) {
      double sum=0.0;
      for(size_t j=0;j<len.size();j++) {
        sum+=-(x1[j]-x2[j])*(x1[j]-x2[j])/len[j]/len[j]/2.0;
      }
      return -exp(sum)/len[ix]/len[ix]*(x1[ix]-x2[ix]);
    }
    
    /** \brief The second derivative of the covariance function with
        respect to the first argument
    */
    template<class vec_t, class vec2_t>
    double deriv2(size_t ic, vec_t &x1, vec2_t &x2, size_t ix,
                  size_t iy) {
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

  /** \brief Covariance function: 1D RBF with a noise term

      \note There's no point making a base class, since there
      aren't really any virtual functions. The covariance functions
      have to be templates, to handle multiple vector types, so 
      no virtual functions are allowed.
   */
  class mcovar_funct_rbf_noise {
    
  public:

    /// Length parameters
    std::vector<double> len;
    
    /// Noise parameter
    double log10_noise;
    
    /// Get the number of parameters
    size_t get_n_params(size_t ic) {
      return len.size()+1;
    }
    
    /// Set the parameters
    template<class vec_t>
    void set_params(size_t ic, vec_t &p) {
      for(size_t j=0;j<len.size();j++) {
        len[j]=p[j];
      }
      log10_noise=p[len.size()];
      return;
    }
    
    /// The covariance function
    template<class vec_t, class vec2_t>
    double operator()(size_t ic, const vec_t &x1, const vec2_t &x2) {
      double sum=0.0;
      bool equal=true;
      for(size_t j=0;j<len.size();j++) {
        if (x1[j]!=x2[j]) equal=false;
        sum+=-(x1[j]-x2[j])*(x1[j]-x2[j])/len[j]/len[j]/2.0;
      }
      if (equal) return exp(sum)+pow(10.0,log10_noise);
      return exp(sum);
    }

    /** \brief The derivative of the covariance function with
        respect to the first argument
    */
    template<class vec_t, class vec2_t>
    double deriv(size_t ic, vec_t &x1, vec2_t &x2, size_t ix) {
      double sum=0.0;
      for(size_t j=0;j<len.size();j++) {
        sum+=-(x1[j]-x2[j])*(x1[j]-x2[j])/len[j]/len[j]/2.0;
      }
      return -exp(sum)/len[ix]/len[ix]*(x1[ix]-x2[ix]);
    }

    /** \brief The second derivative of the covariance function with
        respect to the first argument
    */
    template<class vec_t, class vec2_t>
    double deriv2(size_t ic, vec_t &x1, vec2_t &x2, size_t ix,
                  size_t iy) {
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

  /** \brief Multi-dimensional interpolation by kriging

      \note The set data functions for this class uses a particular
      format, one different format than that in \ref
      o2scl::interpm_idw . This design choice makes it easier to pass
      vector arguments to the covariance function and the linear
      algebra routines. The x and y objects should be of the form
      <tt>x[n_points][n_in]</tt> and <tt>y[n_out][n_points]</tt>.
      A separate covariance function is required for each output.

      \note This class assumes that the function specified in the
      call to set_data() is the same as that passed to the
      eval() functions. If this is not the case, the
      behavior of this class is undefined.

      Because the data structure need not match the vector type of the
      input vectors, there are actually three different covariance
      functions (for the different pairwise combinations of the two
      vector types). In addition, different outputs may have different
      covariance functions. Thus there are three array types used for
      covariance functions. For convenience, these are typedef'd as
      f1_t, f2_t, and f3_t.
      
      \note Experimental.
  */
  template<class vec_t, class mat_x_t, class mat_x_row_t, 
           class mat_y_t, class mat_y_row_t, class mat_inv_kxx_t,
           class mat_inv_t=
           o2scl_linalg::matrix_invert_det_cholesky<mat_inv_kxx_t> >
  class interpm_krige {    
    
  public:

    /// The internal vector type for \f$ K^{-1} f \f$
    typedef boost::numeric::ublas::vector<double> ubvector;

    /// \name The three different ways the covariance function is used
    //@{
    typedef std::vector<std::function<double(mat_x_row_t &,
                                             mat_x_row_t &) > > f1_t;
    typedef std::vector<std::function<double(const vec_t &,
                                             const vec_t &) > > f3_t;
    //@}
    
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

  public:

    interpm_krige() {
      data_set=false;
      verbose=0;
      keep_matrix=true;
      rescaled=false;
    }

    virtual ~interpm_krige() {
    }
    
    /// If true, keep \f$ K^{-1} \f$ (default true)
    bool keep_matrix;
    
    /** \brief Verbosity parameter (default 0)
     */
    int verbose;
    
    /** \brief Initialize the data for the interpolation
    */
    template<class covar_func_t>
    int set_data_internal
    (size_t n_in, size_t n_out, size_t n_points, mat_x_t &user_x, 
     mat_y_t &user_y, covar_func_t &fcovar, bool rescale=false,
     bool err_on_fail=true) {

      if (n_points<2) {
        O2SCL_ERR2("Must provide at least two points in ",
                   "interpm_krige::set_data_internal()",
                   exc_efailed);
      }
      if (n_in<1) {
        O2SCL_ERR2("Must provide at least one input column in ",
                   "interpm_krige::set_data_internal()",
                   exc_efailed);
      }
      if (n_out<1) {
        O2SCL_ERR2("Must provide at least one output column in ",
                   "interpm_krige::set_data_internal()",
                   exc_efailed);
      }
      
      np=n_points;
      nd_in=n_in;
      nd_out=n_out;

      // Check that the data is properly sized
      if (user_x.size1()!=n_points || user_x.size2()!=n_in) {
        std::cout << "Object user_x, function size1() and size2(): "
                  << user_x.size1() << " " << user_x.size2() << std::endl;
        O2SCL_ERR2("Size of x not correct in ",
                   "interpm_krige::set_data_internal().",
                   o2scl::exc_efailed);
      }
    
      if (user_y.size2()!=n_points || user_y.size1()!=n_out) {
        std::cout << "Object user_y, function size1() and size2(): "
                  << user_y.size1() << " " << user_y.size2() << std::endl;
        O2SCL_ERR2("Size of y not correct in ",
                   "interpm_krige::set_data_internal().",
                   o2scl::exc_efailed);
      }

      std::swap(x,user_x);
      std::swap(y,user_y);
      
      rescaled=rescale;
      
      data_set=true;
    
      if (verbose>0) {
        std::cout << "interpm_krige::set_data_internal():\n  "
                  << "Using " << n_points
                  << " points with " << nd_in << " input variables and "
                  << nd_out << " output variables." << std::endl;
      }
      
      if (verbose>0 && rescale==false) {
        std::cout << "interpm_krige::set_data_internal(): "
                  << "No rescaling."
                  << std::endl;
      }
      
      if (rescale==true) {
        if (verbose>0) {
          std::cout << "interpm_krige::set_data_internal(): "
                    << "Rescaling."
                    << std::endl;
        }
        mean_y.resize(n_out);
        std_y.resize(n_out);
        for(size_t j=0;j<n_out;j++) {
          mat_y_row_t vec(y,j);
          mean_y[j]=vector_mean(n_points,vec);
          std_y[j]=vector_stddev(n_points,vec);
          if (verbose>1) {
            std::cout << "  Mean, std. dev. for y " << j << " of "
                      << n_out << " is "
                      << mean_y[j] << " " << std_y[j] << std::endl;
          }
          for(size_t i=0;i<n_points;i++) {
            y(j,i)=(y(j,i)-mean_y[j])/std_y[j];
          }
        }
      }

      Kinvf.resize(n_out);
      size_t n_covar=fcovar.size();

      inv_KXX.resize(n_out);
      
      // Loop over all output functions
      for(size_t iout=0;iout<n_out;iout++) {

        size_t icovar=iout % n_covar;

        // Select the row of the data matrix
        mat_y_row_t yiout(y,iout);

        // Construct the KXX matrix
        mat_inv_kxx_t KXX(n_points,n_points);
        for(size_t irow=0;irow<n_points;irow++) {
          mat_x_row_t xrow(x,irow);
          for(size_t icol=0;icol<n_points;icol++) {
            mat_x_row_t xcol(x,icol);
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else {
              KXX(irow,icol)=fcovar[icovar](xrow,xcol);
            }
          }
        }
        
        inv_KXX[iout].resize(n_points,n_points);
        mi.invert(n_points,KXX,inv_KXX[iout]);
        
        // Inverse covariance matrix times function vector
        Kinvf[iout].resize(n_points);
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           n_points,n_points,1.0,inv_KXX[iout],
                           yiout,0.0,Kinvf[iout]);
	
        if (verbose>1) {
          std::cout << "interpm_krige::set_data_internal() "
                    << "finished " << iout+1
                    << " of " << n_out << "." << std::endl;
        }
        
      }

      if (!keep_matrix) {
        inv_KXX.clear();
      }
      
      if (verbose>1) {
        std::cout << "interpm_krige::set_data_internal() done."
                  << std::endl;
      }
      
      return 0;
    }
    
    /** \brief Remove the rescaling of the data
     */
    void unscale() {
      
      if (rescaled==true) {
        for(size_t j=0;j<nd_out;j++) {
          for(size_t i=0;i<np;i++) {
            y(j,i)=y(j,i)*std_y[j]+mean_y[j];
          }
        }
        if (verbose>1) {
          std::cout << "interpm_krige::unscale(): "
                    << "returned to original values." 
                    << std::endl;
        }
      }
      
      return;
    }      

    /** \brief Restore the data to the user
     */
    template<class covar_func_t>
    void restore_data(mat_x_t &user_x, mat_y_t &user_y) {
      
      unscale();
      
      std::swap(x,user_x);
      std::swap(y,user_y);

      np=0;
      nd_in=0;
      nd_out=0;
      
      return;
    }
    
    /** \brief Initialize the data for the interpolation
    */
    template<class covar_func_t>
    int set_data(size_t n_in, size_t n_out, size_t n_points,
                 mat_x_t &user_x, mat_y_t &user_y,
                 covar_func_t &fcovar, bool rescale=false,
                 bool err_on_fail=true) {
      
      return set_data_internal
        (n_in,n_out,n_points,user_x,user_y,fcovar,
         rescale,err_on_fail);
    }

    /** \brief Given covariance function \c fcovar and input vector \c x0
        store the result of the interpolation in \c y0
    */
    template<class vec2_t, class vec3_t, class covar_func2_t>
    void eval_covar(const vec2_t &x0, vec3_t &y0, covar_func2_t &fcovar) {
      
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_krige::eval_covar().",
                  exc_einval);
      }

      // Evaluate the interpolated result
      for(size_t iout=0;iout<nd_out;iout++) {
        y0[iout]=0.0;
        for(size_t ipoints=0;ipoints<np;ipoints++) {
          mat_x_row_t xrow(x,ipoints);
          double covar_val=fcovar(iout,xrow,x0);
          y0[iout]+=covar_val*Kinvf[iout][ipoints];
        }
        if (rescaled) {
          y0[iout]*=std_y[iout];
          y0[iout]+=mean_y[iout];
        }
      }
      
      return;
      
    }
    
    /** \brief Given covariance function \c fcovar and input vector \c
        x store the result of the interpolation in \c y
    */
    template<class vec2_t, class vec3_t, class covar_func2_t,
      class covar_func3_t>
    void sigma_covar(const vec2_t &x0, vec3_t &dy0, covar_func2_t &f2,
                     covar_func3_t &f3) {
      
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_krige::sigma_covar().",
                  exc_einval);
      }
      if (!keep_matrix) {
        O2SCL_ERR2("Matrix information missing (keep_matrix==false) in ",
                   "interpm_krige::sigma_covar().",o2scl::exc_einval);
      }
      
      // Evaluate the interpolated result
      for(size_t iout=0;iout<nd_out;iout++) {
        double kx0x0=f3(iout,x0,x0);
        
        vec_t kxx0(np), prod(np);
        
        for(size_t ipoints=0;ipoints<np;ipoints++) {
          mat_x_row_t xrow(x,ipoints);
          kxx0[ipoints]=f2(iout,xrow,x0);
        }
        
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           np,np,1.0,inv_KXX[iout],kxx0,0.0,prod);
        dy0[iout]=kx0x0-o2scl_cblas::ddot(np,kxx0,prod);
        
        if (rescaled) {
          dy0[iout]*=std_y[iout];
        }
      }

      return;
      
    }
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
    /// The number of points
    size_t np;
    /// The number of dimensions of the inputs
    size_t nd_in;
    /// The number of dimensions of the outputs
    size_t nd_out;
    /// The data
    mat_x_t x;
    /// The data
    mat_y_t y;
    /// True if the data has been specified
    bool data_set;
    /// Desc
    ubvector mean_y;
    /// Desc
    ubvector std_y;
    /// True if the data needs to be rescaled
    bool rescaled;
  
#endif
  
  };

  /** \brief Multi-dimensional interpolation using an 
      optimized covariance function

      \verbatim embed:rst
      See also the :ref:`Higher-dimensional Interpolation` 
      section of the User's guide. 
      \endverbatim
  */
  template<class func_t, class vec_t, class mat_x_t, class mat_x_row_t, 
           class mat_y_t, class mat_y_row_t, class mat_inv_kxx_t,
           class mat_inv_t=
           o2scl_linalg::matrix_invert_det_cholesky<mat_inv_kxx_t>,
           class vec_vec_t=std::vector<std::vector<double>> >
  class interpm_krige_optim :
    public interpm_krige<vec_t,mat_x_t,mat_x_row_t,mat_y_t,
                         mat_y_row_t,mat_inv_kxx_t,
                         mat_inv_t> {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    
  protected:

    /// Pointer to the covariance function
    func_t *cf;
    
    /// List of parameter values to try
    vec_vec_t plists;
  
     /// The quality factor of the optimization for each output function
    std::vector<double> qual;

    /// Pointer to the user-specified minimizer
    mmin_base<> *mp;

  public:
    
    /** \brief Function to optimize the covariance parameters
     */
    virtual double qual_fun(size_t iout, int &success) {

      // Select the row of the data matrix
      mat_y_row_t yiout2(this->y,iout);
      
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
                inv_KXX2(irow,icol)=(*cf)(iout,xrow,xcol);
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
            kxx0[i]=(*cf)(iout,xi2,xk);
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
            if (irow>icol) {
              this->inv_KXX[iout](irow,icol)=this->inv_KXX[iout](icol,irow);
            } else {
              this->inv_KXX[iout](irow,icol)=(*cf)(iout,xrow,xcol);
            }
          }
        }

        if (verbose>2) {
          std::cout << "Done creating covariance matrix with size "
                    << size << std::endl;
        }

        if (timing) {
          t1=time(0);          
        }
        
        // Construct the inverse of KXX
        if (verbose>2) {
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
          std::cout << "Matrix inversion took "
                    << t2-t1 << " seconds." << std::endl;
        }
        
        if (verbose>2) {
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
          std::cout << "Matrix vector multiply took "
                    << t3-t2 << " seconds." << std::endl;
        }
        
        ret=0.0;
        
        // Select the row of the data matrix
        mat_y_row_t yiout(this->y,iout);
        
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
          std::cout << "Final evaluation took "
                    << t4-t3 << " seconds." << std::endl;
        }
        
        if (verbose>2) {
          std::cout << "ret: " << ret << std::endl;
        }
      
      } else if (mode==mode_max_lml || mode==mode_final) {

        if (verbose>2) {
          std::cout << "Creating covariance matrix with size "
                    << size << std::endl;
        }

        // Construct the KXX matrix
        mat_inv_kxx_t KXX(size,size);
        for(size_t irow=0;irow<size;irow++) {
          mat_x_row_t xrow(this->x,irow);
          for(size_t icol=0;icol<size;icol++) {
            mat_x_row_t xcol(this->x,icol);
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else {
              KXX(irow,icol)=(*cf)(iout,xrow,xcol);
            }
          }
        }

        if (verbose>2) {
          std::cout << "Done creating covariance matrix with size "
                    << size << std::endl;
        }

        // Perform the matrix inversion and compute the determinant

        double lndet;
	
        if (timing) {
          t1=time(0);          
        }
        
        // Construct the inverse of KXX
        if (verbose>2) {
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
	
        if (verbose>2) {
          std::cout << "Done performing matrix inversion with size "
                    << size << std::endl;
        }
        
        if (timing) {
          t2=time(0);          
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
          std::cout << "Final evaluation took "
                    << t4-t3 << " seconds." << std::endl;
        }
        
      }

      if (!isfinite(ret)) success=3;
      
      return ret;
    }

    interpm_krige_optim() {
      full_min=false;
      mp=&def_mmin;
      verbose=0;
      mode=mode_loo_cv;
      loo_npts=100;
      timing=false;
    }

    virtual ~interpm_krige_optim() {
    }
    
    /// \name Function to minimize and various option
    //@{
    /// Leave-one-out cross validation
    static const size_t mode_loo_cv=1;
    /// Minus Log-marginal-likelihood
    static const size_t mode_max_lml=2;
    static const size_t mode_loo_cv_bf=3;
    /// No optimization (for internal use)
    static const size_t mode_final=10;
    /// Function to minimize (default \ref mode_loo_cv)
    size_t mode;
    /// Number of points to test for cross validation (default 100)
    size_t loo_npts;
    ///@}
    
    /// Verbosity parameter
    int verbose;
    
    /// Default minimizer
    mmin_simp2<> def_mmin;

    /** \brief If true, output timing results
     */
    bool timing;
    
    /// If true, use the full minimizer
    bool full_min;
  
    /** \brief Set the covariance function and parameter lists
     */
    int set_covar(func_t &covar, vec_vec_t &param_lists) {
      cf=&covar;
      plists=param_lists;
      return 0;
    }

    /** \brief Initialize interpolation routine
     */
    int set_data_internal
    (size_t n_in, size_t n_out, size_t n_points,
     mat_x_t &user_x, mat_y_t &user_y, 
     bool rescale=false, bool err_on_fail=true) {

      if (n_points<2) {
        O2SCL_ERR2("Must provide at least two points in ",
                   "interpm_krige_optim::set_data_internal()",
                   exc_efailed);
      }
      if (n_in<1) {
        O2SCL_ERR2("Must provide at least one input column in ",
                   "interpm_krige_optim::set_data_internal()",
                   exc_efailed);
      }
      if (n_out<1) {
        O2SCL_ERR2("Must provide at least one output column in ",
                   "interpm_krige_optim::set_data_internal()",
                   exc_efailed);
      }
   
      // Check that the data is properly sized
      if (user_x.size1()!=n_points || user_x.size2()!=n_in) {
        std::cout << "Object user_x, function size1() and size2(): "
                  << user_x.size1() << " " << user_x.size2() << std::endl;
        O2SCL_ERR2("Size of x not correct in ",
                   "interpm_krige::set_data_internal().",
                   o2scl::exc_efailed);
      }
    
      if (user_y.size2()!=n_points || user_y.size1()!=n_out) {
        std::cout << "Object user_y, function size1() and size2(): "
                  << user_y.size1() << " " << user_y.size2() << std::endl;
        O2SCL_ERR2("Size of y not correct in ",
                   "interpm_krige::set_data_internal().",
                   o2scl::exc_efailed);
      }

      // Set parent data members
      this->np=n_points;
      this->nd_in=n_in;
      this->nd_out=n_out;
      this->rescaled=rescale;
      this->data_set=true;
      
      time_t t1=0, t2=0, t3=0, t4=0, t5=0;
      if (timing) {
        t1=time(0);
      }
      
      std::swap(this->x,user_x);
      std::swap(this->y,user_y);

      if (timing) {
        t2=time(0);
        std::cout << "Swap took " << t2-t1 << " seconds." << std::endl;
      }
      
      if (verbose>0) {
        std::cout << "interpm_krige_optim::set_data_internal(): "
                  << "Using " << n_points
                  << " points with\n " << n_in << " input variables and "
                  << this->nd_out << " output variables." << std::endl;
      }

      if (rescale==true) {
        this->mean_y.resize(n_out);
        this->std_y.resize(n_out);
        for(size_t j=0;j<n_out;j++) {
          mat_y_row_t vec(this->y,j);
          this->mean_y[j]=vector_mean(n_points,vec);
          this->std_y[j]=vector_stddev(n_points,vec);
          if (verbose>1) {
            std::cout << "Mean,stddev of y " << j << " of "
                      << n_out << " is " << this->mean_y[j] << " "
                      << this->std_y[j] << std::endl;
          }
          for(size_t i=0;i<n_points;i++) {
            this->y(j,i)=(this->y(j,i)-this->mean_y[j])/this->std_y[j];
          }
        }
        if (verbose>1) {
          std::cout << "interpm_krige_optim::set_data_internal(): "
                    << "data rescaled." << std::endl;
        }
      }

      if (timing) {
        t3=time(0);
        std::cout << "Rescale took " << t3-t2 << " seconds." << std::endl;
      }
      
      int success=0;

      this->Kinvf.resize(n_out);
      this->inv_KXX.resize(n_out);

      qual.resize(n_out);

      // Loop over all output functions
      for(size_t iout=0;iout<n_out;iout++) {
        
        if (timing) {
          t4=time(0);
        }
        
        // Initialize to zero to prevent uninit'ed var. warnings
        double min_qual=0.0;
	
        if (verbose>1) {
          std::cout << "qual fail min_qual" << std::endl;
        }
	
        bool min_set=false, done=false;
        
        size_t np_covar=cf->get_n_params(iout);
        if (verbose>1) {
          std::cout << "interpm_krige_optim::set_data_internal() : "
                    << "simple minimization with " << np_covar
                    << " parameters." << std::endl;
          for(size_t jk=0;jk<plists.size();jk++) {
            std::cout << jk << " ";
            o2scl::vector_out(std::cout,plists[jk],true);
          }
        }
        std::vector<size_t> index_list(np_covar);
        vector_set_all(np_covar,index_list,0);
        std::vector<double> params(np_covar), min_params(np_covar);
        
        while (done==false) {
          
          for(size_t i=0;i<np_covar;i++) {
            params[i]=plists[i][index_list[i]];
          }
          cf->set_params(iout,params);
          
          qual[iout]=qual_fun(iout,success);
          
          if (success==0 && (min_set==false || qual[iout]<min_qual)) {
            min_params=params;
            min_qual=qual[iout];
            min_set=true;
          }
          
          if (verbose>1) {
            std::cout << "interpm_krige_optim: ";
            o2scl::vector_out(std::cout,index_list);
            std::cout << " ";
            o2scl::vector_out(std::cout,params);
            std::cout << " " << qual[iout] << " "
                      << success << " " << min_qual << std::endl;
            if (verbose>2) {
              char ch;
              std::cin >> ch;
            }
          }

          index_list[0]++;
          for(size_t k=0;k<np_covar;k++) {
            if (index_list[k]==plists[k].size()) {
              if (k==np_covar-1) {
                done=true;
              } else {
                index_list[k]=0;
                index_list[k+1]++;
              }
            }
          }
          
        }
        
        if (verbose>1) {
          std::cout << "interpm_krige_optim: ";
          std::cout.width(2);
          std::cout << "   " << min_qual << std::endl;
        }
        cf->set_params(iout,min_params);
        size_t mode_temp=mode;
        mode=mode_final;
        qual[iout]=qual_fun(iout,success);
        mode=mode_temp;
	
        if (verbose>0) {
          std::cout << "interpm_krige_optim::set_data_"
                    << "internal(),\n  "
                    << "optimal parameters: ";
          o2scl::vector_out(std::cout,min_params,true);
        }
        
        if (timing) {
          t5=time(0);
          std::cout << "Optimization of output " << iout << " took "
                    << t5-t4 << " seconds." << std::endl;
        }
        
        // End of loop over iout
      }

      return 0;
    }

    /** \brief Given input vector \c x
        store the result of the interpolation in \c y
    */
    template<class vec2_t, class vec3_t>
    void eval(const vec2_t &x0, vec3_t &y0) {

      return this->eval_covar(x0,y0,*cf);

      return;
    }

    /** \brief Return the interpolation uncertainty from the 
        Gaussian process
    */
    template<class vec2_t, class vec3_t>
    void sigma(const vec2_t &x0, vec3_t &y0) {

      return this->sigma_covar(x0,y0,*cf);

    }
    
    /** \brief Given input vector \c x
        store the result of the interpolation in \c y
    */
    template<class vec2_t, class vec3_t>
    void deriv(const vec2_t &x0, vec3_t &y0, size_t ix) {

      // Evaluate the interpolated result
      for(size_t iout=0;iout<this->nd_out;iout++) {
        y0[iout]=0.0;
        for(size_t ipoints=0;ipoints<this->np;ipoints++) {
          mat_x_row_t xrow(this->x,ipoints);
          double covar_val=cf->deriv(iout,xrow,x0,ix);
          y0[iout]-=covar_val*this->Kinvf[iout][ipoints];
        }
        if (this->rescaled) {
          y0[iout]*=this->std_y[iout];
        }
      }

      return;
    }

    /** \brief Given input vector \c x
        store the result of the interpolation in \c y
    */
    template<class vec2_t, class vec3_t>
    void deriv2(const vec2_t &x0, vec3_t &y0, size_t ix, size_t iy) {

      // Evaluate the interpolated result
      for(size_t iout=0;iout<this->nd_out;iout++) {
        y0[iout]=0.0;
        for(size_t ipoints=0;ipoints<this->np;ipoints++) {
          mat_x_row_t xrow(this->x,ipoints);
          double covar_val=cf->deriv2(iout,xrow,x0,ix,iy);
          y0[iout]+=covar_val*this->Kinvf[iout][ipoints];
        }
        if (this->rescaled) {
          y0[iout]*=this->std_y[iout];
        }
      }

      return;
    }

    /** \brief Initialize the data for the interpolation
    */
    int set_data(size_t n_in, size_t n_out, size_t n_points,
                 mat_x_t &user_x, mat_y_t &user_y, 
                 bool rescale=false, bool err_on_fail=true) {
      
      return set_data_internal
        (n_in,n_out,n_points,user_x,
         user_y,rescale,err_on_fail);
    }

  };
  
}
    
#endif



