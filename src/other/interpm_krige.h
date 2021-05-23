/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2021, Andrew W. Steiner
  
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
    \ref o2scl::interpm_krige_nn and \ref o2scl::interpm_krige_optim
*/

#include <iostream>
#include <string>
#include <cmath>

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
#include <o2scl/cblas.h>
#include <o2scl/invert.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

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

      \note Experimental.
  */
  template<class vec_t, class mat_t, class mat_row_t,
           class mat_col_t,
           class mat2_t, class mat2_row_t, class mat3_t,
           class mat_inv_t=o2scl_linalg::matrix_invert_det_cholesky<mat3_t> >
  class interpm_krige {    
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;
    typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
    
  protected:

    /** \brief Inverse covariance matrix times function vector
     */
    std::vector<ubvector> Kinvf;
    
    /// The inverse of the covariance matrix
    mat3_t inv_KXX;

    /// The matrix inversion object
    mat_inv_t mi;
    
  public:

    interpm_krige() {
      data_set=false;
      verbose=0;
      keep_matrix=true;
    }
    
    /// If true, keep \f$ K^{-1} \f$
    bool keep_matrix;
    
    /** \brief Verbosity parameter (default 0)
     */
    int verbose;
    
    /** \brief Initialize the data for the interpolation

        \note This function works differently than 
        \ref o2scl::interpm_idw::set_data() . See this
        class description for more details.
    */
    template<class func_vec_t>
    int set_data_noise(size_t n_in, size_t n_out, size_t n_points,
                       mat_t &user_x, mat2_t &user_y,
                       func_vec_t &fcovar,
                       const vec_t &noise_var, bool rescale=false,
                       bool err_on_fail=true) {

      if (n_points<2) {
        O2SCL_ERR2("Must provide at least two points in ",
                   "interpm_krige::set_data_noise()",exc_efailed);
      }
      if (n_in<1) {
        O2SCL_ERR2("Must provide at least one input column in ",
                   "interpm_krige::set_data_noise()",exc_efailed);
      }
      if (n_out<1) {
        O2SCL_ERR2("Must provide at least one output column in ",
                   "interpm_krige::set_data_noise()",exc_efailed);
      }
      if (noise_var.size()<1) {
        O2SCL_ERR2("Noise vector empty in ",
                   "interpm_krige::set_data_noise()",exc_efailed);
      }
      np=n_points;
      nd_in=n_in;
      nd_out=n_out;
    
      if (user_x.size1()!=n_points || user_x.size2()!=n_in) {
        
        std::cout << user_x.size1() << std::endl;
        std::cout << user_x.size2() << std::endl;
        
        O2SCL_ERR2("Size of x not correct in ",
                   "interpm_krige::set_data_noise().",o2scl::exc_efailed);
      }
    
      // Check that the data is properly sized
      if (user_y.size2()!=n_points || user_y.size1()!=n_out) {
        std::cout << user_y.size1() << std::endl;
        std::cout << user_y.size2() << std::endl;
      
        O2SCL_ERR2("Size of y not correct in ",
                   "interpm_krige::set_data_noise().",o2scl::exc_efailed);
      }

      x=&user_x;
      y=&user_y;
      
      rescaled=rescale;
    
      data_set=true;
    
      if (verbose>0) {
        std::cout << "interpm_krige::set_data_noise() : Using " << n_points
                  << " points with " << nd_in << " input variables and\n\t"
                  << nd_out << " output variables." << std::endl;
      }
    
      if (rescale==true) {
        mean_x.resize(n_in);
        std_x.resize(n_in);
        for(size_t j=0;j<n_in;j++) {
          mat_col_t vec(*x,j);
          mean_x[j]=vector_mean(n_points,vec);
          std_x[j]=vector_stddev(n_points,vec);
          if (verbose>1) {
            std::cout << "Mean,stddev of x " << j << " of " << n_in << " is "
                      << mean_x[j] << " " << std_x[j] << std::endl;
          }
          for(size_t i=0;i<n_points;i++) {
            user_x(i,j)=(user_x(i,j)-mean_x[j])/std_x[j];
          }
        }
        mean_y.resize(n_out);
        std_y.resize(n_out);
        for(size_t j=0;j<n_out;j++) {
          mat2_row_t vec(*y,j);
          mean_y[j]=vector_mean(n_points,vec);
          std_y[j]=vector_stddev(n_points,vec);
          if (verbose>1) {
            std::cout << "Mean,stddev of y " << j << " of " << n_out << " is "
                      << mean_y[j] << " " << std_y[j] << std::endl;
          }
          for(size_t i=0;i<n_points;i++) {
            user_y(j,i)=(user_y(j,i)-mean_y[j])/std_y[j];
          }
        }
        if (verbose>1) {
          std::cout << "interpm_krige::set_data_noise() rescaled."
                    << std::endl;
        }
      }

      Kinvf.resize(n_out);
      size_t n_covar=fcovar.size();

      // Loop over all output functions
      for(size_t iout=0;iout<n_out;iout++) {

        size_t icovar=iout % n_covar;
        size_t inoise=iout & noise_var.size();

        // Select the row of the data matrix
        mat2_row_t yiout(*y,iout);

        // Construct the KXX matrix
        mat3_t KXX(n_points,n_points);
        for(size_t irow=0;irow<n_points;irow++) {
          mat_row_t xrow(*x,irow);
          for(size_t icol=0;icol<n_points;icol++) {
            mat_row_t xcol(*x,icol);
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else if (irow==icol) {
              KXX(irow,icol)=fcovar[icovar](xrow,xcol)+
                noise_var[inoise];
            } else {
              KXX(irow,icol)=fcovar[icovar](xrow,xcol);
            }
          }
        }

        inv_KXX.resize(n_points,n_points);
        mi.invert(n_points,KXX,inv_KXX);
        
        // Inverse covariance matrix times function vector
        Kinvf[iout].resize(n_points);
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           n_points,n_points,1.0,inv_KXX,
                           yiout,0.0,Kinvf[iout]);
	
        if (verbose>1) {
          std::cout << "interpm_krige::set_data_noise() finished " << iout+1
                    << " of " << n_out << "." << std::endl;
        }
        
      }
      
      if (verbose>1) {
        std::cout << "interpm_krige::set_data_noise() done."
                  << std::endl;
      }
      
      return 0;
    }

    void unscale(size_t n_in, size_t n_out, size_t n_points) {
      if (rescaled==true) {
        for(size_t j=0;j<n_in;j++) {
          for(size_t i=0;i<n_points;i++) {
            (*x)(i,j)=(*x)(i,j)*std_x[j]+mean_x[j];
          }
        }
        for(size_t j=0;j<n_out;j++) {
          for(size_t i=0;i<n_points;i++) {
            (*y)(j,i)=(*y)(j,i)*std_y[j]+mean_y[j];
          }
        }
        if (verbose>1) {
          std::cout << "interpm_krige::set_data_noise() "
                    << "returned to original values." 
                    << std::endl;
        }
      }
      return;
    }      

    /** \brief Initialize the data for the interpolation
      
        \note This function works differently than 
        \ref o2scl::interpm_idw::set_data() . See this
        class description for more details.
    */
    template<class func_vec_t>
    int set_data(size_t n_in, size_t n_out, size_t n_points,
                 mat_t &user_x, mat2_t &user_y,
                 func_vec_t &fcovar, bool rescale=false,
                 bool err_on_fail=true) {
      vec_t noise_vec;
      noise_vec.resize(1);
      noise_vec[0]=0.0;
      return set_data_noise<func_vec_t>
        (n_in,n_out,n_points,user_x,user_y,fcovar,
         noise_vec,rescale,err_on_fail);
    }

    /** \brief Given covariance function \c fcovar and input vector \c x
        store the result of the interpolation in \c y
    */
    template<class vec2_t, class vec3_t, class vec_func_t>
    void eval(const vec2_t &x0, vec3_t &y0, vec_func_t &fcovar) {
    
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_krige::eval().",
                  exc_einval);
      }
      if (fcovar.size()==0) {
        O2SCL_ERR("No covariance functions in interpm_krige::eval().",
                  exc_einval);
      }

      if (rescaled) {
        
        // If necessary, rescale before evaluating the interpolated
        // result
        vec2_t x0p(nd_in);
        for(size_t iin=0;iin<nd_in;iin++) {
          x0p[iin]=(x0[iin]-mean_x[iin])/std_x[iin];
        }

        // Evaluate the interpolated result
        for(size_t iout=0;iout<nd_out;iout++) {
          size_t icovar=iout % fcovar.size();
          y0[iout]=0.0;
          for(size_t ipoints=0;ipoints<np;ipoints++) {
            mat_row_t xrow(*x,ipoints);
            double covar_val=fcovar[icovar](xrow,x0p);
            y0[iout]+=covar_val*Kinvf[iout][ipoints];
          }
          y0[iout]*=std_y[iout];
          y0[iout]+=mean_y[iout];
        }

      } else {
      
        // Evaluate the interpolated result
        for(size_t iout=0;iout<nd_out;iout++) {
          size_t icovar=iout % fcovar.size();
          y0[iout]=0.0;
          for(size_t ipoints=0;ipoints<np;ipoints++) {
            mat_row_t xrow(*x,ipoints);
            y0[iout]+=fcovar[icovar](xrow,x0)*Kinvf[iout][ipoints];
          }
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
    mat_t* x;
    /// The data
    mat2_t* y;
    /// True if the data has been specified
    bool data_set;
    /// Desc
    ubvector mean_x;
    /// Desc
    ubvector std_x;
    /// Desc
    ubvector mean_y;
    /// Desc
    ubvector std_y;
    /// True if the data needs to be rescaled
    bool rescaled;
  
#endif
  
  };
  
#ifdef O2SCL_NEVER_DEFINED

  /** \brief One-dimensional interpolation using an 
      optimized covariance function
      
      \verbatim embed:rst
      See also the :ref:`Higher-dimensional Interpolation` 
      section of the User's guide. 
      \endverbatim

      \note This class is experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_t=boost::numeric::ublas::vector<double>,
           class mat_row_t=boost::numeric::ublas::matrix_row
           <boost::numeric::ublas::vector<double> > >
  class interpm_krige_optim :
    public interpm_krige<vec_t,mat_t,mat_row_t> {    

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

    /// Function objects for the covariance
    std::vector<std::function<double(const mat_row_t &, const vec_t &)> >
    ff2;
  
  protected:

    /// Function objects for the covariance
    std::vector<std::function<double(const mat_row_t &, const mat_row_t &)> >
    ff1;
  
    /// The covariance function length scale for each output function
    std::vector<double> len;
  
    /// The quality factor of the optimization for each output function
    std::vector<double> qual;

    /// If true, min and max has been set for the length parameter
    bool len_guess_set;

    /// Minimum for length parameter range
    double len_min;

    /// Maximum for length parameter range
    double len_max;
  
    /// The covariance function
    template<class vec2_t, class vec3_t>
    double covar(const vec2_t &x1, const vec3_t &x2, size_t sz, double len2) {
      double ret=0.0;
      for(size_t i=0;i<sz;i++) {
        ret+=pow(x1[i]-x2[i],2.0);
      }
      ret=exp(-ret/len2/len2/2.0);
      return ret;
    }

    /// Pointer to the user-specified minimizer
    min_base<> *mp;
  
    /** \brief Function to optimize the covariance parameters
     */
    template<class vec3_t> 
    double qual_fun(double xlen, double noise_var, size_t iout,
                    vec3_t &y, int &success) {

      double ret=0.0;
    
      success=0;

      size_t size=this->x.size1();

      if (mode==mode_loo_cv) {

        for(size_t ell=0;ell<loo_npts;ell++) {

          // Create the new data objects, x_jk and y_jk
          size_t row=ell*size/loo_npts;
          matrix_view_omit_row<mat_t> x_jk(this->x,row);
          ubvector y_jk(size-1);
          vector_copy_jackknife(size,y,row,y_jk);

          // Now perform the matrix analysis with those objects

          // Construct the KXX matrix
          mat3_t KXX(size-1,size-1);
          for(size_t irow=0;irow<size-1;irow++) {
            matrix_row_gen<matrix_view_omit_row<mat_t> > xrow(x_jk,irow);
            for(size_t icol=0;icol<size-1;icol++) {
              matrix_row_gen<matrix_view_omit_row<mat_t> > xcol(x_jk,icol);
              if (irow>icol) {
                KXX(irow,icol)=KXX(icol,irow);
              } else {
                KXX(irow,icol)=covar<
                  matrix_row_gen<matrix_view_omit_row<mat_t> >,
                  matrix_row_gen<matrix_view_omit_row<mat_t> > >
                  (xrow,xcol,this->nd_in,xlen);
                if (irow==icol) KXX(irow,icol)+=noise_var;
              }
            }
          }
          
          // Construct the inverse of KXX
          if (verbose>2) {
            std::cout << "Performing matrix inversion with size "
                      << size-1 << std::endl;
          }
          int cret=mi.invert_inplace(size-1,KXX);
          if (cret!=0) {
            success=1;
            return 1.0e99;
          }
	  
          // Inverse covariance matrix times function vector
          this->Kinvf[iout].resize(size-1);
          o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                             o2scl_cblas::o2cblas_NoTrans,
                             size-1,size-1,1.0,KXX,
                             y,0.0,this->Kinvf[iout]);
          
          double ypred=0.0;
          double yact=y[row];
          for(size_t i=0;i<size-1;i++) {
            matrix_row_gen<matrix_view_omit_row<mat_t> > xrow(x_jk,i);
            mat_row_t xcol(this->x,row);
            ypred+=covar<matrix_row_gen<matrix_view_omit_row<mat_t> >,
                         mat_row_t>(xrow,xcol,this->nd_in,xlen)*
              this->Kinvf[iout][i];
          }

          std::cout << "act,pred: " << yact << " " << ypred << std::endl;
        
          // Measure the quality with a chi-squared like function
          ret+=pow(yact-ypred,2.0);

	
          // Proceed to next point to omit
        }
        std::cout << "ret: " << ret << std::endl;
      
      } else if (mode==mode_max_lml || mode==mode_final) {

        if (verbose>2) {
          std::cout << "Creating covariance matrix with size "
                    << size << std::endl;
        }

        // Construct the KXX matrix
        mat3_t KXX(size,size);
        for(size_t irow=0;irow<size;irow++) {
          mat_row_t xrow(this->x,irow);
          for(size_t icol=0;icol<size;icol++) {
            mat_row_t xcol(this->x,icol);
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else {
              KXX(irow,icol)=covar<mat_row_t,mat_row_t>(xrow,xcol,
                                                        this->nd_in,xlen);
              if (irow==icol) KXX(irow,icol)+=noise_var;
            }
          }
        }
      
        // Note: We have to use LU here because O2scl doesn't yet
        // have a lndet() function for Cholesky decomp

        double lndet;
	
        // Construct the inverse of KXX
        if (verbose>2) {
          std::cout << "Performing matrix inversion with size "
                    << size << std::endl;
        }
        int cret=mi.invert_det(size,KXX,this->inv_KXX,lndet);
        if (cret!=0) {
          success=1;
          return 1.0e99;
        }
	
        lndet=log(lndet);
	
        // Inverse covariance matrix times function vector
        this->Kinvf[iout].resize(size);
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           size,size,1.0,KXX,
                           y,0.0,this->Kinvf[iout]);
	
        if (mode==mode_max_lml) {
          // Compute the log of the marginal likelihood, without
          // the constant term
          for(size_t i=0;i<size;i++) {
            ret+=0.5*y[i]*this->Kinvf[iout][i];
          }
          ret+=0.5*lndet;
        }

      }

      return ret;
    }

  public:

    interpm_krige_optim() {
      nlen=20;
      full_min=false;
      mp=&def_min;
      verbose=0;
      mode=mode_loo_cv;
      loo_npts=100;
      len_guess_set=false;
    }

    /// \name Function to minimize and various option
    //@{
    /// Leave-one-out cross validation
    static const size_t mode_loo_cv=1;
    /// Minus Log-marginal-likelihood
    static const size_t mode_max_lml=2;
    /// No optimization (for internal use)
    static const size_t mode_final=10;
    /// Function to minimize (default \ref mode_loo_cv)
    size_t mode;
    /// Number of points to test for cross validation (default 100)
    size_t loo_npts;
    ///@}
    
    /// Verbosity parameter
    int verbose;
  
    /** \brief Number of length scale points to try when full minimizer 
        is not used (default 20)
    */
    size_t nlen;

    /// Default minimizer
    min_brent_gsl<> def_min;

    /// If true, use the full minimizer
    bool full_min;
  
    /** \brief Set the range for the length parameter
     */
    void set_len_range(double min2, double max2) {
      len_guess_set=true;
      len_min=min2;
      len_max=max2;
      return;
    }
  
    /// Initialize interpolation routine
    template<class mat2_row_t, class mat2_t, class vec2_t,
             class vec3_t>
    int set_data_noise(size_t n_in, size_t n_out, size_t n_points,
                       mat_t &user_x, mat2_t &user_y, 
                       const vec2_t &noise_var, const vec3_t &len_precompute,
                       bool rescale=false, bool err_on_fail=true) {

      if (n_points<2) {
        O2SCL_ERR2("Must provide at least two points in ",
                   "interpm_krige_optim::set_data_noise()",exc_efailed);
      }
      if (n_in<1) {
        O2SCL_ERR2("Must provide at least one input column in ",
                   "interpm_krige_optim::set_data_noise()",exc_efailed);
      }
      if (n_out<1) {
        O2SCL_ERR2("Must provide at least one output column in ",
                   "interpm_krige_optim::set_data_noise()",exc_efailed);
      }
      if (noise_var.size()<1) {
        O2SCL_ERR2("Noise vector empty in ",
                   "interpm_krige::set_data_noise()",exc_efailed);
      }
   
      // Set parent data members
      this->np=n_points;
      this->nd_in=n_in;
      this->nd_out=n_out;
      std::swap(this->x,user_x);
      this->rescaled=rescale;
      this->data_set=true;
       
      if (verbose>0) {
        std::cout << "interpm_krige_optim::set_data_noise() : Using "
                  << n_points
                  << " points with\n " << n_in << " input variables and "
                  << this->nd_out << " output variables." << std::endl;
      }

      // Get max and min of all parameters
      this->min.resize(n_in);
      this->max.resize(n_in);
      for(size_t j=0;j<n_in;j++) {
        this->min[j]=(this->x)(0,j);
        this->max[j]=(this->x)(0,j);
        for(size_t i=1;i<n_points;i++) {
          double val=(this->x)(i,j);
          if (val>this->max[j]) this->max[j]=val;
          if (val<this->min[j]) this->min[j]=val;
        }
      }

      if (rescale==true) {
        for(size_t j=0;j<n_in;j++) {
          for(size_t i=0;i<n_points;i++) {
            (this->x)(i,j)=(((this->x)(i,j)-this->min[j])/
                            (this->max[j]-this->min[j])-0.5)*2.0;
          }
        }
        if (verbose>1) {
          std::cout << "interpm_krige_optim::set_data_noise() rescaled."
                    << std::endl;
        }
      }

      int success=0;

      this->Kinvf.resize(n_out);
      qual.resize(n_out);
      len.resize(n_out);
      ff1.resize(n_out);
      ff2.resize(n_out);

      // Loop over all output functions
      for(size_t iout=0;iout<n_out;iout++) {
      
        // Select the row of the data matrix
        mat2_row_t yiout(user_y,iout);
      
        if (iout<len_precompute.size()) {
	
          if (verbose>1) {
            std::cout << "interp_krige_optim: precomputed length "
                      << len_precompute[iout] << std::endl;
          }
          len[iout]=len_precompute[iout];
          size_t mode_temp=mode;
          mode=mode_final;
          qual[iout]=qual_fun(len[iout],noise_var[iout],iout,yiout,success);
          mode=mode_temp;
	
        } else if (full_min) {
	
          if (verbose>1) {
            std::cout << "interp_krige_optim: full minimization"
                      << std::endl;
          }
	
          double len_opt;
          // Choose average distance for first guess
          if (this->rescaled) {
            len_opt=(this->max[iout]-this->min[iout])/((double)this->np);
          } else {
            len_opt=2.0/((double)this->np);
          }
	
          funct mf=std::bind
            (std::mem_fn<double(double,double,size_t,mat2_row_t &,int &)>
             (&interpm_krige_optim<vec_t,mat_t,
              mat_row_t>::qual_fun<mat2_row_t>),
             this,std::placeholders::_1,noise_var[iout],iout,yiout,
             std::ref(success));
	
          mp->min(len_opt,qual[iout],mf);
          len[iout]=len_opt;
	
          if (success!=0) {
            if (err_on_fail) {
              O2SCL_ERR2("Minimization failed in ",
                         "interp_krige_optim::set_noise().",
                         o2scl::exc_efailed);
            }
          }
	
        } else {
	
          if (verbose>1) {
            std::cout << "interp_krige_optim::set_data_noise() : "
                      << "simple minimization" << std::endl;
          }
	
          if (len_guess_set==false) {
            double len_avg;
            // Choose average distance for first guess
            if (this->rescaled) {
              len_avg=(this->max[iout]-this->min[iout])/((double)this->np);
            } else {
              len_avg=2.0/((double)this->np);
            }

            len_min=len_avg/1.0e2;
            len_max=len_avg*1.0e2;
          }
          double len_ratio=len_max/len_min;

          if (verbose>1) {
            std::cout << "iout, len (min,max,step): " << iout
                      << " " << len_min << " " << len_max << " "
                      << pow(len_ratio,((double)1)/((double)nlen-1))
                      << std::endl;
          }
	
          // Initialize to zero to prevent uninit'ed var. warnings
          double min_qual=0.0, len_opt=0.0;
	
          if (verbose>1) {
            std::cout << "ilen len qual fail min_qual len_opt" << std::endl;
          }
	
          // Loop over the full range, finding the optimum
          bool min_set=false;
          for(size_t j=0;j<nlen;j++) {
            double xlen=len_min*pow(len_ratio,((double)j)/((double)nlen-1));

            success=0;
            qual[iout]=qual_fun(xlen,noise_var[iout],iout,yiout,success);
	
            if (success==0 && (min_set==false || qual[iout]<min_qual)) {
              len_opt=xlen;
              min_qual=qual[iout];
              min_set=true;
            }
	
            if (verbose>1) {
              std::cout << "interp_krige_optim: ";
              std::cout.width(2);
              std::cout << j << " " << xlen << " " << qual[iout] << " "
                        << success << " " << min_qual << " "
                        << len_opt << std::endl;
            }
	  
          }
      
          if (verbose>1) {
            std::cout << "interp_krige_optim: ";
            std::cout.width(2);
            std::cout << "   " << len_opt << " " << min_qual << std::endl;
          }
          qual[iout]=qual_fun(len_opt,noise_var[iout],iout,yiout,success);
	
          len[iout]=len_opt;
	
	
        }
      
        ff1[iout]=std::bind(std::mem_fn<double(const mat_row_t &,
                                               const mat_row_t &,
                                               size_t,double)>
                            (&interpm_krige_optim<vec_t,mat_t,
                             mat_row_t>::covar<mat_row_t,
                             mat_row_t>),this,
                            std::placeholders::_1,std::placeholders::_2,
                            n_in,len[iout]);
        ff2[iout]=std::bind(std::mem_fn<double(const mat_row_t &,
                                               const vec_t &,
                                               size_t,double)>
                            (&interpm_krige_optim<vec_t,mat_t,
                             mat_row_t>::covar<mat_row_t,
                             vec_t>),this,
                            std::placeholders::_1,std::placeholders::_2,
                            n_in,len[iout]);
      }
    
      return 0;
    }

    /** \brief Initialize the data for the interpolation
      
        \note This function works differently than 
        \ref o2scl::interpm_idw::set_data() . See this
        class description for more details.
    */
    template<class mat2_row_t, class mat2_t, class vec2_t>
    int set_data(size_t n_in, size_t n_out, size_t n_points,
                 mat_t &user_x, mat2_t &user_y,
                 const vec2_t &len_precompute,
                 bool rescale=false, bool err_on_fail=true) {
      vec_t noise_vec;
      noise_vec.resize(1);
      noise_vec[0]=0.0;
      return set_data_noise<mat2_row_t,mat2_t,vec_t,vec2_t>
        (n_in,n_out,n_points,user_x,
         user_y,noise_vec,len_precompute,rescale,err_on_fail);
    }
  
  };

  /** \brief Multi-dimensional interpolation by kriging with 
      nearest-neighbor 

      \note This class assumes that the function specified in the
      call to set_data() is the same as that passed to the
      eval() functions. If this is not the case, the
      behavior of this class is undefined.

      \note Experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
           class mat_t=boost::numeric::ublas::vector<double> >
  class interpm_krige_nn : public interpm_krige<vec_t,mat_t> {    
    
  public:
  
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
  
    interpm_krige_nn() {
      data_set=false;
      verbose=0;
    }

    /** \brief Verbosity parameter (default 0)
     */
    int verbose;

    /** \brief Initialize the data for the interpolation
     */
    template<class mat2_t, class vec_func_t>
    void set_data(size_t n_in, size_t n_out, size_t n_points,
                  mat_t &user_x, mat_t &user_y,
                  vec_func_t &fcovar, size_t order) {

      if (n_points<2) {
        O2SCL_ERR2("Must provide at least two points in ",
                   "interpm_krige_nn::set_data()",exc_efailed);
      }
      if (n_in<1) {
        O2SCL_ERR2("Must provide at least one input column in ",
                   "interpm_krige_nn::set_data()",exc_efailed);
      }
      if (n_out<1) {
        O2SCL_ERR2("Must provide at least one output column in ",
                   "interpm_krige_nn::set_data()",exc_efailed);
      }
      np=n_points;
      this->nd_in=n_in;
      nd_out=n_out;
      x.resize(n_points);
      n_order=order;
      for(size_t i=0;i<n_points;i++) {
        if (x[i].size()!=n_in) {
          O2SCL_ERR2("Size of x not correct in ",
                     "interpm_krige_nn::set_data().",o2scl::exc_efailed);
        }
        std::swap(x[i],x[i]);
      }
      y.resize(n_out);
      for(size_t i=0;i<n_out;i++) {
        if (y[i].size()!=n_points) {
          O2SCL_ERR2("Size of y not correct in ",
                     "interpm_krige_nn::set_data().",o2scl::exc_efailed);
        }
        std::swap(y[i],y[i]);
      }
      data_set=true;
      
      if (verbose>0) {
        std::cout << "interpm_krige_nn::set_data() : Using " << n_points
                  << " points with " << nd_in << " input variables and\n\t"
                  << nd_out << " output variables and order "
                  << n_order << " ." << std::endl;
      }

      return;
    }

    /** \brief Given covariance function \c fcovar and input vector \c x
        store the result of the interpolation in \c y
    */
    template<class vec2_t, class vec3_t, class vec_func_t>
    void eval(const vec2_t &x0, vec3_t &y0,
              vec_func_t &fcovar) const {
      
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_krige_nn::eval().",
                  exc_einval);
      }
      
      y0.resize(nd_out);
    
      // Loop over all output functions
      for(size_t iout=0;iout<nd_out;iout++) {

        // Find points closest to requested point, as defined
        // by the negative covariance for this output function
        ubvector dists(np);
        for(size_t ip=0;ip<np;ip++) {
          dists[ip]=-fcovar[iout](x0,x[ip]);
        }
      
        // Empty index vector (resized by the vector_smallest_index
        // function)
        ubvector_size_t index;
        o2scl::vector_smallest_index<ubvector,double,ubvector_size_t>
          (np,dists,n_order,index);
      
        // Construct subset of function values for nearest neighbors
        ubvector func(n_order);
        for(size_t io=0;io<n_order;io++) {
          func[io]=y[iout][index[io]];
        }
      
        // Construct the nearest neighbor KXX matrix
        mat3_t KXX(n_order,n_order);
        for(size_t irow=0;irow<n_order;irow++) {
          for(size_t icol=0;icol<n_order;icol++) {
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else {
              KXX(irow,icol)=fcovar[iout](x[index[irow]],
                                          x[index[icol]]);
            }
          }
        }
      
        // Construct the inverse of KXX
        o2scl_linalg::cholesky_decomp(n_order,KXX);
        mat3_t &inv_KXX=KXX;
        o2scl_linalg::cholesky_invert<mat3_t>(n_order,inv_KXX);
      
        // Inverse covariance matrix times function vector
        ubvector Kinvf(n_order);
        boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);

        // Comput the final result
        y0[iout]=0.0;
        for(size_t ipoints=0;ipoints<n_order;ipoints++) {
          y0[iout]+=-dists[index[ipoints]]*Kinvf[ipoints];
        }
      
      }
    
      return;
    
    }
  
    /** \brief Find a set of linearly independent points 

        Given a point \c x, a covariance function 
        \c fcovar, the index of the output function
        \c iout, and an array specifying the closest points 
        \c index, this function produces a list of 
    */
    template<class vec2_t, class vec_func_t>
    void find_lin_indep(const vec2_t &x2, size_t iout,
                        vec_func_t &fcovar,
                        ubvector_size_t &index,
                        ubvector_size_t &indep) const {
    
      if (x2.size()<nd_in || index.size()<np || indep.size()<n_order
          || iout>=nd_out) {
        std::cout << x2.size() << " " << nd_in << std::endl;
        std::cout << index.size() << " " << np << std::endl;
        std::cout << indep.size() << " " << n_order << std::endl;
        std::cout << iout << " " << nd_out << std::endl;
        O2SCL_ERR("Vectors not of correct size in find_lin_indep().",
                  o2scl::exc_einval);
      }
    
      bool done=false;
      while (done==false) {
      
        // Construct subset of function values for nearest neighbors
        ubvector func(n_order);
        for(size_t io=0;io<n_order;io++) {
          func[io]=y[iout][index[indep[io]]];
        }
      
        // Construct the nearest neighbor KXX matrix
        mat3_t KXX(n_order,n_order);
        for(size_t irow=0;irow<n_order;irow++) {
          for(size_t icol=0;icol<n_order;icol++) {
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else {
              KXX(irow,icol)=fcovar[iout](x2[index[indep[irow]]],
                                          x2[index[indep[icol]]]);
            }
          }
        }
      
        // Construct the inverse of KXX
        int cret=o2scl_linalg::cholesky_decomp(n_order,KXX);
        if (cret==0) {
          done=true;
        } else {
          if (verbose>1) {
            std::cout << "Finding new independent rows." << std::endl;
            for(size_t j=0;j<n_order;j++) {
              std::cout << indep[j] << " " << KXX(j,j) << std::endl;
            }
          }
          size_t max=o2scl::vector_max_value<ubvector_size_t,
                                             double>(indep);
          if (verbose>1) {
            std::cout << "Max is: " << max << std::endl;
          }
          for(size_t j=0;j<n_order;j++) {
            if (KXX(j,j)==0.0) {
              if (max+1<np) {
                if (verbose>1) {
                  std::cout << "Entry " << j << " is zero so replacing "
                            << "entry with " << max+1 << std::endl;
                }
                indep[j]=max+1;
                max++;
              } else {
                O2SCL_ERR3("Failed to find set of independent points ",
                           "in interpm_krige_nn::find_lin_indep",
                           "(const vec2_t &, size_t).",
                           o2scl::exc_efailed);
              }
            }
          }
        }
      
      }
      return;
    }
  
    /** \brief Given covariance function \c fcovar and input vector \c x
        return the result of the interpolation for function with 
        index \c iout
    */
    template<class vec2_t, class vec_func_t>
    double eval(const vec2_t &x2, size_t iout,
                vec_func_t &fcovar) const {
      
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_krige_nn::eval().",
                  exc_einval);
      }

      double ret;
    
      // Find points closest to requested point, as defined
      // by the negative covariance for this output function
      ubvector dists(np);
      for(size_t ip=0;ip<np;ip++) {
        dists[ip]=-fcovar[iout](x2,x2[ip]);
      }
      
      ubvector_size_t index(np);
      o2scl::vector_sort_index<ubvector,ubvector_size_t>(np,dists,index);

      // Vector for storing the indexes in the index array which
      // will store the closest n_order points which are
      // linearly independent
      ubvector_size_t indep(n_order);
      for(size_t io=0;io<n_order;io++) {
        indep[io]=io;
      }

      find_lin_indep(x2,iout,fcovar,index,indep);
    
      // Construct subset of function values for nearest neighbors
      ubvector func(n_order);
      for(size_t io=0;io<n_order;io++) {
        func[io]=y[iout][index[indep[io]]];
      }
      
      // Construct the nearest neighbor KXX matrix
      mat3_t KXX(n_order,n_order);
      for(size_t irow=0;irow<n_order;irow++) {
        for(size_t icol=0;icol<n_order;icol++) {
          if (irow>icol) {
            KXX(irow,icol)=KXX(icol,irow);
          } else {
            KXX(irow,icol)=fcovar[iout](x2[index[indep[irow]]],
                                        x2[index[indep[icol]]]);
          }
        }
      }
	
      // Construct the inverse of KXX
      o2scl_linalg::cholesky_decomp(n_order,KXX);
      mat3_t &inv_KXX=KXX;
      o2scl_linalg::cholesky_invert<mat3_t>(n_order,inv_KXX);
      
      // Inverse covariance matrix times function vector
      ubvector Kinvf(n_order);
      boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);

      // Comput the final result
      ret=0.0;
      for(size_t ipoints=0;ipoints<n_order;ipoints++) {
        ret+=-dists[index[indep[ipoints]]]*Kinvf[ipoints];
      }
      
      return ret;
      
    }

    /** \brief Compute a quality factor for interpolation
        using jackknife resampling
    */
    template<class vec2_t, class func_vec_t>
    double eval_jackknife(const vec2_t &x2, size_t iout,
                          func_vec_t &fcovar) const {
    
      if (data_set==false) {
        O2SCL_ERR("Data not set in interpm_krige_nn::eval_jackknife().",
                  exc_einval);
      }

      double qual=0.0;
    
      // Interpolated function value inside jackknife loop
      double ytmp;
    
      // For a distance measurement, just use the the negative
      // covariance for this output function
      ubvector dists(np);
      for(size_t ip=0;ip<np;ip++) {
        dists[ip]=-fcovar[iout](x2,x2[ip]);
      }
    
      // Create an index array which sorts by distance
      ubvector_size_t index(np);
      o2scl::vector_sort_index<ubvector,ubvector_size_t>(np,dists,index);
    
      // Vector for storing the indexes in the index array which
      // will store the closest n_order points which are
      // linearly independent
      ubvector_size_t indep(n_order);
      for(size_t io=0;io<n_order;io++) {
        indep[io]=io;
      }
    
      // -------------------------------------------------------------
      // Before the jackknife loop, we want to create a full
      // set of n_order linearly independent points

      find_lin_indep(x2,iout,fcovar,index,indep);

      // -------------------------------------------------------------
      // Now, the jackknife loop, removing one point at a time
    
      for(size_t jk=0;jk<n_order;jk++) {

        if (verbose>0) {
          std::cout << "Jackknife: " << jk << " matching function value "
                    << y[iout][index[jk]] << std::endl;
        }
      
        ubvector_size_t indep_jk;
        vector_copy_jackknife(indep,jk,indep_jk);
      
        // Construct subset of function values for nearest neighbors
        ubvector func(n_order-1);
        for(size_t io=0;io<n_order-1;io++) {
          func[io]=y[iout][index[indep_jk[io]]];
        }
      
        // Construct the nearest neighbor KXX matrix
        mat3_t KXX(n_order-1,n_order-1);
        for(size_t irow=0;irow<n_order-1;irow++) {
          for(size_t icol=0;icol<n_order-1;icol++) {
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else {
              KXX(irow,icol)=fcovar[iout](x2[index[indep_jk[irow]]],
                                          x2[index[indep_jk[icol]]]);
            }
          }
        }
	  
        // Construct the inverse of KXX
        o2scl_linalg::cholesky_decomp(n_order-1,KXX);
        mat3_t &inv_KXX=KXX;
        o2scl_linalg::cholesky_invert<mat3_t>(n_order-1,inv_KXX);
      
        // Inverse covariance matrix times function vector
        ubvector Kinvf(n_order-1);
        boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);
      
        // Comput the final result
        ytmp=0.0;
        for(size_t ipoints=0;ipoints<n_order-1;ipoints++) {
          ytmp+=-dists[index[indep_jk[ipoints]]]*Kinvf[ipoints];
        }
      
        // Add the squared deviation to y[iout]
        qual+=pow(y[iout][index[jk]]-ytmp,2.0);

        if (verbose>0) {
          std::cout << "Original value: "
                    << y[iout][index[jk]] << " interpolated: "
                    << ytmp << std::endl;
        }
      
        // End of jackknife loop
      }

      return qual;
    }
  
  
#ifndef DOXYGEN_INTERNAL
  
  protected:
    
    /// The order of the interpolation (specified by \ref set_data() )
    size_t n_order;
    /// The number of points
    size_t np;
    /// The number of dimensions of the inputs
    size_t nd_in;
    /// The number of dimensions of the outputs
    size_t nd_out;
    /// A vector of pointers holding the data
    std::vector<vec_t> x;
    /// A vector of pointers holding the data
    std::vector<vec_t> y;
    /// True if the data has been specified
    bool data_set;
    
#endif
    
  };

#endif
  
#ifndef DOXYGEN_NO_O2NS
}
#endif
    
#endif



