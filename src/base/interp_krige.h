/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2018, Andrew W. Steiner
  
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
#ifndef O2SCL_INTERP_KRIGE_H
#define O2SCL_INTERP_KRIGE_H

/** \file interp_krige.h
    \brief One-dimensional interpolation by Kriging
*/

#include <iostream>
#include <string>

#include <gsl/gsl_sf_erf.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <o2scl/interp.h>
#include <o2scl/lu.h>
#include <o2scl/columnify.h>
#include <o2scl/vector.h>
#include <o2scl/vec_stats.h>
#include <o2scl/mmin_bfgs2.h>
#include <o2scl/constants.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Interpolation by Kriging with a user-specified 
      covariance function
      
      See also the \ref intp_section section of the \o2 User's guide. 

      \note The function \ref set() stores a pointer to the covariance
      function and its derivatives and integrals so they cannot go out
      of scope before any of the interpolation functions are called.

      \note This class is experimental.
  */
  template<class vec_t, class vec2_t=vec_t,
    class covar_func_t=std::function<double(double,double)>,
    class covar_integ_t=std::function<double(double,double,double)> >
    class interp_krige : public interp_base<vec_t,vec2_t> {
    
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;
    
  protected:

    /** \brief Inverse covariance matrix times function vector
     */
    ubvector Kinvf;

    /** \brief Pointer to user-specified covariance function
     */
    covar_func_t *f;
    
    /** \brief Pointer to user-specified derivative
     */
    covar_func_t *df;
    
    /** \brief Pointer to user-specified second derivative
     */
    covar_func_t *df2;
    
    /** \brief Pointer to user-specified second derivative
     */
    covar_integ_t *intp;
    
  public:
    
    interp_krige() {
      this->min_size=2;
    }
    
    virtual ~interp_krige() {}

    /// Initialize interpolation routine
    virtual void set(size_t size, const vec_t &x, const vec2_t &y) {
      O2SCL_ERR2("Function set(size_t,vec_t,vec_t) unimplemented ",
		 "in interp_krige.",o2scl::exc_eunimpl);
      return;
    }

    /** \brief Initialize interpolation routine, specifying derivatives
	and integrals [not yet implemented]
    */
    virtual void set_covar_di_noise(size_t n_dim, const vec_t &x,
				    const vec_t &y, covar_func_t &fcovar,
				    covar_func_t &fderiv,
				    covar_func_t &fderiv2,
				    covar_func_t &finteg,
				    double noise_var) {
      O2SCL_ERR("Function set_covar_di_noise not yet implemented.",
		o2scl::exc_eunimpl);
      return;
    }
    
    /// Initialize interpolation routine
    virtual void set_covar_noise(size_t n_dim, const vec_t &x, const vec_t &y,
				 covar_func_t &fcovar, double noise_var) {
      
      if (n_dim<this->min_size) {
	O2SCL_ERR((((std::string)"Vector size, ")+szttos(n_dim)+", is less"+
		   " than "+szttos(this->min_size)+" in interp_krige::"+
		   "set().").c_str(),exc_einval);
      }

      // Store pointer to covariance function
      f=&fcovar;
      
      // Construct the KXX matrix
      ubmatrix KXX(n_dim,n_dim);
      for(size_t irow=0;irow<n_dim;irow++) {
	for(size_t icol=0;icol<n_dim;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else if (irow==icol) {
	    KXX(irow,icol)=fcovar(x[irow],x[icol])+noise_var;
	  } else {
	    KXX(irow,icol)=fcovar(x[irow],x[icol]);
	  }
	}
      }
      
      // Construct the inverse of KXX
      ubmatrix inv_KXX(n_dim,n_dim);
      o2scl::permutation p(n_dim);
      int signum;
      o2scl_linalg::LU_decomp(n_dim,KXX,p,signum);
      if (o2scl_linalg::diagonal_has_zero(n_dim,KXX)) {
	O2SCL_ERR("KXX matrix is singular in interp_krige::set().",
		  o2scl::exc_efailed);
      }
      o2scl_linalg::LU_invert<ubmatrix,ubmatrix,ubmatrix_column>
	(n_dim,KXX,p,inv_KXX);
      
      // Inverse covariance matrix times function vector
      Kinvf.resize(n_dim);
      boost::numeric::ublas::axpy_prod(inv_KXX,y,Kinvf,true);

      // Set parent data members
      this->px=&x;
      this->py=&y;
      this->sz=n_dim;
      return;
    }
    
    /// Initialize interpolation routine
    virtual void set_covar(size_t n_dim, const vec_t &x, const vec_t &y,
			   covar_func_t &fcovar) {
      set_covar_noise(n_dim,x,y,fcovar,0.0);
      return;
    }

    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(double x0) const {

      // Initialize to zero to prevent uninit'ed var. warnings
      double ret=0.0;
      for(size_t i=0;i<this->sz;i++) {
	ret+=(*f)(x0,(*this->px)[i])*Kinvf[i];
      }

      return ret;
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

    /// Return the type, \c "interp_linear".
    virtual const char *type() const { return "interp_krige"; }

#ifndef DOXYGEN_INTERNAL

  private:

    interp_krige<vec_t,vec2_t,covar_func_t,covar_integ_t>
      (const interp_krige<vec_t,vec2_t,covar_func_t,covar_integ_t> &);
    interp_krige<vec_t,vec2_t,covar_func_t,covar_integ_t>& operator=
      (const interp_krige<vec_t,vec2_t,covar_func_t,covar_integ_t>&);

#endif

  };


  /** \brief One-dimensional interpolation using an 
      optimized covariance function

      See also the \ref intp_section section of the \o2 User's guide. 

      \note This class is experimental.
  */
  template<class vec_t, class vec2_t=vec_t>
    class interp_krige_optim : public interp_krige<vec_t,vec2_t> {

  public:

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

  protected:

  /// The covariance function length scale
  double len;

  /// The covariance function coefficient
  double var;

  /// The quality factor of the optimization
  double qual;
  
  /// The covariance function
  double covar(double x1, double x2) {
    return var*exp(-(x1-x2)*(x1-x2)/len/len);
  }

  /// The derivative of the covariance function
  double deriv(double x1, double x2) {
    return -2.0*var*exp(-(x1-x2)*(x1-x2)/len/len)/len/len*(x1-x2);
  }

  /// The second derivative of the covariance function
  double deriv2(double x1, double x2) {
    return (4.0*(x1-x2)*(x1-x2)-2.0*len*len)*
    var*exp(-(x1-x2)*(x1-x2)/len/len)/len/len/len/len;
  }

  /// The integral of the covariance function
  double integ(double x, double x1, double x2) {
    return 0.5*len*sqrt(o2scl_const::pi)*var*
    (gsl_sf_erf((x2-x)/len)+gsl_sf_erf((x-x1)/len));
  }

  /// Functor for the covariance function \ref covar()
  std::function<double(double,double)> ff;

  /// Pointer to the user-specified minimizer
  mmin_base<> *mp;
  
  /** \brief Function to optimize the covariance
   */
  double qual_fun(size_t nv, const ubvector &x) {

    var=x[0];
    len=x[1];

    size_t size=this->sz;
    
    qual=0.0;
    for(size_t k=0;k<size;k++) {
      
      ubvector x2(size-1);
      o2scl::vector_copy_jackknife((*this->px),k,x2);
      ubvector y2(size-1);
      o2scl::vector_copy_jackknife((*this->py),k,y2);
      
      // Construct the KXX matrix
      ubmatrix KXX(size-1,size-1);
      for(size_t irow=0;irow<size-1;irow++) {
	for(size_t icol=0;icol<size-1;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else {
	    KXX(irow,icol)=exp(-var*pow(((*this->px)[irow]-
					 (*this->px)[icol])/len,2.0));
	  }
	}
      }
      
      // Construct the inverse of KXX
      ubmatrix inv_KXX(size-1,size-1);
      o2scl::permutation p(size-1);
      int signum;
      o2scl_linalg::LU_decomp(size-1,KXX,p,signum);
      if (o2scl_linalg::diagonal_has_zero(size-1,KXX)) {
	O2SCL_ERR("KXX matrix is singular in interp_krige::set().",
		  o2scl::exc_efailed);
      }
      o2scl_linalg::LU_invert<ubmatrix,ubmatrix,ubmatrix_column>
	(size-1,KXX,p,inv_KXX);
      
      // Inverse covariance matrix times function vector
      this->Kinvf.resize(size-1);
      boost::numeric::ublas::axpy_prod(inv_KXX,y2,this->Kinvf,true);
      
      double ypred=0.0;
      double yact=(*this->py)[k];
      for(size_t i=0;i<size-1;i++) {
	ypred+=exp(-var*pow(((*this->px)[k]-x2[i])/len,2.0))*this->Kinvf[i];
      }
      
      qual+=pow(yact-ypred,2.0);
      
    }

    return qual;
  }
  
  public:

  /// Number of variance points to try
  size_t nvar;

  /// Number of length scale points to try
  size_t nlen;

  /// Default minimizer
  mmin_bfgs2<> def_mmin;

  /// If true, use the full minimizer
  bool full_min;
  
  interp_krige_optim() {
    ff=std::bind(std::mem_fn<double(double,double)>
		 (&interp_krige_optim<vec_t,vec2_t>::covar),this,
		 std::placeholders::_1,std::placeholders::_2);
    nvar=20;
    nlen=20;
    full_min=false;
    mp=&def_mmin;
  }

  /// Initialize interpolation routine
  virtual void set_noise(size_t size, const vec_t &x, const vec2_t &y,
			 double noise_var) {

    if (full_min) {

      // Set parent data members
      this->px=&x;
      this->py=&y;
      this->sz=size;

      ubvector p(2);
      p[0]=o2scl::vector_variance(size,y);
      p[1]=x[1]-x[0];
      
      multi_funct mf=std::bind
      (std::mem_fn<double(size_t, const ubvector &)>
       (&interp_krige_optim<vec_t,vec2_t>::qual_fun),this,
       std::placeholders::_1,std::placeholders::_2);
      
      mp->mmin(2,p,qual,mf);

    } else {
      
      // Range of the coefficient parameter
      double var_min=o2scl::vector_variance(size,y);
      std::vector<double> diff(size-1);
      for(size_t i=0;i<size-1;i++) {
	diff[i]=fabs(x[i+1]-x[i]);
      }
      double var_ratio=1.0e2;
      
      // Range of the length parameter
      double len_min=o2scl::vector_min_value
      <std::vector<double>,double>(size-1,diff)/3.0;
      double len_max=fabs(x[size-1]-x[0])*3.0;
      double len_ratio=len_max/len_min;
      
      // Initialize to zero to prevent uninit'ed var. warnings
      double min_qual=0.0, var_opt=0.0, len_opt=0.0;
      
      // Loop over the full range, finding the optimum 
      for(size_t i=0;i<nvar;i++) {
	var=var_min*pow(var_ratio,((double)i)/((double)nvar-1));
	for(size_t j=0;j<nlen;j++) {
	  len=len_min*pow(len_ratio,((double)j)/((double)nlen-1));
	  
	  qual=0.0;
	  for(size_t k=0;k<size;k++) {
	    
	    ubvector x2(size-1);
	    o2scl::vector_copy_jackknife(x,k,x2);
	    ubvector y2(size-1);
	    o2scl::vector_copy_jackknife(y,k,y2);
	    
	    // Construct the KXX matrix
	    ubmatrix KXX(size-1,size-1);
	    for(size_t irow=0;irow<size-1;irow++) {
	      for(size_t icol=0;icol<size-1;icol++) {
		if (irow>icol) {
		  KXX(irow,icol)=KXX(icol,irow);
		} else {
		  KXX(irow,icol)=exp(-var*pow((x[irow]-x[icol])/len,2.0));
		}
	      }
	    }
	    
	    // Construct the inverse of KXX
	    ubmatrix inv_KXX(size-1,size-1);
	    o2scl::permutation p(size-1);
	    int signum;
	    o2scl_linalg::LU_decomp(size-1,KXX,p,signum);
	    if (o2scl_linalg::diagonal_has_zero(size-1,KXX)) {
	      O2SCL_ERR("KXX matrix is singular in interp_krige::set().",
			o2scl::exc_efailed);
	    }
	    o2scl_linalg::LU_invert<ubmatrix,ubmatrix,ubmatrix_column>
	      (size-1,KXX,p,inv_KXX);
	    
	    // Inverse covariance matrix times function vector
	    this->Kinvf.resize(size-1);
	    boost::numeric::ublas::axpy_prod(inv_KXX,y2,this->Kinvf,true);
	    
	    double ypred=0.0;
	    double yact=y[k];
	    for(size_t i=0;i<size-1;i++) {
	      ypred+=exp(-var*pow((x[k]-x2[i])/len,2.0))*this->Kinvf[i];
	    }
	    
	    qual+=pow(yact-ypred,2.0);
	    
	  }
	  
	  if ((i==0 && j==0) || qual<min_qual) {
	    var_opt=var;
	    len_opt=len;
	    min_qual=qual;
	  }
	}
      }
      
      // Now that we've optimized the covariance function,
      // just use the parent class to interpolate
      var=var_opt;
      len=len_opt;

    }

    this->set_covar_noise(size,x,y,ff,noise_var);
      
    return;
  }

  /// Initialize interpolation routine
  virtual void set(size_t size, const vec_t &x, const vec2_t &y) {
    set_noise(size,x,y,0.0);
    return;
  }
  
  
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
