/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017, Andrew W. Steiner
  
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
    \brief One-dimensional interpolation classes and interpolation types
*/

#include <iostream>
#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <o2scl/interp.h>
#include <o2scl/lu.h>
#include <o2scl/columnify.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Interpolation by Kriging
      
      See also the \ref intp_section section of the \o2 User's guide. 

      \note The function \ref set() stores a pointer to the covariance
      function so it cannot go out of scope before any of the
      interpolation functions are called.

      \note This class is experimental.
  */
  template<class vec_t, class vec2_t=vec_t,
    class covar_func_t=std::function<double(double,double)> >
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
    
  public:
    
    interp_krige() {
      this->min_size=2;
      optimize=true;
    }
    
    virtual ~interp_krige() {}
    
    /// Initialize interpolation routine
    virtual void set(size_t size, const vec_t &x, const vec2_t &y) {

      double var_min=o2scl::vector_variance(size,y);
      std::vector<double> diff(size-1);
      for(size_t i=0;i<size-1;i++) {
	diff[i]=fabs(x[i+1]-x[i]);
      }
      double var_ratio=1.0e2;
      double len_min=o2scl::vector_min_value(size-1,diff)/3.0;
      double len_max=fabs(x[size-1]-x[0])*3.0;
      double len_ratio=len_max/len_min;

      double min_diff, var_opt, len_opt;
      size_t nvar=10, nlen=10;
      for(size_t i=0;i<nvar;i++) {
	double var=var_min*pow(var_ratio,((double)i)/((double)nvar-1));
	for(size_t j=0;j<nlen;j++) {
	  double len=len_min*pow(len_ratio,((double)i)/((double)nlen-1));
	}
      }
      
      return;
    }

    /// Initialize interpolation routine
    virtual void set(size_t n_dim, const vec_t &x, const vec_t &y,
		     covar_func_t &fcovar) {
      set_noise(n_dim,x,y,fcovar,0.0);
      return;
    }

    /// Initialize interpolation routine
    virtual void set_noise(size_t n_dim, const vec_t &x, const vec_t &y,
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
      o2scl::permutation p;
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
    
    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(double x0) const {

      double ret;
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

    interp_krige<vec_t,vec2_t,covar_func_t>
      (const interp_krige<vec_t,vec2_t,covar_func_t> &);
    interp_krige<vec_t,vec2_t,covar_func_t>& operator=
      (const interp_krige<vec_t,vec2_t,covar_func_t>&);

#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
