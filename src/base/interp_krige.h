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
#include <o2scl/vector.h>
#include <o2scl/vec_stats.h>
#include <o2scl/min_brent_gsl.h>
#include <o2scl/constants.h>
#include <o2scl/invert.h>
#include <o2scl/prob_dens_func.h>

namespace o2scl {

  /** \brief Desc
   */
  class covar_funct {
    
  public:
    
    /// The function
    virtual double operator()(double x, double y)=0;

    /// The function
    virtual double deriv(double x, double y)=0;
    
    /// The function
    virtual double deriv2(double x, double y)=0;
    
    /// The function
    virtual double integ(double x, double a, double b)=0;
    
  };
  
  /** \brief Desc
   */
  class covar_funct_rbf_noise : public covar_funct {
    
  public:

    /// Desc
    double len;
    
    /// Desc
    double log10_noise;
    
    /// The covariance function
    virtual double operator()(double x1, double x2) {
      double ret=exp(-(x1-x2)*(x1-x2)/len/len/2.0);
      if (x1==x2) ret+=pow(10.0,log10_noise);
      return ret;
    }

    /** \brief The derivative of the covariance function with
        respect to the first argument
    */
    virtual double deriv(double x1, double x2) {
      return -exp(-(x1-x2)*(x1-x2)/len/len/2.0)/len/len*(x1-x2);
    }
    
    /** \brief The second derivative of the covariance function with
        respect to the first argument
    */
    virtual double deriv2(double x1, double x2) {
      return ((x1-x2)*(x1-x2)-len*len)*
        exp(-(x1-x2)*(x1-x2)/len/len/2.0)/len/len/len/len;
    }

    /** \brief The integral of the covariance function over 
        \f$ [a,b] \f$
        
        The integral of the function is
        \f[
        \int_a^b f(x) dx = \sum_i A_i \int_a^b C(x,x_i) dx
        \f]
        where \f$ A_i = (K^{-1})_{ij} f_j \f$. To compute
        the integral we use
        \f[
        \int_a^b C(x,x_i) dx = 
        \int_{a+x_i}^{b+x_i} \exp \left( - \frac{x^2}{2 L^2} \right) dx = 
        \int_{(a+x_i)/(L\sqrt{2})}^{(b+x_i)/(L\sqrt{2})} 
        L \sqrt{2} \exp \left( - y^2 \right) dy 
        \f]
        But 
        \f[
        \mathrm{erf}(x) \equiv \frac{2}{\sqrt{\pi}} \int_0^{x} e^{-t^2}
        \f]
        so
        \f[
        \int_a^b C(x,x_i) dx = 
        L \frac{\sqrt{\pi}}{\sqrt{2}} \left[ 
        \mathrm{erf}\left( \frac{b+x_i}{L \sqrt{2}} \right) - 
        \mathrm{erf}\left( \frac{a+x_i}{L \sqrt{2}} \right) \right]
        \f]
    */
    virtual double integ(double x, double a, double b) {
      double alpha=1.0/(len*len*2.0);
      return sqrt(o2scl_const::pi/alpha)/2.0*
        (gsl_sf_erf(sqrt(alpha)*(b-x))-
         gsl_sf_erf(sqrt(alpha)*(a-x)));
    }

  };

  /** \brief Desc
   */
  class covar_funct_strings : public covar_funct {
    
  public:

    covar_funct_strings() {
    }
    
    /** \brief Desc
     */
    template<class vec_string_t=std::vector<std::string> >
    covar_funct_strings(std::string expr, std::string expr_d,
                        std::string expr_d2, std::string expr_i,
                        vec_string_t &parms, 
                        std::string var, std::string var2,
                        std::string var_lo, std::string var_hi) {
      set(expr,expr_d,expr_d2,expr_i,parms,var,var2,var_lo,var_hi);
    }
    
    /** \brief Desc
     */
    template<class vec_string_t=std::vector<std::string> >
    void set(std::string expr, std::string expr_d,
             std::string expr_d2, std::string expr_i,
             vec_string_t &parms, 
             std::string var, std::string var2,
             std::string var_lo, std::string var_hi) {
      
      calc.compile(expr.c_str(),0);
      st_expr=expr;
      calc_d.compile(expr_d.c_str(),0);
      st_expr_d=expr_d;
      calc_d2.compile(expr_d2.c_str(),0);
      st_expr_d2=expr_d2;
      calc_i.compile(expr_i.c_str(),0);
      st_expr_i=expr_i;
      
      int np=parms.size();
      st_parms.resize(np);
      for (int i=0;i<np;i++) {
	st_parms[i]=parms[i];
      }
      st_var=var;
      st_var2=var;
      st_var_lo=var_lo;
      st_var_hi=var_hi;
      
      return;
    }

    /// The function parser
    o2scl::calc_utf8<> calc;
    o2scl::calc_utf8<> calc_d;
    o2scl::calc_utf8<> calc_d2;
    o2scl::calc_utf8<> calc_i;

    /// The variable values
    std::map<std::string,double> vars;
      
    /// The expression
    std::string st_expr;
    std::string st_expr_d;
    std::string st_expr_d2;
    std::string st_expr_i;
      
    /// The parameters
    std::vector<std::string> st_parms; 

    /// The variable
    std::string st_var; 

    /// The variable
    std::string st_var2; 

    /// The variable
    std::string st_var_lo; 

    /// The variable
    std::string st_var_hi; 

    /// To evaluate the fit quality
    template<class vec_t>
    void set_params(size_t np, vec_t &v) {
      for(size_t i=0;i<np;i++) {
	vars[st_parms[i]]=v[i];
      }
      return;
    }
    
    /// The function
    double operator()(double x, double y) {
      vars[st_var]=x;
      vars[st_var2]=y;
      double z=calc.eval(&vars);
      return z;
    }

    /// The function
    double deriv(double x, double y) {
      vars[st_var]=x;
      vars[st_var2]=y;
      double z=calc_d.eval(&vars);
      return z;
    }
    
    /// The function
    double deriv2(double x, double y) {
      vars[st_var]=x;
      vars[st_var2]=y;
      double z=calc_d2.eval(&vars);
      return z;
    }
    
    /// The function
    double integ(double x, double a, double b) {
      vars[st_var]=x;
      vars[st_var_lo]=a;
      vars[st_var_hi]=b;
      double z=calc_i.eval(&vars);
      return z;
    }
    
  };
  
  /** \brief Interpolation by Kriging with a user-specified 
      covariance function
      
      \verbatim embed:rst
      See also the :ref:`Interpolation` section of the 
      O\ :sub:`2`\ scl User's guide. 
      \endverbatim

      \note The function \ref set() stores a pointer to the covariance
      function and its derivatives and integrals so they cannot go out
      of scope before any of the interpolation functions are called.

      \note This class is experimental.
  */
  template<class vec_t, class vec2_t=vec_t,
           class covar_func_t=std::function<double(double,double)>,
           class covar_integ_t=std::function<double(double,double,double)>,
           class mat_t=boost::numeric::ublas::matrix<double>,
           class mat_inv_t=o2scl_linalg::matrix_invert_det_cholesky<mat_t> >
  class interp_krige : public interp_base<vec_t,vec2_t> {
    
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    
  protected:
    
    /** \brief Inverse covariance matrix times function vector
     */
    ubvector Kinvf;

    /** \brief Pointer to user-specified covariance function
     */
    covar_func_t *f;
    
    /** \brief Pointer to user-specified derivative
     */
    covar_func_t *fd;
    
    /** \brief Pointer to user-specified second derivative
     */
    covar_func_t *fd2;
    
    /** \brief Pointer to user-specified integral
     */
    covar_integ_t *fi;

    /** \brief If true, then the data was rescaled to zero mean and unit
        standard deviation
    */
    bool rescaled;

    /// Mean before rescaling
    double mean_y;

    /// Standard deviation before rescaling
    double std_y;

    // Rescaled y vector
    ubvector y_r;

    /// The inverse of the covariance matrix
    mat_t inv_KXX;

    /// The matrix inversion object
    mat_inv_t mi;

    /** \brief Initialize interpolation routine, internal version
        with pointers for derivative and integral functions
    */
    virtual int set_covar_di_noise_internal
      (size_t n_dim, const vec_t &x, const vec_t &y,
       covar_func_t &fcovar, covar_func_t *fderiv, covar_func_t *fderiv2,
       covar_integ_t *finteg, double noise_var, bool rescale=false) {
      
      f=&fcovar;
      fd=fderiv;
      fd2=fderiv2;
      fi=finteg;
      
      if (n_dim<this->min_size) {
	O2SCL_ERR((((std::string)"Vector size, ")+szttos(n_dim)+
                   ", is less than "+szttos(this->min_size)+
                   " in interp_krige::"+"set().").c_str(),
                  exc_einval);
      }

      if (rescale==true) {
        rescaled=true;
        mean_y=o2scl::vector_mean(n_dim,y);
        std_y=o2scl::vector_stddev(n_dim,y,mean_y);
        y_r.resize(n_dim);
        for (size_t j=0;j<n_dim;j++) {
          y_r[j]=(y[j]-mean_y)/std_y;
        }
      } else {
        rescaled=false;
      }

      // We put the original KXX matrix in inv_KXX, knowing that
      // we will invert it below
      
      // Construct the KXX matrix
      inv_KXX.resize(n_dim,n_dim);
      for(size_t irow=0;irow<n_dim;irow++) {
        for(size_t icol=0;icol<n_dim;icol++) {
          if (irow>icol) {
            inv_KXX(irow,icol)=inv_KXX(icol,irow);
          } else if (irow==icol) {
            inv_KXX(irow,icol)=fcovar(x[irow],x[icol])+noise_var;
          } else {
            inv_KXX(irow,icol)=fcovar(x[irow],x[icol]);
          }
        }
      }

      // Now invert in-place to get inv_KXX
      mi.invert_inplace(n_dim,inv_KXX);

      // Inverse covariance matrix times function vector
      Kinvf.resize(n_dim);
      if (rescaled) {
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           n_dim,n_dim,1.0,inv_KXX,y_r,0.0,Kinvf);
      } else {
        o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                           o2scl_cblas::o2cblas_NoTrans,
                           n_dim,n_dim,1.0,inv_KXX,y,0.0,Kinvf);
      }
      
      if (!keep_matrix) {
        // If not needed, free the associated memory
        inv_KXX.resize(0,0);
      }
      
      // Set parent data members
      this->px=&x;
      this->py=&y;
      this->sz=n_dim;
      
      return 0;
   }
   
  public:
    
    interp_krige() {
      this->min_size=2;
      keep_matrix=true;
      f=0;
      fd=0;
      fd2=0;
      fi=0;
      err_nonconv=true;
    }
    
    virtual ~interp_krige() {}

    /** \brief If true, call the error handler for convergence 
        errors (default true)
     */
    bool err_nonconv;
    
    /** \brief If true, keep \f$ K^{-1} \f$ (default true)

        This must be set to true in order to use \ref sigma() or
        \ref gen_dist() .
    */
    bool keep_matrix;
    
    /// Initialize interpolation routine
    virtual void set(size_t size, const vec_t &x, const vec2_t &y) {
      O2SCL_ERR2("Function set(size_t,vec_t,vec_t) unimplemented ",
		 "in interp_krige.",o2scl::exc_eunimpl);
      return;
    }
    
    /** \brief Initialize interpolation routine with covariance
        function, derivatives and integrals, and noise term
    */
    virtual int set_covar_di_noise
      (size_t n_dim, const vec_t &x,
       const vec_t &y, covar_func_t &fcovar,
       covar_func_t &fderiv, covar_func_t &fderiv2,
       covar_integ_t &finteg, double noise_var, bool rescale=false) {
      return set_covar_di_noise_internal(n_dim,x,y,fcovar,
                                         &fderiv,&fderiv2,&finteg,noise_var,
                                         rescale);
    }
    
    /** \brief Initialize interpolation routine with covariance
        function and noise term, but without derivatives and integrals
    */
    virtual int set_covar_noise(size_t n_dim, const vec_t &x, const vec_t &y,
				covar_func_t &fcovar, double noise_var,
                                bool rescale=false) {
      return set_covar_di_noise_internal(n_dim,x,y,fcovar,0,0,0,noise_var,
                                         rescale);
    }
    
    /** \brief Initialize interpolation routine with covariance
        function, with automatic noise term and without derivatives and
        integrals
    */
    virtual int set_covar(size_t n_dim, const vec_t &x, const vec_t &y,
			  covar_func_t &fcovar, bool rescale=false) {
      
      // Use the mean absolute value to determine noise

      // AWS, 7/30/22, I don't know why, but using mean_abs/1.0e8
      // makes the testing code fail
      double mean_abs=0.0;
      for(size_t j=0;j<n_dim;j++) {
        mean_abs+=fabs(y[j]);
      }
      mean_abs/=n_dim;
      
      return set_covar_noise(n_dim,x,y,fcovar,0.0,rescale);
                             
    }

    /// Give the value of the function \f$ y(x=x_0) \f$ .
    virtual double eval(double x0) const {

      double ret=0.0;
      
      for(size_t i=0;i<this->sz;i++) {
        ret+=(*f)(x0,(*this->px)[i])*Kinvf[i];
      }

      if (rescaled) {
        ret=ret*std_y+mean_y;
      }

      return ret;
    }

    /** \brief Evaluate the interpolation at \f$ x=x_0 \f$ using
        an alternate covariance function
    */
    template<class covar_func2_t>
      double eval_covar(double x0, covar_func2_t &user_f) const {
      
      double ret=0.0;
      
      for(size_t i=0;i<this->sz;i++) {
        ret+=user_f(x0,(*this->px)[i])*Kinvf[i];
      }

      if (rescaled) {
        ret=ret*std_y+mean_y;
      }

      return ret;
    }

    /// Give the value of the derivative \f$ y^{\prime}(x=x_0) \f$ .
    virtual double deriv(double x0) const {

      if (fd==0) {
        O2SCL_ERR("Derivative not specified in interp_krige::deriv().",
                  o2scl::exc_einval);
      }
      
      double ret=0.0;
      
      for(size_t i=0;i<this->sz;i++) {
        ret+=(*fd)(x0,(*this->px)[i])*Kinvf[i];
      }
      
      if (rescaled) {
        ret=ret*std_y+mean_y;
      }

      return ret;
    }
    
    /** \brief Give the value of the second derivative  
	\f$ y^{\prime \prime}(x=x_0) \f$
    */
    virtual double deriv2(double x0) const {
      
      if (fd2==0) {
        O2SCL_ERR("Derivative not specified in interp_krige::deriv().",
                  o2scl::exc_einval);
      }
      
      double ret=0.0;
      
      for(size_t i=0;i<this->sz;i++) {
        ret+=(*fd2)(x0,(*this->px)[i])*Kinvf[i];
      }

      if (rescaled) {
        ret=ret*std_y+mean_y;
      }

      return ret;
    }
    
    /// Give the value of the integral \f$ \int_a^{b}y(x)~dx \f$ .
    virtual double integ(double a, double b) const {
      double ret=0.0;
      
      if (fi==0) {
        O2SCL_ERR("Integral not specified in interp_krige::deriv().",
                  o2scl::exc_einval);
      }
      
      for(size_t i=0;i<this->sz;i++) {
        ret+=(*fi)((*this->px)[i],a,b)*Kinvf[i];
      }
      
      if (rescaled) {
        ret=ret*std_y+mean_y;
      }

      return ret;
    }

    /** \brief Return the interpolation uncertainty from the 
        Gaussian process
    */
    double sigma(double x0) const {
      
      if (!keep_matrix) {
        O2SCL_ERR2("Matrix information missing (keep_matrix==false) in ",
                   "interp_krige::sigma().",o2scl::exc_einval);
      }
      
      double sigma;
      ubvector kxx0(this->sz), prod(this->sz);
      double kx0x0=(*f)(x0,x0);
      
      for(size_t i=0;i<this->sz;i++) {
        kxx0[i]=(*f)(x0,(*this->px)[i]);
      }

      o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                         o2scl_cblas::o2cblas_NoTrans,
                         this->sz,this->sz,1.0,inv_KXX,kxx0,0.0,prod);
      sigma=kx0x0-o2scl_cblas::ddot(this->sz,kxx0,prod);
      //boost::numeric::ublas::axpy_prod(inv_KXX,kxx0,prod,true);
      //sigma=kx0x0-boost::numeric::ublas::inner_prod(kxx0,prod);
      
      if (rescaled) {
        sigma*=std_y;
      }

      return sigma;
    }
    
    /** \brief Generate a probability distribution for the interpolation
        at a specified point

        This function implements Eq. 2.19 of R&W
    */
    prob_dens_gaussian gen_dist(double x0) const {

      if (!keep_matrix) {
        O2SCL_ERR2("Matrix information missing (keep_matrix==false) in ",
                   "interp_krige::gen_dist().",o2scl::exc_einval);
      }
      
      double cent=0.0;
      double sigma;
      ubvector kxx0(this->sz), prod(this->sz);
      double kx0x0=(*f)(x0,x0);
      
      for(size_t i=0;i<this->sz;i++) {
        kxx0[i]=(*f)(x0,(*this->px)[i]);
        cent+=kxx0[i]*Kinvf[i];
      }
      
      if (rescaled) {
        cent=cent*std_y+mean_y;
      }
        
      o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                         o2scl_cblas::o2cblas_NoTrans,
                         this->sz,this->sz,1.0,inv_KXX,kxx0,0.0,prod);
      sigma=kx0x0-o2scl_cblas::ddot(this->sz,kxx0,prod);
      //boost::numeric::ublas::axpy_prod(inv_KXX,kxx0,prod,true);
      //sigma=kx0x0-boost::numeric::ublas::inner_prod(kxx0,prod);

      if (rescaled) {
        sigma*=std_y;
      }

      if (sigma<0.0) sigma=0.0;
      return prob_dens_gaussian(cent,sigma);
    }
    
    /** \brief Sample the probability distribution for the interpolation
        at a specified point

        This creates a new distribution at the point \x c0 and then
        samples that distribution. If one wants to perform several
        samples at the same point, it is much more efficient to
        use gen_dist() instead.
    */
    double sample(double x0) const {
      return gen_dist(x0)();
    }
    
    /** \brief Generate a probability distribution for the interpolation
        at a vector of points
    */
    template<class vec3_t, class vec4_t>
      void sample_vec(vec3_t &x, vec4_t &y) const {
      y.resize(x.size());
      for(size_t j=0;j<x.size();j++) {
        y[j]=sample(x[j]);
      }
      return;
    }
    
    /// Return the type, \c "interp_krige".
    virtual const char *type() const { return "interp_krige"; }

  private:

    interp_krige<vec_t,vec2_t,covar_func_t,
                 covar_integ_t,mat_t,mat_inv_t>
      (const interp_krige<vec_t,vec2_t,covar_func_t,
       covar_integ_t,mat_t,mat_inv_t> &);
    interp_krige<vec_t,vec2_t,covar_func_t,
                 covar_integ_t,mat_t,mat_inv_t>& operator=
      (const interp_krige<vec_t,vec2_t,covar_func_t,
       covar_integ_t,mat_t,mat_inv_t>&);
    
  };


  /** \brief One-dimensional interpolation using an 
      simple optimized covariance function

      \verbatim embed:rst
      See also the :ref:`Interpolation` section of the 
      O\ :sub:`2`\ scl User's guide. 
      \endverbatim

      \note This class is experimental.

      \verbatim embed:rst

      .. todo:

         Future: Class interp_krige_optim only allows one parameter,
         but we could generalize it to allow several hyperparameters,
         similar to the way sklearn works.
         
      \endverbatim
  */
  template<class vec_t, class vec2_t=vec_t,
           class mat_t=boost::numeric::ublas::matrix<double>,
           class mat_inv_t=o2scl_linalg::matrix_invert_det_cholesky<mat_t> >
  class interp_krige_optim :
    public interp_krige<vec_t,vec2_t,
                        std::function<double(double,double)>,
                        std::function<double(double,double,double)>,
                        mat_t,mat_inv_t> {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

  protected:

    /// Function object for the covariance
    std::function<double(double,double)> ff;
  
    /// Function object for the covariance
    std::function<double(double,double)> ffd;
  
    /// Function object for the covariance
    std::function<double(double,double)> ffd2;
  
    /// Function object for the covariance
    std::function<double(double,double,double)> ffi;
  
    /// The covariance function length scale
    double len;
    
    /// The quality factor of the optimization
    double qual;
    
    /// The covariance function
    double covar(double x1, double x2) {
      return exp(-(x1-x2)*(x1-x2)/len/len/2.0);
    }

    /** \brief The derivative of the covariance function with
        respect to the first argument
    */
    double deriv_covar(double x1, double x2) {
      return -exp(-(x1-x2)*(x1-x2)/len/len/2.0)/len/len*(x1-x2);
    }
    
    /** \brief The second derivative of the covariance function with
        respect to the first argument
    */
    double deriv2_covar(double x1, double x2) {
      return ((x1-x2)*(x1-x2)-len*len)*
        exp(-(x1-x2)*(x1-x2)/len/len/2.0)/len/len/len/len;
    }

    /** \brief The integral of the covariance function
        
        The integral of the function is
        \f[
        \int_a^b f(x) dx = \sum_i A_i \int_a^b C(x,x_i) dx
        \f]
        where \f$ A_i = (K^{-1})_{ij} f_j \f$. To compute
        the integral we use
        \f[
        \int_a^b C(x,x_i) dx = 
        \int_{a+x_i}^{b+x_i} \exp \left( - \frac{x^2}{2 L^2} \right) dx = 
        \int_{(a+x_i)/(L\sqrt{2})}^{(b+x_i)/(L\sqrt{2})} 
        L \sqrt{2} \exp \left( - y^2 \right) dy 
        \f]
        But 
        \f[
        \mathrm{erf}(x) \equiv \frac{2}{\sqrt{\pi}} \int_0^{x} e^{-t^2}
        \f]
        so
        \f[
        \int_a^b C(x,x_i) dx = 
        L \frac{\sqrt{\pi}}{\sqrt{2}} \left[ 
        \mathrm{erf}\left( \frac{b+x_i}{L \sqrt{2}} \right) - 
        \mathrm{erf}\left( \frac{a+x_i}{L \sqrt{2}} \right) \right]
        \f]
    */
    double integ_covar(double x, double a, double b) {
      double alpha=1.0/(len*len*2.0);
      return sqrt(o2scl_const::pi/alpha)/2.0*
        (gsl_sf_erf(sqrt(alpha)*(b-x))-
         gsl_sf_erf(sqrt(alpha)*(a-x)));
    }

    /// Pointer to the user-specified minimizer
    min_base<> *mp;

  public:
    
    /** \brief Function to optimize the covariance parameters
     */
    double qual_fun(double x, double noise_var, int &success) {

      len=x;
      success=0;

      size_t size=this->sz;

      if (mode==mode_loo_cv_bf) {
      
        qual=0.0;
        for(size_t k=0;k<size;k++) {
	
          // Leave one observation out
          ubvector x2(size-1);
          ubvector y2(size-1);
          o2scl::vector_copy_jackknife((*this->px),k,x2);
          if (this->rescaled) {
            o2scl::vector_copy_jackknife(this->y_r,k,y2);
          } else {
            o2scl::vector_copy_jackknife((*this->py),k,y2);
          }
	
          // Construct the inverse of the KXX matrix. Note
          // that we don't use the inv_KXX data member
          // because of the size mismatch
          mat_t inv_KXX2(size-1,size-1);
          for(size_t irow=0;irow<size-1;irow++) {
            for(size_t icol=0;icol<size-1;icol++) {
              if (irow>icol) {
                inv_KXX2(irow,icol)=inv_KXX2(icol,irow);
              } else {
                inv_KXX2(irow,icol)=covar(x2[irow],x2[icol]);
                if (irow==icol) inv_KXX2(irow,icol)+=noise_var;
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
          double yact;
          if (this->rescaled) {
            yact=this->y_r[k];
          } else {
            yact=(*this->py)[k];
          }

          // Compute the predicted value
          double ypred=0.0;
          ubvector kxx0(size-1);
          if (this->rescaled) {
            for(size_t i=0;i<size-1;i++) {
              kxx0[i]=covar((*this->px)[k],x2[i]);
              ypred+=kxx0[i]*Kinvf2[i];
              if (false && k==0) {
                std::cout << "  " << Kinvf2[i] << " " << x2[i] << " "
                          << kxx0[i] << " " << ypred
                          << std::endl;
              }
            }
          } else {
            for(size_t i=0;i<size-1;i++) {
              kxx0[i]=covar((*this->px)[k],x2[i]);
              ypred+=kxx0[i]*Kinvf2[i];
            }
          }
          
          if (false && k==0) {
            std::cout.setf(std::ios::showpos);
            std::cout << "k,x,yact,ypred: " << k << " "
                      << (*this->px)[k] << " "
                      << yact << " " << ypred << std::endl;
            std::cout.unsetf(std::ios::showpos);
          }
          
          qual+=pow(yact-ypred,2.0);
	
        }
        if (verbose>0) {
          std::cout << "len,qual (loo_cv_bf): " << len << " "
                    << qual << std::endl;
          //std::cout << "--" << std::endl;
        }

      } else if (mode==mode_loo_cv) {

#ifdef O2SCL_NEVER_DEFINED
        
        if (false) {
          
          // AWS, 7/30/22: This is the loo_cv
          // method to compute mu and sigma and Eq 5.10 in R&W., I
          // think it should work, but it doesn't quite work, however,
          // using Eq. 5.12 to compute mu and sigma (below) seems to
          // work well (and is more efficient), so I'm keeping that
          
          qual=0.0;
          for(size_t k=0;k<size;k++) {
	
            // Leave one observation out
            ubvector x2(size-1);
            ubvector y2(size-1);
            o2scl::vector_copy_jackknife((*this->px),k,x2);
            if (this->rescaled) {
              o2scl::vector_copy_jackknife(this->y_r,k,y2);
            } else {
              o2scl::vector_copy_jackknife((*this->py),k,y2);
            }
	
            // Construct the KXX matrix
            mat_t inv_KXX2(size-1,size-1);
            for(size_t irow=0;irow<size-1;irow++) {
              for(size_t icol=0;icol<size-1;icol++) {
                if (irow>icol) {
                  inv_KXX2(irow,icol)=inv_KXX2(icol,irow);
                } else {
                  inv_KXX2(irow,icol)=covar(x2[irow],x2[icol]);
                  if (irow==icol) inv_KXX2(irow,icol)+=noise_var;
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
	  
            double yact;
            if (this->rescaled) {
              yact=this->y_r[k];
            } else {
              yact=(*this->py)[k];
            }
          
            double ypred=0.0;
            ubvector kxx0(size-1);
            if (this->rescaled) {
              for(size_t i=0;i<size-1;i++) {
                kxx0[i]=covar((*this->px)[k],x2[i]);
                ypred+=kxx0[i]*Kinvf2[i];
              }
            } else {
              for(size_t i=0;i<size-1;i++) {
                kxx0[i]=covar((*this->px)[k],x2[i]);
                ypred+=kxx0[i]*Kinvf2[i];
              }
            }
          
            double kx0x0=covar((*this->px)[k],(*this->px)[k]);
            ubvector prod(size-1);
          
            boost::numeric::ublas::axpy_prod(inv_KXX2,kxx0,prod,true);
            double sigma=kx0x0-
              boost::numeric::ublas::inner_prod(kxx0,prod);

            if (this->rescaled) {
              sigma*=this->std_y;
            }

            if (false && k==0) {
              std::cout.setf(std::ios::showpos);
              std::cout << "k,x,yact,ypred,sigma: " << k << " "
                        << (*this->px)[k] << " "
                        << yact << " " << ypred << " "
                        << sigma << std::endl;
              std::cout.unsetf(std::ios::showpos);
            }

            // We maximize the predictive log probability, Eq 5.10
            // in R&W
            qual+=pow(yact-ypred,2.0)/sigma/sigma/2.0;
            qual+=0.5*log(sigma*sigma);
	
          }

        }
        
#endif

        if (true) {
          
          // Construct the KXX matrix
          mat_t inv_KXX2(size,size);
          for(size_t irow=0;irow<size;irow++) {
            for(size_t icol=0;icol<size;icol++) {
              if (irow>icol) {
                inv_KXX2(irow,icol)=inv_KXX2(icol,irow);
              } else {
                inv_KXX2(irow,icol)=covar((*this->px)[irow],
                                          (*this->px)[icol]);
                if (irow==icol) inv_KXX2(irow,icol)+=noise_var;
              }
            }
          }
          
          // Construct the inverse of KXX
          this->mi.invert_inplace(size-1,inv_KXX2);
	  
          // Inverse covariance matrix times function vector
          ubvector Kinvf2(size);
          if (this->rescaled) {
            o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                               o2scl_cblas::o2cblas_NoTrans,
                               size,size,1.0,inv_KXX2,this->y_r,0.0,
                               Kinvf2);
          } else {
            o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                               o2scl_cblas::o2cblas_NoTrans,
                               size,size,1.0,inv_KXX2,*(this->py),0.0,
                               Kinvf2);
          }

          qual=0.0;
          for(size_t ii=0;ii<size;ii++) {
            
            double yact;
            if (this->rescaled) {
              yact=this->y_r[ii];
            } else {
              yact=(*this->py)[ii];
            }

            // Compute sigma and ypred from Eq. 5.12
            double sigma2=1.0/inv_KXX2(ii,ii);
            double ypred=yact-Kinvf2[ii]*sigma2;

            // Then use Eq. 5.10
            qual+=pow(yact-ypred,2.0)/sigma2/2.0;
            qual+=0.5*log(sigma2);
          }
	  
        }
        
        if (verbose>0) {
          std::cout << "len,qual (loo_cv): " << len << " "
                    << qual << std::endl;
        }
      
      } else if (mode==mode_max_lml) {

        // Construct the KXX matrix
        mat_t KXX(size,size);
        for(size_t irow=0;irow<size;irow++) {
          for(size_t icol=0;icol<size;icol++) {
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else {
              KXX(irow,icol)=covar((*this->px)[irow],
                                   (*this->px)[icol]);
              if (irow==icol) KXX(irow,icol)+=noise_var;
            }
          }
        }

        // Compute the additive inverse of the log of the marginal
        // likelihood, from Eq. 5.8 of R&W, without the constant term
        
        double lndet;
        this->mi.invert_det(size,KXX,this->inv_KXX,lndet);
        lndet=log(lndet);
      
        // Inverse covariance matrix times function vector
        this->Kinvf.resize(size);

        qual=0.0;
        
        if (this->rescaled) {
          o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                             o2scl_cblas::o2cblas_NoTrans,
                             size,size,1.0,this->inv_KXX,this->y_r,0.0,
                             this->Kinvf);
          
          // Compute the log of the marginal likelihood, without
          // the constant term
          for(size_t i=0;i<size;i++) {
            qual+=0.5*this->y_r[i]*this->Kinvf[i];
          }
          qual+=0.5*lndet;
          
        } else {

          o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                             o2scl_cblas::o2cblas_NoTrans,
                             size,size,1.0,this->inv_KXX,*this->py,0.0,
                             this->Kinvf);
          
          // Compute the log of the marginal likelihood, without
          // the constant term
          for(size_t i=0;i<size;i++) {
            qual+=0.5*(*this->py)[i]*this->Kinvf[i];
          }
          qual+=0.5*lndet;

        }

        // Sometimes the calculation fails, and this helps avoid
        // those regions
        if (qual==-std::numeric_limits<double>::infinity()) {
          qual=std::numeric_limits<double>::infinity();
        }
        
        if (verbose>0) {
          std::cout << "len,qual (lml): " << len << " " << qual << std::endl;
        }
      }

      return qual;
    }
  
    interp_krige_optim() {
      nlen=20;
      full_min=false;
      mp=&def_min;
      def_min.ntrial*=10;
      verbose=0;
      mode=mode_loo_cv;
    }

    /// \name Function to minimize and various option
    //@{
    /// Leave-one-out cross validation (brute force version)
    static const size_t mode_loo_cv_bf=1;
    /// Maximize Log-marginal-likelihood
    static const size_t mode_max_lml=2;
    /// New leave-one-out cross validation method (default)
    static const size_t mode_loo_cv=3;
    /// Function to minimize (default \ref mode_loo_cv)
    size_t mode;
    ///@}
    
    /// Verbosity parameter
    int verbose;

    /** \brief Return the current length parameter
     */
    double get_length() {
      return len;
    }

    /** \brief Set the current length parameter
     */
    void set_length(double user_len) {
      if (user_len<=0.0) {
        O2SCL_ERR("Length cannot be less than or equal to zero.",
                  o2scl::exc_einval);
      }
      len=user_len;
      return;
    }
    
    /** \brief Number of length scale points to try when full minimizer 
        is not used (default 20)
    */
    size_t nlen;

    /// Default minimizer
    min_brent_gsl<> def_min;

    /// If true, use the full minimizer
    bool full_min;
    
    /** \brief Initialize interpolation routine with optional
        rescaling, and user-specified noise parameter and length scale

        No optimization of the length scale is performed.
     */
    virtual int set_noise_len(size_t size, const vec_t &x, const vec2_t &y,
                              double noise_var, double user_len,
                              bool rescale=false) {

      len=user_len;
      
      ff=std::bind(std::mem_fn<double(double,double)>
                   (&interp_krige_optim
                    <vec_t,vec2_t,mat_t,mat_inv_t>::covar),this,
                   std::placeholders::_1,std::placeholders::_2);
      ffd=std::bind(std::mem_fn<double(double,double)>
                    (&interp_krige_optim
                     <vec_t,vec2_t,mat_t,mat_inv_t>::deriv_covar),this,
                    std::placeholders::_1,std::placeholders::_2);
      ffd2=std::bind(std::mem_fn<double(double,double)>
                     (&interp_krige_optim
                      <vec_t,vec2_t,mat_t,mat_inv_t>::deriv2_covar),this,
                     std::placeholders::_1,std::placeholders::_2);
      ffi=std::bind(std::mem_fn<double(double,double,double)>
                    (&interp_krige_optim
                     <vec_t,vec2_t,mat_t,mat_inv_t>::integ_covar),this,
                    std::placeholders::_1,std::placeholders::_2,
                    std::placeholders::_3);
      
      this->set_covar_di_noise(size,x,y,ff,ffd,
                               ffd2,ffi,noise_var,this->rescaled);
      
      return 0;
    }
    
    /** \brief Initialize using optimization to find best length scale
        with optional rescaling and user-specified noise parameter
     */
    virtual int set_noise(size_t size, const vec_t &x, const vec2_t &y,
                          double noise_var, bool rescale=false) {

      // Set parent data members
      this->px=&x;
      this->py=&y;
      this->sz=size;

      if (rescale==true) {
        this->rescaled=true;
        this->mean_y=o2scl::vector_mean(size,y);
        this->std_y=o2scl::vector_stddev(size,y,this->mean_y);
        this->y_r.resize(size);
        for (size_t j=0;j<size;j++) {
          this->y_r[j]=(y[j]-this->mean_y)/this->std_y;
        }
      } else {
        this->rescaled=false;
      }

      int success=0;
      
      if (full_min) {

        if (verbose>1) {
          std::cout << "Class interp_krige_optim: full minimization, "
                    << std::endl;
          if (mode==mode_loo_cv_bf) {
            std::cout << "  leave one-out cross validation (brute force). ";
          } else if (mode==mode_loo_cv) {
            std::cout << "  leave one-out cross validation. ";
          } else {
            std::cout << "  log marginal likelihood. ";
          }
          std::cout << std::endl;
        }
      
        // Choose first interval as initial guess
        double len_opt=x[1]-x[0];

        funct mf=std::bind
          (std::mem_fn<double(double,double,int &)>
           (&interp_krige_optim<vec_t,vec2_t,mat_t,mat_inv_t>::qual_fun),
           this,std::placeholders::_1,noise_var,std::ref(success));
           
      
        mp->min(len_opt,qual,mf);
        len=len_opt;

        if (success!=0) {
          O2SCL_CONV2("Minimization failed in ",
                      "interp_krige_optim::set_noise().",
                      o2scl::exc_efailed,this->err_nonconv);
        }

      } else {

        if (verbose>1) {
          std::cout << "Class interp_krige_optim: simple minimization, "
                    << std::endl;
          if (mode==mode_loo_cv_bf) {
            std::cout << "  leave one-out cross validation (brute force). ";
          } else if (mode==mode_loo_cv) {
            std::cout << "  leave one-out cross validation. ";
          } else {
            std::cout << "  maximize marginal likelihood. ";
          }
          std::cout << std::endl;
        }

        // Compute a finite-difference array
        std::vector<double> diff(size-1);
        for(size_t i=0;i<size-1;i++) {
          diff[i]=fabs(x[i+1]-x[i]);
        }
      
        // Range of the length parameter
        double len_min=o2scl::vector_min_value
          <std::vector<double>,double>(size-1,diff)/3.0;
        double len_max;
        if (this->rescaled) {
          len_max=fabs((*this->px)[size-1]-(*this->px)[0])*3.0;
        } else {
          len_max=fabs(x[size-1]-x[0])*3.0;
        }
        double len_ratio=len_max/len_min;

        if (verbose>1) {
          std::cout << "             len (min,max,ratio): "
                    << len_min << " " << len_max << " "
                    << pow(len_ratio,((double)1)/((double)nlen-1))
                    << std::endl;
        }
	
        // Initialize to zero to prevent uninit'ed var. warnings
        double min_qual=0.0, len_opt=0.0;
      
        if (verbose>1) {
          std::cout << "             "
                    << "ilen qual len fail min_qual best_len"
                    << std::endl;
        }

        // Loop over the full range, finding the optimum
        bool min_set=false;
        for(size_t j=0;j<nlen;j++) {
          len=len_min*pow(len_ratio,((double)j)/((double)nlen-1));

          int success=0;
          qual=qual_fun(len,noise_var,success);
	
          if (success==0 && (min_set==false || qual<min_qual)) {
            len_opt=len;
            min_qual=qual;
            min_set=true;
          }
	
          if (verbose>1) {
            std::cout << "krige_optim: ";
            std::cout.width(2);
            std::cout << j << " ";
            std::cout.setf(std::ios::showpos);
            std::cout << qual << " " << len << " "
                      << success << " " << min_qual << " "
                      << len_opt << std::endl;
            std::cout.unsetf(std::ios::showpos);
          }
	  
        }

        // Now that we've optimized the covariance function,
        // just use the parent class to interpolate
        len=len_opt;

      }

      if (verbose>0) {
        std::cout << "interp_krige_optim::set_noise(): optimum length: "
                  << len << std::endl;
      }
      
      ff=std::bind(std::mem_fn<double(double,double)>
                   (&interp_krige_optim
                    <vec_t,vec2_t,mat_t,mat_inv_t>::covar),this,
                   std::placeholders::_1,std::placeholders::_2);
      ffd=std::bind(std::mem_fn<double(double,double)>
                    (&interp_krige_optim
                     <vec_t,vec2_t,mat_t,mat_inv_t>::deriv_covar),this,
                    std::placeholders::_1,std::placeholders::_2);
      ffd2=std::bind(std::mem_fn<double(double,double)>
                     (&interp_krige_optim
                      <vec_t,vec2_t,mat_t,mat_inv_t>::deriv2_covar),this,
                     std::placeholders::_1,std::placeholders::_2);
      ffi=std::bind(std::mem_fn<double(double,double,double)>
                    (&interp_krige_optim
                     <vec_t,vec2_t,mat_t,mat_inv_t>::integ_covar),this,
                    std::placeholders::_1,std::placeholders::_2,
                    std::placeholders::_3);
      
      this->set_covar_di_noise(size,x,y,ff,ffd,
                               ffd2,ffi,noise_var,this->rescaled);
      
      return 0;
    }

    /** \brief Initialize interpolation with default noise parameter
        and no rescaling
     */
    virtual void set(size_t size, const vec_t &x, const vec2_t &y) {

      // Use the mean absolute value to determine noise
      double mean_abs=0.0;
      for(size_t j=0;j<size;j++) {
        mean_abs+=fabs(y[j]);
      }
      mean_abs/=size;

      set_noise(size,x,y,mean_abs/1.0e8,false);
    
      return;
    }
  
    /** \brief Initialize interpolation with default noise parameter
        and optional rescaling
     */
    virtual void set(size_t size, const vec_t &x, const vec2_t &y,
                     bool rescale) {

      // Use the mean absolute value to determine noise
      double mean_abs=0.0;
      for(size_t j=0;j<size;j++) {
        mean_abs+=fabs(y[j]);
      }
      mean_abs/=size;

      set_noise(size,x,y,mean_abs/1.0e8,rescale);
    
      return;
    }
  
  
  };

#ifdef O2SCL_NEVER_DEFINED
  
  template<class vec_t, class vec2_t=vec_t, class func_t=covar_funct,
           class mat_t=boost::numeric::ublas::matrix<double>,
           class mat_inv_t=o2scl_linalg::matrix_invert_det_cholesky<mat_t> >
  class interp_krige_optim_new :
    public interp_krige<vec_t,vec2_t,
                        std::function<double(double,double)>,
                        std::function<double(double,double,double)>,
                        mat_t,mat_inv_t> {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

  protected:

    /// Desc
    func_t cf;
  
    /// The quality factor of the optimization
    double qual;
    
    /// Pointer to the user-specified minimizer
    mmin_base<> *mp;

  public:
    
    /** \brief Function to optimize the covariance parameters
     */
    double qual_fun(double x, double noise_var, int &success) {

      len=x;
      success=0;

      size_t size=this->sz;

      if (mode==mode_loo_cv_bf) {
      
        qual=0.0;
        for(size_t k=0;k<size;k++) {
	
          // Leave one observation out
          ubvector x2(size-1);
          ubvector y2(size-1);
          o2scl::vector_copy_jackknife((*this->px),k,x2);
          if (this->rescaled) {
            o2scl::vector_copy_jackknife(this->y_r,k,y2);
          } else {
            o2scl::vector_copy_jackknife((*this->py),k,y2);
          }
	
          // Construct the inverse of the KXX matrix. Note
          // that we don't use the inv_KXX data member
          // because of the size mismatch
          mat_t inv_KXX2(size-1,size-1);
          for(size_t irow=0;irow<size-1;irow++) {
            for(size_t icol=0;icol<size-1;icol++) {
              if (irow>icol) {
                inv_KXX2(irow,icol)=inv_KXX2(icol,irow);
              } else {
                inv_KXX2(irow,icol)=cf(x2[irow],x2[icol]);
                if (irow==icol) inv_KXX2(irow,icol)+=noise_var;
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
          double yact;
          if (this->rescaled) {
            yact=this->y_r[k];
          } else {
            yact=(*this->py)[k];
          }

          // Compute the predicted value
          double ypred=0.0;
          ubvector kxx0(size-1);
          if (this->rescaled) {
            for(size_t i=0;i<size-1;i++) {
              kxx0[i]=cf((*this->px)[k],x2[i]);
              ypred+=kxx0[i]*Kinvf2[i];
              if (false && k==0) {
                std::cout << "  " << Kinvf2[i] << " " << x2[i] << " "
                          << kxx0[i] << " " << ypred
                          << std::endl;
              }
            }
          } else {
            for(size_t i=0;i<size-1;i++) {
              kxx0[i]=cf((*this->px)[k],x2[i]);
              ypred+=kxx0[i]*Kinvf2[i];
            }
          }
          
          if (false && k==0) {
            std::cout.setf(std::ios::showpos);
            std::cout << "k,x,yact,ypred: " << k << " "
                      << (*this->px)[k] << " "
                      << yact << " " << ypred << std::endl;
            std::cout.unsetf(std::ios::showpos);
          }
          
          qual+=pow(yact-ypred,2.0);
	
        }
        if (verbose>0) {
          std::cout << "len,qual (loo_cv_bf): " << len << " "
                    << qual << std::endl;
          //std::cout << "--" << std::endl;
        }

      } else if (mode==mode_loo_cv) {

#ifdef O2SCL_NEVER_DEFINED
        
        if (false) {
          
          // AWS, 7/30/22: This is the loo_cv
          // method to compute mu and sigma and Eq 5.10 in R&W., I
          // think it should work, but it doesn't quite work, however,
          // using Eq. 5.12 to compute mu and sigma (below) seems to
          // work well (and is more efficient), so I'm keeping that
          
          qual=0.0;
          for(size_t k=0;k<size;k++) {
	
            // Leave one observation out
            ubvector x2(size-1);
            ubvector y2(size-1);
            o2scl::vector_copy_jackknife((*this->px),k,x2);
            if (this->rescaled) {
              o2scl::vector_copy_jackknife(this->y_r,k,y2);
            } else {
              o2scl::vector_copy_jackknife((*this->py),k,y2);
            }
	
            // Construct the KXX matrix
            mat_t inv_KXX2(size-1,size-1);
            for(size_t irow=0;irow<size-1;irow++) {
              for(size_t icol=0;icol<size-1;icol++) {
                if (irow>icol) {
                  inv_KXX2(irow,icol)=inv_KXX2(icol,irow);
                } else {
                  inv_KXX2(irow,icol)=cf(x2[irow],x2[icol]);
                  if (irow==icol) inv_KXX2(irow,icol)+=noise_var;
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
	  
            double yact;
            if (this->rescaled) {
              yact=this->y_r[k];
            } else {
              yact=(*this->py)[k];
            }
          
            double ypred=0.0;
            ubvector kxx0(size-1);
            if (this->rescaled) {
              for(size_t i=0;i<size-1;i++) {
                kxx0[i]=cf((*this->px)[k],x2[i]);
                ypred+=kxx0[i]*Kinvf2[i];
              }
            } else {
              for(size_t i=0;i<size-1;i++) {
                kxx0[i]=cf((*this->px)[k],x2[i]);
                ypred+=kxx0[i]*Kinvf2[i];
              }
            }
          
            double kx0x0=cf((*this->px)[k],(*this->px)[k]);
            ubvector prod(size-1);
          
            boost::numeric::ublas::axpy_prod(inv_KXX2,kxx0,prod,true);
            double sigma=kx0x0-
              boost::numeric::ublas::inner_prod(kxx0,prod);

            if (this->rescaled) {
              sigma*=this->std_y;
            }

            if (false && k==0) {
              std::cout.setf(std::ios::showpos);
              std::cout << "k,x,yact,ypred,sigma: " << k << " "
                        << (*this->px)[k] << " "
                        << yact << " " << ypred << " "
                        << sigma << std::endl;
              std::cout.unsetf(std::ios::showpos);
            }

            // We maximize the predictive log probability, Eq 5.10
            // in R&W
            qual+=pow(yact-ypred,2.0)/sigma/sigma/2.0;
            qual+=0.5*log(sigma*sigma);
	
          }

        }
        
#endif

        if (true) {
          
          // Construct the KXX matrix
          mat_t inv_KXX2(size,size);
          for(size_t irow=0;irow<size;irow++) {
            for(size_t icol=0;icol<size;icol++) {
              if (irow>icol) {
                inv_KXX2(irow,icol)=inv_KXX2(icol,irow);
              } else {
                inv_KXX2(irow,icol)=cf((*this->px)[irow],
                                          (*this->px)[icol]);
                if (irow==icol) inv_KXX2(irow,icol)+=noise_var;
              }
            }
          }
          
          // Construct the inverse of KXX
          this->mi.invert_inplace(size-1,inv_KXX2);
	  
          // Inverse covariance matrix times function vector
          ubvector Kinvf2(size);
          if (this->rescaled) {
            o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                               o2scl_cblas::o2cblas_NoTrans,
                               size,size,1.0,inv_KXX2,this->y_r,0.0,
                               Kinvf2);
          } else {
            o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                               o2scl_cblas::o2cblas_NoTrans,
                               size,size,1.0,inv_KXX2,*(this->py),0.0,
                               Kinvf2);
          }

          qual=0.0;
          for(size_t ii=0;ii<size;ii++) {
            
            double yact;
            if (this->rescaled) {
              yact=this->y_r[ii];
            } else {
              yact=(*this->py)[ii];
            }

            // Compute sigma and ypred from Eq. 5.12
            double sigma2=1.0/inv_KXX2(ii,ii);
            double ypred=yact-Kinvf2[ii]*sigma2;

            // Then use Eq. 5.10
            qual+=pow(yact-ypred,2.0)/sigma2/2.0;
            qual+=0.5*log(sigma2);
          }
	  
        }
        
        if (verbose>0) {
          std::cout << "len,qual (loo_cv): " << len << " "
                    << qual << std::endl;
        }
      
      } else if (mode==mode_max_lml) {

        // Construct the KXX matrix
        mat_t KXX(size,size);
        for(size_t irow=0;irow<size;irow++) {
          for(size_t icol=0;icol<size;icol++) {
            if (irow>icol) {
              KXX(irow,icol)=KXX(icol,irow);
            } else {
              KXX(irow,icol)=cf((*this->px)[irow],
                                   (*this->px)[icol]);
              if (irow==icol) KXX(irow,icol)+=noise_var;
            }
          }
        }

        // Compute the additive inverse of the log of the marginal
        // likelihood, from Eq. 5.8 of R&W, without the constant term
        
        double lndet;
        this->mi.invert_det(size,KXX,this->inv_KXX,lndet);
        lndet=log(lndet);
      
        // Inverse covariance matrix times function vector
        this->Kinvf.resize(size);

        qual=0.0;
        
        if (this->rescaled) {
          o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                             o2scl_cblas::o2cblas_NoTrans,
                             size,size,1.0,this->inv_KXX,this->y_r,0.0,
                             this->Kinvf);
          
          // Compute the log of the marginal likelihood, without
          // the constant term
          for(size_t i=0;i<size;i++) {
            qual+=0.5*this->y_r[i]*this->Kinvf[i];
          }
          qual+=0.5*lndet;
          
        } else {

          o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
                             o2scl_cblas::o2cblas_NoTrans,
                             size,size,1.0,this->inv_KXX,*this->py,0.0,
                             this->Kinvf);
          
          // Compute the log of the marginal likelihood, without
          // the constant term
          for(size_t i=0;i<size;i++) {
            qual+=0.5*(*this->py)[i]*this->Kinvf[i];
          }
          qual+=0.5*lndet;

        }

        // Sometimes the calculation fails, and this helps avoid
        // those regions
        if (qual==-std::numeric_limits<double>::infinity()) {
          qual=std::numeric_limits<double>::infinity();
        }
        
        if (verbose>0) {
          std::cout << "len,qual (lml): " << len << " " << qual << std::endl;
        }
      }

      return qual;
    }
  
    interp_krige_optim_new() {
      nlen=20;
      full_min=false;
      mp=&def_min;
      def_min.ntrial*=10;
      verbose=0;
      mode=mode_loo_cv;
    }

    /// \name Function to minimize and various option
    //@{
    /// Leave-one-out cross validation (brute force version)
    static const size_t mode_loo_cv_bf=1;
    /// Maximize Log-marginal-likelihood
    static const size_t mode_max_lml=2;
    /// New leave-one-out cross validation method (default)
    static const size_t mode_loo_cv=3;
    /// Function to minimize (default \ref mode_loo_cv)
    size_t mode;
    ///@}
    
    /// Verbosity parameter
    int verbose;

    /** \brief Set the current length parameter
     */
    void set_length(double user_len) {
      if (user_len<=0.0) {
        O2SCL_ERR("Length cannot be less than or equal to zero.",
                  o2scl::exc_einval);
      }
      len=user_len;
      return;
    }
    
    /** \brief Number of length scale points to try when full minimizer 
        is not used (default 20)
    */
    size_t nlen;

    /// Default minimizer
    mmin_simp2<> def_mmin;

    /// If true, use the full minimizer
    bool full_min;
    
    /** \brief Initialize interpolation routine with optional
        rescaling, and user-specified noise parameter and length scale

        No optimization of the length scale is performed.
     */
    virtual int set_noise_len(size_t size, const vec_t &x, const vec2_t &y,
                              double noise_var, double user_len,
                              bool rescale=false) {

      len=user_len;
      
      ff=std::bind(std::mem_fn<double(double,double)>
                   (&interp_krige_optim_new
                    <vec_t,vec2_t,mat_t,mat_inv_t>::covar),this,
                   std::placeholders::_1,std::placeholders::_2);
      ffd=std::bind(std::mem_fn<double(double,double)>
                    (&interp_krige_optim_new
                     <vec_t,vec2_t,mat_t,mat_inv_t>::deriv_covar),this,
                    std::placeholders::_1,std::placeholders::_2);
      ffd2=std::bind(std::mem_fn<double(double,double)>
                     (&interp_krige_optim_new
                      <vec_t,vec2_t,mat_t,mat_inv_t>::deriv2_covar),this,
                     std::placeholders::_1,std::placeholders::_2);
      ffi=std::bind(std::mem_fn<double(double,double,double)>
                    (&interp_krige_optim_new
                     <vec_t,vec2_t,mat_t,mat_inv_t>::integ_covar),this,
                    std::placeholders::_1,std::placeholders::_2,
                    std::placeholders::_3);
      
      this->set_covar_di_noise(size,x,y,ff,ffd,
                               ffd2,ffi,noise_var,this->rescaled);
      
      return 0;
    }
    
    /** \brief Initialize using optimization to find best length scale
        with optional rescaling and user-specified noise parameter
     */
    virtual int set_noise(size_t size, const vec_t &x, const vec2_t &y,
                          double noise_var, bool rescale=false) {

      // Set parent data members
      this->px=&x;
      this->py=&y;
      this->sz=size;

      if (rescale==true) {
        this->rescaled=true;
        this->mean_y=o2scl::vector_mean(size,y);
        this->std_y=o2scl::vector_stddev(size,y,this->mean_y);
        this->y_r.resize(size);
        for (size_t j=0;j<size;j++) {
          this->y_r[j]=(y[j]-this->mean_y)/this->std_y;
        }
      } else {
        this->rescaled=false;
      }

      int success=0;
      
      if (full_min) {

        if (verbose>1) {
          std::cout << "Class interp_krige_optim_new: full minimization, "
                    << std::endl;
          if (mode==mode_loo_cv_bf) {
            std::cout << "  leave one-out cross validation (brute force). ";
          } else if (mode==mode_loo_cv) {
            std::cout << "  leave one-out cross validation. ";
          } else {
            std::cout << "  log marginal likelihood. ";
          }
          std::cout << std::endl;
        }
      
        // Choose first interval as initial guess
        double len_opt=x[1]-x[0];

        funct mf=std::bind
          (std::mem_fn<double(double,double,int &)>
           (&interp_krige_optim_new<vec_t,vec2_t,mat_t,mat_inv_t>::qual_fun),
           this,std::placeholders::_1,noise_var,std::ref(success));
           
      
        mp->min(len_opt,qual,mf);
        len=len_opt;

        if (success!=0) {
          O2SCL_CONV2("Minimization failed in ",
                      "interp_krige_optim_new::set_noise().",
                      o2scl::exc_efailed,this->err_nonconv);
        }

      } else {

        if (verbose>1) {
          std::cout << "Class interp_krige_optim_new: simple minimization, "
                    << std::endl;
          if (mode==mode_loo_cv_bf) {
            std::cout << "  leave one-out cross validation (brute force). ";
          } else if (mode==mode_loo_cv) {
            std::cout << "  leave one-out cross validation. ";
          } else {
            std::cout << "  maximize marginal likelihood. ";
          }
          std::cout << std::endl;
        }

        // Compute a finite-difference array
        std::vector<double> diff(size-1);
        for(size_t i=0;i<size-1;i++) {
          diff[i]=fabs(x[i+1]-x[i]);
        }
      
        // Range of the length parameter
        double len_min=o2scl::vector_min_value
          <std::vector<double>,double>(size-1,diff)/3.0;
        double len_max;
        if (this->rescaled) {
          len_max=fabs((*this->px)[size-1]-(*this->px)[0])*3.0;
        } else {
          len_max=fabs(x[size-1]-x[0])*3.0;
        }
        double len_ratio=len_max/len_min;

        if (verbose>1) {
          std::cout << "             len (min,max,ratio): "
                    << len_min << " " << len_max << " "
                    << pow(len_ratio,((double)1)/((double)nlen-1))
                    << std::endl;
        }
	
        // Initialize to zero to prevent uninit'ed var. warnings
        double min_qual=0.0, len_opt=0.0;
      
        if (verbose>1) {
          std::cout << "             "
                    << "ilen qual len fail min_qual best_len"
                    << std::endl;
        }

        // Loop over the full range, finding the optimum
        bool min_set=false;
        for(size_t j=0;j<nlen;j++) {
          len=len_min*pow(len_ratio,((double)j)/((double)nlen-1));

          int success=0;
          qual=qual_fun(len,noise_var,success);
	
          if (success==0 && (min_set==false || qual<min_qual)) {
            len_opt=len;
            min_qual=qual;
            min_set=true;
          }
	
          if (verbose>1) {
            std::cout << "krige_optim: ";
            std::cout.width(2);
            std::cout << j << " ";
            std::cout.setf(std::ios::showpos);
            std::cout << qual << " " << len << " "
                      << success << " " << min_qual << " "
                      << len_opt << std::endl;
            std::cout.unsetf(std::ios::showpos);
          }
	  
        }

        // Now that we've optimized the covariance function,
        // just use the parent class to interpolate
        len=len_opt;

      }

      if (verbose>0) {
        std::cout << "interp_krige_optim_new::set_noise(): optimum length: "
                  << len << std::endl;
      }
      
      ff=std::bind(std::mem_fn<double(double,double)>
                   (&interp_krige_optim_new
                    <vec_t,vec2_t,mat_t,mat_inv_t>::covar),this,
                   std::placeholders::_1,std::placeholders::_2);
      ffd=std::bind(std::mem_fn<double(double,double)>
                    (&interp_krige_optim_new
                     <vec_t,vec2_t,mat_t,mat_inv_t>::deriv_covar),this,
                    std::placeholders::_1,std::placeholders::_2);
      ffd2=std::bind(std::mem_fn<double(double,double)>
                     (&interp_krige_optim_new
                      <vec_t,vec2_t,mat_t,mat_inv_t>::deriv2_covar),this,
                     std::placeholders::_1,std::placeholders::_2);
      ffi=std::bind(std::mem_fn<double(double,double,double)>
                    (&interp_krige_optim_new
                     <vec_t,vec2_t,mat_t,mat_inv_t>::integ_covar),this,
                    std::placeholders::_1,std::placeholders::_2,
                    std::placeholders::_3);
      
      this->set_covar_di_noise(size,x,y,ff,ffd,
                               ffd2,ffi,noise_var,this->rescaled);
      
      return 0;
    }

    /** \brief Initialize interpolation with default noise parameter
        and no rescaling
     */
    virtual void set(size_t size, const vec_t &x, const vec2_t &y) {

      // Use the mean absolute value to determine noise
      double mean_abs=0.0;
      for(size_t j=0;j<size;j++) {
        mean_abs+=fabs(y[j]);
      }
      mean_abs/=size;

      set_noise(size,x,y,mean_abs/1.0e8,false);
    
      return;
    }
  
    /** \brief Initialize interpolation with default noise parameter
        and optional rescaling
     */
    virtual void set(size_t size, const vec_t &x, const vec2_t &y,
                     bool rescale) {

      // Use the mean absolute value to determine noise
      double mean_abs=0.0;
      for(size_t j=0;j<size;j++) {
        mean_abs+=fabs(y[j]);
      }
      mean_abs/=size;

      set_noise(size,x,y,mean_abs/1.0e8,rescale);
    
      return;
    }
  
  
  };
  
#endif
  
}

#endif
