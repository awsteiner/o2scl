/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2015, Andrew W. Steiner

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
#ifndef O2SCL_GSL_MMIN_BFGS2_H
#define O2SCL_GSL_MMIN_BFGS2_H

/** \file mmin_bfgs2.h
    \brief File defining \ref o2scl::mmin_bfgs2
*/
#include <gsl/gsl_blas.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_multimin.h>
#include <o2scl/mmin.h>
#include <o2scl/cblas.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Virtual base for the mmin_bfgs2 wrapper

      This class is useful so that the gsl_mmin_linmin class doesn't
      need to depend on any template parameters, even though it will
      need a wrapping object as an argument for the
      gsl_mmin_linmin::min() function.
  */
  class mmin_wrap_gsl {
  public:
    virtual ~mmin_wrap_gsl() {}
    /// Function
    virtual double wrap_f(double alpha)=0;
    /// Derivative
    virtual double wrap_df(double alpha)=0;
    /// Function and derivative
    virtual void wrap_fdf(double alpha, double *f, double *df)=0;
  };

  /** \brief Wrapper class for the mmin_bfgs2 minimizer

      This is a reimplmentation of the internal GSL wrapper for
      function calls in the BFGS minimizer
  */
  template<class func_t=multi_funct11, 
    class vec_t=boost::numeric::ublas::vector<double> , 
    class dfunc_t=grad_funct11,
    class auto_grad_t=
    gradient<multi_funct11,boost::numeric::ublas::vector<double> > > 
    class mmin_wrapper_gsl : public mmin_wrap_gsl {
    
#ifndef DOXYGEN_INTERNAL
    
  protected:

  /// Function
  func_t *func;

  /// Derivative
  dfunc_t *dfunc;

  /// The automatic gradient object
  auto_grad_t *agrad;
    
  /// True if the gradient was given by the user
  bool grad_given;
    
  /** \name fixed values 
   */
  //@{
  vec_t *x;
  vec_t *g;
  vec_t *p;
  //@}
    
  /** \name cached values, for x(alpha)=x + alpha * p 
   */
  //@{
  double f_alpha;
  double df_alpha;
  //@}
    
  /** \name cache keys
   */
  //@{
  double f_cache_key;
  double df_cache_key;
  double x_cache_key;
  double g_cache_key;
  //@}
    
  /// Move to a new point, using the cached value if possible
  void moveto(double alpha) {

    if (alpha == x_cache_key) {
      /* using previously cached position */
      return;
    }
      
    /* set x_alpha=x + alpha * p */
    
    for(size_t i=0;i<dim;i++) {
      av_x_alpha[i]=(*x)[i];
    }
    o2scl_cblas::daxpy(alpha,dim,*p,av_x_alpha);

    x_cache_key=alpha;
  }

  /// Compute the slope
  double slope() {
    return o2scl_cblas::ddot(dim,av_g_alpha,*p);
  }

  /// Evaluate the function
  virtual double wrap_f(double alpha) {
    if (alpha==f_cache_key) {
      return f_alpha;
    }
    moveto(alpha);
    f_alpha=(*func)(dim,av_x_alpha);

    f_cache_key=alpha;
    return f_alpha;
  }

  /// Evaluate the derivative
  virtual double wrap_df(double alpha) {

    /* using previously cached df(alpha) */
    if (alpha==df_cache_key) return df_alpha;

    moveto(alpha);
    if (alpha!=g_cache_key) {
      if (grad_given) {
	(*dfunc)(dim,av_x_alpha,av_g_alpha);
      } else {
	(*agrad)(dim,av_x_alpha,av_g_alpha);
      }
      g_cache_key=alpha;
    }
    df_alpha=slope();
    df_cache_key=alpha;

    return df_alpha;
  }

  /// Evaluate the function and the derivative
  virtual void wrap_fdf(double alpha, double *f, double *df) {
    
    /* Check for previously cached values */
    if (alpha == f_cache_key && alpha == df_cache_key) {
      *f=f_alpha;
      *df=df_alpha;
      return;
    }
    if (alpha == f_cache_key || alpha == df_cache_key) {
      *f=wrap_f(alpha);
      *df=wrap_df(alpha);
      return;
    }
      
    moveto(alpha);
    f_alpha=(*func)(dim,av_x_alpha);
    if (grad_given) {
      (*dfunc)(dim,av_x_alpha,av_g_alpha);
    } else {
      (*agrad)(dim,av_x_alpha,av_g_alpha);
    }
    f_cache_key=alpha;
    g_cache_key=alpha;
      
    df_alpha=slope();
    df_cache_key=alpha;
      
    *f=f_alpha;
    *df=df_alpha;
      
    return;
  }
    
#endif

  public:

  /// Temporary storage
  vec_t av_x_alpha;

  /// Temporary storage
  vec_t av_g_alpha;

  /// Number of minimization dimensions
  size_t dim;
    
  /// Initialize wrapper
  void prepare_wrapper(func_t &ufunc, dfunc_t *udfunc, vec_t &t_x, 
		       double f, vec_t &t_g, vec_t &t_p, auto_grad_t *ag) {
      
    func=&ufunc;
    dfunc=udfunc;
    if (dfunc!=0) grad_given=true;
    else grad_given=false;
    agrad=ag;
      
    x=&t_x;
    g=&t_g;
    p=&t_p;
      
    x_cache_key=0.0;
    f_cache_key=0.0;
    g_cache_key=0.0;
    df_cache_key=0.0;

    for(size_t i=0;i<dim;i++) {
      av_x_alpha[i]=(*x)[i];
      av_g_alpha[i]=(*g)[i];
    }
      
    f_alpha=f;
    df_alpha=slope();

    return;
  }
    
  /// Update position
  void update_position(double alpha, vec_t &t_x, double *t_f, vec_t &t_g) {
      
    /* ensure that everything is fully cached */
    { 
      double t_f_alpha, t_df_alpha; 
      wrap_fdf(alpha,&t_f_alpha,&t_df_alpha); 
    }
      
    *t_f=f_alpha;
    for(size_t i=0;i<dim;i++) {
      t_x[i]=av_x_alpha[i];
      t_g[i]=av_g_alpha[i];
    }
  }
    
  /** \brief Convert cache values to the new minimizer direction

      Convert the cache values from the end of the current minimisation
      to those needed for the start of the next minimisation, alpha=0
  */
  void change_direction() {
      
    /* The new x_alpha for alpha=0 is the current position and the
       new g_alpha for alpha=0 is the current gradient at the
       endpoint
    */
    for(size_t i=0;i<dim;i++) {
      av_x_alpha[i]=(*x)[i];
      av_g_alpha[i]=(*g)[i];
    }
    x_cache_key=0.0;
    g_cache_key=0.0;
      
    /* The function value does not change */
    f_cache_key=0.0;
      
    /* Calculate the slope along the new direction vector, p */
    df_alpha=slope();
    df_cache_key=0.0;
      
    return;
  }
    
  };

  /** \brief The line minimizer for mmin_bfgs2
   */
  class mmin_linmin_gsl {

#ifndef DOXYGEN_INTERNAL

  protected:

    /** \brief Minimize the interpolating quadratic
	
	Find a minimum in x=[0,1] of the interpolating quadratic through
	(0,f0) (1,f1) with derivative fp0 at x=0.  The interpolating
	polynomial is q(x)=f0 + fp0 * z + (f1-f0-fp0) * z^2
    */
    double interp_quad(double f0, double fp0, double f1, double zl, 
		       double zh);
    
    /** \brief Minimize the interpolating cubic
	
	Find a minimum in x=[0,1] of the interpolating cubic through
	(0,f0) (1,f1) with derivatives fp0 at x=0 and fp1 at x=1.
	
	The interpolating polynomial is:
	
	c(x)=f0 + fp0 * z + eta * z^2 + xi * z^3
	
	where eta=3*(f1-f0)-2*fp0-fp1, xi=fp0+fp1-2*(f1-f0).
    */
    double cubic(double c0, double c1, double c2, double c3, double z);
    
    /// Test to see curvature is positive
    void check_extremum(double c0, double c1, double c2, double c3, double z,
			double *zmin, double *fmin);
    
    /// Interpolate using a cubic
    double interp_cubic(double f0, double fp0, double f1, 
			double fp1, double zl, double zh);
    
    /// Perform the interpolation
    double interpolate(double a, double fa, double fpa, double b, 
		       double fb, double fpb, double xmin, double xmax,
		       int order);
#endif
    
  public:

    /** \brief The line minimization
	
	Recommended values from \ref Fletcher87 are
	<tt>rho=0.01, sigma=0.1, tau1=9, tau2=0.05, 
	tau3=0.5 </tt>
    */
    int minimize(mmin_wrap_gsl &wrap, double rho, double sigma,
		 double tau1, double tau2, double tau3,
		 int order, double alpha1, double *alpha_new);
  };

  /** \brief Multidimensional minimization by the BFGS
      algorithm (GSL)
      
      The functions mmin() and mmin_de() min a given function
      until the gradient is smaller than the value of mmin::tol_rel
      (which defaults to \f$ 10^{-4} \f$ ).

      See an example for the usage of this class in 
      \ref ex_mmin_sect .
      
      This class includes the optimizations from the GSL minimizer \c
      vector_bfgs2.

      Default template arguments
      - \c func_t - \ref multi_funct11
      - \c vec_t - \ref boost::numeric::ublas::vector \<double \> 
      - \c dfunc_t - \ref mm_funct11
      - \c auto_grad_t - \ref gradient\<func_t, 
      \ref boost::numeric::ublas::vector \<double \> \>
      - \c def_auto_grad_t - \ref gradient_gsl\<func_t,
      \ref boost::numeric::ublas::vector \< double \> \>

      \todo While BFGS does well in the \c ex_mmin example with the
      initial guess of \f$ (1,0,7\pi) \f$ it seems to converge more
      poorly for the spring function than the other minimizers with
      other initial guesses, and I think this will happen in the GSL
      versions too. I need to examine this more closely with some code
      designed to clearly show this.
  */
  template<class func_t=multi_funct11, 
    class vec_t=boost::numeric::ublas::vector<double> , 
    class dfunc_t=grad_funct11,
    class auto_grad_t=
    gradient<multi_funct11,boost::numeric::ublas::vector<double> >,
    class def_auto_grad_t=
    gradient_gsl<multi_funct11,boost::numeric::ublas::vector<double> > > 
    class mmin_bfgs2 : public mmin_base<func_t,func_t,vec_t> {
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
  
  /// \name The original variables from the GSL state structure
  //@{
  int iter;
  double step;
  double g0norm;
  double pnorm;
  double delta_f;
  /* f'(0) for f(x-alpha*p) */
  double fp0;                   
  vec_t x0;
  vec_t g0;
  vec_t p;
  /* work space */
  vec_t dx0;
  vec_t dg0;
  /* wrapper function */
  mmin_wrapper_gsl<func_t,vec_t,dfunc_t,auto_grad_t> wrap;
  /* minimization parameters */
  double rho;
  double sigma;
  double tau1;
  double tau2;
  double tau3;
  int order;
  //@}
  
  /// The line minimizer
  mmin_linmin_gsl lm;
      
  /// \name Store the arguments to set() so we can use them for iterate()
  //@{
  vec_t *st_x;
  vec_t st_dx;
  vec_t st_grad;
  double st_f;
  //@}

  /// Memory size
  size_t dim;

  /// Automatic gradient object
  auto_grad_t *agrad;
    
#endif
    
  public:
     
  mmin_bfgs2() {
    // The default values of lmin_tol and step_size from the
    // example in the GSL documentation
    lmin_tol=1.0e-4;
    step_size=0.01;
    agrad=&def_grad;
  }
    
  virtual ~mmin_bfgs2() {}

  /// Perform an iteration
  virtual int iterate() {

    double alpha=0.0, alpha1;

    double pg, dir;
    int status;
      
    double f0=st_f;
      
    if (pnorm == 0.0 || g0norm == 0.0 || fp0 == 0) {
      for(size_t i=0;i<this->dim;i++) st_dx[i]=0.0;
      O2SCL_CONV2("Either pnorm, g0norm, or fp0 vanished in ",
		  "mmin_bfgs2::iterate().",
		  exc_enoprog,this->err_nonconv);
      // The functions mmin() and mmin_de() use this return value, so
      // when err_nonconv is false, we still want to return a non-zero
      // value.
      return exc_enoprog;
    }
      
      double dbl_eps=std::numeric_limits<double>::epsilon();

    if (delta_f < 0) {
      double del;
      if (-delta_f>10.0*dbl_eps*fabs(f0)) del=-delta_f;
      else del=10.0*dbl_eps*fabs(f0);
      if (2.0*del/(-fp0)<1.0) alpha1=2.0*del/(-fp0);
      else alpha1=1.0;
    } else {
      alpha1=fabs(step);
    }
	
    /* line minimisation, with cubic interpolation (order=3) */
      
    status=lm.minimize(wrap,rho,sigma,tau1,tau2,tau3,order,
		       alpha1,&alpha);
      
    if (status != success) {
      O2SCL_CONV2("Variable stepb vanished in ",
		  "mmin_bfgs2::iterate().",
		  exc_enoprog,this->err_nonconv);
      // The functions mmin() and mmin_de() use this return value, so
      // when err_nonconv is false, we still want to return a non-zero
      // value.
      return exc_enoprog;
    }
      
    wrap.update_position(alpha,*st_x,&st_f,st_grad);

    delta_f=st_f - f0;

    /* Choose a new direction for the next step */

    {
      /* This is the BFGS update: */
      /* p'=g1 - A dx - B dg */
      /* A=- (1+ dg.dg/dx.dg) B + dg.g/dx.dg */
      /* B=dx.g/dx.dg */

      double dxg, dgg, dxdg, dgnorm, A, B;

      /* dx0=x - x0 */
      for(size_t i=0;i<dim;i++) dx0[i]=(*st_x)[i];
      o2scl_cblas::daxpy(-1.0,dim,x0,dx0);

      /* keep a copy */
      for(size_t i=0;i<dim;i++) st_dx[i]=dx0[i];

      /* dg0=g - g0 */
      for(size_t i=0;i<dim;i++) dg0[i]=st_grad[i];
      o2scl_cblas::daxpy(-1.0,dim,g0,dg0);

      dxg=o2scl_cblas::ddot(dim,dx0,st_grad);
      dgg=o2scl_cblas::ddot(dim,dg0,st_grad);
      dxdg=o2scl_cblas::ddot(dim,dx0,dg0);
      
      dgnorm=o2scl_cblas::dnrm2(dim,dg0);

      if (dxdg != 0) {
	B=dxg / dxdg;
	A=-(1.0 + dgnorm * dgnorm / dxdg) * B + dgg / dxdg;
      } else {
	B=0;
	A=0;
      }
	
      for(size_t i=0;i<dim;i++) p[i]=st_grad[i];
      o2scl_cblas::daxpy(-A,dim,dx0,p);
      o2scl_cblas::daxpy(-B,dim,dg0,p);
    }
      
    for(size_t i=0;i<dim;i++) {
      g0[i]=st_grad[i];
      x0[i]=(*st_x)[i];
    }

    g0norm=o2scl_cblas::dnrm2(dim,g0);
    pnorm=o2scl_cblas::dnrm2(dim,p);

    /* update direction and fp0 */

    pg=o2scl_cblas::ddot(dim,st_grad,p);
    dir=(pg >= 0.0) ? -1.0 : +1.0;
    o2scl_cblas::dscal(dir/pnorm,dim,p);
    pnorm=o2scl_cblas::dnrm2(dim,p);
    fp0=o2scl_cblas::ddot(dim,p,g0);

    wrap.change_direction();

    return success;
  }
    
  /// Return string denoting type("mmin_bfgs2")
  virtual const char *type() { return "mmin_bfgs2";}

  /// Allocate the memory
  virtual int allocate(size_t n) {

    p.resize(n);
    x0.resize(n);
    g0.resize(n);
    dx0.resize(n);
    dg0.resize(n);

    for(size_t i=0;i<n;i++) {
      p[i]=0.0;
      x0[i]=0.0;
      g0[i]=0.0;
      dx0[i]=0.0;
      dg0[i]=0.0;
    }
    
    st_dx.resize(n);
    st_grad.resize(n);
    
    wrap.av_x_alpha.resize(n);
    wrap.av_g_alpha.resize(n);
    wrap.dim=n;
    dim=n;

    return success;
  }
    
  /// Free the allocated memory
  virtual int free() {
    wrap.av_x_alpha.clear();
    wrap.av_g_alpha.clear();
    st_grad.clear();
    dg0.clear();
    dx0.clear();
    g0.clear();
    x0.clear();
    p.clear();
    st_dx.clear();
    wrap.dim=0;
    dim=0;
    return 0;
  }
      
  /// Reset the minimizer to use the current point as a new starting point
  int restart() {
    iter=0;
    return success;
  }

  /// Set the function and initial guess
  virtual int set(vec_t &x, double u_step_size, double tol_u, 
		  func_t &ufunc) {
      
    iter=0;
    step=u_step_size;
    delta_f=0;
      
    st_x=&x;
      
    st_f=ufunc(dim,x);
    agrad->set_function(ufunc);
    (*agrad)(dim,x,st_grad);
      
    /* Use the gradient as the initial direction */

    for(size_t i=0;i<dim;i++) {
      x0[i]=x[i];
      g0[i]=st_grad[i];
    }
    g0norm=o2scl_cblas::dnrm2(dim,g0);
    for(size_t i=0;i<dim;i++) {
      p[i]=st_grad[i];
    }
    o2scl_cblas::dscal(-1.0/g0norm,dim,p);

    // pnorm should be 1
    pnorm=o2scl_cblas::dnrm2(dim,p);
    fp0=-g0norm;

    /* Prepare the wrapper */
      
    wrap.prepare_wrapper(ufunc,0,x0,st_f,g0,p,agrad);
      
    /* Prepare 1d minimisation parameters */

    rho=0.01;
    sigma=tol_u;
    tau1=9;
    tau2=0.05;
    tau3=0.5;
    // use cubic interpolation where possible 
    order=3;  

    return success;
 
  }

  /// Set the function, the gradient, and the initial guess
  virtual int set_de(vec_t &x, double u_step_size, double tol_u, 
		     func_t &ufunc, dfunc_t &udfunc) {
      
    iter=0;
    step=u_step_size;
    delta_f=0;
      
    st_x=&x;
      
    st_f=ufunc(dim,x);
    udfunc(dim,x,st_grad);
      
    /* Use the gradient as the initial direction */

    for(size_t i=0;i<dim;i++) {
      x0[i]=x[i];
      g0[i]=st_grad[i];
    }
    g0norm=o2scl_cblas::dnrm2(dim,g0);
    for(size_t i=0;i<dim;i++) {
      p[i]=st_grad[i];
    }
    o2scl_cblas::dscal(-1.0/g0norm,dim,p);

    // pnorm should be 1
    pnorm=o2scl_cblas::dnrm2(dim,p);
    fp0=-g0norm;

    /* Prepare the wrapper */
      
    wrap.prepare_wrapper(ufunc,&udfunc,x0,st_f,g0,p,agrad);
      
    /* Prepare 1d minimisation parameters */

    rho=0.01;
    sigma=tol_u;
    tau1=9;
    tau2=0.05;
    tau3=0.5;
    // use cubic interpolation where possible 
    order=3;  

    return success;
 
  }

  /// The size of the first trial step (default 0.01)
  double step_size;
    
  /// The tolerance for the 1-dimensional minimizer
  double lmin_tol;

  /// Default automatic gradient object
  def_auto_grad_t def_grad;

  /** \brief Calculate the minimum \c min of \c func w.r.t the
      array \c x of size \c nn.
  */
  virtual int mmin(size_t nn, vec_t &xx, double &fmin, 
		   func_t &ufunc) {

    if (nn==0) {
      O2SCL_ERR2("Tried to min over zero variables ",
		     " in mmin_bfgs2::mmin().",exc_einval);
    }

    int xiter=0, status;

    allocate(nn);

    set(xx,step_size,lmin_tol,ufunc);

    do {

      xiter++;

      status=iterate();
	
      if (status) {
	break;
      }

      // Equivalent to gsl_multimin_test_gradient with
      // additional code to print out present iteration
	
      double norm=o2scl_cblas::dnrm2(nn,st_grad);
	
      if(this->verbose>0) {
	this->print_iter(nn,*st_x,st_f,xiter,
			 norm,this->tol_rel,"mmin_bfgs2");
      }

      if (norm<this->tol_rel) {
	status=success;
      } else {
	status=gsl_continue;
      }
	
    } while (status == gsl_continue && xiter < this->ntrial);
      
    for(size_t i=0;i<nn;i++) xx[i]=(*st_x)[i];
    fmin=st_f;
      
    free();
    this->last_ntrial=xiter;
      
    if (status==gsl_continue && xiter==this->ntrial) {
      O2SCL_CONV_RET("Too many iterations in mmin_bfgs2::mmin().",
		     exc_emaxiter,this->err_nonconv);
    }
    return status;
  }

  /** \brief Calculate the minimum \c min of \c func w.r.t the
      array \c x of size \c nn.
  */
  virtual int mmin_de(size_t nn, vec_t &xx, double &fmin, 
		      func_t &ufunc, dfunc_t &udfunc) {

    if (nn==0) {
      O2SCL_ERR2("Tried to min over zero variables ",
		     "in mmin_bfgs2::mmin().",exc_einval);
    }

    int xiter=0, status;

    allocate(nn);

    set_de(xx,step_size,lmin_tol,ufunc,udfunc);

    do {
      xiter++;

      status=iterate();
	
      if (status) {
	break;
      }

      // Equivalent to gsl_multimin_test_gradient with
      // additional code to print out present iteration
	
      double norm=o2scl_cblas::dnrm2(nn,st_grad);
	
      if(this->verbose>0) {
	this->print_iter(nn,*st_x,st_f,xiter,
			 norm,this->tol_rel,"mmin_bfgs2");
      }

      if (norm<this->tol_rel) status=success;
      else status=gsl_continue;

    } while (status == gsl_continue && xiter < this->ntrial);
      
    for(size_t i=0;i<nn;i++) xx[i]=(*st_x)[i];
    fmin=st_f;
      
    free();
    this->last_ntrial=xiter;

    if (status==gsl_continue && xiter==this->ntrial) {
      O2SCL_CONV_RET("Too many iterations in mmin_bfgs2::mmin_de().",
		     exc_emaxiter,this->err_nonconv);
    }

    return status;
  }

#ifndef DOXYGEN_INTERNAL

 private:

  mmin_bfgs2<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t>
  (const mmin_bfgs2<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t> &);
  mmin_bfgs2<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t>& operator=
  (const mmin_bfgs2<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t>&);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
