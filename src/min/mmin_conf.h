/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2022, Andrew W. Steiner

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
#ifndef O2SCL_GSL_MMIN_CONF_H
#define O2SCL_GSL_MMIN_CONF_H

/** \file mmin_conf.h
    \brief File defining \ref o2scl::mmin_conf
*/

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/mmin.h>
#include <o2scl/misc.h>
#include <o2scl/cblas.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Base minimization routines for mmin_conf and mmin_conp

      This class is used by the mmin_conf and mmin_conp
      minimizers to perform the line minimization along a specified
      direction. It is not intended for a casual end-user.
      
      Default template arguments
      - \c func_t - \ref multi_funct
      - \c vec_t - \ref boost::numeric::ublas::vector \<double \>
      - \c dfunc_t - \ref mm_funct
      - \c auto_grad_t - \ref gradient \< func_t \>
      - \c def_auto_grad_t - \ref gradient_gsl \<func_t \>
      
      \comment
      Doxygen has had trouble parsing this template, so I have
      simplified it slightly for the documentation below.
      \endcomment
  */
  template<class func_t = multi_funct, 
    class vec_t = boost::numeric::ublas::vector<double>, 
    class dfunc_t = grad_funct,
    class auto_grad_t = gradient<multi_funct,
    boost::numeric::ublas::vector<double> >,
    class def_auto_grad_t = 
    gradient_gsl<multi_funct,boost::numeric::ublas::vector<double> > > 
    class mmin_gsl_base : public mmin_base<func_t,dfunc_t,vec_t> {
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
  
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    /// User-specified function
    func_t *func;
  
    /// User-specified gradient
    dfunc_t *grad;
  
    /// Automatic gradient object
    auto_grad_t *agrad;
      
    /// If true, a gradient has been specified
    bool grad_given;

    /// Memory size
    size_t dim;

    /// Take a step
    void take_step(const vec_t &x, const vec_t &px,
		   double stepx, double lambda, vec_t &x1x, 
		   vec_t &dx) {

      for(size_t i=0;i<this->dim;i++) dx[i]=0.0;
      o2scl_cblas::daxpy(-stepx*lambda,this->dim,px,dx);
    
      for(size_t i=0;i<this->dim;i++) x1x[i]=x[i];
      o2scl_cblas::daxpy(1.0,this->dim,dx,x1x);
	
      return;
    }

    /** \brief Line minimization
	  
	Do a line minimisation in the region (xa,fa) (xc,fc) to find
	an intermediate (xb,fb) satisifying fa > fb < fc.  Choose an
	initial xb based on parabolic interpolation.
    */
    void intermediate_point(const vec_t &x, const vec_t &px,
			    double lambda, double pg, double stepa, 
			    double stepc, double fa, double fc,
			    vec_t &x1x, vec_t &dx, vec_t &gradient, 
			    double *stepx, double *f) {
      double stepb, fb;
    
      bool trial_failed;

      // This looks like an infinite loop, but will eventually 
      // stop because stepb will eventually become zero. This
      // error condition is handled later in iterate().

      do {

	double u = fabs (pg * lambda * stepc);
	stepb = 0.5 * stepc * u / ((fc - fa) + u);
	
	take_step(x, px, stepb, lambda, x1x, dx);

	// Immediately exit if trial point is the same as the
	// initial point (as updated in GSL 1.15)
	bool vector_equal=true;
	for(size_t i=0;i<dim;i++) {
	  if (x[i]!=x1x[i]) vector_equal=false;
	}
	if (vector_equal) {

	  *stepx=0.0;
	  *f=fa;

	  // Evaluate the gradient
	  if (grad_given) {
	    (*grad)(dim,x1x,gradient);
	  } else {
	    (*agrad)(dim,x1x,gradient);
	  }
	  
	  return;
	}

	// Evaluate the function
	fb=(*func)(dim,x1x);
	  
	trial_failed=false;
	if (fb >= fa  && stepb > 0.0) {
	  // [GSL] downhill step failed, reduce step-size and try again 
	  fc = fb;
	  stepc = stepb;
	  trial_failed=true;
	}

      } while (trial_failed);

      *stepx = stepb;
      *f = fb;

      // Evaluate the gradient
      if (grad_given) {
	(*grad)(dim,x1x,gradient);
      } else {
	(*agrad)(dim,x1x,gradient);
      }

      return;
    }
      
    /** \brief Perform the minimization
	
	Starting at (x0, f0) move along the direction p to find a
	minimum f(x0 - lambda * p), returning the new point x1 =
	x0-lambda*p, f1=f(x1) and g1 = grad(f) at x1.
    */
    void min(const vec_t &x, const vec_t &xp, double lambda,
	     double stepa, double stepb, double stepc, double fa, 
	     double fb, double fc, double xtol, vec_t &x1x, 
	     vec_t &dx1x, vec_t &x2x, vec_t &dx2x, vec_t &gradient, 
	     double *xstep, double *f, double *gnorm_u) {
    
      double u = stepb;
      double v = stepa;
      double w = stepc;
      double fu = fb;
      double fv = fa;
      double fw = fc;

      double old2 = fabs(w - v);
      double old1 = fabs(v - u);

      double stepm, fm, pg, gnorm1;

      int xiter = 0;

      for(size_t i=0;i<dim;i++) {
	x2x[i]=x1x[i];
	dx2x[i]=dx1x[i];
      }

      *f = fb;
      *xstep = stepb;
      *gnorm_u = o2scl_cblas::dnrm2(dim,gradient);

      while (true) {

	xiter++;

	if (xiter > nmaxiter) {
	  // Exceeded maximum number of iterations
	  return;  
	}
	
	{
	  double dw = w - u;
	  double dv = v - u;
	  double du = 0.0;
	  double e1 = ((fv - fu) * dw * dw + (fu - fw) * dv * dv);
	  double e2 = 2.0 * ((fv - fu) * dw + (fu - fw) * dv);

	  if (e2 != 0.0) {
	    du = e1 / e2;
	  }

	  if (du > 0.0 && du < (stepc - stepb) && fabs(du) < 0.5 * old2) {
	    stepm = u + du;
	  } else if (du < 0.0 && du > (stepa - stepb) && 
		     fabs(du) < 0.5 * old2) {
	    stepm = u + du;
	  } else if ((stepc - stepb) > (stepb - stepa)) {
	    stepm = 0.38 * (stepc - stepb) + stepb;
	  } else {
	    stepm = stepb - 0.38 * (stepb - stepa);
	  }
	}

	take_step (x, xp, stepm, lambda, x1x, dx1x);

	fm=(*func)(dim,x1x);
	
	if (fm > fb) {

	  if (fm < fv) {
	    w = v;
	    v = stepm;
	    fw = fv;
	    fv = fm;
	  } else if (fm < fw) {
	    w = stepm;
	    fw = fm;
	  }
	  if (stepm < stepb) {
	    stepa = stepm;
	    fa = fm;
	  } else {
	    stepc = stepm;
	    fc = fm;
	  }

	} else if (fm <= fb) {

	  old2 = old1;
	  old1 = fabs(u - stepm);
	  w = v;
	  v = u;
	  u = stepm;
	  fw = fv;
	  fv = fu;
	  fu = fm;

	  for(size_t i=0;i<dim;i++) {
	    x2x[i]=x1x[i];
	    dx2x[i]=dx1x[i];
	  }

	  if (grad_given) {
	    (*grad)(dim,x1x,gradient);
	  } else {
	    (*agrad)(dim,x1x,gradient);
	  }

	  pg=o2scl_cblas::ddot(dim,xp,gradient);
	  gnorm1 = o2scl_cblas::dnrm2(dim,gradient);

	  *f = fm;
	  *xstep = stepm;
	  *gnorm_u = gnorm1;

	  if (fabs (pg * lambda / gnorm1) < xtol) {
	    // [GSL] SUCCESS 
	    return; 
	  }
	    
	  if (stepm < stepb) {
	    stepc = stepb;
	    fc = fb;
	    stepb = stepm;
	    fb = fm;
	  } else {
	    stepa = stepb;
	    fa = fb;
	    stepb = stepm;
	    fb = fm;
	  }
      
	}

      }
    
      return;
    }

#endif

    public:

    /// Stepsize for finite-differencing ( default \f$ 10^{-4} \f$ )
    double deriv_h;

    /// Maximum iterations for line minimization (default 10)
    int nmaxiter;
      
    /// Default automatic \gradient object
    def_auto_grad_t def_grad;
      
    mmin_gsl_base() {
      deriv_h=1.0e-4;
      nmaxiter=10;
      grad_given=false;
      agrad=&def_grad;
    }
      
    /// Set the function
    int base_set(func_t &ufunc, auto_grad_t &u_def_grad) {
      func=&ufunc;
      agrad=&u_def_grad;
      grad_given=false;
      return 0;
    }

    /// Set the function and the \gradient
    int base_set_de(func_t &ufunc, dfunc_t &udfunc) {
      func=&ufunc;
      grad=&udfunc;
      grad_given=true;
      return 0;
    }

    /// Allocate memory
    int base_allocate(size_t nn) {
      dim=nn;
      return 0;
    }

    /// Clear allocated memory
    int base_free() {
      dim=0;
      return 0;
    }

    
    };
  
  /** \brief Multidimensional minimization by the Fletcher-Reeves 
      conjugate \gradient algorithm (GSL)

      This class performs multidimensional minimization by the
      Fletcher-Reeves conjugate \gradient algorithm (GSL). The
      functions mmin() and mmin_de() \minimize a given function until
      the \gradient is smaller than the value of mmin::tol_rel
      (which defaults to \f$ 10^{-4} \f$ ).

      This class has a high-level interface using mmin() or mmin_de()
      which automatically performs the memory allocation and
      minimization, or a GSL-like interface using allocate(), free(),
      interate() and set() or set_simplex().

      \verbatim embed:rst
      See an example for the usage of this class in 
      :ref:`Multidimensional minimizer example`.
      \endverbatim
      
      Default template arguments
      - \c func_t - \ref multi_funct
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \>
      - \c dfunc_t - \ref grad_funct
      - \c auto_grad_t - \ref gradient \< func_t \>
      - \c def_auto_grad_t - \ref gradient_gsl \< func_t \>

      Note that the state variable \c max_iter has not been included
      here, because it was not really used in the original GSL code
      for these minimizers.
  */
  template<class func_t = multi_funct, 
    class vec_t = boost::numeric::ublas::vector<double>, 
    class dfunc_t = grad_funct,
    class auto_grad_t = gradient<multi_funct,
    boost::numeric::ublas::vector<double> >,
    class def_auto_grad_t = 
    gradient_gsl<multi_funct,boost::numeric::ublas::vector<double> > > 
    class mmin_conf : 
    public mmin_gsl_base<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t> {
      
#ifndef DOXYGEN_INTERNAL
      
    protected:
    
    /// \name The original variables from the GSL state structure
    //@{
    /// Iteration number
    int iter;
    /// Stepsize
    double step;
    /// Tolerance
    double tol;
    /// Desc
    vec_t x1;
    /// Desc
    vec_t dx1;
    /// Desc
    vec_t x2;
    /// Desc
    double pnorm;
    /// Desc
    vec_t p;
    /// Desc
    double g0norm;
    /// Desc
    vec_t g0;
    //@}

    /// \name Store the arguments to set() so we can use them for iterate()
    //@{
    /// Proposed minimum
    vec_t ugx;
    /// Gradient
    vec_t ugg;
    /// Proposed step
    vec_t udx;
    /// Desc
    double it_min;
    //@}

#endif
    
    public:
      
    mmin_conf() {
      lmin_tol=1.0e-4;
      this->tol_rel=1.0e-4;
      step_size=0.01;
    }
    
    virtual ~mmin_conf() {}

    /// Tolerance for the line minimization (default \f$ 10^{-4} \f$)
    double lmin_tol;

    /// Size of the initial step (default 0.01)
    double step_size;

    /// \name GSL-like lower level interface 
    //@{
    /// Perform an iteration
    virtual int iterate() {
	
      if (this->dim==0) {
	O2SCL_ERR2("Memory not allocated in ",
		       "mmin_conf::iterate().",exc_efailed);
      }
	
      vec_t &x=ugx;
      vec_t &gradient=ugg;
      vec_t &dx=udx;

      double fa = it_min, fb, fc;
      double dir;
      double stepa = 0.0, stepb, stepc = step;

      double g1norm;
      double pg;

      if (pnorm == 0.0 || g0norm == 0.0) {
	for(size_t i=0;i<this->dim;i++) dx[i]=0.0;
	O2SCL_CONV2_RET("Either pnorm or g0norm vanished in ",
			"in mmin_conf::iterate().",
			exc_enoprog,this->err_nonconv);
      }
    
      /* [GSL] Determine which direction is downhill, +p or -p */

      pg=o2scl_cblas::ddot(this->dim,p,gradient);

      dir = (pg >= 0.0) ? +1.0 : -1.0;

      /* [GSL] Compute new trial point at x_c= x - step * p, where p is the
	 current direction */

      this->take_step (x, p, stepc, dir / pnorm, x1, dx);

      /* [GSL] Evaluate function and gradient at new point xc */
	
      fc=(*this->func)(this->dim,x1);

      if (fc < fa) {

	/* [GSL] Success, reduced the function value */
	step = stepc * 2.0;
	it_min = fc;
	for(size_t i=0;i<this->dim;i++) {
	  x[i]=x1[i];
	}
	  
	if (this->grad_given) {
	  (*this->grad)(this->dim,x1,gradient);
	} else {
	  (*this->agrad)(this->dim,x1,gradient);
	}

	return success;
      }

      /* [GSL] Do a line minimization in the region (xa,fa) (xc,fc) to find an
	 intermediate (xb,fb) satisifying fa > fb < fc.  Choose an initial
	 xb based on parabolic interpolation 
      */
      this->intermediate_point(x,p,dir / pnorm,pg,stepa,stepc,fa,fc,x1,
			       dx1,gradient,&stepb,&fb);
	
      if (stepb == 0.0) {
	O2SCL_CONV2_RET("Variable stepb vanished in ",
			"in mmin_conf::iterate().",
			exc_enoprog,this->err_nonconv);
      }
	
      this->min(x,p,dir / pnorm,stepa,stepb,stepc,fa,fb,fc,tol,
		x1,dx1,x2,dx,gradient,&step,&it_min,&g1norm);
    
      for(size_t i=0;i<this->dim;i++) x[i]=x2[i];
    
      /* [GSL] Choose a new conjugate direction for the next step */
      iter=(iter+1) % this->dim;

      if (iter==0) {

	for(size_t i=0;i<this->dim;i++) p[i]=gradient[i];
	pnorm=g1norm;

      } else {

	/* [GSL] p' = g1 - beta * p */
	double beta=-pow(g1norm/g0norm, 2.0);
	o2scl_cblas::dscal(-beta,this->dim,p);
	o2scl_cblas::daxpy(1.0,this->dim,gradient,p);
	pnorm=o2scl_cblas::dnrm2(this->dim,p);
      }
	
      g0norm=g1norm;
      for(size_t i=0;i<this->dim;i++) {
	g0[i]=gradient[i];
      }

      return success;
    }
    
    /** \brief Allocate the memory
     */
    virtual int allocate(size_t n) {
    
      x1.resize(n);
      dx1.resize(n);
      x2.resize(n);
      p.resize(n);
      g0.resize(n);
      udx.resize(n);
      ugg.resize(n);
      ugx.resize(n);

      this->base_allocate(n);

      return success;

    }
    
    /// Free the allocated memory
    virtual int free() {
    
      this->base_free();

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

      if (this->dim==0) {
	O2SCL_ERR2("Memory not allocated in ",
		       "mmin_conf::set().",exc_efailed);
      }
    
      this->base_set(ufunc,*this->agrad);
      for(size_t i=0;i<this->dim;i++) ugx[i]=x[i];
	
      iter=0;
      step=u_step_size;
      tol=tol_u;
    
      /// Evaluate the function and its gradient
      it_min=ufunc(this->dim,x);
      this->agrad->set_function(ufunc);
      (*this->agrad)(this->dim,x,ugg);
	
      /* [GSL] Use the gradient as the initial direction */
	
      for(size_t i=0;i<this->dim;i++) {
	p[i]=ugg[i];
	g0[i]=ugg[i];
      }
      
      double gnorm=o2scl_cblas::dnrm2(this->dim,ugg);
	
      pnorm=gnorm;
      g0norm=gnorm;
	
      return success;
    }

    /// Set the function and initial guess
    virtual int set_de(vec_t &x, double u_step_size, double tol_u,
		       func_t &ufunc, dfunc_t &udfunc) {

      if (this->dim==0) {
	O2SCL_ERR2("Memory not allocated in ",
		       "mmin_conf::set_de().",exc_efailed);
      }
	
      this->base_set_de(ufunc,udfunc);
      for(size_t i=0;i<this->dim;i++) ugx[i]=x[i];
	
      iter=0;
      step=u_step_size;
      tol=tol_u;

      // Evaluate the function and its gradient
      it_min=ufunc(this->dim,x);
      udfunc(this->dim,x,ugg);
	
      /* Use the gradient as the initial direction */
	
      for(size_t i=0;i<this->dim;i++) {
	p[i]=ugg[i];
	g0[i]=ugg[i];
      }
	
      double gnorm=o2scl_cblas::dnrm2(this->dim,ugg);
	
      pnorm=gnorm;
      g0norm=gnorm;

      return success;
    }
    //@}

    /// \name Basic usage
    //@{
    /** \brief Calculate the minimum \c min of \c func w.r.t the
	array \c x of size \c nvar.
    */
    virtual int mmin(size_t nn, vec_t &xx, double &fmin,
		     func_t &ufunc) {

      if (nn==0) {
	O2SCL_ERR2("Tried to min over zero variables ",
		       "in mmin_conf::mmin().",exc_einval);
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

	double norm=o2scl_cblas::dnrm2(nn,ugg);
	  
	if(this->verbose>0) {
	  this->print_iter(nn,ugx,it_min,xiter,
			   norm,this->tol_rel,type());
	}
    
	if (norm<this->tol_rel) status=success;
	else status=gsl_continue;

      } while (status==gsl_continue && xiter < this->ntrial);

      for(size_t i=0;i<nn;i++) xx[i]=ugx[i];
      fmin=it_min;
	
      free();
      this->last_ntrial=xiter;
	
      if (status==gsl_continue && xiter==this->ntrial) {
	std::string err="Exceeded max number of iterations, "+
	  dtos(this->ntrial)+", in mmin_conf::mmin().";
	O2SCL_CONV_RET(err.c_str(),exc_emaxiter,this->err_nonconv);
      }

      return status;
    }

    /** \brief Calculate the minimum \c min of \c func w.r.t the
	array \c x of size \c nvar.
    */
    virtual int mmin_de(size_t nn, vec_t &xx, double &fmin, 
			func_t &ufunc, dfunc_t &udfunc) {

      if (nn==0) {
	O2SCL_ERR2("Tried to min over zero variables ",
		       "in mmin_conf::mmin_de().",exc_einval);
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

	double norm=o2scl_cblas::dnrm2(nn,ugg);
	  
	if(this->verbose>0) {
	  this->print_iter(nn,ugx,it_min,xiter,
			   norm,this->tol_rel,type());
	}
    
	if (norm<this->tol_rel) status=success;
	else status=gsl_continue;

      }
      while (status==gsl_continue && xiter < this->ntrial);
     
      for(size_t i=0;i<nn;i++) xx[i]=ugx[i];
      fmin=it_min;
	
      free();
      this->last_ntrial=xiter;

      if (status==gsl_continue && xiter==this->ntrial) {
	std::string err="Exceeded max number of iterations, "+
	  dtos(this->ntrial)+", in mmin_conf::mmin().";
	O2SCL_CONV_RET(err.c_str(),exc_emaxiter,this->err_nonconv);
      }

      return status;
    }
    //@}

    /// Return string denoting type("mmin_conf")
    virtual const char *type() { return "mmin_conf";}

#ifndef DOXYGEN_INTERNAL

 private:

  mmin_conf<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t>
  (const mmin_conf<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t> &);
  mmin_conf<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t>& operator=
  (const mmin_conf<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t>&);

#endif

    };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
