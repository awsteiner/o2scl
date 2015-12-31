/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2016, Andrew W. Steiner

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
#ifndef O2SCL_MMIN_CONP_H
#define O2SCL_MMIN_CONP_H

/** \file mmin_conp.h
    \brief File defining \ref o2scl::mmin_conp
*/
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <o2scl/mmin_conf.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multidimensional minimization by the Polak-Ribiere
      conjugate gradient algorithm (GSL)

      The functions mmin() and mmin_de() min a given function
      until the gradient is smaller than the value of mmin::tol_rel
      (which defaults to \f$ 10^{-4} \f$ ).

      See an example for the usage of this class in 
      \ref ex_mmin_sect .
  */
  template<class func_t = multi_funct11, 
    class vec_t = boost::numeric::ublas::vector<double>, 
    class dfunc_t = grad_funct11,
    class auto_grad_t = gradient<multi_funct11,
    boost::numeric::ublas::vector<double> >,
    class def_auto_grad_t = 
    gradient_gsl<multi_funct11,boost::numeric::ublas::vector<double> > > 
    class mmin_conp : 
    public mmin_conf<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t> {
    
    public:
 
    mmin_conp() : mmin_conf<func_t,vec_t,dfunc_t,auto_grad_t,
    def_auto_grad_t>() {
      }
 
    /// Perform an iteration
    virtual int iterate() {
      
      if (this->dim==0) {
	O2SCL_ERR("Memory not allocated in iterate().",exc_efailed);
      }
      
      vec_t &x=this->ugx;
      vec_t &gradient=this->ugg;
      vec_t &dx=this->udx;
      
      double fa = this->it_min, fb, fc;
      double dir;
      double stepa = 0.0, stepb, stepc=this->step;
      
      double g1norm;
      double pg;
      
      if (this->pnorm == 0.0 || this->g0norm == 0.0) {
	for(size_t i=0;i<this->dim;i++) dx[i]=0.0;
	O2SCL_CONV2_RET("Either pnorm or g0norm vanished in ",
			"in mmin_conp::iterate().",
			exc_enoprog,this->err_nonconv);
	return 0;
      }
      
      /* Determine which direction is downhill, +p or -p */

      pg=o2scl_cblas::ddot(this->dim,this->p,gradient);
      
      dir = (pg >= 0.0) ? +1.0 : -1.0;
      
      /* Compute new trial point at x_c= x - step * p, where p is the
	 current direction */
      
      this->take_step(x, this->p, stepc, dir / this->pnorm, this->x1, dx);
      
      /* Evaluate function and gradient at new point xc */
      
      fc=(*this->func)(this->dim,this->x1);
      
      if (fc < fa) {

	/* Success, reduced the function value */
	this->step = stepc * 2.0;
	this->it_min = fc;
	for(size_t i=0;i<this->dim;i++) {
	  x[i]=this->x1[i];
	}

	if (this->grad_given) {
	  (*this->grad)(this->dim,this->x1,gradient);
	} else {
	  (*this->agrad)(this->dim,this->x1,gradient);
	}

	return success;
      }
      
      /* Do a line minimisation in the region (xa,fa) (xc,fc) to find an
	 intermediate (xb,fb) satisifying fa > fb < fc.  Choose an initial
	 xb based on parabolic interpolation */
      
      this->intermediate_point(x, this->p, dir / this->pnorm, pg,
			       stepa, stepc, fa, fc, this->x1, this->dx1,
			       gradient, &stepb, &fb);
      
      if (stepb == 0.0) {
	O2SCL_CONV2_RET("Variable stepb vanished in mmin_conp::",
			"iterate().",exc_enoprog,this->err_nonconv);
      }
      
      this->min(x,this->p,dir / this->pnorm,stepa,stepb,stepc, fa, 
		fb, fc, this->tol, this->x1, this->dx1, this->x2, 
		dx, gradient, &(this->step), &(this->it_min), &g1norm);
      
      for(size_t i=0;i<this->dim;i++) x[i]=this->x2[i];
      
      /* Choose a new conjugate direction for the next step */
      
      this->iter = (this->iter + 1) % this->dim;
      
      if (this->iter == 0) {
	for(size_t i=0;i<this->dim;i++) this->p[i]=gradient[i];
	this->pnorm = g1norm;
      } else {
	/* p' = g1 - beta * p */
	double g0g1, beta;
	
	/* g0' = g0 - g1 */
	o2scl_cblas::daxpy(-1.0,this->dim,gradient,this->g0);
	/* g1g0 = (g0-g1).g1 */
	g0g1=o2scl_cblas::ddot(this->dim,this->g0,gradient);
	/* beta = -((g1 - g0).g1)/(g0.g0) */
	beta = g0g1 / (this->g0norm*this->g0norm);       
	
	o2scl_cblas::dscal(-beta,this->dim,this->p);
	o2scl_cblas::daxpy(1.0,this->dim,gradient,this->p);
	this->pnorm = o2scl_cblas::dnrm2(this->dim,this->p);
      }
      
      this->g0norm = g1norm;
      for(size_t i=0;i<this->dim;i++) {
	this->g0[i]=gradient[i];
      }
      
      return success;
    }
    
    /// Return string denoting type("mmin_conp")
    virtual const char *type() { return "mmin_conp";}

#ifndef DOXYGEN_INTERNAL

    private:

    mmin_conp<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t>
    (const mmin_conp<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t> &);
    mmin_conp<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t>& operator=
    (const mmin_conp<func_t,vec_t,dfunc_t,auto_grad_t,def_auto_grad_t>&);

#endif
    
    };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
