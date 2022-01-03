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
/*----------------------------------------------------------------------------*
 * Open Optimization Library - Constrained Minimization
 * 
 * gencan/gencan.c
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * Sergio Drumond Ventura
 * Luis Alberto D'Afonseca
 * Ricardo Biloti
 * since: Feb 16th, 2004
 *
 * $Id: gencan.c,v 1.16 2005/05/17 19:08:18 biloti Exp $
 *----------------------------------------------------------------------------*/
#ifndef O2SCL_OOL_MMIN_GENCAN_H
#define O2SCL_OOL_MMIN_GENCAN_H

/** \file mmin_constr_gencan.h
    \brief File defining \ref o2scl::mmin_constr_gencan 
*/

#include <o2scl/text_file.h>
#include <o2scl/multi_funct.h>
#include <o2scl/ool_constr_min.h>
#include <gsl/gsl_math.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Constrained minimization by the "GENCAN" method (OOL)
      
      \note Not working yet
  */
  template<class param_t, class func_t, class dfunc_t=func_t, 
    class hfunc_t=func_t, class vec_t=boost::numeric::ublas::vector<double> >
    class mmin_constr_gencan : 
    public ool_constr_min<param_t,func_t,dfunc_t,hfunc_t,vec_t> {
    
#ifndef DOXYGEN_INTERNAL
    
    protected:

    /// Desc (default 1.0)
    double cg_src;

    /// Temporary vector
    vec_t S;
    /// Temporary vector
    vec_t Y;
    /// Temporary vector
    vec_t D;
    /// Temporary vector
    vec_t cg_W;
    /// Temporary vector
    vec_t cg_R;
    /// Temporary vector
    vec_t cg_D;
    /// Temporary vector
    vec_t cg_Sprev;
    /// Temporary vector
    vec_t Xtrial;
    /// Temporary vector
    vec_t tnls_Xtemp;
    /// Temporary vector
    vec_t near_l;
    /// Temporary vector
    vec_t near_u;
    
    /// Desc
    int *Ind;

#ifdef NEVER_DEFINED

    /// Desc
    int spgls() {

      gsl_vector *X        = M->x;
      gsl_vector *gradient = M->gradient;

      size_t nn = X->size;

      /* Direct access to vector data */
      double *l = st->L->data;
      double *u = st->U->data;
      double *d = st->D->data;
      double *x = X->data;

      double *xtrial = st->Xtrial->data;

      /* Internal variables */
      size_t interp;
      size_t imax;

      double alpha;
      double dinfn;
      double gtd;
      double ftrial;

      /* Compute first trial point, spectral projected gradient direction,
       * and directional derivative <g,d> */
      alpha = 1;

      /* Xtrial = min{ U, max[ L, ( X-lambda G ) ] } */
      gsl_vector_memcpy( st->Xtrial, X );
      gsl_blas_daxpy( -(st->lambda), gradient, st->Xtrial );
      conmin_vector_minofmax( st->n, xtrial, u, l, xtrial );

      /* D = Xtrial - X */
      gsl_vector_memcpy( st->D, st->Xtrial );
      gsl_vector_sub( st->D, X );

      /* Inifite norm of D and < G, D > */
      imax  = gsl_blas_idamax( st->D );
      dinfn = fabs( gsl_vector_get( st->D, imax ) );
      gsl_blas_ddot( gradient, st->D, &gtd );

      /* Evaluate objective function */
      OOL_CONMIN_EVAL_F( M, st->Xtrial, ftrial );

      interp = 0;

      /* While Armijo isn't satisfied repeat */
      while (ftrial > M->f + P->gamma*alpha*gtd ) {

	/* Test if the trial point has a function value lower than fmin */
	if (ftrial < M->f ) {

	  M->f = ftrial;
	  gsl_vector_memcpy( X, st->Xtrial );

	  return OOL_UNBOUNDEDF;
	}

	interp++;

	if (alpha < P->sigma1 ) {
	  alpha /= P->nint;
	} else {
	  /* quadratic model */
	  double atemp = ( -gtd*alpha*alpha )/
	    ( 2.0*(ftrial-M->f-alpha*gtd) );

	  if (atemp < P->sigma1 || atemp > P->sigma2*alpha ) {
	    alpha /= P->nint;
	  } else {
	    alpha  = atemp;
	  }
	}

	/* Compute new trial point
	 * Xtrial = X + alpha D */
	gsl_vector_memcpy( st->Xtrial, X );
	gsl_blas_daxpy( alpha, st->D, st->Xtrial );

	/* Evaluate objective function */
	OOL_CONMIN_EVAL_F( M, st->Xtrial, ftrial );

	/* Test whether at least mininterp interpolations were made
	 * and two consecutive iterates are close enough */
	if( interp > P->mininterp &&
	    are_close( nn, alpha, d, x, P->epsrel, P->epsabs )) {

	  M->f = ftrial;
	  gsl_vector_memcpy( X, st->Xtrial );
	    
	  return OOL_FLSEARCH;
	}
      }
      
      /* Set new values of f and X */
      M->f = ftrial;
      gsl_vector_memcpy( X, st->Xtrial );

      return OOL_SUCCESS;
    }

    /// Desc
    int tnls() {

      gsl_vector *X        = M->x;
      gsl_vector *gradient = M->gradient;
      gsl_vector *Xplus    = st->Xtrial;

      /* Direct access to vector data */
      double *x = X->data;
      double *g = gradient->data;
      double *d = st->D->data;
      double *xplus = Xplus->data;

      /* Constant values */
      const size_t nind = st->nind;

      /* Internal variables */
      double fplus;
      double gtd;
      double alpha;
      double gptd;

      /* Compute directional derivative */
      gtd = cblas_ddot( (int)nind, g, 1, d, 1 );

      /* Compute first trial */
      alpha = GSL_MIN( 1.0, st->tnls_amax );

      /* Xplus = X + alpha D */
      conmin_vector_memcpy( nind, xplus, x );
      cblas_daxpy( alpha,  (int)nind,d, 1, xplus, 1 );

      /* Evaluate objective function */
      fplus = conmin_calc_f( M, nind, st->Ind, Xplus, X );

      /* Test Armijo and beta-condition and decide for:
       * 1 - accepting the trial point,
       * 2 - interpolate or
       * 3 - extrapolate. */
      if( st->tnls_amax > 1.0 ) {

	/* X+D belongs to the interior of the feasible set (amax > 1) */

	/* Verify Armijo */
	if( fplus <= M->f + P->gamma * alpha * gtd ) {

	  /* Armijo condition holds */

	  /* Evaluate the gradient of objective function */
	  conmin_calc_g( M, nind, st->Ind, Xplus, X, gradient );
	  /* Eval gptd = < g, d > */
	  gptd = cblas_ddot( (int)nind, g, 1, d, 1 );

	  /* Verify directional derivative (beta condition) */
	  if ( gptd >= P->beta * gtd ) {

	    /* Step = 1 was ok, finish the line search */

	    M->f = fplus;
	    conmin_vector_memcpy( nind, x, xplus );

	    return OOL_SUCCESS;
	  } else {
	    return tnls_extrapolation( M, st, P, alpha, fplus );
	  }
	} else {
	  return tnls_interpolation(M, st, P, alpha, fplus, gtd);
	}
      } else {
	/* x + d does not belong to the feasible set (amax <= 1) */
	if( fplus < M->f ) {
	  return tnls_extrapolation( M, st, P, alpha, fplus );
	} else {
	  return tnls_interpolation(M, st, P, alpha, fplus, gtd);
	}
      }
    }

    /// Desc
    int tnls_extrapolation(double alpha, double fplus) {

      gsl_vector *X        = M->x;
      gsl_vector *gradient = M->gradient;
      gsl_vector *Xplus    = st->Xtrial;

      /* Direct access to vector data */
      double *x = X->data;
      double *d = st->D->data;
      double *l = st->L->data;
      double *u = st->U->data;

      double *xplus = Xplus->data;
      double *xtemp = st->tnls_Xtemp->data;

      /* Constant values */
      const size_t nind = st->nind;

      /* Internal variables */
      double atemp;
      double ftemp;

      size_t ii, extrap;
      short same;

      /* Iterations */
      extrap = 0;
      do {

	extrap++;

	/* Test if maximum number of extrapolation was exceeded */
	if ( extrap > P->maxextrap ) {

	  M->f = fplus;
	  conmin_vector_memcpy( nind, x, xplus );

	  if (extrap > 0 || st->tnls_amax < 1){
	    conmin_calc_g( M, nind, st->Ind, Xplus, X, gradient );
	  }
	  return TNLS_MAXEXTRAP;
	}

	/* Chose new step */
	if (alpha < st->tnls_amax && st->tnls_amax < P->next*alpha ) {
	  atemp = st->tnls_amax;
	} else {
	  atemp = P->next * alpha;
	}

	/* Compute new trial point. Xtemp = X + atemp*D */
	conmin_vector_memcpy( nind, xtemp, x );
	cblas_daxpy( atemp, (int)nind, d, 1, xtemp, 1 );

	/* Project */
	if (atemp > st->tnls_amax ) {
	  conmin_vector_maxofmin( nind, xtemp, l, xtemp, u );
	}

	/* Test if this is not the same point as the previous one.
	 * This test is performed only when alpha >= amax. */
	if( alpha > st->tnls_amax ) {

	  same = 1;
	  for (ii = 0; ii<nind && same; ii++) {

	    double aux;

	    aux = P->epsrel * fabs( xplus[ii] );
	    
	    if ( fabs(xtemp[ii]-xplus[ii]) > GSL_MAX(aux,P->epsabs)) {
	      same = 0;
	    }
	  }

	  if (same) {

	    /* Finish the extrapolation with the current point */
	    M->f = fplus;

	    conmin_vector_memcpy( nind, x, xplus );

	    if (extrap > 0 || st->tnls_amax < 1){
	      conmin_calc_g( M, nind, st->Ind, Xplus, X, gradient );
	    }
	    return OOL_SUCCESS;
	  }
	}

	ftemp = conmin_calc_f( M, nind, st->Ind, st->tnls_Xtemp, X );

	if (ftemp < fplus) {
	  
	  /* If the functional value decreases then set the current
	   * point and continue the extrapolation */
	  
	  alpha = atemp;
	  fplus = ftemp;
	  conmin_vector_memcpy( nind, xplus, xtemp );

	  continue;

	} else {

	  /* If the functional value does not decrease then discard the
	   * last trial and finish the extrapolation with the previous
	   * point */

	  M->f = fplus;

	  conmin_vector_memcpy( nind, x, xplus );
	  if (extrap > 0 || st->tnls_amax < 1) {
	    conmin_calc_g( M, nind, st->Ind, X, X, gradient );
	  }

	  return OOL_SUCCESS;
	}

      } while (1);

      /* Just to make gcc happy */
      return OOL_SUCCESS;
    }

    /// Desc
    int tnls_interpolation(double alpha, double fplus, double gtd) {

      gsl_vector *X        = M->x;
      gsl_vector *gradient = M->gradient;
      gsl_vector *Xplus    = st->Xtrial;

      /* Direct access to vector data */
      double *x = X->data;
      double *d = st->D->data;
      double *xplus = Xplus->data;

      /* Constant values */
      const size_t nind = st->nind;

      /* Internal variables */
      size_t interp;
      double atemp;

      /* Initialization */
      interp = 0;

      /* Iterations */
      do {
	interp++;
	
	/* Test Armijo condition */
	if (fplus <= M->f + P->gamma * alpha * gtd ) {

	  /* Finish the line search */
	  M->f = fplus;

	  /* X = Xplus */
	  conmin_vector_memcpy( nind, x, xplus );

	  /* Evaluate objective function gradient */
	  conmin_calc_g( M, nind, st->Ind, X, X, gradient );

	  return OOL_SUCCESS;
	}
	/* Compute new step */
	if (alpha < P->sigma1 ) {
	  alpha= alpha / P->nint;
	} else {
	  /* quadratic model */
	  atemp = -gtd * alpha*alpha /
	    (2 * (fplus - M->f - alpha * gtd));

	  if( atemp < P->sigma1       ||
	      atemp > P->sigma2*alpha  ) {
	    alpha = alpha / P->nint;
	  } else {
	    alpha = atemp;
	  }
	}

	/* Compute new trial point: xplus = x + alpha d */
	conmin_vector_memcpy( nind, xplus, x );
	cblas_daxpy( alpha, (int)nind, d, 1, xplus, 1 );

	/* Evaluate objective function */
	fplus = conmin_calc_f( M, nind, st->Ind, Xplus, X );

	/* Test whether at least mininterp interpolations were made
	 * and the steplength is soo small */
	if ( are_close( nind, alpha, d, x,  P->epsrel, P->epsabs ) &&
	     interp > P->mininterp ){
	  return OOL_FLSEARCH;
	}
      }
      while( 1 );

      /* Just to make gcc happy */
      return OOL_SUCCESS;
    }

    /** Truncated Newton maximum step*/
    double tnls_maximum_step() {

      /* Direct access to vector data */
      double *x =  M->x->data;
      double *l = st->L->data;
      double *u = st->U->data;
      double *d = st->D->data;

      double step = P->infabs;
      size_t ii;
      
      for( ii = 0; ii < st->nind; ii++ ) {

	if( *d > 0 ) {
	      double aux = ( *u - *x ) / *d;
	      step = GSL_MIN( step, aux );
	    } else if( *d < 0 ) {
	      double aux = ( *l - *x ) / *d;
	      step = GSL_MIN( step, aux );
	    }

	  x++;
	  l++;
	  u++;
	  d++;
	}

      return step;
    }

    /** Spectral step length */
    void spg_steplength() {
      
      if (st->sty <= 0.0) {
	st->lambda = GSL_MAX( 1.0, st->xeucn ) / sqrt( st->gpeucn2 );
      } else {
	double aux;
	double ss   = st->sts / st->sty;
	
	aux = GSL_MAX( P->lspgmi, ss );
	st->lambda = GSL_MIN( P->lspgma, aux );
      }
    }

    /** Iterate */
    int actual_iterate() {
      
      /* Direct access to vector data */
      double *x =  M->x->data;
      double *l = st->L->data;
      double *u = st->U->data;
      /*  double *d = st->D->data; */

      /* Status of internal iterations */
      int lsflag;
      int cgflag;

      /* Internal variables */
      size_t ii, imax;

      /* Saving previous values */
      gsl_vector_memcpy( st->S, M->x        );
      gsl_vector_memcpy( st->Y, M->gradient );

      /* The actual iteration */
      if ( st->gieucn2 <= st->ometa2 * st->gpeucn2 ) {
	/* Compute the new iterate using an SPG iteration */

	/* Perform a line search with spectral continuous 
	   projected gradient */
	lsflag = spgls( M, st, P );

	/* Compute the gradient for the new iterate */
	OOL_CONMIN_EVAL_DF( M, M->x, M->gradient );

      } else {

	/* The new iterate will belong to the closure of the actual face */

	/* Shrink the point, its gradient and the bounds */
	conmin_shrink( st->nind, st->Ind, M->x        );
	conmin_shrink( st->nind, st->Ind, M->gradient );
	conmin_shrink( st->nind, st->Ind, st->L       );
	conmin_shrink( st->nind, st->Ind, st->U       );
    
	/* Compute the descent direction solving the newtonian system */
	cgflag = cg( M, st, P );

	/* Compute maximum step for truncated newton line search */
	if ( cgflag == CG_BOUNDARY ) {
	  st->tnls_amax = 1.0;
	} else {
	  st->tnls_amax = tnls_maximum_step( M, st, P );
	}

	/* Perform the line search */
	lsflag = tnls( M, st, P );

	/* Expand the point, its gradient and the bounds */
	conmin_expand( st->nind, st->Ind, M->x        );
	conmin_expand( st->nind, st->Ind, M->gradient );
	conmin_expand( st->nind, st->Ind, st->L       );
	conmin_expand( st->nind, st->Ind, st->U       );

	/* If the line search in the Truncated Newton direction
	   stopped due to a very small step discard this iteration
	   and force a SPG iteration */
	if ( lsflag == OOL_FLSEARCH ) {

	  /* Perform a line search with spectral projected gradient */
	  lsflag = spgls( M, st, P );

	  /* Compute the gradient for the new iterate */
	  OOL_CONMIN_EVAL_DF( M, M->x, M->gradient );
	}
      }

      /* Prepare for the next iteration */
      /* Adjustment */
      for( ii = 0; ii < st->n; ii++ ) {
	/* In principle, the current point could be slightly changed
	 * here, requiring a new function and gradient
	 * evaluation. But, according to the algorithms authors, this
	 * is done just to account for points that are "numerically"
	 * at faces already. Thus, no additional evaluations are
	 * performed. (May 11th, 2005).
	 */
	if     ( x[ii] <= st->near_l[ii] ) x[ii] = l[ii];
	else if( x[ii] >= st->near_u[ii] ) x[ii] = u[ii];
      }

      /* Compute infinite and Euclidian-norm of X */
      imax      = gsl_blas_idamax( M->x );
      st->xsupn = fabs( gsl_vector_get( M->x, imax ) );
      st->xeucn = gsl_blas_dnrm2 ( M->x ); 

      /* Until now S = X_prev, now S = X - X_prev 
       * Compute s = x_{k+1} - x_k = X - S
       * and     y = g_{k+1} - g_k = G - Y  */
      gsl_vector_sub  ( st->S, M->x        ); /* S = S - X */
      gsl_vector_scale( st->S, -1.0        ); /* S = -S = X - S_prev */
      gsl_vector_sub  ( st->Y, M->gradient ); /* Y = Y - G */
      gsl_vector_scale( st->Y, -1.0        ); /* Y = -Y = G - Y_prev */

      /* Compute sts = s dot s 
       *         sty = s dot y
       * and     sinf = |s|_inf */
      gsl_blas_ddot( st->S, st->S, &(st->sts) );    
      gsl_blas_ddot( st->S, st->Y, &(st->sty) );
      imax      = gsl_blas_idamax( st->S );
      st->sinf  = fabs( gsl_vector_get( st->S, imax ) );

      /* Compute continuous project gradient */
      projected_gradient( st, M->x, M->gradient );

      /* Update spectral steplength */
      spg_steplength( st, P );

      /* Update trust-region radius */
      if ( P->trtype ) st->cg_delta = GSL_MAX( P->delmin, 10*sqrt( st->sts ) );
      else             st->cg_delta = GSL_MAX( P->delmin, 10*    ( st->sinf) );

      return lsflag;
    }

#endif

#endif

    public:

    mmin_constr_gencan() {
      epsgpen=1.0e-5;
      epsgpsn=1.0e-5;
      fmin=-1.0e99;
      udelta0=-1.0;
      ucgmia=-1.0;
      ucgmib=-1.0;
      cg_src=1.0;
      cg_scre=1.0;
      cg_gpnf=1.0e-5;
      cg_epsi=1.0e-1;
      cg_epsf=1.0e-5;
      cg_epsnqmp=1.0e-4;
      cg_maxitnqmp=5;
      nearlyq=0;
      nint=2.0;
      next=2.0;
      mininterp=4;
      maxextrap=100;
      trtype=0;
      eta=0.9;
      delmin=0.1;
      lspgmi=1.0e-10;
      lspgma=1.0e10;
      theta=1.0e-6;
      gamma=1.0e-4;
      beta=0.5;
      sigma1=0.1;
      sigma2=0.9;
      epsrel=1.0e-7;
      epsabs=1.0e-10;
      infrel=1.0e20;
      infabs=1.0e99;
    }

    /// Tolerance on Euclidean norm of projected gradient (default 1.0e-5)
    double epsgpen;
    /// Tolerance on infinite norm of projected gradient (default 1.0e-5)
    double epsgpsn;
    /** \brief Minimum function value (default \f$ 10^{-99} \f$)
	
	If the function value is below this value, then the algorithm
	assumes that the function is not bounded and exits.
    */
    double fmin;
    /// Trust-region radius (default -1.0)
    double udelta0;
    /// Maximum interations per variable (default -1.0)
    double ucgmia;
    /// Extra maximum iterations (default -1.0)
    double ucgmib;
    /// Conjugate gradient condition type (default 1)
    int cg_scre;
    /// Projected gradient norm (default 1.0e-5)
    double cg_gpnf;
    /// Desc (default 1.0e-1)
    double cg_epsi;
    /// Desc (default 1.0e-5)
    double cg_epsf;
    /// Stopping fractional tolerance for conjugate gradient (default 1.0e-4)
    double cg_epsnqmp;
    /// Maximum iterations for conjugate gradient (default 5)
    int cg_maxitnqmp;
    /// Set to 1 if the function is nearly quadratic (default 0)
    int nearlyq;
    /// Interpolation constant (default 2.0)
    double nint;
    /// Extrapolation constant (default 2.0)
    double next;
    /// Minimum interpolation size (default 4)
    int mininterp;
    /// Maximum extrapolations in truncated Newton direction (default 100)
    int maxextrap;
    /// Type of trust region (default 0)
    int trtype;
    /// Threshold for abandoning current face (default 0.9)
    double eta;
    /// Minimum trust region for truncated Newton direction (default 0.1)
    double delmin;
    /// Minimum spectral steplength (default 1.0e-10)
    double lspgmi;
    /// Maximum spectral steplength (default 1.0e10)
    double lspgma;
    /// Constant for the angle condition (default 1.0e-6)
    double theta;
    /// Constant for Armijo condition (default 1.0e-4)
    double gamma;
    /// Constant for beta condition (default 0.5)
    double beta;
    /// Lower bound to the step length reduction (default 0.1)
    double sigma1;
    /// Upper bound to the step length reduction (default 0.9)
    double sigma2;
    /// Relative small number (default 1.0e-7)
    double epsrel;
    /// Absolute small number (default 1.0e-10)
    double epsabs;
    /// Relative infinite number (default 1.0e20)
    double infrel;
    /// Absolute infinite number (default 1.0e99)
    double infabs;
    
    /// Allocate memory
    virtual int alloc(const size_t n) {
      if (this->dim>0) free();
      this->ao.allocate(xx,n);
      this->ao.allocate(d,n);
      this->ao.allocate(s,n);
      this->ao.allocate(y,n);
      return ool_constr_min<param_t,func_t,dfunc_t,hfunc_t,vec_t,vec_t,
	alloc_t>::alloc(n);
    }
    
    /// Free previously allocated memory
    virtual int free() {
      if (this->dim>0) this->ao.free(xx);
      return ool_constr_min<param_t,func_t,dfunc_t,hfunc_t,vec_t,vec_t,
      alloc_t>::free();
    }
    
    /// Set the function, the initial guess, and the parameters
    virtual int set(func_t &fn, dfunc_t &dfn, hfunc_t &hfn,
		    vec_t &init, param_t &par) {
      
      int ret=ool_constr_min<param_t,func_t,dfunc_t,hfunc_t,vec_t,vec_t,
      alloc_t>::set(fn,dfn,hfn,init,par);

#ifdef NEVER_DEFINED
      // Internal variables 
      size_t nn = M->x->size;
      
      // Checking dimensions 
      if( nn != st->n || nn != M->fdf->n || nn != M->con->n  )
	{
	  return OOL_EINVAL;
	}
      
      // Copy boundary vectors 
      gsl_vector_memcpy( st->L, M->con->L );
      gsl_vector_memcpy( st->U, M->con->U );
      
#endif

      prepare_iteration();

      return 0;
    }

#ifdef NEVER_DEFINED

    int prepare_iteration {

      /* Direct access to vector data */
      double *x =  M->x->data;
      double *l = st->L->data;
      double *u = st->U->data;

      /* Internal variables */
      size_t nn = M->x->size;
      size_t ii, imax;

      /* Impose factibility */
      conmin_vector_maxofmin( nn, x, l, u, x );

      /* Eval Euclidean and infinity norms of X */
      st->xeucn = gsl_blas_dnrm2 ( M->x );
      imax      = gsl_blas_idamax( M->x );
      st->xsupn = fabs( gsl_vector_get( M->x, imax ) );

      /* Evaluate objective function and its gradient */
      OOL_CONMIN_EVAL_FDF( M, M->x, &(M->f), M->gradient );

      /* Define near_l and near_u vector */
      for (ii=0; ii < nn; ii++){
	st->near_l[ii] = l[ii] + GSL_MAX( P->epsrel*fabs( l[ii] ), P->epsabs );
	st->near_u[ii] = u[ii] - GSL_MAX( P->epsrel*fabs( u[ii] ), P->epsabs );
      }

      /* Setting constant parameters */
      st->ometa2 = gsl_pow_2( 1.0 - P->eta );
      st->epsgpen2 = gsl_pow_2( P->epsgpen );

      /* Compute continuous project gradient */
      projected_gradient( st, M->x, M->gradient );

      /* Compute a linear relation between gpeucn2 and cgeps2, i.e.,
       * scalars a and b such that
       *
       *     a * log10(||g_P(x_0)||_2^2) + b = log10(cgeps_0^2) and
       *
       *     a * log10(||g_P(x_f)||_2^2) + b = log10(cgeps_f^2),
       *
       *  where cgeps_0 and cgeps_f are provided. Note that if
       *  cgeps_0 is equal to cgeps_f then cgeps will be always
       *  equal to cgeps_0 and cgeps_f.
       *
       *  We introduce now a linear relation between gpsupn and cgeps also.
       */
      if (P->cg_scre == 1 ) {
	st->acgeps = 2 *( log10( P->cg_epsf / P->cg_epsi ) /
			  log10( P->cg_gpnf * P->cg_gpnf / st->gpeucn2 ));

	st->bcgeps = 2 * log10( P->cg_epsi ) -
	  st->acgeps * log10( st->gpeucn2 );
      } else {
	st->acgeps = ( log10( P->cg_epsf / P->cg_epsi ) /
		       log10( P->cg_gpnf / st->gpsupn ) );
	st->bcgeps = log10( P->cg_epsi ) - st->acgeps * log10( st->gpsupn );
      }

      /*     And it will be used for the linear relation of cgmaxit */
      st->gpsupn0  = st->gpsupn;
      st->gpeucn20 = st->gpeucn2;

      /* Initialize the spectral steplength */
      if ( st->gpeucn2 != 0.0 ) {
	st->lambda = GSL_MAX( 1.0, st->xeucn ) / sqrt( st->gpeucn2 );
      }

      /* Initialize the trust-region radius */
      if (P->udelta0 < 0.0 ) {

	double aux;
	if ( P->trtype ) {
	  aux = 0.1 * GSL_MAX( 1.0, st->xeucn );
	} else  {
	  aux = 0.1 * GSL_MAX( 1.0, st->xsupn );
	}
	
	st->cg_delta = GSL_MAX( P->delmin, aux );
	
      } else {
	st->cg_delta = GSL_MAX( P->delmin, P->udelta0 );
      }

      /* Export size */
      M->size = st->gpsupn;

      return OOL_SUCCESS;
    }

#endif

    /// Restart the minimizer
    virtual int restart() {

      /*
      // Restarting dx 
      gsl_vector_set_zero( M->dx );
      
      return prepare_iteration( M );
      */

      return 0;
    }

    /// Perform an iteration
    virtual int iterate() {

#ifdef NEVER_DEFINED

      int status;

      status = actual_iterate( M, st, P );

      /* Export size and dx variables */
      M->size = st->gpsupn;

      /* In the next version does dx replace st->S ? */
      gsl_vector_memcpy( M->dx, st->S );

      return status;

#endif

      return 0;
    }

    /// See if we're finished
    virtual int is_optimal() {

      //return (( st->gpeucn2 <= st->epsgpen2 ||
      //st->gpsupn  <= P->epsgpsn   ||
      //M->f        <= P->fmin      )? OOL_SUCCESS : OOL_CONTINUE );

    }

    /// Return string denoting type ("mmin_constr_gencan")
    const char *type() { return "mmin_constr_gencan"; }

#ifndef DOXYGEN_INTERNAL

  private:
  
  mmin_constr_gencan<func_t,dfunc_t,hfunc_t,vec_t>
  (const mmin_constr_gencan<func_t,dfunc_t,hfunc_t,vec_t> &);
  mmin_constr_gencan<func_t,dfunc_t,hfunc_t,vec_t>& operator=
  (const mmin_constr_gencan<func_t,dfunc_t,hfunc_t,vec_t>&);

#endif
      
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

