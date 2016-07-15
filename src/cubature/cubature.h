/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015-2016, Andrew W. Steiner
  
  This file is part of O2scl. It has been adapted from cubature 
  written by Steven G. Johnson. 
  
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
/* Adaptive multidimensional integration of a vector of integrands.
 *
 * Copyright (c) 2005-2013 Steven G. Johnson
 *
 * Portions (see comments) based on HIntLib (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 2002-2005 Rudolf Schuerer.
 *     (http://www.cosy.sbg.ac.at/~rschuer/hintlib/)
 *
 * Portions (see comments) based on GNU GSL (also distributed under
 * the GNU GPL, v2 or later), copyright (c) 1996-2000 Brian Gough.
 *     (http://www.gnu.org/software/gsl/)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */
/** \file cubature.h
    \brief File for definitions of \ref o2scl::inte_hcubature and 
    \ref o2scl::inte_pcubature
*/
#ifndef O2SCL_CUBATURE_H
#define O2SCL_CUBATURE_H

#include <cmath>
#include <functional>
#include <boost/numeric/ublas/vector.hpp>

#ifndef O2SCL_CLENCURT_H
#define O2SCL_CLENCURT_H
#include <o2scl/clencurt.h>
#endif
#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>

namespace o2scl {

  /** \brief Base class for integration routines from the 
      Cubature library
  */
  class inte_cubature_base {
    
  public:
    
    /** \brief Different ways of measuring the absolute and 
	relative error

	Error estimates given a vector e of error estimates in the
	individual components of a vector v of integrands. These are
	all equivalent when there is only a single integrand.
    */
    typedef enum {
      /** \brief individual relerr criteria in each component */
      ERROR_INDIVIDUAL = 0, 
      /** \brief paired L2 norms of errors in each component,
	  mainly for integrating vectors of complex numbers */
      ERROR_PAIRED, 
      /** \brief abserr is L_2 norm |e|, and relerr is |e|/|v| */
      ERROR_L2, 
      /** \brief abserr is L_1 norm |e|, and relerr is |e|/|v| */
      ERROR_L1, 
      /** \brief abserr is \f$ L_{\infty} \f$ norm |e|, and 
	  relerr is |e|/|v| */
      ERROR_LINF 
    } error_norm;
    
  };

  /** \brief Adaptive multidimensional integration on hyper-rectangles
      using cubature rules from the Cubature library

      This class is experimental.

      \hline
      \b Documentation \b adapted \b from \b Cubature
      
      A cubature rule takes a function and a hypercube and evaluates
      the function at a small number of points, returning an estimate
      of the integral as well as an estimate of the error, and also
      a suggested dimension of the hypercube to subdivide.
      
      Given such a rule, the adaptive integration is simple:
      
      1) Evaluate the cubature rule on the hypercube(s).
      Stop if converged.
      
      2) Pick the hypercube with the largest estimated error,
      and divide it in two along the suggested dimension.
      
      3) Goto (1).
      
      The basic algorithm is based on the adaptive cubature described
      in \ref Genz80 and subsequently extended to integrating a vector
      of integrands in \ref Berntsen91 .

      \hline      

  */
  template<class func_t> class inte_hcubature
    : public inte_cubature_base {
    
  protected:

    /** \brief Desc
     */
    class esterr {
    public:
      /** \brief Desc */
      double val;
      /** \brief Desc */
      double err;
    };

    /** \brief Desc
     */
    double errMax(size_t fdim, const esterr *ee) {

      double errmax = 0;
      for (size_t k = 0; k < fdim; ++k) {
	if (ee[k].err > errmax) errmax = ee[k].err;
      }
      return errmax;
    }

    /** \brief Desc
     */
    class hypercube {
    public:
      /** \brief Desc */
      size_t dim;
      /** \brief length 2*dim = center followed by half-widths */
      double *data; 
      /** \brief cache volume = product of widths */
      double vol;   
    };

    /** \brief Desc
     */
    double compute_vol(const hypercube &h) {
      double vol = 1;
      for (size_t i = 0; i < h.dim; ++i) {
	vol *= 2 * h.data[i + h.dim];
      }
      return vol;
    }

    /** \brief Desc
     */
    void make_hypercube(size_t dim, const std::vector<double> &center,
			const std::vector<double> &halfwidth, hypercube &h) {

      h.dim = dim;
      h.data = (double *) malloc(sizeof(double) * dim * 2);
      h.vol = 0;
      for (size_t i = 0; i < dim; ++i) {
	h.data[i] = center[i];
	h.data[i + dim] = halfwidth[i];
      }
      h.vol = compute_vol(h);
      return;
    }

    /** \brief Desc
     */
    void make_hypercube2(size_t dim, const double *dat, hypercube &h) {

      h.dim = dim;
      h.data = (double *) malloc(sizeof(double) * dim * 2);
      h.vol = 0;
      for (size_t i = 0; i < dim; ++i) {
	h.data[i] = dat[i];
	h.data[i + dim] = dat[i+dim];
      }
      h.vol = compute_vol(h);
      return;
    }
    
    /** \brief Desc
     */
    void make_hypercube_range
      (size_t dim, const std::vector<double> &xmin,
       const std::vector<double> &xmax, hypercube &h) {

      make_hypercube(dim,xmin,xmax,h);
      for (size_t i = 0; i < dim; ++i) {
	h.data[i] = 0.5 * (xmin[i] + xmax[i]);
	h.data[i + dim] = 0.5 * (xmax[i] - xmin[i]);
      }
      h.vol = compute_vol(h);
      return;
    }

    /** \brief Desc
     */
    void destroy_hypercube(hypercube &h) {
      free(h.data);
      h.dim = 0;
      return;
    }

    /** \brief Desc
     */
    class region {
    public:
      /** \brief Desc */
      hypercube h;
      /** \brief Desc */
      size_t splitDim;
      /** \brief dimensionality of vector integrand */
      size_t fdim; 
      /** \brief array of length fdim */
      esterr *ee; 
      /** \brief max ee[k].err */
      double errmax; 
    };

    /** \brief Desc
     */
    void make_region(const hypercube &h, size_t fdim, region &R) {

      make_hypercube2(h.dim, h.data,R.h);
      R.splitDim = 0;
      R.fdim = fdim;
      R.ee=(esterr *) malloc(sizeof(esterr) * fdim);
      R.errmax = HUGE_VAL;

      return;
    }

    /** \brief Desc
     */
    int cut_region(region &R, region &R2) {

      size_t d = R.splitDim, dim = R.h.dim;
      R2=R;
      R.h.data[d + dim] *= 0.5;
      R.h.vol *= 0.5;
      make_hypercube2(dim, R.h.data,R2.h);
      R.h.data[d] -= R.h.data[d + dim];
      R2.h.data[d] += R.h.data[d + dim];
      R2.ee = (esterr *) malloc(sizeof(esterr) * R2.fdim);
      return 0;
    }

    /** \brief Desc
     */
    class rule {

    public:
      
      /** \brief The dimensionality and the number of functions */
      size_t dim, fdim;
      /** \brief The number of evaluation points */
      size_t num_points;
      /** \brief The max number of regions evaluated at once */
      size_t num_regions;
      /** \brief points to eval: num_regions * num_points * dim */
      double *pts;
      /** \brief num_regions * num_points * fdim */
      double *vals;
    };

    /** \brief Desc
     */
    void alloc_rule_pts(rule &r, size_t num_regions) {
      if (num_regions > r.num_regions) {
	free(r.pts);
	r.pts = r.vals = 0;
	r.num_regions = 0;

	/* Allocate extra so that repeatedly calling alloc_rule_pts
	   with growing num_regions only needs a logarithmic number of
	   allocations 
	*/
	num_regions *= 2; 

	r.pts = (double *) malloc(sizeof(double) * 
				  (num_regions
				   * r.num_points * (r.dim + r.fdim)));
	r.vals = r.pts + num_regions * r.num_points * r.dim;
	r.num_regions = num_regions;
      }
      return;
    }

    /** \brief Desc
     */
    rule *make_rule(size_t sz, /* >= sizeof(rule) */
		    size_t dim, size_t fdim, size_t num_points) {
      
      rule *r;

      r = (rule *) malloc(sz);
      r->pts = r->vals = 0;
      r->num_regions = 0;
      r->dim = dim;
      r->fdim = fdim;
      r->num_points = num_points;
      return r;
    }

    /** \brief Desc

	\note All regions must have same fdim 
    */
    int eval_regions(size_t nR, region *R, func_t &f, rule *r) {

      size_t iR;
      if (nR == 0) {
	/* nothing to evaluate */
	return o2scl::success;
      }
      if (r->dim==1) {
	if (rule15gauss_evalError(r, R->fdim, f, nR, R)) {
	  return o2scl::gsl_failure;
	}
      } else {
	if (rule75genzmalik_evalError(r, R->fdim, f, nR, R)) {
	  return o2scl::gsl_failure;
	}
      }
      for (iR = 0; iR < nR; ++iR) {
	R[iR].errmax = errMax(R->fdim, R[iR].ee);
      }
      return o2scl::success;
    }

    /** \brief Functions to loop over points in a hypercube.
    
	Based on orbitrule.cpp in HIntLib-0.0.10
    
	ls0 returns the least-significant 0 bit of n (e.g. it returns
	0 if the LSB is 0, it returns 1 if the 2 LSBs are 01, etcetera).
    */
    static size_t ls0(size_t n) {

#if defined(__GNUC__) &&                                        \
  ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ > 3)
      return __builtin_ctz(~n); /* gcc builtin for version >= 3.4 */
#else
      const size_t bits[256] = {
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 7,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4,
	0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 8,
      };
      size_t bit = 0;
      while ((n & 0xff) == 0xff) {
	n >>= 8;
	bit += 8;
      }
      return bit + bits[n & 0xff];
#endif
    }

    /** \brief Evaluate the integration points for all \f$ 2^n \f$ points 
	(+/-r,...+/-r)
	
	A Gray-code ordering is used to minimize the number of
	coordinate updates in p, although this doesn't matter as much
	now that we are saving all pts.
    */
    void evalR_Rfs(double *pts, size_t dim, double *p,
		   const double *c, const double *r) {
      
      size_t i;
      /* 0/1 bit = +/- for corresponding element of r[] */
      size_t signs = 0; 

      /* We start with the point where r is ADDed in every coordinate
	 (this implies signs=0). */
      for (i = 0; i < dim; ++i) {
	p[i] = c[i] + r[i];
      }

      /* Loop through the points in Gray-code ordering */
      for (i = 0;; ++i) {
	size_t mask, d;

	for(size_t k=0;k<dim;k++) pts[k]=p[k];
	pts += dim;

	d = ls0(i); /* which coordinate to flip */
	if (d >= dim) {
	  break;
	}

	/* flip the d-th bit and add/subtract r[d] */
	mask = 1U << d;
	signs ^= mask;
	p[d] = (signs & mask) ? c[d] - r[d] : c[d] + r[d];
      }
      return;
    }

    /** \brief Desc
     */
    void evalRR0_0fs(double *pts, size_t dim, double *p,
		     const double *c, const double *r) {
      
      for (size_t i = 0; i < dim - 1; ++i) {
	p[i] = c[i] - r[i];
	for (size_t j = i + 1; j < dim; ++j) {
	  p[j] = c[j] - r[j];
	  for(size_t k=0;k<dim;k++) pts[k]=p[k];
	  pts += dim;
	  p[i] = c[i] + r[i];
	  for(size_t k=0;k<dim;k++) pts[k]=p[k];
	  pts += dim;
	  p[j] = c[j] + r[j];
	  for(size_t k=0;k<dim;k++) pts[k]=p[k];
	  pts += dim;
	  p[i] = c[i] - r[i];
	  for(size_t k=0;k<dim;k++) pts[k]=p[k];
	  pts += dim;

	  // Done with j -> Restore p[j]
	  p[j] = c[j];      
	}
	// Done with i -> Restore p[i]
	p[i] = c[i];                
      }
      return;
    }

    /** \brief Desc
     */
    void evalR0_0fs4d
      (double *pts, size_t dim, double *p, const double *c,
       const double *r1, const double *r2) {
      
      for(size_t k=0;k<dim;k++) pts[k]=p[k];
      pts += dim;

      for (size_t i = 0; i < dim; i++) {
	p[i] = c[i] - r1[i];
	for(size_t k=0;k<dim;k++) pts[k]=p[k];
	pts += dim;

	p[i] = c[i] + r1[i];
	for(size_t k=0;k<dim;k++) pts[k]=p[k];
	pts += dim;

	p[i] = c[i] - r2[i];
	for(size_t k=0;k<dim;k++) pts[k]=p[k];
	pts += dim;

	p[i] = c[i] + r2[i];
	for(size_t k=0;k<dim;k++) pts[k]=p[k];
	pts += dim;

	p[i] = c[i];
      }
      return;
    }

    /** \brief Desc
     */
    size_t num0_0(size_t dim) { return 1; }
    /** \brief Desc
     */
    size_t numR0_0fs(size_t dim) { return 2*dim; }
    /** \brief Desc
     */
    size_t numRR0_0fs(size_t dim) { return 2*dim*(dim-1); }
    /** \brief Desc
     */
    size_t numR_Rfs(size_t dim) { return 1U << dim; }
      
    /** \brief Desc

	Based on rule75genzmalik.cpp in HIntLib-0.0.10: An embedded
	cubature rule of degree 7 (embedded rule degree 5) 
	from \ref Genz83.
    */
    class rule75genzmalik : public rule {

    public:

      /** \brief temporary arrays of length dim */
      double *widthLambda;
      /** \brief Desc */
      double *widthLambda2;
      /** \brief Desc */
      double *p;
      
      /** \brief dimension-dependent constants */
      double weight1;
      /** \brief Desc */
      double weight3;
      /** \brief Desc */
      double weight5;
      /** \brief Desc */
      double weightE1;
      /** \brief Desc */
      double weightE3;
    };
    
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
    /** \brief Desc
     */
    int rule75genzmalik_evalError
      (rule *r_, size_t fdim, func_t &f, size_t nR, region *R) {
    
      /* lambda2 = sqrt(9/70), lambda4 = sqrt(9/10), lambda5 = sqrt(9/19) */
      const double lambda2 = 0.3585685828003180919906451539079374954541;
      const double lambda4 = 0.9486832980505137995996680633298155601160;
      const double lambda5 = 0.6882472016116852977216287342936235251269;
      const double weight2 = 980. / 6561.;
      const double weight4 = 200. / 19683.;
      const double weightE2 = 245. / 486.;
      const double weightE4 = 25. / 729.;
      const double ratio = (lambda2 * lambda2) / (lambda4 * lambda4);

      rule75genzmalik *r = (rule75genzmalik *) r_;
      size_t i, j, iR, dim = r_->dim;
      size_t npts = 0;
      double *diff, *pts, *vals;

      alloc_rule_pts(*r_, nR);
      pts = r_->pts; vals = r_->vals;

      for (iR = 0; iR < nR; ++iR) {
	const double *center = R[iR].h.data;
          
	for (i = 0; i < dim; ++i) {
	  r->p[i] = center[i];
	}
          
	for (i = 0; i < dim; ++i) {
	  r->widthLambda2[i] = center[i+dim] * lambda2;
	}
	for (i = 0; i < dim; ++i) {
	  r->widthLambda[i] = center[i+dim] * lambda4;
	}

	/* Evaluate points in the center, in (lambda2,0,...,0) and
	   (lambda3=lambda4, 0,...,0).  */
	evalR0_0fs4d(pts + npts*dim, dim, r->p, center, 
		     r->widthLambda2, r->widthLambda);
	npts += num0_0(dim) + 2 * numR0_0fs(dim);

	/* Calculate points for (lambda4, lambda4, 0, ...,0) */
	evalRR0_0fs(pts + npts*dim, dim, r->p, center, r->widthLambda);
	npts += numRR0_0fs(dim);

	/* Calculate points for (lambda5, lambda5, ..., lambda5) */
	for (i = 0; i < dim; ++i) {
	  r->widthLambda[i] = center[i+dim] * lambda5;
	}
	evalR_Rfs(pts + npts*dim, dim, r->p, center, r->widthLambda);
	npts += numR_Rfs(dim);
      }

      /* Evaluate the integrand function(s) at all the points */
      if (f(dim, npts, pts, fdim, vals)) {
	return o2scl::gsl_failure;
      }

      /* we are done with the points, and so we can re-use the pts
	 array to store the maximum difference diff[i] in each dimension 
	 for each hypercube */
      diff = pts;
      for (i = 0; i < dim * nR; ++i) diff[i] = 0;
      
      for (j = 0; j < fdim; ++j) {
    
	const double *v = vals + j;
    
	for (iR = 0; iR < nR; ++iR) {
	  double result, res5th;
	  double val0, sum2=0, sum3=0, sum4=0, sum5=0;
	  size_t k, k0 = 0;
	  /* accumulate j-th function values into j-th integrals
	     NOTE: this relies on the ordering of the eval functions
	     above, as well as on the internal structure of
	     the evalR0_0fs4d function */

	  val0 = v[fdim*(0)]; /* central point */
	  k0 += 1;

	  for (k = 0; k < dim; ++k) {
	    double v0 = v[fdim*(k0 + 4*k)];
	    double v1 = v[fdim*((k0 + 4*k) + 1)];
	    double v2 = v[fdim*((k0 + 4*k) + 2)];
	    double v3 = v[fdim*((k0 + 4*k) + 3)];
              
	    sum2 += v0 + v1;
	    sum3 += v2 + v3;
    
	    diff[iR * dim + k] += 
	      fabs(v0 + v1 - 2*val0 - ratio * (v2 + v3 - 2*val0));
	  }
#ifdef O2SCL_NEVER_DEFINED
	}{
#endif
	  k0 += 4*k;
    
	  for (k = 0; k < numRR0_0fs(dim); ++k) {
	    sum4 += v[fdim*(k0 + k)];
	  }
	  k0 += k;
    
	  for (k = 0; k < numR_Rfs(dim); ++k) {
	    sum5 += v[fdim*(k0 + k)];
	  }
    
	  /* Calculate fifth and seventh order results */
	  result = R[iR].h.vol * (r->weight1 * val0 + weight2 * sum2 +
				  r->weight3 * sum3 + weight4 * sum4 +
				  r->weight5 * sum5);
	  res5th = R[iR].h.vol * (r->weightE1 * val0 + weightE2 * sum2 +
				  r->weightE3 * sum3 + weightE4 * sum4);
    
	  R[iR].ee[j].val = result;
	  R[iR].ee[j].err = fabs(res5th - result);
    
	  v += r_->num_points * fdim;
	}
      }

      /* figure out dimension to split: */
      for (iR = 0; iR < nR; ++iR) {
	double maxdiff = 0;
	size_t dimDiffMax = 0;
  
	for (i = 0; i < dim; ++i) {
	  if (diff[iR*dim + i] > maxdiff) {
	    maxdiff = diff[iR*dim + i];
	    dimDiffMax = i;
	  }
	}
	R[iR].splitDim = dimDiffMax;
      }
      return o2scl::success;
    }

#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
    /** \brief Desc
     */
    rule *make_rule75genzmalik(size_t dim, size_t fdim) {

      rule75genzmalik *r;
      
      if (dim < 2) return 0; /* this rule does not support 1d integrals */
      
      /* Because of the use of a bit-field in evalR_Rfs, we are limited
	 to be < 32 dimensions (or however many bits are in size_t).
	 This is not a practical limitation...long before you reach
	 32 dimensions, the Genz-Malik cubature becomes excruciatingly
	 slow and is superseded by other methods (e.g. Monte-Carlo). */
      if (dim >= sizeof(size_t) * 8) return 0;
      
      r = (rule75genzmalik *) make_rule(sizeof(rule75genzmalik),
					dim, fdim,
					num0_0(dim) + 2 * numR0_0fs(dim)
					+ numRR0_0fs(dim) + numR_Rfs(dim));

      r->weight1=(12824.0-9120.0*dim+400.0*dim*dim)/19683.0;
      r->weight3=(1820.0-400.0*dim)/19683.0;
      r->weight5=6859.0/19683.0/((double)(1U << dim));
      r->weightE1=(729.0-950.0*dim+50.0*dim*dim)/729.0;
      r->weightE3=(265.0-100.0*dim)/1458.0;

      r->p = (double *) malloc(sizeof(double) * dim * 3);
      r->widthLambda = r->p + dim;
      r->widthLambda2 = r->p + 2 * dim;

      return (rule *) r;
    }

    /** \brief 1d 15-point Gaussian quadrature rule, based on qk15.c
	and qk.c in GNU GSL (which in turn is based on QUADPACK).
    */
    int rule15gauss_evalError
      (rule *r, size_t fdim, func_t &f, size_t nR, region *R) {

      static const double cub_dbl_min=std::numeric_limits<double>::min();
      static const double cub_dbl_eps=std::numeric_limits<double>::epsilon();

      /* Gauss quadrature weights and kronrod quadrature abscissae and
	 weights as evaluated with 80 decimal digit arithmetic by
	 L. W. Fullerton, Bell Labs, Nov. 1981. */
      const size_t n = 8;
      const double xgk[8] = {  /* abscissae of the 15-point kronrod rule */
	0.991455371120812639206854697526329,
	0.949107912342758524526189684047851,
	0.864864423359769072789712788640926,
	0.741531185599394439863864773280788,
	0.586087235467691130294144838258730,
	0.405845151377397166906606412076961,
	0.207784955007898467600689403773245,
	0.000000000000000000000000000000000
	/* xgk[1], xgk[3], ... abscissae of the 7-point gauss rule. 
	   xgk[0], xgk[2], ... to optimally extend the 7-point gauss rule */
      };
      static const double wg[4] = {  /* weights of the 7-point gauss rule */
	0.129484966168869693270611432679082,
	0.279705391489276667901467771423780,
	0.381830050505118944950369775488975,
	0.417959183673469387755102040816327
      };
      static const double wgk[8] = { /* weights of the 15-point kronrod rule */
	0.022935322010529224963732008058970,
	0.063092092629978553290700663189204,
	0.104790010322250183839876322541518,
	0.140653259715525918745189590510238,
	0.169004726639267902826583426598550,
	0.190350578064785409913256402421014,
	0.204432940075298892414161999234649,
	0.209482141084727828012999174891714
      };
      size_t j, k, iR;
      size_t npts = 0;
      double *pts, *vals;

      alloc_rule_pts(*r, nR);
      pts = r->pts; vals = r->vals;

      for (iR = 0; iR < nR; ++iR) {
	const double center = R[iR].h.data[0];
	const double halfwidth = R[iR].h.data[1];

	pts[npts++] = center;

	for (j = 0; j < (n - 1) / 2; ++j) {
	  int j2 = 2*j + 1;
	  double w = halfwidth * xgk[j2];
	  pts[npts++] = center - w;
	  pts[npts++] = center + w;
	}
	for (j = 0; j < n/2; ++j) {
	  int j2 = 2*j;
	  double w = halfwidth * xgk[j2];
	  pts[npts++] = center - w;
	  pts[npts++] = center + w;
	}

	R[iR].splitDim = 0; /* no choice but to divide 0th dimension */
      }

      if (f(1, npts, pts, fdim, vals)) {
	return o2scl::gsl_failure;
      }
     
      for (k = 0; k < fdim; ++k) {
	const double *vk = vals + k;
	for (iR = 0; iR < nR; ++iR) {
	  const double halfwidth = R[iR].h.data[1];
	  double result_gauss = vk[0] * wg[n/2 - 1];
	  double result_kronrod = vk[0] * wgk[n - 1];
	  double result_abs = fabs(result_kronrod);
	  double result_asc, mean, err;

	  /* accumulate integrals */
	  npts = 1;
	  for (j = 0; j < (n - 1) / 2; ++j) {
	    int j2 = 2*j + 1;
	    double v = vk[fdim*npts] + vk[fdim*npts+fdim];
	    result_gauss += wg[j] * v;
	    result_kronrod += wgk[j2] * v;
	    result_abs += wgk[j2] * (fabs(vk[fdim*npts]) 
				     + fabs(vk[fdim*npts+fdim]));
	    npts += 2;
	  }
	  for (j = 0; j < n/2; ++j) {
	    int j2 = 2*j;
	    result_kronrod += wgk[j2] * (vk[fdim*npts] 
					 + vk[fdim*npts+fdim]);
	    result_abs += wgk[j2] * (fabs(vk[fdim*npts]) 
				     + fabs(vk[fdim*npts+fdim]));
	    npts += 2;
	  }
               
	  /* integration result */
	  R[iR].ee[k].val = result_kronrod * halfwidth;

	  /* error estimate 
	     (from GSL, probably dates back to QUADPACK
	     ... not completely clear to me why we don't just use
	     fabs(result_kronrod - result_gauss) * halfwidth */
	  mean = result_kronrod * 0.5;
	  result_asc = wgk[n - 1] * fabs(vk[0] - mean);
	  npts = 1;
	  for (j = 0; j < (n - 1) / 2; ++j) {
	    int j2 = 2*j + 1;
	    result_asc += wgk[j2] * (fabs(vk[fdim*npts]-mean)
				     + fabs(vk[fdim*npts+fdim]-mean));
	    npts += 2;
	  }
	  for (j = 0; j < n/2; ++j) {
	    int j2 = 2*j;
	    result_asc += wgk[j2] * (fabs(vk[fdim*npts]-mean)
				     + fabs(vk[fdim*npts+fdim]-mean));
	    npts += 2;
	  }
	  err = fabs(result_kronrod - result_gauss) * halfwidth;
	  result_abs *= halfwidth;
	  result_asc *= halfwidth;
	  if (result_asc != 0 && err != 0) {
	    double scale = pow((200 * err / result_asc), 1.5);
	    err = (scale < 1) ? result_asc * scale : result_asc;
	  }
	  if (result_abs > cub_dbl_min / (50 * cub_dbl_eps)) {
	    double min_err = 50 * cub_dbl_eps * result_abs;
	    if (min_err > err) err = min_err;
	  }
	  R[iR].ee[k].err = err;
               
	  /* increment vk to point to next batch of results */
	  vk += 15*fdim;
	}
      }
      return o2scl::success;
    }
     
    /** \brief Desc
     */
    rule *make_rule15gauss(size_t dim, size_t fdim) {

      if (dim != 1) {
	O2SCL_ERR("this rule is only for 1d integrals.",o2scl::exc_esanity);
      }
       
      return make_rule(sizeof(rule),dim,fdim,15);
    }

    /** \name Binary heap implementation
	
	Based on \ref Cormen and used as a priority queue of
	regions to integrate.
    */
    //@{
    /** \brief Desc
     */
    typedef region heap_item;

    /** \brief Desc
     */
    class heap {
    public:
      /** \brief Desc */
      size_t n;
      /** \brief Desc */
      size_t nalloc;
      /** \brief Desc */
      heap_item *items;
      /** \brief Desc */
      size_t fdim;
      /** array of length fdim of the total integrand & error */
      std::vector<esterr> ee;
    };

    /** \brief Desc
     */
    void heap_resize(heap &h, size_t nalloc) {

      h.nalloc = nalloc;
      if (nalloc) {
	h.items = (heap_item *) realloc(h.items, sizeof(heap_item)*nalloc);
      } else {
	/* BSD realloc does not free for a zero-sized reallocation */
	free(h.items);
	h.items = 0;
      }
      return;
    }

    /** \brief Desc
     */
    heap heap_alloc(size_t nalloc, size_t fdim) {

      heap h;
      h.n = 0;
      h.nalloc = 0;
      h.items = 0;
      h.fdim = fdim;
      h.ee.resize(fdim);
      for (size_t i = 0; i < fdim; ++i) h.ee[i].val = h.ee[i].err = 0;
      heap_resize(h, nalloc);
      return h;
    }

    /** \brief Note that heap_free does not deallocate anything referenced by
	the items */
    void heap_free(heap &h) {

      h.n = 0;
      heap_resize(h, 0);
      h.fdim = 0;
      h.ee.clear();
      return;
    }

    /** \brief Desc
     */
    int heap_push(heap &h, heap_item hi) {

      int insert;
      size_t fdim = h.fdim;

      for (size_t i = 0; i < fdim; ++i) {
	h.ee[i].val += hi.ee[i].val;
	h.ee[i].err += hi.ee[i].err;
      }
      insert = h.n;
      if (++(h.n) > h.nalloc) {
	heap_resize(h, h.n * 2);
	if (!h.items) return o2scl::gsl_failure;
      }

      while (insert) {
	int parent = (insert - 1) / 2;
	if (hi.errmax <= h.items[parent].errmax) {
	  break;
	}
	h.items[insert] = h.items[parent];
	insert = parent;
      }
      h.items[insert] = hi;
      return o2scl::success;
    }

    /** \brief Desc
     */
    int heap_push_many(heap &h, size_t ni, heap_item *hi) {
      for (size_t i = 0; i < ni; ++i) {
	if (heap_push(h, hi[i])) return o2scl::gsl_failure;
      }
      return o2scl::success;
    }

    /** \brief Desc
     */
    heap_item heap_pop(heap &h) {

      heap_item ret;
      int i, n, child;

      if (!(h.n)) {
	O2SCL_ERR("Attempted to pop an empty heap in cubature.",
		  o2scl::exc_esanity);
      }

      ret = h.items[0];
      h.items[i = 0] = h.items[n = --(h.n)];

      while ((child = i * 2 + 1) < n) {

	int largest;
	heap_item swap;

	if (h.items[child].errmax <= h.items[i].errmax) {
	  largest = i;
	} else {
	  largest = child;
	}
	
	if (++child < n && h.items[largest].errmax <
	    h.items[child].errmax) {
	  largest = child;
	}
	if (largest == i) {
	  break;
	}
	swap = h.items[i];
	h.items[i] = h.items[largest];
	h.items[i = largest] = swap;
      }

      {
	size_t i, fdim = h.fdim;
	for (i = 0; i < fdim; ++i) {
	  h.ee[i].val -= ret.ee[i].val;
	  h.ee[i].err -= ret.ee[i].err;
	}
      }
      return ret;
    }
    //@}

    /** \brief Desc
     */
    int converged(size_t fdim, const std::vector<esterr> &ee,
		  double reqAbsError, double reqRelError,
		  error_norm norm) {

      size_t j;

      switch (norm) {
	
      case ERROR_INDIVIDUAL:
	
	for (j = 0; j < fdim; ++j) {
	  if (ee[j].err > reqAbsError && ee[j].err >
	      fabs(ee[j].val)*reqRelError) {
	    return 0;
	  }
	}
	return 1;
              
      case ERROR_PAIRED:

	for (j = 0; j+1 < fdim; j += 2) {
	  double maxerr, serr, err, maxval, sval, val;
	  /* scale to avoid overflow/underflow */
	  maxerr = ee[j].err > ee[j+1].err ? ee[j].err : ee[j+1].err;
	  maxval = ee[j].val > ee[j+1].val ? ee[j].val : ee[j+1].val;
	  serr = maxerr > 0 ? 1/maxerr : 1;
	  sval = maxval > 0 ? 1/maxval : 1;
	  err=std::hypot(ee[j].err*serr,ee[j+1].err*serr)*maxerr;
	  val=std::hypot(ee[j].val*sval,ee[j+1].val*sval)*maxval;
	  if (err > reqAbsError && err > val*reqRelError) {
	    return 0;
	  }
	}

	/* fdim is odd, do last dimension individually */
	if (j < fdim) {
	  if (ee[j].err > reqAbsError && ee[j].err >
	      fabs(ee[j].val)*reqRelError) {
	    return 0;
	  }
	}
	
	return 1;
	
      case ERROR_L1:
	{
	  double err = 0, val = 0;
	  for (j = 0; j < fdim; ++j) {
	    err += ee[j].err;
	    val += fabs(ee[j].val);
	  }
	  return err <= reqAbsError || err <= val*reqRelError;
	}
	
      case ERROR_LINF:
	{
	  double err = 0, val = 0;
	  for (j = 0; j < fdim; ++j) {
	    double absval = fabs(ee[j].val);
	    if (ee[j].err > err) err = ee[j].err;
	    if (absval > val) val = absval;
	  }
	  return err <= reqAbsError || err <= val*reqRelError;
	}
	
      case ERROR_L2:
	{
	  double maxerr = 0, maxval = 0, serr, sval, err = 0, val = 0;
	  /* scale values by 1/max to avoid overflow/underflow */
	  for (j = 0; j < fdim; ++j) {
	    double absval = fabs(ee[j].val);
	    if (ee[j].err > maxerr) maxerr = ee[j].err;
	    if (absval > maxval) maxval = absval;
	  }
	  serr = maxerr > 0 ? 1/maxerr : 1;
	  sval = maxval > 0 ? 1/maxval : 1;
	  for (j = 0; j < fdim; ++j) {
	    err += (ee[j].err * serr)*(ee[j].err * serr);
	    val += fabs(ee[j].val) * sval*fabs(ee[j].val) * sval;
	  }
	  err = sqrt(err) * maxerr;
	  val = sqrt(val) * maxval;
	  return err <= reqAbsError || err <= val*reqRelError;
	}
	
      }
      return 1; /* unreachable */
    }

    /** \brief Desc
     */
    int rulecubature(rule &r, size_t fdim, func_t &f, 
		     const hypercube &h, size_t maxEval,
		     double reqAbsError, double reqRelError,
		     error_norm norm, double *val, double *err,
		     int parallel) {
      
      size_t numEval = 0;
      heap regions;
      size_t i, j;
      /* array of regions to evaluate */
      region *R = 0; 
      size_t nR_alloc = 0;
      std::vector<esterr> ee(fdim);

      /* norm is irrelevant */
      if (fdim <= 1) norm = ERROR_INDIVIDUAL; 
      /* invalid norm */
      if (norm < 0 || norm > ERROR_LINF) return o2scl::gsl_failure; 

      regions = heap_alloc(1, fdim);
     
      nR_alloc = 2;
      R = (region *) malloc(sizeof(region) * nR_alloc);
      make_region(h, fdim, R[0]);
      if (!R[0].ee || eval_regions(1, R, f, &r) ||
	  heap_push(regions, R[0])) {
	heap_free(regions);
	free(R);
	return o2scl::gsl_failure;
      }
      numEval += r.num_points;
     
      while (numEval < maxEval || !maxEval) {

	if (converged(fdim, regions.ee, reqAbsError, reqRelError, norm)) {
	  break;
	}

	/* maximize potential parallelism */
	if (parallel) { 

	  /* Adapted from I. Gladwell, "Vectorization of one
	     dimensional quadrature codes," pp. 230--238 in _Numerical
	     Integration. Recent Developments, Software and
	     Applications_, G. Fairweather and P. M. Keast, eds., NATO
	     ASI Series C203, Dordrecht (1987), as described in J. M.
	     Bull and T. L. Freeman, "Parallel Globally Adaptive
	     Algorithms for Multi-dimensional Integration,"
	     http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.42.6638
	     (1994).

	     Basically, this evaluates in one shot all regions that
	     *must* be evaluated in order to reduce the error to the
	     requested bound: the minimum set of largest-error regions
	     whose errors push the total error over the bound.

	     [Note: Bull and Freeman claim that the Gladwell approach
	     is intrinsically inefficent because it "requires
	     sorting", and propose an alternative algorithm that
	     "only" requires three passes over the entire set of
	     regions. Apparently, they didn't realize that one could
	     use a heap data structure, in which case the time to pop
	     K biggest-error regions out of N is only O(K log N), much
	     better than the O(N) cost of the Bull and Freeman
	     algorithm if K << N, and it is also much simpler.] 
	  */
	  size_t nR = 0;
	  for (j = 0; j < fdim; ++j) ee[j] = regions.ee[j];
	  do {

	    if (nR + 2 > nR_alloc) {
	      nR_alloc = (nR + 2) * 2;
	      R = (region *) realloc(R, nR_alloc * sizeof(region));
	    }
	    R[nR] = heap_pop(regions);
	    for (j = 0; j < fdim; ++j) ee[j].err -= R[nR].ee[j].err;
	    numEval += r.num_points * 2;
	    nR += 2;
	    if (converged(fdim, ee, reqAbsError, reqRelError, norm)) {
	      /* other regions have small errs */
	      break; 
	    }
	    
	  } while (regions.n > 0 && (numEval < maxEval || !maxEval));

	  if (eval_regions(nR, R, f, &r)
	      || heap_push_many(regions, nR, R)) {
	    heap_free(regions);
	    free(R);
	    return o2scl::gsl_failure;
	  }

	} else { 

	  /* minimize number of function evaluations */
	  
	  /* get worst region */
	  R[0] = heap_pop(regions); 
	  if (cut_region(R[0], R[1]) || eval_regions(2, R, f, &r)
	      || heap_push_many(regions, 2, R)) {
	    heap_free(regions);
	    free(R);
	    return o2scl::gsl_failure;
	  }
	  numEval += r.num_points * 2;
	}
      }

      /* re-sum integral and errors */
      for (j = 0; j < fdim; ++j) val[j] = err[j] = 0;  
      for (i = 0; i < regions.n; ++i) {
	for (j = 0; j < fdim; ++j) { 
	  val[j] += regions.items[i].ee[j].val;
	  err[j] += regions.items[i].ee[j].err;
	}
	{
	  destroy_hypercube(regions.items[i].h);
	  free(regions.items[i].ee);
	  regions.items[i].ee = 0;
	}
      }

      heap_free(regions);
      free(R);

      return o2scl::success;
    }
    
    /** \brief Desc
     */
    int cubature(size_t fdim, func_t &f, size_t dim,
		 const std::vector<double> &xmin,
		 const std::vector<double> &xmax, 
		 size_t maxEval, double reqAbsError, double reqRelError, 
		 error_norm norm, std::vector<double> &val,
		 std::vector<double> &err, int parallel) {
      
      rule *r;
      hypercube h;
      int status;
      size_t i;
      
      if (fdim == 0) {
	/* nothing to do */
	return o2scl::success;
      }
      if (dim == 0) {
	/* trivial integration */
	if (f(0, 1, &(xmin[0]), fdim, &(val[0]))) {
	  return o2scl::gsl_failure;
	}
	for (i = 0; i < fdim; ++i) err[i] = 0;
	return o2scl::success;
      }
      if (dim==1) {
	r=make_rule15gauss(dim,fdim);
      } else {
	r=make_rule75genzmalik(dim,fdim);
      }
      make_hypercube_range(dim,xmin,xmax,h);
      status = rulecubature(*r, fdim, f, h,
			    maxEval, reqAbsError, reqRelError, norm,
			    &(val[0]), &(err[0]), parallel);
      destroy_hypercube(h);

      if (dim==1) {
	rule75genzmalik *r2=(rule75genzmalik *)r;
	free(r2->p);
      }
      free(r->pts);
      free(r);

      return status;
    }
    
  public:
    
    /** \brief Desc
     */
    int integ(size_t fdim, func_t &f, size_t dim,
	      const std::vector<double> &xmin,
	      const std::vector<double> &xmax, size_t maxEval,
	      double reqAbsError, double reqRelError, error_norm norm,
	      std::vector<double> &val, std::vector<double> &err) {
	      
      if (fdim == 0) {
	/* nothing to do */     
	return o2scl::success;
      }
      return cubature(fdim,f,dim,xmin,xmax,
		      maxEval,reqAbsError,reqRelError,norm,val,err,0);
    }
    
  };

#ifdef O2SCL_NEVER_DEFINED
}{
#endif
    
  /** \brief Integration by p-adaptive cubature from the Cubature library

      This class is experimental.

      \hline
      \b Documentation \b adapted \b from \b Cubature
      
      This class performs adaptive integration by increasing the
      degree of the cubature rule rather than subdividing the domain,
      using products of Clenshaw-Curtis rules. This algorithm may be
      superior to Genz-Malik for smooth integrands lacking
      strongly-localized features, in moderate dimensions.
      
      \hline      
  */
  template<class func_t, class vec_t=std::vector<double> >
    class inte_pcubature : public inte_cubature_base {
    
  protected:
  
  /** \brief Maximum integral dimension
   */
  static const size_t MAXDIM=20;
  
  /** \brief Cache of the values for the m[dim] grid.  

      For adaptive cubature, thanks to the nesting of the C-C rules, we
      can re-use the values from coarser grids for finer grids, and the
      coarser grids are also used for error estimation. 
	
      A grid is determined by an m[dim] array, where m[i] denotes
      2^(m[i]+1)+1 points in the i-th dimension.

      If mi < dim, then we only store the values corresponding to
      the difference between the m grid and the grid with m[mi] ->
      m[mi]-1. (m[mi]-1 == -1 corresponds to the trivial grid of one
      point in the center.) 
  */
  class cache {
  public:
    cache() {
      m.resize(MAXDIM);
    }
    /** \brief Desc */
    std::vector<size_t> m;
    /** \brief Desc */
    size_t mi;
    /** \brief Desc */
    std::vector<double> val;
  };

  /** \brief Desc

      recursive loop over all cubature points for the given (m,mi)
      cache entry: add each point to the buffer buf, evaluating all
      at once whenever the buffer is full or when we are done
  */
  int compute_cacheval(const std::vector<size_t> &m, size_t mi, 
		       std::vector<double> &val, size_t &vali,
		       size_t fdim, func_t &f, size_t dim, size_t id,
		       std::vector<double> &p, const vec_t &xmin,
		       const vec_t &xmax, std::vector<double> &buf,
		       size_t nbuf, size_t &ibuf) {

    if (id == dim) {
      /* add point to buffer of points */
      for(size_t k=0;k<dim;k++) {
	buf[k+ibuf*dim]=p[k];
      }
      ibuf++;
      if (ibuf == nbuf) {
	/* flush buffer */
	if (f(dim, nbuf, &(buf[0]), fdim, &(val[0]) + vali)) {
	  return o2scl::gsl_failure;
	}
	vali += ibuf * fdim;
	ibuf = 0;
      }

    } else {
      
      double c = (xmin[id] + xmax[id]) * 0.5;
      double r = (xmax[id] - xmin[id]) * 0.5;
      const double *x = clencurt_x 
	+ ((id == mi) ? (m[id] ? (1 << (m[id] - 1)) : 0) : 0);
      size_t i;
      size_t nx = (id == mi ? (m[id] ? (1 << (m[id] - 1)) : 1)
		   : (1 << (m[id])));
      if (id != mi) {
	p[id] = c;
	if (compute_cacheval(m, mi, val, vali, fdim, f,
			     dim, id + 1, p,
			     xmin, xmax, buf, nbuf, ibuf)) {
	  return o2scl::gsl_failure;
	}
      }
      for (i = 0; i < nx; ++i) {
	p[id] = c + r * x[i];
	if (compute_cacheval(m, mi, val, vali, fdim, f,
			     dim, id + 1, p,
			     xmin, xmax, buf, nbuf, ibuf)) {
	  return o2scl::gsl_failure;
	}
	p[id] = c - r * x[i];
	if (compute_cacheval(m, mi, val, vali, fdim, f,
			     dim, id + 1, p,
			     xmin, xmax, buf, nbuf, ibuf)) {
	  return o2scl::gsl_failure;
	}
      }
    }
    return o2scl::success;
  }

  /** \brief Desc
   */
  size_t num_cacheval(const std::vector<size_t> &m, size_t mi, size_t dim,
		      size_t i_shift) {

    size_t nval = 1;
    for (size_t i = 0; i < dim; ++i) {
      if (i == mi) {
	if (m[i+i_shift]==0) {
	  nval*=2;
	} else {
	  nval*= (1 << (m[i+i_shift]));
	}
      } else {
	nval *= (1 << (m[i+i_shift] + 1)) + 1;
      }
    }
    return nval;
  }
    
  /** \brief Desc
   */
  int add_cacheval(std::vector<cache> &vc, const std::vector<size_t> &m,
		   size_t mi, size_t fdim, func_t &f, size_t dim, 
		   const vec_t &xmin, const vec_t &xmax,
		   std::vector<double> &buf, size_t nbuf) {
      
    size_t ic = vc.size();
    size_t nval, vali = 0, ibuf = 0;
    std::vector<double> p(MAXDIM);

    vc.resize(vc.size()+1);

    vc[ic].mi = mi;
    for(size_t j=0;j<dim;j++) {
      vc[ic].m[j]=m[j];
    }
    nval = fdim * num_cacheval(m,mi,dim,0);
    vc[ic].val.resize(nval);

    if (compute_cacheval(m,mi,vc[ic].val,vali,fdim,f,
			 dim,0,p,xmin,xmax,buf,nbuf,ibuf)) {
      return o2scl::gsl_failure;
    }

    if (ibuf > 0) {
      /* flush remaining buffer */
      return f(dim, ibuf, &(buf[0]), fdim, &((vc[ic].val)[0]) + vali);
    }

    return o2scl::success;
  }
    
  /** \brief Desc
	
      Recursive loop to evaluate the integral contribution from the
      cache entry c, accumulating in val, for the given m[] except
      with m[md] -> m[md] - 1 if md < dim, using the cached values
      (cm,cmi,cval). id is the current loop dimension (from 0 to
      dim-1).
  */
  size_t eval(const std::vector<size_t> &cm, size_t cmi,
	      std::vector<double> &cval,
	      const std::vector<size_t> &m, size_t md,
	      size_t fdim, size_t dim, size_t id,
	      double weight, vec_t &val, size_t voff2) {

    size_t voff = 0; /* amount caller should offset cval array afterwards */
    if (id == dim) {
      size_t i;
      for (i = 0; i < fdim; ++i) {
	val[i] += cval[i+voff2] * weight;
      }
      voff = fdim;

    } else if (m[id] == 0 && id == md) {

      /* using trivial rule for this dim */
      voff = eval(cm, cmi, cval, m, md, fdim, dim, id+1, weight*2, val,voff2);
      voff += fdim * (1 << cm[id]) * 2
      * num_cacheval(cm, cmi - (id+1), dim - (id+1),id+1);

    } else {
      
      size_t i;
      /* order of C-C rule */
      size_t mid = m[id] - (id == md); 
      const double *w = clencurt_w + mid + (1 << mid) - 1
      + (id == cmi ? (cm[id] ? 1 + (1 << (cm[id]-1)) : 1) : 0);
      size_t cnx = (id == cmi ? (cm[id] ? (1 << (cm[id]-1)) : 1)
		    : (1 << (cm[id])));
      size_t nx = cm[id] <= mid ? cnx : (1 << mid);

      if (id != cmi) {
	voff = eval(cm, cmi, cval, m, md, fdim, dim, id + 1,
		    weight * w[0], val,voff2);
	++w;
      }
      for (i = 0; i < nx; ++i) {
	voff += eval(cm, cmi, cval, m, md, fdim, dim, id + 1,
		     weight * w[i], val,voff+voff2);
	voff += eval(cm, cmi, cval, m, md, fdim, dim, id + 1,
		     weight * w[i], val,voff+voff2);
      }

      voff += (cnx - nx) * fdim * 2
      * num_cacheval(cm, cmi - (id+1), dim - (id+1),id+1);
    }
    return voff;
  }

  /** \brief Desc

      Loop over all cache entries that contribute to the integral,
      (with m[md] decremented by 1) 
  */
  void evals(std::vector<cache> &vc, const std::vector<size_t> &m,
	     size_t md,
	     size_t fdim, size_t dim, double V,
	     vec_t &val) {

    for(size_t k=0;k<fdim;k++) {
      val[k]=0.0;
    }
    for (size_t i = 0; i < vc.size(); ++i) {
      if (vc[i].mi >= dim ||
	  vc[i].m[vc[i].mi] + (vc[i].mi == md) <= m[vc[i].mi]) {
	eval(vc[i].m, vc[i].mi, vc[i].val,
	     m, md, fdim, dim, 0, V, val,0);
      }
    }
    return;
  }

  /** \brief Desc

      Evaluate the integrals for the given m[] using the cached
      values in vc, storing the integrals in val[], the error
      estimate in err[], and the dimension to subdivide next (the
      largest error contribution) in *mi
  */
  void eval_integral(std::vector<cache> &vc, const std::vector<size_t> &m, 
		     size_t fdim, size_t dim, double V,
		     size_t &mi, vec_t &val,
		     vec_t &err, vec_t &val1) {

    double maxerr = 0;
    size_t i, j;
     
    evals(vc, m, dim, fdim, dim, V, val);

    /* error estimates along each dimension by comparing val with
       lower-order rule in that dimension; overall (conservative)
       error estimate from maximum error of lower-order rules. */
    for(size_t j=0;j<fdim;j++) {
      err[j]=0.0;
    }
    mi = 0;
    for (i = 0; i < dim; ++i) {
      double emax = 0;
      evals(vc, m, i, fdim, dim, V, val1);
      for (j = 0; j < fdim; ++j) {
	double e = fabs(val[j] - val1[j]);
	if (e > emax) emax = e;
	if (e > err[j]) err[j] = e;
      }
      if (emax > maxerr) {
	maxerr = emax;
	mi = i;
      }
    }
    return;
  }

  /** \brief Desc
   */
  template<class vec2_t>
  int converged(size_t fdim, const vec2_t &vals, const vec2_t &errs,
		double reqAbsError, double reqRelError, error_norm norm) {

    switch (norm) {
	
    case ERROR_INDIVIDUAL:

    for (size_t j = 0; j < fdim; ++j) {
      if (errs[j] > reqAbsError && errs[j] > fabs(vals[j])*reqRelError) {
	return 0;
      }
    }
    return 1;
	      
    case ERROR_PAIRED:

    {
      size_t j;

      for (j = 0; j+1 < fdim; j += 2) {
	double maxerr, serr, err, maxval, sval, val;
	/* scale to avoid overflow/underflow */
	maxerr = errs[j] > errs[j+1] ? errs[j] : errs[j+1];
	maxval = vals[j] > vals[j+1] ? vals[j] : vals[j+1];
	serr = maxerr > 0 ? 1/maxerr : 1;
	sval = maxval > 0 ? 1/maxval : 1;
	err=std::hypot(errs[j]*serr,errs[j+1]*serr)*maxerr;
	val=std::hypot(vals[j]*sval,vals[j+1]*sval)*maxval;
	if (err > reqAbsError && err > val*reqRelError) {
	  return 0;
	}
      }
      /* fdim is odd, do last dimension individually */
      if (j < fdim) {
	if (errs[j] > reqAbsError && errs[j] > fabs(vals[j])*reqRelError) {
	  return 0;
	}
      }
    }
    return 1;

    case ERROR_L1:

    {
      double err = 0, val = 0;
      for (size_t j = 0; j < fdim; ++j) {
	err += errs[j];
	val += fabs(vals[j]);
      }
      return err <= reqAbsError || err <= val*reqRelError;
    }
	
    case ERROR_LINF:

    {
      double err = 0, val = 0;
      for (size_t j = 0; j < fdim; ++j) {
	double absval = fabs(vals[j]);
	if (errs[j] > err) err = errs[j];
	if (absval > val) val = absval;
      }
      return err <= reqAbsError || err <= val*reqRelError;
    }
	
    case ERROR_L2:

    {
      double maxerr = 0, maxval = 0, serr, sval, err = 0, val = 0;
      /* scale values by 1/max to avoid overflow/underflow */
      for (size_t j = 0; j < fdim; ++j) {
	double absval = fabs(vals[j]);
	if (errs[j] > maxerr) maxerr = errs[j];
	if (absval > maxval) maxval = absval;
      }
      serr = maxerr > 0 ? 1/maxerr : 1;
      sval = maxval > 0 ? 1/maxval : 1;
      for (size_t j = 0; j < fdim; ++j) {
	err += (errs[j] * serr)*(errs[j] * serr);
	val += (fabs(vals[j]) * sval)*(fabs(vals[j]) * sval);
      }
      err = sqrt(err) * maxerr;
      val = sqrt(val) * maxval;
      return err <= reqAbsError || err <= val*reqRelError;
    }
    }
      
    O2SCL_ERR("Improper value of 'norm' in cubature::converged().",
	      o2scl::exc_einval);
    return 1;
  }
    
  public:
    
  /** \brief Desc

      Vectorized version with user-supplied buffer to store points
      and values. The buffer *buf should be of length *nbuf * dim on
      entry (these parameters are changed upon return to the final
      buffer and length that was used). The buffer length will be
      kept <= max(max_nbuf, 1) * dim.

      Also allows the caller to specify an array m[dim] of starting degrees
      for the rule, which upon return will hold the final degrees.  The
      number of points in each dimension i is 2^(m[i]+1) + 1. 
  */
  int integ_v_buf(size_t fdim, func_t &f, 
		  size_t dim, const vec_t &xmin, const vec_t &xmax,
		  size_t maxEval, double reqAbsError, double reqRelError,
		  error_norm norm, std::vector<size_t> &m,
		  std::vector<double> &buf, size_t &nbuf, size_t max_nbuf,
		  vec_t &val, vec_t &err) {
      
    int ret = o2scl::gsl_failure;
    double V = 1;
    size_t numEval = 0, new_nbuf;
    size_t i;
    std::vector<cache> vc;

    vec_t val1(fdim);

    /* norm is irrelevant */
    if (fdim <= 1) norm = ERROR_INDIVIDUAL;
    /* invalid norm */
    if (norm < 0 || norm > ERROR_LINF) return o2scl::gsl_failure; 
    /* nothing to do */
    if (fdim == 0) return o2scl::success; 
    /* unsupported */
    if (dim > MAXDIM) return o2scl::gsl_failure; 
    /* trivial case */
    if (dim == 0) {
      // AWS: this is one location where vector types need sync'ing
      if (f(0, 1, &xmin[0], fdim, &(val[0]))) return o2scl::gsl_failure;
      for (i = 0; i < fdim; ++i) err[i] = 0;
      return o2scl::success;
    }

    for (i = 0; i < fdim; ++i) {
      val[i] = 0;
      err[i] = HUGE_VAL;
    }

    for (i = 0; i < dim; ++i) {
      /* scale factor for C-C volume */
      V *= (xmax[i] - xmin[i]) * 0.5; 
    }

    new_nbuf = num_cacheval(m,dim,dim,0);

    if (max_nbuf < 1) max_nbuf = 1;
    if (new_nbuf > max_nbuf) new_nbuf = max_nbuf;
    if (nbuf < new_nbuf) {
      buf.resize((nbuf=new_nbuf)*dim);
    }

    /* start by evaluating the m=0 cubature rule */
    if (add_cacheval(vc,m, dim, fdim, f, dim, xmin, xmax, 
		     buf, nbuf) != o2scl::success) {
      return ret;
    }
    while (1) {
      size_t mi;

      eval_integral(vc,m,fdim,dim,V,mi,val,err,val1);
	
      if (converged(fdim,val,err,reqAbsError, reqRelError, norm) ||
	  (numEval > maxEval && maxEval)) {
	ret = o2scl::success;
	return ret;
      }
      m[mi] += 1;
      if (m[mi] > clencurt_M) {
	/* FAILURE */
	return ret;
      }

      new_nbuf = num_cacheval(m, mi, dim,0);
      if (new_nbuf > nbuf && nbuf < max_nbuf) {
	nbuf = new_nbuf;
	if (nbuf > max_nbuf) nbuf = max_nbuf;
	buf.resize(nbuf*dim);
      }

      if (add_cacheval(vc,m, mi, fdim, f, 
		       dim, xmin, xmax, buf, nbuf) != o2scl::success) {
	/* FAILURE */
	return ret;
      }
      numEval += new_nbuf;
    }

    return ret;
  }

  /** \brief Desc
   */
  static const size_t DEFAULT_MAX_NBUF=(1U << 20);
    
  /** \brief Desc
   */
  int integ(size_t fdim, func_t &f, size_t dim,
	    const vec_t &xmin, const vec_t &xmax, size_t maxEval,
	    double reqAbsError, double reqRelError, error_norm norm,
	    vec_t &val, vec_t &err) {
      
    int ret;
    size_t nbuf = 0;
    std::vector<size_t> m(dim);
    std::vector<double> buf;

    /* max_nbuf > 0 to amortize function overhead */
    ret = integ_v_buf(fdim,f,dim,xmin,xmax,
		      maxEval,reqAbsError,reqRelError,norm,
		      m,buf,nbuf,16,val,err);

    return ret;
  }

  };

}

#endif
