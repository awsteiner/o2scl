/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015, Andrew W. Steiner
  
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
    \brief Desc
*/
#ifndef O2SCL_CUBATURE_H
#define O2SCL_CUBATURE_H

// For size_t
#include <cstdlib>
// For memcpy
#include <cstring>
// For DBL_MIN
#include <cfloat>
// For stderr and fprintf
#include <cstdio>

#include <cmath>

#include <o2scl/clencurt.h>

namespace o2scl {

  /** \brief Desc
   */
  class inte_cubature_base {
    
  public:
    
    static const int SUCCESS=0;
    static const int FAILURE=1;
    
    /** \brief Desc
     */
    int isqr(int x) {
      return x * x;
    }

    /* Different ways of measuring the absolute and relative error when
       we have multiple integrands, given a vector e of error estimates
       in the individual components of a vector v of integrands.  These
       are all equivalent when there is only a single integrand. */
    typedef enum {
      ERROR_INDIVIDUAL = 0, /* individual relerr criteria in each component */
      ERROR_PAIRED, /* paired L2 norms of errors in each component,
		       mainly for integrating vectors of complex numbers */
      ERROR_L2, /* abserr is L_2 norm |e|, and relerr is |e|/|v| */
      ERROR_L1, /* abserr is L_1 norm |e|, and relerr is |e|/|v| */
      ERROR_LINF /* abserr is L_\infty norm |e|, and relerr is |e|/|v| */
    } error_norm;
    
  };

  /** \brief Desc

      This class is experimental.
   */
  template<class func_t, class func_v_t> class inte_hcubature
    : public inte_cubature_base {
    
  protected:

    /** \brief Desc
     */
    typedef struct {
      double val, err;
    } esterr;

    /** \brief Desc
     */
    double errMax(unsigned fdim, const esterr *ee) {

      double errmax = 0;
      unsigned k;
      for (k = 0; k < fdim; ++k) {
	if (ee[k].err > errmax) errmax = ee[k].err;
      }
      return errmax;
    }

    /** \brief Desc
     */
    typedef struct {
      unsigned dim;
      double *data; /* length 2*dim = center followed by half-widths */
      double vol;   /* cache volume = product of widths */
    } hypercube;

    /** \brief Desc
     */
    double compute_vol(const hypercube *h) {
      unsigned i;
      double vol = 1;
      for (i = 0; i < h->dim; ++i)
	vol *= 2 * h->data[i + h->dim];
      return vol;
    }

    /** \brief Desc
     */
    hypercube make_hypercube(unsigned dim, const double *center,
			     const double *halfwidth) {

      unsigned i;
      hypercube h;
      h.dim = dim;
      h.data = (double *) malloc(sizeof(double) * dim * 2);
      h.vol = 0;
      if (h.data) {
	for (i = 0; i < dim; ++i) {
	  h.data[i] = center[i];
	  h.data[i + dim] = halfwidth[i];
	}
	h.vol = compute_vol(&h);
      }
      return h;
    }

    /** \brief Desc
     */
    hypercube make_hypercube_range
      (unsigned dim, const double *xmin, const double *xmax) {

      hypercube h = make_hypercube(dim, xmin, xmax);
      unsigned i;
      if (h.data) {
	for (i = 0; i < dim; ++i) {
	  h.data[i] = 0.5 * (xmin[i] + xmax[i]);
	  h.data[i + dim] = 0.5 * (xmax[i] - xmin[i]);
	}
	h.vol = compute_vol(&h);
      }
      return h;
    }

    /** \brief Desc
     */
    void destroy_hypercube(hypercube *h) {
      free(h->data);
      h->dim = 0;
      return;
    }

    /** \brief Desc
     */
    typedef struct {
      hypercube h;
      unsigned splitDim;
      unsigned fdim; /* dimensionality of vector integrand */
      esterr *ee; /* array of length fdim */
      double errmax; /* max ee[k].err */
    } region;

    /** \brief Desc
     */
    region make_region(const hypercube *h, unsigned fdim) {

      region R;
      R.h = make_hypercube(h->dim, h->data, h->data + h->dim);
      R.splitDim = 0;
      R.fdim = fdim;
      R.ee = R.h.data ? (esterr *) malloc(sizeof(esterr) * fdim) : NULL;
      R.errmax = HUGE_VAL;

      return R;
    }

    /** \brief Desc
     */
    void destroy_region(region *R) {
      destroy_hypercube(&R->h);
      free(R->ee);
      R->ee = 0;
      return;
    }

    /** \brief Desc
     */
    int cut_region(region *R, region *R2) {

      unsigned d = R->splitDim, dim = R->h.dim;
      *R2 = *R;
      R->h.data[d + dim] *= 0.5;
      R->h.vol *= 0.5;
      R2->h = make_hypercube(dim, R->h.data, R->h.data + dim);
      if (!R2->h.data) return FAILURE;
      R->h.data[d] -= R->h.data[d + dim];
      R2->h.data[d] += R->h.data[d + dim];
      R2->ee = (esterr *) malloc(sizeof(esterr) * R2->fdim);
      return R2->ee == NULL;
    }

    struct rule_s; /* forward declaration */

    /** \brief Desc
     */
    typedef int (*evalError_func)(struct rule_s *r,
				  unsigned fdim, func_v_t &f, void *fdata,
				  unsigned nR, region *R);

    /** \brief Desc
     */
    typedef void (*destroy_func)(struct rule_s *r);

    /** \brief Desc
     */
    typedef struct rule_s {
      /* the dimensionality & number of functions */
      unsigned dim, fdim;
      /* number of evaluation points */
      unsigned num_points;
      /* max number of regions evaluated at once */
      unsigned num_regions;
      /* points to eval: num_regions * num_points * dim */
      double *pts;
      /* num_regions * num_points * fdim */
      double *vals; 
      evalError_func evalError;
      destroy_func destroy;
    } rule;

    /** \brief Desc
     */
    void destroy_rule(rule *r) {
      if (r) {
	if (r->destroy) r->destroy(r);
	free(r->pts);
	free(r);
      }
      return;
    }

    /** \brief Desc
     */
    static int alloc_rule_pts(rule *r, unsigned num_regions) {
      if (num_regions > r->num_regions) {
	free(r->pts);
	r->pts = r->vals = NULL;
	r->num_regions = 0;

	/* allocate extra so that repeatedly calling alloc_rule_pts
	   with growing num_regions only needs a logarithmic number of
	   allocations 
	*/
	num_regions *= 2; 

	r->pts = (double *) malloc(sizeof(double) * 
				   (num_regions
				    * r->num_points * (r->dim + r->fdim)));
	if (r->fdim + r->dim > 0 && !r->pts) return FAILURE;
	r->vals = r->pts + num_regions * r->num_points * r->dim;
	r->num_regions = num_regions;
      }
      return SUCCESS;
    }

    /** \brief Desc
     */
    rule *make_rule(size_t sz, /* >= sizeof(rule) */
		    unsigned dim, unsigned fdim, unsigned num_points,
		    evalError_func evalError, destroy_func destroy) {

      rule *r;

      if (sz < sizeof(rule)) return NULL;
      r = (rule *) malloc(sz);
      if (!r) return NULL;
      r->pts = r->vals = NULL;
      r->num_regions = 0;
      r->dim = dim; r->fdim = fdim; r->num_points = num_points;
      r->evalError = evalError;
      r->destroy = destroy;
      return r;
    }

    /** \brief Desc

	note: all regions must have same fdim 
     */
    int eval_regions(unsigned nR, region *R, 
		     func_v_t &f, void *fdata, rule *r) {

      unsigned iR;
      if (nR == 0) {
	/* nothing to evaluate */
	return SUCCESS;
      }
      if (r->evalError(r, R->fdim, f, fdata, nR, R)) return FAILURE;
      for (iR = 0; iR < nR; ++iR) {
	R[iR].errmax = errMax(R->fdim, R[iR].ee);
      }
      return SUCCESS;
    }

    /** \brief Desc

	Functions to loop over points in a hypercube.
    
	Based on orbitrule.cpp in HIntLib-0.0.10
    
	ls0 returns the least-significant 0 bit of n (e.g. it returns
	0 if the LSB is 0, it returns 1 if the 2 LSBs are 01, etcetera).
    */
    static unsigned ls0(unsigned n) {

#if defined(__GNUC__) &&                                        \
  ((__GNUC__ == 3 && __GNUC_MINOR__ >= 4) || __GNUC__ > 3)
      return __builtin_ctz(~n); /* gcc builtin for version >= 3.4 */
#else
      const unsigned bits[256] = {
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
      unsigned bit = 0;
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
    static void evalR_Rfs(double *pts, unsigned dim, double *p,
			  const double *c, const double *r) {
      
      unsigned i;
      unsigned signs = 0; /* 0/1 bit = +/- for corresponding element of r[] */

      /* We start with the point where r is ADDed in every coordinate
	 (this implies signs=0). */
      for (i = 0; i < dim; ++i)
	p[i] = c[i] + r[i];

      /* Loop through the points in Gray-code ordering */
      for (i = 0;; ++i) {
	unsigned mask, d;

	memcpy(pts, p, sizeof(double) * dim); pts += dim;

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
    static void evalRR0_0fs(double *pts, unsigned dim, double *p,
			    const double *c, const double *r) {
      
      unsigned i, j;
      
      for (i = 0; i < dim - 1; ++i) {
	p[i] = c[i] - r[i];
	for (j = i + 1; j < dim; ++j) {
	  p[j] = c[j] - r[j];
	  memcpy(pts, p, sizeof(double) * dim); pts += dim;
	  p[i] = c[i] + r[i];
	  memcpy(pts, p, sizeof(double) * dim); pts += dim;
	  p[j] = c[j] + r[j];
	  memcpy(pts, p, sizeof(double) * dim); pts += dim;
	  p[i] = c[i] - r[i];
	  memcpy(pts, p, sizeof(double) * dim); pts += dim;

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
    static void evalR0_0fs4d
      (double *pts, unsigned dim, double *p, const double *c,
       const double *r1, const double *r2) {
      
      unsigned i;

      memcpy(pts, p, sizeof(double) * dim); pts += dim;

      for (i = 0; i < dim; i++) {
	p[i] = c[i] - r1[i];
	memcpy(pts, p, sizeof(double) * dim); pts += dim;

	p[i] = c[i] + r1[i];
	memcpy(pts, p, sizeof(double) * dim); pts += dim;

	p[i] = c[i] - r2[i];
	memcpy(pts, p, sizeof(double) * dim); pts += dim;

	p[i] = c[i] + r2[i];
	memcpy(pts, p, sizeof(double) * dim); pts += dim;

	p[i] = c[i];
      }
      return;
    }

    /** \brief Desc
     */
    static size_t num0_0(size_t dim) { return 1; }
    /** \brief Desc
     */
    static size_t numR0_0fs(size_t dim) { return 2*dim; }
    /** \brief Desc
     */
    static size_t numRR0_0fs(size_t dim) { return 2*dim*(dim-1); }
    /** \brief Desc
     */
    static size_t numR_Rfs(size_t dim) { return 1U << dim; }
      
    /** \brief Desc

	Based on rule75genzmalik.cpp in HIntLib-0.0.10: An embedded
	cubature rule of degree 7 (embedded rule degree 5) due to A. C. Genz
	and A. A. Malik.  See:
	
	A. C. Genz and A. A. Malik, "An imbedded [sic] family of fully
	symmetric numerical integration rules," SIAM
	J. Numer. Anal. 20 (3), 580-588 (1983).
    */
    typedef struct {
    rule parent;

    /* temporary arrays of length dim */
    double *widthLambda, *widthLambda2, *p;

    /* dimension-dependent constants */
    double weight1, weight3, weight5;
    double weightE1, weightE3;
  } rule75genzmalik;
    
    /** \brief Desc
     */
    static double real(int x) {
    return ((double)(x));
  }
    
    /** \brief Desc
     */
    static int to_int(double n) {
    return ((int)(n));
  }

    /** \brief Desc
     */
    static void destroy_rule75genzmalik(rule *r_)
    {
    rule75genzmalik *r = (rule75genzmalik *) r_;
    free(r->p);
    return;
  }
    
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
    /** \brief Desc
     */
    static int rule75genzmalik_evalError
      (rule *r_, unsigned fdim, func_v_t &f, void *fdata,
      unsigned nR, region *R) {
    
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
    unsigned i, j, iR, dim = r_->dim;
    size_t npts = 0;
    double *diff, *pts, *vals;

    if (alloc_rule_pts(r_, nR)) return FAILURE;
    pts = r_->pts; vals = r_->vals;

    for (iR = 0; iR < nR; ++iR) {
    const double *center = R[iR].h.data;
    const double *halfwidth = R[iR].h.data + dim;
          
    for (i = 0; i < dim; ++i)
      r->p[i] = center[i];
          
    for (i = 0; i < dim; ++i)
      r->widthLambda2[i] = halfwidth[i] * lambda2;
    for (i = 0; i < dim; ++i)
      r->widthLambda[i] = halfwidth[i] * lambda4;

    /* Evaluate points in the center, in (lambda2,0,...,0) and
       (lambda3=lambda4, 0,...,0).  */
    evalR0_0fs4d(pts + npts*dim, dim, r->p, center, 
      r->widthLambda2, r->widthLambda);
    npts += num0_0(dim) + 2 * numR0_0fs(dim);


    /* Calculate points for (lambda4, lambda4, 0, ...,0) */
    evalRR0_0fs(pts + npts*dim, dim, r->p, center, r->widthLambda);
    npts += numRR0_0fs(dim);

    /* Calculate points for (lambda5, lambda5, ..., lambda5) */
    for (i = 0; i < dim; ++i)
      r->widthLambda[i] = halfwidth[i] * lambda5;
    evalR_Rfs(pts + npts*dim, dim, r->p, center, r->widthLambda);
    npts += numR_Rfs(dim);
  }

    /* Evaluate the integrand function(s) at all the points */
    if (f(dim, npts, pts, fdata, fdim, vals))
      return FAILURE;

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
    unsigned k, k0 = 0;
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
      unsigned dimDiffMax = 0;
  
      for (i = 0; i < dim; ++i) {
	if (diff[iR*dim + i] > maxdiff) {
	  maxdiff = diff[iR*dim + i];
	  dimDiffMax = i;
	}
      }
      R[iR].splitDim = dimDiffMax;
    }
return SUCCESS;
}

#ifdef O2SCL_NEVER_DEFINED
}{
#endif
    
    /** \brief Desc
     */
    rule *make_rule75genzmalik(unsigned dim, unsigned fdim)
    {
      rule75genzmalik *r;
      
      if (dim < 2) return NULL; /* this rule does not support 1d integrals */
      
      /* Because of the use of a bit-field in evalR_Rfs, we are limited
	 to be < 32 dimensions (or however many bits are in unsigned).
	 This is not a practical limitation...long before you reach
	 32 dimensions, the Genz-Malik cubature becomes excruciatingly
	 slow and is superseded by other methods (e.g. Monte-Carlo). */
      if (dim >= sizeof(unsigned) * 8) return NULL;
      
      r = (rule75genzmalik *) make_rule(sizeof(rule75genzmalik),
					dim, fdim,
					num0_0(dim) + 2 * numR0_0fs(dim)
					+ numRR0_0fs(dim) + numR_Rfs(dim),
					rule75genzmalik_evalError,
					destroy_rule75genzmalik);
      if (!r) return NULL;

      r->weight1 = (real(12824 - 9120 * to_int(dim) + 400 * isqr(to_int(dim)))
		    / real(19683));
      r->weight3 = real(1820 - 400 * to_int(dim)) / real(19683);
      r->weight5 = real(6859) / real(19683) / real(1U << dim);
      r->weightE1 = (real(729 - 950 * to_int(dim) + 50 * isqr(to_int(dim)))
		     / real(729));
      r->weightE3 = real(265 - 100 * to_int(dim)) / real(1458);

      r->p = (double *) malloc(sizeof(double) * dim * 3);
      if (!r->p) { destroy_rule((rule *) r); return NULL; }
      r->widthLambda = r->p + dim;
      r->widthLambda2 = r->p + 2 * dim;

      return (rule *) r;
    }

    /** \brief 1d 15-point Gaussian quadrature rule, based on qk15.c
	and qk.c in GNU GSL (which in turn is based on QUADPACK).
    */
    static int rule15gauss_evalError(rule *r,
				     unsigned fdim, func_v_t &f, void *fdata,
				     unsigned nR, region *R)
    {
      /* Gauss quadrature weights and kronrod quadrature abscissae and
	 weights as evaluated with 80 decimal digit arithmetic by
	 L. W. Fullerton, Bell Labs, Nov. 1981. */
      const unsigned n = 8;
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
      unsigned j, k, iR;
      size_t npts = 0;
      double *pts, *vals;

      if (alloc_rule_pts(r, nR)) return FAILURE;
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

      if (f(1, npts, pts, fdata, fdim, vals)) {
	return FAILURE;
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
	  if (result_abs > DBL_MIN / (50 * DBL_EPSILON)) {
	    double min_err = 50 * DBL_EPSILON * result_abs;
	    if (min_err > err) err = min_err;
	  }
	  R[iR].ee[k].err = err;
               
	  /* increment vk to point to next batch of results */
	  vk += 15*fdim;
	}
      }
      return SUCCESS;
    }
     
    /** \brief Desc
     */
    rule *make_rule15gauss(unsigned dim, unsigned fdim) {

      if (dim != 1) return NULL; /* this rule is only for 1d integrals */
       
      return make_rule(sizeof(rule),dim,fdim,15,
		       rule15gauss_evalError,0);
    }
     
    /* binary heap implementation (ala _Introduction to Algorithms_ by
       Cormen, Leiserson, and Rivest), for use as a priority queue of
       regions to integrate. */
     
    /** \brief Desc
     */
    typedef region heap_item;

    /** \brief Desc
     */
    typedef struct {
      size_t n, nalloc;
      heap_item *items;
      unsigned fdim;
      esterr *ee; /* array of length fdim of the total integrand & error */
    } heap;

    /** \brief Desc
     */
    void heap_resize(heap *h, size_t nalloc) {

      h->nalloc = nalloc;
      if (nalloc)
	h->items = (heap_item *) realloc(h->items, sizeof(heap_item)*nalloc);
      else {
	/* BSD realloc does not free for a zero-sized reallocation */
	free(h->items);
	h->items = NULL;
      }
    }

    /** \brief Desc
     */
    heap heap_alloc(size_t nalloc, unsigned fdim) {

      heap h;
      unsigned i;
      h.n = 0;
      h.nalloc = 0;
      h.items = 0;
      h.fdim = fdim;
      h.ee = (esterr *) malloc(sizeof(esterr) * fdim);
      if (h.ee) {
	for (i = 0; i < fdim; ++i) h.ee[i].val = h.ee[i].err = 0;
	heap_resize(&h, nalloc);
      }
      return h;
    }

    /** \brief Note that heap_free does not deallocate anything referenced by
	the items */
    void heap_free(heap *h) {

      h->n = 0;
      heap_resize(h, 0);
      h->fdim = 0;
      free(h->ee);
    }

    /** \brief Desc
     */
    int heap_push(heap *h, heap_item hi) {

      int insert;
      unsigned i, fdim = h->fdim;

      for (i = 0; i < fdim; ++i) {
	h->ee[i].val += hi.ee[i].val;
	h->ee[i].err += hi.ee[i].err;
      }
      insert = h->n;
      if (++(h->n) > h->nalloc) {
	heap_resize(h, h->n * 2);
	if (!h->items) return FAILURE;
      }

      while (insert) {
	int parent = (insert - 1) / 2;
	if (hi.errmax <= h->items[parent].errmax) {
	  break;
	}
	h->items[insert] = h->items[parent];
	insert = parent;
      }
      h->items[insert] = hi;
      return SUCCESS;
    }

    /** \brief Desc
     */
    int heap_push_many(heap *h, size_t ni, heap_item *hi) {
      size_t i;
      for (i = 0; i < ni; ++i)
	if (heap_push(h, hi[i])) return FAILURE;
      return SUCCESS;
    }

    /** \brief Desc
     */
    heap_item heap_pop(heap *h) {

      heap_item ret;
      int i, n, child;

      if (!(h->n)) {
	fprintf(stderr, "attempted to pop an empty heap\n");
	exit(EXIT_FAILURE);
      }

      ret = h->items[0];
      h->items[i = 0] = h->items[n = --(h->n)];
      while ((child = i * 2 + 1) < n) {
	int largest;
	heap_item swap;

	if (h->items[child].errmax <= h->items[i].errmax) {
	  largest = i;
	} else {
	  largest = child;
	}
	
	if (++child < n && h->items[largest].errmax <
	    h->items[child].errmax) {
	  largest = child;
	}
	if (largest == i)
	  break;
	swap = h->items[i];
	h->items[i] = h->items[largest];
	h->items[i = largest] = swap;
      }

      {
	unsigned i, fdim = h->fdim;
	for (i = 0; i < fdim; ++i) {
	  h->ee[i].val -= ret.ee[i].val;
	  h->ee[i].err -= ret.ee[i].err;
	}
      }
      return ret;
    }

    /** \brief Desc
     */
    int converged(unsigned fdim, const esterr *ee,
		  double reqAbsError, double reqRelError,
		  error_norm norm) {

      unsigned j;

      switch (norm) {
      case ERROR_INDIVIDUAL:
	for (j = 0; j < fdim; ++j)
	  if (ee[j].err > reqAbsError && ee[j].err >
	      fabs(ee[j].val)*reqRelError)
	    return 0;
	return 1;
              
      case ERROR_PAIRED:

	for (j = 0; j+1 < fdim; j += 2) {
	  double maxerr, serr, err, maxval, sval, val;
	  /* scale to avoid overflow/underflow */
	  maxerr = ee[j].err > ee[j+1].err ? ee[j].err : ee[j+1].err;
	  maxval = ee[j].val > ee[j+1].val ? ee[j].val : ee[j+1].val;
	  serr = maxerr > 0 ? 1/maxerr : 1;
	  sval = maxval > 0 ? 1/maxval : 1;
	  err = sqrt(isqr(ee[j].err*serr) + isqr(ee[j+1].err*serr)) * maxerr;
	  val = sqrt(isqr(ee[j].val*sval) + isqr(ee[j+1].val*sval)) * maxval;
	  if (err > reqAbsError && err > val*reqRelError)
	    return 0;
	}

	if (j < fdim) /* fdim is odd, do last dimension individually */
	  if (ee[j].err > reqAbsError && ee[j].err >
	      fabs(ee[j].val)*reqRelError)
	    return 0;
	return 1;

      case ERROR_L1: {
	double err = 0, val = 0;
	for (j = 0; j < fdim; ++j) {
	  err += ee[j].err;
	  val += fabs(ee[j].val);
	}
	return err <= reqAbsError || err <= val*reqRelError;
      }

      case ERROR_LINF: {
	double err = 0, val = 0;
	for (j = 0; j < fdim; ++j) {
	  double absval = fabs(ee[j].val);
	  if (ee[j].err > err) err = ee[j].err;
	  if (absval > val) val = absval;
	}
	return err <= reqAbsError || err <= val*reqRelError;
      }

      case ERROR_L2: {
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
	  err += isqr(ee[j].err * serr);
	  val += isqr(fabs(ee[j].val) * sval);
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
    int rulecubature(rule *r, unsigned fdim, func_v_t &f, void *fdata, 
		     const hypercube *h, size_t maxEval,
		     double reqAbsError, double reqRelError,
		     error_norm norm, double *val, double *err,
		     int parallel) {

      size_t numEval = 0;
      heap regions;
      unsigned i, j;
      region *R = NULL; /* array of regions to evaluate */
      size_t nR_alloc = 0;
      esterr *ee = NULL;

      if (fdim <= 1) norm = ERROR_INDIVIDUAL; /* norm is irrelevant */
      if (norm < 0 || norm > ERROR_LINF) return FAILURE; /* invalid norm */

      regions = heap_alloc(1, fdim);
      if (!regions.ee || !regions.items) goto bad;

      ee = (esterr *) malloc(sizeof(esterr) * fdim);
      if (!ee) goto bad;
     
      nR_alloc = 2;
      R = (region *) malloc(sizeof(region) * nR_alloc);
      if (!R) goto bad;
      R[0] = make_region(h, fdim);
      if (!R[0].ee || eval_regions(1, R, f, fdata, r) ||
	  heap_push(&regions, R[0])) {
	goto bad;
      }
      numEval += r->num_points;
     
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
	      if (!R) goto bad;
	    }
	    R[nR] = heap_pop(&regions);
	    for (j = 0; j < fdim; ++j) ee[j].err -= R[nR].ee[j].err;
	    if (cut_region(R+nR, R+nR+1)) goto bad;
	    numEval += r->num_points * 2;
	    nR += 2;
	    if (converged(fdim, ee, reqAbsError, reqRelError, norm))
	      break; /* other regions have small errs */
	  } while (regions.n > 0 && (numEval < maxEval || !maxEval));
	  if (eval_regions(nR, R, f, fdata, r)
	      || heap_push_many(&regions, nR, R))
	    goto bad;
	}
	else { /* minimize number of function evaluations */
	  R[0] = heap_pop(&regions); /* get worst region */
	  if (cut_region(R, R+1)
	      || eval_regions(2, R, f, fdata, r)
	      || heap_push_many(&regions, 2, R))
	    goto bad;
	  numEval += r->num_points * 2;
	}
      }

      /* re-sum integral and errors */
      for (j = 0; j < fdim; ++j) val[j] = err[j] = 0;  
      for (i = 0; i < regions.n; ++i) {
	for (j = 0; j < fdim; ++j) { 
	  val[j] += regions.items[i].ee[j].val;
	  err[j] += regions.items[i].ee[j].err;
	}
	destroy_region(&regions.items[i]);
      }

      /* printf("regions.nalloc = %d\n", regions.nalloc); */
      free(ee);
      heap_free(&regions);
      free(R);

      return SUCCESS;

    bad:
      
      free(ee);
      heap_free(&regions);
      free(R);

      return FAILURE;
    }
    
    /** \brief Desc
     */
    int cubature(unsigned fdim, func_v_t f, void *fdata, 
		 unsigned dim, const double *xmin, const double *xmax, 
		 size_t maxEval, double reqAbsError, double reqRelError, 
		 error_norm norm,
		 double *val, double *err, int parallel) {

      rule *r;
      hypercube h;
      int status;
      unsigned i;
      
      if (fdim == 0) /* nothing to do */ return SUCCESS;
      if (dim == 0) { /* trivial integration */
	if (f(0, 1, xmin, fdata, fdim, val)) return FAILURE;
	for (i = 0; i < fdim; ++i) err[i] = 0;
	return SUCCESS;
      }
      r = dim == 1 ? make_rule15gauss(dim, fdim)
	: make_rule75genzmalik(dim, fdim);
      if (!r) { 
	for (i = 0; i < fdim; ++i) {
	  val[i] = 0;
	  err[i] = HUGE_VAL; 
	}
	return FAILURE;
      }
      h = make_hypercube_range(dim, xmin, xmax);
      status = !h.data ? FAILURE
	: rulecubature(r, fdim, f, fdata, &h,
		       maxEval, reqAbsError, reqRelError, norm,
		       val, err, parallel);
      destroy_hypercube(&h);
      destroy_rule(r);
      return status;
    }
    
  public:
    
    /** \brief Desc
     */
    int integ_v(unsigned fdim, func_t f, void *fdata,
		unsigned dim, const double *xmin, const double *xmax, 
		size_t maxEval, double reqAbsError, double reqRelError, 
		error_norm norm, double *val, double *err) {
      return cubature(fdim,f,fdata,dim,xmin,xmax,
		      maxEval,reqAbsError,reqRelError,norm,val,err,1);
    }
    
    /** \brief Desc
     */
    typedef struct fv_data_s { func_t f; void *fdata; } fv_data;
    
    /** \brief Desc
     */
    static int fv(unsigned ndim, size_t npt,
		  const double *x, void *d_,
		  unsigned fdim, double *fval) {
      fv_data *d = (fv_data *) d_;
      func_t f = d->f;
      void *fdata = d->fdata;
      unsigned i;
      /* printf("npt = %u\n", npt); */
      for (i = 0; i < npt; ++i) {
	if (f(ndim, x + i*ndim, fdata, fdim, fval + i*fdim)) {
	  return FAILURE;
	}
      }
      return SUCCESS;
    }

    /** \brief Desc
     */
    int integ(unsigned fdim, func_t f, void *fdata,
	      unsigned dim, const double *xmin,
	      const double *xmax, size_t maxEval, double reqAbsError,
	      double reqRelError, error_norm norm, double *val,
	      double *err) {

      int ret;
      fv_data d;

      if (fdim == 0) return SUCCESS; /* nothing to do */     
     
      d.f = f;
      d.fdata = fdata;
      ret = cubature(fdim,fv,&d,dim,xmin,xmax,
		     maxEval,reqAbsError,reqRelError,norm,val,err,0);
      return ret;

    }
    
  };

#ifdef O2SCL_NEVER_DEFINED
}{
#endif
    
  /** \brief Desc
      
      This class is experimental.
   */
  template<class func_t, class func_v_t> class inte_pcubature
    : public inte_cubature_base {
      
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
    typedef struct cacheval_s {
      unsigned m[MAXDIM];
      unsigned mi;
      double *val;
    } cacheval;

    /** \brief Desc  array of ncache cachevals c[i] 
     */
    typedef struct valcache_s {
      size_t ncache;
      cacheval *c;
    } valcache;

    /** \brief Desc
     */
    void free_cachevals(valcache *v) {
      if (!v) return;
      if (v->c) {
	size_t i;
	for (i = 0; i < v->ncache; ++i) {
	  free(v->c[i].val);
	}
	free(v->c);
	v->c = NULL;
      }
      v->ncache = 0;
      return;
    }

    /** \brief Desc

	recursive loop over all cubature points for the given (m,mi)
	cache entry: add each point to the buffer buf, evaluating all
	at once whenever the buffer is full or when we are done
    */
    int compute_cacheval(const unsigned *m, unsigned mi, 
			 double *val, size_t *vali,
			 unsigned fdim, func_v_t f, void *fdata,
			 unsigned dim, unsigned id, double *p,
			 const double *xmin, const double *xmax,
			 double *buf, size_t nbuf, size_t *ibuf) {

      if (id == dim) { /* add point to buffer of points */
	memcpy(buf + (*ibuf)++ * dim, p, sizeof(double) * dim);
	if (*ibuf == nbuf) { /* flush buffer */
	  if (f(dim, nbuf, buf, fdata, fdim, val + *vali)) {
	    return FAILURE;
	  }
	  *vali += *ibuf * fdim;
	  *ibuf = 0;
	}

      } else {
      
	double c = (xmin[id] + xmax[id]) * 0.5;
	double r = (xmax[id] - xmin[id]) * 0.5;
	const double *x = clencurt_x 
	  + ((id == mi) ? (m[id] ? (1 << (m[id] - 1)) : 0) : 0);
	unsigned i, nx = (id == mi ? (m[id] ? (1 << (m[id] - 1)) : 1)
			  : (1 << (m[id])));
	if (id != mi) {
	  p[id] = c;
	  if (compute_cacheval(m, mi, val, vali, fdim, f, fdata,
			       dim, id + 1, p,
			       xmin, xmax, buf, nbuf, ibuf)) {
	    return FAILURE;
	  }
	}
	for (i = 0; i < nx; ++i) {
	  p[id] = c + r * x[i];
	  if (compute_cacheval(m, mi, val, vali, fdim, f, fdata,
			       dim, id + 1, p,
			       xmin, xmax, buf, nbuf, ibuf)) {
	    return FAILURE;
	  }
	  p[id] = c - r * x[i];
	  if (compute_cacheval(m, mi, val, vali, fdim, f, fdata,
			       dim, id + 1, p,
			       xmin, xmax, buf, nbuf, ibuf)) {
	    return FAILURE;
	  }
	}
      }
      return SUCCESS;
    }

    /** \brief Desc
     */
    size_t num_cacheval(const unsigned *m, unsigned mi, unsigned dim) {

      unsigned i;
      size_t nval = 1;
      for (i = 0; i < dim; ++i) {
	if (i == mi) {
	  nval *= m[i] == 0 ? 2 : (1 << (m[i]));
	} else {
	  nval *= (1 << (m[i] + 1)) + 1;
	}
      }
      return nval;
    }

    /** \brief Desc
     */
    int add_cacheval(valcache *vc, const unsigned *m, unsigned mi,
		     unsigned fdim, func_v_t f, void *fdata,
		     unsigned dim, const double *xmin,
		     const double *xmax, double *buf, size_t nbuf) {
      
      size_t ic = vc->ncache;
      size_t nval, vali = 0, ibuf = 0;
      double p[MAXDIM];

      vc->c = (cacheval *) realloc(vc->c, sizeof(cacheval) * ++(vc->ncache));
      if (!vc->c) return -1;

      vc->c[ic].mi = mi;
      memcpy(vc->c[ic].m, m, sizeof(unsigned) * dim);
      nval = fdim * num_cacheval(m, mi, dim);
      vc->c[ic].val = (double *) malloc(sizeof(double) * nval);
      if (!vc->c[ic].val) return FAILURE;

      if (compute_cacheval(m, mi, vc->c[ic].val, &vali,
			   fdim, f, fdata,
			   dim, 0, p, xmin, xmax,
			   buf, nbuf, &ibuf)) {
	return FAILURE;
      }

      if (ibuf > 0) {
	/* flush remaining buffer */
	return f(dim, ibuf, buf, fdata, fdim, vc->c[ic].val + vali);
      }

      return SUCCESS;
    }
    
    /** \brief Desc
	
	recursive loop to evaluate the integral contribution from the
	cache entry c, accumulating in val, for the given m[] except
	with m[md] -> m[md] - 1 if md < dim, using the cached values
	(cm,cmi,cval). id is the current loop dimension (from 0 to
	dim-1).
    */
    unsigned eval(const unsigned *cm, unsigned cmi, double *cval,
		  const unsigned *m, unsigned md,
		  unsigned fdim, unsigned dim, unsigned id,
		  double weight, double *val) {

      size_t voff = 0; /* amount caller should offset cval array afterwards */
      if (id == dim) {
	unsigned i;
	for (i = 0; i < fdim; ++i) val[i] += cval[i] * weight;
	voff = fdim;

      } else if (m[id] == 0 && id == md) {

	/* using trivial rule for this dim */
	voff = eval(cm, cmi, cval, m, md, fdim, dim, id+1, weight*2, val);
	voff += fdim * (1 << cm[id]) * 2
	  * num_cacheval(cm + id+1, cmi - (id+1), dim - (id+1));

      } else {
      
	unsigned i;
	unsigned mid = m[id] - (id == md); /* order of C-C rule */
	const double *w = clencurt_w + mid + (1 << mid) - 1
	  + (id == cmi ? (cm[id] ? 1 + (1 << (cm[id]-1)) : 1) : 0);
	unsigned cnx = (id == cmi ? (cm[id] ? (1 << (cm[id]-1)) : 1)
			: (1 << (cm[id])));
	unsigned nx = cm[id] <= mid ? cnx : (1 << mid);

	if (id != cmi) {
	  voff = eval(cm, cmi, cval, m, md, fdim, dim, id + 1,
		      weight * w[0], val);
	  ++w;
	}
	for (i = 0; i < nx; ++i) {
	  voff += eval(cm, cmi, cval + voff, m, md, fdim, dim, id + 1,
		       weight * w[i], val);
	  voff += eval(cm, cmi, cval + voff, m, md, fdim, dim, id + 1,
		       weight * w[i], val);
	}

	voff += (cnx - nx) * fdim * 2
	  * num_cacheval(cm + id+1, cmi - (id+1), dim - (id+1));
      }
      return voff;
    }

    /** \brief Desc

	loop over all cache entries that contribute to the integral,
	(with m[md] decremented by 1) 
    */
    void evals(valcache vc, const unsigned *m, unsigned md,
	       unsigned fdim, unsigned dim, double V, double *val) {

      size_t i;

      memset(val, 0, sizeof(double) * fdim);
      for (i = 0; i < vc.ncache; ++i) {
	if (vc.c[i].mi >= dim ||
	    vc.c[i].m[vc.c[i].mi] + (vc.c[i].mi == md) <= m[vc.c[i].mi]) {
	  eval(vc.c[i].m, vc.c[i].mi, vc.c[i].val,
	       m, md, fdim, dim, 0, V, val);
	}
      }
      return;
    }

    /** \brief Desc

	evaluate the integrals for the given m[] using the cached
	values in vc, storing the integrals in val[], the error
	estimate in err[], and the dimension to subdivide next (the
	largest error contribution) in *mi
    */
    void eval_integral(valcache vc, const unsigned *m, 
		       unsigned fdim, unsigned dim, double V,
		       unsigned *mi, double *val, double *err,
		       double *val1) {

      double maxerr = 0;
      unsigned i, j;
     
      evals(vc, m, dim, fdim, dim, V, val);

      /* error estimates along each dimension by comparing val with
	 lower-order rule in that dimension; overall (conservative)
	 error estimate from maximum error of lower-order rules. */
      memset(err, 0, sizeof(double) * fdim);
      *mi = 0;
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
	  *mi = i;
	}
      }
      /* printf("eval: %g +/- %g (dim %u)\n", val[0], err[0], *mi); */
      return;
    }

    /** \brief Desc
     */
    int converged(unsigned fdim, const double *vals, const double *errs,
		  double reqAbsError, double reqRelError,
		  error_norm norm) {

      unsigned j;
      switch (norm) {
      case ERROR_INDIVIDUAL:
	for (j = 0; j < fdim; ++j)
	  if (errs[j] > reqAbsError && errs[j] > fabs(vals[j])*reqRelError)
	    return 0;
	return 1;
	      
      case ERROR_PAIRED:
	for (j = 0; j+1 < fdim; j += 2) {
	  double maxerr, serr, err, maxval, sval, val;
	  /* scale to avoid overflow/underflow */
	  maxerr = errs[j] > errs[j+1] ? errs[j] : errs[j+1];
	  maxval = vals[j] > vals[j+1] ? vals[j] : vals[j+1];
	  serr = maxerr > 0 ? 1/maxerr : 1;
	  sval = maxval > 0 ? 1/maxval : 1;
	  err = sqrt(isqr(errs[j]*serr) + isqr(errs[j+1]*serr)) * maxerr;
	  val = sqrt(isqr(vals[j]*sval) + isqr(vals[j+1]*sval)) * maxval;
	  if (err > reqAbsError && err > val*reqRelError)
	    return 0;
	}
	if (j < fdim) /* fdim is odd, do last dimension individually */
	  if (errs[j] > reqAbsError && errs[j] > fabs(vals[j])*reqRelError)
	    return 0;
	return 1;

      case ERROR_L1: {
	double err = 0, val = 0;
	for (j = 0; j < fdim; ++j) {
	  err += errs[j];
	  val += fabs(vals[j]);
	}
	return err <= reqAbsError || err <= val*reqRelError;
      }

      case ERROR_LINF: {
	double err = 0, val = 0;
	for (j = 0; j < fdim; ++j) {
	  double absval = fabs(vals[j]);
	  if (errs[j] > err) err = errs[j];
	  if (absval > val) val = absval;
	}
	return err <= reqAbsError || err <= val*reqRelError;
      }

      case ERROR_L2: {
	double maxerr = 0, maxval = 0, serr, sval, err = 0, val = 0;
	/* scale values by 1/max to avoid overflow/underflow */
	for (j = 0; j < fdim; ++j) {
	  double absval = fabs(vals[j]);
	  if (errs[j] > maxerr) maxerr = errs[j];
	  if (absval > maxval) maxval = absval;
	}
	serr = maxerr > 0 ? 1/maxerr : 1;
	sval = maxval > 0 ? 1/maxval : 1;
	for (j = 0; j < fdim; ++j) {
	  err += isqr(errs[j] * serr);
	  val += isqr(fabs(vals[j]) * sval);
	}
	err = sqrt(err) * maxerr;
	val = sqrt(val) * maxval;
	return err <= reqAbsError || err <= val*reqRelError;
      }
      }
      return 1; /* unreachable */
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
    int integ_v_buf(unsigned fdim, func_v_t f, void *fdata,
		    unsigned dim, const double *xmin, const double *xmax,
		    size_t maxEval, double reqAbsError, double reqRelError,
		    error_norm norm, unsigned *m,
		    double **buf, size_t *nbuf, size_t max_nbuf,
		    double *val, double *err) {

      int ret = FAILURE;
      double V = 1;
      size_t numEval = 0, new_nbuf;
      unsigned i;
      valcache vc = {0, NULL};
      double *val1 = NULL;

      if (fdim <= 1) norm = ERROR_INDIVIDUAL; /* norm is irrelevant */
      if (norm < 0 || norm > ERROR_LINF) return FAILURE; /* invalid norm */

      if (fdim == 0) return SUCCESS; /* nothing to do */
      if (dim > MAXDIM) return FAILURE; /* unsupported */
      if (dim == 0) { /* trivial case */
	if (f(0, 1, xmin, fdata, fdim, val)) return FAILURE;
	for (i = 0; i < fdim; ++i) err[i] = 0;
	return SUCCESS;
      }

      for (i = 0; i < fdim; ++i) {
	val[i] = 0;
	err[i] = HUGE_VAL;
      }

      for (i = 0; i < dim; ++i) {
	V *= (xmax[i] - xmin[i]) * 0.5; /* scale factor for C-C volume */
      }

      new_nbuf = num_cacheval(m, dim, dim);

      if (max_nbuf < 1) max_nbuf = 1;
      if (new_nbuf > max_nbuf) new_nbuf = max_nbuf;
      if (*nbuf < new_nbuf) {
	free(*buf);
	*buf = (double *) malloc(sizeof(double) 
				 * (*nbuf = new_nbuf) * dim);
	if (!*buf) goto done;
      }

      /* start by evaluating the m=0 cubature rule */
      if (add_cacheval(&vc, m, dim, fdim, f, fdata, dim, xmin, xmax, 
		       *buf, *nbuf) != SUCCESS)
	goto done;

      val1 = (double *) malloc(sizeof(double) * fdim);

      while (1) {
	unsigned mi;

	eval_integral(vc, m, fdim, dim, V, &mi, val, err, val1);
	if (converged(fdim, val, err, reqAbsError, reqRelError, norm)
	    || (numEval > maxEval && maxEval)) {
	  ret = SUCCESS;
	  goto done;
	}
	m[mi] += 1;
	if (m[mi] > clencurt_M) goto done; /* FAILURE */

	new_nbuf = num_cacheval(m, mi, dim);
	if (new_nbuf > *nbuf && *nbuf < max_nbuf) {
	  *nbuf = new_nbuf;
	  if (*nbuf > max_nbuf) *nbuf = max_nbuf;
	  free(*buf);
	  *buf = (double *) malloc(sizeof(double) * *nbuf * dim);
	  if (!*buf) goto done; /* FAILURE */
	}

	if (add_cacheval(&vc, m, mi, fdim, f, fdata, 
			 dim, xmin, xmax, *buf, *nbuf) != SUCCESS)
	  goto done; /* FAILURE */
	numEval += new_nbuf;
      }

    done:

      free(val1);
      free_cachevals(&vc);
      
      return ret;
    }

    /** \brief Desc
     */
    static const size_t DEFAULT_MAX_NBUF=(1U << 20);
    
    /** \brief Desc
     */
    int integ_v(unsigned fdim, func_v_t f, void *fdata,
		unsigned dim, const double *xmin, const double *xmax,
		size_t maxEval, double reqAbsError, double reqRelError,
		error_norm norm, double *val, double *err) {

      int ret;
      size_t nbuf = 0;
      unsigned m[MAXDIM];
      double *buf = NULL;
      memset(m, 0, sizeof(unsigned) * dim);
      ret = integ_v_buf(fdim, f, fdata, dim, xmin, xmax,
			maxEval, reqAbsError, reqRelError, norm,
			m, &buf, &nbuf, DEFAULT_MAX_NBUF, val, err);
      free(buf);

      return ret;
    }

    /** \brief Desc
     */
    typedef struct fv_data_s { func_t f; void *fdata; } fv_data;

    /** \brief Desc
     */
    static int fv(unsigned ndim, size_t npt, const double *x, void *d_,
		  unsigned fdim, double *fval) {

      fv_data *d = (fv_data *) d_;
      func_t f = d->f;
      void *fdata = d->fdata;
      unsigned i;
      /* printf("npt = %u\n", npt); */
      for (i = 0; i < npt; ++i) 
	if (f(ndim, x + i*ndim, fdata, fdim, fval + i*fdim))
	  return FAILURE;
      return SUCCESS;
    }

    /** \brief Desc
     */
    int integ(unsigned fdim, func_t f, void *fdata,
	      unsigned dim, const double *xmin, const double *xmax,
	      size_t maxEval, double reqAbsError, double reqRelError,
	      error_norm norm, double *val, double *err) {
      
      int ret;
      size_t nbuf = 0;
      unsigned m[MAXDIM];
      double *buf = NULL;
      fv_data d;

      d.f = f;
      d.fdata = fdata;
      memset(m, 0, sizeof(unsigned) * dim);
      /* max_nbuf > 0 to amortize function overhead */
      ret = integ_v_buf(fdim, fv, &d, dim, xmin, xmax, 
			maxEval, reqAbsError, reqRelError, norm,
			m, &buf, &nbuf, 16, val, err);
      free(buf);
      return ret;
    }

  };

  /* USAGE: Call hcubature or pcubature with your function as described
     in the README file. */
  
  /* a vector integrand - evaluates the function at the given point x
     (an array of length ndim) and returns the result in fval (an
     array of length fdim). The void* parameter is there in case you
     have to pass any additional data through to your function (it
     corresponds to the fdata parameter you pass to cubature). Return
     0 on success or nonzero to terminate the integration. */
  typedef int (*integrand)(unsigned ndim, const double *x, void *,
			   unsigned fdim, double *fval);
  
  /* a vector integrand of a vector of npt points: x[i*ndim + j] is the
     j-th coordinate of the i-th point, and the k-th function evaluation
     for the i-th point is returned in fval[i*fdim + k].  Return 0 on success
     or nonzero to terminate the integration. */
  typedef int (*integrand_v)(unsigned ndim, size_t npt,
			     const double *x, void *,
			     unsigned fdim, double *fval);
  
  /* Different ways of measuring the absolute and relative error when
     we have multiple integrands, given a vector e of error estimates
     in the individual components of a vector v of integrands.  These
     are all equivalent when there is only a single integrand. */
  typedef enum {
    ERROR_INDIVIDUAL = 0, /* individual relerr criteria in each component */
    ERROR_PAIRED, /* paired L2 norms of errors in each component,
		     mainly for integrating vectors of complex numbers */
    ERROR_L2, /* abserr is L_2 norm |e|, and relerr is |e|/|v| */
    ERROR_L1, /* abserr is L_1 norm |e|, and relerr is |e|/|v| */
    ERROR_LINF /* abserr is L_\infty norm |e|, and relerr is |e|/|v| */
  } error_norm;

  /* Integrate the function f from xmin[dim] to xmax[dim], with at most
     maxEval function evaluations (0 for no limit), until the given
     absolute or relative error is achieved.  val returns the integral,
     and err returns the estimate for the absolute error in val; both
     of these are arrays of length fdim, the dimension of the vector
     integrand f(x). The return value of the function is 0 on success
     and non-zero if there  was an error. */
  
  /* adapative integration by partitioning the integration domain
     ("h-adaptive") and using the same fixed-degree quadrature in each
     subdomain, recursively, until convergence is achieved. */
  int hcubature(unsigned fdim, integrand f, void *fdata,
		unsigned dim, const double *xmin, const double *xmax, 
		size_t maxEval, double reqAbsError, double reqRelError, 
		error_norm norm,
		double *val, double *err);

  /* as hcubature, but vectorized integrand */
  int hcubature_v(unsigned fdim, integrand_v f, void *fdata,
		  unsigned dim, const double *xmin, const double *xmax, 
		  size_t maxEval, double reqAbsError, double reqRelError, 
		  error_norm norm,
		  double *val, double *err);

  /* adaptive integration by increasing the degree of (tensor-product
     Clenshaw-Curtis) quadrature rules ("p-adaptive"), rather than
     subdividing the domain ("h-adaptive").  Possibly better for
     smooth integrands in low dimensions. */
  int pcubature_v_buf(unsigned fdim, integrand_v f, void *fdata,
		      unsigned dim, const double *xmin, const double *xmax,
		      size_t maxEval, 
		      double reqAbsError, double reqRelError,
		      error_norm norm,
		      unsigned *m,
		      double **buf, size_t *nbuf, size_t max_nbuf,
		      double *val, double *err);
  
  int pcubature_v(unsigned fdim, integrand_v f, void *fdata,
		  unsigned dim, const double *xmin, const double *xmax, 
		  size_t maxEval, double reqAbsError, double reqRelError, 
		  error_norm norm,
		  double *val, double *err);

  int pcubature(unsigned fdim, integrand f, void *fdata,
		unsigned dim, const double *xmin, const double *xmax, 
		size_t maxEval, double reqAbsError, double reqRelError, 
		error_norm norm,
		double *val, double *err);

}

#endif
