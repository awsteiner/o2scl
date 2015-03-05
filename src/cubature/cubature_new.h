/*
  -------------------------------------------------------------------
  
  Copyright (C) 2015, Andrew W. Steiner
  
  This file is part of O2scl. It has been adapted from cubature_new 
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
/** \file cubature_new.h
    \brief File for definitions of \ref o2scl::inte_hcubature_new and 
    \ref o2scl::inte_pcubature_new
*/
#ifndef O2SCL_CUBATURE_NEW_H
#define O2SCL_CUBATURE_NEW_H

// For memcpy
#include <cstring>

#include <cmath>
#include <functional>
#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/clencurt.h>
#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>

namespace o2scl {

#ifndef O2SCL_NEVER_DEFINED
  
  /** \brief Desc
   */
  typedef std::function<
    int(size_t,const boost::numeric::ublas::vector<double> &,
	size_t,boost::numeric::ublas::vector<double> &) > cub_funct11;
  
  /** \brief Desc
   */
  typedef std::function<
    int(size_t,size_t,const boost::numeric::ublas::vector<double> &,
	size_t,boost::numeric::ublas::vector<double> &) > cub_vec_funct11;

  /** \brief Desc
   */
  template<class vec_t, class func_t> class cub_wrapper {

  protected:

    /** \brief Desc
     */
    func_t *fp;
    
  public:

    /** \brief Desc
     */
    cub_wrapper(func_t &f) {
      fp=&f;
    }
    
    /** \brief Desc
     */
    int operator()(size_t ndim, size_t npts, const vec_t &x,
		   size_t fdim, vec_t &fval) {
      for (size_t i = 0; i < npts; i++) {
	vec_t s=o2scl::vector_range<double>(x,i*ndim,(i+1)*ndim);
	vec_t sf=o2scl::vector_range<double>(fval,i*fdim,(i+1)*fdim);
	f(ndim,s,fdim,sf);
      }
      return o2scl::success;
    }
    
  };

#endif
  
  /** \brief Base class for integration routines from the 
      Cubature_New library
  */
  class inte_cubature_new_base {
    
  public:
    
    /** \brief Square an integer
     */
    int isqr(int x) {
      return x*x;
    }
    
    /** \brief Square an double-precision floating-point number
     */
    double dsqr(double x) {
      return x*x;
    }

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
      using cubature_new rules from the Cubature_New library

      This class is experimental.

      \hline
      \b Documentation \b adapted \b from \b Cubature_New
      
      A cubature_new rule takes a function and a hypercube and evaluates
      the function at a small number of points, returning an estimate
      of the integral as well as an estimate of the error, and also
      a suggested dimension of the hypercube to subdivide.
      
      Given such a rule, the adaptive integration is simple:
      
      1) Evaluate the cubature_new rule on the hypercube(s).
      Stop if converged.
      
      2) Pick the hypercube with the largest estimated error,
      and divide it in two along the suggested dimension.
      
      3) Goto (1).
      
      The basic algorithm is based on the adaptive cubature_new described
      in \ref Genz80 and subsequently extended to integrating a vector
      of integrands in \ref Berntsen91 .

      \hline      

  */
  template<class func_t, class vec_t> class inte_hcubature_new
    : public inte_cubature_new_base {
    
  protected:

    /** \brief Simple class for a value and its uncertainty
	\comment
	Default copy constructors are ok for this class
	\endcomment
     */
    class esterr {

    public:

      /** \brief Value */
      double val;
      /** \brief Error */
      double err;

    };
    
    /** \brief Return the maximum error from the array \c ee
     */
    double errMax(unsigned fdim, const std::vector<esterr> &ee) {
      
      double errmax = 0;
      for (unsigned k = 0; k < fdim; ++k) {
	if (ee[k].err > errmax) errmax = ee[k].err;
      }
      return errmax;
    }

    /** \brief Specification of the hypercubic region over which 
	one wants to integrate
	
	\comment
	Default copy constructors are not ok for this class
	\endcomment
     */
    class hypercube {
      
    public:

      /** \brief Desc */
      unsigned dim;
      /** \brief length 2*dim = center followed by half-widths */
      vec_t data; 
      /** \brief cache volume = product of widths */
      double vol;   
    };

    /** \brief Desc
     */
    double compute_vol(const hypercube &h) {
      double vol = 1.0;
      for (unsigned i = 0; i < h.dim; ++i) {
	vol *= 2 * h.data[i + h.dim];
      }
      return vol;
    }

    /** \brief Desc
     */
    hypercube make_hypercube(unsigned dim, const vec_t &center,
			     const vec_t &halfwidth) {
      
      hypercube h;
      h.dim = dim;
      h.data.resize(dim*2);
      h.vol = 0;
      if (h.data.size()>0) {
	for (unsigned i = 0; i < dim; ++i) {
	  h.data[i] = center[i];
	  h.data[i + dim] = halfwidth[i];
	}
	h.vol = compute_vol(h);
      }
      return h;
    }

    /** \brief Desc
     */
    hypercube make_hypercube_range
      (unsigned dim, const vec_t &xmin, const vec_t &xmax) {

      hypercube h = make_hypercube(dim, xmin, xmax);
      if (h.data.size()>0) {
	for (unsigned i = 0; i < dim; ++i) {
	  h.data[i] = 0.5 * (xmin[i] + xmax[i]);
	  h.data[i + dim] = 0.5 * (xmax[i] - xmin[i]);
	}
	h.vol = compute_vol(h);
      }
      return h;
    }

    /** \brief Desc
     */
    void destroy_hypercube(hypercube &h) {
      h.data.clear();
      h.dim = 0;
      return;
    }

    /** \brief Desc
	\comment
	Default copy constructors are not ok for this class
	\endcomment
     */
    class region {

    public:

      /** \brief Desc */
      hypercube h;
      /** \brief Desc */
      unsigned splitDim;
      /** \brief dimensionality of vector integrand */
      unsigned fdim; 
      /** \brief array of length fdim */
      std::vector<esterr> ee;
      /** \brief max ee[k].err */
      double errmax; 

    };

    /** \brief Desc
     */
    region make_region(const hypercube &h, unsigned fdim) {

      region R;
      vec_t htmp=o2scl::vector_range(h.data,h.dim,h.data.size());
      R.h = make_hypercube(h.dim, h.data, htmp);
      R.splitDim = 0;
      R.fdim = fdim;
      if (R.h.data.size()>0) {
	R.ee.resize(fdim);
      }
      R.errmax = HUGE_VAL;

      return R;
    }

    /** \brief Desc
     */
    void destroy_region(region &R) {
      destroy_hypercube(R.h);
      R.ee.clear();
      return;
    }

    /** \brief Desc
     */
    int cut_region(region *R, region *R2) {

      unsigned d = R->splitDim, dim = R->h.dim;
      *R2 = *R;
      R->h.data[d + dim] *= 0.5;
      R->h.vol *= 0.5;
      vec_t vtmp=o2scl::vector_range(R->h.data,dim,R->h.data.size());
      R2->h = make_hypercube(dim, R->h.data, vtmp);
      if (R2->h.data.size()==0) return o2scl::gsl_failure;
      R->h.data[d] -= R->h.data[d + dim];
      R2->h.data[d] += R->h.data[d + dim];
      R2->ee.resize(R2->fdim);
      return 0;
    }

    struct rule_s; /* forward declaration */

    /** \brief Desc
     */
    typedef int (*evalError_func)(struct rule_s *r, unsigned fdim,
				  func_t &f, unsigned nR, region *R);
    
    /** \brief Desc
     */
    typedef void (*destroy_func)(struct rule_s *r);

    /** \brief Desc
     */
    typedef struct rule_s {
      /** \brief The dimensionality and the number of functions */
      unsigned dim, fdim;
      /** \brief The number of evaluation points */
      unsigned num_points;
      /** \brief The max number of regions evaluated at once */
      unsigned num_regions;
      /** \brief points to eval: num_regions * num_points * dim */
      vec_t pts;
      /** \brief num_regions * num_points * fdim */
      vec_t vals;
      /** \brief Desc */
      evalError_func evalError;
      /** \brief Desc */
      destroy_func destroy;
    } rule;

    /** \brief Desc
     */
    void destroy_rule(rule *r) {
      if (r) {
	if (r->destroy) r->destroy(r);
	free(r);
      }
      return;
    }

    /** \brief Desc
     */
    static int alloc_rule_pts(rule *r, unsigned num_regions) {
      if (num_regions > r->num_regions) {
	r->pts.clear();
	r->vals.clear();
	r->num_regions = 0;

	/* Allocate extra so that repeatedly calling alloc_rule_pts
	   with growing num_regions only needs a logarithmic number of
	   allocations 
	*/
	num_regions *= 2; 

	r->pts.resize((num_regions* r->num_points * (r->dim + r->fdim)));
	r->vals=o2scl::vector_range(r->pts,num_regions * r->num_points *
				    r->dim,r->vals.size());
	r->num_regions = num_regions;
      }
      return o2scl::success;
    }

    /** \brief Desc
     */
    rule *make_rule(size_t sz, /* >= sizeof(rule) */
		    unsigned dim, unsigned fdim, unsigned num_points,
		    evalError_func evalError, destroy_func destroy) {
      
      rule *r;

      if (sz < sizeof(rule)) return 0;
      r = (rule *) malloc(sz);
      if (!r) return 0;
      r->pts.clear();
      r->vals.clear();
      r->num_regions = 0;
      r->dim = dim;
      r->fdim = fdim;
      r->num_points = num_points;
      r->evalError = evalError;
      r->destroy = destroy;
      return r;
    }

    /** \brief Desc

	\note All regions must have same fdim 
    */
    int eval_regions(unsigned nR, region *R, 
		     func_t &f, rule *r) {

      unsigned iR;
      if (nR == 0) {
	/* nothing to evaluate */
	return o2scl::success;
      }
      if (r->evalError(r, R->fdim, f, nR, R)) return o2scl::gsl_failure;
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
    static void evalR_Rfs(vec_t &pts, unsigned dim, double *p,
			  const vec_t &c, const double *r) {
      
      unsigned i;
      /* 0/1 bit = +/- for corresponding element of r[] */
      unsigned signs = 0; 

      size_t istart=0;
      vec_t vtmp=o2scl::vector_range(pts,istart,pts.size());

      /* We start with the point where r is ADDed in every coordinate
	 (this implies signs=0). */
      for (i = 0; i < dim; ++i)
	p[i] = c[i] + r[i];

      /* Loop through the points in Gray-code ordering */
      for (i = 0;; ++i) {
	unsigned mask, d;
	
	for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
	istart+=dim;
	vtmp=o2scl::vector_range(pts,istart,pts.size());
	
	/* which coordinate to flip */
	d = ls0(i); 
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
    static void evalRR0_0fs(vec_t &pts, unsigned dim, double *p,
			    const vec_t &c, const double *r) {
      
      size_t istart=0;
      vec_t vtmp=o2scl::vector_range(pts,istart,pts.size());

      for (unsigned i = 0; i < dim - 1; ++i) {
	p[i] = c[i] - r[i];
	for (unsigned j = i + 1; j < dim; ++j) {
	  p[j] = c[j] - r[j];

	  for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
	  istart+=dim;
	  vtmp=o2scl::vector_range(pts,istart,pts.size());

	  p[i] = c[i] + r[i];

	  for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
	  istart+=dim;
	  vtmp=o2scl::vector_range(pts,istart,pts.size());

	  p[j] = c[j] + r[j];

	  for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
	  istart+=dim;
	  vtmp=o2scl::vector_range(pts,istart,pts.size());

	  p[i] = c[i] - r[i];

	  for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
	  istart+=dim;
	  vtmp=o2scl::vector_range(pts,istart,pts.size());

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
      (vec_t &pts, unsigned dim, double *p, const vec_t &c,
       const double *r1, const double *r2) {

      size_t istart=0;
      vec_t vtmp=o2scl::vector_range(pts,istart,pts.size());
      
      for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
      istart+=dim;
      vtmp=o2scl::vector_range(pts,istart,pts.size());

      for (unsigned i = 0; i < dim; i++) {
	p[i] = c[i] - r1[i];
	for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
	istart+=dim;
	vtmp=o2scl::vector_range(pts,istart,pts.size());
	
	p[i] = c[i] + r1[i];
	for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
	istart+=dim;
	vtmp=o2scl::vector_range(pts,istart,pts.size());
	
	p[i] = c[i] - r2[i];
	for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
	istart+=dim;
	vtmp=o2scl::vector_range(pts,istart,pts.size());
	
	p[i] = c[i] + r2[i];
	for(size_t kk=0;kk<dim;kk++) vtmp[kk]=p[kk];
	istart+=dim;
	vtmp=o2scl::vector_range(pts,istart,pts.size());

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
	cubature_new rule of degree 7 (embedded rule degree 5) 
	from \ref Genz83.
    */
    typedef struct {

      /** \brief Desc */
      rule parent;
      
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
    } rule75genzmalik;
    
    /** \brief Convert integer to double
     */
    static double real(int x) {
      return ((double)(x));
    }
    
    /** \brief Convert double to intger
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
      (rule *r_, unsigned fdim, func_t &f, 
       unsigned nR, region *R) {
    
      // lambda2 = sqrt(9/70)
      // lambda4 = sqrt(9/10)
      // lambda5 = sqrt(9/19)

      const double lambda2 = 0.3585685828003180919906451539079374954541;
      const double lambda4 = 0.9486832980505137995996680633298155601160;
      const double lambda5 = 0.6882472016116852977216287342936235251269;
      const double weight2 = 980.0/6561.0;
      const double weight4 = 200.0/19683.0;
      const double weightE2 = 245.0/486.0;
      const double weightE4 = 25.0/729.0;
      const double ratio = (lambda2 * lambda2) / (lambda4 * lambda4);

      rule75genzmalik *r = (rule75genzmalik *) r_;
      unsigned i, j, iR, dim = r_->dim;
      size_t npts = 0;

      if (alloc_rule_pts(r_, nR)) return o2scl::gsl_failure;
      vec_t &pts=r_->pts;
      vec_t &vals=r_->vals;

      for (iR = 0; iR < nR; ++iR) {
	const vec_t &center = R[iR].h.data;
	const vec_t halfwidth=o2scl::vector_range(R[iR].h.data,dim,
						  R[iR].h.data.size());
	
	for (i = 0; i < dim; ++i) {
	  r->p[i] = center[i];
	}
          
	for (i = 0; i < dim; ++i) {
	  r->widthLambda2[i] = halfwidth[i] * lambda2;
	}
	for (i = 0; i < dim; ++i) {
	  r->widthLambda[i] = halfwidth[i] * lambda4;
	}

	/* Evaluate points in the center, in (lambda2,0,...,0) and
	   (lambda3=lambda4, 0,...,0).  */
	vec_t ptmp=o2scl::vector_range(pts,npts*dim,pts.size());
	evalR0_0fs4d(ptmp,dim,r->p,center, 
		     r->widthLambda2,r->widthLambda);
	npts += num0_0(dim) + 2 * numR0_0fs(dim);

	/* Calculate points for (lambda4, lambda4, 0, ...,0) */
	evalRR0_0fs(ptmp, dim, r->p, center, r->widthLambda);
	npts += numRR0_0fs(dim);

	/* Calculate points for (lambda5, lambda5, ..., lambda5) */
	for (i = 0; i < dim; ++i) {
	  r->widthLambda[i] = halfwidth[i] * lambda5;
	}
	evalR_Rfs(ptmp, dim, r->p, center, r->widthLambda);
	npts += numR_Rfs(dim);
      }

      /* Evaluate the integrand function(s) at all the points */
      if (f(dim, npts, pts, fdim, vals)) {
	return o2scl::gsl_failure;
      }

      /* We are done with the points, and so we can re-use the pts
	 array to store the maximum difference diff[i] in each dimension 
	 for each hypercube. */
      vec_t &diff = pts;
      for (i = 0; i < dim * nR; ++i) diff[i] = 0;
      
      for (j = 0; j < fdim; ++j) {

	unsigned vj=j;
	vec_t v=o2scl::vector_range(vals,vj,vals.size());
    
	for (iR = 0; iR < nR; ++iR) {
	  double result, res5th;
	  double val0, sum2=0, sum3=0, sum4=0, sum5=0;
	  unsigned k, k0 = 0;

	  /* accumulate j-th function values into j-th integrals
	     NOTE: this relies on the ordering of the eval functions
	     above, as well as on the internal structure of
	     the evalR0_0fs4d function */

	  /* central point */
	  val0 = v[fdim*(0)]; 
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
    
	  vj += r_->num_points * fdim;
	  v=o2scl::vector_range(vals,vj,vals.size());
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
      return o2scl::success;
    }

#ifdef O2SCL_NEVER_DEFINED
  }{
#endif
    
    /** \brief Desc
     */
    rule *make_rule75genzmalik(unsigned dim, unsigned fdim) {

      rule75genzmalik *r;
      
      /* this rule does not support 1d integrals */
      if (dim < 2) return 0; 
      
      /* Because of the use of a bit-field in evalR_Rfs, we are
	 limited to be < 32 dimensions (or however many bits are in
	 unsigned). This is not a practical limitation...long before
	 you reach 32 dimensions, the Genz-Malik cubature_new becomes
	 excruciatingly slow and is superseded by other methods (e.g.
	 Monte-Carlo). */
      if (dim >= sizeof(unsigned) * 8) return 0;
      
      r = (rule75genzmalik *) make_rule(sizeof(rule75genzmalik),
					dim, fdim,
					num0_0(dim) + 2 * numR0_0fs(dim)
					+ numRR0_0fs(dim) + numR_Rfs(dim),
					rule75genzmalik_evalError,
					destroy_rule75genzmalik);
      if (!r) return 0;

      r->weight1 = (real(12824 - 9120 * to_int(dim) + 400 * isqr(to_int(dim)))
		    / real(19683));
      r->weight3 = real(1820 - 400 * to_int(dim)) / real(19683);
      r->weight5 = real(6859) / real(19683) / real(1U << dim);
      r->weightE1 = (real(729 - 950 * to_int(dim) + 50 * isqr(to_int(dim)))
		     / real(729));
      r->weightE3 = real(265 - 100 * to_int(dim)) / real(1458);

      r->p = (double *) malloc(sizeof(double) * dim * 3);
      if (!r->p) { destroy_rule((rule *) r); return 0; }
      r->widthLambda = r->p + dim;
      r->widthLambda2 = r->p + 2 * dim;

      return (rule *) r;
    }

    /** \brief 1d 15-point Gaussian quadrature rule
	
	Based on qk15.c and qk.c in GNU GSL (which in turn is based on
	QUADPACK).
    */
    static int rule15gauss_evalError(rule *r,
				     unsigned fdim, func_t &f, 
				     unsigned nR, region *R) {

      static const double DBL_MIN=std::numeric_limits<double>::min();
      static const double DBL_EPSILON=std::numeric_limits<double>::epsilon();

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
      vec_t pts, vals;

      if (alloc_rule_pts(r, nR)) return o2scl::gsl_failure;
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

	unsigned vk0=k;
	vec_t vk=o2scl::vector_range(vals,vk0,vals.size());

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
	  vk0+=15*fdim;
	  vk=o2scl::vector_range(vals,vk0,vals.size());
	}
      }
      return o2scl::success;
    }
     
    /** \brief Desc
     */
    rule *make_rule15gauss(unsigned dim, unsigned fdim) {

      if (dim != 1) return 0; /* this rule is only for 1d integrals */
       
      return make_rule(sizeof(rule),dim,fdim,15,
		       rule15gauss_evalError,0);
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
    typedef struct {
      /** \brief Desc */
      size_t n;
      /** \brief Desc */
      size_t nalloc;
      /** \brief Desc */
      heap_item *items;
      /** \brief Desc */
      unsigned fdim;
      /** array of length fdim of the total integrand & error */
      std::vector<esterr> ee; 
    } heap;

    /** \brief Desc
     */
    void heap_resize(heap *h, size_t nalloc) {

      h->nalloc = nalloc;
      if (nalloc) {
	h->items = (heap_item *) realloc(h->items, sizeof(heap_item)*nalloc);
      } else {
	/* BSD realloc does not free for a zero-sized reallocation */
	free(h->items);
	h->items = 0;
      }
      return;
    }

    /** \brief Desc
     */
    heap heap_alloc(size_t nalloc, unsigned fdim) {

      heap h;
      h.n = 0;
      h.nalloc = 0;
      h.items = 0;
      h.fdim = fdim;
      h.ee.resize(fdim);
      for (unsigned i = 0; i < fdim; ++i) h.ee[i].val = h.ee[i].err = 0;
      heap_resize(&h, nalloc);
      return h;
    }

    /** \brief Note that heap_free does not deallocate anything referenced by
	the items */
    void heap_free(heap *h) {

      h->n = 0;
      heap_resize(h, 0);
      h->fdim = 0;
      h->ee.clear();
      return;
    }

    /** \brief Desc
     */
    int heap_push(heap *h, heap_item hi) {

      int insert;
      unsigned fdim = h->fdim;

      for (unsigned i = 0; i < fdim; ++i) {
	h->ee[i].val += hi.ee[i].val;
	h->ee[i].err += hi.ee[i].err;
      }
      insert = h->n;
      if (++(h->n) > h->nalloc) {
	heap_resize(h, h->n * 2);
	if (!h->items) return o2scl::gsl_failure;
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
      return o2scl::success;
    }

    /** \brief Desc
     */
    int heap_push_many(heap *h, size_t ni, heap_item *hi) {
      for (size_t i = 0; i < ni; ++i) {
	if (heap_push(h, hi[i])) return o2scl::gsl_failure;
      }
      return o2scl::success;
    }

    /** \brief Desc
     */
    heap_item heap_pop(heap *h) {

      heap_item ret;
      int i, n, child;

      if (!(h->n)) {
	O2SCL_ERR("Attempted to pop an empty heap in cubature_new.",
		  o2scl::exc_esanity);
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
	if (largest == i) {
	  break;
	}
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
    //@}

    /** \brief Desc
     */
    int converged(unsigned fdim, const std::vector<esterr> &ee,
		  double reqAbsError, double reqRelError,
		  error_norm norm) {

      unsigned j;

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
	  err = sqrt(dsqr(ee[j].err*serr) + dsqr(ee[j+1].err*serr)) * maxerr;
	  val = sqrt(dsqr(ee[j].val*sval) + dsqr(ee[j+1].val*sval)) * maxval;
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
	    err += dsqr(ee[j].err * serr);
	    val += dsqr(fabs(ee[j].val) * sval);
	  }
	  err = sqrt(err) * maxerr;
	  val = sqrt(val) * maxval;
	  return err <= reqAbsError || err <= val*reqRelError;
	}
	
      }
      
      O2SCL_ERR2("Invalid value of 'norm' in ",
		 "cubature::converged().",o2scl::exc_einval);
      return o2scl::exc_einval; 
    }

    /** \brief Desc
     */
    int rulecubature_new(rule *r, unsigned fdim, func_t &f, 
			 const hypercube *h, size_t maxEval,
			 double reqAbsError, double reqRelError,
			 error_norm norm, vec_t &val, vec_t &err,
			 int parallel) {
      
      size_t numEval = 0;
      heap regions;
      unsigned i, j;
      /* array of regions to evaluate */
      region *R = 0; 
      size_t nR_alloc = 0;
      std::vector<esterr> ee;

      /* norm is irrelevant */
      if (fdim <= 1) norm = ERROR_INDIVIDUAL; 
      /* invalid norm */
      if (norm < 0 || norm > ERROR_LINF) return o2scl::gsl_failure; 

      regions = heap_alloc(1, fdim);
      
      ee.resize(fdim);
     
      nR_alloc = 2;
      R = (region *) malloc(sizeof(region) * nR_alloc);

      R[0] = make_region(*h, fdim);
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
	    if (converged(fdim, ee, reqAbsError, reqRelError, norm)) {
	      /* other regions have small errs */
	      break; 
	    }
	    
	  } while (regions.n > 0 && (numEval < maxEval || !maxEval));

	  if (eval_regions(nR, R, f, r)
	      || heap_push_many(&regions, nR, R)) {
	    goto bad;
	  }

	} else { 

	  /* minimize number of function evaluations */
	  
	  /* get worst region */
	  R[0] = heap_pop(&regions); 
	  if (cut_region(R, R+1)
	      || eval_regions(2, R, f, r)
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
	destroy_region(regions.items[i]);
      }
      
      ee.clear();
      heap_free(&regions);
      free(R);

      return o2scl::success;

    bad:
      
      ee.clear();
      heap_free(&regions);
      free(R);

      return o2scl::gsl_failure;
    }
    
    /** \brief Desc
     */
    int cubature_new(unsigned fdim, func_t &f, 
		     unsigned dim, const vec_t &xmin, const vec_t &xmax, 
		     size_t maxEval, double reqAbsError, double reqRelError, 
		     error_norm norm, vec_t &val, vec_t &err, int parallel) {

      rule *r;
      hypercube h;
      int status;
      unsigned i;
      
      if (fdim == 0) {
	/* nothing to do */
	return o2scl::success;
      }
      if (dim == 0) {
	/* trivial integration */
	if (f(0, 1, xmin, fdim, val)) return o2scl::gsl_failure;
	for (i = 0; i < fdim; ++i) err[i] = 0;
	return o2scl::success;
      }
      r = dim == 1 ? make_rule15gauss(dim, fdim)
	: make_rule75genzmalik(dim, fdim);
      if (!r) { 
	for (i = 0; i < fdim; ++i) {
	  val[i] = 0;
	  err[i] = HUGE_VAL; 
	}
	return o2scl::gsl_failure;
      }
      h = make_hypercube_range(dim, xmin, xmax);
      status=rulecubature_new(r, fdim, f, &h, maxEval, reqAbsError, 
			      reqRelError, norm, val, err, parallel);
			      
      destroy_hypercube(h);
      destroy_rule(r);
      return status;
    }
    
  public:
    
    /** \brief Desc
     */
    int integ(unsigned fdim, func_t &f,
	      unsigned dim, const vec_t &xmin,
	      const vec_t &xmax, size_t maxEval, double reqAbsError,
	      double reqRelError, error_norm norm, vec_t &val,
	      vec_t &err) {
      if (fdim == 0) {
	/* nothing to do */     
	return o2scl::success;
      }
      return cubature_new(fdim,f,dim,xmin,xmax,
			  maxEval,reqAbsError,reqRelError,norm,val,err,0);
    }
    
  };

#ifdef O2SCL_NEVER_DEFINED
}{
#endif
  
  /** \brief Integration by p-adaptive cubature_new from the Cubature
      library

      This class is experimental.

      \hline
      \b Documentation \b adapted \b from \b Cubature_New
      
      This class performs adaptive integration by increasing the
      degree of the cubature_new rule rather than subdividing the domain,
      using products of Clenshaw-Curtis rules. This algorithm may be
      superior to Genz-Malik for smooth integrands lacking
      strongly-localized features, in moderate dimensions.

      \hline      
  */
  template<class func_t, class vec_t> class inte_pcubature_new
    : public inte_cubature_new_base {
      
  protected:
      
    /** \brief Maximum integral dimension
     */
    static const size_t MAXDIM=20;
    
    /** \brief Cache of the values for the m[dim] grid.  

	For adaptive cubature_new, thanks to the nesting of the C-C rules, we
	can re-use the values from coarser grids for finer grids, and the
	coarser grids are also used for error estimation. 
	
	A grid is determined by an m[dim] array, where m[i] denotes
	2^(m[i]+1)+1 points in the i-th dimension.

	If mi < dim, then we only store the values corresponding to
	the difference between the m grid and the grid with m[mi] ->
	m[mi]-1. (m[mi]-1 == -1 corresponds to the trivial grid of one
	point in the center.) 
    */
    class cacheval {
    public:
      /** \brief Desc */
      unsigned m[MAXDIM];
      /** \brief Desc */
      unsigned mi;
      /** \brief Desc */
      vec_t &val;
    };

    /** \brief Desc  array of ncache cachevals c[i] 
     */
    class valcache {
    public:
      /** \brief Desc */
      size_t ncache;
      /** \brief Desc */
      cacheval *c;
    };

    /** \brief Desc
     */
    void free_cachevals(valcache *v) {
      if (!v) return;
      if (v->c) {
	free(v->c);
	v->c = 0;
      }
      v->ncache = 0;
      return;
    }

    /** \brief Desc

	recursive loop over all cubature_new points for the given (m,mi)
	cache entry: add each point to the buffer buf, evaluating all
	at once whenever the buffer is full or when we are done
    */
    int compute_cacheval(const unsigned *m, unsigned mi, 
			 vec_t &val, size_t *vali,
			 unsigned fdim, func_t &f, 
			 unsigned dim, unsigned id, double *p,
			 const vec_t &xmin, const vec_t &xmax,
			 vec_t &buf, size_t nbuf, size_t *ibuf) {

      if (id == dim) {

	/* add point to buffer of points */
	for (size_t kk=0;kk<dim;kk++) {
	  buf[kk+(*ibuf)*dim]=p[kk];
	}
	*ibuf++;

	if (*ibuf == nbuf) {
	  /* flush buffer */
	  vec_t vtmp=o2scl::vector_range(val,*vali,val.size());
	  if (f(dim, nbuf, buf, fdim, vtmp)) {
	    return o2scl::gsl_failure;
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
    size_t num_cacheval(const unsigned *m, unsigned mi, unsigned dim) {

      size_t nval = 1;
      for (unsigned i = 0; i < dim; ++i) {
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
		     unsigned fdim, func_t &f,
		     unsigned dim, const vec_t &xmin,
		     const vec_t &xmax, vec_t &buf, size_t nbuf) {
      
      size_t ic = vc->ncache;
      size_t nval, vali = 0, ibuf = 0;
      double p[MAXDIM];

      vc->c = (cacheval *) realloc(vc->c, sizeof(cacheval) * ++(vc->ncache));
      if (!vc->c) return -1;

      vc->c[ic].mi = mi;
      memcpy(vc->c[ic].m, m, sizeof(unsigned) * dim);
      nval = fdim * num_cacheval(m, mi, dim);
      vc->c[ic].val.resize(nval);

      if (compute_cacheval(m, mi, vc->c[ic].val, &vali, fdim, f,
			   dim, 0, p, xmin, xmax, buf, nbuf, &ibuf)) {
	return o2scl::gsl_failure;
      }

      if (ibuf > 0) {
	/* flush remaining buffer */
	std::vector<double> val2=o2scl::vector_range
	  (vc->c[ic].val,vali,vc->c[ic].val.size());
	return f(dim, ibuf, buf, fdim, val2);
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
    unsigned eval(const unsigned *cm, unsigned cmi, vec_t &cval,
		  const unsigned *m, unsigned md,
		  unsigned fdim, unsigned dim, unsigned id,
		  double weight, vec_t &val) {

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
	/* order of C-C rule */
	unsigned mid = m[id] - (id == md); 
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
	  vec_t vtmp=o2scl::vector_range(cval,voff,cval.size());
	  voff += eval(cm, cmi, vtmp, m, md, fdim, dim, id + 1,
		       weight * w[i], val);
	  voff += eval(cm, cmi, vtmp, m, md, fdim, dim, id + 1,
		       weight * w[i], val);
	}

	voff += (cnx - nx) * fdim * 2
	  * num_cacheval(cm + id+1, cmi - (id+1), dim - (id+1));
      }
      return voff;
    }

    /** \brief Desc

	Loop over all cache entries that contribute to the integral,
	(with m[md] decremented by 1) 
    */
    void evals(valcache vc, const unsigned *m, unsigned md,
	       unsigned fdim, unsigned dim, double V, vec_t &val) {

      for(size_t kk=0;kk<fdim;kk++) val[kk]=0.0;

      for (size_t i = 0; i < vc.ncache; ++i) {
	if (vc.c[i].mi >= dim ||
	    vc.c[i].m[vc.c[i].mi] + (vc.c[i].mi == md) <= m[vc.c[i].mi]) {
	  eval(vc.c[i].m, vc.c[i].mi, vc.c[i].val,
	       m, md, fdim, dim, 0, V, val);
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
    void eval_integral(valcache vc, const unsigned *m, 
		       unsigned fdim, unsigned dim, double V,
		       unsigned *mi, vec_t &val, vec_t &err,
		       vec_t &val1) {

      double maxerr = 0;
      unsigned i, j;
     
      evals(vc, m, dim, fdim, dim, V, val);

      /* error estimates along each dimension by comparing val with
	 lower-order rule in that dimension; overall (conservative)
	 error estimate from maximum error of lower-order rules. */
      for(size_t kk=0;kk<fdim;kk++) err[kk]=0.0;

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
    int converged(unsigned fdim, const vec_t &vals, const vec_t &errs,
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
	  err = sqrt(dsqr(errs[j]*serr) + dsqr(errs[j+1]*serr)) * maxerr;
	  val = sqrt(dsqr(vals[j]*sval) + dsqr(vals[j+1]*sval)) * maxval;
	  if (err > reqAbsError && err > val*reqRelError)
	    return 0;
	}
	/* fdim is odd, do last dimension individually */
	if (j < fdim) {
	  if (errs[j] > reqAbsError && errs[j] > fabs(vals[j])*reqRelError) {
	    return 0;
	  }
	}
	return 1;

      case ERROR_L1:

	{
	  double err = 0, val = 0;
	  for (j = 0; j < fdim; ++j) {
	    err += errs[j];
	    val += fabs(vals[j]);
	  }
	  return err <= reqAbsError || err <= val*reqRelError;
	}
	
      case ERROR_LINF:

	{
	  double err = 0, val = 0;
	  for (j = 0; j < fdim; ++j) {
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
	  for (j = 0; j < fdim; ++j) {
	    double absval = fabs(vals[j]);
	    if (errs[j] > maxerr) maxerr = errs[j];
	    if (absval > maxval) maxval = absval;
	  }
	  serr = maxerr > 0 ? 1/maxerr : 1;
	  sval = maxval > 0 ? 1/maxval : 1;
	  for (j = 0; j < fdim; ++j) {
	    err += dsqr(errs[j] * serr);
	    val += dsqr(fabs(vals[j]) * sval);
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
    int integ_v_buf(unsigned fdim, func_t &f, 
		    unsigned dim, const vec_t &xmin, const vec_t &xmax,
		    size_t maxEval, double reqAbsError, double reqRelError,
		    error_norm norm, unsigned *m,
		    vec_t &buf, size_t *nbuf, size_t max_nbuf,
		    vec_t &val, vec_t &err) {

      int ret = o2scl::gsl_failure;
      double V = 1;
      size_t numEval = 0, new_nbuf;
      unsigned i;
      valcache vc;
      vc.ncache=0;
      vc.c=0;
      vec_t val1;

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
	if (f(0, 1, xmin, fdim, val)) return o2scl::gsl_failure;
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

      new_nbuf = num_cacheval(m, dim, dim);

      if (max_nbuf < 1) max_nbuf = 1;
      if (new_nbuf > max_nbuf) new_nbuf = max_nbuf;
      if (*nbuf < new_nbuf) {
	buf.resize(new_nbuf*dim);
      }

      /* start by evaluating the m=0 cubature_new rule */
      if (add_cacheval(&vc, m, dim, fdim, f, dim, xmin, xmax, 
		       buf, *nbuf) != o2scl::success) {
	goto done;
      }

      val1.resize(fdim);

      while (1) {
	unsigned mi;

	eval_integral(vc, m, fdim, dim, V, &mi, val, err, val1);
	if (converged(fdim, val, err, reqAbsError, reqRelError, norm)
	    || (numEval > maxEval && maxEval)) {
	  ret = o2scl::success;
	  goto done;
	}
	m[mi] += 1;
	if (m[mi] > clencurt_M) {
	  /* FAILURE */
	  goto done; 
	}

	new_nbuf = num_cacheval(m, mi, dim);
	if (new_nbuf > *nbuf && *nbuf < max_nbuf) {
	  *nbuf = new_nbuf;
	  if (*nbuf > max_nbuf) *nbuf = max_nbuf;
	  buf.resize(*nbuf*dim);
	}

	if (add_cacheval(&vc, m, mi, fdim, f, 
			 dim, xmin, xmax, buf, *nbuf) != o2scl::success) {
	  /* FAILURE */
	  goto done; 
	}
	numEval += new_nbuf;
      }

    done:

      free_cachevals(&vc);
      
      return ret;
    }

    /** \brief Desc
     */
    static const size_t DEFAULT_MAX_NBUF=(1U << 20);
    
    /** \brief Desc
     */
    int integ(unsigned fdim, func_t &f,
	      unsigned dim, const vec_t &xmin, const vec_t &xmax,
	      size_t maxEval, double reqAbsError, double reqRelError,
	      error_norm norm, vec_t &val, vec_t &err) {
      
      int ret;
      size_t nbuf = 0;
      unsigned m[MAXDIM];
      vec_t buf;

      memset(m, 0, sizeof(unsigned) * dim);
      /* max_nbuf > 0 to amortize function overhead */
      ret = integ_v_buf(fdim, f, dim, xmin, xmax, 
			maxEval, reqAbsError, reqRelError, norm,
			m, buf, &nbuf, 16, val, err);
      return ret;
    }

  };

}

#endif
