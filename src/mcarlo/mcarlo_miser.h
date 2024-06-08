/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2024, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
/* monte/miser.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
#ifndef O2SCL_MCARLO_MISER_H
#define O2SCL_MCARLO_MISER_H

/** \file mcarlo_miser.h
    \brief File defining \ref o2scl::mcarlo_miser
*/
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/misc.h>
#include <o2scl/mcarlo.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_monte_miser.h>

namespace o2scl {
  
  /** \brief Multidimensional integration using the MISER Monte Carlo 
      algorithm (GSL)

      This class uses recursive stratified sampling to estimate the
      value of an integral over a hypercubic region.

      By default the minimum number of calls to estimate the variance
      is 16 times the number of dimensions. This ratio is
      stored in \ref calls_per_dim. By default the minimum number of
      calls per bisection is 32 times \ref calls_per_dim times
      the number of dimensions. This ratio is stored in 
      \ref bisection_ratio. These ratios are employed by
      \ref minteg_err().

      Alternatively, the user can directly set these minimums by \ref
      set_min_calls() and \ref set_min_calls_per_bisection() and use
      \ref miser_minteg_err() which ignores \ref calls_per_dim and
      \ref bisection_ratio. 

      If \ref mcarlo::verbose is greater than zero, then the
      status of the integration is output at every level of bisection
      less than \ref n_levels_out. If it is greater than 1, then the
      boundaries of the current region are also output. Finally, if it
      is greater than 2, a keypress is required after each output.

      \verbatim embed:rst
      Based on [Press90]_.
      \endverbatim
  */
  template<class func_t=multi_funct, 
    class vec_t=boost::numeric::ublas::vector<double>,
           class rng_t=rng<> >
    class mcarlo_miser : public mcarlo<func_t,vec_t,rng_t> {
    
    public:
  
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
  
    /** \brief Number of calls per dimension (default 16)
     */
    size_t calls_per_dim;
  
    /** \brief Factor to set \ref min_calls_per_bisection (default 32)
     */
    size_t bisection_ratio;
  
    /** \brief Introduce random variation into bisection (default 0.0)
      
	From GSL documentation:
	\verbatim
	This parameter introduces a random fractional variation of size
	DITHER into each bisection, which can be used to break the
	symmetry of integrands which are concentrated near the exact
	center of the hypercubic integration region.  The default value of
	dither is zero, so no variation is introduced. If needed, a
	typical value of DITHER is 0.1.
	\endverbatim
    */
    double dither;

    /** \brief Specify fraction of function calls for estimating variance
	(default 0.1)
	
	From GSL documentation:
	\verbatim
	This parameter specifies the fraction of the currently available
	number of function calls which are allocated to estimating the
	variance at each recursive step. The default value is 0.1.
	\endverbatim
    */
    double estimate_frac;

    /** \brief How estimated variances for two sub-regions are combined
	(default 2.0)

	The error handler will be called if this is less than zero.

	From GSL documentation:
	\verbatim
	This parameter controls how the estimated variances for the two
	sub-regions of a bisection are combined when allocating points.
	With recursive sampling the overall variance should scale better
	than 1/N, since the values from the sub-regions will be obtained
	using a procedure which explicitly minimizes their variance.  To
	accommodate this behavior the MISER algorithm allows the total
	variance to depend on a scaling parameter \alpha,
	
	\Var(f) = {\sigma_a \over N_a^\alpha} + {\sigma_b \over N_b^\alpha}.
	
	The authors of the original paper describing MISER recommend the
	value \alpha = 2 as a good choice, obtained from numerical
	experiments, and this is used as the default value in this
	implementation.
	\endverbatim
    */
    double alpha;

    /** \brief Minimum number of calls to estimate the variance

	This is set by minteg() and minteg_err() to be \ref
	calls_per_dim times the number of dimensions in the problem. The
	default value of calls_per_dim is 16 (which is the GSL default).

	From GSL documentation:
	\verbatim
	This parameter specifies the minimum number of function calls
	required for each estimate of the variance. If the number of
	function calls allocated to the estimate using ESTIMATE_FRAC falls
	below MIN_CALLS then MIN_CALLS are used instead.  This ensures
	that each estimate maintains a reasonable level of accuracy. 
	\endverbatim
    */
    int set_min_calls(size_t &mc) {
      min_calls=mc;
      return 0;
    }

    /** \brief Minimum number of calls required to proceed with bisection
	
	This is set by minteg() and minteg_err() to be \ref
	calls_per_dim times \ref bisection_ratio times the number of
	dimensions in the problem. The default values give 512 times the
	number of dimensions (also the GSL default).

	From GSL documentation:
	\verbatim
	This parameter specifies the minimum number of function calls
	required to proceed with a bisection step.  When a recursive step
	has fewer calls available than MIN_CALLS_PER_BISECTION it performs
	a plain Monte Carlo estimate of the current sub-region and
	terminates its branch of the recursion. 
	\endverbatim
    */
    int set_min_calls_per_bisection(size_t &mcb) {
      min_calls_per_bisection=mcb;
      return 0;
    }

    /** \brief The number of recursive levels to output when verbose is 
	greater than zero (default 5)
    */
    size_t n_levels_out;
  
#ifndef DOXYGEN_INTERNAL

    protected:

    /// Minimum number of calls to estimate the variance
    size_t min_calls;

    /// Minimum number of calls required to proceed with bisection
    size_t min_calls_per_bisection;

    /// The number of dimensions of memory allocated
    size_t dim;

    /// \name Arrays which contain a value for each dimension 
    //@{
    /// The current midpoint 
    ubvector xmid;
    /// The left variance
    ubvector sigma_l;
    /// The right variance
    ubvector sigma_r;
    /// The maximum function value in the left half 
    ubvector fmax_l;
    /// The maximum function value in the right half
    ubvector fmax_r;
    /// The minimum function value in the left half
    ubvector fmin_l;
    /// The minimum function value in the right half
    ubvector fmin_r;
    /// The sum in the left half
    ubvector fsum_l;
    /// The sum in the right half 
    ubvector fsum_r;
    /// The sum of the squares in the left half
    ubvector fsum2_l;
    /// The sum of the squares in the right half 
    ubvector fsum2_r;
    /// The number of evaluation points in the left half
    ubvector_size_t hits_l;
    /// The number of evaluation points in the right half
    ubvector_size_t hits_r;
    //@}
    
    /** \brief Estimate the variance

	\future Remove the reference to GSL_POSINF and replace with a
	function parameter.
    */
    virtual int estimate_corrmc(func_t &func, size_t ndim,
				const vec_t &xl, const vec_t &xu,
				size_t calls, double &res,
				double &err, const ubvector &lxmid,
				ubvector &lsigma_l, ubvector &lsigma_r) {
      size_t i, n;
      
      double m=0.0, q=0.0;
      double vol=1.0;

      for (i=0;i<dim;i++) {
	vol*=xu[i]-xl[i];
	hits_l[i]=hits_r[i]=0;
	fsum_l[i]=fsum_r[i]=0.0;
	fsum2_l[i]=fsum2_r[i]=0.0;
	lsigma_l[i]=lsigma_r[i]=-1;
      }

      for (n=0;n<calls;n++) {
	double fval;

	unsigned int j=(n/2) % dim;
	unsigned int side=(n % 2);

	for (i=0;i<dim;i++) {

	  // The equivalent of gsl_rng_uniform_pos()
	  double z;
	  do { 
	    z=this->rng.random();
	  } while (z==0);
	  
	  if (i != j) {
	    x[i]=xl[i]+z*(xu[i]-xl[i]);
	  } else {
	    if (side == 0) {
	      x[i]=lxmid[i]+z*(xu[i]-lxmid[i]);
	    } else {
	      x[i]=xl[i]+z*(lxmid[i]-xl[i]);
	    }
	  }
	}
	
	fval=func(ndim,x);
	
	/* recurrence for mean and variance */
	{
	  double d=fval-m;
	  m+=d/(n+1.0);
	  q+=d*d*(n/(n+1.0));
	}

	/* compute the variances on each side of the bisection */
	for (i=0;i<dim;i++) {
	  if (x[i] <= lxmid[i]) {
	    fsum_l[i]+=fval;
	    fsum2_l[i]+=fval*fval;
	    hits_l[i]++;
	  } else {
	    fsum_r[i]+=fval;
	    fsum2_r[i]+=fval*fval;
	    hits_r[i]++;
	  }
	}
      }

      for (i=0;i<dim;i++) {
	double fraction_l=(lxmid[i]-xl[i])/(xu[i]-xl[i]);

	if (hits_l[i] > 0) {
	  fsum_l[i] /= hits_l[i];
	  lsigma_l[i]=sqrt(fsum2_l[i]-fsum_l[i]*fsum_l[i]/hits_l[i]);
	  lsigma_l[i]*=fraction_l*vol/hits_l[i];
	}
	
	if (hits_r[i] > 0) {
	  fsum_r[i] /= hits_r[i];
	  lsigma_r[i]=sqrt(fsum2_r[i]-fsum_r[i]*fsum_r[i]/hits_r[i]);
	  lsigma_r[i]*=(1-fraction_l)*vol/hits_r[i];
	}
      }
      
      res=vol*m;
      
      if (calls<2) {
	err=GSL_POSINF;
      } else {
	err=vol*sqrt(q/(calls*(calls-1.0)));
      }
      
      return success;
    }

    /// The most recent integration point
    vec_t x;

#endif
    
    public:
    
    mcarlo_miser() {
      estimate_frac=0.1;
      alpha=2.0;
      dither=0.0;
      min_calls=100;
      min_calls_per_bisection=3000;
      dim=0;
      n_levels_out=5;
      calls_per_dim=16;
      bisection_ratio=32;
    }

    virtual ~mcarlo_miser() {
    }
  
    /// Allocate memory
    virtual int allocate(size_t ldim) {

      if (ldim==0) {
	O2SCL_ERR2("Can't allocate zero memory in ",
		       "mcarlo_miser::allocate().",exc_efailed);
      }

      if (ldim!=dim) {
	x.resize(ldim);
	xmid.resize(ldim);
	sigma_l.resize(ldim);
	sigma_r.resize(ldim);
	fmax_l.resize(ldim);
	fmax_r.resize(ldim);
	fmin_l.resize(ldim);
	fmin_r.resize(ldim);
	fsum_l.resize(ldim);
	fsum_r.resize(ldim);
	fsum2_l.resize(ldim);
	fsum2_r.resize(ldim);
	hits_l.resize(ldim);
	hits_r.resize(ldim);
      }
    
      dim=ldim;

      return 0;
    }

    /** \brief Integrate function \c func over the hypercube from
	\f$ x_i=\mathrm{xl}_i \f$ to \f$ x_i=\mathrm{xu}_i \f$ for
	\f$ 0<i< \f$ ndim-1

	\note The values of \ref min_calls and \ref
	min_calls_per_bisection should be set before calling this
	function. The default values if not set are 100 and 3000
	respectively, which correspond to the GSL default setting
	for a 6 dimensional problem. 
    */
    virtual int miser_minteg_err(func_t &func, size_t ndim, const vec_t &xl, 
				 const vec_t &xu, size_t calls, size_t level,
				 double &res, double &err) {

      if (min_calls==0 || min_calls_per_bisection==0) {
	O2SCL_ERR2("Variables min_calls or min_calls_per_bisection ",
		       "are zero in mcarlo_miser::miser_minteg_err().",
		       exc_einval);
      }
      
      size_t n, estimate_calls, calls_l, calls_r;
      size_t i;
      size_t i_bisect;
      int found_best;

      double res_est=0, err_est=0;
      double res_r=0, err_r=0, res_l=0, err_l=0;
      double xbi_l, xbi_m, xbi_r, s;

      double vol;
      double weight_l, weight_r;

      for (i=0;i<dim;i++) {
	if (xu[i] <= xl[i]) {
	  std::string str="Upper limit, "+dtos(xu[i])+", must be greater "+
	    "than lower limit, "+dtos(xl[i])+", in mcarlo_miser::"+
	    "miser_minteg_err().";
	  O2SCL_ERR(str.c_str(),exc_einval);
	}
	if (xu[i]-xl[i]>GSL_DBL_MAX) {
	  O2SCL_ERR2("Range of integration is too large ",
			 "in mcarlo_miser::miser_minteg_err().",exc_einval);
	}
      }

      if (alpha<0) {
	std::string str="Parameter alpha, "+dtos(alpha)+", must be non-"+
	"negative in mcarlo_miser::mister_minteg_err().";
	O2SCL_ERR(str.c_str(),exc_einval);
      }

      /* [GSL] Compute volume */

      vol=1;

      for (i=0;i<dim;i++) {
	vol*=xu[i]-xl[i];
      }

      if (calls<min_calls_per_bisection) {
	double m=0.0, q=0.0;

	if (calls<2) {
	  O2SCL_ERR2("Insufficient calls for subvolume ", 
			 "in mcarlo_miser::miser_minteg_err().",exc_einval);
	}

	for (n=0;n<calls;n++) {
	  /* Choose a random point in the integration region */

	  for (i=0;i<dim;i++) {
	    
	    // The equivalent of gsl_rng_uniform_pos()
	    double rdn;
	    do { 
	      rdn=this->rng.random();
	    } while (rdn==0);
	  
	    x[i]=xl[i]+rdn*(xu[i]-xl[i]);
	  }

	  {
	    double fval;
	    fval=func(ndim,x);
	    
	    /* [GSL] recurrence for mean and variance */

	    double d=fval-m;
	    m+=d/(n+1.0);
	    q+=d*d*(n/(n+1.0));
	  }
	}

	res=vol*m;

	err=vol*sqrt(q/(calls*(calls-1.0)));

	return success;
      }
    
      // The equivalent of what's in GSL but rewritten with explicit
      // typecasts
      size_t prod=(size_t)(((double)calls)*estimate_frac);
      estimate_calls=(min_calls > prod ? min_calls : prod);
      
      if (estimate_calls<4*dim) {
	O2SCL_ERR2("Insufficient calls to sample all halfspaces ", 
		       "in mcarlo_miser::miser_minteg_err().",exc_esanity);
      }

      // [GSL] Flip coins to bisect the integration region with some fuzz 

      for (i=0;i<dim;i++) {
	s=((this->rng.random())-0.5) >= 0.0 ? dither : -dither;
	xmid[i]=(0.5+s)*xl[i]+(0.5-s)*xu[i];
      }

      /* 
	 [GSL] The idea is to chose the direction to bisect based on
	 which will give the smallest total variance. We could (and may
	 do so later) use MC to compute these variances. But the NR guys
	 simply estimate the variances by finding the min and max
	 function values for each half-region for each bisection.
      */
      estimate_corrmc(func,dim,xl,xu,estimate_calls,
		      res_est,err_est,xmid,sigma_l,sigma_r);

      // [GSL] We have now used up some calls for the estimation 

      calls -= estimate_calls;

      // [GSL] Now find direction with the smallest total "variance"

      {
	double best_var=GSL_DBL_MAX;
	double beta=2.0/(1.0+alpha);
	found_best=0;
	i_bisect=0;
	weight_l=weight_r=1.0;

	for (i=0;i<dim;i++) {
	  if (sigma_l[i] >= 0 && sigma_r[i] >= 0) {
	    // [GSL] Estimates are okay 
	    double var=pow (sigma_l[i], beta)+pow (sigma_r[i], beta);
	    
	    if (var <= best_var) {
	      found_best=1;
	      best_var=var;
	      i_bisect=i;
	      weight_l=pow (sigma_l[i], beta);
	      weight_r=pow (sigma_r[i], beta);
	      if (weight_l==0 && weight_r==0) {
		weight_l=1;
		weight_r=1;
	      }
	    }
	  } else {
	    if (sigma_l[i]<0) {
	      O2SCL_ERR2("No points in left-half space ", 
			     "in mcarlo_miser::miser_minteg_err().",
			     exc_esanity);
	    }
	    if (sigma_r[i]<0) {
	      O2SCL_ERR2("No points in right-half space ", 
                         "in mcarlo_miser::miser_minteg_err().",
                         exc_esanity);
	    }
	  }
	}
      }
      
      if (!found_best) {
	// [GSL] All estimates were the same, so chose a direction at
	// random
	i_bisect=((int)((this->rng.random())*(dim-1.0e-10)));
	//std::uniform_int_distribution<int> int_dist(0,dim-1);
	//i_bisect=int_dist(this->rng);
	//gsl_rnga gr;
	//i_bisect=gr.random_int(dim);
      }

      xbi_l=xl[i_bisect];
      xbi_m=xmid[i_bisect];
      xbi_r=xu[i_bisect];

      // [GSL] Get the actual fractional sizes of the two "halves", and
      // distribute the remaining calls among them 
      {
	double fraction_l=fabs ((xbi_m-xbi_l)/(xbi_r-xbi_l));
	double fraction_r=1-fraction_l;

	double a=fraction_l*weight_l;
	double b=fraction_r*weight_r;

	calls_l=(size_t)(min_calls+(calls-2*min_calls)*a/(a+b));
	calls_r=(size_t)(min_calls+(calls-2*min_calls)*b/(a+b));
      }

      if (this->verbose>0 && level<n_levels_out) {
	std::cout << "mcarlo_miser: level,calls_l,calls_r,calls,min_calls_pb: " 
	<< level << " " << calls_l << " " << calls_r << " " << calls << " "
	<< min_calls_per_bisection << std::endl;
	std::cout << "\tres,err: " << res_est << " " << err_est << std::endl;
	if (this->verbose>1) {
	  std::cout << "\ti,left,mid,right: " << i_bisect << " "
		    << xbi_l << " " << xbi_m << " " << xbi_r << std::endl;
	  for(size_t j=0;j<dim;j++) {
	    std::cout << "\t\ti,low,high: " << j << " " << xl[j] << " "
                      << xu[j] << std::endl;
	  }
	}
	if (this->verbose>2) {
	  char ch;
	  std::cin >> ch;
	}
      }

      /* [GSL] Compute the integral for the left hand side of the
	 bisection. Due to the recursive nature of the algorithm we must
	 allocate some new memory for the integration limits for each
	 recursive call
      */
      {
	int status;
	
	vec_t xu_tmp(dim);
	
	for (i=0;i<dim;i++) {
	  xu_tmp[i]=xu[i];
	}
	
	xu_tmp[i_bisect]=xbi_m;
	
	status=miser_minteg_err(func,dim,xl,xu_tmp,calls_l,level+1,
				res_l,err_l);

	if (status != success) {
	  return status;
	}
      }

      // [GSL] Compute the integral for the right hand side of the
      // bisection 
      {
	int status;

	vec_t xl_tmp(dim);

	for (i=0;i<dim;i++) {
	  xl_tmp[i]=xl[i];
	}

	xl_tmp[i_bisect]=xbi_m;

	status=miser_minteg_err(func,dim,xl_tmp,xu,calls_r,level+1,
				res_r,err_r);
	
	if (status != success) {
	  return status;
	}
      }

      res=res_l+res_r;
      err=sqrt(err_l*err_l+err_r*err_r);

      return 0;
    }
    
    /** \brief Integrate function \c func from x=a to x=b.

        The result of the integral is stored in \c res and the
        error estimate in \c err.
        
	This function is just a wrapper to miser_minteg_err() which
	allocates the memory if necessary, sets \c min_calls and \c
	min_calls_per_bisection, calls \ref miser_minteg_err(), and then
	frees the previously allocated memory.
    */
    virtual int minteg_err(func_t &func, size_t ndim, const vec_t &a, 
			   const vec_t &b, double &res, double &err) {
      if (ndim!=dim) allocate(ndim);
      min_calls=calls_per_dim*ndim;
      min_calls_per_bisection=bisection_ratio*min_calls;
      int ret=miser_minteg_err(func,ndim,a,b,this->n_points,0,res,err);
      // Set these to back to zero to ensure that if 
      // miser_minteg_err() is called directly, the user sets 
      // min_calls and min_calls_per_bisection first.
      min_calls=0;
      min_calls_per_bisection=0;
      return ret;
    }
    
    /** \brief Integrate function \c func over the hypercube from
	\f$ x_i=a_i \f$ to \f$ x_i=b_i \f$ for
	\f$ 0<i< \f$ ndim-1
	
	This function is just a wrapper to minteg_err() which allocates
	the memory, sets min_calls and min_calls_per_bisection, calls
	miser_minteg_err(), and then frees the previously allocated
	memory.
    */
    virtual double minteg(func_t &func, size_t ndim, const vec_t &a, 
			  const vec_t &b) {
      double res;
      minteg_err(func,ndim,a,b,res,this->interror);
      return res;
    }

    /// Return string denoting type ("mcarlo_miser")
    virtual const char *type() { return "mcarlo_miser"; }
    
    };
  
}

#endif

