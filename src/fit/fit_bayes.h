/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2014, Andrew W. Steiner
  
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
/** \file fit_bayes.h
    \brief File defining \ref o2scl::fit_bayes and related classes
*/
#ifndef O2SCL_FIT_BAYES_H
#define O2SCL_FIT_BAYES_H

#include <vector>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/hist.h>
#include <o2scl/rng_gsl.h>
#include <o2scl/interp.h>
#include <o2scl/fit_base.h>
#include <o2scl/expval.h>
#include <o2scl/mcarlo.h>
#include <o2scl/mcarlo_vegas.h>
#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief An unnormalized uniform prior distribution for several 
      variables
  */
  template<class vec_t=boost::numeric::ublas::vector<double> > 
    class uniform_prior : public multi_funct<vec_t> {
    
  public:
  
  /** \brief Return the value of the prior at the specified
      point in parameter space
  */
  virtual double operator()(size_t npar, const vec_t &par) {
    return 1.0;
  }
 
  };

  /** \brief Fit a function to data using Bayesian methods

      This class is experimental.

      This class uses Markov Chain Monte Carlo (MCMC) and
      marginal estimation to give a probability distribution
      for parameters in a fit.

      \future Also make weight_fun() an object of type multi_func_t?
      \future Offer two ways to do the evidence: direct MC or the 
      interpolation method from SLB13
      \future Build upon gen_fit_funct instead of fit_funct?
  */
  template<class fit_func_t=fit_funct<>, class multi_func_t=uniform_prior<>, 
    class vec_t=boost::numeric::ublas::vector<double> > class fit_bayes {

  public:

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef boost::numeric::ublas::vector<int> ubvector_int;

  protected:

  /** \brief The integrand for the evidence
   */
  virtual double integrand(size_t npar, const vec_t &par) {
    return weight_fun(lndat,*lxdat,*lydat,*lyerr,npar,par)*(*pri)(npar,par);
  }

  /// User-specified function
  fit_func_t *ff;

  /// Number of data points
  size_t lndat;

  /// X-values
  vec_t *lxdat;

  /// Y-values
  vec_t *lydat;

  /// Y-errors
  vec_t *lyerr;

  public:
  
  fit_bayes() {
    n_warm_up=100;
    n_iter=10000;
    hsize=20;
    nmeas=20;
  }
 
  /// Number of warmup iterations (default 100)
  size_t n_warm_up;
 
  /// Number of total iterations (default 1000)
  size_t n_iter;
  
  /// Histogram size (default 20)
  size_t hsize;

  /// Number of measurements (default 20)
  size_t nmeas;

  /// Prior distribution
  multi_func_t *pri;

  /// Random number generator
  rng_gsl gr;

  /// Default Monte Carlo integrator
  mcarlo_vegas<> def_inte;

  /** \brief Compute the evidence
   */
  virtual void evidence(size_t ndat, vec_t &xdat, vec_t &ydat,
			vec_t &yerr, size_t npar, vec_t &plo2, 
			vec_t &phi2, multi_func_t &prior_fun,
			double &evi, double &err) {
    lndat=ndat;
    lxdat=&xdat;
    lydat=&ydat;
    lyerr=&yerr;
    pri=&prior_fun;
#ifndef O2SCL_NO_CPP11
    multi_funct11 mfm=
      std::bind(std::mem_fn<double(size_t,const ubvector &)>
		(&fit_bayes::integrand),
		this,std::placeholders::_1,std::placeholders::_2);
#else
    multi_funct_mfptr<fit_bayes,vec_t> mfm(this,&fit_bayes::integrand);
#endif
    def_inte.minteg_err(mfm,npar,plo2,phi2,evi,err);
    return;
  }

  /** \brief The weight function (based on a \f$ \chi^2 \f$
      distribution)
  */
  virtual double weight_fun(size_t ndat, const vec_t &xdat, 
			    const vec_t &ydat, const vec_t &yerr, 
			    size_t npar, const vec_t &par) {
			    
    double weight=1.0;
    for(size_t i=0;i<ndat;i++) {
      weight*=exp(-pow((*ff)(npar,par,xdat[i])-ydat[i],2.0)/
		  2.0/yerr[i]/yerr[i]);
    }
    return weight;
  }

  /** \brief Fit \c ndat data points in \c xdat and \c ydat with 
      errors \c yerr to function \c fitfun with \c npar parameters.

      The initial values of the parameters should be specified in
      \c par. 
  */
  virtual int fit(size_t ndat, vec_t &xdat, vec_t &ydat, vec_t &yerr,
		  size_t npar, vec_t &plo2, vec_t &pmax, vec_t &phi2, 
		  vec_t &plo_err, vec_t &pmax_err, vec_t &phi_err, 
		  fit_func_t &fitfun, multi_func_t &prior_fun) {
    
    ff=&fitfun;
    pri=&prior_fun;

    double weight, weight_old;

    // Create par_old and set up initial weight in weight_old
    vec_t par_old, par;
    par_old.resize(npar);
    par.resize(npar);

    // Choose the first point as the midpoint of the parameter space
    for(size_t i=0;i<npar;i++) {
      par_old[i]=(plo2[i]+phi2[i])/2.0;
    }
    weight_old=weight_fun(ndat,xdat,ydat,yerr,npar,par_old)*
    (*pri)(npar,par_old);
    if (weight_old==0.0) {
      O2SCL_ERR2_RET("Initial weight zero in ",
		     "fit_bayes::fit().",exc_einval);
    }

    // Set up the histograms for marginal estimation
    std::vector<hist> harr(npar);
    for(size_t i=0;i<npar;i++) {
      harr[i].set_bin_edges(uniform_grid_end<>(plo2[i],phi2[i],hsize));
      harr[i].clear_wgts();
    }
    
    // Create ev objects for lower and upper limits of each parameter
    size_t meas_size=((size_t)(((double)n_iter)/((double)nmeas)));
    size_t imeas=meas_size-1;
    expval_vector low(npar,nmeas,1);
    expval_vector high(npar,nmeas,1);
    expval_vector max(npar,nmeas,1);
    ubvector_int nonzero_bins(npar);
    
    // Main iteration
    for(size_t it=0;it<n_iter+n_warm_up;it++) {
      
      // Select new point in parameter space and evaluate it's weight
      for (size_t i=0;i<npar;i++) {
	par[i]=gr.random()*(phi2[i]-plo2[i])+plo2[i];
      }
      weight=weight_fun(ndat,xdat,ydat,yerr,npar,par)*(*pri)(npar,par);

      // Metropolis sampling
      double r=gr.random();
      if (weight>0.0 && (weight>weight_old || r<weight/weight_old)) {
	for (size_t i=0;i<npar;i++) {
	  par_old[i]=par[i];
	}
	weight_old=weight;
      }

      // If we're not in the warm up phase
      if (it>=n_warm_up) {

	for(size_t i=0;i<npar;i++) {
	  // Only update if its in range
	  if (par_old[i]>plo2[i] && par_old[i]<phi2[i]) {
	    harr[i].update(par_old[i]);
	  }
	}
	
	if (it-n_warm_up==imeas) {

	  // From histogram, determine parameter value and error by finding
	  // 68% confidence region
	  std::vector<double> xlo, xhi, xmax;
	  for(size_t i=0;i<npar;i++) {

	    // Correct for non-zero values at the edges
	    std::vector<double> edge_x, edge_y;
	    edge_x.push_back(2.0*harr[i].get_rep_i(0)-
			     harr[i].get_rep_i(1));
	    edge_y.push_back(0.0);
	    for(size_t j=0;j<hsize;j++) {
	      if (harr[i].get_wgt_i(j)>0.0) {
		nonzero_bins[i]++;
	      }
	      edge_x.push_back(harr[i].get_rep_i(j));
	      edge_y.push_back(harr[i].get_wgt_i(j));
	    }
	    edge_x.push_back(2.0*harr[i].get_rep_i(hsize-1)-
			     harr[i].get_rep_i(hsize-2));
	    edge_y.push_back(0.0);

	    // Ensure the integral is unchanged
	    edge_y[1]/=2.0;
	    edge_y[edge_y.size()-2]/=2.0;

	    double total=vector_integ_linear(hsize+2,edge_x,edge_y);
	    
	    double lev;
	    std::vector<double> locs;
	    vector_invert_enclosed_sum(0.68*total,hsize+2,edge_x,edge_y,lev);
	    vector_find_level(lev,hsize+2,edge_x,edge_y,locs);

 	    double lmin=*min_element(locs.begin(),locs.end());
	    double lmax=*max_element(locs.begin(),locs.end());
	    xlo.push_back(lmin);
	    xhi.push_back(lmax);

	    size_t imax=vector_max_index<std::vector<double>,double>
	      (hsize+2,edge_y);
	    xmax.push_back(quadratic_extremum_x
			   (edge_x[imax-1],edge_x[imax],edge_x[imax+1],
			    edge_y[imax-1],edge_y[imax],edge_y[imax+1]));
	  }

	  // Add the limits to the expval_vector object
	  low.add(xlo);
	  high.add(xhi);
	  max.add(xmax);

	  // Clear the histogram data but leave the grid
	  for(size_t i=0;i<npar;i++) {
	    harr[i].clear_wgts();
	  }

	  // Setup for the next measurement
	  imeas+=meas_size;
	}
	
      }
      
    }

    // Now get parameter values and errors from the expval_vector objects
    ubvector sd1(npar), sd2(npar), sd3(npar);
    low.current_avg(plo2,sd1,plo_err);
    high.current_avg(phi2,sd2,phi_err);
    max.current_avg(pmax,sd3,pmax_err);

    for(size_t i=0;i<npar;i++) {
      if (((double)nonzero_bins[i])<((double)nmeas)*2.5) {
	O2SCL_ERR2_RET("Not enough data for accurate parameter ",
		       "estimation in fit_bayes::fit().",exc_efailed);
      }
    }

    return 0;
  }

  /** \brief Desc

      For each measurement block, this function collects the data for
      all the parameters into 1d histogram objects. Then, at the end
      of the block, the histogram information is added to a hist
      object for each parameter.
  */
  virtual int fit_hist(size_t ndat, vec_t &xdat, vec_t &ydat, vec_t &yerr,
		       size_t npar, vec_t &plo2, vec_t &phi2, 
		       std::vector<hist> &par_hist, fit_func_t &fitfun, 
		       multi_func_t &prior_fun) {

    ff=&fitfun;
    pri=&prior_fun;

    double weight, weight_old;

    // Create par_old and set up initial weight in weight_old
    vec_t par_old, par;
    par_old.resize(npar);
    par.resize(npar);

    // Choose the first point as the midpoint of the parameter space
    for(size_t i=0;i<npar;i++) {
      par_old[i]=(plo2[i]+phi2[i])/2.0;
    }
    weight_old=weight_fun(ndat,xdat,ydat,yerr,npar,par_old)*
    (*pri)(npar,par_old);
    if (weight_old==0.0) {
      O2SCL_ERR2_RET("Initial weight zero in ",
		     "fit_bayes::fit().",exc_einval);
    }

    // Set up the histograms for marginal estimation
    par_hist.resize(npar);
    std::vector<hist> harr(npar);
    for(size_t i=0;i<npar;i++) {
      harr[i].set_bin_edges(uniform_grid_end<>(plo2[i],phi2[i],hsize));
      harr[i].clear_wgts();
      par_hist[i].set_bin_edges(uniform_grid_end<>(plo2[i],phi2[i],hsize));
      par_hist[i].clear_wgts();
    }
    
    // Create ev objects for lower and upper limits of each parameter
    size_t meas_size=((size_t)(((double)n_iter)/((double)nmeas)));
    size_t imeas=meas_size-1;
    expval_matrix hist_ev(npar,hsize,nmeas,1);
    ubmatrix mat_tmp(npar,hsize);
    
    // Main iteration
    for(size_t it=0;it<n_iter+n_warm_up;it++) {
      
      // Select new point in parameter space and evaluate it's weight
      for (size_t i=0;i<npar;i++) {
	par[i]=gr.random()*(phi2[i]-plo2[i])+plo2[i];
      }
      weight=weight_fun(ndat,xdat,ydat,yerr,npar,par)*(*pri)(npar,par);

      // Metropolis sampling
      double r=gr.random();
      if (weight>0.0 && (weight>weight_old || r<weight/weight_old)) {
	for (size_t i=0;i<npar;i++) {
	  par_old[i]=par[i];
	}
	weight_old=weight;
      }

      // If we're not in the warm up phase
      if (it>=n_warm_up) {

	for(size_t i=0;i<npar;i++) {
	  // Only update if its in range
	  if (par_old[i]>plo2[i] && par_old[i]<phi2[i]) {
	    harr[i].update(par_old[i]);
	  }
	}
	
	if (it-n_warm_up==imeas) {

	  for(size_t i=0;i<npar;i++) {
	    for(size_t j=0;j<hsize;j++) {
	      mat_tmp(i,j)=harr[i].get_wgt_i(j);
	    }
	  }
	  hist_ev.add(mat_tmp);

	  // Clear the histogram data but leave the grid
	  for(size_t i=0;i<npar;i++) {
	    harr[i].clear_wgts();
	  }

	  // Setup for the next measurement
	  imeas+=meas_size;
	}
	
      }
      
    }

    ubmatrix avg(npar,hsize), std_dev(npar,hsize), avg_err(npar,hsize);
    hist_ev.current_avg(avg,std_dev,avg_err);
    for(size_t i=0;i<npar;i++) {
      for(size_t j=0;j<hsize;j++) {
	par_hist[i].set_wgt_i(j,avg(i,j));
      }
    }

    return 0;
  }
  
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
