/*
   ───────────────────────────────────────────────────────────────────

   Copyright (C) 2026, Andrew W. Steiner

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
#ifndef O2SCL_CMA_ES_H
#define O2SCL_CMA_ES_H

/** \file cma_es.h
    \brief File defining \ref o2scl::cma_es and \ref o2scl::cma_es_eigen
*/

#include <cmath>
#include <vector>
#include <deque>
#include <algorithm>
#include <limits>
#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <gsl/gsl_eigen.h>

#include <o2scl/rng.h>
#include <o2scl/mmin.h>
#include <o2scl/multi_funct.h>

#ifdef O2SCL_SET_EIGEN
#include <eigen3/Eigen/Dense>
#endif

namespace o2scl {

  /** \brief Multidimensional minimization by covariance matrix
      adaptation evolution strategy (CMA-ES)

      This class minimizes a function using the covariance matrix
      adaptation evolution strategy, a derivative-free stochastic
      optimizer well-suited to non-convex, non-separable, and
      moderately noisy objective functions. At each generation,
      \ref lambda candidate points are sampled from a multivariate
      normal distribution centered at the current mean with
      covariance \f$ \sigma^2 C \f$, the best \ref mu of them (by
      function value) are used to update the mean, and the
      covariance matrix \f$ C \f$ and the global step size \f$
      \sigma \f$ are adapted using the standard cumulative step-size
      adaptation and rank-one/rank-\f$ \mu \f$ update rules.

      \verbatim embed:rst
      This implements the basic :math:`(\mu/\mu_W,\lambda)`-CMA-ES
      of [Hansen01]_ with weighted recombination and, optionally,
      restarts using either the IPOP strategy of [Auger05]_ (which
      restarts with a doubled population size until a function
      evaluation budget is exhausted) or the BIPOP strategy of
      [Hansen09]_ (which alternates restarts between a large,
      IPOP-like population regime and a small population regime with
      a randomly chosen population size and initial step size,
      allocating each restart to whichever regime has consumed the
      smaller share of the evaluation budget so far). The restart
      strategy is selected with \ref restart_mode.
      \endverbatim

      By default, the per-generation eigendecomposition of the
      covariance matrix (used both to sample new points efficiently
      and to compute \f$ C^{-1/2} \f$ for the evolution path update)
      is performed with GSL's <tt>gsl_eigen_symmv()</tt>. If O2scl
      was compiled with Eigen support, \ref cma_es_eigen may be used
      instead, which performs the same eigendecomposition with
      <tt>Eigen::SelfAdjointEigenSolver</tt>.

      A single run (i.e. without restarts) is stopped early using
      the two tolerances inherited from \ref o2scl::mmin_base: \ref
      o2scl::mmin_base::tol_abs is a stopping tolerance on the
      spread of the search distribution (the run stops once \f$
      \sigma \f$ times the largest coordinate of the square root of
      the eigenvalues of \f$ C \f$ falls below this value), and \ref
      o2scl::mmin_base::tol_rel is a stopping tolerance on
      fitness stagnation (the run stops once the range of the best
      function value found in each of the most recent \f$ 10+\lceil
      30 n/\lambda\rceil \f$ generations falls below this value,
      following the "TolFun" criterion of the reference CMA-ES
      implementations). The constructor overrides the generic \ref
      o2scl::mmin_base defaults for both of these (\f$ 10^{-4} \f$)
      with values more appropriate for CMA-ES (\f$ 10^{-11} \f$ for
      \ref o2scl::mmin_base::tol_abs and \f$ 10^{-12} \f$ for \ref
      o2scl::mmin_base::tol_rel).

      \note This class always minimizes; there is no \ref mmin_de()
      override since gradient information is not used.
  */
  template<class func_t=multi_funct,
    class vec_t=boost::numeric::ublas::vector<double> >
    class cma_es : public mmin_base<func_t,func_t,vec_t> {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    /// \name Restart strategy selection
    //@{
    /// No restarts (default)
    static const int cma_es_restart_none=0;
    /// IPOP restarts (population doubling)
    static const int cma_es_restart_ipop=1;
    /// BIPOP restarts (alternating large/small population regimes)
    static const int cma_es_restart_bipop=2;
    //@}

    /** \brief The restart strategy, one of \ref cma_es_restart_none
        (default), \ref cma_es_restart_ipop, or \ref
        cma_es_restart_bipop
    */
    int restart_mode;

    /** \brief The initial global step size (default 0.3)

        This is interpreted relative to the scale of the initial
        guess, i.e. offspring are initially sampled with standard
        deviation \ref sigma0 in each coordinate.
    */
    double sigma0;

    /** \brief The population size (default 0)

        If zero (the default), then \f$ 4+\lfloor 3\ln n\rfloor \f$
        is used, where \f$ n \f$ is the number of parameters, as
        suggested in [Hansen01]_ .
    */
    size_t lambda;

    /** \brief The maximum total number of function evaluations,
        summed over all restarts (default 0)

        If zero (the default), then \f$ 1000 n^2 \f$ is used, where
        \f$ n \f$ is the number of parameters. This limit is only
        relevant when \ref restart_mode is not \ref
        cma_es_restart_none, since it determines when the restart
        loop stops; a single run without restarts is instead limited
        by \ref o2scl::mmin_base::ntrial generations.
    */
    size_t max_evals;

    /** \brief The total number of function evaluations used by the
        most recent call to \ref mmin(), summed over all restarts
        (default 0)

        This is analogous to \ref o2scl::mmin_base::last_ntrial,
        but counts function evaluations rather than generations, and
        follows the convention described in the documentation of
        \ref o2scl::mmin_base.
    */
    size_t last_n_evals;

    /** \brief The population size multiplier used for each IPOP or
        BIPOP large-regime restart (default 2.0)
    */
    double ipop_factor;

    /** \brief The number of restarts performed by the last call to
        \ref mmin() (default 0)
    */
    size_t n_restarts;

    /** \brief Maximum allowed condition number of the covariance
        matrix (default \f$ 10^{14} \f$)

        A single run is stopped early if the ratio of the largest to
        the smallest eigenvalue of \f$ C \f$ exceeds this value.
    */
    double cond_max;

    /// The random number generator
    mutable rng<> rg;

    cma_es() {
      restart_mode=cma_es_restart_none;
      sigma0=0.3;
      lambda=0;
      max_evals=0;
      last_n_evals=0;
      ipop_factor=2.0;
      n_restarts=0;
      cond_max=1.0e14;
      this->ntrial=10000;
      // Override the generic o2scl::mmin_base defaults (both
      // 1.0e-4) with values more appropriate for CMA-ES; see the
      // class documentation for the meaning of these two
      // tolerances in this context.
      this->tol_abs=1.0e-11;
      this->tol_rel=1.0e-12;
    }

    virtual ~cma_es() {
    }

    /// Return string denoting type ("cma_es")
    virtual const char *type() { return "cma_es"; }

    /** \brief Calculate the minimum \c fmin of \c func w.r.t the
        array \c x of size \c nvar

        The initial point \c x is used as the starting mean of the
        search distribution, and the final point found (over all
        restarts, if \ref restart_mode is not \ref
        cma_es_restart_none) is stored back in \c x on exit.
    */
    virtual int mmin(size_t nvar, vec_t &x, double &fmin, func_t &func) {

      if (nvar==0) {
        O2SCL_ERR2("Tried to min. over zero variables ",
                   "in cma_es::mmin().",exc_einval);
      }

      n_restarts=0;

      size_t budget=(max_evals>0) ? max_evals : 1000*nvar*nvar;
      size_t lambda_def=(lambda>0) ? lambda : default_lambda(nvar);

      vec_t xbest(x);
      double fbest;
      size_t evals_used=0;

      // Initial (first, always large-regime) run
      run_single(nvar,x,sigma0,lambda_def,fbest,xbest,func,
                budget,evals_used);

      if (restart_mode==cma_es_restart_ipop) {

        size_t lambda_cur=lambda_def;

        while (evals_used<budget) {
          lambda_cur=(size_t)(((double)lambda_cur)*ipop_factor);
          if (lambda_cur<4) lambda_cur=4;

          vec_t x0(x), xr(nvar);
          double fr;
          size_t used=0;
          run_single(nvar,x0,sigma0,lambda_cur,fr,xr,func,
                    budget-evals_used,used);
          evals_used+=used;
          n_restarts++;

          if (fr<fbest) {
            fbest=fr;
            xbest=xr;
          }
        }

      } else if (restart_mode==cma_es_restart_bipop) {

        size_t lambda_large=lambda_def;
        size_t evals_large=0, evals_small=0;

        while (evals_used<budget) {

          vec_t x0(x), xr(nvar);
          double fr, sigma_r;
          size_t lambda_r;
          size_t used=0;

          if (evals_small==0 || evals_large<=evals_small) {

            // Large-population regime: double lambda as in IPOP
            lambda_large=(size_t)(((double)lambda_large)*ipop_factor);
            if (lambda_large<4) lambda_large=4;
            lambda_r=lambda_large;
            sigma_r=sigma0;

            run_single(nvar,x0,sigma_r,lambda_r,fr,xr,func,
                      budget-evals_used,used);
            evals_large+=used;

          } else {

            // Small-population regime: randomized population size
            // and initial step size, following [Hansen09]_
            double u1=rg.random(), u2=rg.random();
            double lam_small_d=((double)lambda_def)*
              std::pow(0.5*((double)lambda_large)/((double)lambda_def),
                      u1*u1);
            lambda_r=(size_t)lam_small_d;
            if (lambda_r<4) lambda_r=4;
            sigma_r=sigma0*std::pow(10.0,-2.0*u2);

            run_single(nvar,x0,sigma_r,lambda_r,fr,xr,func,
                      budget-evals_used,used);
            evals_small+=used;

          }

          evals_used=evals_large+evals_small;
          n_restarts++;

          if (fr<fbest) {
            fbest=fr;
            xbest=xr;
          }
        }
      }

      x=xbest;
      fmin=fbest;
      last_n_evals=evals_used;

      return 0;
    }

  protected:

    /// \name Strategy state for the run currently in progress
    //@{
    vec_t xmean_, pc_, ps_;
    ubmatrix Cmat_, Bmat_;
    ubvector Dvec_;
    ubvector weights_;
    double sigma_;
    size_t mu_;
    double mueff_, cc_, cs_, c1_, cmu_, damps_, chiN_;
    size_t eigen_interval_, gen_since_eigen_;
    //@}

    /** \brief Return the default population size for a problem of
        dimension \c nvar
    */
    virtual size_t default_lambda(size_t nvar) {
      size_t lam=4+((size_t)(3.0*std::log((double)nvar)));
      if (lam<4) lam=4;
      return lam;
    }

    /** \brief Set up the (fixed, for the duration of a single run)
        strategy parameters for a problem of dimension \c nvar and
        population size \c lam

        This sets \ref mu_, \ref mueff_, \ref weights_, \ref cc_,
        \ref cs_, \ref c1_, \ref cmu_, \ref damps_, \ref chiN_, and
        \ref eigen_interval_ following the standard formulas of
        [Hansen01]_ with only positive recombination weights.
    */
    virtual void setup_strategy_params(size_t nvar, size_t lam) {

      double n=(double)nvar;

      mu_=lam/2;
      if (mu_<1) mu_=1;

      weights_.resize(mu_);
      double wsum=0.0, wsqsum=0.0;
      for(size_t i=0;i<mu_;i++) {
        double w=std::log(((double)mu_)+0.5)-std::log((double)(i+1));
        weights_[i]=w;
        wsum+=w;
      }
      for(size_t i=0;i<mu_;i++) {
        weights_[i]/=wsum;
        wsqsum+=weights_[i]*weights_[i];
      }
      mueff_=1.0/wsqsum;

      cc_=(4.0+mueff_/n)/(n+4.0+2.0*mueff_/n);
      cs_=(mueff_+2.0)/(n+mueff_+5.0);
      c1_=2.0/((n+1.3)*(n+1.3)+mueff_);
      double cmu_alt=2.0*(mueff_-2.0+1.0/mueff_)/
        ((n+2.0)*(n+2.0)+mueff_);
      cmu_=std::min(1.0-c1_,cmu_alt);
      double t=std::sqrt((mueff_-1.0)/(n+1.0))-1.0;
      damps_=1.0+2.0*std::max(0.0,t)+cs_;

      chiN_=std::sqrt(n)*(1.0-1.0/(4.0*n)+1.0/(21.0*n*n));

      eigen_interval_=(size_t)(1.0/((c1_+cmu_)*n*10.0));
      if (eigen_interval_<1) eigen_interval_=1;
    }

    /** \brief Initialize the search distribution state for a new
        run
    */
    virtual void init_state(size_t nvar, vec_t &x0, double sigma_in) {

      xmean_.resize(nvar);
      pc_.resize(nvar);
      ps_.resize(nvar);
      for(size_t i=0;i<nvar;i++) {
        xmean_[i]=x0[i];
        pc_[i]=0.0;
        ps_[i]=0.0;
      }

      Cmat_.resize(nvar,nvar);
      Bmat_.resize(nvar,nvar);
      Dvec_.resize(nvar);
      for(size_t i=0;i<nvar;i++) {
        for(size_t j=0;j<nvar;j++) {
          Cmat_(i,j)=(i==j) ? 1.0 : 0.0;
          Bmat_(i,j)=(i==j) ? 1.0 : 0.0;
        }
        Dvec_[i]=1.0;
      }

      sigma_=sigma_in;
      gen_since_eigen_=0;
    }

    /** \brief Update \ref Bmat_ and \ref Dvec_ from an
        eigendecomposition of \ref Cmat_

        \ref Dvec_ is set to the square roots of the eigenvalues of
        \ref Cmat_ (i.e. the standard deviations along the principal
        axes), and \ref Bmat_ is set to the corresponding matrix of
        orthonormal eigenvectors (as columns), so that \f$ C = B
        D^2 B^T \f$.

        This default implementation uses GSL's
        <tt>gsl_eigen_symmv()</tt>. The Eigen-based alternative is
        \ref cma_es_eigen::eigen_decomp_cov() .
    */
    virtual int eigen_decomp_cov(size_t nvar) {

      gsl_matrix *gC=gsl_matrix_alloc(nvar,nvar);
      for(size_t i=0;i<nvar;i++) {
        for(size_t j=0;j<nvar;j++) {
          gsl_matrix_set(gC,i,j,Cmat_(i,j));
        }
      }

      gsl_vector *geval=gsl_vector_alloc(nvar);
      gsl_matrix *gevec=gsl_matrix_alloc(nvar,nvar);
      gsl_eigen_symmv_workspace *gw=gsl_eigen_symmv_alloc(nvar);

      gsl_eigen_symmv(gC,geval,gevec,gw);

      gsl_eigen_symmv_free(gw);
      gsl_matrix_free(gC);

      for(size_t i=0;i<nvar;i++) {
        double ev=gsl_vector_get(geval,i);
        if (ev<0.0) ev=0.0;
        Dvec_[i]=std::sqrt(ev);
        for(size_t j=0;j<nvar;j++) {
          Bmat_(j,i)=gsl_matrix_get(gevec,j,i);
        }
      }

      gsl_vector_free(geval);
      gsl_matrix_free(gevec);

      return 0;
    }

    /** \brief Perform a single CMA-ES run (no restarts) with the
        specified initial mean, initial step size, and population
        size

        This runs generations until one of several stopping criteria
        is met (the function-evaluation budget \c eval_budget is
        exhausted, \ref o2scl::mmin_base::ntrial generations have
        elapsed, the spread of the distribution falls below \ref
        o2scl::mmin_base::tol_abs, the range of the best function
        value found in each of the most recent \f$ 10+\lceil 30
        n/\lambda\rceil \f$ generations falls below \ref
        o2scl::mmin_base::tol_rel, or the condition number of the
        covariance matrix exceeds \ref cond_max), and stores the
        best point found (not necessarily the final generation's
        mean) in \c xbest with function value \c fmin. The actual
        number of function evaluations used is returned in \c
        evals_used.
    */
    virtual int run_single(size_t nvar, vec_t &x0, double sigma_in,
                           size_t lam, double &fmin, vec_t &xbest,
                           func_t &func, size_t eval_budget,
                           size_t &evals_used) {

      setup_strategy_params(nvar,lam);
      init_state(nvar,x0,sigma_in);

      xbest.resize(nvar);
      for(size_t i=0;i<nvar;i++) xbest[i]=x0[i];
      fmin=func(nvar,xbest);
      evals_used=1;

      std::vector<ubvector> zpts(lam), ypts(lam);
      for(size_t k=0;k<lam;k++) {
        zpts[k].resize(nvar);
        ypts[k].resize(nvar);
      }
      ubvector fvals(lam);
      std::vector<size_t> idx(lam);

      std::normal_distribution<double> ndist(0.0,1.0);

      // History of each generation's best function value, used by
      // the fitness-stagnation stopping criterion below (the
      // "TolFun" criterion of the reference CMA-ES
      // implementations). The window length follows the same
      // formula used there, based on the population size at the
      // start of this run.
      size_t fit_hist_len=10+((size_t)std::ceil
                              (30.0*((double)nvar)/((double)lam)));
      std::deque<double> fit_hist;

      size_t gen=0;

      while (gen<this->ntrial && evals_used<eval_budget) {

        gen++;

        // Sample lambda offspring
        for(size_t k=0;k<lam;k++) {
          for(size_t i=0;i<nvar;i++) {
            zpts[k][i]=ndist(rg.def_engine);
          }
          for(size_t i=0;i<nvar;i++) {
            double sum=0.0;
            for(size_t j=0;j<nvar;j++) {
              sum+=Bmat_(i,j)*Dvec_[j]*zpts[k][j];
            }
            ypts[k][i]=sum;
          }
          vec_t xk(nvar);
          for(size_t i=0;i<nvar;i++) {
            xk[i]=xmean_[i]+sigma_*ypts[k][i];
          }
          fvals[k]=func(nvar,xk);
          evals_used++;
          if (fvals[k]<fmin) {
            fmin=fvals[k];
            for(size_t i=0;i<nvar;i++) xbest[i]=xk[i];
          }
          if (evals_used>=eval_budget) {
            lam=k+1;
            break;
          }
        }

        // Sort by function value (ascending)
        for(size_t k=0;k<lam;k++) idx[k]=k;
        std::sort(idx.begin(),idx.begin()+lam,
                 [&fvals](size_t a, size_t b)
                 { return fvals[a]<fvals[b]; });

        if (lam<mu_) {
          // Ran out of budget mid-generation with too few
          // offspring evaluated to recombine; stop here.
          break;
        }

        // Weighted recombination in y-space (and z-space, for the
        // evolution path updates)
        ubvector ymean(nvar), zmean(nvar);
        for(size_t i=0;i<nvar;i++) {
          ymean[i]=0.0;
          zmean[i]=0.0;
        }
        for(size_t m=0;m<mu_;m++) {
          size_t k=idx[m];
          for(size_t i=0;i<nvar;i++) {
            ymean[i]+=weights_[m]*ypts[k][i];
            zmean[i]+=weights_[m]*zpts[k][i];
          }
        }

        for(size_t i=0;i<nvar;i++) {
          xmean_[i]+=sigma_*ymean[i];
        }

        // Evolution path for step-size control: ps += B*zmean,
        // scaled appropriately
        ubvector Bz(nvar);
        for(size_t i=0;i<nvar;i++) {
          double sum=0.0;
          for(size_t j=0;j<nvar;j++) {
            sum+=Bmat_(i,j)*zmean[j];
          }
          Bz[i]=sum;
        }
        double csfac=std::sqrt(cs_*(2.0-cs_)*mueff_);
        for(size_t i=0;i<nvar;i++) {
          ps_[i]=(1.0-cs_)*ps_[i]+csfac*Bz[i];
        }

        double psnorm=0.0;
        for(size_t i=0;i<nvar;i++) psnorm+=ps_[i]*ps_[i];
        psnorm=std::sqrt(psnorm);

        double hsig_rhs=(1.4+2.0/(nvar+1.0))*chiN_*
          std::sqrt(1.0-std::pow(1.0-cs_,2.0*((double)gen)));
        bool hsig=(psnorm<hsig_rhs);

        double ccfac=std::sqrt(cc_*(2.0-cc_)*mueff_);
        for(size_t i=0;i<nvar;i++) {
          pc_[i]=(1.0-cc_)*pc_[i]+
            (hsig ? ccfac*ymean[i] : 0.0);
        }

        // Rank-one and rank-mu covariance updates
        for(size_t i=0;i<nvar;i++) {
          for(size_t j=0;j<nvar;j++) {
            double cij=(1.0-c1_-cmu_)*Cmat_(i,j)+
              c1_*pc_[i]*pc_[j];
            if (!hsig) {
              cij+=c1_*cc_*(2.0-cc_)*Cmat_(i,j);
            }
            for(size_t m=0;m<mu_;m++) {
              size_t k=idx[m];
              cij+=cmu_*weights_[m]*ypts[k][i]*ypts[k][j];
            }
            Cmat_(i,j)=cij;
          }
        }

        // Step-size update
        sigma_*=std::exp((cs_/damps_)*(psnorm/chiN_-1.0));

        // Periodically refresh the eigendecomposition
        gen_since_eigen_++;
        if (gen_since_eigen_>=eigen_interval_) {
          eigen_decomp_cov(nvar);
          gen_since_eigen_=0;
        }

        // Stopping criteria
        double dmax=0.0, dmin=std::numeric_limits<double>::max();
        for(size_t i=0;i<nvar;i++) {
          if (Dvec_[i]>dmax) dmax=Dvec_[i];
          if (Dvec_[i]<dmin) dmin=Dvec_[i];
        }

        if (this->verbose>0) {
          print_iter(nvar,gen,evals_used,lam,fmin,sigma_,xmean_,
                    dmax,dmin);
        }

        if (sigma_*dmax<this->tol_abs) {
          break;
        }
        if (dmin>0.0 && (dmax*dmax)/(dmin*dmin)>cond_max) {
          break;
        }

        // Fitness-stagnation stopping criterion: track the best
        // function value found among this generation's offspring,
        // and stop once its range over the trailing window of
        // fit_hist_len generations falls below tol_rel.
        fit_hist.push_back(fvals[idx[0]]);
        if (fit_hist.size()>fit_hist_len) {
          fit_hist.pop_front();
        }
        if (this->tol_rel>0.0 && fit_hist.size()>=fit_hist_len) {
          double fit_hist_max=*std::max_element
            (fit_hist.begin(),fit_hist.end());
          double fit_hist_min=*std::min_element
            (fit_hist.begin(),fit_hist.end());
          if (fit_hist_max-fit_hist_min<this->tol_rel) {
            break;
          }
        }
      }

      this->last_ntrial=gen;

      return 0;
    }

    /** \brief Print out iteration information

        Depending on the value of \ref o2scl::mmin_base::verbose,
        this prints out information about the generation just
        completed, following the same convention as \ref
        o2scl::diff_evo::print_iter() and \ref
        o2scl::mmin_simp2::print_iter(): if verbose=0, nothing is
        printed; if verbose>=1, a small, single-generation summary
        (generation number, function evaluations used so far,
        current population size, the best function value found so
        far, and the current global step size \ref sigma_) is
        printed after each generation; if verbose>=2, this is
        followed by the current mean point \c xmean, the standard
        deviations along the principal axes of the search
        distribution (the square roots of the eigenvalues of the
        covariance matrix, i.e. \ref Dvec_), and the resulting
        condition number, prints a prompt, and then waits for a
        keypress (read from \ref o2scl::mmin_base::ins) before
        continuing. All output (including the prompt) goes to \ref
        o2scl::mmin_base::outs, so both streams respect any prior
        call to \ref o2scl::mmin_base::set_verbose_stream().
    */
    virtual void print_iter(size_t nvar, int gen, size_t evals_used,
                            size_t lam, double fmin, double sigma,
                            vec_t &xmean, double dmax, double dmin) {

      (*this->outs) << "cma_es Generation: " << gen << " of "
                    << this->ntrial << ", evals: " << evals_used
                    << ", lambda: " << lam << std::endl;
      (*this->outs) << "  fmin: " << fmin << " sigma: " << sigma
                    << std::endl;

      if (this->verbose>1) {
        (*this->outs) << "  mean: ";
        for(size_t i=0;i<nvar;i++) (*this->outs) << xmean[i] << " ";
        (*this->outs) << std::endl;
        (*this->outs) << "  D (std. dev. along principal axes): ";
        for(size_t i=0;i<nvar;i++) (*this->outs) << Dvec_[i] << " ";
        (*this->outs) << std::endl;
        double cond=(dmin>0.0) ? (dmax*dmax)/(dmin*dmin) :
          std::numeric_limits<double>::infinity();
        (*this->outs) << "  cond(C): " << cond << std::endl;
        (*this->outs) << "Press a key and type enter to continue. ";
        char ch;
        (*this->ins) >> ch;
      }

      return;
    }

  private:

    cma_es<func_t,vec_t>(const cma_es<func_t,vec_t> &);
    cma_es<func_t,vec_t> &operator=(const cma_es<func_t,vec_t> &);

  };

#if defined(O2SCL_SET_EIGEN) || defined(DOXYGEN)

  /** \brief CMA-ES minimization with the covariance matrix
      eigendecomposition performed by Eigen rather than GSL

      This class is identical to \ref cma_es, except that \ref
      eigen_decomp_cov() uses
      <tt>Eigen::SelfAdjointEigenSolver</tt> rather than GSL's
      <tt>gsl_eigen_symmv()</tt>. It is only defined if O2scl was
      compiled with Eigen support (i.e. if \c O2SCL_SET_EIGEN is
      defined).
  */
  template<class func_t=multi_funct,
    class vec_t=boost::numeric::ublas::vector<double> >
    class cma_es_eigen : public cma_es<func_t,vec_t> {

  public:

    virtual ~cma_es_eigen() {
    }

    /// Return string denoting type ("cma_es_eigen")
    virtual const char *type() { return "cma_es_eigen"; }

  protected:

    /** \brief Update \ref cma_es::Bmat_ and \ref cma_es::Dvec_ from
        an eigendecomposition of \ref cma_es::Cmat_, using
        <tt>Eigen::SelfAdjointEigenSolver</tt>

        See \ref cma_es::eigen_decomp_cov() for the meaning of the
        outputs.
    */
    virtual int eigen_decomp_cov(size_t nvar) {

      Eigen::MatrixXd eC(nvar,nvar);
      for(size_t i=0;i<nvar;i++) {
        for(size_t j=0;j<nvar;j++) {
          eC(i,j)=this->Cmat_(i,j);
        }
      }

      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(eC);

      const Eigen::VectorXd &evals=es.eigenvalues();
      const Eigen::MatrixXd &evecs=es.eigenvectors();

      for(size_t i=0;i<nvar;i++) {
        double ev=evals(i);
        if (ev<0.0) ev=0.0;
        this->Dvec_[i]=std::sqrt(ev);
        for(size_t j=0;j<nvar;j++) {
          this->Bmat_(j,i)=evecs(j,i);
        }
      }

      return 0;
    }

  };

#endif

}

#endif
