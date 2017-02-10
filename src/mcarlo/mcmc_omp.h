/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2017, Andrew W. Steiner
  
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
/** \file mcmc_omp.h
    \brief File for definition of \ref o2scl::mcmc_omp_base and 
    \ref o2scl::mcmc_omp_table
*/
#ifndef O2SCL_MCMC_OMP_H
#define O2SCL_MCMC_OMP_H

#include <iostream>
#include <random>

#ifdef O2SCL_OPENMP
#include <omp.h>
#endif
#ifdef O2SCL_MPI
#include <mpi.h>
#endif

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/uniform_grid.h>
#include <o2scl/table3d.h>
#include <o2scl/hdf_file.h>
#include <o2scl/exception.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/cholesky.h>
#include <o2scl/vector.h>
#include <o2scl/multi_funct.h>

/** \brief Main namespace
    
    This file is documented in mcmc_omp.h .
*/
namespace o2scl {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  
  /** \brief A generic MCMC simulation class

      This class generates a Markov chain of a user-specified
      function. The chain can be generated using the
      Metropolis-Hastings algorithm with a user-specified proposal
      distribution or using the affine-invariant sampling method of
      Goodman and Weare.

      By default, the Metropolis-Hastings algorithm is executed with a
      simple walk, with steps in each dimension of size \f$
      (\mathrm{high} - \mathrm{low})/\mathrm{step\_fac} \f$ with the
      denominator specified in \ref step_fac.

      The function type is a template type, \c func_t, which should
      be of the form 
      \code
      int f(size_t num_of_parameters, const vec_t &parameters,
      double &log_pdf, data_t &dat)
      \endcode
      which computes \c log_pdf, the natural logarithm of the function
      value, for any point in parameter space (any point between \c
      low and \c high ).

      If the function being simulated returns \ref mcmc_skip then the
      point is automatically rejected. After each acceptance or
      rejection, a user-specified "measurement" function (of type \c
      measure_t ) is called, which can be used to store the results.
      In order to stop the simulation, either this function or the
      probability distribution being simulated should return the value
      \ref mcmc_done . 
      
      A generic proposal distribution can be specified in \ref
      set_proposal(). To go back to the default random walk method,
      one can call the function \ref unset_proposal().

      If \ref aff_inv is set to true and the number of walkers, \ref
      n_walk is set to a number larger than 1, then affine-invariant
      sampling is used. For affine-invariant sampling, the variable 
      \ref step_fac represents the value of \f$ a \f$, the 
      limits of the distribution for \f$ z \f$.

      In order to store data at each point, the user can store this
      data in any object of type \c data_t . If affine-invariant
      sampling is used, then each chain has it's own data object. The
      class keeps twice as many copies of these data object as would
      otherwise be required, in order to avoid copying of data objects
      in the case that the steps are accepted or rejected.

      \note This class is experimental.
  */
  template<class func_t, class measure_t,
    class data_t, class vec_t=ubvector> class mcmc_omp_base {
    
  protected:
  
  /// \name MPI properties
  //@{
  /// The MPI processor rank
  int mpi_rank;

  /// The MPI number of processors
  int mpi_nprocs;

  /// The MPI starting time
  double mpi_start_time;
  //@}
  
  /// Random number generators
  std::vector<rng_gsl> rg;
  
  /// Proposal distribution
  o2scl::prob_cond_mdim<vec_t> *prop_dist;
  
  /// If true, then use the user-specified proposal distribution
  bool pd_mode;

  /// If true, we are in the warm up phase
  bool warm_up;

  /// Current points in parameter space
  std::vector<vec_t> current;

  /// Data array
  std::vector<data_t> data_arr;

  /// Data switch array
  std::vector<bool> switch_arr;
  
  /// Return value counters
  std::vector<std::vector<size_t> > ret_value_counts;
  
  /// \name Interface customization
  //@{
  /** \brief Initializations before the MCMC 
   */
  virtual int mcmc_init() {
    return 0;
  }

  /** \brief Cleanup after the MCMC
   */
  virtual void mcmc_cleanup() {
    return;
  }

  /** \brief Function to run for the best point 
   */
  virtual void best_point(vec_t &best, double w_best, data_t &dat) {
    return;
  }
  //@}

  /// Index of the current walker
  size_t curr_walker;

  /// Number of initial points specified by the user;
  size_t n_init_points;
  
  public:

  /// Integer to indicate completion
  static const int mcmc_done=-10;

  /// Integer to indicate rejection
  static const int mcmc_skip=-20;

  /// \name Output quantities
  //@{
  /// The number of Metropolis steps which were accepted in each thread
  std::vector<size_t> n_accept;
  
  /// The number of Metropolis steps which were rejected in each thread
  std::vector<size_t> n_reject;
  //@}

  /// \name Settings
  //@{
  /// If true, use affine-invariant Monte Carlo
  bool aff_inv;
  
  /// Stepsize factor (default 10.0)
  double step_fac;

  /** \brief Number of warm up steps (successful steps not
      iterations) (default 0)
	
      \note Not to be confused with <tt>warm_up</tt>, which is 
      a boolean local variable in some functions not an int.
  */
  size_t n_warm_up;

  /** \brief If non-zero, use as the seed for the random number 
      generator (default 0)
  */
  int user_seed;

  /// Output control (default 0)
  int verbose;

  /** \brief Maximum number of failed steps when generating initial points
      with affine-invariant sampling (default 1000)
  */
  size_t max_bad_steps;

  /** \brief Number of walkers for affine-invariant MC or 1 
      otherwise (default 1)
  */
  size_t n_walk;

  /** \brief If true, call the error handler if msolve() or
      msolve_de() does not converge (default true)
  */
  bool err_nonconv;

  /** \brief If true, accept all steps
   */
  bool always_accept;

  /** \brief Initial step fraction for affine-invariance sampling walkers
   */
  double ai_initial_step;
  //@}
  
  mcmc_omp_base() {
    user_seed=0;
    n_warm_up=0;

    // MC step parameters
    aff_inv=false;
    pd_mode=false;
    step_fac=10.0;
    n_walk=1;
    err_nonconv=true;
    verbose=0;
    warm_up=false;
    max_bad_steps=1000;

    prop_dist=0;
    always_accept=false;
    ai_initial_step=0.1;

    n_init_points=0;
    n_threads=1;

    // Initial values for MPI paramers
    mpi_nprocs=1;
    mpi_rank=0;
    mpi_start_time=0.0;

#ifdef O2SCL_MPI
    // Get MPI rank, etc.
    MPI_Comm_rank(MPI_COMM_WORLD,&this->mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&this->mpi_nprocs);
#endif
    
  }

  /// Requested number of threads
  size_t n_threads;
  
  /// \name Basic usage
  //@{
  /** \brief Perform an MCMC simulation

      Perform an MCMC simulation over \c nparams parameters starting
      at initial point \c init, limiting the parameters to be between
      \c low and \c high, using \c func as the objective function and
      calling the measurement function \c meas at each MC point.
  */
  virtual int mcmc(size_t nparams, vec_t &init,
		   vec_t &low, vec_t &high, std::vector<func_t> &func,
		   std::vector<measure_t> &meas) {

    bool valid_guess=true;
    for(size_t k=0;k<nparams;k++) {
      if (init[k]<low[k] || init[k]>high[k]) valid_guess=false;
    }
    if (!valid_guess) {
      O2SCL_ERR2("Initial guess outside of user-specified ",
		 "lower and upper bounds in mcmc_base::mcmc_omp().",
		 o2scl::exc_einval);
    }
    
    // Set number of threads
#ifdef O2SCL_OPENMP
    omp_set_num_threads(n_threads);
#else
    n_threads=1;
#endif

    // Set starting time
#ifdef O2SCL_MPI
    mpi_start_time=MPI_Wtime();
#else
    mpi_start_time=time(0);
#endif

    // Storage for return values from each thread
    std::vector<int> func_ret(n_threads), meas_ret(n_threads);
      
    // Fix the number of walkers if it is too small
    if (aff_inv) {
      if (n_walk<=1) n_walk=2;
      if (n_walk<n_threads) n_walk=n_threads;
    }
    
    // Fix 'step_fac' if it's less than or equal to zero
    if (step_fac<=0.0) {
      if (aff_inv) step_fac=2.0;
      else step_fac=10.0;
    }
    
    // Set RNGs with a different seed for each thread and rank
    rg.resize(n_threads);
    unsigned long int seed=time(0);
    if (this->user_seed!=0) {
      seed=this->user_seed;
    }
    for(size_t it=0;it<n_threads;it++) {
      seed*=(mpi_rank*n_threads+it+1);
      rg[it].set_seed(seed);
    }
    
    // Keep track of successful and failed MH moves
    n_accept.resize(n_threads);
    n_reject.resize(n_threads);
    for(size_t it=0;it<n_threads;it++) {
      n_accept[it]=0;
      n_reject[it]=0;
    }

    // Warm-up flag, not to be confused with 'n_warm_up', i.e. the
    // number of warm_up iterations.
    warm_up=true;
    if (n_warm_up==0) warm_up=false;

    // Storage size required
    size_t ssize=n_walk*n_threads;

    // Allocate current point and current weight
    current.resize(ssize);
    std::vector<double> w_current(ssize);
    for(size_t i=0;i<ssize;i++) {
      current[i].resize(nparams);
      w_current[i]=0.0;
    }

    // Allocate ret_value_counts
    ret_value_counts.resize(n_threads);

    // Initialize data and switch arrays
    data_arr.resize(2*ssize);
    switch_arr.resize(ssize);
    for(size_t i=0;i<switch_arr.size();i++) switch_arr[i]=false;
    
    // Next point and next weight for each thread
    std::vector<vec_t> next(n_threads);
    for(size_t it=0;it<n_threads;it++) {
      next[it].resize(nparams);
    }
    std::vector<double> w_next(n_threads);

    // Best point over all threads
    vec_t best(nparams);
    double w_best;

    // Generally, these flags are are true for any thread if func_ret
    // or meas_ret is equal to mcmc_done.
    std::vector<bool> mcmc_done_flag(n_threads);
    for(size_t it=0;it<n_threads;it++) {
      mcmc_done_flag[it]=false;
    }
	  
    // Proposal weight
    double q_prop;
    
    // Set current[0] to user-specified initial point
    for(size_t i=0;i<nparams;i++) {
      current[0][i]=init[i];
    }

    // Run mcmc_init() function. The initial point, stored in
    // current[0] can be modified by this function and the local
    // variable 'init' is not accessible to the mcmc_init() function.
    int init_ret=mcmc_init();
    if (init_ret!=0) {
      O2SCL_ERR("Function mcmc_init() failed in mcmc_base::mcmc().",
		o2scl::exc_einval);
      return init_ret;
    }

    // Initial verbose output
    if (verbose>=1) {
      if (aff_inv) {
	std::cout << "mcmc: Affine-invariant step, n_params="
		  << nparams << " n_walk=" << n_walk
		  << ", n_threads=" << n_threads << ", n_ranks="
		  << mpi_nprocs << std::endl;
      } else if (pd_mode==true) {
	std::cout << "mcmc: With proposal distribution, n_params="
		  << nparams << ", n_threads=" << n_threads << ", n_ranks="
		  << mpi_nprocs << std::endl;
      } else {
	std::cout << "mcmc: Random-walk w/uniform dist., n_params="
		  << nparams << ", n_threads=" << n_threads << ", n_ranks="
		  << mpi_nprocs << std::endl;
      }
    }
    
    // ---------------------------------------------------
    // Compute initial point and initial weights
    
    // Initial point and weights for stretch-move algorithm
    if (aff_inv) {
      
      // The mcmc_init() function may have changed the
      // initial point, so we copy back to init here
      for(size_t ipar=0;ipar<nparams;ipar++) {
	init[ipar]=current[0][ipar];
      }

#pragma omp parallel default(shared)
      {
#pragma omp for
	for(size_t it=0;it<n_threads;it++) {
	  
	  // Initialize each walker in turn
	  for(curr_walker=0;curr_walker<n_walk &&
		mcmc_done_flag[it]==false;curr_walker++) {

	    // Index in storage
	    size_t sindex=n_walk*it+curr_walker;
	    
	    size_t init_iters=0;
	    bool done=false;
	    
	    // If we already have a guess, try to use that
	    if (sindex<n_init_points) {
	      
	      // Compute the weight
	      func_ret[it]=func[it](nparams,current[sindex],
				    w_current[sindex],data_arr[sindex]);
	      
	      if (func_ret[it]==mcmc_done) {
		mcmc_done_flag[it]=true;
	      } else if (func_ret[it]==o2scl::success) {

		// If we have a good point, update ret_value_counts
		// and call the measurement function 
		if (func_ret[it]>=0 &&
		    func_ret[it]<((int)ret_value_counts[it].size())) {
		  ret_value_counts[it][func_ret[it]]++;
		}
		meas_ret[it]=meas[it](current[sindex],w_current[sindex],
				      curr_walker,true,data_arr[sindex]);
		if (meas_ret[it]==mcmc_done) {
		  mcmc_done_flag[it]=true;
		}
		done=true;
	      }
	    }
	    
	    while (!done && !mcmc_done_flag[it]) {
	      
	      // Make a perturbation from the initial point
	      for(size_t ipar=0;ipar<nparams;ipar++) {
		if (init[ipar]<low[ipar] || init[ipar]>high[ipar]) {
		  O2SCL_ERR((((std::string)"Parameter ")+o2scl::szttos(ipar)+
			     " of "+o2scl::szttos(nparams)+
			     " out of range (value="+o2scl::dtos(init[ipar])+
			     ") in mcmc_base::mcmc().").c_str(),
			    o2scl::exc_einval);
		}
		do {
		  current[sindex][ipar]=init[ipar]+(rg[it].random()*2.0-1.0)*
		    (high[ipar]-low[ipar])*ai_initial_step;
		} while (current[sindex][ipar]>high[ipar] ||
			 current[sindex][ipar]<low[ipar]);
	      }
	      
	      // Compute the weight
	      func_ret[it]=func[it](nparams,current[sindex],
				    w_current[sindex],data_arr[sindex]);
		
	      // ------------------------------------------------
	      
	      // Increment iteration count
	      init_iters++;
	      
	      if (func_ret[it]==mcmc_done) {
		mcmc_done_flag[it]=true;
	      } else {
		// If we have a good point, update ret_value_counts,
		// call the measurement function and stop the loop
		if (func_ret[it]==o2scl::success) {
		  if (func_ret[it]>=0 &&
		      func_ret[it]<((int)ret_value_counts[it].size())) {
		    ret_value_counts[it][func_ret[it]]++;
		  }
		  if (meas_ret[it]!=mcmc_done) {
		    meas_ret[it]=meas[it](current[sindex],w_current[sindex],
					  curr_walker,true,data_arr[sindex]);
		  } else {
		    mcmc_done_flag[it]=true;
		  }
		  done=true;
		} else if (init_iters>max_bad_steps) {
		  O2SCL_ERR2("Initial walkers failed in ",
			     "mcmc_omp_base::mcmc().",o2scl::exc_einval);
		}
	      }
	    }
	  }
	}
      }
      // End of parallel region

      // Stop early if mcmc_done was returned
      bool stop_early=false;
      for(size_t it=0;it<n_threads;it++) {
	if (mcmc_done_flag[it]==true) {
	  if (verbose>=1) {
	    std::cout << "mcmc (" << it << "): Returned mcmc_done "
		      << "(initial; ai)." << std::endl;
	  }
	  stop_early=true;
	}
      }
      if (stop_early) {
	mcmc_cleanup();
	return 0;
      }

      // Set initial values for best point
      w_best=w_current[0];
      size_t best_index=0;
      for(size_t it=0;it<n_threads;it++) {
	for(curr_walker=0;curr_walker<n_walk;curr_walker++) {
	  size_t sindex=n_walk*it+curr_walker;
	  if (w_current[sindex]>w_current[0]) {
	    best_index=sindex;
	    w_best=w_current[sindex];
	  }
	}
      }
      best=current[best_index];
      best_point(best,w_best,data_arr[best_index]);

      // Verbose output
      if (verbose>=2) {
	for(size_t it=0;it<n_threads;it++) {
	  for(curr_walker=0;curr_walker<n_walk;curr_walker++) {
	    size_t sindex=n_walk*it+curr_walker;
	    std::cout.precision(4);
	    std::cout << "mcmc: ";
	    std::cout.width((int)(1.0+log10((double)(n_walk-1))));
	    std::cout << it << " " << curr_walker << " "
		      << w_current[sindex] << " (initial; ai)" << std::endl;
	    std::cout.precision(6);
	  }
	}
      }

    } else {

      // Initial point and weights without stretch-move

      func_ret[0]=func[0](nparams,current[0],w_current[0],data_arr[0]);
      if (func_ret[0]==mcmc_done) {
	std::cout << "mcmc: Initial point returned mcmc_done."
		  << std::endl;
	mcmc_cleanup();
	return 0;
      }
      if (func_ret[0]!=o2scl::success) {
	if (err_nonconv) {
	  O2SCL_ERR("Initial weight vanished in mcmc_omp_base::mcmc().",
		    o2scl::exc_einval);
	}
	return 2;
      }
      
#pragma omp parallel default(shared)
      {
#pragma omp for
	for(size_t it=0;it<n_threads;it++) {
	  // Copy the results over from the initial point
	  if (it!=0) {
	    func_ret[it]=func_ret[0];
	    current[it]=current[0];
	    w_current[it]=w_current[0];
	    data_arr[it]=data_arr[0];
	  }
	  // Update the return value count
	  if (func_ret[it]>=0 &&
	      func_ret[it]<((int)ret_value_counts[it].size())) {
	    ret_value_counts[it][func_ret[it]]++;
	  }
	  // Call the measurement function	  
	  meas_ret[it]=meas[it](current[it],w_current[it],0,
				true,data_arr[it]);
	  if (meas_ret[it]==mcmc_done) {
	    mcmc_done_flag[it]=true;
	  }
	}
      }
      // End of parallel region

      // Stop early if mcmc_done was returned
      bool stop_early=false;
      for(size_t it=0;it<n_threads;it++) {
	if (mcmc_done_flag[it]==true) {
	  if (verbose>=1) {
	    std::cout << "mcmc (" << it << "): Returned mcmc_done "
		      << "(initial)." << std::endl;
	  }
	  stop_early=true;
	}
      }
      if (stop_early) {
	mcmc_cleanup();
	return 0;
      }

      // Set initial values for best point
      best=current[0];
      w_best=w_current[0];
      best_point(best,w_best,data_arr[0]);

      if (verbose>=2) {
	std::cout.precision(4);
	std::cout << "mcmc (0): "
		  << w_current[0] << " (initial)" << std::endl;
	std::cout.precision(6);
      }
      
    }

    if (verbose>=3) {
      std::cout << "Press a key and type enter to continue. ";
      char ch;
      std::cin >> ch;
    }

    // End of initial point and weight section
    // ---------------------------------------------------

    bool main_done=false;
    size_t mcmc_iters=0;
    
    while (!main_done) {

      // Walker to move (or zero when aff_inv is false)
      curr_walker=0;
      std::vector<double> smove_z(n_threads);
      for(size_t it=0;it<n_threads;it++) {
	smove_z[it]=0.0;
      }
      
      // Choose walker to move (same for all threads)
      if (aff_inv) {
	curr_walker=mcmc_iters % n_walk;
      }
      
#pragma omp parallel default(shared)
      {
#pragma omp for
	for(size_t it=0;it<n_threads;it++) {
	  
	  // ---------------------------------------------------
	  // Select next point
	  
	  if (aff_inv) {
	    // Choose jth walker
	    size_t ij;
	    do {
	      ij=((size_t)(rg[it].random()*((double)n_walk)));
	    } while (ij==curr_walker || ij>=n_walk);
	    
	    // Select z 
	    double p=rg[it].random();
	    double a=step_fac;
	    smove_z[it]=(1.0-2.0*p+2.0*a*p+p*p-2.0*a*p*p+a*a*p*p)/a;
	    
	    // Create new trial point
	    for(size_t i=0;i<nparams;i++) {
	      next[it][i]=current[n_walk*it+ij][i]+
		smove_z[it]*(current[n_walk*it+curr_walker][i]-
			     current[n_walk*it+ij][i]);
	    }
	    
	  } else if (pd_mode) {
	    
	    // Use proposal distribution and compute associated weight
	    (*prop_dist)(current[it],next[it]);
	    q_prop=prop_dist->log_pdf(current[it],next[it])-
	      prop_dist->log_pdf(next[it],current[it]);
	    if (!std::isfinite(q_prop)) {
	      O2SCL_ERR2("Proposal distribution not finite in ",
			 "mcmc_omp_base::mcmc().",o2scl::exc_efailed);
	    }
	    
	  } else {
	    
	    // Uniform random-walk step
	    for(size_t k=0;k<nparams;k++) {
	      next[it][k]=current[it][k]+(rg[it].random()*2.0-1.0)*
		(high[k]-low[k])/step_fac;
	    }
	    
	  }
	  
	  // ---------------------------------------------------
	  // Compute next weight
      
	  func_ret[it]=o2scl::success;
	  // If the next point out of bounds, ensure that the point is rejected
	  for(size_t k=0;k<nparams;k++) {
	    if (next[it][k]<low[k] || next[it][k]>high[k]) {
	      func_ret[it]=mcmc_skip;
	    }
	  }
	  if (func_ret[it]!=mcmc_skip) {
	    if (switch_arr[it]==false) {
	      func_ret[it]=func[it](nparams,next[it],w_next[it],
				    data_arr[it+n_threads*n_walk]);
	    } else {
	      func_ret[it]=func[it](nparams,next[it],w_next[it],data_arr[it]);
	    }
	    if (func_ret[it]==mcmc_done) {
	      mcmc_done_flag[it]=true;
	    } else {
	      if (func_ret[it]>=0 &&
		  func_ret[it]<((int)ret_value_counts[it].size())) {
		ret_value_counts[it][func_ret[it]]++;
	      }
	    }

	  }
	}
      }
      // End of parallel region

      // Handle verbose output serially
      if (verbose>=1) {
	for(size_t it=0;it<n_threads;it++) {
	  if (func_ret[it]==mcmc_done) {
	    std::cout << "mcmc (" << it << "): Returned mcmc_done." 
		      << std::endl;
	  } else if (func_ret[it]==mcmc_skip && verbose>=3) {
	    std::cout << "mcmc (" << it
		      << "): Parameter(s) out of range: " << std::endl;
	    std::cout.setf(std::ios::showpos);
	    for(size_t k=0;k<nparams;k++) {
	      std::cout << k << " " << low[k] << " "
			<< next[it][k] << " " << high[k];
	      if (next[it][k]<low[k] || next[it][k]>high[k]) {
		std::cout << " <-";
	      }
	      std::cout << std::endl;
	    }
	    std::cout.unsetf(std::ios::showpos);
	  } else if (func_ret[it]!=o2scl::success &&
		     func_ret[it]!=mcmc_skip) {
	    if (verbose>=2) {
	      std::cout << "mcmc (" << it << "): Function returned failure " 
			<< func_ret[it] << " at point ";
	      for(size_t k=0;k<nparams;k++) {
		std::cout << next[it][k] << " ";
	      }
	      std::cout << std::endl;
	    }
	  }
	}
	if (verbose>=3) {
	  std::cout << "Press a key and type enter to continue. ";
	  char ch;
	  std::cin >> ch;
	}
      }

#pragma omp parallel default(shared)
      {
#pragma omp for
	for(size_t it=0;it<n_threads;it++) {
	  
	  // Index in storage
	  size_t sindex=n_walk*it+curr_walker;
	  
	  // ---------------------------------------------------
	  // Accept or reject
    
	  bool accept=false;
	  if (always_accept && func_ret[it]==success) accept=true;

	  if (func_ret[it]==o2scl::success) {
	    double r=rg[it].random();
	    
	    if (aff_inv) {
	      double ai_ratio=pow(smove_z[it],((double)nparams)-1.0)*
		exp(w_next[it]-w_current[sindex]);
	      if (r<ai_ratio) {
		accept=true;
	      }
	    } else if (pd_mode) {
	      if (r<exp(w_next[it]-w_current[sindex]+q_prop)) {
		accept=true;
	      }
	    } else {
	      // Metropolis algorithm
	      if (r<exp(w_next[it]-w_current[sindex])) {
		accept=true;
	      }
	    }

	    // End of 'if (func_ret[it]==o2scl::success)'
	  }

	  if (accept) {
	  
	    n_accept[it]++;
	  
	    // Store results from new point
	    if (!warm_up) {
	      if (switch_arr[sindex]==false) {
		meas_ret[it]=meas[it](next[it],w_next[it],curr_walker,true,
				      data_arr[sindex+n_threads*n_walk]);
	      } else {
		meas_ret[it]=meas[it](next[it],w_next[it],curr_walker,true,
				      data_arr[sindex]);
	      }
	    }

	    // Prepare for next point
	    current[sindex]=next[it];
	    w_current[sindex]=w_next[it];
	    switch_arr[sindex]=!(switch_arr[sindex]);
	  
	  } else {
	    
	    // Point was rejected
	    n_reject[it]++;

	    // Repeat measurement of old point
	    if (!warm_up) {
	      if (switch_arr[sindex]==false) {
		meas_ret[it]=meas[it](current[sindex],
				      w_current[sindex],
				      curr_walker,false,data_arr[sindex]);
	      } else {
		meas_ret[it]=meas[it](current[sindex],
				      w_current[sindex],
				      curr_walker,false,
				      data_arr[sindex+n_walk*n_threads]);
	      }
	    }

	  }

	}
      }
      // End of parallel region

      // Verbose output
      if (verbose>=2) {
	for(size_t it=0;it<n_threads;it++) {
	  size_t sindex=n_walk*it+curr_walker;
	  std::cout.precision(4);
	  std::cout << "mcmc (" << it << "): ";
	  std::cout.width((int)(1.0+log10((double)(nparams-1))));
	  std::cout << mcmc_iters << " "
		    << curr_walker << " " << w_current[sindex]
		    << std::endl;
	  std::cout.precision(6);
	}
      }
      
      // Collect best point
      for(size_t it=0;it<n_threads;it++) {
	if (func_ret[it]==o2scl::success && w_best>w_next[it]) {
	  best=next[it];
	  w_best=w_next[it];
	  if (switch_arr[curr_walker]==false) {
	    best_point(best,w_best,data_arr[curr_walker+n_walk]);
	  } else {
	    best_point(best,w_best,data_arr[curr_walker]);
	  }
	}
      }

      // Check to see if mcmc_done was returned or if meas_ret
      // returned an error
      for(size_t it=0;it<n_threads;it++) {
	if (meas_ret[it]==mcmc_done || func_ret[it]==mcmc_done) {
	  main_done=true;
	}
	if (meas_ret[it]!=mcmc_done && meas_ret[it]!=o2scl::success) {
	  if (err_nonconv) {
	    O2SCL_ERR((((std::string)"Measurement function returned ")+
		       o2scl::dtos(meas_ret[it])+
		       " in mcmc_omp_base::mcmc().").c_str(),
		      o2scl::exc_efailed);
	  }
	  main_done=true;
	}
      }

      if (main_done==false) {
	
	mcmc_iters++;
	
	if (warm_up && mcmc_iters==n_warm_up) {
	  warm_up=false;
	  mcmc_iters=0;
	  for(size_t it=0;it<n_threads;it++) {
	    n_accept[it]=0;
	    n_reject[it]=0;
	  }
	  if (verbose>=1) {
	    std::cout << "Finished warmup." << std::endl;
	  }
	  
	}
      }
      
      if (verbose>=3) {
	std::cout << "Press a key and type enter to continue. ";
	char ch;
	std::cin >> ch;
      }

      // --------------------------------------------------------------
      // End of main loop
    }
    
    // --------------------------------------------------------------
    
    mcmc_cleanup();
    
    return 0;
  }
  //@}

  /// \name Proposal distribution
  //@{
  /** \brief Set the proposal distribution
   */
  virtual void set_proposal(o2scl::prob_cond_mdim<vec_t> &p) {
    prop_dist=&p;
    pd_mode=true;
    aff_inv=false;
    n_walk=1;
    return;
  }

  /** \brief Go back to random-walk Metropolis with a uniform distribution
   */
  virtual void unset_proposal() {
    if (pd_mode) {
      prop_dist=0;
      pd_mode=false;
    }
    aff_inv=false;
    n_walk=1;
    return;
  }
  //@}

  };

  // End of namespace
}

#endif
