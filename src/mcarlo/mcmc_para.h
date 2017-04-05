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
/** \file mcmc_para.h
    \brief File for definition of \ref o2scl::mcmc_para_base and 
    \ref o2scl::mcmc_para_table
*/
#ifndef O2SCL_MCMC_PARA_H
#define O2SCL_MCMC_PARA_H

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
#include <o2scl/vector.h>
#include <o2scl/multi_funct.h>
#include <o2scl/interpm_idw.h>
#include <o2scl/vec_stats.h>

namespace o2scl {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  
  /** \brief A generic MCMC simulation class

      This class performs a Markov chain Monte Carlo simulation of a
      user-specified function using OpenMP and/or MPI. Either the
      Metropolis-Hastings algorithm with a user-specified proposal
      distribution or the affine-invariant sampling method of Goodman
      and Weare can be used.

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
    class data_t, class vec_t=ubvector> class mcmc_para_base {
    
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
  
  /** \brief Time in seconds (default is 0)
   */
  double max_time;

  /// The screen output file
  std::ofstream scr_out;
  
  /// Random number generators
  std::vector<rng_gsl> rg;
  
  /// Proposal distribution
  o2scl::prob_cond_mdim<vec_t> *prop_dist;
  
  /// If true, then use the user-specified proposal distribution
  bool pd_mode;

  /// If true, we are in the warm up phase
  bool warm_up;

  /** \brief Current points in parameter space for each walker and 
      each OpenMP thread

      This is an array of size \ref n_threads times \ref n_walk initial
      guesses, indexed by <tt>thread_index*n_walk+walker_index</tt> .
  */
  std::vector<vec_t> current;

  /** \brief Data array

      This is an array of size 2 times \ref n_threads times \ref
      n_walk . The two copies of data objects are indexed by
      <tt>i_copy*n_walk*n_threads+thread_index*n_walk+walker_index</tt>
   */
  std::vector<data_t> data_arr;

  /** \brief Data switch array for each walker and each OpenMP thread

      This is an array of size \ref n_threads times \ref n_walk initial
      guesses, indexed by <tt>thread_index*n_walk+walker_index</tt> .
   */
  std::vector<bool> switch_arr;
  
  /** \brief Return value counters, one vector for each OpenMP thread
   */
  std::vector<std::vector<size_t> > ret_value_counts;
  
  /// \name Interface customization
  //@{
  /** \brief Initializations before the MCMC 
   */
  virtual int mcmc_init() {
    if (verbose>0) {
      // Open main output file for this rank
      scr_out.open((prefix+"_"+
		    o2scl::itos(mpi_rank)+"_scr").c_str());
      scr_out.setf(std::ios::scientific);
    }
    return 0;
  }

  /** \brief Cleanup after the MCMC
   */
  virtual void mcmc_cleanup() {
    if (verbose>0) {
      for(size_t it=0;it<n_threads;it++) {
	scr_out << "mcmc (" << it << "): accept=" << n_accept[it]
	<< " reject=" << n_reject[it] << std::endl;
      }
      scr_out.close();
    }
    return;
  }

  /** \brief Function to run for the best point 
   */
  virtual void best_point(vec_t &best, double w_best, data_t &dat) {
    return;
  }
  //@}

  /** \brief Index of the current walker
      
      This quantity has to be a vector because different threads
      may have different values for the current walker during
      the initialization phase for the affine sampling algorithm.
   */
  std::vector<size_t> curr_walker;

  public:

  /// If non-zero, the maximum number of MCMC iterations (default 0)
  size_t max_iters;
  
  /** \brief Prefix for output filenames
   */
  std::string prefix;

  /// Integer to indicate completion
  static const int mcmc_done=-10;

  /// Integer to indicate rejection
  static const int mcmc_skip=-20;

  /// \name Output quantities
  //@{
  /** \brief The number of Metropolis steps which were accepted in 
      each thread (summed over all walkers)
  */
  std::vector<size_t> n_accept;
  
  /** \brief The number of Metropolis steps which were rejected in 
      each thread (summed over all walkers)
  */
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
      (default 0.1)
   */
  double ai_initial_step;
  //@}
  
  mcmc_para_base() {
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
    
    prefix="mcmc";
    max_time=0.0;
    max_iters=0;
  }

  /// Number of OpenMP threads
  size_t n_threads;
  
  /** \brief Initial points in parameter space

      To fully specify all of the initial points, this should be 
      a vector of size \ref n_walk times \ref n_threads .
   */
  std::vector<ubvector> initial_points;
    
  /// \name Basic usage
  //@{
  /** \brief Perform an MCMC simulation

      Perform an MCMC simulation over \c nparams parameters starting
      at initial point \c init, limiting the parameters to be between
      \c low and \c high, using \c func as the objective function and
      calling the measurement function \c meas at each MC point.
  */
  virtual int mcmc(size_t nparams, vec_t &low, vec_t &high,
		   std::vector<func_t> &func, std::vector<measure_t> &meas) {

    if (initial_points.size()==0) {
      // Setup initial guess if not specified
      initial_points.resize(1);
      for(size_t k=0;k<nparams;k++) {
	initial_points[0][k]=(low[k]+high[k])/2.0;
      }
    } else {
      // If initial points are specified, make sure they're within
      // the user-specified limits
      for(size_t iip=0;iip<initial_points.size();iip++) {
	for(size_t ipar=0;ipar<nparams;ipar++) {
	  if (initial_points[iip][ipar]<low[ipar] ||
	      initial_points[iip][ipar]>high[ipar]) {
	    O2SCL_ERR((((std::string)"Parameter ")+o2scl::szttos(ipar)+
		       " of "+o2scl::szttos(nparams)+" out of range (value="+
		       o2scl::dtos(initial_points[iip][ipar])+
		       ") in mcmc_base::mcmc().").c_str(),
		      o2scl::exc_einval);
	  }
	}
      }
    }
    
    // Set number of threads
#ifdef O2SCL_OPENMP
    omp_set_num_threads(n_threads);
    n_threads=omp_get_num_threads();
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

    // Allocate ret_value_counts and curr_walker
    ret_value_counts.resize(n_threads);
    curr_walker.resize(n_threads);

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
	scr_out << "mcmc: Affine-invariant step, n_params="
		<< nparams << ", n_walk=" << n_walk
		<< ", n_threads=" << n_threads << ", n_ranks="
		<< mpi_nprocs << std::endl;
      } else if (pd_mode==true) {
	scr_out << "mcmc: With proposal distribution, n_params="
		<< nparams << ", n_threads=" << n_threads << ", n_ranks="
		<< mpi_nprocs << std::endl;
      } else {
	scr_out << "mcmc: Random-walk w/uniform dist., n_params="
		<< nparams << ", n_threads=" << n_threads << ", n_ranks="
		<< mpi_nprocs << std::endl;
      }
    }
    
    // ---------------------------------------------------
    // Compute initial point and initial weights
    
    // Initial point and weights for stretch-move algorithm
    if (aff_inv) {
      
#ifdef O2SCL_OPENMP
#pragma omp parallel default(shared)
#endif
      {
#ifdef O2SCL_OPENMP
#pragma omp for
#endif
	for(size_t it=0;it<n_threads;it++) {
	  
	  // Initialize each walker in turn
	  for(curr_walker[it]=0;curr_walker[it]<n_walk &&
		mcmc_done_flag[it]==false;curr_walker[it]++) {

	    // Index in storage
	    size_t sindex=n_walk*it+curr_walker[it];
	    
	    size_t init_iters=0;
	    bool done=false;

	    // If we already have a guess, try to use that
	    if (sindex<initial_points.size()) {

	      // Copy from the initial points array
	      for(size_t ipar=0;ipar<nparams;ipar++) {
		current[sindex][ipar]=initial_points[sindex][ipar];
	      }
	      
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
				      curr_walker[it],true,data_arr[sindex]);
		if (meas_ret[it]==mcmc_done) {
		  mcmc_done_flag[it]=true;
		}
		done=true;
	      }
	    }

	    while (!done && !mcmc_done_flag[it]) {
	      
	      // Make a perturbation from the initial point
	      for(size_t ipar=0;ipar<nparams;ipar++) {
		do {
		  current[sindex][ipar]=
		    initial_points[sindex % initial_points.size()][ipar]+
		    (rg[it].random()*2.0-1.0)*
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
					  curr_walker[it],true,
					  data_arr[sindex]);
		  } else {
		    mcmc_done_flag[it]=true;
		  }
		  done=true;
		} else if (init_iters>max_bad_steps) {
		  O2SCL_ERR2("Initial walkers failed in ",
			     "mcmc_para_base::mcmc().",o2scl::exc_einval);
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
	    scr_out << "mcmc (" << it << "): Returned mcmc_done "
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
	for(curr_walker[it]=0;curr_walker[it]<n_walk;curr_walker[it]++) {
	  size_t sindex=n_walk*it+curr_walker[it];
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
	  for(curr_walker[it]=0;curr_walker[it]<n_walk;curr_walker[it]++) {
	    size_t sindex=n_walk*it+curr_walker[it];
	    scr_out.precision(4);
	    scr_out << "mcmc (" << it << "): i_walk: ";
	    scr_out.width((int)(1.0+log10((double)(n_walk-1))));
	    scr_out << curr_walker[it] << " log_wgt: "
		    << w_current[sindex] << " (initial; ai)" << std::endl;
	    scr_out.precision(6);
	  }
	}
      }

    } else {

      // Note that this value is used (e.g. in
      // mcmc_para_table::add_line() ) even if aff_inv is false, so we
      // set it to zero here.
      for(size_t it=0;it<n_threads;it++) {
	curr_walker[it]=0;
      }
      
      // Copy from the initial points array
      for(size_t ipar=0;ipar<nparams;ipar++) {
	current[0][ipar]=initial_points[0][ipar];
      }
      
      // Initial point and weights without stretch-move

      func_ret[0]=func[0](nparams,current[0],w_current[0],data_arr[0]);
      if (func_ret[0]==mcmc_done) {
	if (verbose>=1) {
	  scr_out << "mcmc: Initial point returned mcmc_done."
		  << std::endl;
	}
	mcmc_cleanup();
	return 0;
      }
      if (func_ret[0]!=o2scl::success) {
	if (err_nonconv) {
	  O2SCL_ERR("Initial weight vanished in mcmc_para_base::mcmc().",
		    o2scl::exc_einval);
	}
	return 2;
      }
      
#ifdef O2SCL_OPENMP
#pragma omp parallel default(shared)
#endif
      {
#ifdef O2SCL_OPENMP
#pragma omp for
#endif
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
	    scr_out << "mcmc (" << it << "): Returned mcmc_done "
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
	scr_out.precision(4);
	scr_out << "mcmc (0): "
		<< w_current[0] << " (initial)" << std::endl;
	scr_out.precision(6);
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
      std::vector<double> smove_z(n_threads);
      for(size_t it=0;it<n_threads;it++) {
	curr_walker[it]=0;
	smove_z[it]=0.0;
	// Choose walker to move (same for all threads)
	if (aff_inv) {
	  curr_walker[it]=mcmc_iters % n_walk;
	}
      }
      
#ifdef O2SCL_OPENMP
#pragma omp parallel default(shared)
#endif
      {
#ifdef O2SCL_OPENMP
#pragma omp for
#endif
	for(size_t it=0;it<n_threads;it++) {
	  
	  // ---------------------------------------------------
	  // Select next point
	  
	  if (aff_inv) {
	    // Choose jth walker
	    size_t ij;
	    do {
	      ij=((size_t)(rg[it].random()*((double)n_walk)));
	    } while (ij==curr_walker[it] || ij>=n_walk);
	    
	    // Select z 
	    double p=rg[it].random();
	    double a=step_fac;
	    smove_z[it]=(1.0-2.0*p+2.0*a*p+p*p-2.0*a*p*p+a*a*p*p)/a;
	    
	    // Create new trial point
	    for(size_t i=0;i<nparams;i++) {
	      next[it][i]=current[n_walk*it+ij][i]+
		smove_z[it]*(current[n_walk*it+curr_walker[it]][i]-
			     current[n_walk*it+ij][i]);
	    }
	    
	  } else if (pd_mode) {
	    
	    // Use proposal distribution and compute associated weight
	    (*prop_dist)(current[it],next[it]);
	    q_prop=prop_dist->log_pdf(current[it],next[it])-
	      prop_dist->log_pdf(next[it],current[it]);
	    if (!std::isfinite(q_prop)) {
	      O2SCL_ERR2("Proposal distribution not finite in ",
			 "mcmc_para_base::mcmc().",o2scl::exc_efailed);
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
	    if (switch_arr[n_walk*it+curr_walker[it]]==false) {
	      func_ret[it]=func[it](nparams,next[it],w_next[it],
				    data_arr[it*n_walk+curr_walker[it]+
					     n_walk*n_threads]);
	    } else {
	      func_ret[it]=func[it](nparams,next[it],w_next[it],
				    data_arr[it*n_walk+curr_walker[it]]);
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
	    scr_out << "mcmc (" << it << "): Returned mcmc_done." 
		    << std::endl;
	  } else if (func_ret[it]==mcmc_skip && verbose>=3) {
	    scr_out << "mcmc (" << it
		    << "): Parameter(s) out of range: " << std::endl;
	    scr_out.setf(std::ios::showpos);
	    for(size_t k=0;k<nparams;k++) {
	      scr_out << k << " " << low[k] << " "
		      << next[it][k] << " " << high[k];
	      if (next[it][k]<low[k] || next[it][k]>high[k]) {
		scr_out << " <-";
	      }
	      scr_out << std::endl;
	    }
	    scr_out.unsetf(std::ios::showpos);
	  } else if (func_ret[it]!=o2scl::success &&
		     func_ret[it]!=mcmc_skip) {
	    if (verbose>=2) {
	      scr_out << "mcmc (" << it << "): Function returned failure " 
		      << func_ret[it] << " at point ";
	      for(size_t k=0;k<nparams;k++) {
		scr_out << next[it][k] << " ";
	      }
	      scr_out << std::endl;
	    }
	  }
	}
	if (verbose>=3) {
	  std::cout << "Press a key and type enter to continue. ";
	  char ch;
	  std::cin >> ch;
	}
      }

#ifdef O2SCL_OPENMP
#pragma omp parallel default(shared)
#endif
      {
#ifdef O2SCL_OPENMP
#pragma omp for
#endif
	for(size_t it=0;it<n_threads;it++) {
	  
	  // Index in storage
	  size_t sindex=n_walk*it+curr_walker[it];
	  
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
		meas_ret[it]=meas[it](next[it],w_next[it],curr_walker[it],true,
				      data_arr[sindex+n_threads*n_walk]);
	      } else {
		meas_ret[it]=meas[it](next[it],w_next[it],curr_walker[it],true,
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
				      curr_walker[it],false,data_arr[sindex]);
	      } else {
		meas_ret[it]=meas[it](current[sindex],
				      w_current[sindex],
				      curr_walker[it],false,
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
	  size_t sindex=n_walk*it+curr_walker[it];
	  scr_out.precision(4);
	  scr_out << "mcmc (" << it << "): iter: ";
	  scr_out.width((int)(1.0+log10((double)(nparams-1))));
	  scr_out << mcmc_iters << " i_walk: "
		  << curr_walker[it] << " log_wgt: "
		  << w_current[sindex] << std::endl;
	  scr_out.precision(6);
	}
      }
      
      // Collect best point
      for(size_t it=0;it<n_threads;it++) {
	if (func_ret[it]==o2scl::success && w_best>w_next[it]) {
	  best=next[it];
	  w_best=w_next[it];
	  if (switch_arr[n_walk*it+curr_walker[it]]==false) {
	    best_point(best,w_best,data_arr[curr_walker[it]+n_walk*it+
					    n_threads*n_walk]);
	  } else {
	    best_point(best,w_best,data_arr[curr_walker[it]+n_walk*it]);
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
		       " in mcmc_para_base::mcmc().").c_str(),
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
	    scr_out << "Finished warmup." << std::endl;
	  }
	  
	}
      }

      if (verbose>=3) {
	std::cout << "Press a key and type enter to continue. ";
	char ch;
	std::cin >> ch;
      }

      if (main_done==false && warm_up==false && max_iters>0 &&
	  mcmc_iters==max_iters) {
	scr_out << "mcmc (0): Stopping because number of iterations "
		<< "equal to 'max_iters'." << std::endl;
	main_done=true;
      }
      
      if (main_done==false) {
	// Check to see if we're out of time
#ifdef O2SCL_MPI
	double elapsed=MPI_Wtime()-mpi_start_time;
#else
	double elapsed=time(0)-mpi_start_time;
#endif
	if (max_time>0.0 && elapsed>max_time) {
	  if (verbose>=0) {
	    scr_out << "mcmc (0): Stopping because elapsed > max_time."
		    << std::endl;
	  }
	  main_done=true;
	}
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

  /** \brief A generic MCMC simulation class writing data to a 
      \ref o2scl::table_units object

      This class performs a MCMC simulation and stores the 
      results in a \ref o2scl::table_units object. The
      user must specify the column names and units in 
      \ref set_names_units() before \ref mcmc() is called.

      The function \ref add_line is the measurement function of type
      \c measure_t in the parent. The overloaded function \ref mcmc()
      in this class works a bit differently in that it takes a
      function object (type \c fill_t) of the form
      \code
      int fill_func(const vec_t &pars, double log_weight, 
      std::vector<double> &line, data_t &dat);
      \endcode
      which should store any auxillary values stored in the data
      object to \c line, in order to be added to the table.

      The output table will contain the parameters, the logarithm of
      the function (called "log_wgt") and a multiplying factor called
      "mult". This "fill" function is called only when a step is
      accepted and the multiplier for that row is set to 1. If a
      future step is rejected, then the multiplier is increased by
      one, rather than adding the same row to the table again.

      This class forms the basis of the MCMC used in the Bayesian
      analysis of neutron star mass and radius in
      http://github.com/awsteiner/bamr .

      \note This class is experimental.
  */
  template<class func_t, class fill_t, class data_t, class vec_t=ubvector>
    class mcmc_para_table : public mcmc_para_base<func_t,
    std::function<int(const vec_t &,double,size_t,bool,data_t &)>,
    data_t,vec_t> {

  protected:

  /// Measurement functor type for the parent
  typedef std::function<int(const vec_t &,double,size_t,bool,data_t &)>
  internal_measure_t;
  
  /// Type of parent class
  typedef mcmc_para_base<func_t,internal_measure_t,data_t,vec_t> parent_t;

  /// Column names
  std::vector<std::string> col_names;
    
  /// Column units
  std::vector<std::string> col_units;
    
  /// Main data table for Markov chain
  std::shared_ptr<o2scl::table_units<> > table;

  bool first_write;
  
  /** \brief MCMC initialization function

      This function sets the column names and units.
  */
  virtual int mcmc_init() {

    // -----------------------------------------------------------
    // Init tables
    
    table=std::shared_ptr<o2scl::table_units<> >(new o2scl::table_units<>);
    table->new_column("thread");
    table->new_column("walker");
    table->new_column("mult");
    table->new_column("log_wgt");
    for(size_t i=0;i<col_names.size();i++) {
      table->new_column(col_names[i]);
      if (col_units[i].length()>0) {
	table->set_unit(col_names[i],col_units[i]);
      }
    }
    
    walker_rows.resize(this->n_walk*this->n_threads);
    for(size_t i=0;i<this->n_walk*this->n_threads;i++) {
      walker_rows[i]=-1;
    }

    /*
      if (this->verbose>=2) {
      std::cout << "mcmc: Table column names and units: " << std::endl;
      for(size_t i=0;i<tab->get_ncolumns();i++) {
      std::cout << tab->get_column_name(i) << " "
      << tab->get_unit(tab->get_column_name(i)) << std::endl;
      }
      }
      
      if (this->verbose>=2) {
      std::cout << "End mcmc_para_table::mcmc_init()." << std::endl;
      }
    */
    
    return parent_t::mcmc_init();
  }
  
  /** \brief Fill \c line with data for insertion into the table
   */
  virtual int fill_line(const vec_t &pars, double log_weight, 
			std::vector<double> &line, data_t &dat,
			size_t i_walker, fill_t &fill) {

#ifdef O2SCL_OPENMP
    size_t i_thread=omp_get_thread_num();
#else
    size_t i_thread=0;
#endif

    // Thread
    line.push_back(i_thread);
    // Walker (set later)
    line.push_back(i_walker);
    // Initial multiplier
    line.push_back(1.0);
    line.push_back(log_weight);
    for(size_t i=0;i<pars.size();i++) {
      line.push_back(pars[i]);
    }
    int tempi=fill(pars,log_weight,line,dat);
    return tempi;
  }
  
  /** \brief Record the last row in the table which corresponds
      to each walker
  */
  std::vector<int> walker_rows;

  /// Likelihood estimator
  interpm_idw<double *> esti;
  
  public:

  /** \brief Write MCMC tables to files
   */
  virtual void write_files() {

#ifdef O2SCL_MPI
    // Ensure that multiple threads aren't writing to the
    // filesystem at the same time
    int tag=0, buffer=0;
    if (this->mpi_nprocs>1 && this->mpi_rank>0) {
      MPI_Recv(&buffer,1,MPI_INT,this->mpi_rank-1,tag,MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
    }
#endif
    
    o2scl_hdf::hdf_file hf;
    std::string fname=this->prefix+"_"+o2scl::itos(this->mpi_rank)+"_out";
    hf.open_or_create(fname);

    if (first_write==false) {
      hf.set_szt("max_iters",this->max_iters);
      hf.sets("prefix",this->prefix);
      hf.seti("aff_inv",this->aff_inv);
      hf.setd("step_fac",this->step_fac);
      hf.setd("ai_initial_step",this->ai_initial_step);
      hf.set_szt("n_warm_up",this->n_warm_up);
      hf.seti("user_seed",this->user_seed);
      hf.seti("mpi_rank",this->mpi_rank);
      hf.seti("mpi_nprocs",this->mpi_nprocs);
      hf.seti("verbose",this->verbose);
      hf.set_szt("max_bad_steps",this->max_bad_steps);
      hf.set_szt("n_walk",this->n_walk);
      hf.set_szt("n_threads",this->n_threads);
      hf.seti("always_accept",this->always_accept);
      hf.seti("allow_estimates",allow_estimates);
      first_write=true;
    }
    
    hf.set_szt_vec("n_accept",this->n_accept);
    hf.set_szt_vec("n_reject",this->n_reject);
    hf.set_szt_arr2d_copy("ret_value_counts",this->ret_value_counts.size(),
			  this->ret_value_counts[0].size(),
			  this->ret_value_counts);
    hf.setd_arr2d_copy("initial_points",this->initial_points.size(),
		       this->initial_points[0].size(),
		       this->initial_points);
    
    hdf_output(hf,*table,"markov_chain0");
    
    hf.close();
    
#ifdef O2SCL_MPI
    if (this->mpi_nprocs>1 && this->mpi_rank>0) {
      MPI_Send(&buffer,1,MPI_INT,this->mpi_rank+1,tag,MPI_COMM_WORLD);
    }
#endif
    
    return;
  }
  
  /// If true, allow estimates of the weight (default false)
  bool allow_estimates;
  
  mcmc_para_table() {
    allow_estimates=false;
    first_write;
  }
  
  /// \name Basic usage
  //@{
  /** \brief Set the table names and units
   */
  virtual void set_names_units(std::vector<std::string> names,
			       std::vector<std::string> units) {
    col_names=names;
    col_units=units;
    return;
  }
  
  /** \brief Perform an MCMC simulation
      
      Perform an MCMC simulation over \c nparams parameters starting
      at initial point \c init, limiting the parameters to be between
      \c low and \c high, using \c func as the objective function and
      calling the measurement function \c meas at each MC point.
  */
  virtual int mcmc(size_t nparams, 
		   vec_t &low, vec_t &high, std::vector<func_t> &func,
		   std::vector<fill_t> &fill) {
    
    // Set number of threads (this is done in the child as well, but
    // we need this number to set up the vector of measure functions
    // below).
#ifdef O2SCL_OPENMP
    omp_set_num_threads(this->n_threads);
    this->n_threads=omp_get_num_threads();
#else
    this->n_threads=1;
#endif

    // Setup the vector of measure functions
    std::vector<internal_measure_t> meas(this->n_threads);
    for(size_t it=0;it<this->n_threads;it++) {
      meas[it]=std::bind
	(std::mem_fn<int(const vec_t &,double,size_t,bool,
			 data_t &, size_t, fill_t &)>
	 (&mcmc_para_table::add_line),this,std::placeholders::_1,
	 std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,
	 std::placeholders::_5,it,std::ref(fill[it]));
    }
    
    return parent_t::mcmc(nparams,low,high,func,meas);
  }
  
  /** \brief Get the output table
   */
  std::shared_ptr<o2scl::table_units<> > get_table() {
    return table;
  }
  
  /** \brief Set the output table
   */
  void set_table(std::shared_ptr<o2scl::table_units<> > &t) {
    table=t;
    return;
  }
  
  /** \brief Determine the chain sizes
   */
  void get_chain_sizes(std::vector<size_t> &csizes) {

    size_t ntot=this->n_threads*this->n_walk;
    csizes.resize(ntot);
    
    for(size_t it=0;it<this->n_threads;it++) {
      for(size_t iw=0;iw<this->n_walk;iw++) {
	size_t ix=it*this->n_walk+iw;
	size_t istart=ix;
	csizes[ix]=0;
	for(size_t j=istart;j<table->get_nlines();j+=ntot) {
	  if (table->get("mult",j)>0.5) csizes[ix]++;
	}
      }
    }
    
    return;
  }
  
  /** \brief A measurement function which adds the point to the
      table
  */
  virtual int add_line(const vec_t &pars, double log_weight,
		       size_t walker_ix, bool new_meas, data_t &dat,
		       size_t i_thread, fill_t &fill) {

    // The combined walker/thread index 
    size_t windex=i_thread*this->n_walk+walker_ix;

    // The total number of walkers * threads
    size_t ntot=this->n_threads*this->n_walk;

    int ret_value=o2scl::success;
    
    // Make sure only one thread writes to the table at a time by
    // making this next region 'critical'. This is required because
    // writes to the table may require resizing the table structure.
#ifdef O2SCL_OPENMP
#pragma omp critical (o2scl_mcmc_para_table_add_line)
#endif
    {

      // If there's not enough space in the table for this iteration,
      // create it
      if ((walker_rows[windex]<0 && table->get_nlines()<ntot) ||
	  table->get_nlines()<=walker_rows[windex]+ntot) {
	size_t istart=table->get_nlines();
	// Create enough space
	table->set_nlines(table->get_nlines()+ntot);
	// Now additionally initialize the first four colums
	for(size_t j=0;j<this->n_threads;j++) {
	  for(size_t i=0;i<this->n_walk;i++) {
	    table->set("thread",istart+j*this->n_walk+i,j);
	    table->set("walker",istart+j*this->n_walk+i,i);
	    table->set("mult",istart+j*this->n_walk+i,0.0);
	    table->set("log_wgt",istart+j*this->n_walk+i,-1.0);
	  }
	}
      }

      // Test to see if we need to add a new line of data or increment
      // the weight on the previous line. If the fill function has reset
      // the table data, then the walker_rows will refer to a row which
      // doesn't currently exist, so we have to add a new line of data.

      if (new_meas==true || walker_rows[windex]<0) {

	// We need to create a new measurement, so update the
	// row corresponding to this index
	if (walker_rows[windex]>=0) {
	  walker_rows[windex]+=ntot;
	} else {
	  walker_rows[windex]=windex;
	}
      
	if (walker_rows[windex]>=((int)(table->get_nlines()))) {
	  O2SCL_ERR("Not enough space in table.",o2scl::exc_esanity);
	}

	std::vector<double> line;
	int fret=fill_line(pars,log_weight,line,dat,walker_ix,fill);

	// The fill_line() function doesn't set the walker index,
	// so we do this here
	//line[1]=walker_ix;
      
	if (fret!=o2scl::success) {
	  // If we're done, we stop before adding the last point to the
	  // table. This is important because otherwise the last line in
	  // the table will always only have unit multiplicity, which
	  // may or may not be correct.
	  ret_value=this->mcmc_done;
	} else {
      
	  if (line.size()!=table->get_ncolumns()) {
	    std::cout << "line: " << line.size() << " columns: "
		      << table->get_ncolumns() << std::endl;
	    O2SCL_ERR("Table misalignment in mcmc_para_table::add_line().",
		      exc_einval);
	  }
	  
	  table->set_row(((size_t)walker_rows[windex]),line);
	  
	}
      
      } else {
	
	// Otherwise, just increment the multiplier on the previous line
      
	double mult_old=table->get("mult",walker_rows[windex]);
	table->set("mult",walker_rows[windex],mult_old+1.0);
      
      }

    }
    // End of parallel region
    
    return ret_value;
  }
  //@}
  
  /** \brief Desc
   */
  virtual void mcmc_cleanup() {

    // This section removes empty rows at the end of the
    // table that were allocated but not used
    int i;
    bool done=false;
    for(i=table->get_nlines()-1;i>=0 && done==false;i--) {
      done=true;
      if (table->get("mult",i)<1.0e-10) {
	done=false;
      } 
    }
    if (i+2<((int)table->get_nlines())) {
      table->set_nlines(i+2);
    }
    return parent_t::mcmc_cleanup();
  }

  /** \brief Desc
   */
  virtual void ac_coeffs(size_t ncols, ubmatrix &ac_coeffs) {
    std::vector<size_t> csizes;
    get_chain_sizes(csizes);
    size_t min_size=csizes[0];
    for(size_t i=1;i<csizes.size();i++) {
      if (csizes[i]<min_size) min_size=csizes[i];
    }
    size_t N_max=min_size/2;
    ac_coeffs.resize(ncols,N_max-1);
    for(size_t i=0;i<ncols;i++) {
      for(size_t ell=1;ell<N_max;ell++) {
	ac_coeffs(i,ell-1)=0.0;
      }
    }
    size_t n_tot=n_threads*n_walk;
    size_t table_row=0;
    size_t cstart=table->lookup_column("log_wgt")+1;
    for(size_t i=0;i<ncols;i++) {
      for(size_t j=0;j<n_threads;j++) {
	for(size_t k=0;k<n_walk;k++) {
	  size_t tindex=j*n_walk+k;
	  for(size_t ell=1;ell<N_max;ell++) {
	    double mean=o2scl::vector_mean_double
	      (csizes[tindex]+1,&((*table)[cstart+i][table_row]));
	    ac_coeffs(i,ell-1)+=o2scl::vector_lagk_autocorr
	      (csizes[tindex]+1,&((*table)[cstart+i][table_row]),
	       ell,mean);
	  }
	  table_row+=csizes[tindex]+1;
	}
      }
      for(size_t ell=1;ell<N_max;ell++) {
	ac_coeffs(i,ell-1)/=((double)n_tot);
      }
    }
    return;
  }

  /** \brief Desc
   */
  virtual void ac_lengths(size_t ncols, ubmatrix &ac_coeffs_cols,
			  ubvector &ac_lengths) {
    size_t N_max=ac_coeffs_cols.size2();
    ac_lengths.resize(ncols);
    for(size_t icol=0;icol<ncols;icol++) {
      std::vector<double> tau(N_max);
      for(size_t i=5;i<N_max;i++) {
	double sum=0.0;
	for(size_t j=0;j<i;j++) {
	  sum+=ac_coeffs_cols(icol,j);
	}
	tau[i]=1.0+2.0*sum;
	std::cout << tau[i] << " " << ((double)i)/5.0 << std::endl;
      }
      std::cout << std::endl;
    }
    return;
  }
  
  /** \brief Reorder the table by thread and walker index
   */
  virtual void reorder_table() {

    // Create a new table
    std::shared_ptr<o2scl::table_units<> > table2=
    std::shared_ptr<o2scl::table_units<> >(new o2scl::table_units<>);

    for(size_t i=0;i<n_threads;i++) {
      for(size_t j=0;j<n_walk;j++) {
	std::string func=std::string("abs(walker-")+o2scl::szttos(j)+
	  ")<0.1 && abs(thread-"+o2scl::szttos(i)+")<0.1";
	table->copy_rows(func,*table2);
      }
    }
    
    return;
  }
  
  /** \brief Reaverage the data into blocks of a fixed
      size in order to avoid autocorrelations
      
      \note The number of blocks \c n_blocks must be larger than the
      current table size. This function expects to find a column named
      "mult" which contains the multiplicity of each column, as is the
      case after a call to \ref mcmc_para_base::mcmc().
      
      This function is useful to remove autocorrelations to the table
      so long as the autocorrelation length is shorter than the block
      size. This function does not compute the autocorrelation length
      to check that this is the case.
  */
  void reblock(size_t n_blocks) {
    
    for(size_t it=0;it<this->n_threads;it++) {
      
      size_t n=table->get_nlines();
      if (n_blocks>n) {
	O2SCL_ERR2("Cannot reblock. Not enough data in ",
		   "mcmc_para_table::reblock().",o2scl::exc_einval);
      }
      size_t n_block=n/n_blocks;
      size_t m=table->get_ncolumns();
      for(size_t j=0;j<n_blocks;j++) {
	double mult=0.0;
	ubvector dat(m);
	for(size_t i=0;i<m;i++) {
	  dat[i]=0.0;
	}
	for(size_t k=j*n_block;k<(j+1)*n_block;k++) {
	  mult+=(*table)["mult"][k];
	  for(size_t i=1;i<m;i++) {
	    dat[i]+=(*table)[i][k]*(*table)["mult"][k];
	  }
	}
	table->set("mult",j,mult);
	for(size_t i=1;i<m;i++) {
	  dat[i]/=mult;
	  table->set(i,j,dat[i]);
	}
      }
      table->set_nlines(n_blocks);

    }
    
    return;
  }
  
  };
  
  // End of namespace
}

#endif
