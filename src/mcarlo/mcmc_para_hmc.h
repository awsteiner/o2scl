/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2012-2023, Andrew W. Steiner
  
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
/** \file mcmc_para.h
    \brief File for definition of \ref o2scl::mcmc_para_base,
    \ref o2scl::mcmc_para_table and \ref o2scl::mcmc_para_cli
*/
#ifndef O2SCL_MCMC_PARA_HMC_H
#define O2SCL_MCMC_PARA_HMC_H

#include <iostream>
#include <random>

#include <o2scl/set_openmp.h>

#ifdef O2SCL_SET_OPENMP
#include <omp.h>
#endif
#ifdef O2SCL_MPI
#include <mpi.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/hdf_file.h>
#include <o2scl/exception.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/vector.h>
#include <o2scl/multi_funct.h>
#include <o2scl/vec_stats.h>
#include <o2scl/cli.h>
#include <o2scl/mmin.h>

namespace o2scl {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

  /** \brief A generic MCMC simulation class

      Significant changes:
      
      * There is no data vector in the class, it's moved to a
      parameter of the mcmc() function.

      * New outside_parallel() function: the idea is that 
      an emulator can be retrained in this function.

      * New steps_in_parallel variable

      * The best point mechanism has been reworked, there is now a
      best point over all threads and a best point array which stores
      the best point from each thread

      * The initial point evaluation is essentially the same, 
      but the main loop is reorganized.

      * When aff_inv is false, there are three loops, a main loop,
      a parallel loop over threads, and then an inner loop
      of size steps_in_parallel.

      * When aff_inv is true, there is a main loop and two sequential
      parallel loops over the number of threads.

      Todos:

      * Figure out what to do with the message vector which is
      commented out in both versions

      * The main loop with the affine-invariant sampling could be
      modified with a new inner loop to do many function evaluations
      for each thread. However, I think this would demand combining
      the two sequential parallel loops.

      ---------------------------------------------------------

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

      If \ref aff_inv is set to true, then affine-invariant sampling
      is used. For affine-invariant sampling, the variable \ref
      step_fac represents the value of \f$ a \f$, the limits of the
      distribution for \f$ z \f$ (which defaults to 2). If \ref
      aff_inv is true and an initial point fails, then \ref mcmc()
      chooses random points inside the hypercube to attempt to find
      enough initial points to proceed. This method of finding initial
      points, however, is often too slow for large parameter spaces.

      Affine-invariant sampling works best when the number of walkers
      is much larger than the number of parameters. If \ref n_walk is
      0 or 1, then this class automatically sets \ref n_walk to three
      times the number of parameters. This class will otherwise allow
      the user to set a smaller number of walkers than parameters
      without notifying the user.

      In order to store data at each point, the user can store this
      data in any object of type \c data_t . If affine-invariant
      sampling is used, then each chain has it's own data object. The
      class keeps twice as many copies of these data object as would
      otherwise be required, in order to avoid copying of data objects
      in the case that the steps are accepted or rejected.

      Whether or not \ref aff_inv is true, there is a virtual function
      called \ref outside_parallel() which is called during the MCMC.
      Class descendants can replace this function with code which must
      be run outside of an OpenMP parallel region. Note that this is
      not implemented via, e.g. an OpenMP CRITICAL region so that no
      loss of performance is expected. If \ref aff_inv is false, then
      \ref outside_parallel() is called every \ref steps_in_parallel
      MCMC steps (for each OpenMP thread). If \ref aff_inv is true,
      then \ref outside_parallel() is called after all the walkers
      have completed for each thread.

      <b>Verbose output:</b> If verbose is 0, no output is generated
      (the default). If verbose is 1, then output to <tt>cout</tt>
      occurs only if the settings are somehow misconfigured and the
      class attempts to recover from them, for example if not enough
      functions are specified for the requested number of OpenMP
      threads, or if more than one thread was requested but
      O2SCL_SET_OPENMP was not defined, or if a negative value for \ref
      step_fac was requested. When verbose is 1, a couple messages are
      written to \ref scr_out : a summary of the number
      of walkers, chains, and threads at the beginning of the MCMC
      simulation, a message indicating why the MCMC simulation
      stopped, a message when the warm up iterations are completed, a
      message every time files are written to disk, and a message at
      the end counting the number of acceptances and rejections.
      If verbose is 2, then the file prefix is output to <tt>cout</tt>
      during initialization.

      \note This class is experimental.

      \future There is a little code in mcmc_init() and mcmc_cleanup()
      and I should document why that code needs to be there.

      \note Currently, this class requires that the data_t 
      has a good copy constructor. 

      \future The copy constructor for the data_t type is used when
      the user doesn't specify enough initial points for the
      corresponding number of threads and walkers. 
  */
  template<class func_t, class measure_t,
           class data_t, class vec_t=ubvector> class mcmc_para_base {
    
    protected:
  
    /// \name MPI properties
    //@{
    /// The MPI processor rank
    int mpi_rank;

    /// The MPI number of processors
    int mpi_size;
    //@}
  
    /// The screen output file
    std::ofstream scr_out;
  
    /// Random number generators, one for each thread
    std::vector<rng<> > rg;
  
    /// Pointer to proposal distribution for each thread
    std::vector<o2scl::prob_cond_mdim<vec_t> *> prop_dist;
  
    /** \brief If true, then use the user-specified proposal
        distribution (default false) */
    bool pd_mode;

    /// If true, we are in the warm up phase (default false)
    bool warm_up;

    /** \brief Current points in parameter space for each walker and 
        each OpenMP thread

        This is an array of size \ref n_threads times \ref n_walk initial
        guesses, indexed by <tt>thread_index*n_walk+walker_index</tt> .
    */
    std::vector<vec_t> current;

    /** \brief Data switch array for each walker and each OpenMP thread

        This is an array of size \ref n_threads times \ref n_walk initial
        guesses, indexed by <tt>thread_index*n_walk+walker_index</tt> .
    */
    std::vector<bool> switch_arr;
  
    /** \brief Return value counters, one vector for each independent
        chain
     */
    std::vector<std::vector<size_t> > ret_value_counts;
  
    /// \name Interface customization
    //@{
    /** \brief Initializations before the MCMC 
     */
    virtual int mcmc_init() {

      if (verbose>1) {
        std::cout << "Prefix is: " << prefix << std::endl;
      }

      if (verbose>0) {
        // Open main output file for this rank
        scr_out.open((prefix+"_"+
                      o2scl::itos(mpi_rank)+"_scr").c_str());
        scr_out.setf(std::ios::scientific);
      }
    
      // End of mcmc_init()
      return 0;
    }

    /** \brief Cleanup after the MCMC
     */
    virtual void mcmc_cleanup() {
      if (verbose>0) {
        for(size_t it=0;it<n_threads;it++) {
          scr_out << "mcmc (" << it << "," << mpi_rank
                  << "): accept=" << n_accept[it]
                  << " reject=" << n_reject[it] << std::endl;
        }
        scr_out.close();
      }
      return;
    }
    //@}

    /** \brief Index of the current walker, one for each OpenMP thread
      
        This quantity has to be a vector because different threads
        may have different values for the current walker during
        the initialization phase for the affine sampling algorithm.
    */
    std::vector<size_t> curr_walker;

    public:

    /** \brief If true, call the measurement function for the
        initial point
    */
    bool meas_for_initial;
  
    /// Integer to indicate completion
    static const int mcmc_done=-10;

    /// Integer to indicate rejection
    static const int mcmc_skip=-20;

    /// \name Output quantities
    //@{
    /** \brief The number of Metropolis steps which were accepted in 
        each independent chain (summed over all walkers)

        This vector has a size equal to \ref n_threads .
    */
    std::vector<size_t> n_accept;
  
    /** \brief The number of Metropolis steps which were rejected in 
        each independent chain (summed over all walkers)

        This vector has a size equal to \ref n_threads .
    */
    std::vector<size_t> n_reject;
    //@}

    /// \name Settings
    //@{
    /** \brief The MPI starting time (defaults to 0.0)

        This can be set by the user before mcmc() is called, so
        that the time required for initializations before
        the MCMC starts can be counted.
    */
    double mpi_start_time;

    /** \brief If non-zero, the maximum number of MCMC iterations 
        (default 0)

        If both \ref max_iters and \ref max_time are nonzero, the
        MCMC will stop when either the number of iterations 
        exceeds \ref max_iters or the time exceeds \ref max_time,
        whichever occurs first.
    */
    size_t max_iters;
  
    /** \brief Time in seconds (default is 0)

        If both \ref max_iters and \ref max_time are nonzero, the
        MCMC will stop when either the number of iterations 
        exceeds \ref max_iters or the time exceeds \ref max_time,
        whichever occurs first.
    */
    double max_time;

    /** \brief Prefix for output filenames (default "mcmc")
     */
    std::string prefix;

    /// If true, use affine-invariant Monte Carlo (default false)
    bool aff_inv;
  
    /// Stepsize factor (default 10.0)
    double step_fac;
  
    /// Optionally specify step sizes for each parameter
    std::vector<double> step_vec;

    /** \brief If true, couple the walkers across threads during
        affine-invariant sampling (default false)
    */
    bool couple_threads;

    /** \brief Number of warm up steps (successful steps not
        iterations) (default 0)
        
        \note Not to be confused with <tt>warm_up</tt>, which is 
        a protected boolean local variable in some functions which
        indicates whether we're in warm up mode or not.
    */
    size_t n_warm_up;

    /** \brief If non-zero, use as the seed for the random number 
        generator (default 0)

        The random number generator is modified so that each thread and
        each rank has a different set of random numbers.

        If this value is zero, then the random number generators are
        seeded by the clock time in seconds, so that if two separate
        simulations begin at the same time (to within 1 second) they
        will produce identical results. This can be avoided simply by
        ensuring that user_seed is different between the two
        simulations.
    */
    int user_seed;

    /// Output control (default 0)
    int verbose;

    /** \brief Maximum number of failed steps when generating initial points
        with affine-invariant sampling (default 1000)
    */
    size_t max_bad_steps;

    /** \brief Number of walkers (per openMP thread) for
        affine-invariant MC or 1 otherwise (default 1)
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
      verbose=1;
      warm_up=false;
      max_bad_steps=1000;

      always_accept=false;
      ai_initial_step=0.1;

      n_threads=1;
      n_walk=1;

      // Initial values for MPI paramers
      mpi_size=1;
      mpi_rank=0;
      mpi_start_time=0.0;

      #ifdef O2SCL_MPI
        // Get MPI rank, etc.
        MPI_Comm_rank(MPI_COMM_WORLD,&this->mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD,&this->mpi_size);
      #endif
    
      prefix="mcmc";
      max_time=0.0;
      max_iters=0;
      meas_for_initial=true;
      couple_threads=false;
      steps_in_parallel=100;
    }

    /// Number of OpenMP threads
    size_t n_threads;
  
    /** \brief Initial points in parameter space

        To fully specify all of the initial points, this should be 
        a vector of size \ref n_walk times \ref n_threads . Initial
        points are used for multiple threads and/or walkers if the
        full number of initial points is not specified.
    */
    std::vector<ubvector> initial_points;

    /** \brief The number of steps in parallel (default 100)
     */
    size_t steps_in_parallel;

    /** \brief Function outside parallel region
     */
    virtual void outside_parallel() {
      return;
    }

    // HMC Trajectory length
    int traj_length;

    /// \name Basic usage
    //@{
    /** \brief Perform a MCMC simulation

        Perform a MCMC simulation over \c n_params parameters starting
        at initial point \c init, limiting the parameters to be between
        \c low and \c high, using \c func as the objective function and
        calling the measurement function \c meas at each MC point.

        The vector \c data should be of size
        <tt>2*n_walk*n_threads</tt>.
    */
    virtual int mcmc(size_t n_params, vec_t &low, vec_t &high,
                  std::vector<func_t> &func,
                  std::vector<measure_t> &meas,
                  std::vector<data_t> &data) {
      
      // ----------------------------------------------------------------
      // Check inputs and settings, and fix them if necessary
      if (func.size()==0 || meas.size()==0) {
        O2SCL_ERR2("Size of 'func' or 'meas' array is zero in ",
                   "mcmc_para::mcmc().",o2scl::exc_einval);
      }
      if (func.size()<n_threads) {
        if (verbose>0) {
          std::cout << "mcmc_para::mcmc(): Not enough functions for "
                    << n_threads << " threads. Setting n_threads to "
                    << func.size() << "." << std::endl;
        }
        n_threads=func.size();
      }
      if (meas.size()<n_threads) {
        if (verbose>0) {
          std::cout << "mcmc_para::mcmc(): Not enough measurement "
                    << "objects for "
                    << n_threads << " threads. Setting n_threads to "
                    << meas.size() << "." << std::endl;
        }
        n_threads=meas.size();
      } // End of checking inputs and settings
      
      // ----------------------------------------------------------------
      // Set start time if necessary
      if (mpi_start_time==0.0) {
        #ifdef O2SCL_MPI
          mpi_start_time=MPI_Wtime();
        #else
          mpi_start_time=time(0);
        #endif
      } // End of setting start time

      // ----------------------------------------------------------------
      // Check intitial points
      
      // If not specified, setup initial guess from midpoint
      if (initial_points.size()==0) {
        initial_points.resize(1);
        initial_points[0].resize(n_params);
        for(size_t k=0;k<n_params;k++) {
          initial_points[0][k]=(low[k]+high[k])/2.0;
        } 
      } // End of not specified initial points
      
      // If specified, check existence and validity of initial points
      else {
        for(size_t iip=0;iip<initial_points.size();iip++) {
          
          // Check existence: Initial guess exists for all points
          if (initial_points[iip].size()<n_params) {
            std::cerr << "Not enough parameters." << std::endl;
            std::cerr << "On initial point " << iip << " of "
                      << initial_points.size() << "." << std::endl;
            std::cerr << "Expecting size " << n_params
                      << " and got size "
                      << initial_points[iip].size() << "." << std::endl;
            O2SCL_ERR2("Initial point vector not correctly sized ",
                       "in mcmc_base::mcmc().",o2scl::exc_efailed);
          } // End of check existence

          // Check validity (a): Initial points are within their limits
          for(size_t ipar=0;ipar<n_params;ipar++) {
            if (initial_points[iip][ipar]<low[ipar] ||
                initial_points[iip][ipar]>high[ipar]) {
              std::cout << "Parameters for point " << iip+1 << " of "
                        << initial_points.size() << "." << std::endl;
              for(size_t iki=0;iki<n_params;iki++) {
                std::cout << iki << " " << low[iki] << " "
                          << initial_points[iip][iki] << " "
                          << high[iki];
                if (initial_points[iip][ipar]<low[ipar]) {
                  std::cout << " L";
                } else if (initial_points[iip][ipar]<high[ipar]) {
                  std::cout << " H";
                }
                std::cout << std::endl;
              }
              O2SCL_ERR((((std::string)"Parameter ")+o2scl::szttos(ipar)+
                         " of "+o2scl::szttos(n_params)+
                         " out of range (value="+
                         o2scl::dtos(initial_points[iip][ipar])+
                         " low="+o2scl::dtos(low[ipar])+" high="+
                         o2scl::dtos(high[ipar])+
                         ") in mcmc_base::mcmc().").c_str(),
                        o2scl::exc_einval);
            }
          } // End of check check validity (a)

          // Check validity (b): Initial points are finite
          for (size_t ipar=0;ipar<initial_points[iip].size();ipar++) {
            if (!std::isfinite(initial_points[iip][ipar])) {
              O2SCL_ERR2("Initial point not finite in ",
                         "mcmc_para::mcmc().",o2scl::exc_einval);
            }
          } // End of check validity (b)
        } // End of if specified initial points
      } // End of check initial points

      // ----------------------------------------------------------------
      // Set number of OpenMP threads
      #ifdef O2SCL_SET_OPENMP
        omp_set_num_threads(n_threads);
      #else
        if (n_threads>1) {
          std::cout << "mcmc_para::mcmc(): "
                    << n_threads << " threads were requested but the "
                    << "-DO2SCL_SET_OPENMP flag was not used during "
                    << "compilation. Setting n_threads to 1."
                    << std::endl;
          n_threads=1;
        }
      #endif // End of setting OpenMP threads

      // ----------------------------------------------------------------
      // Fix 'step_fac' if less than or equal to zero
      if (step_fac<=0.0) {
        if (aff_inv) {
          std::cout << "mcmc_para::mcmc(): Requested negative or zero "
                    << "step_fac with aff_inv=true.\nSetting to 2.0."
                    << std::endl;
          step_fac=2.0;
        } else {
          std::cout << "mcmc_para::mcmc(): Requested negative or zero "
                    << "step_fac. Setting to 10.0." << std::endl;
          step_fac=10.0;
        }
      } // End of fixing 'step_fac'

      // ----------------------------------------------------------------
      // Set RNGs with a different seed for each thread and rank. 
      rg.resize(n_threads);
      unsigned long int seed=time(0);
      if (this->user_seed!=0) {
        seed=this->user_seed;
      }
      for(size_t it=0;it<n_threads;it++) {
        seed*=(mpi_rank*n_threads+it+1);
        rg[it].set_seed(seed);
      } // End of setting RNGs

      // ----------------------------------------------------------------
      // Keep track of successful and failed MH moves in each
      // independent chain
      n_accept.resize(n_threads);
      n_reject.resize(n_threads);
      for(size_t it=0;it<n_threads;it++) {
        n_accept[it]=0;
        n_reject[it]=0;
      } // End of tracking successful and failed MH moves

      // ----------------------------------------------------------------
      // Warm-up flag, not to be confused with 'n_warm_up', the
      // number of warm_up iterations.
      warm_up=true;
      if (n_warm_up==0) warm_up=false;

      // ----------------------------------------------------------------
      // Allocate space for storage
      
      // Set storage of return values from each thread
      std::vector<int> func_ret(n_threads), meas_ret(n_threads);

      // Set required storage size
      size_t ssize=n_walk*n_threads;

      // Allocate current point and current weight
      current.resize(ssize);
      std::vector<double> w_current(ssize);
      for(size_t i=0;i<ssize;i++) {
        current[i].resize(n_params);
        w_current[i]=0.0;
      }

      // Allocate curr_walker
      curr_walker.resize(n_threads);

      // Allocation of ret_value_counts should be handled by the user
      // in mcmc_init(), because this class can't determine what the
      // possible and interesting return values are.
      // ret_value_counts.resize(n_threads);

      // Initialize data and switch arrays
      switch_arr.resize(ssize);
      for(size_t i=0;i<switch_arr.size();i++) switch_arr[i]=false;
    
      // Next point and next weight for each thread
      std::vector<vec_t> next(n_threads);
      for(size_t it=0;it<n_threads;it++) {
        next[it].resize(n_params);
      }
      std::vector<double> w_next(n_threads);

      // Best point and best weight for each thread (only used when
      // aff_inv=false and not used until after the initial points are
      // computed)
      std::vector<vec_t> best_t(n_threads);
      for(size_t it=0;it<n_threads;it++) {
        best_t[it].resize(n_params);
      }
      std::vector<double> w_best_t(n_threads);
      
      // Best point over all threads
      vec_t best(n_params);
      double w_best;

      // Generally, these flags are are true for any thread if func_ret
      // or meas_ret is equal to mcmc_done.
      std::vector<bool> mcmc_done_flag(n_threads);
      for(size_t it=0;it<n_threads;it++) {
        mcmc_done_flag[it]=false;
      }
          
      // Proposal weight
      std::vector<double> q_prop(n_threads);

      // ----------------------------------------------------------------
      // Run the mcmc_init() function. 
      int init_ret=mcmc_init();
      if (init_ret!=0) {
        O2SCL_ERR("Function mcmc_init() failed in mcmc_base::mcmc().",
                  o2scl::exc_einval);
        return init_ret;
      }

      // ----------------------------------------------------------------
      // ----------------------------------------------------------------

      std::vector<double> q=position;
      
      std::default_random_engine gen;
      std::normal_distribution<double> dist(0.0, 1.0);
      
      std::vector<double> p(q.size(), dist(gen));
      momentum=p;

      // Make a half step for momentum at the beginning
      for (size_t i=0; i<p.size(); i++) {
        p[i]=p[i]-0.5*step_size[i]*grad_potential(q)[i];
      }

      // Alternate full steps for position and momentum
      for (int j=1; j<=traj_length; j++) {
        
        // Make a full step for the position
        for (size_t i=0; i<q.size(); i++) {
          q[i]=q[i]+step_size[i]*p[i];
          
          // Make a full step for the momentum except at the end of trajectory
          if (j!=traj_length) {
            p[i]=p[i]-step_size[i]*grad_potential(q)[i];
          }
        }
      }
      // Make a half step for momentum at the end
      for (size_t i=0; i<p.size(); i++) {
        p[i]=p[i]-0.5*step_size[i]*grad_potential(q)[i];
      }

      // Negate momentum at end of trajectory to make the proposal symmetric
      for (size_t i=0; i<p.size(); i++) {
        p[i]=-p[i];
      }

      // Evaluate potential and kinetic energies at start and end of trajectory
      double u_curr=potential_en(position);
      double u_prop=potential_en(q);
      double k_curr=kinetic_en(momentum);
      double k_prop=kinetic_en(p);

      // Accept or reject the state at end of trajectory, returning either
      // the position at the end of the trajectory or the initial position
      if (dist(gen) < exp(u_curr-u_prop+k_curr-k_prop)) {

        // Accept: Update current positions
        position=q;
        accept++;
      }
      else {

        // Reject: Keep initial positions
        reject++;
      }

      return 0;
    }
  }; // End of class mcmc_para_base
} // End of namespace o2scl

#endif // O2SCL_MCMC_PARA_H