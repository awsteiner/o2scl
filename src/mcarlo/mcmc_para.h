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
#ifndef O2SCL_MCMC_PARA_H
#define O2SCL_MCMC_PARA_H

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

#include <o2scl/hdf_file.h>
#include <o2scl/exception.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/vector.h>
#include <o2scl/multi_funct.h>
#include <o2scl/vec_stats.h>
#include <o2scl/cli.h>

namespace o2scl {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

  // Do we need one stepper for each thread? maybe not.

  /** \brief Desc
   */
  template<class func_t, class data_t, class vec_t>
  class mcmc_stepper_base {
    
  public:

    /** \brief Desc
     */
    rng<> rg;

    mcmc_stepper_base() {
      unsigned long int seed=time(0);
      rg.set_seed(seed);
    }

    /** \brief Desc
     */
    virtual void step(size_t n_params, func_t &f, vec_t &current, vec_t &next, 
                      double w_current, double &w_next, vec_t &low,
                      vec_t &high, int &func_ret, bool &accept,
                      data_t &dat) {
      O2SCL_ERR("Parent class.",o2scl::exc_eunimpl);
      return;
    }
    
    virtual ~mcmc_stepper_base() {
    }
    
  };

#ifdef O2SCL_NEVER_DEFINED

  /** \brief Desc
   */
  template<class func_t, class data_t, class vec_t> mcmc_stepper_mh :
  public<func_t,data_t,vec_t> {
    
  public:

    /** \brief Desc
     */
    double step_fac;

    mcmc_stepper_mh() {
      step_fac=2.0;
    }
    
    virtual ~mcmc_stepper_mh() {
    }
    
    /** \brief Desc
     */
    virtual void step(size_t n_params, func_t &f, vec_t &current, vec_t &next,
                      double w_curr, double &w_next, vec_t &low, vec_t &high, 
                      int &func_ret, bool &accept, data_t &dat) {
      
      for(size_t k=0;k<n_params;k++) {
        next[k]=current[k]+(rg.random()*2.0-1.0)*
          (high[k]-low[k])/step_fac;
      }

      func_ret=success;
      
      for(size_t k=0;k<n_params;k++) {
        if (next[k]<low[k] || next[k]>high[k]) {
          func_ret=mcmc_skip;
        }
      }

      func_ret=f(n_params,next,w_next,dat);

      double r=rg.random();
            
      // Metropolis algorithm
      accept=false;
      if (r<exp(w_next-w_current)) {
        accept=true;
      }
      
      return;
    }
    
  };

  /** \brief Desc
   */
  template<class grad_t, class func_t, class data_t,
           class vec_t> mcmc_stepper_hmc :
  public<func_t,data_t,vec_t> {
    
  public:

    /** \brief Desc
     */
    int traj_length;

    /** \brief Desc
     */
    double inv_mass;

    /** \brief Desc
     */
    prob_dens_gaussian pdg;

    /** \brief Desc
     */
    vec_t mom;

    /** \brief Desc
     */
    double mom_step;

    /** \brief Desc
     */
    grad_t g;
    
    mcmc_stepper_hmc() {
      
      unsigned long int seed=time(0);
      rg.set_seed(seed);
      inv_mass=0.1;
      traj_length=20;
      mom_step=0.18;
    }

    virtual ~mcmc_stepper_hmc() {
    }
    
    /** \brief Desc
     */
    virtual void step(size_t n_params, func_t &f, vec_t &current, vec_t &next,
                      double w_curr, double &w_next, vec_t &low, vec_t &high, 
                      int &func_ret, bool &accept, data_t &dat) {

      vec_t mom, grad, mom_next;
      
      for(size_t k=0;k<n_params;k++) {
        mom[k]=pdg()*mom_step;
      }

      g(n_params,current,grad);
      
      for(size_t k=0;k<n_params;k++) {
        mom_next[k]=mom[k]-0.5*mom_step*grad[k];
      }

      for(size_t i=0;i<traj_length;i++) {

        for(size_t k=0;k<n_params;k++) {
          next[k]=current[k]+mom_step*mom_next[k];
        }
        
        func_ret=f(n_params,next,w_next,dat);
        g(n_params,next,grad);

        if (i<traj_length-1) {
          for(size_t k=0;k<n_params;k++) {
            mom_next[k]=mom_next[k]-mom_step*grad[k];
          }

        }
        
      }
      
      for(size_t k=0;k<n_params;k++) {
        mom_next[k]-=0.5*mom_step*grad[k];
      }

      double pot_curr=-log(0.5*w_current);
      double pot_next=-log(0.5*w_next);

      double kin_curr=0.0, kin_next=0.0;
      for(size_t k=0;k<n_params;k++) {
        kin_curr+=mom[k]*mass_inv*mom[k]/2.0;
        kin_next+=mom_next[k]*mass_inv*mom_next[k]/2.0;
      }
        
      double r=rg.random();
            
      // Metropolis algorithm
      accept=false;
      if (r<exp(pot_curr-pot_next+kin_curr-kin_next)) {
        accept=true;
      }

      return;
    }
    
  };

#endif
  
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
           class data_t, class vec_t=ubvector,
           class stepper_t=
           mcmc_stepper_base<func_t,data_t,vec_t>
           > class mcmc_para_base {
    
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

    /// If true, use new initial point code
    bool new_ip;

    /// If true, use the new stepper
    bool new_step;

    /// The stepper
    stepper_t stepper;
    
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
      new_ip=false;
      new_step=false;
      
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

      // Verify that the input and settings make sense and fix
      // them if we can

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
      }
      if (data.size()<2*n_walk*n_threads) {
        std::cout << "mcmc_para::mcmc() data.size(): " << data.size()
                  << " n_walk: " << n_walk << " threads: "
                  << n_threads << std::endl;
        O2SCL_ERR2("Not enough data objects for walkers and threads in ",
                   "mcmc_para::mcmc()",o2scl::exc_einval);
      }
      if (aff_inv==true && n_walk<=1) {
        if (verbose>0) {
          std::cout << "mcmc_para::mcmc(): Affine-invariant "
                    << "sampling selected "
                    << "but n_walk is <= 1.\n  Setting n_walk to 3 * n_params."
                    << std::endl;
        }
        n_walk=3*n_params;
      }
      if (pd_mode && aff_inv) {
        O2SCL_ERR2("Using a proposal distribution with affine-invariant ",
                   "sampling not implemented in mcmc_para::mcmc().",
                   o2scl::exc_eunimpl);
      }

      // Set start time if necessary
      if (mpi_start_time==0.0) {
#ifdef O2SCL_MPI
        mpi_start_time=MPI_Wtime();
#else
        mpi_start_time=time(0);
#endif
      }

      if (initial_points.size()==0) {
      
        // Setup initial guess from midpoint if not specified
        initial_points.resize(1);
        initial_points[0].resize(n_params);
        for(size_t k=0;k<n_params;k++) {
          initial_points[0][k]=(low[k]+high[k])/2.0;
        }
      
      } else {
      
        // If initial points are specified, make sure they're within
        // the user-specified limits
        for(size_t iip=0;iip<initial_points.size();iip++) {
          if (initial_points[iip].size()<n_params) {
            std::cerr << "Not enough parameters." << std::endl;
            std::cerr << "On initial point " << iip << " of "
                      << initial_points.size() << "." << std::endl;
            std::cerr << "Expecting size " << n_params
                      << " and got size "
                      << initial_points[iip].size() << "." << std::endl;
            O2SCL_ERR2("Initial point vector not correctly sized ",
                       "in mcmc_base::mcmc().",o2scl::exc_efailed);
          }
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
          }
        }
      
        // Double check that the initial points are distinct and finite
        for(size_t i=0;i<initial_points.size();i++) {
          for(size_t k=0;k<initial_points[i].size();k++) {
            if (!std::isfinite(initial_points[i][k])) {
              O2SCL_ERR2("Initial point not finite in ",
                         "mcmc_para::mcmc().",o2scl::exc_einval);
            }
          }
          // 2/14/19 - AWS: I waffle on whether or not this ought to be
          // included, but it's too confusing and I'm having too much
          // trouble with it right now so I'm taking it out for now.
          if (false) {
            for(size_t j=i+1;j<initial_points.size();j++) {
              bool vec_equal=true;
              for(size_t k=0;k<initial_points[i].size();k++) {
                if (initial_points[i][k]!=initial_points[j][k]) {
                  vec_equal=false;
                }
              }
              if (vec_equal) {
                std::cout.setf(std::ios::scientific);
                std::cout << i << " ";
                o2scl::vector_out(std::cout,initial_points[i],true);
                std::cout << j << " ";
                o2scl::vector_out(std::cout,initial_points[j],true);
                O2SCL_ERR2("Initial points not distinct in ",
                           "mcmc_para::mcmc().",o2scl::exc_einval);
              }
            }
          }
        }
      
      }

      // Set number of threads
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
#endif

      // Storage for return values from each thread
      std::vector<int> func_ret(n_threads), meas_ret(n_threads);
      
      // Fix 'step_fac' if it's less than or equal to zero
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
      }
    
      // Set RNGs with a different seed for each thread and rank. 
      rg.resize(n_threads);
      unsigned long int seed=time(0);
      if (this->user_seed!=0) {
        seed=this->user_seed;
      }
      for(size_t it=0;it<n_threads;it++) {
        seed*=(mpi_rank*n_threads+it+1);
        rg[it].set_seed(seed);
      }
    
      // Keep track of successful and failed MH moves in each
      // independent chain
      n_accept.resize(n_threads);
      n_reject.resize(n_threads);
      for(size_t it=0;it<n_threads;it++) {
        n_accept[it]=0;
        n_reject[it]=0;
      }

      // Warm-up flag, not to be confused with 'n_warm_up', the
      // number of warm_up iterations.
      warm_up=true;
      if (n_warm_up==0) warm_up=false;

      // Storage size required
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
    
      // --------------------------------------------------------------
      // Run the mcmc_init() function. 
    
      int init_ret=mcmc_init();
      if (init_ret!=0) {
        O2SCL_ERR("Function mcmc_init() failed in mcmc_base::mcmc().",
                  o2scl::exc_einval);
        return init_ret;
      }

      // --------------------------------------------------------------
      // Initial verbose output (note that scr_out isn't created until
      // the mcmc_init() function call above.
      
      if (verbose>=1) {
        if (aff_inv) {
          scr_out << "mcmc_para_base::mcmc(): "
                  << "Affine-invariant step, n_params="
                  << n_params << ", n_walk=" << n_walk
                  << ", n_threads=" << n_threads << ",\n  rank="
                  << mpi_rank << ", n_ranks="
                  << mpi_size << std::endl;
        } else if (pd_mode==true) {
          scr_out << "mcmc_para_base::mcmc(): "
                  << "With proposal distribution, n_params="
                  << n_params << ", n_threads=" << n_threads << ", rank="
                  << mpi_rank << ", n_ranks="
                  << mpi_size << std::endl;
        } else {
          scr_out << "mcmc_para_base::mcmc(): "
                  << "Random-walk w/uniform dist., n_params="
                  << n_params << ", n_threads=" << n_threads << ", rank="
                  << mpi_rank << ", n_ranks="
                  << mpi_size << std::endl;
        }
        scr_out << "Set start time to: " << mpi_start_time << std::endl;
      }

      // --------------------------------------------------------
      // Initial point and weights for affine-invariant sampling
    
      if (aff_inv) {
      
#ifdef O2SCL_SET_OPENMP
#pragma omp parallel default(shared)
#endif
        {
#ifdef O2SCL_SET_OPENMP
#pragma omp for
#endif
          for(size_t it=0;it<n_threads;it++) {

            // Initialize each walker in turn
            for(curr_walker[it]=0;curr_walker[it]<n_walk &&
                  mcmc_done_flag[it]==false;curr_walker[it]++) {

              // Index in storage
              size_t sindex=n_walk*it+curr_walker[it];

              // Index in initial_point specification
              size_t ip_index=sindex % initial_points.size();
            
              size_t init_iters=0;
              bool done=false;

              // If we already have a unique guess for this
              // walker/thread, try to use that
            
              if (sindex<initial_points.size()) {

                // Copy from the initial points array
                for(size_t ipar=0;ipar<n_params;ipar++) {
                  current[sindex][ipar]=initial_points[ip_index][ipar];
                }
              
                // Compute the weight
                func_ret[it]=func[it](n_params,current[sindex],
                                      w_current[sindex],data[sindex]);

                if (func_ret[it]==mcmc_done) {
                  mcmc_done_flag[it]=true;
                } else if (func_ret[it]==o2scl::success) {

                  // If we have a good point, update ret_value_counts
                  // and call the measurement function
                  if (func_ret[it]>=0 && ret_value_counts.size()>it && 
                      func_ret[it]<((int)ret_value_counts[it].size())) {
                    ret_value_counts[it][func_ret[it]]++;
                  }
                  if (meas_for_initial) {
                    meas_ret[it]=meas[it](current[sindex],w_current[sindex],
                                          curr_walker[it],func_ret[it],
                                          true,data[sindex]);
                  } else {
                    meas_ret[it]=0;
                  }
                  if (meas_ret[it]==mcmc_done) {
                    mcmc_done_flag[it]=true;
                  }
                  done=true;
                }
              }
            
              // Otherwise, if the initial guess wasn't provided or
              // failed for some reason, try generating a new point
            
              while (!done && !mcmc_done_flag[it]) {

                // Make a perturbation from the initial point
                for(size_t ipar=0;ipar<n_params;ipar++) {
                  do {
                    current[sindex][ipar]=
                      initial_points[ip_index][ipar]+
                      (rg[it].random()*2.0-1.0)*
                      (high[ipar]-low[ipar])*ai_initial_step;
                  } while (current[sindex][ipar]>high[ipar] ||
                           current[sindex][ipar]<low[ipar]);
                }
              
                // Compute the weight
                func_ret[it]=func[it](n_params,current[sindex],
                                      w_current[sindex],data[sindex]);
                
                // ------------------------------------------------
              
                // Increment iteration count
                init_iters++;
              
                if (func_ret[it]==mcmc_done) {
                  mcmc_done_flag[it]=true;
                } else {
                  // If we have a good point, update ret_value_counts,
                  // call the measurement function and stop the loop
                  if (func_ret[it]==o2scl::success) {
                    if (verbose>=2) {
                      scr_out << "Found initial guess for thread "
                              << it << ". func_ret,weight,params=\n  "
                              << func_ret[it] << " "
                              << w_current[sindex] << " ";
                      for(size_t iji=0;iji<n_params;iji++) {
                        scr_out << current[sindex][iji] << " ";
                      }
                      scr_out << std::endl;
                    }
                    if (func_ret[it]>=0 && ret_value_counts.size()>it && 
                        func_ret[it]<((int)ret_value_counts[it].size())) {
                      ret_value_counts[it][func_ret[it]]++;
                    }
                    if (meas_ret[it]!=mcmc_done) {
                      if (meas_for_initial) {
                        meas_ret[it]=meas[it](current[sindex],w_current[sindex],
                                              curr_walker[it],func_ret[it],true,
                                              data[sindex]);
                      } else {
                        meas_ret[it]=0;
                      }
                    } else {
                      mcmc_done_flag[it]=true;
                    }
                    done=true;
                  } else if (init_iters>max_bad_steps) {
                    std::string err=((std::string)"In loop with thread ")+
                      o2scl::szttos(it)+" iterations required to obtain an "+
                      "initial point exceeded "+o2scl::szttos(max_bad_steps)+
                      " in mcmc_para_base::mcmc().";
                    O2SCL_ERR(err.c_str(),o2scl::exc_einval);
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
              scr_out << "mcmc (" << it << "," << mpi_rank
                      << "): Returned mcmc_done "
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
          for(curr_walker[it]=0;curr_walker[it]<n_walk;
              curr_walker[it]++) {
            size_t sindex=n_walk*it+curr_walker[it];
            if (w_current[sindex]>w_current[0]) {
              best_index=sindex;
              w_best=w_current[sindex];
            }
          }
        }
        best=current[best_index];

        // Verbose output
        if (verbose>=2) {
          for(size_t it=0;it<n_threads;it++) {
            for(curr_walker[it]=0;curr_walker[it]<n_walk;curr_walker[it]++) {
              size_t sindex=n_walk*it+curr_walker[it];
              scr_out.precision(4);
              scr_out << "mcmc (" << it << "," << mpi_rank << "): i_walk: ";
              scr_out.width((int)(1.0+log10((double)(n_walk-1))));
              scr_out << curr_walker[it] << " log_wgt: "
                      << w_current[sindex] << " (initial; ai)" << std::endl;
              scr_out.precision(6);
            }
          }
        }

        // End of 'if (aff_inv)' for initial point evaluation
      
      } else {

        // --------------------------------------------------------
        // Initial point evaluation when aff_inv is false.

        if (new_ip) {
          
          size_t ip_size=initial_points.size();
          
          if (ip_size<n_threads) {
            O2SCL_ERR("mcmc_para, insufficient number of initial points.",
                      o2scl::exc_einval);
          }
          
#ifdef O2SCL_SET_OPENMP
#pragma omp parallel default(shared)
#endif
          {
#ifdef O2SCL_SET_OPENMP
#pragma omp for
#endif
            for(size_t it=0;it<n_threads;it++) {
              
              // Note that this value is used (e.g. in
              // mcmc_para_table::add_line() ) even if aff_inv is
              // false, so we set it to zero here.
              curr_walker[it]=0;
              
              // Copy from the initial points array into current point
              for(size_t ipar=0;ipar<n_params;ipar++) {
                current[it][ipar]=initial_points[it % ip_size][ipar];
              }
              
              // If we have a new unique initial point, then
              // perform a function evaluation
              func_ret[it]=func[it](n_params,current[it],w_current[it],
                                    data[it]);
            }
            
          }
          // End of parallel region
          
        } else {
          
#ifdef O2SCL_SET_OPENMP
#pragma omp parallel default(shared)
#endif
          {
#ifdef O2SCL_SET_OPENMP
#pragma omp for
#endif
            for(size_t it=0;it<n_threads;it++) {
              
              // Note that this value is used (e.g. in
              // mcmc_para_table::add_line() ) even if aff_inv is
              // false, so we set it to zero here.
              curr_walker[it]=0;
              
              // Copy from the initial points array into current point
              size_t ip_size=initial_points.size();
              for(size_t ipar=0;ipar<n_params;ipar++) {
                current[it][ipar]=initial_points[it % ip_size][ipar];
              }
              
              if (it<ip_size) {
                // If we have a new unique initial point, then
                // perform a function evaluation
                func_ret[it]=func[it](n_params,current[it],w_current[it],
                                      data[it]);
              } else {
                // Otherwise copy the result already computed
                func_ret[it]=func_ret[it % ip_size];
                w_current[it]=w_current[it % ip_size];
                // This loop requires the data to have a valid
                // copy constructor
                for(size_t j=0;j<data.size();j++) {
                  data[it]=data[it % ip_size];
                }
              }
              
            }
            
          }
          // End of parallel region
          
        }
        
        // Check return values from initial point function evaluations
        for(size_t it=0;it<n_threads;it++) {
          if (func_ret[it]==mcmc_done) {
            if (verbose>=1) {
              scr_out << "mcmc (" << it
                      << "): Initial point returned mcmc_done."
                      << std::endl;
            }
            mcmc_cleanup();
            return 0;
          }
          if (func_ret[it]!=o2scl::success) {
            if (err_nonconv) {
              O2SCL_ERR((((std::string)"Initial weight from thread ")+
                         o2scl::szttos(it)+
                         " vanished in mcmc_para_base::mcmc().").c_str(),
                        o2scl::exc_einval);
            }
            return 2;
          }
        }
        
        // --------------------------------------------------------
        // Post-processing initial point when aff_inv is false.
        
#ifdef O2SCL_SET_OPENMP
#pragma omp parallel default(shared)
#endif
        {
#ifdef O2SCL_SET_OPENMP
#pragma omp for
#endif
          for(size_t it=0;it<n_threads;it++) {
            
            size_t ip_size=initial_points.size();
            if (it>=ip_size) {
              // If no initial point was specified, copy one of
              // the other initial points
              func_ret[it]=func_ret[it % ip_size];
              current[it]=current[it % ip_size];
              w_current[it]=w_current[it % ip_size];
              data[it]=data[it % ip_size];
            }
            
            // Update the return value count
            if (func_ret[it]>=0 && ret_value_counts.size()>it && 
                func_ret[it]<((int)ret_value_counts[it].size())) {
              ret_value_counts[it][func_ret[it]]++;
            }
            if (meas_for_initial) {
              // Call the measurement function    
              meas_ret[it]=meas[it](current[it],w_current[it],0,
                                    func_ret[it],true,data[it]);
            } else {
              meas_ret[it]=0;
            }
            if (meas_ret[it]==mcmc_done) {
              mcmc_done_flag[it]=true;
            }
            
            // End of loop over threads
          }
          
        }
        // End of parallel region
        
        // Stop early if mcmc_done was returned from one of the
        // measurement function calls
        bool stop_early=false;
        for(size_t it=0;it<n_threads;it++) {
          if (mcmc_done_flag[it]==true) {
            if (verbose>=1) {
              scr_out << "mcmc (" << it << "," << mpi_rank
                      << "): Returned mcmc_done "
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
        for(size_t it=1;it<n_threads;it++) {
          if (w_current[it]<w_best) {
            best=current[it];
            w_best=w_current[it];
          }
        }
        
        if (verbose>=2) {
          scr_out.precision(4);
          scr_out << "mcmc: ";
          for(size_t it=0;it<n_threads;it++) {
            scr_out << w_current[it] << " ";
          }
          scr_out << " (initial)" << std::endl;
          scr_out.precision(6);
        }
        
        // End of initial point region for 'aff_inv=false'
      }
      
      // Set meas_for_initial back to true if necessary
      meas_for_initial=true;
      
      // --------------------------------------------------------
      // Require keypress after initial point if verbose is
      // sufficiently large.
      
      if (verbose>=4) {
        std::cout << "Initial point(s) done. "
                  << "Press a key and type enter to continue. ";
        char ch;
        std::cin >> ch;
      }
      
      // End of initial point and weight section
      // --------------------------------------------------------
      
      // The main section split into two parts, aff_inv=false and
      // aff_inv=true.
      
      if (aff_inv==false) {
        
        // ---------------------------------------------------
        // Start of main loop over threads for aff_inv=false
        
        bool main_done=false;
        
        // Initialize the number of iterations for each thread
        std::vector<size_t> mcmc_iters(n_threads);
        for(size_t it=0;it<n_threads;it++) {
          mcmc_iters[it]=0;
        }
        
        while (!main_done) {
          
#ifdef O2SCL_SET_OPENMP
#pragma omp parallel default(shared)
#endif
          {
#ifdef O2SCL_SET_OPENMP
#pragma omp for
#endif
            for(size_t it=0;it<n_threads;it++) {
              
              bool inner_done=false;
              while (!inner_done && !main_done) {
                
                if (verbose>=2) {
                  scr_out << "Iteration: " << mcmc_iters[it] << " of "
                          << max_iters << " thread " << it << " accept: "
                          << n_accept[it] << " " << n_reject[it]
                          << std::endl;
                }
                
                // ---------------------------------------------------
                // Select next point for aff_inv=false
                
                if (new_step) {
                  
                  bool accept_new;
                  size_t sindex2=n_walk*it+curr_walker[it];
                  
                  if (switch_arr[sindex2]==false) {
                    stepper.step(n_params,func[it],current[it],
                                 next[it],w_current[sindex2],w_next[it],
                                 low,high,func_ret[it],accept_new,
                                 data[sindex2+n_walk*n_threads]);
                  } else {
                    stepper.step(n_params,func[it],current[it],
                                 next[it],w_current[sindex2],w_next[it],
                                 low,high,func_ret[it],accept_new,
                                 data[sindex2]);
                  }
                  
                } else if (pd_mode) {
                  
                  // Use proposal distribution and compute associated weight
                  q_prop[it]=prop_dist[it]->log_metrop_hast(current[it],
                                                            next[it]);
              
                  if (!std::isfinite(q_prop[it])) {
                    O2SCL_ERR2("Proposal distribution not finite in ",
                               "mcmc_para_base::mcmc().",
                               o2scl::exc_efailed);
                  }
            
                } else {
            
                  // Uniform random-walk step
                  for(size_t k=0;k<n_params;k++) {
                    if (step_vec.size()>0) {
                      next[it][k]=current[it][k]+(rg[it].random()*2.0-1.0)*
                        step_vec[k%step_vec.size()];
                    } else {
                      next[it][k]=current[it][k]+(rg[it].random()*2.0-1.0)*
                        (high[k]-low[k])/step_fac;
                    }
                  }
            
                }   
          
                // ---------------------------------------------------
                // Compute next weight for aff_inv=false
          
                func_ret[it]=o2scl::success;

                // If the next point out of bounds, ensure that the point is
                // rejected without attempting to evaluate the function
                for(size_t k=0;k<n_params;k++) {
                  if (next[it][k]<low[k] || next[it][k]>high[k]) {
                    func_ret[it]=mcmc_skip;
                    if (verbose>=3) {
                      if (next[it][k]<low[k]) {
                        std::cout << "mcmc (" << it << ","
                                  << mpi_rank << "): Parameter with index "
                                  << k << " and value " << next[it][k]
                                  << " smaller than limit " << low[k]
                                  << std::endl;
                        scr_out << "mcmc (" << it << ","
                                << mpi_rank << "): Parameter with index " << k
                                << " and value " << next[it][k]
                                << " smaller than limit " << low[k]
                                << std::endl;
                      } else {
                        std::cout << "mcmc (" << it << "," << mpi_rank
                                  << "): Parameter with index " << k
                                  << " and value " << next[it][k]
                                  << " larger than limit " << high[k]
                                  << std::endl;
                        scr_out << "mcmc (" << it << "," << mpi_rank
                                << "): Parameter with index " << k
                                << " and value " << next[it][k]
                                << " larger than limit " << high[k]
                                << std::endl;
                      }
                    }
                  }
                }

                // Evaluate the function, set the 'done' flag if
                // necessary, and update the return value array
                if (func_ret[it]!=mcmc_skip) {
                  if (switch_arr[n_walk*it+curr_walker[it]]==false) {
                    func_ret[it]=func[it](n_params,next[it],w_next[it],
                                          data[it*n_walk+curr_walker[it]+
                                               n_walk*n_threads]);
                  } else {
                    func_ret[it]=func[it](n_params,next[it],w_next[it],
                                          data[it*n_walk+curr_walker[it]]);
                  }
                  if (func_ret[it]==mcmc_done) {
                    mcmc_done_flag[it]=true;
                  } else {
                    if (func_ret[it]>=0 && ret_value_counts.size()>it && 
                        func_ret[it]<((int)ret_value_counts[it].size())) {
                      ret_value_counts[it][func_ret[it]]++;
                    }
                  }
                }
            
                // ------------------------------------------------------
                // Accept or reject and call the measurement function for
                // aff_inv=false
            
                // Index in storage
                size_t sindex=n_walk*it+curr_walker[it];
            
                bool accept=false;
                if (always_accept && func_ret[it]==success) accept=true;
            
                if (func_ret[it]==o2scl::success) {
                  double r=rg[it].random();
            
                  if (pd_mode) {
                    if (r<exp(w_next[it]-w_current[sindex]+q_prop[it])) {
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
                      meas_ret[it]=meas[it](next[it],w_next[it],
                                            curr_walker[it],func_ret[it],true,
                                            data[sindex+n_threads*n_walk]);
                    } else {
                      meas_ret[it]=meas[it](next[it],w_next[it],
                                            curr_walker[it],func_ret[it],true,
                                            data[sindex]);
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
                    {
                      if (switch_arr[sindex]==false) {
                        meas_ret[it]=meas[it](next[it],w_next[it],
                                              curr_walker[it],
                                              func_ret[it],false,
                                              data[sindex+
                                                   n_threads*n_walk]);
                      } else {
                        meas_ret[it]=meas[it](next[it],w_next[it],
                                              curr_walker[it],
                                              func_ret[it],false,
                                              data[sindex]);
                      }
                    }
                  }

                }

                // ---------------------------------------------------
                // Best point, update iteration counts, and check if done

                // Collect best point
                if (func_ret[it]==o2scl::success && w_best_t[it]>w_next[it]) {
                  best_t[it]=next[it];
                  w_best_t[it]=w_next[it];
                }
              
                // Check to see if mcmc_done was returned or if meas_ret
                // returned an error
                if (meas_ret[it]==mcmc_done || func_ret[it]==mcmc_done) {
                  main_done=true;
                }
                if (meas_ret[it]!=mcmc_done &&
                    meas_ret[it]!=o2scl::success) {
                  if (err_nonconv) {
                    O2SCL_ERR((((std::string)"Measurement function ")+
                               "returned "+o2scl::dtos(meas_ret[it])+
                               " in mcmc_para_base::mcmc().").c_str(),
                              o2scl::exc_efailed);
                  }
                  main_done=true;
                }
            
                // Update iteration count and reset counters for
                // warm up iterations if necessary
                if (main_done==false) {

                  scr_out << "Incrementing mcmc_iters." << std::endl;
                  mcmc_iters[it]++;
              
                  if (warm_up && mcmc_iters[it]==n_warm_up) {
                    warm_up=false;
                    mcmc_iters[it]=0;
                    n_accept[it]=0;
                    n_reject[it]=0;
                    if (verbose>=1) {
                      scr_out << "o2scl::mcmc_para: Thread " << it
                              << " finished warmup." << std::endl;
                    }
                
                  }
                }
            
                // Stop this thread if iterations greater than max_iters
                if (inner_done==false && warm_up==false && max_iters>0 &&
                    mcmc_iters[it]==max_iters) {
                  if (verbose>=1) {
                    scr_out << "o2scl::mcmc_para: Thread " << it
                            << " stopping because number of iterations ("
                            << mcmc_iters[it] << ") equal to max_iters ("
                            << max_iters << ")." << std::endl;
                  }
                  inner_done=true;
                }
            
                // If we're out of time, stop all threads
                if (main_done==false) {
#ifdef O2SCL_MPI
                  double elapsed=MPI_Wtime()-mpi_start_time;
#else
                  double elapsed=time(0)-mpi_start_time;
#endif
                  if (max_time>0.0 && elapsed>max_time) {
                    if (verbose>=1) {
                      scr_out << "o2scl::mcmc_para: Thread " << it
                              << " stopping because elapsed (" << elapsed
                              << ") > max_time (" << max_time << ")."
                              << std::endl;
                    }
                    main_done=true;
                  }
                }

                // If we have requested a particular number of steps
                // in parallel, then end the inner loop and continue
                // later
                if (steps_in_parallel>0 &&
                    mcmc_iters[it]%steps_in_parallel==0) {
                  inner_done=true;
                }

                // End of while loop for inner_done==false and
                // main_done==false
              }

              // Loop over threads for aff_inv=false
            }

            // End of parallel region for aff_inv=false
          }
          
          // Collect best point over all threads
          for(size_t it=0;it<n_threads;it++) {
            if (w_best_t[it]>w_best) {
              w_best=w_best_t[it];
              best=best_t[it];
            }
          }

          if (main_done==false && max_iters>0) {
            main_done=true;
            for(size_t it=0;it<n_threads;it++) {
              if (mcmc_iters[it]<max_iters) main_done=false;
            }
          }
            
          // Call function outside parallel region 
          outside_parallel();
          
          // End of 'main_done' while loop for aff_inv=false
        }
        
      } else {
    
        // ---------------------------------------------------
        // Start of main loop for aff_inv=true

        bool main_done=false;
        size_t mcmc_iters=0;

        //std::vector<std::string> message(n_threads);
      
        while (!main_done) {

          std::vector<double> smove_z(n_threads);
      
          // ----------------------------------------------------------
          // First parallel region to make the stretch move and 
          // call the object function
          
#ifdef O2SCL_SET_OPENMP
#pragma omp parallel default(shared)
#endif
          {
#ifdef O2SCL_SET_OPENMP
#pragma omp for
#endif
            for(size_t it=0;it<n_threads;it++) {
              
              // Choose walker to move. If the threads are not coupled,
              // then each thread maintains its own ensemble, and we
              // just loop over all of the walkers
              curr_walker[it]=mcmc_iters % n_walk;
            
              // ---------------------------------------------------
              // Select next point
          
              // Total number of walkers
              size_t n_tot;
              if (couple_threads) {
                n_tot=n_walk*n_threads;
              } else {
                n_tot=n_walk;
              }

              // Choose jth walker, between 0 and n_tot-1, making
              // sure ij != curr_walker[it]
              size_t ij;
              do {
                ij=((size_t)(rg[it].random()*((double)n_tot)));
              } while (ij==curr_walker[it] || ij>=n_tot);
            
              // Select z 
              double p=rg[it].random();
              double a=step_fac;
              smove_z[it]=(1.0-2.0*p+2.0*a*p+p*p-2.0*a*p*p+a*a*p*p)/a;

              // The variable jt is the thread index of the walker with
              // index ij. If couple_threads is true and ij is greater
              // than n_walk, then we modify ij and jt to refer to the
              // correct walker in the next thread, wrapping around if
              // necessary
              size_t jt=it;
              if (couple_threads && ij>=n_walk) {
                jt=(jt+ij/n_walk)%n_threads;
                ij=ij%n_walk;
              }

              if (couple_threads) {
                if (jt>=n_threads || ij>=n_walk) {
                  O2SCL_ERR("Variable jt or ij wrong in mcmc_para",
                            o2scl::exc_esanity);
                }
              }

              // At this point, the index of the current walker
              //
              // is n_walk*it+curr_walker[it]
              // 
              // and the index of the jth walker for the stretch-move
              // 
              // is n_walk*jt+ij
            
              // Create new trial point
              for(size_t i=0;i<n_params;i++) {
                if (n_walk*jt+ij>=current.size() ||
                    n_walk*it+curr_walker[it]>=current.size()) {
                  O2SCL_ERR("Walker arithmetic wrong in mcmc_para",
                            o2scl::exc_esanity);
                }
                next[it][i]=current[n_walk*jt+ij][i]+
                  smove_z[it]*(current[n_walk*it+curr_walker[it]][i]-
                               current[n_walk*jt+ij][i]);
              }
              /*
                std::ostringstream os;
                os << n_walk << " "
                << it << " " << curr_walker[it] << " "
                << jt << " " << ij << " "
                << smove_z[it] << std::endl;
                o2scl::vector_out(os,current[n_walk*it+curr_walker[it]],true);
                o2scl::vector_out(os,next[it],true);
                message[it]=os.str();
              */
            
              // ---------------------------------------------------
              // Compute next weight
      
              func_ret[it]=o2scl::success;
              // If the next point out of bounds, ensure that the point is
              // rejected without attempting to evaluate the function
              for(size_t k=0;k<n_params;k++) {
                if (next[it][k]<low[k] || next[it][k]>high[k]) {
                  func_ret[it]=mcmc_skip;
                  if (verbose>=3) {
                    if (next[it][k]<low[k]) {
                      std::cout << "mcmc (" << it << ","
                                << mpi_rank << "): Parameter with index "
                                << k << " and value " << next[it][k]
                                << " smaller than limit " << low[k]
                                << std::endl;
                      scr_out << "mcmc (" << it << ","
                              << mpi_rank << "): Parameter with index "
                              << k << " and value " << next[it][k]
                              << " smaller than limit " << low[k]
                              << std::endl;
                    } else {
                      std::cout << "mcmc (" << it << "," << mpi_rank
                                << "): Parameter with index " << k
                                << " and value " << next[it][k]
                                << " larger than limit " << high[k]
                                << std::endl;
                      scr_out << "mcmc (" << it << "," << mpi_rank
                              << "): Parameter with index " << k
                              << " and value " << next[it][k]
                              << " larger than limit " << high[k]
                              << std::endl;
                    }
                  }
                }
              }

              // Evaluate the function, set the 'done' flag if
              // necessary, and update the return value array
              if (func_ret[it]!=mcmc_skip) {
                if (switch_arr[n_walk*it+curr_walker[it]]==false) {
                  func_ret[it]=func[it](n_params,next[it],w_next[it],
                                        data[it*n_walk+curr_walker[it]+
                                             n_walk*n_threads]);
                } else {
                  func_ret[it]=func[it](n_params,next[it],w_next[it],
                                        data[it*n_walk+curr_walker[it]]);
                }
                if (func_ret[it]==mcmc_done) {
                  mcmc_done_flag[it]=true;
                } else {
                  if (func_ret[it]>=0 && ret_value_counts.size()>it && 
                      func_ret[it]<((int)ret_value_counts[it].size())) {
                    ret_value_counts[it][func_ret[it]]++;
                  }
                }

              }
            }
          }
          // End of first parallel region for aff_inv=true

          // ---------------------------------------------------------
          // Post-function verbose output in case parameter was out of
          // range, function returned "done" or a failure. More
          // verbose output is performed below after the possible call
          // to the measurement function.

          if (verbose>=1) {
            for(size_t it=0;it<n_threads;it++) {
              //std::cout << "thread: " << it << " message: "
              //<< message[it] << std::endl;
              if (pd_mode) {
                scr_out << "it: " << it << " q_prop[it]: "
                        << q_prop[it] << std::endl;
              }
              if (func_ret[it]==mcmc_done) {
                scr_out << "mcmc (" << it << "," << mpi_rank
                        << "): Returned mcmc_done." 
                        << std::endl;
              } else if (func_ret[it]==mcmc_skip && verbose>=3) {
                scr_out << "mcmc (" << it
                        << "): Parameter(s) out of range: " << std::endl;
                scr_out.setf(std::ios::showpos);
                for(size_t k=0;k<n_params;k++) {
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
                  scr_out << "mcmc (" << it << "," << mpi_rank
                          << "): Function returned failure " 
                          << func_ret[it] << " at point ";
                  for(size_t k=0;k<n_params;k++) {
                    scr_out << next[it][k] << " ";
                  }
                  scr_out << std::endl;
                }
              }
            }
          }

          // ----------------------------------------------------------
          // Second parallel region to accept or reject, and call
          // measurement function
      
#ifdef O2SCL_SET_OPENMP
#pragma omp parallel default(shared)
#endif
          {
#ifdef O2SCL_SET_OPENMP
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
            
                double ai_ratio=pow(smove_z[it],((double)n_params)-1.0)*
                  exp(w_next[it]-w_current[sindex]);
                if (r<ai_ratio) {
                  accept=true;
                }

                // End of 'if (func_ret[it]==o2scl::success)'
              }

              if (accept) {
          
                n_accept[it]++;
          
                // Store results from new point
                if (!warm_up) {
                  if (switch_arr[sindex]==false) {
                    meas_ret[it]=meas[it](next[it],w_next[it],
                                          curr_walker[it],func_ret[it],true,
                                          data[sindex+n_threads*n_walk]);
                  } else {
                    meas_ret[it]=meas[it](next[it],w_next[it],
                                          curr_walker[it],func_ret[it],true,
                                          data[sindex]);
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
                    meas_ret[it]=meas[it](next[it],w_next[it],
                                          curr_walker[it],func_ret[it],false,
                                          data[sindex+n_threads*n_walk]);
                  } else {
                    meas_ret[it]=meas[it](next[it],w_next[it],
                                          curr_walker[it],func_ret[it],false,
                                          data[sindex]);
                  }
                }

              }

            }
          }
          // End of second parallel region for aff_inv=true

          // -----------------------------------------------------------
          // Post-measurement verbose output of iteration count, weight,
          // and walker index for each thread
      
          if (verbose>=2) {
            for(size_t it=0;it<n_threads;it++) {
              size_t sindex=n_walk*it+curr_walker[it];
              scr_out.precision(4);
              scr_out << "mcmc (" << it << "," << mpi_rank << "): iter: ";
              scr_out.width((int)(1.0+log10((double)(n_params-1))));
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

          // Update iteration count and reset counters for
          // warm up iterations if necessary
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
                scr_out << "mcmc: Finished warmup." << std::endl;
              }
          
            }
          }

          if (verbose>=3) {
            std::cout << "Press a key and type enter to continue. ";
            char ch;
            std::cin >> ch;
          }

          // Stop if iterations greater than max
          if (main_done==false && warm_up==false && max_iters>0 &&
              mcmc_iters==max_iters) {
            if (verbose>=1) {
              scr_out << "mcmc: Stopping because number of iterations ("
                      << mcmc_iters << ") equal to max_iters (" << max_iters
                      << ")." << std::endl;
            }
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
              if (verbose>=1) {
                scr_out << "mcmc: Stopping because elapsed (" << elapsed
                        << ") > max_time (" << max_time << ")."
                        << std::endl;
              }
              main_done=true;
            }
          }

          outside_parallel();

          // --------------------------------------------------------------
          // End of main loop for aff_inv=true
        }

        // End of conditional for aff_inv=true
      }
    
      // --------------------------------------------------------------
    
      mcmc_cleanup();

      return 0;
    }
    
    /** \brief Perform a MCMC simulation with a thread-safe function
        or with only one OpenMP thread
    */
    virtual int mcmc(size_t n_params, vec_t &low, vec_t &high,
                     func_t &func, measure_t &meas, data_t data) {
    
#ifdef O2SCL_SET_OPENMP
      omp_set_num_threads(n_threads);
#else
      n_threads=1;
#endif
      std::vector<func_t> vf(n_threads);
      std::vector<measure_t> vm(n_threads);
      std::vector<data_t> vd(2*n_threads*n_walk);
      for(size_t i=0;i<n_threads;i++) {
        vf[i]=func;
        vm[i]=meas;
      }
      for(size_t i=0;i<2*n_threads*n_walk;i++) {
        vd[i]=data;
      }
      return mcmc(n_params,low,high,vf,vm,vd);
    }
    //@}

    /// \name Proposal distribution
    //@{
    /** \brief Set the proposal distribution

        \note This function automatically sets \ref aff_inv to false
        and \ref n_walk to 1. The vector of proposal distributions needs
        to have an element for each thread.
    */
    template<class prob_vec_t> void set_proposal(prob_vec_t &pv) {
      prop_dist.resize(pv.size());
      for(size_t i=0;i<pv.size();i++) {
        prop_dist[i]=&pv[i];
      }
      pd_mode=true;
      aff_inv=false;
      n_walk=1;
      return;
    }

    /** \brief Set pointers to proposal distributions

        \note This function automatically sets \ref aff_inv to false
        and \ref n_walk to 1. The vector of proposal distributions
        needs to have an element for each thread.
    */
    template<class prob_ptr_vec_t> void set_proposal_ptrs(prob_ptr_vec_t &pv) {
      prop_dist.resize(pv.size());
      for(size_t i=0;i<pv.size();i++) {
        prop_dist[i]=pv[i];
      }
      pd_mode=true;
      aff_inv=false;
      n_walk=1;
      return;
    }

    /** \brief Go back to random-walk Metropolis with a uniform distribution
     */
    virtual void unset_proposal() {
      if (pd_mode) {
        prop_dist.clear();
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

      There is some output which occurs in addition to the output
      from \ref o2scl::mcmc_para_base depending on the value 
      of \ref o2scl::mcmc_para_base::verbose . If there is 
      a misalignment between the number of columns in the 
      table and the number of data points in any line, then
      some debugging information is sent to <tt>cout</tt>.
      When verbose is 2 or larger, ... (FIXME)

      \note This class is experimental.

      \future Verbose output may need improvement
      \future Use reorder_table() and possibly reblock()
      to create a full post-processing function.
  */
  template<class func_t, class fill_t, class data_t, class vec_t=ubvector,
           class stepper_t=mcmc_stepper_base<func_t,data_t,vec_t>>
  class mcmc_para_table :
    public mcmc_para_base<func_t,
                          std::function<int(const vec_t &,
                                            double,size_t,
                                            int,bool,data_t &)>,
                          data_t,vec_t,stepper_t> {
    
  protected:
  
    /// Measurement functor type for the parent
    typedef std::function<int(const vec_t &,double,size_t,int,bool,data_t &)>
    internal_measure_t;
  
    /// Type of parent class
    typedef mcmc_para_base<func_t,internal_measure_t,data_t,vec_t,
                           stepper_t> parent_t;

    /// Column names
    std::vector<std::string> col_names;
    
    /// Column units
    std::vector<std::string> col_units;

    /// Number of parameters
    size_t n_params;
  
    /// Main data table for Markov chain
    std::shared_ptr<o2scl::table_units<> > table;

    /** \brief If true, the HDF5 I/O initial info has been written
        to the file (set by \ref mcmc() )
    */
    bool first_write;

    /** \brief MCMC initialization function

        This function sets the column names and units.
    */
    virtual int mcmc_init() {
    
      if (!prev_read) {
      
        // -----------------------------------------------------------
        // Initialize table, walker_accept_rows, and walker_reject_rows

        if (table==0) {
          table=std::shared_ptr<o2scl::table_units<> >
            (new o2scl::table_units<>);
        } else {
          table->clear();
        }
        if (table_prealloc>0) {
          table->set_maxlines(table_prealloc);
        }
        table->new_column("rank");
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
      
        walker_accept_rows.resize(this->n_walk*this->n_threads);
        for(size_t i=0;i<this->n_walk*this->n_threads;i++) {
          walker_accept_rows[i]=-1;
        }
        walker_reject_rows.resize(this->n_walk*this->n_threads);
        for(size_t i=0;i<this->n_walk*this->n_threads;i++) {
          walker_reject_rows[i]=-1;
        }
      
        if (false && this->verbose>=3) {
          // AWS, 8/19/23: I took this out because it sends too much
          // output to cout
          std::cout << "mcmc: Table column names and units: " << std::endl;
          for(size_t i=0;i<table->get_ncolumns();i++) {
            std::cout << table->get_column_name(i) << " "
                      << table->get_unit(table->get_column_name(i))
                      << std::endl;
          }
        }
      
      } else {

        if (table==0) {
          O2SCL_ERR2("Flag 'prev_read' is true but table pointer is 0 ",
                     "in mcmc_para_table::mcmc_init().",o2scl::exc_esanity);
        }
      
        // -----------------------------------------------------------
        // Previous results are already present

        if (table->get_ncolumns()!=5+col_names.size()) {
          std::string str=((std::string)"Table does not have correct ")+
            "number of columns in mcmc_para_table::mcmc_init()."+
            o2scl::szttos(table->get_ncolumns())+" columns and "+
            o2scl::szttos(col_names.size())+" entries in col_names.";
          O2SCL_ERR(str.c_str(),o2scl::exc_einval);
        }
        if (!table->is_column("rank") ||
            !table->is_column("thread") ||
            !table->is_column("walker") ||
            !table->is_column("mult") ||
            !table->is_column("log_wgt")) {
          O2SCL_ERR2("Table does not have the correct internal columns ",
                     "in mcmc_para_table::mcmc_init().",o2scl::exc_einval);
        }
        if (walker_accept_rows.size()!=this->n_walk*this->n_threads) {
          O2SCL_ERR2("Array walker_accept_rows does not have correct size ",
                     "in mcmc_para_table::mcmc_init().",o2scl::exc_einval);
        }
        if (walker_reject_rows.size()!=this->n_walk*this->n_threads) {
          O2SCL_ERR2("Array walker_reject_rows does not have correct size ",
                     "in mcmc_para_table::mcmc_init().",o2scl::exc_einval);
        }

        // Set prev_read to false so that next call to mcmc()
        // doesn't use the previously read results.
        prev_read=false;
      
      }
    
      last_write_iters=0;
#ifdef O2SCL_MPI
      last_write_time=MPI_Wtime();
#else
      last_write_time=time(0);
#endif

      return parent_t::mcmc_init();
    }
  
    /** \brief Fill \c line with data for insertion into the table
     */
    virtual int fill_line(const vec_t &pars, double log_weight, 
                          std::vector<double> &line, data_t &dat,
                          size_t i_walker, fill_t &fill) {

#ifdef O2SCL_SET_OPENMP
      size_t i_thread=omp_get_thread_num();
#else
      size_t i_thread=0;
#endif

      // Rank
      line.push_back(this->mpi_rank);
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
  
    /** \brief For each walker and thread, record the last row in the
        table which corresponds to an accept
    */
    std::vector<int> walker_accept_rows;

    /** \brief For each walker and thread, record the last row in the
        table which corresponds to an reject
    */
    std::vector<int> walker_reject_rows;

    /** \brief Initial write to HDF5 file 
     */
    virtual void file_header(o2scl_hdf::hdf_file &hf) {
      return;
    }
  
    /** \brief A copy of the lower limits for HDF5 output
     */
    vec_t low_copy;
  
    /** \brief A copy of the upper limits for HDF5 output
     */
    vec_t high_copy;
  
    /** \brief Total number of MCMC acceptances over all threads at last
        file write() (default 0)
    */
    size_t last_write_iters;
  
    /** \brief Time at last
        file write() (default 0.0)
    */
    double last_write_time;

    /** \brief If true, previous results have been read
      
        This is set to <tt>true</tt> by \ref read_prev_results()
        and set back to <tt>false</tt> after mcmc_init() is called.
    */
    bool prev_read;
  
  public:

    /// \name Settings
    //@{
    /** \brief If true, ensure sure walkers and OpenMP threads are
        written to the table with equal spacing between rows (default
        true)
    */
    bool table_sequence;
  
    /** \brief Iterations between file updates (default 0 for no file updates)
     */
    size_t file_update_iters;
  
    /** \brief Time between file updates (default 0.0 for no file updates)
     */
    double file_update_time;

    /** \brief Number of rows to allocate for the table before the MCMC run
     */
    size_t table_prealloc;
  
    /** \brief The number of tables to combine before I/O (default 1)
     */
    int table_io_chunk;

    /** \brief If true, store MCMC rejections in the table (default false)
     */
    bool store_rejects;

    /** \brief If true, check rows (default true)
     */
    bool check_rows;
    //@}
  
    /** \brief Write MCMC tables to files
     */
    virtual void write_files(bool sync_write=false) {

      if (this->verbose>=2) {
        this->scr_out << "mcmc: Start write_files(). mpi_rank: "
                      << this->mpi_rank << " mpi_size: "
                      << this->mpi_size <<  " table_io_chunk: "
                      << table_io_chunk << std::endl;
      }
    
      std::vector<o2scl::table_units<> > tab_arr;
      bool rank_sent=false;
    
#ifdef O2SCL_MPI
      if (table_io_chunk>1) {
        if (this->mpi_rank%table_io_chunk==0) {
          // Parent ranks
          for(int i=0;i<table_io_chunk-1;i++) {
            int child=this->mpi_rank+i+1;
            if (child<this->mpi_size) {
              table_units<> t;
              tab_arr.push_back(t);
              o2scl_table_mpi_recv(child,tab_arr[tab_arr.size()-1]);
            }
          }
        } else {
          // Child ranks
          size_t parent=this->mpi_rank-(this->mpi_rank%table_io_chunk);
          o2scl_table_mpi_send(*table,parent);
          rank_sent=true;
        }
      }
#endif

#ifdef O2SCL_MPI
      // Ensure that multiple threads aren't writing to the
      // filesystem at the same time
      int tag=0, buffer=0;
      if (sync_write && this->mpi_size>1 &&
          this->mpi_rank>=table_io_chunk) {
        MPI_Recv(&buffer,1,MPI_INT,this->mpi_rank-table_io_chunk,
                 tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
#endif
    
      o2scl_hdf::hdf_file hf;
      std::string fname=this->prefix+"_"+o2scl::itos(this->mpi_rank)+"_out";
      hf.open_or_create(fname);

      if (first_write==false) {
        hf.setd("ai_initial_step",this->ai_initial_step);
        hf.seti("aff_inv",this->aff_inv);
        hf.seti("always_accept",this->always_accept);
        hf.setd_vec_copy("high",this->high_copy);
        hf.setd_vec_copy("low",this->low_copy);
        hf.set_szt("max_bad_steps",this->max_bad_steps);
        hf.set_szt("max_iters",this->max_iters);
        hf.set_szt("max_time",this->max_time);
        hf.set_szt("file_update_iters",this->file_update_iters);
        hf.setd("file_update_time",this->file_update_time);
        hf.seti("mpi_rank",this->mpi_rank);
        hf.seti("mpi_size",this->mpi_size);
        hf.set_szt("n_params",this->n_params);
        hf.set_szt("n_threads",this->n_threads);
        hf.set_szt("n_walk",this->n_walk);
        hf.set_szt("n_warm_up",this->n_warm_up);
        hf.seti("pd_mode",this->pd_mode);
        hf.sets("prefix",this->prefix);
        hf.setd("step_fac",this->step_fac);
        hf.seti("store_rejects",this->store_rejects);
        hf.seti("table_sequence",this->table_sequence);
        hf.seti("user_seed",this->user_seed);
        hf.seti("verbose",this->verbose);
        file_header(hf);
        first_write=true;
      }
    
      hf.set_szt_vec("n_accept",this->n_accept);
      hf.set_szt_vec("n_reject",this->n_reject);
      if (this->ret_value_counts.size()>0) {
        hf.set_szt_arr2d_copy("ret_value_counts",this->ret_value_counts.size(),
                              this->ret_value_counts[0].size(),
                              this->ret_value_counts);
      }
      hf.setd_arr2d_copy("initial_points",this->initial_points.size(),
                         this->initial_points[0].size(),
                         this->initial_points);

      hf.seti("n_tables",tab_arr.size()+1);
      if (rank_sent==false) {
        hdf_output(hf,*table,"markov_chain_0");
      }
      for(size_t i=0;i<tab_arr.size();i++) {
        std::string name=((std::string)"markov_chain_")+szttos(i+1);
        hdf_output(hf,tab_arr[i],name);
      }
    
      hf.close();
    
#ifdef O2SCL_MPI
      if (sync_write && this->mpi_size>1 &&
          this->mpi_rank<this->mpi_size-1) {
        MPI_Send(&buffer,1,MPI_INT,this->mpi_rank+table_io_chunk,
                 tag,MPI_COMM_WORLD);
      }
#endif
    
      if (this->verbose>=2) {
        this->scr_out << "mcmc: Done write_files()." << std::endl;
      }

      return;
    }
  
    mcmc_para_table() {
      table_io_chunk=1;
      file_update_iters=0;
      file_update_time=0.0;
      last_write_iters=0;
      store_rejects=false;
      table_sequence=true;
      prev_read=false;
      table_prealloc=0;
      check_rows=true;
    }
  
    /// \name Basic usage
    //@{
    /** \brief Set the table names and units
     */
    virtual void set_names_units(std::vector<std::string> names,
                                 std::vector<std::string> units) {
      if (names.size()!=units.size()) {
        O2SCL_ERR2("Size of names and units arrays don't match in ",
                   "mcmc_para_table::set_names_units().",
                   o2scl::exc_einval);
      }
      col_names=names;
      col_units=units;
      return;
    }

    /** \brief Read initial points from the last points recorded in file
        named \c fname

        The values of \ref o2scl::mcmc_para_base::n_walk and \ref
        o2scl::mcmc_para_base::n_threads, must be set to their correct
        values before calling this function. This function requires that
        a table is present in \c fname which stores parameters in a
        block of columns and has columns named \c mult, \c thread, \c
        walker, and \c log_wgt. This function does not double check
        that the columns in the file associated with the parameters 
        have the correct names.
    */
    virtual void initial_points_file_last(std::string fname,
                                          size_t n_param_loc,
                                          size_t offset=5) {

      o2scl::table_units<> tip;

#ifdef O2SCL_MPI
      // Ensure that multiple threads aren't reading from the
      // filesystem at the same time
      int tag=0, buffer=0;
      if (this->mpi_size>1 && this->mpi_rank>0) {
        MPI_Recv(&buffer,1,MPI_INT,this->mpi_rank-1,
                 tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
#endif
    
      o2scl_hdf::hdf_file hf;
      hf.open(fname);
      std::string tname;
      hdf_input(hf,tip,tname);
      hf.close();
    
#ifdef O2SCL_MPI
      if (this->mpi_size>1 && this->mpi_rank<this->mpi_size-1) {
        MPI_Send(&buffer,1,MPI_INT,this->mpi_rank+1,
                 tag,MPI_COMM_WORLD);
      }
#endif

      // Determine number of points
      size_t n_points=this->n_walk*this->n_threads;

      if (this->verbose>0) {
        std::cout << "Initial points: Finding last " << n_points
                  << " points from file named "
                  << fname << " ." << std::endl;
      }
    
      this->initial_points.resize(n_points);
        
      for(size_t it=0;it<this->n_threads;it++) {
        for(size_t iw=0;iw<this->n_walk;iw++) {
        
          // The combined walker/thread index 
          size_t windex=it*this->n_walk+iw;

          bool found=false;
          for(int row=tip.get_nlines()-1;row>=0 && found==false;row--) {
            if (tip.get("walker",row)==iw &&
                tip.get("thread",row)==it &&
                tip.get("mult",row)>0.5) {

              found=true;
            
              std::cout << "Function initial_point_file_last():\n\tit: "
                        << it << " rank: " << this->mpi_rank
                        << " iw: " << iw << " row: "
                        << row << " log_wgt: " << tip.get("log_wgt",row)
                        << std::endl;
            
              // Copy the entries from this row into the initial_points object
              this->initial_points[windex].resize(n_param_loc);
              for(size_t ip=0;ip<n_param_loc;ip++) {
                this->initial_points[windex][ip]=tip.get(ip+offset,row);
              }
            }
          }

          // If we can't find a row with the proper thread and walker
          // index, then just use one of the points from the end of
          // the file
          if (found==false && tip.get_nlines()>this->n_walk*
              this->n_threads) {
            std::cout << "Could not find guess for rank " << this->mpi_rank
                      << " thread " << it
                      << " and walker " << iw
                      << ". Trying to find a point at the end of "
                      << "the file to use." << std::endl;
            int row=tip.get_nlines()-this->n_walk*this->n_threads+windex;
            this->initial_points[windex].resize(n_param_loc);
            for(size_t ip=0;ip<n_param_loc;ip++) {
              this->initial_points[windex][ip]=tip.get(ip+offset,row);
            }
            found=true;
          }
        
          if (found==false) {
            std::cout << "No initial guess found for rank "
                      << this->mpi_rank
                      << " thread " << it << " and walker " << iw
                      << std::endl;
            O2SCL_ERR("Function initial_points_file_last() failed.",
                      o2scl::exc_einval);
          }
        }
      }
    
      return;
    }
  
    /** \brief Read initial points from file named \c fname,
        distributing across the chain if necessary

        The values of \ref o2scl::mcmc_para_base::n_walk and \ref
        o2scl::mcmc_para_base::n_threads, must be set to their
        correct values before calling this function. This function
        requires that a table is present in \c fname which stores
        parameters in a block of columns. This function does not
        double check that the columns in the file associated with the
        parameters have the correct names.
    */
    virtual void initial_points_file_dist(std::string fname,
                                          size_t n_param_loc,
                                          size_t offset=5) {

      o2scl::table_units<> tip;

#ifdef O2SCL_MPI
      // Ensure that multiple threads aren't reading from the
      // filesystem at the same time
      int tag=0, buffer=0;
      if (this->mpi_size>1 && this->mpi_rank>0) {
        MPI_Recv(&buffer,1,MPI_INT,this->mpi_rank-1,
                 tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
#endif
    
      o2scl_hdf::hdf_file hf;
      hf.open(fname);
      std::string tname;
      hdf_input(hf,*table,tname);
      hf.close();
    
#ifdef O2SCL_MPI
      if (this->mpi_size>1 && this->mpi_rank<this->mpi_size-1) {
        MPI_Send(&buffer,1,MPI_INT,this->mpi_rank+1,
                 tag,MPI_COMM_WORLD);
      }
#endif

      // Determine number of points
      size_t n_points=this->n_walk*this->n_threads;

      if (this->verbose>0) {
        std::cout << "Initial points: Finding last " << n_points
                  << " points from file named "
                  << fname << " ." << std::endl;
      }
    
      this->initial_points.resize(n_points);
        
      size_t nlines=tip.get_nlines();
      size_t decrement=nlines/n_points;
      if (decrement<1) decrement=1;

      int row=nlines-1;
      for(size_t k=0;k<n_points;k++) {
        row-=decrement;
        if (row<0) row=0;
      
        std::cout << "Function initial_point_file_dist():\n\trow: "
                  << row << " log_wgt: " << tip.get("log_wgt",row)
                  << std::endl;
    
        // Copy the entries from this row into the initial_points object
        this->initial_points[k].resize(n_param_loc);
        for(size_t ip=0;ip<n_param_loc;ip++) {
          this->initial_points[k][ip]=tip.get(ip+offset,row);
        }
      
      }
    
      return;
    }
  
    /** \brief Read initial points from the best points recorded in file
        named \c fname

        Before calling this function, the values the values of \ref
        o2scl::mcmc_para_base::n_walk and \ref
        o2scl::mcmc_para_base::n_threads should have been set by the
        user to the correct values for the future call to \ref
        mcmc_para_base::mcmc() .

        In order for this function to succeed, a table must be present
        (the function just reads the first \ref o2scl::table_units or
        \ref o2scl::table_units object it can find) in the HDF5 file
        named \c fname. This table should store parameters in a block
        of columns beginning with column \c offset. It should also
        contain a separate column named \c log_wgt for the log
        likelihood. The table must have at least as unique rows as
        input points required by the \ref mcmc_para_base::mcmc()
        function, i.e. the product of \ref
        o2scl::mcmc_para_base::n_walk and
        o2scl::mcmc_para_base::n_threads. Rows are presumed to be
        identical if all of their values differ by a value less than
        \c thresh, which defaults to \f$ 10^{-6} \f$.

        This function does not double check that the columns in the
        file associated with the parameters have the correct names.
        The values in the "walker", "thread", and "rank" columns in
        the table, if present, are ignored by this function.
    */
    virtual void initial_points_file_best(std::string fname,
                                          size_t n_param_loc,
                                          double thresh=1.0e-6,
                                          size_t offset=5) {

      o2scl::table_units<> tip;
    
#ifdef O2SCL_MPI
      // Ensure that multiple threads aren't reading from the
      // filesystem at the same time
      int tag=0, buffer=0;
      if (this->mpi_size>1 && this->mpi_rank>0) {
        MPI_Recv(&buffer,1,MPI_INT,this->mpi_rank-1,
                 tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      }
#endif
    
      o2scl_hdf::hdf_file hf;
      hf.open(fname);
      std::string tname;
      hdf_input(hf,tip,tname);
      hf.close();
    
#ifdef O2SCL_MPI
      if (this->mpi_size>1 && this->mpi_rank<this->mpi_size-1) {
        MPI_Send(&buffer,1,MPI_INT,this->mpi_rank+1,
                 tag,MPI_COMM_WORLD);
      }
#endif

      // Determine number of points
      size_t n_points=this->n_walk*this->n_threads;

      if (this->verbose>0) {
        std::cout << "Initial points: Finding best " << n_points
                  << " unique points from file named "
                  << fname << " using threshold " << thresh
                  << " and offset " << offset << std::endl;
      }

      typedef std::map<double,int,std::greater<double> > map_t;
      map_t m;

      // Sort by inserting into a map
      for(size_t k=0;k<tip.get_nlines();k++) {
        m.insert(std::make_pair(tip.get("log_wgt",k),k));
      }

      // Remove near duplicates. The map insert function will 
      // just overwrite duplicate key entries, but we also
      // want to avoid near duplicates, so we have to search
      // manually for those. 
      bool found;
      do {
        found=false;
        for(map_t::iterator mit=m.begin();mit!=m.end();mit++) {
          map_t::iterator mit2=mit;
          mit2++;
          if (mit2!=m.end()) {
            if (fabs(mit->first-mit2->first)<thresh) {
              if (this->verbose>0) {
                std::cout << "mcmc_para::initial_points_file_best():\n  "
                          << "Removing duplicate weights: "
                          << mit->first << " " << mit2->first << std::endl;
                
              }
              m.erase(mit2);
              mit=m.begin();
              found=true;
            }
          }
        }
      } while (found==true);

      // Check to see if we have enough
      if (m.size()<n_points) {
	std::cerr << "m.size() is " << m.size() << " and n_points is "
		  << n_points << std::endl;
        O2SCL_ERR2("Could not find enough points in file in ",
                   "mcmc_para::initial_points_file_best().",
                   o2scl::exc_efailed);
      }

      // Copy the entries from this row into the initial_points object
      this->initial_points.resize(n_points);
      map_t::iterator mit=m.begin();
      for(size_t k=0;k<n_points;k++) {
        int row=mit->second;
        if (this->verbose>0) {
          std::cout << "Initial point " << k << " at row "
                    << row << " has log_weight= "
                    << tip.get("log_wgt",row) << std::endl;
        }
        this->initial_points[k].resize(n_param_loc);
        for(size_t ip=0;ip<n_param_loc;ip++) {
          this->initial_points[k][ip]=tip.get(ip+offset,row);
        }
        mit++;
      }

      return;
    }
    
    /** \brief Perform an MCMC simulation
        
        Perform an MCMC simulation over \c n_params parameters starting
        at initial point \c init, limiting the parameters to be between
        \c low and \c high, using \c func as the objective function and
        calling the measurement function \c meas at each MC point.
    */
    virtual int mcmc_fill(size_t n_params_local, 
                          vec_t &low, vec_t &high,
                          std::vector<func_t> &func,
                          std::vector<fill_t> &fill,
                          std::vector<data_t> &data) {
    
      n_params=n_params_local;
      low_copy=low;
      high_copy=high;
    
      first_write=false;
    
      // Set number of threads (this is done in the child as well, but
      // we need this number to set up the vector of measure functions
      // below).
#ifdef O2SCL_SET_OPENMP
      omp_set_num_threads(this->n_threads);
#else
      this->n_threads=1;
#endif

      // Setup the vector of measure functions
      std::vector<internal_measure_t> meas(this->n_threads);
      for(size_t it=0;it<this->n_threads;it++) {
        meas[it]=std::bind
          (std::mem_fn<int(const vec_t &,double,size_t,int,bool,
                           data_t &, size_t, fill_t &)>
           (&mcmc_para_table::add_line),this,std::placeholders::_1,
           std::placeholders::_2,std::placeholders::_3,
           std::placeholders::_4,std::placeholders::_5,
           std::placeholders::_6,it,std::ref(fill[it]));
      }
    
      return parent_t::mcmc(n_params,low,high,func,meas,data);
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

        \future This algorithm could be improved by started from the end
        of the table and going backwards instead of starting from the
        front of the table and going forwards.
    */
    void get_chain_sizes(std::vector<size_t> &chain_sizes) {

      size_t ntot=this->n_threads*this->n_walk;
      chain_sizes.resize(ntot);
    
      for(size_t it=0;it<this->n_threads;it++) {
        for(size_t iw=0;iw<this->n_walk;iw++) {
          size_t ix=it*this->n_walk+iw;
          size_t istart=ix;
          chain_sizes[ix]=0;
          for(size_t j=istart;j<table->get_nlines();j+=ntot) {
            if (table->get("mult",j)>0.5) chain_sizes[ix]++;
          }
        }
      }
    
      return;
    }

    /** \brief Read previous results (number of threads and 
        walkers must be set first)

        \note By default, this tries to obtain the initial points
        for the next call to \ref mcmc() by the previously 
        accepted point in the table. 

        \note This function requires a table correctly stored with
        the right column order
    */
    virtual void read_prev_results(o2scl_hdf::hdf_file &hf,
                                   size_t n_param_loc,
                                   std::string name="") {

      // Create the table object
      table=std::shared_ptr<o2scl::table_units<> >(new o2scl::table_units<>);

      // Read the table data from the HDF5 file
      hdf_input(hf,*table,name);

      if (!table->is_column("rank") ||
          !table->is_column("thread") ||
          !table->is_column("walker") ||
          !table->is_column("mult") ||
          !table->is_column("log_wgt")) {
        O2SCL_ERR2("Table does not have the correct internal columns ",
                   "in mcmc_para_table::read_prev_results().",
                   o2scl::exc_einval);
      }
    
      // -----------------------------------------------------------
      // Set the values of walker_accept_rows and walker_reject_rows
      // by finding the last accepted row and the last rejected row
      // for each walker and each thread.
    
      // The total number of walkers * threads
      size_t ntot=this->n_threads*this->n_walk;
    
      walker_accept_rows.resize(ntot);
      walker_reject_rows.resize(ntot);

      for(size_t j=0;j<ntot;j++) {
        walker_accept_rows[j]=-1;
        walker_reject_rows[j]=-1;
      }
    
      for(size_t j=0;j<table->get_nlines();j++) {
      
        size_t i_thread=((size_t)(table->get("thread",j)+1.0e-12));
        size_t i_walker=((size_t)(table->get("walker",j)+1.0e-12));

        // The combined walker/thread index 
        size_t windex=i_thread*this->n_walk+i_walker;

        if (table->get("mult",j)>0.5) {
          walker_accept_rows[windex]=j;
        } else if (table->get("mult",j)<-0.5) {
          walker_reject_rows[windex]=j;
        }

      }

      // Only set initial points if we found an acceptance for
      // all walkers and threads
      bool found=true;
      for(size_t j=0;j<ntot;j++) {
        if (walker_accept_rows[j]<0) found=false;
      }

      if (found) {
        // Set up initial points
        this->initial_points.clear();
        this->initial_points.resize(ntot);
        for(size_t j=0;j<ntot;j++) {
          this->initial_points[j].resize(n_param_loc);
          for(size_t k=0;k<n_param_loc;k++) {
            this->initial_points[j][k]=table->get(k+5,walker_accept_rows[j]);
          }
        }
      } else {
        std::cout << "Previous table was read, but initial points not set."
                  << std::endl;
      }

      if (this->verbose>0) {
        std::cout << "mcmc_para_table::read_prev_results():" << std::endl;
        std::cout << "  index walker_accept_rows walker_reject_rows"
                  << std::endl;
        for(size_t j=0;j<ntot;j++) {
          std::cout << "  ";
          std::cout.width(3);
          std::cout << j << " ";
          std::cout.width(5);
          std::cout << walker_accept_rows[j] << " ";
          std::cout.width(5);
          std::cout << walker_reject_rows[j] << std::endl;
        }
      }
    
      prev_read=true;
      this->meas_for_initial=false;

      return;
    }

    /** \brief Additional code to execute inside the
        OpenMP critical section
    */
    virtual void critical_extra(size_t i_thread) {

      if (i_thread==0) {
      
        // If necessary, output to files (only thread 0)
        bool updated=false;
        if (file_update_iters>0) {
          size_t total_iters=0;
          for(size_t it=0;it<this->n_threads;it++) {
            total_iters+=this->n_accept[it]+this->n_reject[it];
          }
          if (total_iters>=last_write_iters+file_update_iters) {
            if (this->verbose>=1) {
              this->scr_out << "mcmc: Writing to file. total_iters: "
                            << total_iters << " file_update_iters: "
                            << file_update_iters << " last_write_iters: "
                            << last_write_iters << std::endl;
            }
            write_files(false);
            last_write_iters=total_iters;
            updated=true;
          }
        }
        if (updated==false && file_update_time>0.0) {
#ifdef O2SCL_MPI
          double elapsed=MPI_Wtime()-last_write_time;
#else
          double elapsed=time(0)-last_write_time;
#endif
          if (elapsed>file_update_time) {
            if (this->verbose>=1) {
              this->scr_out << "mcmc: Writing to file. elapsed: "
                            << elapsed << " file_update_time: "
                            << file_update_time << " last_write_time: "
                            << last_write_time << std::endl;
            }
            write_files(false);
#ifdef O2SCL_MPI
            last_write_time=MPI_Wtime();
#else
            last_write_time=time(0);
#endif
          }
        }
      }

      return;
    }
  
    /** \brief A measurement function which adds the point to the
        table
    */
    virtual int add_line(const vec_t &pars, double log_weight,
                         size_t walker_ix, int func_ret,
                         bool mcmc_accept, data_t &dat,
                         size_t i_thread, fill_t &fill) {

      // The combined walker/thread index 
      size_t windex=i_thread*this->n_walk+walker_ix;

      // The total number of walkers * threads
      size_t ntot=this->n_threads*this->n_walk;
    
      int ret_value=o2scl::success;

      // Determine the next row
      int next_row;
      if ((mcmc_accept || store_rejects) && walker_accept_rows[windex]<0) {
        next_row=windex;
      } else {
        if (table_sequence) {
          if (walker_accept_rows[windex]>walker_reject_rows[windex]) {
            next_row=walker_accept_rows[windex]+ntot;
          } else {
            next_row=walker_reject_rows[windex]+ntot;
          }
        } else {
          if (walker_accept_rows[windex]>walker_reject_rows[windex]) {
            next_row=walker_accept_rows[windex]+1;
          } else {
            next_row=walker_reject_rows[windex]+1;
          }
        }
      }
    
#ifdef O2SCL_SET_OPENMP
#pragma omp critical (o2scl_mcmc_para_table_add_line)
#endif
      {

        while (next_row<((int)table->get_nlines()) &&
               fabs(table->get("mult",next_row))>0.1) {
          next_row++;
        }
      
        // If there's not enough space in the table for this iteration,
        // then create it. There is not enough space if any of the
        // walker_accept_rows array entries is -1, if we have an
        // acceptance but there isn't room to store it, or if
        // we have a rejection and there isn't room to store it.
        if (next_row>=((int)table->get_nlines())) {
          size_t istart=table->get_nlines();
          // Create enough space
          table->set_nlines(table->get_nlines()+ntot);
          // Now additionally initialize the first four colums
          for(size_t j=0;j<this->n_threads;j++) {
            for(size_t i=0;i<this->n_walk;i++) {
              table->set("rank",istart+j*this->n_walk+i,this->mpi_rank);
              table->set("thread",istart+j*this->n_walk+i,j);
              table->set("walker",istart+j*this->n_walk+i,i);
              table->set("mult",istart+j*this->n_walk+i,0.0);
              table->set("log_wgt",istart+j*this->n_walk+i,0.0);
            }
          }
        }
      
        // If needed, add the line to the next row
        if (func_ret==0 && (mcmc_accept || store_rejects)) {
        
          if (next_row>=((int)(table->get_nlines()))) {
            O2SCL_ERR("Not enough space in table.",o2scl::exc_esanity);
          }
        
          std::vector<double> line;
          int fret=fill_line(pars,log_weight,line,dat,walker_ix,fill);
        
          // For rejections, set the multiplier to -1.0 (it was set to
          // 1.0 in the fill_line() call above)
          if (store_rejects && mcmc_accept==false) {
            line[3]=-1.0;
          }
        
          if (fret!=o2scl::success) {
          
            // If we're done, we stop before adding the last point to the
            // table. This is important because otherwise the last line in
            // the table will always only have unit multiplicity, which
            // may or may not be correct.
            ret_value=this->mcmc_done;
          
          } else if (fret!=this->mcmc_skip) {

            // First, double check that the table has the right
            // number of columns
            if (line.size()!=table->get_ncolumns()) {
              std::cout << "line: " << line.size() << " columns: "
                        << table->get_ncolumns() << std::endl;
              for(size_t k=0;k<table->get_ncolumns() || k<line.size();k++) {
                std::cout << k << ". ";
                if (k<table->get_ncolumns()) {
                  std::cout << table->get_column_name(k) << " ";
                  std::cout << table->get_unit(table->get_column_name(k))
                            << " ";
                }
                if (k<line.size()) std::cout << line[k] << " ";
                std::cout << std::endl;
              }
              O2SCL_ERR2("Table misalignment in ",
                         "mcmc_para_table::add_line().",
                         exc_einval);
            }
          
            if (this->verbose>=4) {
              for(size_t i=0;i<line.size();i++) {
                std::cout << table->get_column_name(i) << " "
                          << line[i] << std::endl;
              }
            }
          
            // Verbose output
            if (this->verbose>=2) {
              this->scr_out << "mcmc: Thread " << i_thread
                            << " setting data at row " << next_row
                            << std::endl;
              this->scr_out << "  func_ret: " << func_ret << " mcmc_accept: "
                            << mcmc_accept << " walker_ix: "
                            << walker_ix << " store_rejects: "
                            << store_rejects << std::endl;
            }
            if (this->verbose>=3) {
              for(size_t k=0;k<line.size();k++) {
                this->scr_out << k << ". ";
                this->scr_out << table->get_column_name(k) << " ";
                this->scr_out << table->get_unit(table->get_column_name(k));
                this->scr_out << " " << line[k] << std::endl;
              }
            }

            if (check_rows && fabs(line[3])>0.5 &&
                fabs(table->get(3,((size_t)next_row)))>0.5) {
              std::cout << "mult for line is " << line[3]
                        << " and mult at next_row (" << next_row
                        << ") is "
                        << table->get(3,((size_t)next_row)) << std::endl;
              std::cout << "  Walker and thread at next row is "
                        << table->get(2,((size_t)next_row)) << " and "
                        << table->get(1,((size_t)next_row)) << std::endl;
              std::cout << "  " << table->get_nlines() << " "
                        << table->get_maxlines() << std::endl;
              O2SCL_ERR("Row arithmetic problem in mcmc_para_table.",
                        o2scl::exc_esanity);
            }

            // AWS, 9/30/23: This code was used for testing in 2023
            // but is probably unnecessary now
            //std::cout << "Setting " << next_row << " " << line[3]
            //<< " " << line.size() << " "
            //<< table->get_ncolumns() << std::endl;
            
            // Set the row
            table->set_row(((size_t)next_row),line);
            
          }
        
          // End of 'if (mcmc_accept || store_rejects)'
        }

        critical_extra(i_thread);
      
        // End of critical region
      }
    
      // If necessary, increment the multiplier on the previous point
      if (ret_value==o2scl::success && mcmc_accept==false) {
        if (walker_accept_rows[windex]<0 ||
            walker_accept_rows[windex]>=((int)table->get_nlines())) {
          O2SCL_ERR2("Invalid row for incrementing multiplier in ",
                     "mcmc_para_table::add_line().",o2scl::exc_efailed);
        }
        double mult_old=table->get("mult",walker_accept_rows[windex]);
        if (mult_old<0.5) {
          O2SCL_ERR2("Old multiplier less than 1 in ",
                     "mcmc_para_table::add_line().",o2scl::exc_efailed);
        }
        table->set("mult",walker_accept_rows[windex],mult_old+1.0);
        if (this->verbose>=2) {
          this->scr_out << "mcmc: Updating mult of row "
                        << walker_accept_rows[windex]
                        << " from " << mult_old << " to "
                        << mult_old+1.0 << std::endl;
        }
      
      }
    
      // Increment row counters if necessary
      if (ret_value==o2scl::success) {
        if (mcmc_accept) {
          walker_accept_rows[windex]=next_row;
        } else if (store_rejects && func_ret==0) {
          walker_reject_rows[windex]=next_row;
        }
      }

      
      return ret_value;
    }
    //@}
  
    /** \brief Perform cleanup after an MCMC simulation
     */
    virtual void mcmc_cleanup() {

      // This section removes empty rows at the end of the
      // table that were allocated but not used.
      int i;
      bool done=false;
      for(i=table->get_nlines()-1;i>=0 && done==false;i--) {
        done=true;
        if (fabs(table->get("mult",i))<0.1) {
          done=false;
        } 
      }
      if (i+2<((int)table->get_nlines())) {
        table->set_nlines(i+2);
      }

      write_files(true);

      return parent_t::mcmc_cleanup();
    }

    /** \brief Compute autocorrelation coefficient for column
        with index \c icol averaging over all walkers and 
        all threads
    */
    virtual void ac_coeffs(size_t icol,
                           std::vector<double> &ac_coeff_avg,
                           int loc_verbose=0) {

      ac_coeff_avg.resize(0);

      size_t n_tot=this->n_threads*this->n_walk;
    
      std::vector<std::vector<double> > ac_coeff_temp(n_tot);
      size_t max_size=0;
      for(size_t j=0;j<this->n_threads;j++) {
        for(size_t k=0;k<this->n_walk;k++) {
          size_t tindex=j*this->n_walk+k;
          std::vector<double> vec, mult;
          for(size_t ell=0;ell<table->get_nlines();ell++) {
            if (fabs(table->get("walker",ell)-k)<0.1 &&
                fabs(table->get("thread",ell)-j)<0.1 &&
                table->get("mult",ell)>0.5) {
              mult.push_back(table->get("mult",ell));
              vec.push_back(table->get(icol,ell));
            }
          }
          if (mult.size()>1) {
            o2scl::vector_autocorr_vector_mult
              (vec.size(),vec,mult,ac_coeff_temp[tindex]);
            if (ac_coeff_temp[tindex].size()>max_size) {
              max_size=ac_coeff_temp[tindex].size();
            }
          }
          if (loc_verbose>0) {
            std::cout << "column index, thread, walker, length, mult. sum: "
                      << icol << " " << j << " " << k << " "
                      << mult.size() << " "
                      << o2scl::vector_sum<std::vector<double>,double>
              (mult.size(),mult) << " " << max_size << std::endl;
          }
        }
      }
      for(size_t j=0;j<max_size;j++) {
        double ac_sum=0.0;
        int ac_count=0;
        for(size_t k=0;k<n_tot;k++) {
          if (j<ac_coeff_temp[k].size()) {
            ac_sum+=ac_coeff_temp[k][j];
            ac_count++;
          }
        }
        if (ac_count>0) {
          ac_coeff_avg.push_back(ac_sum/((double)ac_count));
        }
      }
    
      return;
    }

    /** \brief Reorder the table by thread and walker index
     */
    virtual void reorder_table() {

      // Create a new table
      std::shared_ptr<o2scl::table_units<> > table2=
        std::shared_ptr<o2scl::table_units<> >(new o2scl::table_units<>);

      for(size_t i=0;i<this->n_threads;i++) {
        for(size_t j=0;j<this->n_walk;j++) {
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

  /** \brief MCMC class with a command-line interface

      This class forms the basis of the MCMC used in the Bayesian
      analysis of neutron star mass and radius in
      http://github.com/awsteiner/bamr .

  */
  template<class func_t, class fill_t, class data_t, class vec_t=ubvector,
           class stepper_t=mcmc_stepper_base<func_t,data_t,vec_t>>
  class mcmc_para_cli : public mcmc_para_table<func_t,fill_t,
                                               data_t,vec_t,stepper_t> {
    
  protected:
  
    /** \brief The parent typedef
     */
    typedef o2scl::mcmc_para_table<func_t,fill_t,data_t,vec_t,
                                   stepper_t> parent_t;

    /// \name Parameter objects for the 'set' command
    //@{
    o2scl::cli::parameter_double p_step_fac;
    o2scl::cli::parameter_size_t p_n_warm_up;
    o2scl::cli::parameter_int p_user_seed;
    o2scl::cli::parameter_size_t p_max_bad_steps;
    o2scl::cli::parameter_size_t p_n_walk;
    o2scl::cli::parameter_bool p_aff_inv;
    o2scl::cli::parameter_bool p_table_sequence;
    o2scl::cli::parameter_bool p_store_rejects;
    o2scl::cli::parameter_bool p_check_rows;
    o2scl::cli::parameter_bool p_couple_threads;
    o2scl::cli::parameter_double p_max_time;
    o2scl::cli::parameter_size_t p_max_iters;
    //o2scl::cli::parameter_int p_max_chain_size;
    o2scl::cli::parameter_size_t p_file_update_iters;
    o2scl::cli::parameter_double p_file_update_time;
    //o2scl::cli::parameter_bool p_output_meas;
    o2scl::cli::parameter_string p_prefix;
    o2scl::cli::parameter_int p_verbose;
    //@}
  
    /** \brief Initial write to HDF5 file 
     */
    virtual void file_header(o2scl_hdf::hdf_file &hf) {
      hf.sets_vec_copy("cl_args",this->cl_args);
      return;
    }
  
  public:
  
    /** \brief The arguments sent to the command-line
     */
    std::vector<std::string> cl_args;

    /// \name Customization functions
    //@{
    /** \brief Set up the 'cli' object
      
        This function just adds the four commands and the 'set' parameters
    */
    virtual void setup_cli(cli &cl) {
      
      // ---------------------------------------
      // Set commands/options

      /*
        static const size_t nopt=1;
        o2scl::comm_option_s options[nopt]={
        {'i',"initial-point","Set the starting point in the parameter space",
        1,-1,"<mode> [...]",
        ((std::string)"Mode can be one of 'best', 'last', 'N', or ")+
        "'values'. If mode is 'best', then it uses the point with the "+
        "largest weight and the second argument specifies the file. If "+
        "mode is 'last' then it uses the last point and the second "+
        "argument specifies the file. If mode is 'N' then it uses the Nth "+
        "point, the second argument specifies the value of N and the third "+
        "argument specifies the file. If mode is 'values', then the "+
        "remaining arguments specify all the parameter values. On the "+
        "command-line, enclose negative values in quotes and parentheses, "+
        "i.e. \"(-1.00)\" to ensure they do not get confused with other "+
        "options.",new o2scl::comm_option_mfptr<mcmc_para_cli>
        (this,&mcmc_para_cli::set_initial_point),
        o2scl::cli::comm_option_both}
        {'s',"hastings","Specify distribution for M-H step",
        1,1,"<filename>",
        ((string)"Desc. ")+"Desc2.",
        new comm_option_mfptr<mcmc_mpi>(this,&mcmc_mpi::hastings),
        cli::comm_option_both}
        };
        this->cl.set_comm_option_vec(nopt,options);
      */
    
      p_file_update_iters.s=&this->file_update_iters;
      p_file_update_iters.help=((std::string)"Number of MCMC successes ")+
        "between file updates (default 0 for no file updates).";
      cl.par_list.insert(std::make_pair("file_update_iters",
                                        &p_file_update_iters));
    
      p_file_update_time.d=&this->file_update_time;
      p_file_update_time.help=((std::string)"Time ")+
        "between file updates (default 0.0 for no file updates).";
      cl.par_list.insert(std::make_pair("file_update_time",
                                        &p_file_update_time));
    
      /*
        p_max_chain_size.i=&this->max_chain_size;
        p_max_chain_size.help=((std::string)"Maximum Markov chain size ")+
        "(default 10000).";
        cl.par_list.insert(std::make_pair("max_chain_size",
        &p_max_chain_size));
      */
    
      p_max_time.d=&this->max_time;
      p_max_time.help=((std::string)"Maximum run time in seconds ")+
        "(default 86400 sec or 1 day).";
      cl.par_list.insert(std::make_pair("max_time",&p_max_time));
    
      p_max_iters.s=&this->max_iters;
      p_max_iters.help=((std::string)"If non-zero, limit the number of ")+
        "iterations to be less than the specified number (default zero).";
      cl.par_list.insert(std::make_pair("max_iters",&p_max_iters));
    
      p_prefix.str=&this->prefix;
      p_prefix.help="Output file prefix (default 'mcmc\').";
      cl.par_list.insert(std::make_pair("prefix",&p_prefix));

      /*
        p_output_meas.b=&this->output_meas;
        p_output_meas.help=((std::string)"If true, output next point ")+
        "to the '_scr' file before calling TOV solver (default true).";
        cl.par_list.insert(std::make_pair("output_meas",&p_output_meas));
      */
    
      p_step_fac.d=&this->step_fac;
      p_step_fac.help=((std::string)"MCMC step factor. The step size for ")+
        "each variable is taken as the difference between the high and low "+
        "limits divided by this factor (default 10.0). This factor can "+
        "be increased if the acceptance rate is too small, but care must "+
        "be taken, e.g. if the conditional probability is multimodal. If "+
        "this step size is smaller than 1.0, it is reset to 1.0 .";
      cl.par_list.insert(std::make_pair("step_fac",&p_step_fac));

      p_n_warm_up.s=&this->n_warm_up;
      p_n_warm_up.help=((std::string)"Minimum number of warm up iterations ")+
        "(default 0).";
      cl.par_list.insert(std::make_pair("n_warm_up",&p_n_warm_up));

      p_verbose.i=&this->verbose;
      p_verbose.help=((std::string)"MCMC verbosity parameter ")+
        "(default 0).";
      cl.par_list.insert(std::make_pair("mcmc_verbose",&p_verbose));

      p_max_bad_steps.s=&this->max_bad_steps;
      p_max_bad_steps.help=((std::string)"Maximum number of bad steps ")+
        "(default 1000).";
      cl.par_list.insert(std::make_pair("max_bad_steps",&p_max_bad_steps));

      p_n_walk.s=&this->n_walk;
      p_n_walk.help=((std::string)"Number of walkers ")+
        "(default 1).";
      cl.par_list.insert(std::make_pair("n_walk",&p_n_walk));

      p_user_seed.i=&this->user_seed;
      p_user_seed.help=((std::string)"Seed for multiplier for random ")+
        "number generator. If zero is given (the default), then mcmc() "+
        "uses time(0) to generate a random seed.";
      cl.par_list.insert(std::make_pair("user_seed",&p_user_seed));
    
      p_aff_inv.b=&this->aff_inv;
      p_aff_inv.help=((std::string)"If true, then use affine-invariant ")+
        "sampling (default false).";
      cl.par_list.insert(std::make_pair("aff_inv",&p_aff_inv));
    
      p_table_sequence.b=&this->table_sequence;
      p_table_sequence.help=((std::string)"If true, then ensure equal ")+
        "spacing between walkers (default true).";
      cl.par_list.insert(std::make_pair("table_sequence",&p_table_sequence));
    
      p_store_rejects.b=&this->store_rejects;
      p_store_rejects.help=((std::string)"If true, then store MCMC ")+
        "rejections (default false).";
      cl.par_list.insert(std::make_pair("store_rejects",&p_store_rejects));
      
      p_check_rows.b=&this->check_rows;
      p_check_rows.help="If true, then check rows";
      cl.par_list.insert(std::make_pair("check_rows",&p_check_rows));

      p_couple_threads.b=&this->couple_threads;
      p_couple_threads.help="help";
      cl.par_list.insert(std::make_pair("couple_threads",&p_couple_threads));
    
      return;
    }
  
  };
  
  // End of namespace
}

#endif
