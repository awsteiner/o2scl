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
    
    // Position variables
    std::vector<double> position;

    // Momentum variables
    std::vector<double> momentum;

    // Scaling parameters
    std::vector<double> scale;

    // Step size
    std::vector<double> step_size;

    // Trajectory length
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
    virtual int hmc(size_t n_params, vec_t &low, vec_t &high,
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

      // ----------------------------------------------------------------
      // Set storage of return values from each thread
      std::vector<int> func_ret(n_threads), meas_ret(n_threads);



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