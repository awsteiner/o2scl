/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#ifndef O2SCL_ANNEAL_MT_H
#define O2SCL_ANNEAL_MT_H

#include <o2scl/rng_gsl.h>
#include <o2scl/multi_funct.h>
#include <o2scl/anneal.h>
#include <o2scl/vector.h>

// for boost::bind()
#include <boost/bind.hpp>
// for boost::thread
#include <boost/thread/thread.hpp>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Multidimensional minimization by simulated annealing 
      (Boost multi-threaded version)

      This header-only class additionally requires the Boost
      libraries. It performs simulated annealing using an arbitrary
      number of processors using <tt>boost::thread</tt>, which is
      closely related to the standard Unix pthread library. It works
      very similarly to \ref anneal_gsl, it performs \ref ntrial
      evaluations over each processor, then applying the metropolis
      algorithm to the results from all of the processors at the end.
     
      Because <tt>np</tt> function calls are happening simultaneously,
      where <tt>np</tt> is the number of processors, <tt>np</tt>
      copies of the function parameters of type <tt>param_t</tt> must
      also be specified. The user-defined function to minimize must
      also be thread-safe, allowing multiple threads to call the
      function at once (albeit given different parameters). The
      default type for these <tt>np</tt> copies of the parameters of
      type <tt>param_t</tt> is <tt>std::vector<param_t></tt>.

      This works particularly well for functions which are not trivial
      to evaluate, i.e. functions where the execution time is more
      longer than the bookkeeping that \ref anneal_mt performs between
      trials. For functions which satisfy this requirement, this
      algorithm scales nearly linearly with the number of processors.

      Verbose I/O for this class happens only outside the theads
      unless the user places I/O in the streams in the function that
      is specified.

      \future There may be a good way to remove the function indirection
      here to make this class a bit faster.
  */
  template<class func_t=multi_funct,
    class vec_t=boost::numeric::ublas::vector<double>,
    class rng_t=int, class rng_dist_t=rng_gsl>
    class anneal_mt : public anneal_base<func_t,vec_t,rng_t,rng_dist_t> {

  public:
  
  anneal_mt() {
    ntrial=100;
    nproc=1;
    outs=&std::cout;
    ins=&std::cin;
    tolx=1.0e-6;
    T_start=1.0;
    T_dec=1.5;
    step_dec=1.5;
    min_step_ratio=100.0;
    step_vec.resize(1);
    step_vec[0]=1.0;
    out_best=false;
    out_step_changes=false;
  }

  virtual ~anneal_mt() {
  }

  /// \name Basic usage
  //@{
  /** \brief Calculate the minimum \c fmin of \c func w.r.t the 
      array \c x0 of size \c nv using \c np threads.
  */
  virtual int mmin(size_t nv, vec_t &x0, double &fmin, 
		   func_t &func, size_t np) {
    
    if (nv==0) {
      O2SCL_ERR2("Tried to minimize over zero variables ",
		     " in anneal_mt::mmin().",o2scl::exc_einval);
    }
    if (np==0) {
      O2SCL_ERR2("Tried to use zero threads in ",
		     "anneal_mt::mmin().",o2scl::exc_einval);
    }

    allocate(nv,np);
    f=&func;
    fmin=0.0;

    double E, best_E, T, old_E;
    int i, iter=0;
    size_t j;

    for(j=0;j<nv;j++) {
      x[j]=x0[j];
      best_x[j]=x0[j];
    }
    
    E=func(nv,x);
    best_E=E;

    // Setup initial temperature and step sizes
    start(nv,T);

    bool done=false;
    
    boost::thread **thrd;
    
    while (!done) {

      // Copy old value of x for next() function

      for(j=0;j<nv;j++) old_x[j]=x[j];
      old_E=E;
	  
      size_t nmoves=0;

      for (i=0;i<ntrial;++i) {

	// Determine the stepsize, create and execute the new threads
	thrd=new boost::thread *[np];
	for(size_t ip=0;ip<np;ip++) {
	  new_x[ip]=x;
	  thrd[ip]=new boost::thread
	    (boost::bind(&anneal_mt::func_wrapper,this,ip));
	}
	// Wait until all the threads are done
	for(size_t ip=0;ip<np;ip++) {
	  thrd[ip]->join();
	}
	// Delete the threads and continue
	for(size_t ip=0;ip<np;ip++) {
	  delete thrd[ip];
	}
	delete[] thrd;
	// Process the results from each thread
	for(size_t ip=0;ip<np;ip++) {

	  // Store best value obtained so far
	  if(new_E[ip]<=best_E){
	    for(j=0;j<nv;j++) best_x[j]=new_x[ip][j];
	    best_E=new_E[ip];
	    if (this->verbose>0 && out_best) {
	      std::cout << "Best: ";
	      o2scl::vector_out(std::cout,best_x);
	      std::cout << " " << best_E << std::endl;
	    }
	  }

	  // Take the crucial step: see if the new point is accepted
	  // or not, as determined by the boltzmann probability
	  if (new_E[ip]<E) {
	    for(j=0;j<nv;j++) x[j]=new_x[ip][j];
	    E=new_E[ip];
	    nmoves++;
	  } else if (this->rng_dist(this->rng)
		     < exp(-(new_E[ip]-E)/(boltz*T)) ) {
	    for(j=0;j<nv;j++) x[j]=new_x[ip][j];
	    E=new_E[ip];
	    nmoves++;
	  }
	}
	
      }
      
      if (this->verbose>0) {
	print_iter(nv,best_x,best_E,iter,T,"anneal_mt");
	iter++;
      }
      
      // See if we're finished and proceed to the next step
      next(nv,old_x,old_E,x,E,T,nmoves,done);
      
    }
    
    for(j=0;j<nv;j++) x0[j]=best_x[j];
    fmin=best_E;

    return 0;
  }
  //@}

  /** \brief Desc
   */
  virtual int mmin(size_t nv, vec_t &x0, double &fmin, 
		   func_t &func) {
    return mmin(nv,x0,fmin,func,1);
  }

  /// \name Iteration control
  //@{
  /// Determine how to change the minimization for the next iteration
  virtual int next(size_t nv, vec_t &x_old, double min_old, vec_t &x_new, 
		   double min_new, double &T, size_t n_moves,
		   bool &finished) {
	
    if (T/T_dec<this->tolx) {
      finished=true;
      return o2scl::success;
    }
    if (n_moves==0) {
      bool changed=false;
      for(size_t i=0;i<nv;i++) {
	if (i<step_vec.size() && step_vec[i]>this->tolx*min_step_ratio) {
	  step_vec[i]/=step_dec;
	  changed=true;
	}
      }
      if (changed && verbose>0 && out_step_changes) {
	std::cout << "Step sizes changed: ";
	o2scl::vector_out(std::cout,step_vec,true);
      }
    }
    T/=T_dec;
    return o2scl::success;
  }

  /// Setup initial temperature and stepsize
  virtual int start(size_t nv, double &T) {
    T=T_start;
    return o2scl::success;
  }
  //@}

  /// \name Parameters
  //@{
  /// Boltzmann factor (default 1.0).
  double boltz;

  /// Number of iterations
  int ntrial;

  /// Output control
  int verbose;

  /// Output step size changes (default false)
  bool out_step_changes;

  /// Output best point (default false)
  bool out_best;
  
  /// The independent variable tolerance (default \f$ 10^{-6} \f$ )
  double tolx;

  /// Initial temperature (default 1.0)
  double T_start;

  /// Factor to decrease temperature by (default 1.5)
  double T_dec;

  /// Factor to decrease step size by (default 1.5)
  double step_dec;

  /// Ratio between minimum step size and \ref tolx (default 100.0)
  double min_step_ratio;
  //@}

  /// The default random number generator
  rng_t def_rng;

  /// Return string denoting type ("anneal_mt")
  virtual const char *type() { return "anneal_mt"; }

  /// Set streams for verbose I/O
  int set_verbose_stream(std::ostream &out, std::istream &in) {
    outs=&out;
    ins=&in;
    return 0;
  }
    
  /** \brief Print out iteration information.
      
      Depending on the value of the variable verbose, this prints out
      the iteration information. If verbose=0, then no information is
      printed. If verbose>0, then after each iteration, the present
      values of x and y are output to std::cout along with the
      iteration number. Also, if verbose>0, every time a new smallest
      function value is found, the location and the function value is
      output. If verbose>=2 then each iteration waits for a character
      between each trial.
  */
  virtual int print_iter(size_t nv, vec_t &xx, double y, int iter,
			 double tptr, std::string comment) 
  {
    if (this->verbose<=0) return 0;
    
    size_t i;
    char ch;

    (*this->outs) << comment << " Iteration: " << iter << std::endl;
    std::cout << "x: ";
    for(i=0;i<nv;i++) std::cout << x[i] << " ";
    std::cout << std::endl;
    (*this->outs) << "y: " << y << " Tptr: " << tptr << std::endl;
    if (this->verbose>1) {
      (*this->outs) << "Press a key and type enter to continue. ";
      (*this->ins) >> ch;
    }
    
    return 0;
  }
  
  /// Set the step sizes 
  template<class vec2_t> int set_step(size_t nv, vec2_t &stepv) {
    if (nv>0) {
      step_vec.resize(nv);
      for(size_t i=0;i<nv;i++) step_vec[i]=stepv[i];
    }
    return 0;
  }

#ifndef DOXYGEN_INTERNAL
      
  protected:

  /// The function wrapper executed by thread with index \c ip
  void func_wrapper(size_t ip) {
    step(nvar,new_x[ip]);
    new_E[ip]=(*f)(nvar,new_x[ip]);
  }
  
  /// Stream for verbose output
  std::ostream *outs;
  
  /// Stream for verbose input
  std::istream *ins;

  /// The number of threads to run
  size_t nproc;

  /// The number of variables over which we minimize
  size_t nvar;

  /// The function to minimize
  func_t *f;
  
  /// \name Storage for present, next, and best vectors
  //@{
  vec_t x, best_x, new_E, old_x;
  std::vector<vec_t> new_x;
  //@}

  /// Vector of step sizes
  vec_t step_vec;
      
  /** \brief Allocate memory for a minimizer over \c n dimensions
      with stepsize \c step
  */
  virtual int allocate(size_t nv, size_t np) {
    nvar=nv;
    nproc=np;
    x.resize(nvar);
    new_x.resize(np);
    old_x.resize(nvar);
    for(size_t i=0;i<np;i++) {
      new_x[i].resize(nvar);
    }
    new_E.resize(np);
    best_x.resize(nvar);
    return 0;
  }

  /// Make a step to a new attempted minimum
  virtual int step(size_t nv, vec_t &sx) {
    size_t nstep=step_vec.size();
    for(size_t i=0;i<nv;i++) {
      double u=this->rng_dist(this->rng);
      sx[i]=(2.0*u-1.0)*step_vec[i%nstep]+sx[i];
    }
    return 0;
  }
  
#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
