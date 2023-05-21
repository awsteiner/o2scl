/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2017-2023, Andrew W. Steiner
  
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
#ifndef O2SCL_ANNEAL_PARA_H
#define O2SCL_ANNEAL_PARA_H

#ifdef O2SCL_OPENMP
#include <omp.h>
#endif

#include <o2scl/anneal_gsl.h>
#include <o2scl/vector.h>

namespace o2scl {
  
  /** \brief Multidimensional minimization by simulated annealing 
      (OpenMP/MPI version)

      This class works particularly well for functions which are not
      trivial to evaluate, i.e. functions where the execution time is
      more longer than the bookkeeping that \ref anneal_para performs
      between trials. For functions which satisfy this requirement,
      this algorithm scales nearly linearly with the number of
      processors.

      Verbose I/O for this class happens only outside the parallel
      regions unless the user places I/O in the streams in the
      function that is specified.
  */
  template<class func_t=multi_funct,
    class vec_t=boost::numeric::ublas::vector<double> >
    class anneal_para : public anneal_gsl<func_t,vec_t> {

  public:

  typedef boost::numeric::ublas::vector<double> ubvector;

  /// The number of OpenMP threads
  size_t n_threads;
  
  /// The MPI starting time
  double start_time;

  /// The maximum time
  double max_time;

  /** \brief If true, obtain the global minimum over all MPI ranks
      (default true)
  */
  bool collect_all_ranks;
  
  anneal_para() {
    n_threads=1;
    start_time=0.0;
    max_time=0.0;
    collect_all_ranks=true;
  }

  virtual ~anneal_para() {
  }

    std::vector<rng<> > vrng;
  
  /** \brief Make a step to a new attempted minimum
   */
  virtual int step_para(vec_t &x, vec_t &sx, int nvar, size_t ithread) {
    size_t nstep=this->step_vec.size();
    for(int i=0;i<nvar;i++) {
      double u=vrng[ithread].random();

      // Construct the step in the ith direction
      double step_i=this->step_norm*this->step_vec[i%nstep];
      // Fix if the step is too small
      if (step_i<this->tol_abs*this->min_step_ratio) {
	step_i=this->tol_abs*this->min_step_ratio;
      }
      
      sx[i]=(2.0*u-1.0)*step_i+x[i];
    }
    return 0;
  }

  /// \name Basic usage
  //@{
  /** \brief Desc
   */
  virtual int mmin(size_t nv, vec_t &x0, double &fmin, 
		   std::vector<func_t> &func) {

    // Check that we have at least one variable
    if (nv==0) {
      O2SCL_ERR2("Tried to minimize over zero variables ",
		 " in anneal_para::mmin().",exc_einval);
    }
    
    // Check that enough function objects were passed
    if (func.size()<n_threads) {
      if (this->verbose>0) {
	std::cout << "mcmc_para::mcmc(): Not enough functions for "
	<< n_threads << " threads. Setting n_threads to "
	<< func.size() << "." << std::endl;
      }
      n_threads=func.size();
    }

    // Set number of threads
#ifdef O2SCL_OPENMP
    omp_set_num_threads(n_threads);
#pragma omp parallel
    {
      n_threads=omp_get_num_threads();
    }
#else
    n_threads=1;
#endif
    
    // Set starting time
#ifdef O2SCL_MPI
    start_time=MPI_Wtime();
#else
    start_time=time(0);
#endif

    // Initial value of step normalization
    this->step_norm=1.0;

    std::vector<vec_t> x(n_threads), new_x(n_threads);
    std::vector<vec_t> old_x(n_threads), best_x(n_threads);
    
    ubvector E(n_threads), new_E(n_threads), best_E(n_threads);
    ubvector old_E(n_threads), r(n_threads);
    
    double T;
    
    int iter=0;

    /// A different random number generator for each OpenMP thread
    vrng.resize(n_threads);

    // Seed the random number generators
    unsigned long int s=time(0);
    for(size_t it=1;it<n_threads;it++) {
      vrng[it].set_seed(s+it);
    }

    // Setup initial temperature and step sizes
    this->start(nv,T);
    
#ifdef O2SCL_OPENMP
#pragma omp parallel default(shared)
#endif
    {
#ifdef O2SCL_OPENMP
#pragma omp for
#endif
      for(size_t it=0;it<n_threads;it++) {

	// Allocate space
	x[it].resize(nv);
	new_x[it].resize(nv);
	best_x[it].resize(nv);
	old_x[it].resize(nv);
	
	// Copy initial point
	for(size_t j=0;j<nv;j++) {
	  x[it][j]=x0[j];
	  best_x[it][j]=x0[j];
	}

	// Perform first function evaluation
	E[it]=func[it](nv,x[it]);
	best_E[it]=E[it];
      }
    }
    // End of parallel region

    bool done=false;
    while (!done) {

      // Count the number of moves accepted for thread 0
      size_t nmoves=0;

#ifdef O2SCL_OPENMP
#pragma omp parallel default(shared)
#endif
      {
#ifdef O2SCL_OPENMP
#pragma omp for
#endif
	for(size_t it=0;it<n_threads;it++) {

	  // Copy old value of x for next() function
	  for(size_t j=0;j<nv;j++) old_x[it][j]=x[it][j];
	  old_E[it]=E[it];
	  
	  for (int i=0;i<this->ntrial;++i) {
	    
	    step_para(x[it],new_x[it],nv,it);
	    
	    new_E[it]=func[it](nv,new_x[it]);
#ifdef O2SCL_MPI
	    double elapsed=MPI_Wtime()-start_time;
#else
	    double elapsed=time(0)-start_time;
#endif
	    if (max_time>0.0 && elapsed>max_time) {
	      i=this->ntrial;
	      done=true;
	    }

	    // Store best value obtained so far
	    if(new_E[it]<=best_E[it]) {
	      for(size_t j=0;j<nv;j++) best_x[it][j]=new_x[it][j];
	      best_E[it]=new_E[it];
	    }
	
	    if (done==false) {
	      // Take the crucial step: see if the new point is accepted
	      // or not, as determined by the Boltzmann probability
	      if (new_E[it]<E[it]) {
		for(size_t j=0;j<nv;j++) x[it][j]=new_x[it][j];
		E[it]=new_E[it];
		if (it==0) nmoves++;
	      } else {
		r[it]=vrng[it].random();
		if (r[it] < exp(-(new_E[it]-E[it])/(this->boltz*T))) {
		  for(size_t j=0;j<nv;j++) x[it][j]=new_x[it][j];
		  E[it]=new_E[it];
		  if (it==0) nmoves++;
		}
	      }
	    }
	  }
	  
	}
      }
      // End of parallel region

      // Find best point over all threads
      fmin=best_E[0];
      for(size_t iv=0;iv<nv;iv++) {
	x0[iv]=best_x[0][iv];
      }
      for(size_t i=1;i<n_threads;i++) {
	if (best_E[i]<fmin) {
	  fmin=best_E[i];
	  for(size_t iv=0;iv<nv;iv++) {
	    x0[iv]=best_x[i][iv];
	  }
	}
      }
      
      if (this->verbose>0) {
	this->print_iter(nv,x0,fmin,iter,T,"anneal_gsl");
	iter++;
	if (this->verbose>1) {
	  char ch;
	  std::cin >> ch;
	}
      }
	  
      // See if we're finished and proceed to the next step
      if (done==false) {
	
	// std::vector<bool> doesn't work here so we use an
	// integer type instead
	std::vector<int> done_arr(n_threads);
	
#ifdef O2SCL_OPENMP
#pragma omp parallel default(shared)
#endif
	{
#ifdef O2SCL_OPENMP
#pragma omp for
#endif
	  for(size_t it=0;it<n_threads;it++) {
	    done_arr[it]=0;
	    next_para(nv,old_x[it],old_E[it],x[it],E[it],T,nmoves,
		 x0,fmin,done_arr[it]);
	  }
	}
	// End of parallel region

	// Decrease the temperature
	T/=this->T_dec;
	
	// We're done when all threads report done
	done=true;
	for(size_t it=0;it<n_threads;it++) {
	  if (done_arr[it]==0) done=false;
	}
      }
      
    }

#ifdef O2SCL_MPI

    if (collect_all_ranks) {
      
      int mpi_rank, mpi_size;
      // Get MPI rank and size
      MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
      MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
      if (mpi_size>1) {

	double fmin_new;
	vector<double> x0_new(nv);

	if (mpi_rank<size-1) {

	  // Get the minimum from the higher rank
	  int tag=0;
	  MPI_Recv(&fmin_new,1,MPI_DOUBLE,mpi_rank+1,tag,MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  tag=1;
	  MPI_Recv(&(x0_new[0]),nv,MPI_DOUBLE,mpi_rank+1,tag,MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);

	  // Update x0 and fmin if necessary
	  if (fmin_new<fmin) {
	    for(size_t ik=0;ik<nv;ik++) {
	      x0[ik]=x0_new[ik];
	    }
	    fmin=fmin_new;
	  } else {
	    for(size_t ik=0;ik<nv;ik++) {
	      x0_new[ik]=x0[ik];
	    }
	  }
	}
	
	if (mpi_size>1 && mpi_rank>0) {
	
	  // Copy the minimum into the buffer
	  for(size_t ik=0;ik<nv;ik++) x0_new[ik]=x0[ik];

	  // Send the minimum down to the lower rank
	  MPI_Send(&fmin,1,MPI_DOUBLE,mpi_rank-1,0,MPI_COMM_WORLD);
	  MPI_Send(&(x0_new[0]),nv,MPI_DOUBLE,mpi_rank-1,1,MPI_COMM_WORLD);
	}

	if (mpi_size>1 && mpi_rank>0) {
	  // Obtain the global minimum from the lower rank
	  int tag=2;
	  MPI_Recv(&(fmin_new),1,MPI_DOUBLE,mpi_rank-1,tag,MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);
	  tag=3;
	  MPI_Recv(&(x0_new[0]),nv,MPI_DOUBLE,mpi_rank-1,tag,MPI_COMM_WORLD,
		   MPI_STATUS_IGNORE);

	  // Update x0 and fmin if necessary
	  if (fmin_new<fmin) {
	    for(size_t ik=0;ik<nv;ik++) {
	      x0[ik]=x0_new[ik];
	    }
	    fmin=fmin_new;
	  }
	}
      
	if (mpi_rank<size-1) {
	  // Send the minimum back to the higher rank
	  MPI_Send(&fmin,1,MPI_DOUBLE,mpi_rank+1,2,MPI_COMM_WORLD);
	  MPI_Send(&(x0_new[0]),nv,MPI_DOUBLE,mpi_rank+1,3,MPI_COMM_WORLD);
	}
      
	// Ensure that x0_new isn't deallocated before the communication is
	// done
	MPI_Barrier(MPI_COMM_WORLD);
      
      }

    }
    
#endif
    
    return 0;
  }

  /// Determine how to change the minimization for the next iteration
  virtual int next_para(size_t nvar, vec_t &x_old, double min_old, 
		   vec_t &x_new, double min_new, double &T, 
		   size_t n_moves, vec_t &best_x, double best_E, 
		   int &finished) {
    
    if (T/this->T_dec<this->tol_abs) {
      finished=1;
      return success;
    }
    if (n_moves==0) {
      // If we haven't made any progress, shrink the step by
      // decreasing step_norm
      this->step_norm/=this->step_dec;
      // Also reset x to best value so far
      for(size_t i=0;i<nvar;i++) {
	x_new[i]=best_x[i];
      }
      min_new=best_E;
    }
    return success;
  }

  /** \brief Desc
   */
  virtual int mmin(size_t nv, vec_t &x0, double &fmin, 
		   func_t &func) {
#ifdef O2SCL_OPENMP
    omp_set_num_threads(n_threads);
#pragma omp parallel
    {
      n_threads=omp_get_num_threads();
    }
#else
    n_threads=1;
#endif
    std::vector<func_t> vf(n_threads);
    for(size_t i=0;i<n_threads;i++) {
      vf[i]=func;
    }
    return mmin(nv,x0,fmin,vf);

  }

  /// Return string denoting type ("anneal_para")
  virtual const char *type() { return "anneal_para"; }

  };

}

#endif
