/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2016, Andrew W. Steiner
  
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
/** \file mcmc_mpi.h
    \brief Definition of \ref o2scl::mcmc_mpi and \ref o2scl::mcmc_cli 
    classes
*/
#ifndef MCMC_MPI_H
#define MCMC_MPI_H

#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#ifdef O2SCL_MPI
#include <mpi.h>
#endif

#include <o2scl/rng_gsl.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/table3d.h>
#include <o2scl/hdf_file.h>
#include <o2scl/hdf_io.h>
#include <o2scl/exception.h>
#include <o2scl/prob_dens_func.h>
#include <o2scl/cholesky.h>
#include <o2scl/vector.h>
#include <o2scl/multi_funct.h>
#include <o2scl/mcmc.h>

#ifdef O2SCL_READLINE
#include <o2scl/cli_readline.h>
#else
#include <o2scl/cli.h>
#endif

namespace o2scl {
  
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  
  /** \brief MCMC with MPI and HDF5 table I/O 

      The parent, \ref o2scl::mcmc_table assembles the MCMC data into
      objects of type \ref o2scl::table_units . This class
      additionally writes the MCMC data to a HDF5 table and controls
      run time on several processors via MPI.
   */
  template<class func_t, class fill_t, class data_t,
    class vec_t=ubvector> class mcmc_mpi :
    public o2scl::mcmc_table<func_t,fill_t,data_t,vec_t> {
    
  protected:

  typedef o2scl::mcmc_table<func_t,fill_t,data_t,vec_t> parent_t;
  
  /// \name MPI properties
  //@{
  /// The MPI processor rank
  int mpi_rank;

  /// The MPI number of processors
  int mpi_nprocs;

  /// The MPI starting time
  double mpi_start_time;
  //@}
  
  /** \brief Error handler for each thread
   */
  o2scl::err_hnd_cpp error_handler;

  /** \brief The number of parameters
   */
  size_t n_params;

  /** \brief A copy of the lower limits for HDF5 output
   */
  vec_t low_copy;

  /** \brief A copy of the upper limits for HDF5 output
   */
  vec_t high_copy;
  
  public:

  /** \brief Default method for setting the random seed
   */
  virtual void set_seed() {
    // Set RNG seed
    unsigned long int seed=time(0);
    if (this->user_seed!=0) {
      seed=this->user_seed;
    }
    seed*=(mpi_rank+1);
    this->rg.set_seed(seed);

    return;
  }

  /** \brief Perform an MCMC simulation
   */
  virtual int mcmc(size_t np, vec_t &init,
		   vec_t &low, vec_t &high, func_t &func,
		   fill_t &fill) {

    low_copy=low;
    high_copy=high;
    
    // The function mcmc_init() needs to know the number of
    // parameters, so we store it here for later use
    n_params=np;
    return parent_t::mcmc(np,init,low,high,func,fill);
  }
    
  public:
  
  /// The screen output file
  std::ofstream scr_out;
  
  /// If true, scr_out has been opened
  bool file_opened;

  /// \name Command-line settings
  //@{
  /** \brief Maximum number of iterations (default 0)
   */
  size_t max_iters;
  
  /** \brief Time in seconds (default is 86400 seconds or 1 day)
   */
  double max_time;

  /** \brief Output each measurement
   */
  bool output_meas;

  /** \brief Prefix for output filenames
   */
  std::string prefix;

  /// If true, output MC accepts and rejects (default true)
  //@}

  /// The first point in the parameter space
  ubvector initial_point;
    
  /// The file containing the initial point
  std::string initial_point_file;

  /// \name Integer designating how to set the initial point
  //@{
  int initial_point_type;
  static const int fp_unspecified=-1;
  static const int fp_last=-2;
  static const int fp_best=-3;
  //@}

  /// If true, then \ref first_update() has been called
  bool first_file_update;
  
  /// \name Command-line settings
  //@{
  /** \brief The number of MCMC successes between file updates
      (default 40)
  */
  int file_update_iters;

  /** \brief Maximum size of Markov chain (default 10000)
   */
  int max_chain_size;
  //@}
    
  /// Number of complete Markov chain segments
  size_t chain_index;
  
  /** \brief Update files with current table
   */
  virtual void update_files() {

    o2scl_hdf::hdf_file hf;
    
    // Open main update file
    hf.open_or_create(prefix+"_"+o2scl::itos(mpi_rank)+"_out");
    
    // First time, output some initial quantities
    if (first_file_update==false) {
      first_update(hf);
      first_file_update=true;
    }

    hf.set_szt("n_accept",this->n_accept);
    hf.set_szt("n_reject",this->n_reject);
    hf.set_szt("n_chains",chain_index+1);
    hf.set_szt_vec("ret_value_counts",this->ret_value_counts);
    
    std::string ch_name="markov_chain"+o2scl::szttos(chain_index);
    hdf_output(hf,*this->tab,ch_name);

    hf.close();
    
    return;
  }
  
  /** \brief Add a measurement to the table
   */
  virtual int add_line(const ubvector &pars, double weight,
		       size_t ix, bool new_meas, data_t &dat,
		       fill_t &fill) {

    if (output_meas) {
      scr_out << "Line: ";
      o2scl::vector_out(scr_out,pars);
      scr_out << " " << weight << " " << ix << " " << new_meas << std::endl;
    }
    
    o2scl::mcmc_table<func_t,fill_t,data_t,ubvector>::add_line
    (pars,weight,ix,new_meas,dat,fill);
      
    bool files_updated=false;
    if (((int)this->tab->get_nlines())==max_chain_size) {
      scr_out << "Writing files and starting new chain." << std::endl;
      update_files();
      chain_index++;
      this->tab->clear_data();
      files_updated=true;
    } else if (this->n_accept>0 && this->n_reject>0 &&
	       (this->n_accept+this->n_reject)%file_update_iters==0) {
      update_files();
      files_updated=true;
    }
    
    if (this->max_iters>0 &&
	(this->n_accept+this->n_reject)==this->max_iters) {
      if (files_updated==false) {
	update_files();
      }
      if (this->verbose>=1) {
	std::cout << "mcmc: Stopping because n_accept+n_reject, "
	<< this->n_accept+this->n_reject
	<< ", is equal to max_iters." << std::endl;
      }
      return this->mcmc_done;
    } else {
      // Determine time elapsed
#ifdef O2SCL_MPI
      double elapsed=MPI_Wtime()-mpi_start_time;
#else
      double elapsed=time(0)-mpi_start_time;
#endif
      if (elapsed>max_time) {
	if (files_updated==false) {
	  update_files();
	}
	if (this->verbose>=1) {
	  std::cout << "mcmc: Stopping because elapsed > max_time."
		    << std::endl;
	}
	return this->mcmc_done;
      }
    }

    return 0;
  }
    
  /// \name Customization functions
  //@{
  /** \brief User-defined initialization function
   */
  virtual int mcmc_init() {
    
    if (this->verbose>=2) {
      std::cout << "Start mcmc_mpi::mcmc_init()." << std::endl;
    }
    
    o2scl::mcmc_table<func_t,fill_t,data_t,vec_t>::mcmc_init();
    
#ifdef O2SCL_MPI
    mpi_start_time=MPI_Wtime();
#else
    mpi_start_time=time(0);
#endif
    
    chain_index=0;

    if (this->file_opened==false) {
      // Open main output file
      this->scr_out.open((this->prefix+"_"+
			  o2scl::itos(this->mpi_rank)+"_scr").c_str());
      this->scr_out.setf(std::ios::scientific);
      this->file_opened=true;
      this->scr_out << "Opened main file in command 'mcmc'." << std::endl;
    }
      
    // Fix file_update_iters if necessary
    if (file_update_iters<1) {
      this->scr_out << "Parameter 'file_update_iters' less "
      << "than 1. Set equal to 1." << std::endl;
      file_update_iters=1;
    }
      
    if (max_chain_size<1) {
      O2SCL_ERR("Parameter 'max_chain_size' must be larger than 1.",
		o2scl::exc_einval);
    }

#ifdef O2SCL_MPI
    int buffer=0, tag=0;
    if (this->mpi_nprocs>1 && this->mpi_rank>0) {
      MPI_Recv(&buffer,1,MPI_INT,this->mpi_rank-1,tag,MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
    }
#endif

    if (initial_point_file.length()>0) {
  
      if (initial_point_type==fp_last) {

	// Read file 
	this->scr_out << "Reading last point from file '" << initial_point_file
	<< "'." << std::endl;
	o2scl_hdf::hdf_file hf;
	hf.open(initial_point_file);
      
	// Read table
	size_t file_n_chains;
	hf.get_szt("n_chains",file_n_chains);
	std::string chain_name=std::string("markov_chain")+
	o2scl::szttos(file_n_chains-1);
	o2scl::table_units<> file_tab;
	hdf_input(hf,file_tab,chain_name);
	size_t last_line=file_tab.get_nlines()-1;

	// Get parameters
	for(size_t i=0;i<n_params;i++) {
	  std::string pname=((std::string)"param_")+this->col_names[i];
	  this->current[0][i]=file_tab.get(pname,last_line);
	  this->scr_out << "Parameter named "
			<< this->col_names[i] << " " 
			<< this->current[0][i] << std::endl;
	}
      
	// Finish up
	this->scr_out << std::endl;
	hf.close();

      } else if (initial_point_type==fp_best) {
	
	std::vector<double> best_point;
	o2scl_hdf::hdf_file hf;
	hf.open(initial_point_file);
	hf.getd_vec("best_point",best_point);
	hf.close();
	this->scr_out << "Reading best point from file '" << initial_point_file
	<< "'." << std::endl;
	for(size_t i=0;i<n_params;i++) {
	  this->current[0][i]=best_point[i];
	  this->scr_out << "Parameter " << i << " : "
			<< this->current[0][i] << std::endl;
	}
	this->scr_out << "Best weight: "
	<< best_point[n_params] << std::endl;
	this->scr_out << std::endl;

      } else {

	// Read file 
	this->scr_out << "Reading "
	<< initial_point_type << "th point from file '" 
	<< initial_point_file
	<< "'." << std::endl;
	o2scl_hdf::hdf_file hf;
	hf.open(initial_point_file);
      
	// Read table
	size_t file_n_chains, row=initial_point_type;
	hf.get_szt("n_chains",file_n_chains);
      
	o2scl::table_units<> file_tab;
	for(size_t k=0;k<file_n_chains;k++) {
	  std::string chain_name=std::string("markov_chain")+o2scl::szttos(k);
	  hdf_input(hf,file_tab,chain_name);
	  if (file_tab.get_nlines()>row) {
	    k=file_n_chains;
	  } else {
	    row-=file_tab.get_nlines();
	  }
	}
	if (row>=file_tab.get_nlines()) {
	  this->scr_out << "Couldn't find point " << initial_point_type 
	  << " in file. Using last point." << std::endl;
	  row=file_tab.get_nlines()-1;
	}
      
	// Get parameters
	for(size_t i=0;i<n_params;i++) {
	  std::string pname=((std::string)"param_")+this->col_names[i];
	  this->current[0][i]=file_tab.get(pname,row);
	  this->scr_out << "Parameter named "
	  << this->col_names[i] << " " 
	  << this->current[0][i] << std::endl;
	}
      
	// Finish up
	this->scr_out << std::endl;
	hf.close();
      }

    } else if (initial_point.size()>0) {
    
      this->scr_out << "First point from command-line." << std::endl;
      for(size_t i=0;i<n_params;i++) {
	this->current[0][i]=initial_point[i];
	this->scr_out << this->current[0][i] << std::endl;
      }
      this->scr_out << std::endl;

    } else {
      
      this->scr_out << "First point from default." << std::endl;
      
    }

#ifdef O2SCL_MPI
    if (this->mpi_nprocs>1 && this->mpi_rank<this->mpi_nprocs-1) {
      MPI_Send(&buffer,1,MPI_INT,this->mpi_rank+1,tag,MPI_COMM_WORLD);
    }
#endif

    this->scr_out << "First point: ";
    o2scl::vector_out(this->scr_out,this->current[0],true);

    if (this->verbose>=2) {
      std::cout << "End mcmc_mpi::mcmc_init()." << std::endl;
    }
    
    return 0;
  };

  /** \brief Output the best point so far
   */
  virtual void best_point(ubvector &best, double w_best) {
    if (this->file_opened==false) {
      // Open main output file
      this->scr_out.open((this->prefix+"_"+
			  o2scl::itos(this->mpi_rank)+"_scr").c_str());
      this->scr_out.setf(std::ios::scientific);
      this->file_opened=true;
      this->scr_out << "Opened main file in function 'best_point()'."
		    << std::endl;
    }
    this->scr_out << "Best: ";
    o2scl::vector_out(this->scr_out,best);
    this->scr_out << " " << w_best << std::endl;
    return;
  }

  /** \brief Initial write to HDF5 file 
   */
  virtual void first_update(o2scl_hdf::hdf_file &hf) {
    
    hf.sets_vec("param_names",this->col_names);
    
    hf.set_szt("n_params",n_params);
    hf.setd("max_time",this->max_time);
    hf.seti("user_seed",this->user_seed);
    hf.seti("n_warm_up",this->n_warm_up);
    hf.setd("step_fac",this->step_fac);
    hf.seti("max_iters",this->max_iters);
    hf.seti("file_update_iters",file_update_iters);
    hf.seti("output_meas",this->output_meas);
    hf.seti("initial_point_type",initial_point_type);
    hf.sets("initial_point_file",initial_point_file);
    hf.setd_vec_copy("initial_point",initial_point);
    hf.set_szt("n_chains",chain_index+1);
    hf.set_szt_vec("ret_value_counts",this->ret_value_counts);
    hf.setd_vec_copy("low",this->low_copy);
    hf.setd_vec_copy("high",this->high_copy);
    
    return;
  }
  
#ifdef NEVER_DEFINED
  /** \brief Choose a Metropolis-Hastings step
   */
  virtual int hastings(std::vector<std::string> &sv, 
		       bool itive_com) {

    bool debug=true;

    if (file_opened==false) {
      // Open main output file
      scr_out.open((prefix+"_"+o2scl::itos(mpi_rank)+"_scr").c_str());
      scr_out.setf(ios::scientific);
      file_opened=true;
      scr_out << "Opened main file in command 'hastings'." << endl;
    }

    if (sv.size()<2) {
      cout << "No arguments given to 'hastings'." << endl;
      return exc_efailed;
    }

    if (model_type.length()==0) {
      cout << "No model selected in 'hastings'." << endl;
      return exc_efailed;
    }

#ifdef O2SCL_MPI
    int buffer=0, tag=0;
    if (mpi_nprocs>1 && mpi_rank>0) {
      MPI_Recv(&buffer,1,MPI_INT,mpi_rank-1,tag,MPI_COMM_WORLD,
	       MPI_STATUS_IGNORE);
    }
#endif
  
    // Read the data file
    std::string fname=sv[1];
    scr_out << "Opening file " << fname << " for hastings." << endl;
    hdf_file hf;
    hf.open(fname);
    table_units<> file_tab;
    hdf_input(hf,file_tab,"markov_chain0");
    hf.close();
    scr_out << "Done opening file " << fname << " for hastings." << endl;

#ifdef O2SCL_MPI
    if (mpi_nprocs>1 && mpi_rank<mpi_nprocs-1) {
      MPI_Send(&buffer,1,MPI_INT,mpi_rank+1,tag,MPI_COMM_WORLD);
    }
#endif

    // Create a new column equal to mult times weight
    file_tab.function_column("mult*weight","mwgt");
  
    // Remove
    double max_mwgt=file_tab.max("mwgt");
    if (debug) scr_out << "lines: " << file_tab.get_nlines() << endl;
    file_tab.add_constant("max_mwgt",max_mwgt);
    file_tab.delete_rows("mwgt<0.1*max_mwgt");
    if (debug) scr_out << "lines: " << file_tab.get_nlines() << endl;
  
    // The total number of variables
    size_t nv=n_params+n_sources;
    if (debug) {
      scr_out << n_params << " parameters and " << n_sources << " sources."
	      << endl;
    }
    hg_best.resize(nv);
  
    // Find the average values
    for(size_t i=0;i<n_params;i++) {
      string str_i=((string)"param_")+data_arr[0].modp->param_name(i);
      hg_best[i]=wvector_mean(file_tab.get_nlines(),file_tab[str_i],
			      file_tab["mwgt"]);
    }
    for(size_t i=0;i<n_sources;i++) {
      string str_i=((string)"Mns_")+source_names[i];
      hg_best[i+n_params]=wvector_mean(file_tab.get_nlines(),file_tab[str_i],
				       file_tab["mwgt"]);
    }
  
    // Construct the covariance matrix
    ubmatrix covar(nv,nv);
    for(size_t i=0;i<n_params;i++) {
      string str_i=((string)"param_")+data_arr[0].modp->param_name(i);
      for(size_t j=i;j<n_params;j++) {
	string str_j=((string)"param_")+data_arr[0].modp->param_name(j);
	covar(i,j)=wvector_covariance(file_tab.get_nlines(),
				      file_tab[str_i],file_tab[str_j],
				      file_tab["mult"]);
	if (debug) {
	  scr_out << "Covar: " << i << " " << j << " "
		  << covar(i,j) << endl;
	}
	covar(j,i)=covar(i,j);
      }
      for(size_t j=0;j<n_sources;j++) {
	string str_j=((string)"Mns_")+source_names[j];
	covar(i,j+n_params)=wvector_covariance(file_tab.get_nlines(),
					       file_tab[str_i],file_tab[str_j],
					       file_tab["mult"]);
	if (debug) {
	  scr_out << "Covar: " << i << " " << j+n_params << " "
		  << covar(i,j+n_params) << endl;
	}
	covar(j+n_params,i)=covar(i,j+n_params);
      }
    }
    for(size_t i=0;i<n_sources;i++) {
      string str_i=((string)"Mns_")+source_names[i];
      for(size_t j=i;j<n_sources;j++) {
	string str_j=((string)"Mns_")+source_names[j];
	covar(i+n_params,j+n_params)=
	  wvector_covariance(file_tab.get_nlines(),
			     file_tab[str_i],file_tab[str_j],
			     file_tab["mult"]);
	if (debug) {
	  scr_out << "Covar: " << i+n_params << " " << j+n_params << " "
		  << covar(i+n_params,j+n_params) << endl;
	}
	covar(j+n_params,i+n_params)=covar(i+n_params,j+n_params);
      }
    }

    // Perform the Cholesky decomposition
    hg_chol=covar;
    o2scl_linalg::cholesky_decomp(nv,hg_chol);

    // Find the inverse
    hg_covar_inv=hg_chol;
    o2scl_linalg::cholesky_invert<ubmatrix>(nv,hg_covar_inv);
  
    // Force hg_chol to be lower triangular
    for(size_t i=0;i<nv;i++) {
      for(size_t j=0;j<nv;j++) {
	if (i<j) hg_chol(i,j)=0.0;
      }
    }

    // Compute the normalization, weighted by the likelihood function
    hg_norm=1.0;
    size_t step=file_tab.get_nlines()/20;
    if (step<1) step=1;
    double renorm=0.0;
    double wgt_sum=0.0;
    for(size_t i=0;i<file_tab.get_nlines();i+=step) {
      ubvector e(n_params,n_sources);
      for(size_t j=0;j<n_params;j++) {
	string str_j=((string)"param_")+data_arr[0].modp->param_name(j);
	e.params[j]=file_tab.get(str_j,i);
      }
      for(size_t j=0;j<n_sources;j++) {
	string str_j=((string)"Mns_")+source_names[j];
	e.mass[j]=file_tab.get(str_j,i);
      }
      double wgt=file_tab.get("mult",i)*file_tab.get("weight",i);
      double rat=wgt/approx_like(e);
      renorm+=wgt*wgt/approx_like(e);
      if (debug) {
	scr_out << wgt << " " << approx_like(e) << " " << rat << endl;
      }
      wgt_sum+=wgt;
    }
    renorm/=((double)wgt_sum);
    hg_norm*=renorm;
    if (debug) {
      scr_out << "New normalization: " << hg_norm << endl;
    }

    step=file_tab.get_nlines()/20;
    if (step<1) step=1;
    for(size_t i=0;i<file_tab.get_nlines();i+=step) {
      ubvector e(n_params,n_sources);
      for(size_t j=0;j<n_params;j++) {
	string str_j=((string)"param_")+data_arr[0].modp->param_name(j);
	e.params[j]=file_tab.get(str_j,i);
      }
      for(size_t j=0;j<n_sources;j++) {
	string str_j=((string)"Mns_")+source_names[j];
	e.mass[j]=file_tab.get(str_j,i);
      }
      double wgt=file_tab.get("mult",i)*file_tab.get("weight",i);
      double rat=wgt/approx_like(e);
      if (debug) {
	scr_out << wgt << " " << approx_like(e) << " " << rat << endl;
      }
    }
    hg_mode=1;

    return 0;
  }

#endif
  
  /** \brief Set the first point
   */
  int set_initial_point(std::vector<std::string> &sv, 
			bool itive_com) {

    if (sv.size()<2) {
      std::cout << "No arguments given to 'initial-point'." << std::endl;
      return o2scl::exc_efailed;
    }

    if (sv[1]==((std::string)"values")) {

      initial_point.resize(sv.size()-1);
      for(size_t i=2;i<sv.size();i++) {
	initial_point[i-2]=o2scl::function_to_double(sv[i]);
      }
      initial_point_type=fp_unspecified;

    } else if (sv[1]==((std::string)"prefix")) {
  
      initial_point_type=fp_last;
      initial_point_file=sv[2]+((std::string)"_")+
      o2scl::itos(this->mpi_rank)+"_out";
      
    } else if (sv[1]==((std::string)"last")) {
      initial_point_type=fp_last;
      initial_point_file=sv[2];
    } else if (sv[1]==((std::string)"best")) {
      initial_point_type=fp_best;
      initial_point_file=sv[2];
    } else if (sv[1]==((std::string)"N")) {
      initial_point_type=o2scl::stoi(sv[2]);
      initial_point_file=sv[3];
    }

    return 0;
  }

  /** \brief Create an MCMC object with model \c m
   */
  mcmc_mpi() {

    // Initial values for MPI paramers
    mpi_nprocs=1;
    mpi_rank=0;
    mpi_start_time=0.0;

    // Parameters
    prefix="mcmc";
    output_meas=true;

    // Default to 24 hours
    max_time=3.6e3*24;
    max_iters=0;

    // True if scr_out has been opened
    file_opened=false;
    first_file_update=false;
    
    initial_point_file="";
    initial_point_type=fp_unspecified;

    file_update_iters=40;
    max_chain_size=10000;

    chain_index=0;
    
    // ---------------------------------------
    // Set error handler for this thread
    
    o2scl::err_hnd=&this->error_handler;
      
  }

  };

  /** \brief MCMC class with a command-line interface
   */
  template<class func_t, class fill_t, class data_t,
    class vec_t=ubvector> class mcmc_cli :
    public o2scl::mcmc_mpi<func_t,fill_t,data_t,vec_t> {

  protected:
    
  typedef o2scl::mcmc_mpi<func_t,fill_t,data_t,vec_t> parent_t;

#ifdef O2SCL_READLINE
  /// Command-line interface
  o2scl::cli_readline cl;
#else
  /// Command-line interface
  o2scl::cli cl;
#endif

  /** \brief The arguments sent to the command-line
   */
  std::vector<std::string> cl_args;

  /// \name Parameter objects for the 'set' command
  //@{
  o2scl::cli::parameter_double p_step_fac;
  o2scl::cli::parameter_size_t p_n_warm_up;
  o2scl::cli::parameter_int p_user_seed;
  o2scl::cli::parameter_size_t p_max_bad_steps;
  o2scl::cli::parameter_size_t p_n_walk;
  o2scl::cli::parameter_bool p_aff_inv;
  o2scl::cli::parameter_double p_max_time;
  o2scl::cli::parameter_size_t p_max_iters;
  o2scl::cli::parameter_int p_max_chain_size;
  o2scl::cli::parameter_int p_file_update_iters;
  o2scl::cli::parameter_bool p_output_meas;
  o2scl::cli::parameter_string p_prefix;
  o2scl::cli::parameter_int p_verbose;
  //@}
  
  /** \brief Initial write to HDF5 file 
   */
  virtual void first_update(o2scl_hdf::hdf_file &hf) {
    parent_t::first_update(hf);
    hf.sets_vec("cl_args",this->cl_args);
    return;
  }
  
  public:
  
  /// Main wrapper for parsing command-line arguments
  virtual void run(int argc, char *argv[]) {
      
    // ---------------------------------------
    // Process command-line arguments and run
      
    setup_cli();
      
#ifdef O2SCL_MPI
    // Get MPI rank, etc.
    MPI_Comm_rank(MPI_COMM_WORLD,&this->mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&this->mpi_nprocs);
#endif
      
    // Process arguments
    for(int i=0;i<argc;i++) {
      this->cl_args.push_back(argv[i]);
    }
      
    this->cl.prompt="mcmc> ";
    this->cl.run_auto(argc,argv);
      
    if (this->file_opened) {
      //Close main output file
      this->scr_out.close();
    }
      
    return;
  }    
    
  /// \name Customization functions
  //@{
  /** \brief Set up the 'cli' object
      
      This function just adds the four commands and the 'set' parameters
  */
  virtual void setup_cli() {

    // ---------------------------------------
    // Set commands/options

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
       "options.",new o2scl::comm_option_mfptr<mcmc_cli>
       (this,&mcmc_cli::set_initial_point),
       o2scl::cli::comm_option_both}
      /*
	{'s',"hastings","Specify distribution for M-H step",
	1,1,"<filename>",
	((string)"Desc. ")+"Desc2.",
	new comm_option_mfptr<mcmc_mpi>(this,&mcmc_mpi::hastings),
	cli::comm_option_both}
      */
    };
    this->cl.set_comm_option_vec(nopt,options);

    p_file_update_iters.i=&this->file_update_iters;
    p_file_update_iters.help=((std::string)"Number of MCMC successes ")+
    "between file upates (default 40, minimum value 1).";
    this->cl.par_list.insert(std::make_pair("file_update_iters",
					    &p_file_update_iters));
    
    p_max_chain_size.i=&this->max_chain_size;
    p_max_chain_size.help=((std::string)"Maximum Markov chain size ")+
    "(default 10000).";
    this->cl.par_list.insert(std::make_pair("max_chain_size",
					    &p_max_chain_size));
    
    p_max_time.d=&this->max_time;
    p_max_time.help=((std::string)"Maximum run time in seconds ")+
    "(default 86400 sec or 1 day).";
    this->cl.par_list.insert(std::make_pair("max_time",&p_max_time));
    
    p_max_iters.s=&this->max_iters;
    p_max_iters.help=((std::string)"If non-zero, limit the number of ")+
    "iterations to be less than the specified number (default zero).";
    this->cl.par_list.insert(std::make_pair("max_iters",&p_max_iters));
    
    p_prefix.str=&this->prefix;
    p_prefix.help="Output file prefix (default 'mcmc\').";
    this->cl.par_list.insert(std::make_pair("prefix",&p_prefix));
    
    p_output_meas.b=&this->output_meas;
    p_output_meas.help=((std::string)"If true, output next point ")+
    "to the '_scr' file before calling TOV solver (default true).";
    this->cl.par_list.insert(std::make_pair("output_meas",&p_output_meas));
    
    p_step_fac.d=&this->step_fac;
    p_step_fac.help=((std::string)"MCMC step factor. The step size for ")+
    "each variable is taken as the difference between the high and low "+
    "limits divided by this factor (default 10.0). This factor can "+
    "be increased if the acceptance rate is too small, but care must "+
    "be taken, e.g. if the conditional probability is multimodal. If "+
    "this step size is smaller than 1.0, it is reset to 1.0 .";
    this->cl.par_list.insert(std::make_pair("step_fac",&p_step_fac));

    p_n_warm_up.s=&this->n_warm_up;
    p_n_warm_up.help=((std::string)"Minimum number of warm up iterations ")+
    "(default 0).";
    this->cl.par_list.insert(std::make_pair("n_warm_up",&p_n_warm_up));

    p_verbose.i=&this->verbose;
    p_verbose.help=((std::string)"Verbosity parameter ")+
    "(default 0).";
    this->cl.par_list.insert(std::make_pair("verbose",&p_verbose));

    p_max_bad_steps.s=&this->max_bad_steps;
    p_max_bad_steps.help=((std::string)"Maximum number of bad steps ")+
    "(default 1000).";
    this->cl.par_list.insert(std::make_pair("max_bad_steps",&p_max_bad_steps));

    p_n_walk.s=&this->n_walk;
    p_n_walk.help=((std::string)"Number of walkers ")+
    "(default 1).";
    this->cl.par_list.insert(std::make_pair("n_walk",&p_n_walk));

    p_user_seed.i=&this->user_seed;
    p_user_seed.help=((std::string)"Seed for multiplier for random ")+
    "number generator. If zero is given (the default), then mcmc() "+
    "uses time(0) to generate a random seed.";
    this->cl.par_list.insert(std::make_pair("user_seed",&p_user_seed));
    
    p_aff_inv.b=&this->aff_inv;
    p_aff_inv.help=((std::string)"If true, then use affine-invariant ")+
    "sampling (default false).";
    this->cl.par_list.insert(std::make_pair("aff_inv",&p_aff_inv));
    
    return;
  }
  
  };
  
  // End of bamr namespace
}

#endif
