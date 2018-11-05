/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017-2018, Andrew W. Steiner
  
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
#ifndef O2SCL_INTERPM_KRIGE_H
#define O2SCL_INTERPM_KRIGE_H

/** \file interpm_krige.h
    \brief File defining \ref o2scl::interpm_krige and 
    \ref o2scl::interpm_krige_nn
*/

#include <iostream>
#include <string>
#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <gsl/gsl_combination.h>

#include <o2scl/err_hnd.h>
#include <o2scl/vector.h>
#include <o2scl/vec_stats.h>
#include <o2scl/linear_solver.h>
#include <o2scl/columnify.h>
#include <o2scl/cholesky.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multi-dimensional interpolation by kriging

      \note The set data functions for this class uses a particular
      format, one different format than that in \ref
      o2scl::interpm_idw . This design choice makes it easier to pass
      vector arguments to the covariance function and the linear
      algebra routines. The x and y objects should be of the form
      <tt>x[n_points][n_in]</tt> and <tt>y[n_out][n_points]</tt>.

      \note This class assumes that the function specified in the
      call to set_data() is the same as that passed to the
      eval() functions. If this is not the case, the
      behavior of this class is undefined.

      \note Experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_col_t=boost::numeric::ublas::matrix_column<
    boost::numeric::ublas::matrix<double> >,
    class covar_func_t=std::function<double(const vec_t &,const vec_t &)> >
    class interpm_krige {    
    
  public:
    
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
  typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;
    
  protected:

  /** \brief Inverse covariance matrix times function vector
   */
  std::vector<ubvector> Kinvf;
  
  public:

  interpm_krige() {
    data_set=false;
    verbose=0;
    matrix_mode=matrix_cholesky;
  }

  /// \name Select matrix inversion method
  //@{
  /** \brief Method for matrix inversion 
   */
  size_t matrix_mode;
  /// Use Cholesky decomposition
  static const size_t matrix_cholesky=0;
  /// Use LU decomposition
  static const size_t matrix_LU=1;
  //@}
    
  /** \brief Verbosity parameter (default 0)
   */
  int verbose;

  /** \brief Initialize the data for the interpolation

      \note This function works differently than 
      \ref o2scl::interpm_idw::set_data() . See this
      class description for more details.
   */
  template<class vec_vec_t, class vec_vec2_t>
  int set_data_noise(size_t n_in, size_t n_out, size_t n_points,
		     vec_vec_t &x, vec_vec2_t &y, 
		     std::vector<covar_func_t> &fcovar,
		     const vec_t &noise_var, bool rescale=false,
		     bool err_on_fail=true) {
    
    if (n_points<2) {
      O2SCL_ERR2("Must provide at least two points in ",
		 "interpm_krige::set_data_noise()",exc_efailed);
    }
    if (n_in<1) {
      O2SCL_ERR2("Must provide at least one input column in ",
		 "interpm_krige::set_data_noise()",exc_efailed);
    }
    if (n_out<1) {
      O2SCL_ERR2("Must provide at least one output column in ",
		 "interpm_krige::set_data_noise()",exc_efailed);
    }
    if (noise_var.size()<1) {
      O2SCL_ERR2("Noise vector empty in ",
		 "interpm_krige::set_data_noise()",exc_efailed);
    }
    np=n_points;
    nd_in=n_in;
    nd_out=n_out;
    ptrs_x.resize(n_points);
    for(size_t i=0;i<n_points;i++) {
      if (x[i].size()!=n_in) {
	O2SCL_ERR2("Size of x not correct in ",
		   "interpm_krige::set_data_noise().",o2scl::exc_efailed);
      }
      std::swap(ptrs_x[i],x[i]);
    }
    // We don't need to copy the 'y' data, but we double check
    // that it is properly sized
    for(size_t i=0;i<n_out;i++) {
      if (y[i].size()!=n_points) {
	O2SCL_ERR2("Size of y not correct in ",
		   "interpm_krige::set_data_noise().",o2scl::exc_efailed);
      }
    }
    data_set=true;
    
    if (verbose>0) {
      std::cout << "interpm_krige::set_data_noise() : Using " << n_points
		<< " points with " << nd_in << " input variables and\n\t"
		<< nd_out << " output variables." << std::endl;
    }
    
    if (rescale==true) {
      min.resize(n_in);
      max.resize(n_in);
      for(size_t j=0;j<n_in;j++) {
	min[j]=ptrs_x[0][j];
	max[j]=ptrs_x[0][j];
	for(size_t i=1;i<n_points;i++) {
	  double val=ptrs_x[i][j];
	  if (val>max[j]) max[j]=val;
	  if (val<min[j]) min[j]=val;
	}
	for(size_t i=1;i<n_points;i++) {
	  ptrs_x[i][j]=(ptrs_x[i][j]-min[j])/(max[j]-min[j])*2.0-0.5;
	}
      }
      if (verbose>1) {
	std::cout << "interpm_krige::set_data_noise() rescaled."
		  << std::endl;
      }
    }

    Kinvf.resize(n_out);

    // Loop over all output functions
    for(size_t iout=0;iout<n_out;iout++) {
	
      // Construct the KXX matrix
      ubmatrix KXX(n_points,n_points);
      for(size_t irow=0;irow<n_points;irow++) {
	for(size_t icol=0;icol<n_points;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else if (irow==icol) {
	    KXX(irow,icol)=fcovar[iout](ptrs_x[irow],ptrs_x[icol])+
	      noise_var[iout % noise_var.size()];
	  } else {
	    KXX(irow,icol)=fcovar[iout](ptrs_x[irow],ptrs_x[icol]);
	  }
	}
      }

      if (matrix_mode==matrix_LU) {
	
	// Construct the inverse of KXX
	ubmatrix inv_KXX(n_points,n_points);
	o2scl::permutation p(n_points);
	int signum;
	o2scl_linalg::LU_decomp(n_points,KXX,p,signum);
	if (o2scl_linalg::diagonal_has_zero(n_points,KXX)) {
	  if (err_on_fail) {
	    O2SCL_ERR2("Matrix singular (LU method) ",
		       "in interpm_krige::set_data_noise().",
		       o2scl::exc_esing);
	  }
	  return 1;
	}
	o2scl_linalg::LU_invert<ubmatrix,ubmatrix,ubmatrix_column>
	  (n_points,KXX,p,inv_KXX);
	
	// Inverse covariance matrix times function vector
	Kinvf[iout].resize(n_points);
	boost::numeric::ublas::axpy_prod(inv_KXX,y[iout],Kinvf[iout],true);
	
      } else {
	
	// Construct the inverse of KXX
	int cret=o2scl_linalg::cholesky_decomp(n_points,KXX,false);
	if (cret!=0) {
	  if (err_on_fail) {
	    O2SCL_ERR2("Matrix singular (Cholesky method) ",
		       "in interpm_krige::set_data_noise().",
		       o2scl::exc_esing);
	  }
	  return 2;
	}
	ubmatrix &inv_KXX=KXX;
	o2scl_linalg::cholesky_invert<ubmatrix>(n_points,inv_KXX);
	
	// Inverse covariance matrix times function vector
	Kinvf[iout].resize(n_points);
	boost::numeric::ublas::axpy_prod(inv_KXX,y[iout],Kinvf[iout],true);
	
      }
      
      if (verbose>1) {
	std::cout << "interpm_krige::set_data_noise() finished " << iout
		  << " of " << n_out << "." << std::endl;
      }
	
      if (verbose>1) {
	std::cout << "interpm_krige::set_data_noise() done."
		  << std::endl;
      }
    }
      
    return 0;
  }

  /** \brief Given covariance function \c fcovar and input vector \c x
      store the result of the interpolation in \c y
  */
  template<class vec2_t, class vec3_t>
  void eval(const vec2_t &x, vec3_t &y,
	    std::vector<covar_func_t> &fcovar) const {
    
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_krige::eval().",
		exc_einval);
    }

    if (min.size()>0) {

      // If necessary, rescale before evaluating the interpolated
      // result
      vec_t x2(nd_in);
      for(size_t iin=0;iin<nd_in;iin++) {
	x2[iin]=(x[iin]-min[iin])/(max[iin]-min[iin])*2.0-0.5;
      }

      // Evaluate the interpolated result
      for(size_t iout=0;iout<nd_out;iout++) {
	y[iout]=0.0;
	for(size_t ipoints=0;ipoints<np;ipoints++) {
	  y[iout]+=fcovar[iout](x2,ptrs_x[ipoints])*Kinvf[iout][ipoints];
	}
      }

    } else {
      
      // Evaluate the interpolated result
      for(size_t iout=0;iout<nd_out;iout++) {
	y[iout]=0.0;
	for(size_t ipoints=0;ipoints<np;ipoints++) {
	  y[iout]+=fcovar[iout](x,ptrs_x[ipoints])*Kinvf[iout][ipoints];
	}
      }

    }

    return;
      
  }
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
  /// The number of points
  size_t np;
  /// The number of dimensions of the inputs
  size_t nd_in;
  /// The number of dimensions of the outputs
  size_t nd_out;
  /// The user-specified data
  std::vector<vec_t> ptrs_x;
  /// The user-specified data
  std::vector<vec_t> ptrs_y;
  /// True if the data has been specified
  bool data_set;
  /// Minimum values for rescaling
  ubvector min;
  /// Maximum values for rescaling
  ubvector max;
  
#endif
    
    };

#ifdef O2SCL_NEVER_DEFINED  

  /** \brief One-dimensional interpolation using an 
      optimized covariance function
      
      See also the \ref intp_section section of the \o2 User's guide. 

      \note This class is experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_col_t=boost::numeric::ublas::matrix_column<
    boost::numeric::ublas::matrix<double> >,
    class interpm_krige_optim :
    public interpm_krige<vec_t,mat_col_t> {    
    
  public:

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;

  protected:

  /// Function object for the covariance
  std::function<double(double,double)> ff;
  
  /// The covariance function length scale
  double len;

  /// The quality factor of the optimization
  double qual;

  /// The covariance function
  double covar(double x1, double x2) {
    return exp(-(x1-x2)*(x1-x2)/len/len/2.0);
  }

  /// If true, data was rescaled
  bool rescaled;
  
  /// Pointer to the user-specified minimizer
  min_base<> *mp;
  
  /** \brief Function to optimize the covariance parameters
   */
  double qual_fun(double x, double noise_var, size_t &iout, int &success) {

    len=x;
    success=0;

    size_t size=this->sz;

    if (mode==mode_loo_cv) {
      
    } else if (mode==mode_max_lml) {
      
      // Construct the KXX matrix
      ubmatrix KXX(size,size);
      for(size_t irow=0;irow<size;irow++) {
	for(size_t icol=0;icol<size;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else {
	    KXX(irow,icol)=exp(-pow((ptrs_x[iout][irow]-
				     ptrs_x[iout][icol])/len,2.0)/2.0);
	    if (irow==icol) KXX(irow,icol)+=noise_var;
	  }
	}
      }

      // Note: We have to use LU here because O2scl doesn't yet
      // have a lndet() function for Cholesky decomp
      
      // Construct the inverse of KXX
      ubmatrix inv_KXX(size,size);
      o2scl::permutation p(size);
      int signum;
      o2scl_linalg::LU_decomp(size,KXX,p,signum);
      if (o2scl_linalg::diagonal_has_zero(size,KXX)) {
	success=1;
	return 1.0e99;
      }
      o2scl_linalg::LU_invert<ubmatrix,ubmatrix,ubmatrix_column>
	(size,KXX,p,inv_KXX);
      
      double lndet=o2scl_linalg::LU_lndet<ubmatrix>(size,KXX);
      
      // Inverse covariance matrix times function vector
      this->Kinvf.resize(size);
      boost::numeric::ublas::axpy_prod(inv_KXX,*this->py,this->Kinvf,true);

      // Compute the log of the marginal likelihood, without
      // the constant term
      for(size_t i=0;i<size;i++) {
	qual+=0.5*(*this->py)[i]*this->Kinvf[i];
      }
      qual+=0.5*lndet;
    }

    return qual;
  }
  
  public:

  interpm_krige_optim() {
    nlen=20;
    full_min=false;
    mp=&def_min;
    verbose=0;
    mode=mode_loo_cv;
    loo_npts=100;
  }

  /// \name Function to minimize and various option
  //@{
  /// Leave-one-out cross validation
  static const size_t mode_loo_cv=1;
  /// Minus Log-marginal-likelihood
  static const size_t mode_max_lml=2;
  /// Function to minimize (default \ref mode_loo_cv)
  size_t mode;
  /// Number of points to test for cross validation (default 100)
  size_t loo_npts;
  ///@}
    
  /// Verbosity parameter
  int verbose;
  
  /** \brief Number of length scale points to try when full minimizer 
      is not used (default 20)
  */
  size_t nlen;

  /// Default minimizer
  min_brent_gsl<> def_min;

  /// If true, use the full minimizer
  bool full_min;

  /// Initialize interpolation routine
  template<class vec_vec_t, class vec_vec2_t>
  void set_data_noise(size_t n_in, size_t n_out, size_t n_points,
		      vec_vec_t &x, vec_vec2_t &y, 
		      const vec_t &noise_var, bool rescale=false,
		      bool err_on_fail=true) {

    if (n_points<2) {
      O2SCL_ERR2("Must provide at least two points in ",
		 "interpm_krige_optim::set_data_noise()",exc_efailed);
    }
    if (n_in<1) {
      O2SCL_ERR2("Must provide at least one input column in ",
		 "interpm_krige_optim::set_data_noise()",exc_efailed);
    }
    if (n_out<1) {
      O2SCL_ERR2("Must provide at least one output column in ",
		 "interpm_krige_optim::set_data_noise()",exc_efailed);
    }
    if (noise_var.size()<1) {
      O2SCL_ERR2("Noise vector empty in ",
		 "interpm_krige::set_data_noise()",exc_efailed);
    }
   
    // Set parent data members
    np=n_points;
    nd_in=n_in;
    nd_out=n_out;
    ptrs_x.resize(n_points);
    rescaled=rescale;
    for(size_t i=0;i<n_points;i++) {
      if (x[i].size()!=n_in) {
	O2SCL_ERR2("Size of x not correct in ",
		   "interpm_krige_optim::set_data_noise().",
		   o2scl::exc_efailed);
      }
      std::swap(ptrs_x[i],x[i]);
    }
    data_set=true;
       
    if (verbose>0) {
      std::cout << "interpm_krige_optim::set_data_noise() : Using "
		<< n_points
		<< " points with " << nd_in << " input variables and\n\t"
		<< nd_out << " output variables." << std::endl;
    }

    // Get max and min of all parameters
    min.resize(n_in);
    max.resize(n_in);
    for(size_t j=0;j<n_in;j++) {
      min[j]=ptrs_x[j][0];
      max[j]=ptrs_x[j][0];
      for(size_t i=1;i<n_points;i++) {
	double val=ptrs_x[j][i];
	if (val>max[j]) max[j]=val;
	if (val<min[j]) min[j]=val;
      }
    }
    
    if (rescale==true) {
      for(size_t j=0;j<n_in;j++) {
	for(size_t i=1;i<n_points;i++) {
	  ptrs_x[j][i]=(ptrs_x[j][i]-min[j])/(max[j]-min[j])*2.0-0.5;
	}
      }
      if (verbose>1) {
	std::cout << "interpm_krige_optim::set_data_noise() rescaled."
		  << std::endl;
      }
    }

    int success=0;

    Kinvf.resize(n_out);

    if (full_min) {
     
      if (verbose>1) {
	std::cout << "interp_krige_optim: full minimization"
		  << std::endl;
      }
     
      // Loop over all output functions
      for(size_t iout=0;iout<n_out;iout++) {
	
	double len_opt;
	// Choose average distance for first guess
	if (rescaled) {
	  len_opt=(max[iout]-min[iout])/((double)np);
	} else {
	  len_opt=2.0/((double)np);
	}
	
	funct mf=std::bind
	  (std::mem_fn<double(double,double,int &)>
	   (&interp_krige_optim<vec_t,vec2_t>::qual_fun),this,
	   std::placeholders::_1,noise_var,std::ref(success));
	
	mp->min(len_opt,qual,mf);
	len=len_opt;
	
	if (success!=0) {
	  if (err_on_fail) {
	    O2SCL_ERR2("Minimization failed in ",
		       "interp_krige_optim::set_noise().",
		       o2scl::exc_efailed);
	  }
	}
	
      } else {
	
	if (verbose>1) {
	  std::cout << "interp_krige_optim: simple minimization"
		    << std::endl;
	}
	
	double len_avg;
	// Choose average distance for first guess
	if (rescaled) {
	  len_avg=(max[iout]-min[iout])/((double)np);
	} else {
	  len_avg=2.0/((double)np);
	}

	double len_min=len_avg/100.0;
	double len_max=len_avg*100.0;
	double len_ratio=len_max/len_min;
	
	if (verbose>1) {
	  std::cout << "len (min,max,ratio): " << len_min << " "
		    << len_max << " "
		    << pow(len_ratio,((double)1)/((double)nlen-1))
		    << std::endl;
	}
	
	// Initialize to zero to prevent uninit'ed var. warnings
	double min_qual=0.0, len_opt=0.0;
	
	if (verbose>1) {
	  std::cout << "ilen len qual fail min_qual" << std::endl;
	}
	
	// Loop over the full range, finding the optimum
	bool min_set=false;
	for(size_t j=0;j<nlen;j++) {
	  len=len_min*pow(len_ratio,((double)j)/((double)nlen-1));
	  
	  int success=0;
	  qual=qual_fun(len,noise_var,success);
	
	if (success==0 && (min_set==false || qual<min_qual)) {
	  len_opt=len;
	  min_qual=qual;
	  min_set=true;
	}
	
	if (verbose>1) {
	  std::cout << "interp_krige_optim: ";
	  std::cout.width(2);
	  std::cout << j << " " << len << " " << qual << " "
		    << success << " " << min_qual << " "
		    << len_opt << std::endl;
	}
	  
      }
      
      // Now that we've optimized the covariance function,
      // just use the parent class to interpolate
      len=len_opt;

    }

    ff=std::bind(std::mem_fn<double(double,double)>
		 (&interp_krige_optim<vec_t,vec2_t>::covar),this,
		 std::placeholders::_1,std::placeholders::_2);
    
    this->set_covar_noise(size,x,y,ff,noise_var);
      
    return 0;
  }

  /// Initialize interpolation routine
  virtual void set(size_t size, const vec_t &x, const vec2_t &y) {

    // Use the mean absolute value to determine noise
    double mean_abs=0.0;
    for(size_t j=0;j<size;j++) {
      mean_abs+=fabs(y[j]);
    }
    mean_abs/=size;

    set_noise(size,x,y,mean_abs/1.0e8,true);
    
    return;
  }
  
  
  };

#endif  

  /** \brief Multi-dimensional interpolation by kriging with 
      nearest-neighbor 

      \note This class assumes that the function specified in the
      call to set_data() is the same as that passed to the
      eval() functions. If this is not the case, the
      behavior of this class is undefined.

      \note Experimental.
  */
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_col_t=boost::numeric::ublas::matrix_column<
    boost::numeric::ublas::matrix<double> >,
    class covar_func_t=std::function<double(const vec_t &,const vec_t &)> >
    class interpm_krige_nn {
    
  public:
    
  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
  
  interpm_krige_nn() {
    data_set=false;
    verbose=0;
  }

  /** \brief Verbosity parameter (default 0)
   */
  int verbose;

  /** \brief Initialize the data for the interpolation
   */
  template<class vec_vec_t, class vec_vec2_t>
  void set_data(size_t n_in, size_t n_out, size_t n_points,
		vec_vec_t &x, vec_vec2_t &y, 
		std::vector<covar_func_t> &fcovar, size_t order) {

    if (n_points<2) {
      O2SCL_ERR2("Must provide at least two points in ",
		 "interpm_krige_nn::set_data()",exc_efailed);
    }
    if (n_in<1) {
      O2SCL_ERR2("Must provide at least one input column in ",
		 "interpm_krige_nn::set_data()",exc_efailed);
    }
    if (n_out<1) {
      O2SCL_ERR2("Must provide at least one output column in ",
		 "interpm_krige_nn::set_data()",exc_efailed);
    }
    np=n_points;
    nd_in=n_in;
    nd_out=n_out;
    ptrs_x.resize(n_points);
    n_order=order;
    for(size_t i=0;i<n_points;i++) {
      if (x[i].size()!=n_in) {
	O2SCL_ERR2("Size of x not correct in ",
		   "interpm_krige_nn::set_data().",o2scl::exc_efailed);
      }
      std::swap(ptrs_x[i],x[i]);
    }
    ptrs_y.resize(n_out);
    for(size_t i=0;i<n_out;i++) {
      if (y[i].size()!=n_points) {
	O2SCL_ERR2("Size of y not correct in ",
		   "interpm_krige_nn::set_data().",o2scl::exc_efailed);
      }
      std::swap(ptrs_y[i],y[i]);
    }
    data_set=true;
      
    if (verbose>0) {
      std::cout << "interpm_krige_nn::set_data() : Using " << n_points
		<< " points with " << nd_in << " input variables and\n\t"
		<< nd_out << " output variables and order "
		<< n_order << " ." << std::endl;
    }

    return;
  }

  /** \brief Given covariance function \c fcovar and input vector \c x
      store the result of the interpolation in \c y
  */
  template<class vec2_t, class vec3_t>
  void eval(const vec2_t &x, vec3_t &y,
	    std::vector<covar_func_t> &fcovar) const {
      
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_krige_nn::eval().",
		exc_einval);
    }
      
    y.resize(nd_out);
    
    // Loop over all output functions
    for(size_t iout=0;iout<nd_out;iout++) {

      // Find points closest to requested point, as defined
      // by the negative covariance for this output function
      ubvector dists(np);
      for(size_t ip=0;ip<np;ip++) {
	dists[ip]=-fcovar[iout](x,ptrs_x[ip]);
      }
      
      // Empty index vector (resized by the vector_smallest_index
      // function)
      ubvector_size_t index;
      o2scl::vector_smallest_index<ubvector,double,ubvector_size_t>
	(np,dists,n_order,index);
      
      // Construct subset of function values for nearest neighbors
      ubvector func(n_order);
      for(size_t io=0;io<n_order;io++) {
	func[io]=ptrs_y[iout][index[io]];
      }
      
      // Construct the nearest neighbor KXX matrix
      ubmatrix KXX(n_order,n_order);
      for(size_t irow=0;irow<n_order;irow++) {
	for(size_t icol=0;icol<n_order;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else {
	    KXX(irow,icol)=fcovar[iout](ptrs_x[index[irow]],
					ptrs_x[index[icol]]);
	  }
	}
      }
      
      // Construct the inverse of KXX
      o2scl_linalg::cholesky_decomp(n_order,KXX);
      ubmatrix &inv_KXX=KXX;
      o2scl_linalg::cholesky_invert<ubmatrix>(n_order,inv_KXX);
      
      // Inverse covariance matrix times function vector
      ubvector Kinvf(n_order);
      boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);

      // Comput the final result
      y[iout]=0.0;
      for(size_t ipoints=0;ipoints<n_order;ipoints++) {
	y[iout]+=-dists[index[ipoints]]*Kinvf[ipoints];
      }
      
    }
    
    return;
    
  }
  
  /** \brief Find a set of linearly independent points 

      Given a point \c x, a covariance function 
      \c fcovar, the index of the output function
      \c iout, and an array specifying the closest points 
      \c index, this function produces a list of 
  */
  template<class vec2_t>
  void find_lin_indep(const vec2_t &x, size_t iout,
		      std::vector<covar_func_t> &fcovar,
		      ubvector_size_t &index,
		      ubvector_size_t &indep) const {
    
    if (x.size()<nd_in || index.size()<np || indep.size()<n_order
	|| iout>=nd_out) {
      std::cout << x.size() << " " << nd_in << std::endl;
      std::cout << index.size() << " " << np << std::endl;
      std::cout << indep.size() << " " << n_order << std::endl;
      std::cout << iout << " " << nd_out << std::endl;
      O2SCL_ERR("Vectors not of correct size in find_lin_indep().",
		o2scl::exc_einval);
    }
    
    bool done=false;
    while (done==false) {
      
      // Construct subset of function values for nearest neighbors
      ubvector func(n_order);
      for(size_t io=0;io<n_order;io++) {
	func[io]=ptrs_y[iout][index[indep[io]]];
      }
      
      // Construct the nearest neighbor KXX matrix
      ubmatrix KXX(n_order,n_order);
      for(size_t irow=0;irow<n_order;irow++) {
	for(size_t icol=0;icol<n_order;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else {
	    KXX(irow,icol)=fcovar[iout](ptrs_x[index[indep[irow]]],
					ptrs_x[index[indep[icol]]]);
	  }
	}
      }
      
      // Construct the inverse of KXX
      int cret=o2scl_linalg::cholesky_decomp(n_order,KXX);
      if (cret==0) {
	done=true;
      } else {
	if (verbose>1) {
	  std::cout << "Finding new independent rows." << std::endl;
	  for(size_t j=0;j<n_order;j++) {
	    std::cout << indep[j] << " " << KXX(j,j) << std::endl;
	  }
	}
	size_t max=o2scl::vector_max_value<ubvector_size_t,
	  double>(indep);
	if (verbose>1) {
	  std::cout << "Max is: " << max << std::endl;
	}
	for(size_t j=0;j<n_order;j++) {
	  if (KXX(j,j)==0.0) {
	    if (max+1<np) {
	      if (verbose>1) {
		std::cout << "Entry " << j << " is zero so replacing "
			  << "entry with " << max+1 << std::endl;
	      }
	      indep[j]=max+1;
	      max++;
	    } else {
	      O2SCL_ERR3("Failed to find set of independent points ",
			 "in interpm_krige_nn::find_lin_indep",
			 "(const vec2_t &, size_t).",
			 o2scl::exc_efailed);
	    }
	  }
	}
      }
      
    }
    return;
  }
  
  /** \brief Given covariance function \c fcovar and input vector \c x
      return the result of the interpolation for function with 
      index \c iout
  */
  template<class vec2_t>
  double eval(const vec2_t &x, size_t iout,
	      std::vector<covar_func_t> &fcovar) const {
      
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_krige_nn::eval().",
		exc_einval);
    }

    double ret;
    
    // Find points closest to requested point, as defined
    // by the negative covariance for this output function
    ubvector dists(np);
    for(size_t ip=0;ip<np;ip++) {
      dists[ip]=-fcovar[iout](x,ptrs_x[ip]);
    }
      
    ubvector_size_t index(np);
    o2scl::vector_sort_index<ubvector,ubvector_size_t>(np,dists,index);

    // Vector for storing the indexes in the index array which
    // will store the closest n_order points which are
    // linearly independent
    ubvector_size_t indep(n_order);
    for(size_t io=0;io<n_order;io++) {
      indep[io]=io;
    }

    find_lin_indep(x,iout,fcovar,index,indep);
    
    // Construct subset of function values for nearest neighbors
    ubvector func(n_order);
    for(size_t io=0;io<n_order;io++) {
      func[io]=ptrs_y[iout][index[indep[io]]];
    }
      
    // Construct the nearest neighbor KXX matrix
    ubmatrix KXX(n_order,n_order);
    for(size_t irow=0;irow<n_order;irow++) {
      for(size_t icol=0;icol<n_order;icol++) {
	if (irow>icol) {
	  KXX(irow,icol)=KXX(icol,irow);
	} else {
	  KXX(irow,icol)=fcovar[iout](ptrs_x[index[indep[irow]]],
				      ptrs_x[index[indep[icol]]]);
	}
      }
    }
	
    // Construct the inverse of KXX
    o2scl_linalg::cholesky_decomp(n_order,KXX);
    ubmatrix &inv_KXX=KXX;
    o2scl_linalg::cholesky_invert<ubmatrix>(n_order,inv_KXX);
      
    // Inverse covariance matrix times function vector
    ubvector Kinvf(n_order);
    boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);

    // Comput the final result
    ret=0.0;
    for(size_t ipoints=0;ipoints<n_order;ipoints++) {
      ret+=-dists[index[indep[ipoints]]]*Kinvf[ipoints];
    }
      
    return ret;
      
  }

  /** \brief Compute a quality factor for interpolation
      using jackknife resampling
  */
  template<class vec2_t>
  double eval_jackknife(const vec2_t &x, size_t iout,
			std::vector<covar_func_t> &fcovar) const {
    
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_krige_nn::eval_jackknife().",
		exc_einval);
    }

    double qual=0.0;
    
    // Interpolated function value inside jackknife loop
    double ytmp;
    
    // For a distance measurement, just use the the negative
    // covariance for this output function
    ubvector dists(np);
    for(size_t ip=0;ip<np;ip++) {
      dists[ip]=-fcovar[iout](x,ptrs_x[ip]);
    }
    
    // Create an index array which sorts by distance
    ubvector_size_t index(np);
    o2scl::vector_sort_index<ubvector,ubvector_size_t>(np,dists,index);
    
    // Vector for storing the indexes in the index array which
    // will store the closest n_order points which are
    // linearly independent
    ubvector_size_t indep(n_order);
    for(size_t io=0;io<n_order;io++) {
      indep[io]=io;
    }
    
    // -------------------------------------------------------------
    // Before the jackknife loop, we want to create a full
    // set of n_order linearly independent points

    find_lin_indep(x,iout,fcovar,index,indep);

    // -------------------------------------------------------------
    // Now, the jackknife loop, removing one point at a time
    
    for(size_t jk=0;jk<n_order;jk++) {

      if (verbose>0) {
	std::cout << "Jackknife: " << jk << " matching function value "
		  << ptrs_y[iout][index[jk]] << std::endl;
      }
      
      ubvector_size_t indep_jk;
      vector_copy_jackknife(indep,jk,indep_jk);
      
      // Construct subset of function values for nearest neighbors
      ubvector func(n_order-1);
      for(size_t io=0;io<n_order-1;io++) {
	func[io]=ptrs_y[iout][index[indep_jk[io]]];
      }
      
      // Construct the nearest neighbor KXX matrix
      ubmatrix KXX(n_order-1,n_order-1);
      for(size_t irow=0;irow<n_order-1;irow++) {
	for(size_t icol=0;icol<n_order-1;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else {
	    KXX(irow,icol)=fcovar[iout](ptrs_x[index[indep_jk[irow]]],
					ptrs_x[index[indep_jk[icol]]]);
	  }
	}
      }
	  
      // Construct the inverse of KXX
      o2scl_linalg::cholesky_decomp(n_order-1,KXX);
      ubmatrix &inv_KXX=KXX;
      o2scl_linalg::cholesky_invert<ubmatrix>(n_order-1,inv_KXX);
      
      // Inverse covariance matrix times function vector
      ubvector Kinvf(n_order-1);
      boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);
      
      // Comput the final result
      ytmp=0.0;
      for(size_t ipoints=0;ipoints<n_order-1;ipoints++) {
	ytmp+=-dists[index[indep_jk[ipoints]]]*Kinvf[ipoints];
      }
      
      // Add the squared deviation to y[iout]
      qual+=pow(ptrs_y[iout][index[jk]]-ytmp,2.0);

      if (verbose>0) {
	std::cout << "Original value: "
		  << ptrs_y[iout][index[jk]] << " interpolated: "
		  << ytmp << std::endl;
      }
      
      // End of jackknife loop
    }

    return qual;
  }
  
  
#ifndef DOXYGEN_INTERNAL
  
  protected:
    
  /// The order of the interpolation (specified by \ref set_data() )
  size_t n_order;
  /// The number of points
  size_t np;
  /// The number of dimensions of the inputs
  size_t nd_in;
  /// The number of dimensions of the outputs
  size_t nd_out;
  /// A vector of pointers holding the data
  std::vector<vec_t> ptrs_x;
  /// A vector of pointers holding the data
  std::vector<vec_t> ptrs_y;
  /// True if the data has been specified
  bool data_set;
    
#endif
    
    };
    
#ifndef DOXYGEN_NO_O2NS
}
#endif
    
#endif



