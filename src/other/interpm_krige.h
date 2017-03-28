/*
  -------------------------------------------------------------------
  
  Copyright (C) 2017, Andrew W. Steiner
  
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

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Multi-dimensional interpolation by kriging

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
    
  protected:

  /** \brief Inverse covariance matrix times function vector
   */
  std::vector<ubvector> Kinvf;
  
  /** \brief Pointer to user-specified covariance function array
   */
  std::vector<covar_func_t> *f;
    
  public:

  interpm_krige() {
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
		std::vector<covar_func_t> &fcovar) {

    if (n_points<3) {
      O2SCL_ERR2("Must provide at least three points in ",
		 "interpm_krige::set_data()",exc_efailed);
    }
    if (n_in<1) {
      O2SCL_ERR2("Must provide at least one input column in ",
		 "interpm_krige::set_data()",exc_efailed);
    }
    if (n_out<1) {
      O2SCL_ERR2("Must provide at least one output column in ",
		 "interpm_krige::set_data()",exc_efailed);
    }
    np=n_points;
    nd_in=n_in;
    nd_out=n_out;
    ptrs_x.resize(n_points);
    for(size_t i=0;i<n_points;i++) {
      if (x[i].size()!=n_in) {
	O2SCL_ERR2("Size of x not correct in ",
		   "interpm_krige::set_data().",o2scl::exc_efailed);
      }
      std::swap(ptrs_x[i],x[i]);
    }
    data_set=true;
    
    if (verbose>0) {
      std::cout << "interpm_krige::set_data() : Using " << n_points
		<< " points with " << nd_in << " input variables and\n\t"
		<< nd_out << " output variables." << std::endl;
    }

    Kinvf.resize(n_out);

    // Store pointer to covariance function
    f=&fcovar;

    // Loop over all output functions
    for(size_t iout=0;iout<n_out;iout++) {
	
      // Construct the KXX matrix
      ubmatrix KXX(n_points,n_points);
      for(size_t irow=0;irow<n_points;irow++) {
	for(size_t icol=0;icol<n_points;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else {
	    KXX(irow,icol)=(*f)[iout](ptrs_x[irow],ptrs_x[icol]);
	  }
	}
      }
	
      // Construct the inverse of KXX
      ubmatrix inv_KXX(n_points,n_points);
      o2scl::permutation p(n_points);
      int signum;
      if (verbose>0) {
	std::cout << "interpm_krige::set_data() : "
		  << "LU decompose and invert " << iout+1 << " of " << n_out
		  << std::endl;
      }
      o2scl_linalg::LU_decomp(n_points,KXX,p,signum);
      if (o2scl_linalg::diagonal_has_zero(n_points,KXX)) {
	O2SCL_ERR2("KXX matrix is singular in ",
		   "interpm_krige::set_data().",
		   o2scl::exc_efailed);
      }
      o2scl_linalg::LU_invert<ubmatrix,ubmatrix,mat_col_t>
	(n_points,KXX,p,inv_KXX);

      // Inverse covariance matrix times function vector
      Kinvf[iout].resize(n_points);
      boost::numeric::ublas::axpy_prod(inv_KXX,y[iout],Kinvf[iout],true);
	
    }
      
    return;
  }

  /** \brief Perform the interpolation
   */
  template<class vec2_t, class vec3_t>
  void eval(const vec2_t &x, vec3_t &y) const {
    
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_krige::eval().",
		exc_einval);
    }
    
    y.resize(nd_out);
    for(size_t iout=0;iout<nd_out;iout++) {
      y[iout]=0.0;
      for(size_t ipoints=0;ipoints<np;ipoints++) {
	y[iout]+=(*f)[iout](x,ptrs_x[ipoints])*Kinvf[iout][ipoints];
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
  /// A vector of pointers holding the data
  std::vector<vec_t> ptrs_x;
  /// True if the data has been specified
  bool data_set;
  /// Number of points to include in each interpolation (default 3)
  size_t order;
    
#endif
    
  };
    
  /** \brief Multi-dimensional interpolation by kriging with 
      nearest-neighbor 

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

    if (n_points<3) {
      O2SCL_ERR2("Must provide at least three points in ",
		 "interpm_krige::set_data()",exc_efailed);
    }
    if (n_in<1) {
      O2SCL_ERR2("Must provide at least one input column in ",
		 "interpm_krige::set_data()",exc_efailed);
    }
    if (n_out<1) {
      O2SCL_ERR2("Must provide at least one output column in ",
		 "interpm_krige::set_data()",exc_efailed);
    }
    np=n_points;
    nd_in=n_in;
    nd_out=n_out;
    ptrs_x.resize(n_points);
    norder=order;
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

    // Store pointer to covariance function array
    f=&fcovar;
      
    if (verbose>0) {
      std::cout << "interpm_krige_nn::set_data() : Using " << n_points
		<< " points with " << nd_in << " input variables and\n\t"
		<< nd_out << " output variables and order "
		<< norder << " ." << std::endl;
    }

    return;
  }

  /** \brief Perform the interpolation
   */
  template<class vec2_t, class vec3_t>
  void eval(const vec2_t &x, vec3_t &y) const {
      
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_krige::eval().",
		exc_einval);
    }
      
    y.resize(nd_out);
    
    // Loop over all output functions
    for(size_t iout=0;iout<nd_out;iout++) {
	
      // Find points closest to requested point, as defined
      // by the negative covariance for this output function
      ubvector dists(np);
      for(size_t ip=0;ip<np;ip++) {
	dists[ip]=-(*f)[iout](x,ptrs_x[ip]);
      }
      
      // Empty index vector (resized by the vector_smallest_index
      // function)
      ubvector_size_t index;
      o2scl::vector_smallest_index<ubvector,double,ubvector_size_t>
	(np,dists,norder,index);

      // Construct subset of function values for nearest neighbors
      ubvector func(norder);
      for(size_t io=0;io<norder;io++) {
	func[io]=ptrs_y[iout][index[io]];
      }
      
      // Construct the nearest neighbor KXX matrix
      ubmatrix KXX(norder,norder);
      for(size_t irow=0;irow<norder;irow++) {
	for(size_t icol=0;icol<norder;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else {
	    KXX(irow,icol)=(*f)[iout](ptrs_x[index[irow]],
				      ptrs_x[index[icol]]);
	  }
	}
      }
	
      // Construct the inverse of KXX
      ubmatrix inv_KXX(norder,norder);
      o2scl::permutation p(norder);
      int signum;
      o2scl_linalg::LU_decomp(norder,KXX,p,signum);
      if (o2scl_linalg::diagonal_has_zero(norder,KXX)) {
	O2SCL_ERR2("KXX matrix is singular in ",
		   "interpm_krige_nn::eval().",
		   o2scl::exc_efailed);
      }
      o2scl_linalg::LU_invert<ubmatrix,ubmatrix,mat_col_t>
	(norder,KXX,p,inv_KXX);
      
      // Inverse covariance matrix times function vector
      ubvector Kinvf(norder);
      boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);

      // Comput the final result
      y[iout]=0.0;
      for(size_t ipoints=0;ipoints<norder;ipoints++) {
	y[iout]+=-dists[index[ipoints]]*Kinvf[ipoints];
      }
      
    }

    return;
      
  }
  
  /** \brief Perform the interpolation
   */
  template<class vec2_t>
  double eval(const vec2_t &x, size_t iout) const {
      
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_krige::eval().",
		exc_einval);
    }

    double ret;
    
      // Find points closest to requested point, as defined
      // by the negative covariance for this output function
      ubvector dists(np);
      for(size_t ip=0;ip<np;ip++) {
	dists[ip]=-(*f)[iout](x,ptrs_x[ip]);
      }
      
      // Empty index vector (resized by the vector_smallest_index
      // function)
      ubvector_size_t index;
      o2scl::vector_smallest_index<ubvector,double,ubvector_size_t>
	(np,dists,norder,index);

      // Construct subset of function values for nearest neighbors
      ubvector func(norder);
      for(size_t io=0;io<norder;io++) {
	func[io]=ptrs_y[iout][index[io]];
      }
      
      // Construct the nearest neighbor KXX matrix
      ubmatrix KXX(norder,norder);
      for(size_t irow=0;irow<norder;irow++) {
	for(size_t icol=0;icol<norder;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else {
	    KXX(irow,icol)=(*f)[iout](ptrs_x[index[irow]],
				      ptrs_x[index[icol]]);
	  }
	}
      }
	
      // Construct the inverse of KXX
      ubmatrix inv_KXX(norder,norder);
      o2scl::permutation p(norder);
      int signum;
      o2scl_linalg::LU_decomp(norder,KXX,p,signum);
      if (o2scl_linalg::diagonal_has_zero(norder,KXX)) {
	O2SCL_ERR2("KXX matrix is singular in ",
		   "interpm_krige_nn::eval().",
		   o2scl::exc_efailed);
      }
      o2scl_linalg::LU_invert<ubmatrix,ubmatrix,mat_col_t>
	(norder,KXX,p,inv_KXX);
      
      // Inverse covariance matrix times function vector
      ubvector Kinvf(norder);
      boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);

      // Comput the final result
      ret=0.0;
      for(size_t ipoints=0;ipoints<norder;ipoints++) {
	ret+=-dists[index[ipoints]]*Kinvf[ipoints];
      }
      
    return ret;
      
  }
  
  /** \brief Use jackknife to match interpolated to real function 
      values
  */
  template<class vec2_t, class vec3_t>
  void eval_jackknife(const vec2_t &x, vec3_t &qual) const {
    
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_krige::eval().",
		exc_einval);
    }
      
    qual.resize(nd_out);
    
    // Interpolated function value inside jackknife loop
    double ytmp=0.0;
      
    // Loop over all output functions
    for(size_t iout=0;iout<nd_out;iout++) {
	
      qual[iout]=0.0;
	
      // Find points closest to requested point, as defined
      // by the negative covariance for this output function
      ubvector dists(np);
      for(size_t ip=0;ip<np;ip++) {
	dists[ip]=-(*f)[iout](x,ptrs_x[ip]);
      }
    
      // Empty index vector (resized by the vector_smallest_index
      // function)
      ubvector_size_t index;
      o2scl::vector_smallest_index<ubvector,double,ubvector_size_t>
	(np,dists,norder,index);
	
      // The jackknife loop
      for(size_t jk=0;jk<norder;jk++) {
	  
	ubvector_size_t index_jk;
	vector_copy_jackknife(index,jk,index_jk);
	  
	// Construct subset of function values for nearest neighbors
	ubvector func(norder-1);
	for(size_t io=0;io<norder-1;io++) {
	  func[io]=ptrs_y[iout][index_jk[io]];
	}
	  
	// Construct the nearest neighbor KXX matrix
	ubmatrix KXX(norder-1,norder-1);
	for(size_t irow=0;irow<norder-1;irow++) {
	  for(size_t icol=0;icol<norder-1;icol++) {
	    if (irow>icol) {
	      KXX(irow,icol)=KXX(icol,irow);
	    } else {
	      KXX(irow,icol)=(*f)[iout](ptrs_x[index_jk[irow]],
					ptrs_x[index_jk[icol]]);
	    }
	  }
	}
	  
	// Construct the inverse of KXX
	ubmatrix inv_KXX(norder-1,norder-1);
	o2scl::permutation p(norder-1);
	int signum;
	o2scl_linalg::LU_decomp(norder-1,KXX,p,signum);
	o2scl_linalg::LU_invert<ubmatrix,ubmatrix,mat_col_t>
	  (norder-1,KXX,p,inv_KXX);
	  
	// Inverse covariance matrix times function vector
	ubvector Kinvf(norder-1);
	boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);
	  
	// Comput the final result
	ytmp=0.0;
	for(size_t ipoints=0;ipoints<norder-1;ipoints++) {
	  ytmp+=-dists[index_jk[ipoints]]*Kinvf[ipoints];
	}

	// Add the squared deviation to y[iout]
	qual[iout]+=pow(ptrs_y[iout][index[jk]]-ytmp,2.0);

      }
	  
    }

    return;
      
  }
    
  /** \brief Use jackknife to match interpolated to real function 
      values
  */
  template<class vec2_t>
  double eval_jackknife(const vec2_t &x, size_t iout) const {
    
    double qual;
    
    if (data_set==false) {
      O2SCL_ERR("Data not set in interpm_krige::eval().",
		exc_einval);
    }
  
    // Interpolated function value inside jackknife loop
    double ytmp=0.0;
      
    qual=0.0;
	
    // Find points closest to requested point, as defined
    // by the negative covariance for this output function
    ubvector dists(np);
    for(size_t ip=0;ip<np;ip++) {
      dists[ip]=-(*f)[iout](x,ptrs_x[ip]);
    }
      
    // Empty index vector (resized by the vector_smallest_index
    // function)
    ubvector_size_t index;
    o2scl::vector_smallest_index<ubvector,double,ubvector_size_t>
    (np,dists,norder,index);
	
    // The jackknife loop
    for(size_t jk=0;jk<norder;jk++) {

      ubvector_size_t index_jk;
      vector_copy_jackknife(index,jk,index_jk);
	  
      // Construct subset of function values for nearest neighbors
      ubvector func(norder-1);
      for(size_t io=0;io<norder-1;io++) {
	func[io]=ptrs_y[iout][index_jk[io]];
      }
	  
      // Construct the nearest neighbor KXX matrix
      ubmatrix KXX(norder-1,norder-1);
      for(size_t irow=0;irow<norder-1;irow++) {
	for(size_t icol=0;icol<norder-1;icol++) {
	  if (irow>icol) {
	    KXX(irow,icol)=KXX(icol,irow);
	  } else {
	    KXX(irow,icol)=(*f)[iout](ptrs_x[index_jk[irow]],
				      ptrs_x[index_jk[icol]]);
	  }
	}
      }

      // Construct the inverse of KXX
      ubmatrix inv_KXX(norder-1,norder-1);
      o2scl::permutation p(norder-1);
      int signum;
      if (verbose>0) {
	std::cout << "interpm_krige_nn::eval_jackknife() "
		  << "LU decompose and invert " << jk+1 << " of "
		  << norder << std::endl;
      }
      o2scl_linalg::LU_decomp(norder-1,KXX,p,signum);
      o2scl_linalg::LU_invert<ubmatrix,ubmatrix,mat_col_t>
	(norder-1,KXX,p,inv_KXX);
	  
      // Inverse covariance matrix times function vector
      ubvector Kinvf(norder-1);
      boost::numeric::ublas::axpy_prod(inv_KXX,func,Kinvf,true);
	  
      // Comput the final result
      ytmp=0.0;
      for(size_t ipoints=0;ipoints<norder-1;ipoints++) {
	ytmp+=-dists[index_jk[ipoints]]*Kinvf[ipoints];
      }

      // Add the squared deviation to y[iout]
      qual+=pow(ptrs_y[iout][index[jk]]-ytmp,2.0);

    }
	  
      

    return qual;
      
  }
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
  
  /** \brief Pointer to user-specified covariance function array
   */
  std::vector<covar_func_t> *f;
    
  /// Desc
  size_t norder;
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
  /// Number of points to include in each interpolation (default 3)
  size_t order;
    
#endif
    
  };
    
#ifndef DOXYGEN_NO_O2NS
}
#endif
    
#endif



