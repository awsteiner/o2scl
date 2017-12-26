/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
#ifndef O2SCL_TENSOR_GRID_H
#define O2SCL_TENSOR_GRID_H

/** \file tensor_grid.h
    \brief File defining \ref o2scl::tensor_grid and rank-specific children
*/

#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_ieee_utils.h>

#include <o2scl/err_hnd.h>
#include <o2scl/interp.h>
#include <o2scl/tensor.h>
#include <o2scl/table3d.h>

// Forward definition of the tensor_grid class for HDF I/O
namespace o2scl {
  template<class vec_t, class vec_size_t> class tensor_grid;
}

// Forward definition of HDF I/O to extend friendship
namespace o2scl_hdf { 
  class hdf_file; 
  template<class vec_t, class vec_size_t>
    void hdf_input(hdf_file &hf, o2scl::tensor_grid<vec_t,vec_size_t> &t, 
		   std::string name);
  template<class vec_t, class vec_size_t>
    void hdf_output(hdf_file &hf, o2scl::tensor_grid<vec_t,vec_size_t> &t, 
		    std::string name);
}

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  typedef boost::numeric::ublas::range ub_range;
  typedef boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector<double> > ubvector_range;
  typedef boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector<size_t> > ubvector_size_t_range;
  
  /** \brief Tensor class with arbitrary dimensions with a grid
      
      This tensor class allows one to assign the indexes to numerical
      scales, effectively defining a data set on an n-dimensional
      grid. To set the grid, use \ref set_grid() or \ref
      set_grid_packed().
      
      By convention, member functions ending in the <tt>_val</tt>
      suffix return the closest grid-point to some user-specified
      values. 

      \comment
      I like how hist_new only allows you to set the
      grid if it's the same size as the old grid or the tensor
      is empty. This same practice could be applied here, to
      force the user to clear the tensor before resetting the grid.
      (10/24/11 - Actually, maybe this is a bad idea, because
      this class is more analogous to ubvector, which doesn't
      behave this way.)
      \endcomment

      \note Currently, HDF5 I/O is only allowed if the tensor is
      allocated with std::vector-based types, and the \ref
      interpolate() function only works with ublas-based vector types.

      \todo It is possible for the user to create a tensor_grid
      object, upcast it to a tensor object, and then use
      tensor::resize() to resize the tensor, failing to resize the
      grid. This probably needs fixing.

      \future Is it really necessary that get_data() is public?
      \future Only allocate space for grid if it is set
      \future Consider creating a new set_grid() function which
      takes grids from an object like uniform_grid. Maybe make a 
      constructor for a tensor_grid object which just takes 
      as input a set of grids?
  */
  template<class vec_t=std::vector<double>, 
    class vec_size_t=std::vector<size_t> > class tensor_grid :
    public tensor<vec_t,vec_size_t> {
    
  public:
  
#ifndef DOXYGEN_INTERNAL
  
  protected:

#ifdef O2SCL_NEVER_DEFINED
  }{
#endif  
    
    /// A rank-sized set of arrays for the grid points
    vec_t grid;
    
    /// If true, the grid has been set by the user
    bool grid_set;

    /// Interpolation type
    size_t itype;
    
#endif
    
  public:
    
  /// Create an empty tensor with zero rank
  tensor_grid() : tensor<vec_t,vec_size_t>() {
      grid_set=false;
      itype=itp_linear;
    }
    
    /** \brief Create a tensor of rank \c rank with sizes given in \c dim
	
	The parameter \c dim must be a vector of sizes with length \c
	rank. If the user requests any of the sizes to be zero, this
	constructor will call the error handler.
    */
    template<class size_vec_t> 
      tensor_grid(size_t rank, const size_vec_t &dim) : 
    tensor<vec_t,vec_size_t>(rank,dim) {
      grid_set=false;
      itype=itp_linear;
      // Note that the parent method sets rk to be equal to rank
      for(size_t i=0;i<this->rk;i++) {
	if (dim[i]==0) {
	  O2SCL_ERR((((std::string)"Requested zero size with non-zero ")+
		     "rank for index "+szttos(i)+" in tensor_grid::"+
		     "tensor_grid(size_t,size_vec_t)").c_str(),
		    exc_einval);
	}
      }
    }

    virtual ~tensor_grid() {
    }

    /// \name Set functions
    //@{
    /// Set the element closest to grid point \c grdp to value \c val
    template<class vec2_t> 
      void set_val(const vec2_t &grdp, double val) {
      
      // Find indices
      vec_size_t index(this->rk);
      for(size_t i=0;i<this->rk;i++) {
	index[i]=lookup_grid(i,grdp[i]);
      }
      
      // Pack
      size_t ix=index[0];
      for(size_t i=1;i<this->rk;i++) {
	ix*=this->size[i];
	ix+=index[i];
      }

      // Set value
      this->data[ix]=val;

      return;
    }

    /** \brief Set the element closest to grid point \c grdp to value 
	\c val
	
	The parameters \c closest and \c grdp may be identical,
	allowing one to update a vector \c grdp with the
	closest grid point.
    */
    template<class vec2_t, class vec3_t> 
      void set_val(const vec2_t &grdp, double val, vec3_t &closest) {

      // Find indices
      vec_size_t index(this->rk);
      for(size_t i=0;i<this->rk;i++) {
	index[i]=lookup_grid_val(i,grdp[i],closest[i]);
      }
      
      // Pack
      size_t ix=index[0];
      for(size_t i=1;i<this->rk;i++) {
	ix*=this->size[i];
	ix+=index[i];
      }

      // Set value
      this->data[ix]=val;

      return;
    }
    //@}

    /// \name Get functions
    //@{
    /// Get the element closest to grid point \c gridp 
    template<class vec2_t> double get_val(const vec2_t &gridp) {
      
      // Find indices
      vec_size_t index(this->rk);
      for(size_t i=0;i<this->rk;i++) {
	index[i]=lookup_grid(i,gridp[i]);
      }

      // Pack
      size_t ix=index[0];
      for(size_t i=1;i<this->rk;i++) {
	ix*=this->size[i];
	ix+=index[i];
      }

      // Set value
      return this->data[ix];
    }

    /** \brief Get the element closest to grid point \c gridp, 
	store grid values in \c closest and return value

	The parameters \c gridp and \c closest may refer to the
	same object. 
    */
    template<class vec2_t, class vec3_t> 
      double get_val(const vec2_t &gridp, vec3_t &closest) {
      
      // Find indices
      vec_size_t index(this->rk);
      for(size_t i=0;i<this->rk;i++) {
	index[i]=lookup_grid_val(i,gridp[i],closest[i]);
      }
      
      // Pack
      size_t ix=index[0];
      for(size_t i=1;i<this->rk;i++) {
	ix*=this->size[i];
	ix+=index[i];
      }

      // Set value
      return this->data[ix];
    }
    //@}
    
    /// \name Resize method
    //@{
    /** \brief Resize the tensor to rank \c rank with sizes
	given in \c dim
	
	The parameter \c dim must be a vector of sizes with a length
	equal to \c rank. This resize method is always destructive,
	and the grid is always reset. 
	
	If the user requests any of the sizes to be zero, this
	function will call the error handler.
    */
    template<class size_vec2_t>
      void resize(size_t rank, const size_vec2_t &dim) {
      // Double check that none of the sizes that the user
      // specified are zero
      for(size_t i=0;i<rank;i++) {
	if (dim[i]==0) {
	  O2SCL_ERR((((std::string)"Requested zero size with non-zero ")+
		     "rank for index "+szttos(i)+" in tensor_grid::"+
		     "resize().").c_str(),exc_einval);
	}
      }
      // Set the new rank
      this->rk=rank;
      // Resize the size vector
      this->size.resize(this->rk);
      // Reset the grid
      grid_set=false;
      grid.resize(0);
      // If the new rank is zero, reset the data, otherwise,
      // update the size vector and reset the data vector
      if (rank==0) {
	this->data.resize(0);
	return;
      } else {
	size_t tot=1;
	for(size_t i=0;i<this->rk;i++) {
	  this->size[i]=dim[i];
	  tot*=this->size[i];
	}
	this->data.resize(tot);
      }
      return;
    }

    //@}

    /// \name Grid manipulation
    //@{
    /// Return true if the grid has been set
    bool is_grid_set() const { return grid_set; }

    /** \brief Set the grid

	The grid must be specified for all of the dimensions at
	once. Denote \f$ (\mathrm{size})_0 \f$ as the size
	of the first dimension, \f$ (\mathrm{size})_1 \f$ as the
	size of the second dimesion, and so on. Then the 
	first \f$ (\mathrm{size})_0 \f$ entries in \c grid must
	be the grid for the first dimension, the next 
	\f$ (\mathrm{size})_1 \f$ entries must be the grid
	for the second dimension, and so on. Thus \c grid
	must be a vector of size
	\f[
	\sum_{i=0}^{\mathrm{rank}} (\mathrm{size})_i
	\f]

	Note that the grid is copied so the function argument may
	be destroyed by the user after calling set_grid_packed() without
	affecting the tensor grid. 

	\future Define a more generic interface for matrix types
    */
    template<class vec2_t>
      void set_grid_packed(const vec2_t &grid_vec) {
      if (this->rk==0) {
	O2SCL_ERR2("Tried to set grid for empty tensor in ",
		   "tensor_grid::set_grid_packed().",exc_einval);
      }
      size_t ngrid=0;
      for(size_t i=0;i<this->rk;i++) ngrid+=this->size[i];
      grid.resize(ngrid);
      for(size_t i=0;i<ngrid;i++) {
	grid[i]=grid_vec[i];
      }
      grid_set=true;
      return;
    }

    /** \brief Set grid from a vector of vectors of grid points
     */
    template<class vec_vec_t>
      void set_grid(const vec_vec_t &grid_vecs) {
      if (this->rk==0) {
	O2SCL_ERR2("Tried to set grid for empty tensor in ",
		   "tensor_grid::set_grid().",exc_einval);
      }
      size_t ngrid=0;
      for(size_t i=0;i<this->rk;i++) ngrid+=this->size[i];
      grid.resize(ngrid);
      size_t k=0;
      for(size_t i=0;i<this->rk;i++) {
	for(size_t j=0;j<this->size[i];j++) {
	  grid[k]=grid_vecs[i][j];
	  k++;
	}
      }
      grid_set=true;
      return;
    }

    /** \brief Copy grid for index \c i to vector \c v
	
	The type \c rvec_t must be a vector with a resize
	method. 
    */
    template<class rvec_t> void copy_grid(size_t i, rvec_t &v) {
      v.resize(this->size[i]);
      size_t istart=0;
      for(size_t k=0;k<i;k++) istart+=this->size[k];
      for(size_t k=0;k<this->size[i];k++) {
	v[k]=grid[istart+k];
      }
      return;
    }
    
    /// Lookup jth value on the ith grid
    double get_grid(size_t i, size_t j) const {
      if (!grid_set) {
	O2SCL_ERR("Grid not set in tensor_grid::get_grid().",
		  exc_einval);
      }
      if (i>=this->rk) {
	O2SCL_ERR((((std::string)"Index ")+szttos(i)+
		   " greater than or equal to rank, "+szttos(this->rk)+
		   ", in tensor_grid::get_grid().").c_str(),
		  exc_einval);
      }
      size_t istart=0;
      for(size_t k=0;k<i;k++) istart+=this->size[k];
      return grid[istart+j];
    }

    /// Set the jth value on the ith grid
    void set_grid(size_t i, size_t j, double val) {
      if (!grid_set) {
	O2SCL_ERR("Grid not set in tensor_grid::get_grid().",
		  exc_einval);
      }
      if (i>=this->rk) {
	O2SCL_ERR((((std::string)"Index ")+szttos(i)+
		   " greater than or equal to rank, "+szttos(this->rk)+
		   ", in tensor_grid::get_grid().").c_str(),
		  exc_einval);
      }
      size_t istart=0;
      for(size_t k=0;k<i;k++) istart+=this->size[k];
      grid[istart+j]=val;
      return;
    }

    /// Lookup index for grid closest to \c val
    size_t lookup_grid(size_t i, double val) {
      if (!grid_set) {
	O2SCL_ERR("Grid not set in tensor_grid::lookup_grid().",
		  exc_einval);
      }
      if (i>=this->rk) {
	O2SCL_ERR((((std::string)"Index ")+szttos(i)+
		   " greater than or equal to rank, "+szttos(this->rk)+
		   ", in tensor_grid::lookup_grid().").c_str(),
		  exc_einval);
      }
      size_t istart=0;
      
      for(size_t j=0;j<i;j++) {
	istart+=this->size[j];
      }
      size_t best=istart;
      double min=fabs(grid[istart]-val);
      for(size_t j=istart;j<istart+this->size[i];j++) {
	if (fabs(grid[j]-val)<min) {
	  best=j;
	  min=fabs(grid[j]-val);
	}
      }
      return best-istart;
    }

    /** \brief Lookup indices for grid closest point to \c vals

        The values in \c vals are not modified by this function.
	
	\comment
	This function must have a different name than 
	lookup_grid() because the template types cause
	confusion between the two functions.
	\endcomment
    */
    template<class vec2_t, class size_vec2_t>
      void lookup_grid_vec(const vec2_t &vals, size_vec2_t &indices) const {
      for(size_t k=0;k<this->rk;k++) {
	indices[k]=lookup_grid(k,vals[k]);
      }
      return;
    }

    /** \brief Lookup index for grid closest to \c val, returning the 
	grid point

	The parameters \c val and \c val2 may refer to the
	same object. 
    */
    size_t lookup_grid_val(size_t i, const double &val, double &val2) {
      if (i>=this->rk) {
	O2SCL_ERR((((std::string)"Index ")+szttos(i)+
		   " greater than or equal to rank, "+szttos(this->rk)+
		   ", in tensor_grid::lookup_grid_val().").c_str(),
		  exc_einval);
      }
      if (grid_set==false) {
	O2SCL_ERR("Grid not set in tensor_grid::lookup_grid_val().",
		  exc_einval);
      }

      // We have to store the grid point in a temporary in case 
      // the parameters gridp and closest refer to the same object
      double temp=val;

      size_t istart=0;
      for(size_t j=0;j<i;j++) istart+=this->size[j];
      size_t best=istart;
      double min=fabs(grid[istart]-temp);
      val2=grid[istart];
      for(size_t j=istart;j<istart+this->size[i];j++) {
	if (fabs(grid[j]-temp)<min) {
	  best=j;
	  min=fabs(grid[j]-temp);
	  val2=grid[j];
	}
      }
      return best-istart;
    }

    /// Lookup index for grid closest to \c val
    size_t lookup_grid_packed(size_t i, double val) {
      if (!grid_set) {
	O2SCL_ERR("Grid not set in tensor_grid::lookup_grid_packed().",
		  exc_einval);
      }
      if (i>=this->rk) {
	O2SCL_ERR((((std::string)"Index ")+szttos(i)+" greater than rank, "+
		   szttos(this->rk)+
		   ", in tensor_grid::lookup_grid_packed().").c_str(),
		  exc_einval);
      }
      size_t istart=0;
      for(size_t j=0;j<i;j++) istart+=this->size[j];
      size_t best=istart;
      double min=fabs(grid[istart]-val);
      for(size_t j=istart;j<istart+this->size[i];j++) {
	if (fabs(grid[j]-val)<min) {
	  best=j;
	  min=fabs(grid[j]-val);
	}
      }
      return best;
    }

    /// Lookup index for grid closest to \c val
    size_t lookup_grid_packed_val(size_t i, double val, double &val2) {
      if (!grid_set) {
	O2SCL_ERR("Grid not set in tensor_grid::lookup_grid_packed().",
		  exc_einval);
      }
      if (i>=this->rk) {
	O2SCL_ERR((((std::string)"Index ")+szttos(i)+" greater than rank, "+
		   szttos(this->rk)+
		   ", in tensor_grid::lookup_grid_packed().").c_str(),
		  exc_einval);
      }
      size_t istart=0;
      for(size_t j=0;j<i;j++) istart+=this->size[j];
      size_t best=istart;
      double min=fabs(grid[istart]-val);
      val2=grid[istart];
      for(size_t j=istart;j<istart+this->size[i];j++) {
	if (fabs(grid[j]-val)<min) {
	  best=j;
	  min=fabs(grid[j]-val);
	  val2=grid[j];
	}
      }
      return best;
    }
    //@}

    /// Return a reference to the data (for HDF I/O)
    vec_t &get_data() {
      return this->data;
    }
    
    /// \name Slicing
    //@{
    /** \brief Create a slice in a table3d object with an aligned
	grid

	This function uses the grid associated with indices \c ix_x
	and \c ix_y, and the tensor interpolation function to copy to
	the slice named \c slice_name in the table3d object \c tab .

	If the table3d object does not currently have a grid set, then
	the grid is automatically set to be the same as that stored in
	the tensor_grid object associated with ranks \c ix_x and \c
	iy_y. If the \ref o2scl::table3d object does have a grid set,
	then the values returned by \ref o2scl::table3d::get_nx() and
	\ref o2scl::table3d::get_ny() must be equal to the size of the
	tensor in indices \c ix_x and ix_y, respectively.

	This currently requires a copy of all the tensor data
	into the table3d object.
    */
    template<class size_vec2_t> 
      void copy_slice_align(size_t ix_x, size_t ix_y, size_vec2_t &index, 
			    table3d &tab, std::string slice_name) {
      
      if (ix_x>=this->rk || ix_y>=this->rk || ix_x==ix_y) {
	O2SCL_ERR2("Either indices greater than rank or x and y ind",
		   "ices equal in tensor_grid::copy_slice_align().",
		   exc_efailed);
      }

      // Get current table3d grid
      size_t nx, ny;
      tab.get_size(nx,ny);

      if (nx==0 && ny==0) {

	// If there's no grid, just create it
	std::vector<double> gx, gy;
	for(size_t i=0;i<this->size[ix_x];i++) {
	  gx.push_back(this->get_grid(ix_x,i));
	}
	for(size_t i=0;i<this->size[ix_y];i++) {
	  gy.push_back(this->get_grid(ix_y,i));
	}
	nx=gx.size();
	ny=gy.size();
	tab.set_xy("x",nx,gx,"y",ny,gy);
      }

      // Check that the grids are commensurate
      if (nx!=this->size[ix_x] || ny!=this->size[ix_y]) {
	O2SCL_ERR2("Grids not commensurate in ",
		   "tensor_grid::copy_slice_align().",exc_einval);
      }

      // Create slice if not already present
      size_t is;
      if (!tab.is_slice(slice_name,is)) tab.new_slice(slice_name);
      
      // Copy over data
      for(size_t i=0;i<nx;i++) {
	for(size_t j=0;j<ny;j++) {
	  index[ix_x]=i;
	  index[ix_y]=j;
	  double val=this->get(index);
	  tab.set(i,j,slice_name,val);
	}
      }
      
      return;
    }

    /** \brief Copy to a slice in a table3d object using interpolation

	This function uses the grid associated with indices \c ix_x
	and \c ix_y, and the tensor interpolation function to copy the
	tensor information to the slice named \c slice_name in the
	table3d object \c tab .
	
	\note This function uses the \ref tensor_grid::interp_linear() 
	for the interpolation.
    */
    template<class size_vec2_t> 
      void copy_slice_interp(size_t ix_x, size_t ix_y, size_vec2_t &index, 
			     table3d &tab, std::string slice_name) {

      if (ix_x>=this->rk || ix_y>=this->rk || ix_x==ix_y) {
	O2SCL_ERR2("Either indices greater than rank or x and y ",
		   "indices equal in tensor_grid::copy_slice_interp().",
		   exc_efailed);
      }

      // Get current table3d grid
      size_t nx, ny;
      tab.get_size(nx,ny);

      if (nx==0 && ny==0) {
	// If there's no grid, then just use the aligned version
	return copy_slice_align(ix_x,ix_y,index,tab,slice_name);
      }

      // Create vector of values to interpolate with
      std::vector<double> vals(this->rk);
      for(size_t i=0;i<this->rk;i++) {
	if (i!=ix_x && i!=ix_y) vals[i]=this->get_grid(i,index[i]);
      }

      // Create slice if not already present
      size_t is;
      if (!tab.is_slice(slice_name,is)) tab.new_slice(slice_name);

      // Loop through the table grid to perform the interpolation
      for(size_t i=0;i<nx;i++) {
	for(size_t j=0;j<ny;j++) {
	  vals[ix_x]=tab.get_grid_x(i);
	  vals[ix_y]=tab.get_grid_y(j);
	  tab.set(i,j,slice_name,this->interp_linear(vals));
	}
      }

      return;
    }
    //@}

    /// \name Clear method
    //@{
    /// Clear the tensor of all data and free allocated memory
    void clear() {
      grid.resize(0);
      grid_set=false;
      tensor<vec_t,vec_size_t>::clear();
      return;
    }
    //@}

    /// \name Interpolation
    //@{
    /// Set interpolation type for \ref interpolate()
    void set_interp_type(size_t interp_type) {
      itype=interp_type;
      return;
    }
  
    /** \brief Interpolate values \c vals into the tensor, 
	returning the result

	This is a quick and dirty implementation of n-dimensional
	interpolation by recursive application of the 1-dimensional
	routine from \ref interp_vec, using the base interpolation
	object specified in the template parameter \c base_interp_t.
	This will be very slow for sufficiently large data sets.

 	\note This function requires a range objects to obtain ranges
 	of vector objects. In ublas, this is done with
 	<tt>ublas::vector_range</tt> objects, so this function will
 	certainly work for \ref tensor_grid objects built on ublas
 	vector types. There is no corresponding <tt>std::range</tt>,
 	but you may be able to use either <tt>ublas::vector_range</tt>
 	or <tt>Boost.Range</tt> in order to call this function with
 	\ref tensor_grid objects built on <tt>std::vector</tt>.
 	However, this is not fully tested at the moment.

	\future It should be straightforward to improve the scaling of
	this algorithm significantly by creating a "window" of local
	points around the point of interest. This could be done easily
	by constructing an initial subtensor. However, this should
	probably be superceded by a more generic alternative which
	avoids explicit use of the 1-d interpolation types.
    */
    template<class range_t=ub_range,
      class data_range_t=ubvector_range, 
      class index_range_t=ubvector_size_t_range> 
      double interpolate(double *vals) {

      typedef interp_vec<vec_t> interp_t;
      
      if (this->rk==1) {
	
	interp_t si(this->size[0],grid,this->data,itype);
	return si.eval(vals[0]);

      } else {
	
	// Get total number of interpolations at this level
	size_t ss=1;
	for(size_t i=1;i<this->rk;i++) ss*=this->size[i];

	// Create space for y vectors and interpolators
	std::vector<vec_t> yvec(ss);
	std::vector<interp_t *> si(ss);
	for(size_t i=0;i<ss;i++) {
	  yvec[i].resize(this->size[0]);
	}
	
	// Create space for interpolation results
	tensor_grid tdat;
	index_range_t size_new(this->size,ub_range(1,this->rk));
	tdat.resize(this->rk-1,size_new);

	// Set grid for temporary tensor
	data_range_t grid_new(grid,ub_range(this->size[0],grid.size()));
	tdat.set_grid_packed(grid_new);
	
	// Create starting coordinate and counter
	vec_size_t co(this->rk);
	for(size_t i=0;i<this->rk;i++) co[i]=0;
	size_t cnt=0;

	// Loop over every interpolation
	bool done=false;
	while(done==false) {

	  // Fill yvector with the appropriate data
	  for(size_t i=0;i<this->size[0];i++) {
	    co[0]=i;
	    yvec[cnt][i]=this->get(co);
	  }
	  
	  si[cnt]=new interp_t(this->size[0],grid,yvec[cnt],itype);
	  
	  index_range_t co2(co,ub_range(1,this->rk));
	  tdat.set(co2,si[cnt]->eval(vals[0]));

	  // Go to next interpolation
	  cnt++;
	  co[this->rk-1]++;
	  // carry if necessary
	  for(int j=((int)this->rk)-1;j>0;j--) {
	    if (co[j]>=this->size[j]) {
	      co[j]=0;
	      co[j-1]++;
	    }
	  }

	  // Test if done
	  if (cnt==ss) done=true;

	  // End of while loop
	}
      
	// Now call the next level of interpolation
	double res=tdat.interpolate(vals+1);
	
	for(size_t i=0;i<ss;i++) {
	  delete si[i];
	}

	return res;
      }
    }

    /** \brief Perform a linear interpolation of \c v into the 
	function implied by the tensor and grid
	
	This performs multi-dimensional linear interpolation (or
	extrapolation) It works by first using \ref o2scl::search_vec
	to find the interval containing (or closest to) the specified
	point in each direction and constructing the corresponding
	hypercube of size \f$ 2^{\mathrm{rank}} \f$ containing \c v.
	It then calls \ref interp_linear_power_two() to perform the
	interpolation in that hypercube.

	\future This starts with a small copy, which can be eliminated
	by creating a new version of interp_linear_power_two
	which accepts an offset vector parameter so that the 
	first interpolation is offset. Remaining interpolations
	don't need to be offset because the tensor has to be
	created from the previous interpolation round.
    */
    template<class vec2_size_t> double interp_linear(vec2_size_t &v) {

      // Find the the corner of the hypercube containing v
      size_t rgs=0;
      std::vector<size_t> loc(this->rk);
      std::vector<double> gnew;
      for(size_t i=0;i<this->rk;i++) {
	std::vector<double> grid_unpacked(this->size[i]);
	for(size_t j=0;j<this->size[i];j++) {
	  grid_unpacked[j]=grid[j+rgs];
	}
	search_vec<std::vector<double> > sv(this->size[i],grid_unpacked);
	loc[i]=sv.find(v[i]);
	gnew.push_back(grid_unpacked[loc[i]]);
	gnew.push_back(grid_unpacked[loc[i]+1]);
	rgs+=this->size[i];
      }

      // Now construct a 2^{rk}-sized tensor containing only that 
      // hypercube
      std::vector<size_t> snew(this->rk);
      for(size_t i=0;i<this->rk;i++) {
	snew[i]=2;
      }
      tensor_grid tnew(this->rk,snew);
      tnew.set_grid_packed(gnew);
      
      // Copy over the relevant data
      for(size_t i=0;i<tnew.total_size();i++) {
	std::vector<size_t> index_new(this->rk), index_old(this->rk);
	tnew.unpack_indices(i,index_new);
	for(size_t j=0;j<this->rk;j++) index_old[j]=index_new[j]+loc[j];
	tnew.set(index_new,this->get(index_old));
      }
      
      // Now use interp_power_two()
      return tnew.interp_linear_power_two(v);
    }
    
    /** \brief Perform linear interpolation assuming that all
	indices can take only two values

	This function works by recursively slicing the hypercube of
	size \f$ 2^{\mathrm{rank}} \f$ into a hypercube of size \f$
	2^{\mathrm{rank-1}} \f$ performing linear interpolation for
	each pair of points.

	\note This is principally a function for internal use
	by \ref interp_linear().
    */
    template<class vec2_size_t>
      double interp_linear_power_two(vec2_size_t &v) {

      if (this->rk==1) {
	return this->data[0]+(this->data[1]-this->data[0])/
	  (grid[1]-grid[0])*(v[0]-grid[0]);
      }

      size_t last=this->rk-1;
      double frac=(v[last]-get_grid(last,0))/
	(get_grid(last,1)-get_grid(last,0));

      // Create new size vector and grid
      tensor_grid tnew(this->rk-1,this->size);
      tnew.set_grid_packed(grid);

      // Create data in new tensor, removing the last index through
      // linear interpolation
      for(size_t i=0;i<tnew.total_size();i++) {
	std::vector<size_t> index(this->rk);
	tnew.unpack_indices(i,index);
	index[this->rk-1]=0;
	double val_lo=this->get(index);
	index[this->rk-1]=1;
	double val_hi=this->get(index);
	tnew.set(index,val_lo+frac*(val_hi-val_lo));
      }

      // Recursively interpolate the smaller tensor
      return tnew.interp_linear_power_two(v);
    }

    /** \brief Perform a linear interpolation of <tt>v[1]</tt>
	to <tt>v[n-1]</tt> resulting in a vector

	\note The type <tt>vec2_t</tt> for the vector <tt>res</tt>
	must have a <tt>resize()</tt> method.
	
	This performs multi-dimensional linear interpolation (or
	extrapolation) in the last <tt>n-1</tt> indices of the
	rank-<tt>n</tt> tensor leaving the first index free and places
	the results in the vector \c res.
    */
    template<class vec2_size_t, class vec2_t>
      void interp_linear_vec0(vec2_size_t &v, vec2_t &res) {

      // Find the the corner of the hypercube containing v
      size_t rgs=0;
      std::vector<size_t> loc(this->rk);
      std::vector<double> gnew;
      for(size_t i=0;i<this->size[0];i++) {
	gnew.push_back(grid[i]);
      }
      rgs=this->size[0];
      loc[0]=0;
      for(size_t i=1;i<this->rk;i++) {
	std::vector<double> grid_unpacked(this->size[i]);
	for(size_t j=0;j<this->size[i];j++) {
	  grid_unpacked[j]=grid[j+rgs];
	}
	search_vec<std::vector<double> > sv(this->size[i],grid_unpacked);
	loc[i]=sv.find(v[i]);
	gnew.push_back(grid_unpacked[loc[i]]);
	gnew.push_back(grid_unpacked[loc[i]+1]);
	rgs+=this->size[i];
      }
      
      // Now construct a n*2^{rk-1}-sized tensor containing only that
      // hypercube
      std::vector<size_t> snew(this->rk);
      snew[0]=this->size[0];
      for(size_t i=1;i<this->rk;i++) {
	snew[i]=2;
      }
      tensor_grid tnew(this->rk,snew);
      tnew.set_grid_packed(gnew);
      
      // Copy over the relevant data
      for(size_t i=0;i<tnew.total_size();i++) {
	std::vector<size_t> index_new(this->rk), index_old(this->rk);
	tnew.unpack_indices(i,index_new);
	for(size_t j=0;j<this->rk;j++) {
	  index_old[j]=index_new[j]+loc[j];
	}
	tnew.set(index_new,this->get(index_old));
      }

      // Now use interp_power_two_vec0()
      tnew.interp_linear_power_two_vec0(v,res);

      return;
    }

    /** \brief Perform linear interpolation assuming that the last
	<tt>n-1</tt> indices can take only two values
	
	This function performs linear interpolation assuming that the
	last <tt>n-1</tt> indices can take only two values and placing
	the result into <tt>res</tt>.

	\note The type <tt>vec2_t</tt> for the vector <tt>res</tt>
 	must have a <tt>resize()</tt> method. This is principally a
 	function for internal use by \ref interp_linear_vec0().
    */
    template<class vec2_size_t, class vec2_t>
      void interp_linear_power_two_vec0(vec2_size_t &v, vec2_t &res) {
      
      if (this->rk==2) {
	size_t n=this->size[0];
	res.resize(n);
	vec_size_t ix0(2), ix1(2);
	ix0[1]=0;
	ix1[1]=1;
	for(size_t i=0;i<n;i++) {
	  ix0[0]=i;
	  ix1[0]=i;
	  res[i]=this->get(ix0)+(this->get(ix1)-this->get(ix0))/
	    (grid[n+1]-grid[n])*(v[1]-grid[n]);
	}
	return;
      }

      size_t last=this->rk-1;
      double frac=(v[last]-get_grid(last,0))/
	(get_grid(last,1)-get_grid(last,0));

      // Create new size vector and grid
      tensor_grid tnew(this->rk-1,this->size);
      tnew.set_grid_packed(grid);
      
      // Create data in new tensor, removing the last index through
      // linear interpolation
      for(size_t i=0;i<tnew.total_size();i++) {
	std::vector<size_t> index(this->rk);
	tnew.unpack_indices(i,index);
	index[this->rk-1]=0;
	double val_lo=this->get(index);
	index[this->rk-1]=1;
	double val_hi=this->get(index);
	tnew.set(index,val_lo+frac*(val_hi-val_lo));
      }
      
      // Recursively interpolate the smaller tensor
      tnew.interp_linear_power_two_vec0(v,res);

      return;
    }

    /** \brief Perform a linear interpolation of <tt>v</tt>
	into the tensor leaving one index free resulting in a vector
	
	This performs multi-dimensional linear interpolation (or
	extrapolation) in the last <tt>n-1</tt> indices of the
	rank-<tt>n</tt> tensor leaving the first index free and places
	the results in the vector \c res.

	\future This function could be more efficient.
    */
    template<class vec2_size_t, class vec2_t>
      void interp_linear_vec(vec2_size_t &v, size_t ifree, vec2_t &res) {

      size_t n=this->size[ifree];

      // This function uses interp_linear_power_two_vec0(), so it
      // works by remapping the indices. This defines the remapping.
      std::vector<size_t> map;
      map.push_back(ifree);
      for(size_t i=0;i<this->rk;i++) {
	if (i!=ifree) {
	  map.push_back(i);
	}
      }

      // Find the the corner of the hypercube containing v
      size_t rgs=0;
      vec_size_t loc(this->rk);
      loc[ifree]=0;
      for(size_t i=0;i<this->rk;i++) {
	vec_t grid_unpacked(this->size[i]);
	for(size_t j=0;j<this->size[i];j++) {
          grid_unpacked[j]=grid[j+rgs];
        }
	search_vec<vec_t> sv(this->size[i],grid_unpacked);
	if (i!=ifree) {
	  loc[i]=sv.find(v[i]);
	}
	rgs+=this->size[i];
      }

      // Compute the remapped grid and interpolating vector
      std::vector<double> gnew, vnew;
      for(size_t new_ix=0;new_ix<this->rk;new_ix++) {
	for(size_t old_ix=0;old_ix<this->rk;old_ix++) {
	  if (map[new_ix]==old_ix) {
	    vnew.push_back(v[old_ix]);
	    if (old_ix==ifree) {
	      for(size_t j=0;j<this->size[old_ix];j++) {
		gnew.push_back(this->get_grid(old_ix,j));
	      }
	    } else {
	      gnew.push_back(this->get_grid(old_ix,loc[old_ix]));
	      gnew.push_back(this->get_grid(old_ix,loc[old_ix]+1));
	    }
	  }
	}
      }

      // Now construct a n*2^{rk-1}-sized tensor containing only the
      // hypercube needed to do the interpolation

      // Specify the size of each rank
      std::vector<size_t> snew;
      snew.push_back(n);
      for(size_t i=0;i<this->rk;i++) {
	if (i!=ifree) {
	  snew.push_back(2);
	}
      }

      // Create the tensor and set the grid
      tensor_grid tnew(this->rk,snew);
      tnew.set_grid_packed(gnew);

      // Copy over the relevant data
      for(size_t i=0;i<tnew.total_size();i++) {
	std::vector<size_t> index_new(this->rk), index_old(this->rk);
	tnew.unpack_indices(i,index_new);
	for(size_t j=0;j<this->rk;j++) {
	  index_old[map[j]]=index_new[j]+loc[map[j]];
	}
	tnew.set(index_new,this->get(index_old));
      }

      // Now use interp_power_two_vec()
      tnew.interp_linear_power_two_vec0(vnew,res);

      return;
    }
    //@}
    
    template<class vecf_t, class vecf_size_t> friend void o2scl_hdf::hdf_output
      (o2scl_hdf::hdf_file &hf, tensor_grid<vecf_t,vecf_size_t> &t, 
       std::string name);
    
    template<class vecf_t, class vecf_size_t> friend void o2scl_hdf::hdf_input
      (o2scl_hdf::hdf_file &hf, tensor_grid<vecf_t,vecf_size_t> &t, 
       std::string name);

  };

  /** \brief Rank 1 tensor with a grid
      
      \future Make rank-specific get_val and set_val functions?
  */
  template<class vec_t=std::vector<double>, 
    class vec_size_t=std::vector<size_t> > class tensor_grid1 : 
    public tensor_grid<vec_t,vec_size_t> {
     
  public:
     
  /// Create an empty tensor
  tensor_grid1() : tensor_grid<vec_t,vec_size_t>() {}
      
  /// Create a rank 2 tensor of size \c (sz,sz2,sz3)
  tensor_grid1(size_t sz) : tensor_grid<vec_t,vec_size_t>() {
      this->rk=1;
      this->size.resize(1);
      this->size[0]=sz;
      this->data.resize(sz);
      this->grid_set=false;
    }
  
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif  
    
    virtual ~tensor_grid1() {
    }
    
    /// Get the element indexed by \c (ix1)
    double &get(size_t ix1) { 
      size_t sz[1]={ix1};
      return tensor_grid<vec_t,vec_size_t>::get(sz); 
    }
    
    /// Get the element indexed by \c (ix1)
    const double &get(size_t ix1) const { 
      size_t sz[1]={ix1};
      return tensor_grid<vec_t,vec_size_t>::get(sz); 
    }
    
    /// Set the element indexed by \c (ix1) to value \c val
    void set(size_t ix1, double val) {
      size_t sz[1]={ix1};
      tensor_grid<vec_t,vec_size_t>::set(sz,val); 
    }
    
    /// Interpolate \c x and return the results
    template<class range_t=ub_range, class data_range_t=ubvector_range, 
      class index_range_t=ubvector_size_t_range> 
      double interp(double x) {
      return tensor_grid<vec_t,vec_size_t>::template interpolate
      <range_t,data_range_t,index_range_t>(&x);
    }
    
    /// Interpolate \c x and return the results
    double interp_linear(double x) {
      double arr[1]={x};
      return tensor_grid<vec_t,vec_size_t>::interp_linear(arr);
    }

  };
  
  /** \brief Rank 2 tensor with a grid
   */
  template<class vec_t=std::vector<double>, 
    class vec_size_t=std::vector<size_t> > class tensor_grid2 : 
    public tensor_grid<vec_t,vec_size_t> {
     
  public:
     
  /// Create an empty tensor
  tensor_grid2() : tensor_grid<vec_t,vec_size_t>() {}

  /// Create a rank 2 tensor of size \c (sz,sz2)
  tensor_grid2(size_t sz, size_t sz2) : tensor_grid<vec_t,vec_size_t>() {
      this->rk=2;
      this->size.resize(2);
      this->size[0]=sz;
      this->size[1]=sz2;
      size_t tot=sz*sz2;
      this->data.resize(tot);
      this->grid_set=false;
    }
   
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif  

    virtual ~tensor_grid2() {
    }
    
    /// Get the element indexed by \c (ix1,ix2)
    double &get(size_t ix1, size_t ix2) { 
      size_t sz[2]={ix1,ix2};
      return tensor_grid<vec_t,vec_size_t>::get(sz); 
    }
    
    /// Get the element indexed by \c (ix1,ix2)
    const double &get(size_t ix1, size_t ix2) const { 
      size_t sz[2]={ix1,ix2};
      return tensor_grid<vec_t,vec_size_t>::get(sz); 
    }
    
    /// Set the element indexed by \c (ix1,ix2) to value \c val
    void set(size_t ix1, size_t ix2, double val) {
      size_t sz[2]={ix1,ix2};
      tensor_grid<vec_t,vec_size_t>::set(sz,val); 
      return;
    }
    
    /// Interpolate \c (x,y) and return the results
    template<class range_t=ub_range, class data_range_t=ubvector_range, 
      class index_range_t=ubvector_size_t_range> 
      double interp(double x, double y) {
      double arr[2]={x,y};
      return tensor_grid<vec_t,vec_size_t>::template interpolate
      <range_t,data_range_t,index_range_t>(arr);
    }
    
    /// Interpolate \c (x,y) and return the results
    double interp_linear(double x, double y) {
      double arr[2]={x,y};
      return tensor_grid<vec_t,vec_size_t>::interp_linear(arr);
    }

  };
  
  /** \brief Rank 3 tensor with a grid
   */
  template<class vec_t=std::vector<double>, 
    class vec_size_t=std::vector<size_t> > class tensor_grid3 : 
    public tensor_grid<vec_t,vec_size_t> {
     
  public:
     
  /// Create an empty tensor
  tensor_grid3() : tensor_grid<vec_t,vec_size_t>() {}

  /// Create a rank 3 tensor of size \c (sz,sz2,sz3)
  tensor_grid3(size_t sz, size_t sz2, size_t sz3) : 
  tensor_grid<vec_t,vec_size_t>() {
    this->rk=3;
    this->size.resize(3);
    this->size[0]=sz;
    this->size[1]=sz2;
    this->size[2]=sz3;
    size_t tot=sz*sz2*sz3;
    this->data.resize(tot);
    this->grid_set=false;
  }
   
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif  

    virtual ~tensor_grid3() {
    }
   
    /// Get the element indexed by \c (ix1,ix2,ix3)
    double &get(size_t ix1, size_t ix2, size_t ix3) { 
      size_t sz[3]={ix1,ix2,ix3};
      return tensor_grid<vec_t,vec_size_t>::get(sz); 
    }
 
    /// Get the element indexed by \c (ix1,ix2,ix3)
    const double &get(size_t ix1, size_t ix2, size_t ix3) const { 
      size_t sz[3]={ix1,ix2,ix3};
      return tensor_grid<vec_t,vec_size_t>::get(sz); 
    }
 
    /// Set the element indexed by \c (ix1,ix2,ix3) to value \c val
    void set(size_t ix1, size_t ix2, size_t ix3, double val) {
      size_t sz[3]={ix1,ix2, ix3};
      tensor_grid<vec_t,vec_size_t>::set(sz,val); 
      return;
    }
    
    /// Interpolate \c (x,y,z) and return the results
    template<class range_t=ub_range, class data_range_t=ubvector_range, 
      class index_range_t=ubvector_size_t_range> 
      double interp(double x, double y, double z) {
      double arr[3]={x,y,z};
      return tensor_grid<vec_t,vec_size_t>::template interpolate
      <range_t,data_range_t,index_range_t>(arr);
    }
    
    /// Interpolate \c (x,y,z) and return the results
    double interp_linear(double x, double y, double z) {
      double arr[3]={x,y,z};
      return tensor_grid<vec_t,vec_size_t>::interp_linear(arr);
    }

  };
  
  /** \brief Rank 4 tensor with a grid
   */
  template<class vec_t=std::vector<double>, 
    class vec_size_t=std::vector<size_t> > class tensor_grid4 : 
    public tensor_grid<vec_t,vec_size_t> {
     
  public:
     
  /// Create an empty tensor
  tensor_grid4() : tensor_grid<vec_t,vec_size_t>() {}

  /// Create a rank 4 tensor of size \c (sz,sz2,sz3,sz4)
  tensor_grid4(size_t sz, size_t sz2, size_t sz3,
	       size_t sz4) : tensor_grid<vec_t,vec_size_t>() {
      this->rk=4;
      this->size.resize(4);
      this->size[0]=sz;
      this->size[1]=sz2;
      this->size[2]=sz3;
      this->size[3]=sz4;
      size_t tot=sz*sz2*sz3*sz4;
      this->data.resize(tot);
      this->grid_set=false;
    }
   
#ifdef O2SCL_NEVER_DEFINED
  }{
#endif  

    virtual ~tensor_grid4() {
    }
   
    /// Get the element indexed by \c (ix1,ix2,ix3,ix4)
    double &get(size_t ix1, size_t ix2, size_t ix3, size_t ix4) { 
      size_t sz[4]={ix1,ix2,ix3,ix4};
      return tensor_grid<vec_t,vec_size_t>::get(sz); 
    }
    
    /// Get the element indexed by \c (ix1,ix2,ix3,ix4)
    const double &get(size_t ix1, size_t ix2, size_t ix3,
		      size_t ix4) const { 
      size_t sz[4]={ix1,ix2,ix3,ix4};
      return tensor_grid<vec_t,vec_size_t>::get(sz); 
    }
    
    /// Set the element indexed by \c (ix1,ix2,ix3,ix4) to value \c val
    void set(size_t ix1, size_t ix2, size_t ix3, size_t ix4,
	     double val) {
      size_t sz[4]={ix1,ix2,ix3,ix4};
      tensor_grid<vec_t,vec_size_t>::set(sz,val); 
      return;
    }

    /// Interpolate \c (x,y,z,a) and return the results
    template<class range_t=ub_range, class data_range_t=ubvector_range, 
      class index_range_t=ubvector_size_t_range> 
      double interp(double x, double y, double z, double a) {
      double arr[4]={x,y,z,a};
      return tensor_grid<vec_t,vec_size_t>::template interpolate
      <range_t,data_range_t,index_range_t>(arr);
    }

    /// Interpolate \c (x,y,z,a) and return the results
    double interp_linear(double x, double y, double z, double a) {
      double arr[4]={x,y,z,a};
      return tensor_grid<vec_t,vec_size_t>::interp_linear(arr);
    }

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif



