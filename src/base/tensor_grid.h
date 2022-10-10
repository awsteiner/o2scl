/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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

namespace o2scl {

  /** \brief A <tt>ublas::range</tt> typedef for \ref 
      o2scl::tensor_grid and related classes in src/base/tensor_grid.h
   */
  typedef boost::numeric::ublas::range ub_range;

  /** \brief A <tt>ublas::vector_range</tt> typedef for \ref 
      o2scl::tensor_grid and related classes in src/base/tensor_grid.h
   */
  typedef boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector<double> > ubvector_range;

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;
  
  /** \brief A <tt>ublas::vector_range</tt> typedef (size_t version)
      for \ref o2scl::tensor_grid and related classes
      in src/base/tensor_grid.h
   */
  typedef boost::numeric::ublas::vector_range
    <boost::numeric::ublas::vector<size_t> > ubvector_size_t_range;
  
  /** \brief Tensor class with arbitrary dimensions with a grid
      
      This tensor class allows one to assign the indexes to numerical
      scales, effectively defining a data set on an n-dimensional
      grid. To set the grid, use \ref default_grid(), \ref set_grid()
      or \ref set_grid_packed().
      
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

      \b Slicing

      New \ref o2scl::tensor_grid objects can be obtained
      by fixing any set of indices using \ref copy_slice_interp().

      Fixing all but two indices also results in a \ref o2scl::table3d
      object, and five functions perform this task in different ways.
      The function \ref copy_table3d_align() copies a two-dimensional
      slice to a \ref o2scl::table3d object presuming that the grid in
      the \ref o2scl::table3d object has already been set and exactly
      matches the corresponding sizes for the selected tensor indices.
      This function does not check that the grids between the two
      objects match, it only ensures that they have the same size. In
      order to copy to a \ref o2scl::table3d object and set its grid
      to match that from the unfixed indices in the \ref
      o2scl::tensor_grid object, the function \ref
      copy_table3d_align_setxy() can be used. The function \ref
      copy_table3d_interp() uses interpolation to extract values from
      the \ref o2scl::tensor_grid object. It allows the user to select
      indices to be fixed and then uses the values in the grid in the
      \ref o2scl::table3d object for the indices which vary.
      Alternatively \ref copy_table3d_interp_values() allows the user
      to specify values on the grid for the indices to be fixed and
      uses the grid in the \ref o2scl::table3d object for the indices
      which vary. Finally, \ref copy_table3d_interp_values_setxy() acts
      like \ref copy_table3d_interp_values() except that it sets the
      \ref o2scl::table3d grid to be the same as the grid in the \ref
      o2scl::tensor_grid object which corresponds to the indices which
      are being varied.

      \b Notes and Todos

      \note Currently, HDF5 I/O is only allowed if the tensor is
      allocated with std::vector-based types, and the \ref
      interpolate() function only works with ublas-based vector types.

      \verbatim embed:rst

      .. todo:: 

         Class tensor_grid: It is possible for the user to create a
         tensor_grid object, upcast it to a tensor object, and then
         use tensor::resize() to resize the tensor, failing to resize
         the grid. Following this, grid access functions will access
         random parts of memory or segfault. This can be fixed by
         ensuring that resize functions are virtual and have a version
         in tensor_grid which ensure that the grid and tensor data are
         matched. The problem is that the resize functions are
         templates, so they cannot be virtual.

         - Future: Create a swap function for the grid similar to the
           data swap function in the parent \ref o2scl::tensor class?

         - Future: Only allocate space for grid if it is set.

         - Future: as with \ref o2scl::tensor, generalize to other
           base data types.

         - Future: The function \ref interp_linear_partial() appears
           to be a generalization of \ref
           copy_table3d_interp_values_setxy(), so there may be some code
           duplication between the two that can be avoided.

      \endverbatim

  */
  template<class vec_t=std::vector<double>, 
    class vec_size_t=std::vector<size_t> > class tensor_grid :
    public tensor<double,vec_t,vec_size_t> {
    
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

    /// \name Constructors and Destructors
    //@{
    /// Create an empty tensor with zero rank
  tensor_grid() : tensor<double,vec_t,vec_size_t>() {
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
    tensor<double,vec_t,vec_size_t>(rank,dim) {
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
    
    /** \brief Create a tensor with a grid defined by a set
        of \ref o2scl::uniform_grid objects
    */
  tensor_grid(std::vector<uniform_grid<double> > &ugs) :
    tensor<double,vec_t,vec_size_t>() {    
      this->rk=ugs.size();
      itype=itp_linear;
      size_t tot=1;
      for(size_t j=0;j<this->rk;j++) {
        this->size.push_back(ugs[j].get_npoints());
        tot*=ugs[j].get_npoints();
      }
      this->data.resize(tot);
      set_grid(ugs);
    }

    /** \brief Destructor
     */
    virtual ~tensor_grid() {
    }
    //@}

    /// \name Method to check for valid object
    //@{
    /** \brief Check that the \ref o2scl::tensor_grid object is valid
     */
    void is_valid() const {
      
      tensor<double,vec_t,vec_size_t>::is_valid();
      
      if (this->rk>0 && grid_set) {
        size_t tot2=0;
        for(size_t i=0;i<this->rk;i++) {
          tot2+=this->size[i];
        }
        
        if (tot2!=grid.size()) {
          O2SCL_ERR2("Value grid_set is true but grid vector size ",
                     "is wrong in tensor_grid::is_valid().",
                     o2scl::exc_esanity);
        }
      }

      if (!grid_set && grid.size()>0) {
        O2SCL_ERR2("Value grid_set is false but grid vector size ",
                   "is not zero in tensor_grid::is_valid().",
                   o2scl::exc_esanity);
      }
      
      return;
    }
    //@}
    
    /// \name Copy constructors
    //@{
    /** \brief Copy using <tt>operator()</tt>
     */
    tensor_grid<vec_t,vec_size_t>
      (const tensor_grid<vec_t,vec_size_t> &t) {
      this->rk=t.rk;
      this->data=t.data;
      this->size=t.size;
      grid=t.grid;
      grid_set=t.grid_set;
      itype=t.itype;
    }
    
    /** \brief Copy using <tt>operator=()</tt>
     */
    tensor_grid<vec_t,vec_size_t> &operator=
      (const tensor_grid<vec_t,vec_size_t> &t) {
      if (this!=&t) {
        this->rk=t.rk;
        this->data=t.data;
        this->size=t.size;
        grid=t.grid;
        grid_set=t.grid_set;
        itype=t.itype;
      }
      return *this;
    }
    //@}
  
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

    /// Get grid
    const vec_t &get_grid() const {
      return grid;
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

        \verbatim embed:rst
        
        .. todo:: 

           In set_grid_packed(): Define a more generic interface for
           matrix types

        \endverbatim
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

    /** \brief Use a default grid which just uses the index
     */
    void default_grid() {
      size_t ngrid=0;
      for(size_t i=0;i<this->rk;i++) ngrid+=this->size[i];
      grid.resize(ngrid);
      size_t k=0;
      for(size_t i=0;i<this->rk;i++) {
        for(size_t j=0;j<this->size[i];j++) {
          grid[k]=((double)j);
          k++;
        }
      }
      grid_set=true;
      return;
    }
    
    /** \brief Set grid for one index from a vector
     */
    template<class vec2_t>
      void set_grid_i_vec(size_t ix, const vec2_t &grid_vec) {
      if (grid_set==false) {
        O2SCL_ERR2("Grid not already set in ",
                   "tensor_grid::set_grid_i_vec().",exc_einval);
      }
      if (this->rk==0) {
        O2SCL_ERR2("Tried to set grid for empty tensor in ",
                   "tensor_grid::set_grid_i_vec().",exc_einval);
      }
      size_t k=0;
      for(size_t i=0;i<this->rk;i++) {
        for(size_t j=0;j<this->size[i];j++) {
          if (j==ix) {
            grid[k]=grid_vec[j];
          }
          k++;
        }
      }
      return;
    }

    /** \brief Set grid for one index from a function
     */
    void set_grid_i_func(size_t ix, std::string func) {
      if (grid_set==false) {
        O2SCL_ERR2("Grid not already set in ",
                   "tensor_grid::set_grid_i_func().",exc_einval);
      }
      if (this->rk==0) {
        O2SCL_ERR2("Tried to set grid for empty tensor in ",
                   "tensor_grid::set_grid_i_func().",exc_einval);
      }
      
      calc_utf8<> calc;
      std::map<std::string,double> vars;
      calc.compile(func.c_str(),&vars);
      
      size_t k=0;
      for(size_t i=0;i<this->rk;i++) {
        for(size_t j=0;j<this->size[i];j++) {
          if (j==ix) {
            vars["i"]=((double)j);
            grid[k]=calc.eval(&vars);
          }
          k++;
        }
      }
      
      return;
    }

    /** \brief Set grid from a vector of uniform grid objects

        \note This is called by one of the constructors.
    */
    void set_grid(std::vector<uniform_grid<double> > &ugs) {
      if (this->rk==0) {
        O2SCL_ERR2("Tried to set grid for empty tensor in ",
                   "tensor_grid::set_grid().",exc_einval);
      }
      size_t ngrid=0;
      for(size_t i=0;i<this->rk;i++) ngrid+=ugs[i].get_npoints();
      grid.resize(ngrid);
      size_t k=0;
      for(size_t i=0;i<this->rk;i++) {
        for(size_t j=0;j<ugs[i].get_npoints();j++) {
          grid[k]=ugs[i][j];
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
    template<class rvec_t> void copy_grid(size_t i, rvec_t &v) const {
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

    /** \brief Lookup index for grid closest to \c val, returning the 
        grid point

        The parameters \c val and \c val2 may refer to the
        same object. 
    */
    size_t lookup_grid_val(size_t i, const double &val, double &val2) const {
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
    size_t lookup_grid(size_t i, double val) const {
      double val2;
      return lookup_grid_val(i,val,val2);
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

    /** \brief Lookup internal packed grid index for point closest 
        to \c val and store closest value in \c val2

        This version, rather than \ref
        o2scl::tensor_grid::lookup_grid_val() can be useful because it
        gives the index of the grid point in the internal grid vector
        object.
    */
    size_t lookup_grid_packed_val(size_t i, double val, double &val2) const {
      if (!grid_set) {
        O2SCL_ERR("Grid not set in tensor_grid::lookup_grid_packed_val().",
                  exc_einval);
      }
      if (i>=this->rk) {
        O2SCL_ERR((((std::string)"Index ")+szttos(i)+" greater than rank, "+
                   szttos(this->rk)+
                   ", in tensor_grid::lookup_grid_packed_val().").c_str(),
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
    
    /** \brief Lookup internal packed grid index for point closest 
        to \c val
    */
    size_t lookup_grid_packed(size_t i, double val) const {
      double val2;
      return lookup_grid_packed_val(i,val,val2);
    }
    //@}

    /// \name Slicing to tensor_grid objects
    //@{
    /** \brief Copy an abitrary slice by fixing 1 or more indices and
        use interpolation to return a new \ref tensor_grid object
    */
    template<class size_vec2_t, class vec2_t> 
      tensor_grid<> copy_slice_interp(size_vec2_t &ifix, vec2_t &vals) const {

      if (this->rk<1+ifix.size()) {
        O2SCL_ERR2("Fixed too many indices in ",
                   "tensor_grid::copy_slice_interp().",
                   o2scl::exc_einval);
      }
      if (ifix.size()!=vals.size()) {
        O2SCL_ERR2("Mismatch between indices and values in ",
                   "tensor_grid::copy_slice_interp().",
                   o2scl::exc_einval);
      }
      
      // Determine new rank
      size_t rank_new=this->rk-ifix.size();

      // Determine the new sizes and new grid
      std::vector<size_t> sz_new;
      std::vector<std::vector<double> > grid_new;
      for(size_t i=0;i<this->rk;i++) {
        bool found=false;
        for(size_t j=0;j<ifix.size();j++) {
          if (ifix[j]==i) found=true;
        }
        if (found==false) {
          sz_new.push_back(this->get_size(i));
          std::vector<double> grid_temp;
          for(size_t j=0;j<this->get_size(i);j++) {
            grid_temp.push_back(this->get_grid(i,j));
          }
          grid_new.push_back(grid_temp);
        }
      }

      // Create the new tensor_grid object and set the new grid
      tensor_grid<> tg_new(rank_new,sz_new);
      tg_new.set_grid(grid_new);

      // Interpolate the data into the new tensor_grid object
      std::vector<size_t> ix_new(rank_new);
      std::vector<double> point_old(this->rk);

      // Loop over the new tensor_grid object
      for(size_t i=0;i<tg_new.total_size();i++) {

        // Find the location in the new tensor_grid object
        tg_new.unpack_index(i,ix_new);

        // Find the point in the old tensor object to interpolate
        for(size_t j=0;j<this->rk;j++) {
          int ix_found=-1;
          for(size_t k=0;k<ifix.size();k++) {
            if (ifix[k]==j) ix_found=k;
          }
          if (ix_found==-1) {
            point_old[j]=this->get_grid(j,ix_new[j]);
          } else {
            point_old[j]=vals[ix_found];
          }
        }
        
        // Set the new point by performing the linear interpolation
        tg_new.set(ix_new,this->interp_linear(point_old));
      }

      return tg_new;
    }
    //@}

    /// \name Slicing and converting to table3d objects
    //@{

    /** \brief Convert to a \ref o2scl::table3d object by
        summing over all but two indices
    */
    void convert_table3d_sum
      (size_t ix_x, size_t ix_y, table3d &tab, std::string x_name="x",
       std::string y_name="y", std::string slice_name="z") const {
      
      // Get current table3d grid
      size_t nx, ny;
      tab.get_size(nx,ny);
      
      if (nx==0 && ny==0) {
        
        if (x_name.length()==0) x_name="x";
        if (y_name.length()==0) y_name="y";
        
        // If there's no grid, then create a grid in the table3d
        // object that is the same as that in the tensor_grid object
        std::vector<double> grid_x, grid_y;
        copy_grid(ix_x,grid_x);
        copy_grid(ix_y,grid_y);
        tab.set_xy("x",grid_x.size(),grid_x,
                   "y",grid_y.size(),grid_y);
        // Now that the grid is set, get nx and ny
        tab.get_size(nx,ny);
      }
      
      // Check that the grids are commensurate
      if (nx!=this->size[ix_x] || ny!=this->size[ix_y]) {
        O2SCL_ERR2("Grids not commensurate in ",
                   "tensor_grid::convert_table3d_sum().",exc_einval);
      }
      
      tab.set_slice_all(slice_name,0.0);
      
      std::vector<size_t> ix;
      for(size_t i=0;i<this->total_size();i++) {
        this->unpack_index(i,ix);
        tab.set(ix[ix_x],ix[ix_y],slice_name,
                tab.get(ix[ix_x],ix[ix_y],slice_name)+
                this->data[i]);
      }
      
      return;
    }

    /** \brief Create a slice in a \ref o2scl::table3d object with an
        aligned grid

        This function uses the grid associated with indices \c ix_x
        and \c ix_y, to copy data to a slice named \c slice_name in
        the table3d object \c tab . All other indices are fixed
        to values specified by the user in \c index and the values
        of <tt>index[ix_x]</tt> and <tt>index[ix_y]</tt> are
        used for temporary storage.

        If the table3d object does not currently have a grid set, then
        the grid is automatically set to be the same as that stored in
        the tensor_grid object associated with ranks \c ix_x and \c
        iy_y. If the \ref o2scl::table3d object does have a grid set,
        then the values returned by \ref o2scl::table3d::get_nx() and
        \ref o2scl::table3d::get_ny() must be equal to the size of the
        tensor in indices \c ix_x and ix_y, respectively. If a
        slice named \c slice_name is not already present in 
        \c tab, then a new slice with that name is created.

        The error handler is called if \c ix_x is the same as 
        \c ix_y, or if either of these two values is greater
        than or equal to the tensor rank.
    */
    template<class size_vec2_t> 
      void copy_table3d_align(size_t ix_x, size_t ix_y, size_vec2_t &index, 
                              table3d &tab, std::string slice_name="z") const {
      
      if (ix_x>=this->rk || ix_y>=this->rk || ix_x==ix_y) {
        O2SCL_ERR2("Either indices greater than rank or x and y ind",
                   "ices equal in tensor_grid::copy_table3d_align().",
                   exc_efailed);
      }

      // Get current table3d grid
      size_t nx, ny;
      tab.get_size(nx,ny);

      // Check that the grids are commensurate
      if (nx!=this->size[ix_x] || ny!=this->size[ix_y]) {
        O2SCL_ERR2("Grids not commensurate in ",
                   "tensor_grid::copy_table3d_align().",exc_einval);
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
    
    /** \brief Create a slice in a table3d object with a new aligned
        grid
    */
    template<class size_vec2_t> 
      void copy_table3d_align_setxy
      (size_t ix_x, size_t ix_y, size_vec2_t &index, 
       table3d &tab, std::string x_name="x", std::string y_name="y",
       std::string slice_name="z") const {

      // Get current table3d grid
      size_t nx, ny;
      tab.get_size(nx,ny);

      if (nx==0 && ny==0) {

        if (x_name.length()==0) x_name="x";
        if (y_name.length()==0) y_name="y";
        
        // If there's no grid, then create a grid in the table3d
        // object that is the same as that in the tensor_grid object
        std::vector<double> grid_x, grid_y;
        copy_grid(ix_x,grid_x);
        copy_grid(ix_y,grid_y);
        tab.set_xy("x",grid_x.size(),grid_x,
                   "y",grid_y.size(),grid_y);
        // Now that the grid is set, get nx and ny
        tab.get_size(nx,ny);
      }
      
      copy_table3d_align(ix_x,ix_y,index,tab,slice_name);
      
      return;
    }
    
    /** \brief Copy to a slice in a table3d object using interpolation

        This function uses the grid associated with indices \c ix_x
        and \c ix_y, and the tensor interpolation function to copy the
        tensor information to the slice named \c slice_name in the
        table3d object \c tab . All other indices are fixed
        to values specified by the user in \c index and the values
        of <tt>index[ix_x]</tt> and <tt>index[ix_y]</tt> are
        used for temporary storage.
        
        If a slice named \c slice_name is not already present in \c
        tab, then a new slice with that name is created.

        The error handler is called if \c ix_x is the same as 
        \c ix_y, or if either of these two values is greater
        than or equal to the tensor rank.

        \note This function uses the \ref tensor_grid::interp_linear() 
        for the interpolation.
    */
    template<class size_vec2_t> 
      void copy_table3d_interp(size_t ix_x, size_t ix_y, size_vec2_t &index, 
                               table3d &tab, std::string slice_name="z") const {

      if (ix_x>=this->rk || ix_y>=this->rk || ix_x==ix_y) {
        O2SCL_ERR2("Either indices greater than rank or x and y ",
                   "indices equal in tensor_grid::copy_table3d_interp().",
                   exc_efailed);
      }

      // Get current table3d grid
      size_t nx, ny;
      tab.get_size(nx,ny);

      if (nx==0 && ny==0) {
        // If there's no grid, then just use the aligned version
        return copy_table3d_align(ix_x,ix_y,index,tab,slice_name);
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
    
    /** \brief Copy to a slice in a table3d object using interpolation
     */
    template<class vec2_t> 
      void copy_table3d_interp_values(size_t ix_x, size_t ix_y,
                                      vec2_t &values, table3d &tab,
                                      std::string slice_name="z",
                                      int verbose=0) const {
      
      if (ix_x>=this->rk || ix_y>=this->rk || ix_x==ix_y) {
        O2SCL_ERR2("Either indices greater than rank or x and y ",
                   "indices equal in tensor_grid::copy_table3d_interp().",
                   exc_efailed);
      }
      if (values.size()!=this->rk) {
        O2SCL_ERR2("Values array not equal to rank ",
                   "in tensor_grid::copy_table3d_interp_values().",
                   exc_efailed);
      }

      if (tab.is_size_set()==false || tab.is_xy_set()==false) {
        O2SCL_ERR2("Grid not set in tensor_grid::",
                   "copy_table3d_interp_value().",o2scl::exc_einval);
      }

      // Get current table3d grid
      size_t nx, ny;
      tab.get_size(nx,ny);

      // Create slice if not already present
      size_t is;
      if (!tab.is_slice(slice_name,is)) tab.new_slice(slice_name);

      // Loop through the table grid to perform the interpolation
      for(size_t i=0;i<nx;i++) {
        for(size_t j=0;j<ny;j++) {
          values[ix_x]=tab.get_grid_x(i);
          values[ix_y]=tab.get_grid_y(j);
          tab.set(i,j,slice_name,this->interp_linear(values));
          if (verbose>0) {
            std::cout << "At location values: ";
            for(size_t k=0;k<values.size();k++) {
              std::cout << values[k] << " ";
            }
            std::cout << "Interpolated to get: "
                      << i << " " << j << " " << slice_name << " "
                      << this->interp_linear(values) << std::endl;
            if (verbose>1) {
              char ch;
              std::cin >> ch;
            }
          }
        }
      }

      return;
    }

    /** \brief Copy to a slice in a table3d object using interpolation
        creating a new table3d grid
    */
    template<class vec2_t> 
      void copy_table3d_interp_values_setxy
      (size_t ix_x, size_t ix_y, vec2_t &values, table3d &tab,
       std::string x_name="x", std::string y_name="y",
       std::string slice_name="z") const {

      // Get current table3d grid
      size_t nx, ny;
      tab.get_size(nx,ny);

      if (nx==0 && ny==0) {

        if (x_name.length()==0) x_name="x";
        if (y_name.length()==0) y_name="y";
        
        // If there's no grid, then create a grid in the table3d
        // object that is the same as that in the tensor_grid object
        std::vector<double> grid_x, grid_y;
        copy_grid(ix_x,grid_x);
        copy_grid(ix_y,grid_y);
        tab.set_xy("x",grid_x.size(),grid_x,
                   "y",grid_y.size(),grid_y);
        // Now that the grid is set, get nx and ny
        tab.get_size(nx,ny);
      }
      
      copy_table3d_interp_values(ix_x,ix_y,values,tab,slice_name);
      return;
    }
    //@}

    /** \brief Create from a table3d object
     */
    void from_table3d_fermi
      (const table3d &t3d, std::string slice, size_t n_points,
       double low=0.0, double high=0.0, double width=0.0) {
      
      if (n_points<2) {
        O2SCL_ERR("Number of points too small.",o2scl::exc_efailed);
      }

      const ubmatrix &sl=t3d.get_slice(slice);

      if (low>=high) {
        low=matrix_min_value<ubmatrix,double>(sl);
        high=matrix_max_value<ubmatrix,double>(sl);
      }
      uniform_grid_end<double> ug(low,high,n_points-1);

      const ubvector &gx=t3d.get_x_data();
      const ubvector &gy=t3d.get_y_data();
      std::vector<size_t> sz={gx.size(),gy.size(),ug.get_npoints()};
      resize(3,sz);

      if (width<=0.0) {
        width=fabs(high-low)/100.0;
      }
      
      std::vector<double> tgrid;
      for(size_t i=0;i<gx.size();i++) {
        tgrid.push_back(gx[i]);
      }
      for(size_t i=0;i<gy.size();i++) {
        tgrid.push_back(gy[i]);
      }
      for(size_t i=0;i<ug.get_npoints();i++) {
        tgrid.push_back(ug[i]);
      }
      this->set_grid_packed(tgrid);

      for(size_t i=0;i<gx.size();i++) {
        for(size_t j=0;j<gy.size();j++) {
          double val=t3d.get(i,j,slice);
          for(size_t k=0;k<ug.get_npoints();k++) {
            double vk=ug[k];
            double v2=1.0/(1.0+exp((vk-val)/width));
            std::vector<size_t> ix={i,j,k};
            this->set(ix,v2);
          }
        }
      }
      
      return;
    }
    

    
    /// \name Clear method
    //@{
    /// Clear the tensor of all data and free allocated memory
    void clear() {
      grid.resize(0);
      grid_set=false;
      tensor<double,vec_t,vec_size_t>::clear();
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

    /** \brief Obtain a value by looking up some indices and
        interpolating the others

        To call this function, the arguments should be
        of the following form
        - The vector \c ix_to_interp should be a list of indices to
        interpolate. The size of \c ix_to_interp must be at least 1
        or larger but smaller than or equal to the full tensor rank.
        All entries in \c ix_to_interp should be smaller than the
        full tensor rank.
        - The vector \c ix should have a size equal to the tensor 
        rank, but values stored in entries corresponding to the
        indices in \c ix_to_interp will be ignored
        - The vector \c val should be a list of values to be
        interpolated and should have a size equal to that of 
        \c ix_to_interp .

        \verbatim embed:rst
        .. todo:: 

           In tensor_grid::interp_linear_partial(): Double check and
           document if the vector "ix_to_interp" needs to be ordered.
           I'm pretty sure it doesn't, so long as the ordering in \c
           val and \c ix_to_interp are consistent.

        \endverbatim
    */
    template<class vec2_size_t, class vec3_size_t, class vec2_t>
      double interp_linear_partial
      (const vec2_size_t &ix_to_interp,
       vec3_size_t &ix, const vec2_t &val) const {

      if (val.size()!=ix_to_interp.size()) {
        O2SCL_ERR2("Index and value list don't match in ",
                   "tensor_grid::interp_linear_partial().",
                   o2scl::exc_einval);
      }
      if (ix_to_interp.size()>this->get_rank() ||
          ix_to_interp.size()==0) {
        O2SCL_ERR2("Index list too large or too small in ",
                   "tensor_grid::interp_linear_partial().",
                   o2scl::exc_einval);
      }

      // Find the the corner of the hypercube containing val for all
      // the indices to be interpolated
      std::vector<size_t> loc(ix_to_interp.size());
      std::vector<double> gnew;
      for(size_t i=0;i<ix_to_interp.size();i++) {
        size_t ixi=ix_to_interp[i];
        if (ixi>=this->get_rank()) {
          O2SCL_ERR2("Index to interpolate larger than tensor rank in ",
                     "tensor_grid::interp_linear_partial().",
                     o2scl::exc_einval);
        }
        std::vector<double> grid_one(this->size[ixi]);
        for(size_t j=0;j<this->size[ixi];j++) {
          grid_one[j]=this->get_grid(ixi,j);
        }
        search_vec<std::vector<double> > sv(this->size[ixi],grid_one);
        loc[i]=sv.find(val[i]);
        gnew.push_back(grid_one[loc[i]]);
        gnew.push_back(grid_one[loc[i]+1]);
      }

      // Now construct a 2^{rk}-sized tensor containing only that 
      // hypercube
      std::vector<size_t> snew(ix_to_interp.size());
      for(size_t i=0;i<ix_to_interp.size();i++) {
        snew[i]=2;
      }
      tensor_grid tnew(ix_to_interp.size(),snew);
      tnew.set_grid_packed(gnew);
      
      // Copy over the relevant data
      for(size_t i=0;i<tnew.total_size();i++) {
        std::vector<size_t> index_new(ix_to_interp.size());
        tnew.unpack_index(i,index_new);
        for(size_t j=0;j<ix_to_interp.size();j++) {
          ix[ix_to_interp[j]]=index_new[j]+loc[j];
        }
        tnew.set(index_new,this->get(ix));
      }
      
      // Now use interp_power_two()
      return tnew.interp_linear_power_two(val);
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

        This function calls the error handler if the user
        tries to interpolate an empty tensor.

        \future This starts with a small copy, which can be eliminated
        by creating a new version of interp_linear_power_two
        which accepts an offset vector parameter so that the 
        first interpolation is offset. Remaining interpolations
        don't need to be offset because the tensor has to be
        created from the previous interpolation round.
    */
    template<class vec2_t> double interp_linear(vec2_t &v) const {

      if (this->rk==0) {
        O2SCL_ERR2("Tried to interpolate in empty tensor in ",
                   "tensor_grid::interp_linear().",o2scl::exc_einval);
      }
      
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
        tnew.unpack_index(i,index_new);
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
    template<class vec2_t>
      double interp_linear_power_two(vec2_t &v) const {

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
        tnew.unpack_index(i,index);
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
    template<class vec2_t, class vec3_t>
      void interp_linear_vec0(vec2_t &v, vec3_t &res) const {

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
        tnew.unpack_index(i,index_new);
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
    template<class vec2_t, class vec3_t>
      void interp_linear_power_two_vec0(vec2_t &v, vec3_t &res) const {
      
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
        tnew.unpack_index(i,index);
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
    template<class vec2_t, class vec3_t>
      void interp_linear_vec(vec2_t &v, size_t ifree, vec3_t &res) const {

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
        tnew.unpack_index(i,index_new);
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
  
  /** \brief Rearrange, sum and copy current tensor to a new tensor
      
      \note Experimental

      \future Some code duplication between this function
      and the one in the tensor class.
  */
  template<class tensor_t, class data_t>
  tensor_t grid_rearrange_and_copy(const tensor_t &t,
                                    std::vector<index_spec> spec,
                               int verbose=0, bool err_on_fail=true) {
    
    // Old rank and new rank (computed later)
    size_t rank_old=t.get_rank();
    size_t rank_new=0;
      
    // Size of the new indices
    std::vector<size_t> size_new;

    // Reorganize the index specifications by both the old
    // and new indexing system
    std::vector<index_spec> spec_old(rank_old);
    std::vector<index_spec> spec_new;
  
    // Number of elements to sum over
    size_t n_sum_loop=1;
    // Number of interpolations to perform
    size_t n_interps=0;

    // Size of sums
    std::vector<size_t> sum_sizes;

    // List of indexes to interpolate
    std::vector<size_t> ix_to_interp;

    // New grid
    std::vector<double> new_grid;
  
    // Loop through the index specifications and add them to
    // spec_old and spec_new (if necessary). This loop also
    // determines the rank of the new tensor, "rank_new", and the
    // sizes of the indices "size_new". Finally, this loop counts
    // the number of interpolations and modifies n_sum_loop as
    // needed. No actual copying or summing is done in this loop
    // yet.
      
    for(size_t i=0;i<spec.size();i++) {
      if (spec[i].type==index_spec::index ||
          spec[i].type==index_spec::reverse) {
        if (spec[i].ix1>=rank_old) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (index,reverse) in ",
                       "tensor_grid::rearrange_and_copy().",
                       o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index too large (index,reverse) in "
                        << "tensor_grid::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_grid<>();
          }
        }
        size_new.push_back(t.get_size(spec[i].ix1));
        // Use ix1 to store the destination index (which is
        // at this point equal to rank_new)
        spec_old[spec[i].ix1]=index_spec(spec[i].type,
                                         rank_new,
                                         spec[i].ix2,0,
                                         spec[i].val1);
        spec_new.push_back(index_spec(spec[i].type,
                                      spec[i].ix1,
                                      spec[i].ix2,0,
                                      spec[i].val1));
        rank_new++;
        // Update the new grid
        if (spec[i].type==index_spec::index) {
          for(size_t k=0;k<t.get_size(spec[i].ix1);k++) {
            new_grid.push_back(t.get_grid(spec[i].ix1,k));
          }
        } else {
          size_t nt=t.get_size(spec[i].ix1);
          for(size_t k=0;k<nt;k++) {
            new_grid.push_back(t.get_grid(spec[i].ix1,nt-1-k));
          }
        }
      } else if (spec[i].type==index_spec::range) {
        if (spec[i].ix1>=rank_old ||
            spec[i].ix2>=t.get_size(spec[i].ix1) ||
            spec[i].ix3>=t.get_size(spec[i].ix1)) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (range) in ",
                       "tensor::rearrange_and_copy().",o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index too large (range) in "
                        << "tensor in tensor::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_grid<>();
          }
        }
        if (spec[i].ix3>spec[i].ix2) {
          size_new.push_back(spec[i].ix3-spec[i].ix2+1);
        } else {
          size_new.push_back(spec[i].ix2-spec[i].ix3+1);
        }
        // Use ix1 to store the destination index (which is
        // at this point equal to rank_new)
        spec_old[spec[i].ix1]=
          index_spec(spec[i].type,rank_new,spec[i].ix2,
                     spec[i].ix3,spec[i].val1);
        spec_new.push_back
          (index_spec(spec[i].type,spec[i].ix1,
                      spec[i].ix2,spec[i].ix3,spec[i].val1));
        rank_new++;
        // Update the new grid
        if (spec[i].ix3>spec[i].ix2) {
          for(size_t k=0;k<spec[i].ix3-spec[i].ix2+1;k++) {
            new_grid.push_back(t.get_grid(spec[i].ix1,k+spec[i].ix2));
          }
        } else {
          for(size_t k=0;k<spec[i].ix2-spec[i].ix3+1;k++) {
            new_grid.push_back(t.get_grid(spec[i].ix1,k+spec[i].ix3));
          }
        }
      } else if (spec[i].type==index_spec::trace) {
        if (spec[i].ix1>=rank_old || spec[i].ix2>=rank_old) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (trace) in ",
                       "tensor_grid::rearrange_and_copy().",
                       o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index too large (trace) in "
                        << "tensor_grid::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_grid<>();
          }
        }
        if (t.get_size(spec[i].ix1)<t.get_size(spec[i].ix2)) {
          n_sum_loop*=t.get_size(spec[i].ix1);
          sum_sizes.push_back(t.get_size(spec[i].ix1));
        } else {
          n_sum_loop*=t.get_size(spec[i].ix2);
          sum_sizes.push_back(t.get_size(spec[i].ix2));
        }
        // We set the values of ix1 and ix2 so that ix2
        // always refers to the other index being traced over
        spec_old[spec[i].ix1]=index_spec(spec[i].type,spec[i].ix1,
                                         spec[i].ix2);
        spec_old[spec[i].ix2]=index_spec(spec[i].type,spec[i].ix2,
                                         spec[i].ix1);
      } else if (spec[i].type==index_spec::sum) {
        if (spec[i].ix1>=rank_old) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (sum) in ",
                       "tensor_grid::rearrange_and_copy().",
                       o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index " << spec[i].ix1
                        << " too large (sum) in "
                        << "tensor_grid::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_grid<>();
          }
        }
        n_sum_loop*=t.get_size(spec[i].ix1);
        sum_sizes.push_back(t.get_size(spec[i].ix1));
        spec_old[spec[i].ix1]=index_spec(spec[i].type,
                                         spec[i].ix1,
                                         spec[i].ix2,0,
                                         spec[i].val1);
      } else if (spec[i].type==index_spec::fixed) {
        if (spec[i].ix1>=rank_old ||
            spec[i].ix2>=t.get_size(spec[i].ix1)) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (fixed) in ",
                       "tensor_grid::rearrange_and_copy().",
                       o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index " << spec[i].ix1 << " or "
                        << spec[i].ix2 << " too large (fixed) in "
                        << "tensor_grid::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_grid<>();
          }
        }
        // Use ix1 to store the destination index (which is
        // at this point equal to rank_new)
        spec_old[spec[i].ix1]=index_spec(spec[i].type,rank_new,
                                         spec[i].ix2);
      } else if (spec[i].type==index_spec::interp) {
        if (spec[i].ix1>=rank_old) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (interp) in ",
                       "tensor_grid::rearrange_and_copy().",
                       o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index " << spec[i].ix1
                        << " too large (interp) in "
                        << "tensor_grid::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_grid<>();
          }
        }
        // Use ix1 to store the destination index (which is
        // at this point equal to rank_new)
        spec_old[spec[i].ix1]=index_spec(spec[i].type,
                                         rank_new,
                                         spec[i].ix2,0,
                                         spec[i].val1);
        n_interps++;
        ix_to_interp.push_back(spec[i].ix1);
      } else if (spec[i].type==index_spec::grid ||
                 spec[i].type==index_spec::gridw) {
        if (spec[i].ix1>=rank_old) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (grid) in ",
                       "tensor_grid::rearrange_and_copy().",
                       o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index " << spec[i].ix1
                        << " too large (grid) in "
                        << "tensor_grid::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_grid<>();
          }
        }

        // Setup new grid in log or linear mode and determine
        // spacing between grid points for later use
        if (spec[i].ix3==1) {
          // Log mode
          double rat;
          int curr_grid_size;
          if (spec[i].type==index_spec::grid) {
            rat=pow(spec[i].val2/spec[i].val1,1.0/((double)spec[i].ix2));
            spec[i].val3=rat;
            curr_grid_size=spec[i].ix2+1;
          } else {
            rat=spec[i].val3;
            curr_grid_size=((size_t)(log(spec[i].val2/spec[i].val1)/
                                     log(spec[i].val3)));
          }
          size_new.push_back(curr_grid_size);
          for(int k=0;k<curr_grid_size;k++) {
            if (k==0) {
              new_grid.push_back(spec[i].val1);
            } else if (curr_grid_size>1 && k==curr_grid_size-1) {
              new_grid.push_back(spec[i].val2);
            } else {
              new_grid.push_back(spec[i].val1*pow(rat,((double)k)));
            }
          }
        } else {
          // Linear mode
          double width;
          int curr_grid_size;
          if (spec[i].type==index_spec::grid) {
            width=(spec[i].val2-spec[i].val1)/((double)spec[i].ix2);
            spec[i].val3=width;
            curr_grid_size=spec[i].ix2+1;
          } else {
            width=spec[i].val3;
            curr_grid_size=((size_t)(1+(spec[i].val2-
                                        spec[i].val1)/spec[i].val3));
          }
          size_new.push_back(curr_grid_size);
          for(int k=0;k<curr_grid_size;k++) {
            if (k==0) {
              new_grid.push_back(spec[i].val1);
            } else if (curr_grid_size>1 && k==curr_grid_size-1) {
              new_grid.push_back(spec[i].val2);
            } else {
              new_grid.push_back(spec[i].val1+width*((double)k));
            }
          }
        }
          
        // Use ix1 to store the destination index (which is
        // at this point equal to rank_new)
        spec_old[spec[i].ix1]=index_spec
          (spec[i].type,rank_new,spec[i].ix2,spec[i].ix3,
           spec[i].val1,spec[i].val2,spec[i].val3);
        spec_new.push_back(index_spec
                           (spec[i].type,spec[i].ix1,spec[i].ix2,
                            spec[i].ix3,spec[i].val1,
                            spec[i].val2,spec[i].val3));
          
        n_interps++;
        ix_to_interp.push_back(spec[i].ix1);
          
        rank_new++;
          
      } else {
        if (err_on_fail) {
          O2SCL_ERR2("Index specification type not allowed in ",
                     "tensor_grid::rearrange_and_copy()",
                     o2scl::exc_einval);
        } else {
          if (verbose>0) {
            std::cout << "Index specification type not allowed in "
                      << "tensor_grid::rearrange_and_copy()." << std::endl;
          }
          return tensor_grid<>();
        }
      }
        
      // End of loop over index specifications
    }

    // Compute the total number of sums to be performed
    size_t n_sums=sum_sizes.size();

    // If the new rank is zero, then return an error
    if (rank_new==0) {
      if (err_on_fail) {
        O2SCL_ERR2("Zero new indices in ",
                   "tensor_grid::rearrange_and_copy()",
                   o2scl::exc_einval);
      } else {
        if (verbose>0) {
          std::cout << "Zero new indices in " 
                    << "tensor_grid::rearrange_and_copy()." << std::endl;
        }
        return tensor_grid<>();
      }
    }

    // If not all indices in the old tensor were accounted for,
    // then return an error
    for(size_t i=0;i<rank_old;i++) {
      if (spec_old[i].type==index_spec::empty) {
        if (err_on_fail) {
          O2SCL_ERR2("Not all indices accounted for in ",
                     "tensor_grid::rearrange_and_copy()",
                     o2scl::exc_einval);
        } else {
          if (verbose>0) {
            std::cout << "Index " << i << " not accounted for in "
                      << "tensor_grid::rearrange_and_copy()." << std::endl;
          }
          return tensor_grid<>();
        }
      }
    }

    // Verbose output of change in rank
    if (verbose>0) {
      std::cout << "Using a " << rank_old
                << " rank tensor to create a new "
                << rank_new << " rank tensor." << std::endl;
    }

    // Verbose output looping over the spec_old list
    if (verbose>1) {
      for(size_t i=0;i<rank_old;i++) {
        std::cout << "Old index " << i;
        if (spec_old[i].type==index_spec::index) {
          std::cout << " is being remapped to new index "
                    << spec_old[i].ix1 << "." << std::endl;
        } else if (spec_old[i].type==index_spec::range) {
          std::cout << " is being remapped to new index "
                    << spec_old[i].ix1 << " with a range from "
                    << spec_old[i].ix2 << " to " << spec_old[i].ix3
                    << "." << std::endl;
        } else if (spec_old[i].type==index_spec::reverse) {
          std::cout << " is being reversed and remapped to new index "
                    << spec_old[i].ix1 << "." << std::endl;
        } else if (spec_old[i].type==index_spec::trace) {
          std::cout << " is being traced with index "
                    << spec_old[i].ix2 << "." << std::endl;
        } else if (spec_old[i].type==index_spec::sum) {
          std::cout << " is being summed." << std::endl;
        } else if (spec_old[i].type==index_spec::fixed) {
          std::cout << " is being fixed to " << spec_old[i].ix2
                    << "." << std::endl;
        } else if (spec_old[i].type==index_spec::interp) {
          std::cout << " is being interpolated from value "
                    << spec_old[i].val1 << "." << std::endl;
        } else if (spec_old[i].type==index_spec::grid ||
                   spec_old[i].type==index_spec::gridw) {
          std::cout << " is being reinterpolated based on grid\n  "
                    << spec_old[i].val1 << " "
                    << spec_old[i].val2 << " ";
          if (spec_old[i].type==index_spec::grid) {
            std::cout << spec_old[i].ix2;
          } else {
            std::cout << spec_old[i].val3;
          }
          if (spec_old[i].ix3==1) {
            std::cout << " (log)." << std::endl;
          } else {
            std::cout << "." << std::endl;
          }
          std::cout << "  and placed in new index " << spec_old[i].ix1
                    << "." << std::endl;
        }
      }

      // Verbose output looping over the spec_new list
      for(size_t i=0;i<rank_new;i++) {
        std::cout << "New index " << i;
        if (spec_new[i].type==index_spec::index) {
          std::cout << " was remapped from old index " << spec_new[i].ix1
                    << "." << std::endl;
        } else if (spec_new[i].type==index_spec::range) {
          std::cout << " was remapped from old index " << spec_new[i].ix1
                    << " using range from " << spec_new[i].ix2 << " to "
                    << spec_new[i].ix3 << "." << std::endl;
        } else if (spec_new[i].type==index_spec::reverse) {
          std::cout << " was reversed and remapped from old index "
                    << spec_new[i].ix1 << "." << std::endl;
        } else if (spec_new[i].type==index_spec::grid ||
                   spec_new[i].type==index_spec::gridw) {
          std::cout << " was obtained from grid\n  "
                    << spec_new[i].val1 << " "
                    << spec_new[i].val2 << " ";
          if (spec_new[i].type==index_spec::grid) {
            std::cout << spec_new[i].ix2;
          } else {
            std::cout << spec_new[i].val3;
          }
          if (spec_new[i].ix3==1) {
            std::cout << " (log)." << std::endl;
          } else {
            std::cout << "." << std::endl;
          }
          std::cout << "  originally stored in index "
                    << spec_new[i].ix1 << "." << std::endl;
        }
      }

      // Output sum information
      std::cout << "n_sums, n_sum_loop: " << n_sums << " "
                << n_sum_loop << std::endl;
        
      // End of 'if (verbose>1)'
    }
    
    // Create the new tensor object
    tensor_grid<> t_new(rank_new,size_new);
    t_new.set_grid_packed(new_grid);
    if (verbose>1) {
      std::cout << "New grid is: " << std::endl;
      for(size_t k=0;k<rank_new;k++) {
        std::cout << k << " (" << t_new.get_size(k) << "): ";
        if (t_new.get_size(k)>3) {
          std::cout << t_new.get_grid(k,0) << ", ";
          std::cout << t_new.get_grid(k,1) << " ... ";
        } else {
          for(int ell=0;ell<((int)t_new.get_size(k))-1;ell++) {
            std::cout << t_new.get_grid(k,ell) << ", ";
          }
        }
        std::cout << t_new.get_grid(k,t_new.get_size(k)-1) << std::endl;
      }
      std::cout << "Number of interpolations: " << n_interps
                << ", indexes to interpolate: ";
      o2scl::vector_out(std::cout,ix_to_interp,true);
    }
    
    // Index arrays. For indices in the old tensor which we are
    // interpolating, the value of ix_old is not used, so it
    // is not set.
    std::vector<size_t> ix_new(rank_new);
    std::vector<size_t> ix_old(rank_old);
    std::vector<size_t> sum_ix(n_sums);

    // The vector "interp_vals" below is ordered according to the
    // indices in the old vector, but ix_to_interp isn't always
    // ordered that way, so we sort ix_to_interp here. This sorting
    // is important for the call interp_linear_partial() below.
    //o2scl::vector_sort<std::vector<size_t>,size_t>(ix_to_interp.size(),
    //ix_to_interp);
      
    // Loop over the new tensor object
    for(size_t i=0;i<t_new.total_size();i++) {

      // Find the location in the new tensor object
      t_new.unpack_index(i,ix_new);

      // List of interpolated values (vector of size n_interps)
      std::vector<double> interp_vals;

      // Determine the location in the old tensor object
      for(size_t j=0;j<rank_old;j++) {
          
        if (spec_old[j].type==index_spec::index) {
          ix_old[j]=ix_new[spec_old[j].ix1];
        } else if (spec_old[j].type==index_spec::range) {
          if (spec_old[j].ix2<spec_old[j].ix3) {
            ix_old[j]=ix_new[spec_old[j].ix1]+spec_old[j].ix2;
          } else {
            ix_old[j]=spec_old[j].ix2-ix_new[spec_old[j].ix1];
          }
        } else if (spec_old[j].type==index_spec::reverse) {
          ix_old[j]=t.get_size(j)-1-ix_new[spec_old[j].ix1];
        } else if (spec_old[j].type==index_spec::fixed) {
          ix_old[j]=spec_old[j].ix2;
        } else if (spec_old[j].type==index_spec::interp) {
          interp_vals.push_back(spec_old[j].val1);
        } else if (spec_old[j].type==index_spec::grid) {
          if (spec_old[j].ix3==1) {
            double val=spec_old[j].val1*
              pow(spec_old[j].val3,ix_new[spec_old[j].ix1]);
            interp_vals.push_back(val);
          } else {
            double val=spec_old[j].val1+
              ix_new[spec_old[j].ix1]*spec_old[j].val3;
            interp_vals.push_back(val);
          }
        } else if (spec_old[j].type==index_spec::gridw) {
          if (spec_old[j].ix3==1) {
            double val=spec_old[j].val1*
              pow(spec_old[j].val3,ix_new[spec_old[j].ix1]);
            interp_vals.push_back(val);
          } else {
            double val=spec_old[j].val1+
              ix_new[spec_old[j].ix1]*spec_old[j].val3;
            interp_vals.push_back(val);
          }
        }

      }
        
      size_t ntmp=t_new.total_size()/10;
      if (ntmp==0) ntmp++;

      if (verbose>1 && i%ntmp==0) {
        for(size_t j2=0;j2<rank_old;j2++) {
          size_t k;
          if (vector_search(ix_to_interp,j2,k)==false) {
            if (j2>=ix_old.size()) {
              std::cout << "\n\nj2,ix_old.size(): "
                        << j2 << " " << ix_old.size() << std::endl;
              O2SCL_ERR2("Indexing problem 1 in",
                         "tensor_grid::rearrange_and_copy().",
                         o2scl::exc_esanity);
            } else {
              std::cout << ix_old[j2] << " ";
            }
          } else {
            if (k>=interp_vals.size()) {
              std::cout << "\n\nk,interp_vals.size(): "
                        << k << " " << interp_vals.size()
                        << std::endl;
              o2scl::vector_out(std::cout,interp_vals,true);
              O2SCL_ERR2("Indexing problem w in",
                         "tensor_grid::rearrange_and_copy().",
                         o2scl::exc_esanity);
            } else {
              std::cout << "(" << interp_vals[k] << ") ";
            }
          }
        }
        std::cout << " ix_new: ";
        vector_out(std::cout,ix_new,true);
      }
        
      double val=0;

      for(size_t j=0;j<n_sum_loop;j++) {

        // This code is similar to tensor::unpack_index(), it unpacks
        // the index j to the indices which we are summing over.
        size_t j2=j, sub_size;
        for(size_t k=0;k<n_sums;k++) {
          if (k==n_sums-1) {
            sum_ix[k]=j2;
          } else {
            sub_size=1;
            for(size_t kk=k+1;kk<n_sums;kk++) sub_size*=sum_sizes[kk];
            sum_ix[k]=j2/sub_size;
            j2-=sub_size*(j2/sub_size);
          }
        }
          
        if (verbose>2) {
          std::cout << "n_sum_loop: " << n_sum_loop << " n_sums: "
                    << n_sums << " sum_sizes: ";
          vector_out(std::cout,sum_sizes,true);
          std::cout << "j: " << j << " sum_ix: ";
          vector_out(std::cout,sum_ix,true);
        }

        // Remap from sum_ix to ix_old
        size_t cnt=0;
        for(size_t k=0;k<rank_old;k++) {
          if (spec_old[k].type==index_spec::sum) {
            if (cnt>=sum_ix.size()) {
              std::cout << "X: " << cnt << " " << sum_ix.size()
                        << std::endl;
              O2SCL_ERR2("Bad sync 1 in sum_ix in ",
                         "tensor_grid::rearrange_and_copy()",
                         o2scl::exc_esanity);
            }
            ix_old[k]=sum_ix[cnt];
            cnt++;
          } else if (spec_old[k].type==index_spec::trace &&
                     spec_old[k].ix1<spec_old[k].ix2) {
            if (cnt>=sum_ix.size()) {
              std::cout << "X: " << cnt << " " << sum_ix.size()
                        << std::endl;
              O2SCL_ERR2("Bad sync 2 in sum_ix in ",
                         "tensor_grid::rearrange_and_copy()",
                         o2scl::exc_esanity);
            }
            ix_old[spec_old[k].ix1]=sum_ix[cnt];
            ix_old[spec_old[k].ix2]=sum_ix[cnt];
            cnt++;
          }
        }

        if (verbose>2) {
          std::cout << "Old: ";
          for(size_t j3=0;j3<rank_old;j3++) {
            std::cout << ix_old[j3] << " ";
          }
          std::cout << "New: ";
          vector_out(std::cout,ix_new,true);
        }
        if (n_interps>0) {
          val+=t.interp_linear_partial(ix_to_interp,ix_old,interp_vals);
        } else {
          val+=t.get(ix_old);
        }
      
      }
      
      // Set the new point by performing the linear interpolation
      t_new.set(ix_new,val);
    }
    
    return t_new;
  }

  /** \brief Rearrange, sum and copy current tensor to a new tensor
      (string input version)
  */
  template<class tensor_t, class data_t>
  tensor_t grid_rearrange_and_copy(const tensor_t &t, std::string spec,
                               int verbose=0, bool err_on_fail=true) {
    
    std::vector<std::string> sv2;
    index_spec_preprocess2(spec,sv2);
    std::vector<o2scl::index_spec> vis;
    strings_to_indexes2(sv2,vis,verbose);
    return grid_rearrange_and_copy<tensor_t,data_t>(t,vis,verbose,err_on_fail);
  }
  
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

  /** \brief Desc
   */
  template<class tensor_grid_t>
  void tensor_grid_out(std::ostream &os, tensor_grid_t &t,
                       bool pretty=true) {

    if (pretty) {
      
      size_t rk=t.get_rank();
      os << "rank: " << rk << " sizes: ";
      for(size_t i=0;i<rk;i++) {
	os << t.get_size(i) << " ";
      }
      os << std::endl;
      os << "grids:" << std::endl;
      auto &grid=t.get_grid();
      size_t jtmp=0;
      for(size_t i=0;i<rk;i++) {
	os << i << ": ";
        for(size_t j=0;j<t.get_size(i);j++) {
          os << grid[jtmp];
          if (j!=t.get_size(i)-1) os << " ";
          jtmp++;
        }
        os << std::endl;
      }
      auto &data=t.get_data();
      std::vector<size_t> ix(rk);
      std::vector<std::string> sv, sv_out;
      for(size_t i=0;i<t.total_size();i++) {
	t.unpack_index(i,ix);
	std::string tmp="(";
	for(size_t j=0;j<rk;j++) {
	  if (j!=rk-1) {
	    tmp+=o2scl::szttos(ix[j])+",";
	  } else {
	    tmp+=o2scl::szttos(ix[j]);
	  }
	}
	tmp+="): "+o2scl::dtos(data[i]);
	sv.push_back(tmp);
      }
      screenify(sv.size(),sv,sv_out);
      for(size_t k=0;k<sv_out.size();k++) {
	os << sv_out[k] << std::endl;
      }

    } else {
      
      size_t rk=t.get_rank();
      os << rk << " ";
      for(size_t i=0;i<rk;i++) {
	os << t.get_size(i) << " ";
      }
      os << std::endl;
      for(size_t i=0;i<rk;i++) {
        auto &grid=t.get_grid();
        vector_out(os,grid,true);
      }
      auto &data=t.get_data();
      for(size_t i=0;i<t.total_size();i++) {
	os << data[i] << " ";
	if (i%10==9) os << std::endl;
      }
      os << std::endl;
      
    }
    
    return;
  }

  
  
}

#endif




