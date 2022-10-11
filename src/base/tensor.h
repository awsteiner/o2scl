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
#ifndef O2SCL_TENSOR_H
#define O2SCL_TENSOR_H

/** \file tensor.h
    \brief File defining \ref o2scl::tensor and rank-specific children
*/
#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_ieee_utils.h>

#include <o2scl/err_hnd.h>
#include <o2scl/interp.h>
#include <o2scl/table3d.h>
#include <o2scl/misc.h>

namespace o2scl {
  
  /** \brief Index specification 

      \comment

      It would be great to add derivatives to the index specifications
      and the rearrange function, but I'm not sure how best to do it
      yet because there are a large number of ways to compute the
      derivative. The interpolation is easy because it's just done
      linearly, but using linear interpolation to do the derivative
      isn't as useful. For now, 'acol -deriv' is a way to compute
      derivatives.

      \endcomment
   */
  class index_spec {
    
  public:
  
    /// Type of specification
    size_t type;
    /// First argument
    size_t ix1;
    /// Second argument
    size_t ix2;
    /// Third argument
    size_t ix3;
    /// First double argument
    double val1;
    /// Second double argument
    double val2;
    /// Third double argument
    double val3;

    /// \name Possible values for type
    //@{
    /// Empty specification
    static const size_t empty=0;
    /// Retain an index
    static const size_t index=1;
    /// Fix the value of an index
    static const size_t fixed=2;
    /// Sum over an index
    static const size_t sum=3;
    /// Perform a trace (sum over two indices)
    static const size_t trace=4;
    /// Reverse an index
    static const size_t reverse=5;
    /// Choose a new range for an index
    static const size_t range=6;
    /// Interpolate a value to fix an index
    static const size_t interp=7;
    /// Interpolate a value to set a new grid (fixed bin number)
    static const size_t grid=8;
    /// Interpolate a value to set a new grid (fixed bin width)
    static const size_t gridw=9;
    /// Interpolate a value to fix an index and take the derivative
    static const size_t deriv=10;
    //@}

    /// Default constructor
    index_spec() {
      type=0;
      ix1=0;
      ix2=0;
      ix3=0;
      val1=0.0;
      val2=0.0;
      val3=0.0;
    }

    /// Explicit full constructor
    index_spec(size_t typ, size_t i1, size_t i2=0, size_t i3=0,
               double v1=0.0, double v2=0.0, double v3=0.0) {
      type=typ;
      ix1=i1;
      ix2=i2;
      ix3=i3;
      val1=v1;
      val2=v2;
      val3=v3;
    }

    /// If true, then two index specifications are equal
    bool equal(index_spec &is) {
      bool ret;
      ret=(type==is.type && ix1==is.ix1 && ix2==is.ix2 && ix3==is.ix3 &&
           val1==is.val1 && val2==is.val2 && val3==is.val3);
      if (ret==false) {
        std::cout << type << " " << ix1 << " " << ix2 << " "
                  << ix3 << " " << val1 << " " << val2 << " " << val3
                  << std::endl;
        std::cout << is.type << " " << is.ix1 << " " << is.ix2 << " "
                  << is.ix3 << " " << is.val1 << " " << is.val2 << " "
                  << is.val3 << std::endl;
      }
      return ret;
    }
    
  };

  /** \brief Unmodified index for tensors
   */
  class ix_index : public index_spec {

  public:

    /// Create an ix_index object from index \c ix
    ix_index(size_t ix);
    
    /// Create an ix_index object from a index_spec object
    ix_index(index_spec &is) {
      if (is.type!=index_spec::index) {
	O2SCL_ERR("Invalid index_spec in ix_index",
		  o2scl::exc_einval);
      }
      this->type=is.type;
      this->ix1=is.ix1;
      this->ix2=0;
      this->ix3=0;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
  };
  
  /** \brief Fix an index for tensors
   */
  class ix_fixed : public index_spec {

  public:

    /// The value at which the specified index is to be fixed
    size_t &fixed_value;

    /// Create an ix_fixed object for index \c ix at value \c fix
    ix_fixed(size_t ix, size_t fix) : fixed_value(this->ix2) {
      this->type=index_spec::fixed;
      this->ix1=ix;
      this->ix2=fix;
      this->ix3=0;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
    /// Create an ix_fixed object from a index_spec object
    ix_fixed(index_spec &is) : fixed_value(this->ix2) {
      if (is.type!=index_spec::fixed) {
	O2SCL_ERR("Invalid index_spec in ix_fixed",
		  o2scl::exc_einval);
      }
      this->type=is.type;
      this->ix1=is.ix1;
      this->ix2=is.ix2;
      this->ix3=0;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
  };
  
  /** \brief Sum over an index for tensors
   */
  class ix_sum : public index_spec {

  public:

    /// Create an ix_sum object for index \c ix 
    ix_sum(size_t ix) {
      this->type=index_spec::sum;
      this->ix1=ix;
      this->ix2=0;
      this->ix3=0;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
    /// Create an ix_sum object from a index_spec object
    ix_sum(index_spec &is) {
      if (is.type!=index_spec::sum) {
	O2SCL_ERR("Invalid index_spec in ix_sum",
		  o2scl::exc_einval);
      }
      this->type=is.type;
      this->ix1=is.ix1;
      this->ix2=0;
      this->ix3=0;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
  };

  /** \brief Perform a trace over two tensor indices
   */
  class ix_trace : public index_spec {

  public:

    /// The second index to trace over
    size_t &second_index;

    /// Create an ix_trace object for indices \c ix and \c jx
    ix_trace(size_t ix, size_t jx) : second_index(this->ix2) {
      this->type=index_spec::trace;
      this->ix1=ix;
      this->ix2=jx;
      this->ix3=0;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
    /// Create an ix_trace object from a index_spec object
    ix_trace(index_spec &is) : second_index(is.ix2) {
      if (is.type!=index_spec::trace) {
	O2SCL_ERR("Invalid index_spec in ix_trace",
		  o2scl::exc_einval);
      }
      this->type=is.type;
      this->ix1=is.ix1;
      this->ix2=is.ix2;
      this->ix3=0;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
  };

  /** \brief Reverse the order of an index
   */
  class ix_reverse : public index_spec {

  public:

    /// Create an ix_reverse object for index \c ix
    ix_reverse(size_t ix) {
      this->type=index_spec::reverse;
      this->ix1=ix;
      this->ix2=0;
      this->ix3=0;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
    /// Create an ix_reverse object from a index_spec object
    ix_reverse(index_spec &is) {
      if (is.type!=index_spec::reverse) {
	O2SCL_ERR("Invalid index_spec in ix_reverse",
		  o2scl::exc_einval);
      }
      this->type=is.type;
      this->ix1=is.ix1;
      this->ix2=0;
      this->ix3=0;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
  };

  /** \brief Restrict an index to a range
   */
  class ix_range : public index_spec {

  public:

    /// The beginning of the range
    size_t &begin;

    /// The end of the range
    size_t &end;

    /** \brief Create an ix_range object for index \c ix beginning at \c 
        start and ending at \c finish
    */
    ix_range(size_t ix, size_t start, size_t finish) :
      begin(this->ix2), end(this->ix3) {
      this->type=index_spec::range;
      this->ix1=ix;
      this->ix2=start;
      this->ix3=finish;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
    /// Create an ix_range object from a index_spec object
    ix_range(index_spec &is) : begin(this->ix2), end(this->ix3) {
      if (is.type!=index_spec::range) {
	O2SCL_ERR("Invalid index_spec in ix_range",
		  o2scl::exc_einval);
      }
      this->type=is.type;
      this->ix1=is.ix1;
      this->ix2=is.ix2;
      this->ix3=is.ix3;
      this->val1=0.0;
      this->val2=0.0;
      this->val3=0.0;
    }
    
  };

  /** \brief For a tensor_grid object, interpolate a value into 
      the grid for one of the indices
  */
  class ix_interp : public index_spec {

  public:

    /// The value to be interpolated into index \c ix
    double &val;

    /// Create an ix_interp object for index \c ix and value \c v
    ix_interp(size_t ix, double v) : val(this->val1) {
      this->type=index_spec::interp;
      this->ix1=ix;
      this->ix2=0;
      this->ix3=0;
      this->val1=v;
      this->val2=0.0;
      this->val3=0.0;
    }
    
    /// Create an ix_interp object from a index_spec object
    ix_interp(index_spec &is) : val(this->val1) {
      if (is.type!=index_spec::interp) {
	O2SCL_ERR("Invalid index_spec in ix_interp",
		  o2scl::exc_einval);
      }
      this->type=is.type;
      this->ix1=is.ix1;
      this->ix2=0;
      this->ix3=0;
      this->val1=is.val1;
      this->val2=0.0;
      this->val3=0.0;
    }
    
  };

  /** \brief For a tensor_grid object, interpolate a grid
  */
  class ix_grid : public index_spec {

  public:

    /// The first grid point
    double &begin;

    /// The last grid point
    double &end;

    /// The number of intervals between grid points
    size_t &n_bins;

    /// True for a logarithmic grid
    size_t &log_flag;
    
    /// Create an ix_grid object from the specified inputs
    ix_grid(size_t ix, double start, double finish, size_t bins,
            bool log=false) : begin(this->val1), end(this->val2),
                        n_bins(this->ix2), log_flag(this->ix3) {
      this->type=index_spec::grid;
      this->ix1=ix;
      this->ix2=bins;
      if (log==true) {
        this->ix3=1;
      } else {
        this->ix3=0;
      }
      this->val1=start;
      this->val2=finish;
      this->val3=0.0;
    }

    /// Create an ix_grid object from a index_spec object
    ix_grid(index_spec &is) : begin(this->val1), end(this->val2),
                              n_bins(this->ix2), log_flag(this->ix3) {
      if (is.type!=index_spec::grid) {
	O2SCL_ERR("Invalid index_spec in ix_grid",
		  o2scl::exc_einval);
      }
      this->type=is.type;
      this->ix1=is.ix1;
      this->ix2=is.ix2;
      this->ix3=is.ix3;
      this->val1=is.val1;
      this->val2=is.val2;
      this->val3=0.0;
    }
    
  };

  /** \brief For a tensor_grid object, interpolate a grid with
      a fixed width
  */
  class ix_gridw : public index_spec {

  public:
    
    /// The first grid point
    double &begin;

    /// The last grid point
    double &end;

    /// The size of the interval between grid points
    double &width;

    /// True for a logarithmic grid
    size_t &log_flag;
    
    /// Create an ix_gridw object from the specified inputs
    ix_gridw(size_t ix, double start, double finish, double wid,
             bool log=false) : begin(this->val1), end(this->val2),
                         width(this->val3), log_flag(this->ix3) {
      this->type=index_spec::gridw;
      this->ix1=ix;
      this->ix2=0;
      if (log==true) {
        this->ix3=1;
      } else {
        this->ix3=0;
      }
      this->val1=start;
      this->val2=finish;
      this->val3=wid;
    }

    /// Create an ix_gridw object from a index_spec object
    ix_gridw(index_spec &is) : begin(this->val1), end(this->val2),
                               width(this->val3), log_flag(this->ix3) {
      if (is.type!=index_spec::gridw) {
	O2SCL_ERR("Invalid index_spec in ix_gridw",
		  o2scl::exc_einval);
      }
      this->type=is.type;
      this->ix1=is.ix1;
      this->ix2=0;
      this->ix3=is.ix3;
      this->val1=is.val1;
      this->val2=is.val2;
      this->val3=is.val3;
    }
    
  };
  
  /** \brief Tensor class with arbitrary dimensions

      The elements of a tensor are typically specified as a list of
      <tt>size_t</tt> numbers with length equal to the tensor rank.
      For a rank-4 tensor named \c t, the element 
      <tt>t[1][2][0][3]</tt> can be obtained with something similar to 
      \code
      size_t ix[4]={1,2,0,3};
      double x=t.get(ix);
      \endcode

      Empty tensors have zero rank.

      The type <tt>vec_t</tt> can be any vector type with
      <tt>operator[]</tt>, <tt>size()</tt> and <tt>resize()</tt>
      methods. The type <tt>vec_size_t</tt> can be any integer-like
      vector type with <tt>operator[]</tt>, <tt>size()</tt> and
      <tt>resize()</tt> methods.

      For I/O with tensors, see \ref o2scl_hdf::hdf_file::setd_ten()
      and \ref o2scl_hdf::hdf_file::getd_ten() . 
      \verbatim embed:rst
      See the the discussion in the sections :ref:`Tensors` and
      :ref:`I/O and contiguous storage` of the User's Guide for more
      details.
      \endverbatim

      The storage pattern is a generalization of row-major order.
      In the case of a 4-rank tensor, the location of a generic 
      element is
      \f[
      \left( \left( i_0 s_1 + i_1 \right) s_2 + i_2 \right) s_3 + i_3 \, .
      \f]
      In this case the distance between two elements \f$(i_0,i_1,
      i_2,i_3)\f$ and \f$ (i_0+1,i_1,i_2,i_3) \f$ is \f$ s_1 s_2 s_3
      \f$, the distance between two elements \f$(i_0,i_1,i_2, i_3)\f$
      and \f$ (i_0,i_1+1,i_2,i_3) \f$ is \f$ s_2 s_3 \f$, and the
      elements \f$(i_0,i_1,i_2,i_3)\f$ and \f$ (i_0,i_1,i_2,i_3+1) \f$
      are adjacent.

      \note Slices of tensors are subsets obtained from fixing the
      index of several dimensions, while letting other dimensions
      vary. For a slice with one dimension not fixed, see \ref
      vector_slice(). The \ref o2scl::tensor::vector_slice() function
      should clearly work for uBlas vectors, and seems to work with
      std::vector objects also, but latter use has not been fully
      tested.
      
      \verbatim embed:rst
      
      .. todo:: 

         In class tensor:
        
         - Future: Create an operator[] for tensor and not just tensor1?

         - Future: Could implement arithmetic operators + and - and some
           different products. 

         - Future: Implement copies to and from vector
           and matrices 

         - Future: Implement tensor contractions, i.e. tensor
           = tensor * tensor 

         - Future: Could be interesting to write an iterator for this class.

      \endverbatim

  */
  template<class data_t=double, class vec_t=std::vector<data_t>, 
           class vec_size_t=std::vector<size_t> > class tensor {

  public:
  
#ifndef DOXYGEN_INTERNAL
  
  protected:
  
    /// The data
    vec_t data;
  
    /// A rank-sized vector of the sizes of each dimension
    vec_size_t size;
  
    /// Rank
    size_t rk;
  
#endif
  
  public:

    /// Create an empty tensor with zero rank
    tensor() {
      rk=0;
    }

    /** \brief Create a tensor of rank \c rank with sizes given in \c dim
	
        The parameter \c dim must be a pointer to a vector of sizes with
        length \c rank. If the user requests any of the sizes to be
        zero, this constructor will call the error handler, create an
        empty tensor, and will allocate no memory.
    */
    template<class size_vec_t> 
    tensor(size_t rank, const size_vec_t &dim) {
      if (rank==0) {
        rk=0;
      } else {
        rk=rank;
        for(size_t i=0;i<rk;i++) {
          if (dim[i]==0) {	  
            rk=0;
            O2SCL_ERR((((std::string)"Requested zero size with non-zero ")+
                       "rank for index "+szttos(i)+
                       " in tensor::tensor(size_t,size_vec_t)").c_str(),
                      exc_einval);
          }
        }
        size.resize(rk);
        size_t tot=1;
        for(size_t i=0;i<rk;i++) {
          size[i]=dim[i];
          tot*=size[i];
        }
        data.resize(tot);
      }
    }

    virtual ~tensor() {
    }

    /// \name Method to check for valid object
    //@{
    /** \brief Check that the \ref o2scl::tensor object is valid
     */
    void is_valid() const {
      if (rk==0) {
        if (data.size()!=0) {
          O2SCL_ERR2("Rank is zero but the data vector has non-zero size ",
                     "in tensor::is_valid().",o2scl::exc_esanity);
        }
      }
    
      if (rk!=size.size()) {
        O2SCL_ERR2("Rank does not match size vector size ",
                   "in tensor::is_valid().",o2scl::exc_esanity);
      }

      if (rk>0) {
        size_t tot=1;
        for(size_t i=0;i<rk;i++) tot*=size[i];
        if (tot==0) {
          O2SCL_ERR2("One entry in the size vector is zero ",
                     "in tensor::is_valid().",o2scl::exc_esanity);
        }
        if (tot!=data.size()) {
          O2SCL_ERR2("Product of size vector entries does not match data ",
                     "vector size in tensor::is_valid().",o2scl::exc_esanity);
        }
      }
    
      return;
    }
    //@}
  
    /// \name Copy constructors
    //@{
    /** \brief Copy using <tt>operator()</tt>
     */
    tensor<data_t,vec_t,vec_size_t>
    (const tensor<data_t,vec_t,vec_size_t> &t) {
      rk=t.rk;
      data=t.data;
      size=t.size;
    }

    /** \brief Copy using <tt>operator=()</tt>
     */
    tensor<data_t,vec_t,vec_size_t> &operator=
    (const tensor<data_t,vec_t,vec_size_t> &t) {
      if (this!=&t) {
        rk=t.rk;
        data=t.data;
        size=t.size;
      }
      return *this;
    }
    //@}
  
    /// \name Clear method
    //@{
    /// Clear the tensor of all data and free allocated memory
    void clear() {
      rk=0;
      data.resize(0);
      size.resize(0);
      return;
    }
    //@}
  
    /// \name Set functions
    //@{
    /// Set the element indexed by \c index to value \c val
    template<class size_vec_t>
    void set(const size_vec_t &index, data_t val) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (rk==0) {
        O2SCL_ERR("Empty tensor in tensor::set().",exc_einval);
      }
      if (index[0]>=size[0]) {
        O2SCL_ERR((((std::string)"Value of index[0]=")+szttos(index[0])+
                   " greater than or equal to size[0]="+
                   szttos(size[0])+" in tensor::set().").c_str(),
                  exc_eindex);
      }
#endif
      size_t ix=index[0];
      for(size_t i=1;i<rk;i++) {
#if O2SCL_NO_RANGE_CHECK
#else
        if (index[i]>=size[i]) {
          O2SCL_ERR((((std::string)"Value of index[")+szttos(i)+"]="+
                     szttos(index[i])+" greater than or equal to size "+
                     szttos(size[i])+" in tensor::set().").c_str(),
                    exc_eindex);
        }
#endif
        ix*=size[i];
        ix+=index[i];
      }
      data[ix]=val;
      return;
    }

    /// Set all elements in a tensor to some fixed value
    void set_all(data_t x) {
      for(size_t i=0;i<total_size();i++) data[i]=x;
      return;
    }
  
    /** \brief Swap the data vector
     */
    void swap_data(vec_t &dat) {
      if (data.size()!=dat.size()) {
        O2SCL_ERR2("Size of new vector does not equal tensor size in ",
                   "tensor::swap_data().",o2scl::exc_einval);
      }
      std::swap(dat,data);
      return;
    }
    //@}

    /// \name Get functions
    //@{
    /// Get the element indexed by \c index
    template<class size_vec_t> data_t &get(const size_vec_t &index) {
#if O2SCL_NO_RANGE_CHECK
#else
      if (rk==0) {
        O2SCL_ERR("Empty tensor in tensor::get().",exc_einval);
      }
      if (index[0]>=size[0]) {
        O2SCL_ERR((((std::string)"Value of index[0]=")+szttos(index[0])+
                   " greater than or equal to size[0]="+
                   szttos(size[0])+" in tensor::get().").c_str(),
                  exc_eindex);
      }
#endif
      size_t ix=index[0];
      for(size_t i=1;i<rk;i++) {
#if O2SCL_NO_RANGE_CHECK
#else
        if (index[i]>=size[i]) {
          O2SCL_ERR((((std::string)"Value of index[")+szttos(i)+"]="+
                     szttos(index[i])+" greater than or equal to size "+
                     szttos(size[i])+" in tensor::get().").c_str(),
                    exc_eindex);
        }
#endif
        ix*=size[i];
        ix+=index[i];
      }
      return data[ix];
    }

    /// Get a const reference to the element indexed by \c index
    template<class size_vec_t> 
    data_t const &get(const size_vec_t &index) const {
#if O2SCL_NO_RANGE_CHECK
#else
      if (rk==0) {
        O2SCL_ERR("Empty tensor in tensor::get() (const).",exc_einval);
      }
      if (index[0]>=size[0]) {
        O2SCL_ERR((((std::string)"Value of index[0]=")+szttos(index[0])+
                   " greater than or equal to size[0]="+
                   szttos(size[0])+" in tensor::get() (const).").c_str(),
                  exc_eindex);
      }
#endif
      size_t ix=index[0];
      for(size_t i=1;i<rk;i++) {
#if O2SCL_NO_RANGE_CHECK
#else
        if (index[i]>=size[i]) {
          O2SCL_ERR((((std::string)"Value of index[")+szttos(i)+"]="+
                     szttos(index[i])+" greater than or equal to size "+
                     szttos(size[i])+" in tensor::get() (const).").c_str(),
                    exc_eindex);
        }
#endif
        ix*=size[i];
        ix+=index[i];
      }
      return data[ix];
    }
    //@}
    
    typedef boost::numeric::ublas::vector_slice<
      boost::numeric::ublas::vector<data_t> > ubvector_slice;
    typedef boost::numeric::ublas::slice slice;

    /// \name Slice function
    //@{
    /** \brief Fix all but one index to create a vector

        This fixes all of the indices to the values given in \c index
        except for the index number \c ix, and returns the
        corresponding vector, whose length is equal to the size of the
        tensor in that index. The value <tt>index[ix]</tt> is ignored.

        For example, for a rank 3 tensor allocated with
        \code
        tensor t;
        size_t dim[3]={3,4,5};
        t.resize(3,dim);
        \endcode
        the following code
        \code
        size_t index[3]={1,0,3};
        ubvector_view v=t.vector_slice(1,index);
        \endcode
        Gives a vector \c v of length 4 which refers to the values
        <tt>t(1,0,3)</tt>, <tt>t(1,1,3)</tt>, <tt>t(1,2,3)</tt>, and
        <tt>t(1,3,3)</tt>.
    */
    template<class size_vec_t> 
    ubvector_slice vector_slice(size_t ix, const size_vec_t &index) {
      if (ix+1>rk) {
        O2SCL_ERR((((std::string)"Specified index ")+szttos(ix)+
                   " greater than or equal to rank "+szttos(rk)+
                   " in tensor::vector_slice()").c_str(),
                  exc_eindex);
      }
      size_t start;
      if (ix==0) start=0;
      else start=index[0];
      for(size_t i=1;i<rk;i++) {
        start*=size[i];
        if (i!=ix) start+=index[i];
      }
      size_t stride=1;
      for(size_t i=ix+1;i<rk;i++) stride*=size[i];
      return ubvector_slice(data,slice(start,stride,size[ix]));
    }
    //@}

#ifdef O2SCL_NEVER_DEFINED

    /** \brief Fix all but two indices to create a matrix

        One could use ublas::make_matrix_from_pointer() and
        then ublas matrix slicing to perform this operation, but
        only if the template types are built on ublas objects.
    */
    template<class size_vec_t> ubmatrix_array matrix_slice() {
    }

#endif

    /// \name Resize method
    //@{
    /** \brief Resize the tensor to rank \c rank with sizes
        given in \c dim
	
        The parameter \c dim must be a vector of sizes with a length
        equal to \c rank. This resize method is always destructive.
	
        If the user requests any of the sizes to be zero, this
        function will call the error handler.
    */
    template<class size_vec_t> 
    void resize(size_t rank, const size_vec_t &dim) {
      for(size_t i=0;i<rank;i++) {
        if (dim[i]==0) {
          O2SCL_ERR((((std::string)"Requested zero size with non-zero ")+
                     "rank for index "+szttos(i)+" in tensor::"+
                     "resize().").c_str(),exc_einval);
        }
      }
      rk=rank;
      size.resize(rk);
      if (rk==0) {
        data.resize(0);
      } else {
        size_t tot=1;
        for(size_t i=0;i<rk;i++) {
          size[i]=dim[i];
          tot*=size[i];
        }
        data.resize(tot);
      }
      return;
    }
    //@}

    /// \name Size functions
    //@{
    /// Return the rank of the tensor
    size_t get_rank() const { return rk; }

    /// Returns the size of the ith index
    size_t get_size(size_t i) const { 
      if (i>=rk) {
        O2SCL_ERR((((std::string)"Specified index ")+szttos(i)+
                   " greater than or equal to rank "+szttos(rk)+
                   " in tensor::get_size()").c_str(),
                  exc_einval);
      }
      return size[i]; 
    }
    
    /// Return the full vector of sizes
    const vec_size_t &get_size_arr() const {
      return size;
    }

    /// Return the full data vector
    const vec_t &get_data() const {
      return data;
    }

    /** \brief Returns the size of the tensor (the product of 
        the sizes over every index)
    */
    size_t total_size() const { 
      if (rk==0) return 0;
      size_t tot=1;
      for(size_t i=0;i<rk;i++) tot*=size[i];
      return tot;
    }
    //@}

    /// \name Index manipulation
    //@{
    /// Pack the indices into a single vector index
    template<class size_vec_t> 
    size_t pack_indices(const size_vec_t &index) {
      if (rk==0) {
        O2SCL_ERR("Empty tensor in tensor::pack_indices().",exc_einval);
        return 0;
      }
      if (index[0]>=size[0]) {
        O2SCL_ERR((((std::string)"Value of index[0]=")+szttos(index[0])+
                   " greater than or equal to size[0]="+
                   szttos(size[0])+" in tensor::pack_indices().").c_str(),
                  exc_eindex);
      }
      size_t ix=index[0];
      for(size_t i=1;i<rk;i++) {
        if (index[i]>=size[i]) {
          O2SCL_ERR((((std::string)"Value of index[")+szttos(i)+"]="+
                     szttos(index[i])+" greater than or equal to size "+
                     szttos(size[i])+" in tensor::pack_indices().").c_str(),
                    exc_eindex);
        }
        ix*=size[i];
        ix+=index[i];
      }
      return ix;
    }
    
    /// Unpack the single vector index into indices
    template<class size_vec_t> 
    void unpack_index(size_t ix, size_vec_t &index) {
      if (ix>=total_size()) {
        O2SCL_ERR((((std::string)"Value of index ")+szttos(ix)+
                   " greater than or equal to total size "+
                   szttos(total_size())+
                   " in tensor::unpack_index().").c_str(),
                  exc_eindex);
        return;
      }
      size_t ix2, sub_size;
      for(size_t i=0;i<rk;i++) {
        if (i==rk-1) {
          index[i]=ix;
        } else {
          sub_size=1;
          for(size_t j=i+1;j<rk;j++) sub_size*=size[j];
          index[i]=ix/sub_size;
          // (Remember we're doing integer arithmetic here.)
          ix-=sub_size*(ix/sub_size);
        }
      }
      return;
    }
    //@}

    /// \name Minimum, maximum, and sum
    //@{
    /** \brief Compute the minimum value in the tensor
     */
    data_t min_value() {
      return o2scl::vector_min_value<vec_t,data_t>(total_size(),data);
    }

    /** \brief Compute the index of the minimum value in the tensor
     */
    void min_index(vec_size_t &index) {
      size_t ix=o2scl::vector_min_index<vec_t,data_t>(total_size(),data);
      unpack_index(ix,index);
      return;
    }

    /** \brief Compute the index of the minimum value in the tensor
        and return the minimum
    */
    void min(vec_size_t &index, data_t &val) {
      size_t ix;
      o2scl::vector_min<vec_t,data_t>(total_size(),data,ix,val);
      unpack_index(ix,index);
      return; 
    }

    /** \brief Compute the maximum value in the tensor
     */
    data_t max_value() {
      return o2scl::vector_max_value<vec_t,data_t>(total_size(),data);
    }

    /** \brief Compute the index of the maximum value in the tensor
     */
    void max_index(vec_size_t &index) {
      size_t ix=o2scl::vector_max_index<vec_t,data_t>(total_size(),data);
      unpack_index(ix,index);
      return; 
    }

    /** \brief Compute the index and value of the maximum value in the tensor
        and return the maximum
    */
    void max(vec_size_t &index, data_t &val) {
      size_t ix;
      o2scl::vector_max<vec_t,data_t>(total_size(),data,ix,val);
      unpack_index(ix,index);
      return;
    }

    /** \brief Compute the minimum and maximum values in the tensor
     */
    void minmax_value(data_t &min, data_t &max) {
      return o2scl::vector_minmax_value<vec_t,data_t>(total_size(),
                                                      data,min,max);
    }

    /** \brief Compute the indices of the minimum and maximum values 
        in the tensor
    */
    void minmax_index(vec_size_t &index_min, vec_size_t &index_max) {
      size_t ix_min, ix_max;
      o2scl::vector_minmax_index<vec_t,data_t>
        (total_size(),data,ix_min,ix_max);
      unpack_index(ix_min,index_min);
      unpack_index(ix_max,index_max);
      return;
    }

    /** \brief Compute the indices and values of the maximum and minimum
        in the tensor
    */
    void minmax(vec_size_t &index_min, data_t &min,
                vec_size_t &index_max, data_t &max) {
      size_t ix_min, ix_max;
      o2scl::vector_minmax<vec_t,data_t>(total_size(),data,ix_min,min,
                                         ix_max,max);
      unpack_index(ix_min,index_min);
      unpack_index(ix_max,index_max);
      return;
    }

    /** \brief Return the sum over every element in the tensor
     */
    double total_sum() const { 
      if (rk==0) return 0.0;
      double tot=0.0;
      for(size_t i=0;i<data.size();i++) {
        tot+=data[i];
      }
      return tot;
    }
    //@}
  
    /// \name Slicing and converting to table3d objects
    //@{
    /** \brief Copy to a \ref o2scl::table3d object by
        summing over all but two indices
    */
    void copy_table3d_sum
    (size_t ix_x, size_t ix_y, table3d &tab, std::string x_name="x",
     std::string y_name="y", std::string slice_name="z") {

      tab.clear();
      
      // Get current table3d grid
      size_t nx, ny;
      tab.get_size(nx,ny);
    
      if (nx==0 && ny==0) {
      
        if (x_name.length()==0) x_name="x";
        if (y_name.length()==0) y_name="y";
      
        // If there's no grid, then create a grid in the table3d
        // object which just enumerates the indices
        std::vector<double> grid_x(size[ix_x]), grid_y(size[ix_y]);
        for(size_t i=0;i<size[ix_x];i++) {
          grid_x[i]=((double)i);
        }
        for(size_t i=0;i<size[ix_y];i++) {
          grid_y[i]=((double)i);
        }
        tab.set_xy("x",grid_x.size(),grid_x,
                   "y",grid_y.size(),grid_y);
        // Now that the grid is set, get nx and ny
        tab.get_size(nx,ny);
      }
    
      // Check that the grids are commensurate
      if (nx!=this->size[ix_x] || ny!=this->size[ix_y]) {
        O2SCL_ERR2("Grids not commensurate in ",
                   "tensor_grid::copy_table3d_sum().",exc_einval);
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

    /** \brief Copy to a table3d object by fixing two indices

        \warning The vector input \c index must be initialized
        before calling this function so that all elements in the
        vector (except for those at index \c ix_x and \c ix_y)
        are specified. If this is not the case, then this function
        will return unpredictable results.
     */
    template<class size_vec2_t> 
    void copy_table3d(size_t ix_x, size_t ix_y, size_vec2_t &index, 
                      table3d &tab, std::string x_name="x",
                      std::string y_name="y",
                      std::string slice_name="z") const {
      
      if (ix_x>=this->rk || ix_y>=this->rk || ix_x==ix_y) {
        O2SCL_ERR2("Either indices greater than rank or x and y ind",
                   "ices equal in tensor_grid::copy_table3d_align().",
                   exc_efailed);
      }

      size_t nx, ny;
      nx=get_size(ix_x);
      ny=get_size(ix_y);

      // Set the table3d grid using the indexes
      tab.clear();
      tab.set_xy(x_name,uniform_grid_end_width<double>(0,nx-1,1),
                 y_name,uniform_grid_end_width<double>(0,ny-1,1));
      
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
    
    /** \brief Copy to a table3d object by fixing two indices

        In this function, all other indices are set to 0. The indices
        \c ix_x and \c ix_y may be the same, in which case this
        function effectively performs a trace over those two indices.
     */
    void copy_table3d(size_t ix_x, size_t ix_y, 
                      table3d &tab, std::string x_name="x",
                      std::string y_name="y",
                      std::string slice_name="z") const {
      std::vector<size_t> ix(rk);
      for(size_t i=0;i<rk;i++) ix[i]=0;
      copy_table3d(ix_x,ix_y,ix,tab,x_name,y_name,slice_name);
      return;
    }
    //@}
    
  };

  /** Given a set of index specifications specified in a list of 
      strings, reformat them into a list of \ref o2scl::index_spec
      objects
      
      The tensor rearrange commands use index specifications to
      specify how the tensor should be rearranged. Index
      specifications may be specified as separate arguments e.g.
      "index(1)" "fixed(2,10)" or multiple index specifications may
      be given in a single argument separated by spaces or commas,
      e.g. "index(1) fixed(2,10)" or "index(1),fixed(2,10)". The
      indices begin with 0, the first index so that index 1 is the
      second index. The list of index specification is:
      
      - index(ix): Retain index ix in the new tensor.
      
      - fixed(ix): Fix the value of index ix.
      
      - sum(ix): Sum over the value of index ix
      
      - trace(ix1,ix2): Trace (sum) over indices ix and ix2. If the
      number of entries in either index is smaller than the other,
      then the remaining entries are ignored in the sum.
      
      - reverse(ix): Retain index ix but reverse the order.
      
      - range(ix,start,end): Retain index ix but modify range.
      Ranges include both of their endpoints.
      
      - interp(ix,value) (for tensor_grid): fix index ix by
      interpolating 'value' into the grid for index ix.
      
      - grid(ix,begin,end,n_bins,log) (for tensor_grid): interpolate
      the specified index on a grid to create a new index. If the
      value of log is 1, then the grid is logarithmic.
      
      - gridw(ix,begin,end,bin_width,log) (for tensor_grid):
      interpolate the specified index on a grid with a fixed bin
      width to create a new index. If the value of log is 1, then
      the grid is logarithmic and the bin_width is the
      multiplicative factor between bin edges.
      
      Note that the index specifications which result in a tensor
      index (all except 'fixed', 'sum', 'trace' and 'interp') must
      be given in the order they should appear in the tensor which
      results. Also, the 'rearrange' commands require that the
      result of the rearrangement must have at least one index left.
      
      Examples:
      
      index(1),index(0) - take the transpose of a rank 2 tensor
      (i.e. a matrix)
      
      index(1),fixed(2,0),index(0) - fix the value of index 2 (i.e.
      the third index) to zero and transpose the other two indices
      
      fixed(2,0),index(1),index(0) - same as above
      
  */
  int strings_to_indexes(std::vector<std::string> sv2,
                          std::vector<o2scl::index_spec> &vis,
                          int verbose=0, bool err_on_fail=false);
  
  /** \brief Take a set of index specifications contained in a
      single string \c str and arrange them in \c sv
      
      \verbatim embed:rst
      
      .. todo:: 
      
         In tensor::index_spec_preprocess():

         Future: Improve this to be more intelligent about whitespace
         and other characters between index specifications. Right now,
         this function fails if there are, e.g. two spaces between
         index specs.
         
      \endverbatim
  */
  void index_spec_preprocess(std::string str,
                              std::vector<std::string> &sv, int verbose=0);
  
  /** \brief Rearrange, sum and copy current tensor to a new tensor
      
      \todo We need to check all of the degenerate cases, for
      example, a range spec. with only one element, a grid specification
      with only one element in the result, etc.
      
      \verbatim embed:rst
        
      .. todo:: 

         In tensor::rarrange_and_copy(): 

         Future: Return a scalar if possible as a rank 1 tensor with
         1 element.
           
      \endverbatim
  */
  template<class tensor_t, class data_t>
  tensor_t rearrange_and_copy(const tensor_t &t,
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
    
    // Size of sums
    std::vector<size_t> sum_sizes;
    
    // Loop through the index specifications and add them to spec_old
    // and spec_new (if necessary). This loop also determines the rank
    // of the new tensor, "rank_new", and the sizes of the indices
    // "size_new". No actual copying or summing is done in this loop
    // yet.
      
    for(size_t i=0;i<spec.size();i++) {
      
      if (spec[i].type==index_spec::index ||
          spec[i].type==index_spec::reverse) {
        if (spec[i].ix1>=rank_old) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (index,reverse) in ",
                       "tensor::rearrange_and_copy().",o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index " << spec[i].ix1
                        << " too large (index,reverse) in "
                        << "tensor in tensor::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_t();
          }
        }
        size_new.push_back(t.get_size(spec[i].ix1));
        // Use ix1 to store the destination index (which is
        // at this point equal to rank_new)
        spec_old[spec[i].ix1]=index_spec(spec[i].type,rank_new);
        spec_new.push_back(index_spec(spec[i].type,spec[i].ix1));
        rank_new++;
        
      } else if (spec[i].type==index_spec::range) {
        
        if (verbose>2) {
          std::cout << "In range " << spec[i].ix1 << " "
                    << spec[i].ix2 << " " << spec[i].ix3 << std::endl;
        }
        if (spec[i].ix1>=rank_old ||
            spec[i].ix2>=t.get_size(spec[i].ix1) ||
            spec[i].ix3>=t.get_size(spec[i].ix1)) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (range) in ",
                       "tensor::rearrange_and_copy().",o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index " << spec[i].ix1 << " "
                        << spec[i].ix2 << " " << spec[i].ix3
                        << " too large (range) in "
                        << "tensor in tensor::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_t();
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
          index_spec(spec[i].type,rank_new,spec[i].ix2,spec[i].ix3);
        spec_new.push_back
          (index_spec(spec[i].type,spec[i].ix1,
                      spec[i].ix2,spec[i].ix3));
        rank_new++;
        
        if (verbose>2) {
          std::cout << "Out range " << size_new[size_new.size()-1]
                    << std::endl;
        }
        
      } else if (spec[i].type==index_spec::trace) {
        
        if (spec[i].ix1>=rank_old || spec[i].ix2>=rank_old) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (trace) in ",
                       "tensor::rearrange_and_copy().",o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Indices " << spec[i].ix1 << " or "
                        << spec[i].ix2 << " too large (trace) in "
                        << "tensor in tensor::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_t();
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
        spec_old[spec[i].ix1]=index_spec(spec[i].type,
                                         spec[i].ix1,spec[i].ix2);
        spec_old[spec[i].ix2]=index_spec(spec[i].type,
                                         spec[i].ix2,spec[i].ix1);
        
      } else if (spec[i].type==index_spec::sum) {
        
        if (spec[i].ix1>=rank_old) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (sum) in ",
                       "tensor::rearrange_and_copy().",o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index " << spec[i].ix1
                        << " too large (sum) in "
                        << "tensor in tensor::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_t();
          }
        }
        n_sum_loop*=t.get_size(spec[i].ix1);
        sum_sizes.push_back(t.get_size(spec[i].ix1));
        spec_old[spec[i].ix1]=index_spec(spec[i].type,
                                         spec[i].ix1,spec[i].ix2,0);
        
      } else if (spec[i].type==index_spec::fixed) {
        
        if (spec[i].ix1>=rank_old ||
            spec[i].ix2>=t.get_size(spec[i].ix1)) {
          if (err_on_fail) {
            O2SCL_ERR2("Index too large (fixed) in ",
                       "tensor::rearrange_and_copy().",o2scl::exc_einval);
          } else {
            if (verbose>0) {
              std::cout << "Index too large (fixed) in "
                        << "tensor in tensor::rearrange_and_copy()."
                        << std::endl;
            }
            return tensor_t();
          }
        }
        
        // Use ix1 to store the destination index (which is
        // at this point equal to rank_new)
        spec_old[spec[i].ix1]=index_spec(spec[i].type,
                                         rank_new,spec[i].ix2);
        
      } else {
        
        if (err_on_fail) {
          O2SCL_ERR2("Index specification type not allowed in ",
                     "tensor::rearrange_and_copy().",o2scl::exc_einval);
        } else {
          if (verbose>0) {
            std::cout << "Index specification type not allowed in "
                      << "tensor::rearrange_and_copy()." << std::endl;
          }
          return tensor_t();
        }
        
      }
    }
    
    size_t n_sums=sum_sizes.size();
    
    // Call the error handler if the input is invalid
    if (rank_new==0) {
      if (err_on_fail) {
        O2SCL_ERR2("Zero new indices in ",
                   "tensor::rearrange_and_copy().",o2scl::exc_einval);
      } else {
        if (verbose>0) {
          std::cout << "Zero new indices in "
                    << "tensor::rearrange_and_copy()." << std::endl;
        }
        return tensor_t();
      }
    }
    
    for(size_t i=0;i<rank_old;i++) {
      if (spec_old[i].type==index_spec::empty) {
        if (err_on_fail) {
          O2SCL_ERR2("Not all indices accounted for in ",
                     "tensor::rearrange_and_copy().",o2scl::exc_einval);
        } else {
          if (verbose>0) {
            std::cout << "Index " << i << " not accounted for in "
                      << "tensor::rearrange_and_copy()." << std::endl;
          }
          return tensor_t();
        }
      }
    }
    
    // Verbose output if necessary
    if (verbose>0) {
      std::cout << "rearrange_and_copy(): using a " << rank_old
                << " rank tensor to create\n  a new "
                << rank_new << " rank tensor." << std::endl;
    }
    if (verbose>1) {
      for(size_t i=0;i<rank_old;i++) {
        std::cout << "  Old index " << i;
        if (spec_old[i].type==index_spec::index) {
          std::cout << " is being remapped to new index " << spec_old[i].ix1
                    << "." << std::endl;
        } else if (spec_old[i].type==index_spec::range) {
          std::cout << " is being remapped to new index " << spec_old[i].ix1
                    << " with a range from " << spec_old[i].ix2
                    << " to " << spec_old[i].ix3 << "." << std::endl;
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
        }
      }
      for(size_t i=0;i<rank_new;i++) {
        std::cout << "  New index " << i;
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
        }
      }
    }
    
    // Create the new tensor object
    tensor_t t_new(rank_new,size_new);
    
    // Index arrays
    std::vector<size_t> ix_new(rank_new);
    std::vector<size_t> ix_old(rank_old);
    std::vector<size_t> sum_ix(n_sums);
    
    // Loop over the new tensor object
    for(size_t i=0;i<t_new.total_size();i++) {
      
      // Find the location in the new tensor object
      t_new.unpack_index(i,ix_new);
      
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
        }
      }
      
      data_t val=0;
      
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
          std::cout << "rearrange_and_copy(): n_sum_loop: "
                    << n_sum_loop << " n_sums: "
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
              std::cout << "X: " << cnt << " " << sum_ix.size() << std::endl;
              O2SCL_ERR2("Bad sync 1 in sum_ix in ",
                         "tensor::rearrange_and_copy()",o2scl::exc_esanity);
            }
            ix_old[k]=sum_ix[cnt];
            cnt++;
          } else if (spec_old[k].type==index_spec::trace &&
                     spec_old[k].ix1<spec_old[k].ix2) {
            if (cnt>=sum_ix.size()) {
              std::cout << "X: " << cnt << " " << sum_ix.size() << std::endl;
              O2SCL_ERR2("Bad sync 2 in sum_ix in ",
                         "tensor::rearrange_and_copy()",o2scl::exc_esanity);
            }
            ix_old[spec_old[k].ix1]=sum_ix[cnt];
            ix_old[spec_old[k].ix2]=sum_ix[cnt];
            cnt++;
          }
        }
        
        if (verbose>2) {
          std::cout << "Here old: ";
          vector_out(std::cout,ix_old,true);
          std::cout << "Here new: ";
          vector_out(std::cout,ix_new,true);
        }
        val+=t.get(ix_old);
        
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
  tensor_t rearrange_and_copy(const tensor_t &t, std::string spec,
                               int verbose=0, bool err_on_fail=true) {
    
    std::vector<std::string> sv2;
    index_spec_preprocess(spec,sv2);
    std::vector<o2scl::index_spec> vis;
    strings_to_indexes(sv2,vis,verbose);
    return rearrange_and_copy<tensor_t,data_t>(t,vis,verbose,err_on_fail);
  }
  
  /** \brief Rank 1 tensor
   */
  template<class data_t=double, class vec_t=std::vector<data_t>, 
           class vec_size_t=std::vector<size_t> > class tensor1 : 
    public tensor<data_t,vec_t,vec_size_t> {
    
  public:
  
    /// Create an empty tensor
    tensor1() : tensor<data_t,vec_t,vec_size_t>() {}
  
    /// Create a rank 1 tensory of size \c sz
    tensor1(size_t sz) : tensor<data_t,vec_t,vec_size_t>() {
      vec_size_t sizex(1);
      sizex[0]=sz;
      this->resize(1,sizex);
    }
  
    /// \name Method to check for valid object
    //@{
    /** \brief Check that the \ref o2scl::tensor1 object is valid
     */
    void is_valid() const {
      tensor<double,vec_t,vec_size_t>::is_valid();
      if (this->rk>1) {
        O2SCL_ERR2("Rank is neither 0 nor 1 in ",
                   "tensor1::is_valid().",
                   o2scl::exc_esanity);
      
      }
      return;
    }
    //@}

    /// \name Specialized get and set functions
    //@{
    /// Get the element indexed by \c ix
    data_t &get(size_t ix) { 
      return tensor<data_t,vec_t,vec_size_t>::get(&ix); 
    }
  
    /// Get the element indexed by \c ix
    const data_t &get(size_t ix) const { 
      return tensor<data_t,vec_t,vec_size_t>::get(&ix); 
    }
  
    /// Set the element indexed by \c index to value \c val
    void set(size_t index, data_t val) 
    { tensor<data_t,vec_t,vec_size_t>::set(&index,val); }
  
    /** \brief Set the element indexed by \c index to value \c val
      
        (We have to explicitly provide this version since the set()
        function is overloaded in this child of \ref tensor.)
    */
    template<class size_vec_t>
    void set(const size_vec_t &index, data_t val) {
      tensor<data_t,vec_t,vec_size_t>::set(index,val);
    }
    //@}
  
    /// \name Specialized operator functions
    //@{
    /// Get an element using array-like indexing
    data_t &operator[](size_t ix) { return this->data[ix]; }
  
    /// Get an element using array-like indexing (const version)
    const data_t &operator[](size_t ix) const { return this->data[ix]; }
  
    /// Get an element using operator()
    data_t &operator()(size_t ix) { return this->data[ix]; }
  
    /// Get an element using operator() (const version)
    const data_t &operator()(size_t ix) const { return this->data[ix]; }
    //@}
  
  };

  /** \brief Rank 2 tensor
   */
  template<class data_t=double, class vec_t=std::vector<data_t>, 
           class vec_size_t=std::vector<size_t> > class tensor2 : 
    public tensor<data_t,vec_t,vec_size_t> {

  public:

    /// Create an empty tensor
    tensor2() : tensor<data_t,vec_t,vec_size_t>() {}

    /// Create a rank 2 tensor of size \c (sz,sz2)
    tensor2(size_t sz, size_t sz2) : tensor<data_t,vec_t,vec_size_t>() {
      this->rk=2;
      this->size.resize(2);
      this->size[0]=sz;
      this->size[1]=sz2;
      size_t tot=sz*sz2;
      this->data.resize(tot);
    }
	
    /// \name Method to check for valid object
    //@{
    /** \brief Check that the \ref o2scl::tensor2 object is valid
     */
    void is_valid() const {
      tensor<double,vec_t,vec_size_t>::is_valid();
      if (this->rk!=0 && this->rk!=2) {
        O2SCL_ERR2("Rank is neither 0 nor 2 in ",
                   "tensor2::is_valid().",
                   o2scl::exc_esanity);
      
      }
      return;
    }
    //@}
  
    /// \name Specialized get and set functions
    //@{
    /// Get the element indexed by \c (ix1,ix2)
    data_t &get(size_t ix1, size_t ix2) { 
      size_t sz[2]={ix1,ix2};
      return tensor<data_t,vec_t,vec_size_t>::get(sz); 
    }

    /// Get the element indexed by \c (ix1,ix2)
    const data_t &get(size_t ix1, size_t ix2) const { 
      size_t sz[2]={ix1,ix2};
      return tensor<data_t,vec_t,vec_size_t>::get(sz); 
    }

    /// Set the element indexed by \c (ix1,ix2) to value \c val
    void set(size_t ix1, size_t ix2, data_t val) {
      size_t sz[2]={ix1,ix2};
      tensor<data_t,vec_t,vec_size_t>::set(sz,val); 
      return;
    }

    /** \brief Set the element indexed by \c index to value \c val

        (We have to explicitly provide this version since the set()
        function is overloaded in this child of \ref tensor.)
    */
    template<class size_vec_t>
    void set(const size_vec_t &index, data_t val) {
      tensor<data_t,vec_t,vec_size_t>::set(index,val);
      return;
    }

    /// Get the element indexed by \c (ix1,ix2)
    data_t &operator()(size_t ix, size_t iy) 
    { return this->data[ix*this->size[1]+iy]; }

    /// Get the element indexed by \c (ix1,ix2) (const version)
    const data_t &operator()(size_t ix, size_t iy) const
    { return this->data[ix*this->size[1]+iy]; }
    //@}
  };
  
  /** \brief Rank 3 tensor
   */
  template<class data_t=double, class vec_t=std::vector<data_t>, 
           class vec_size_t=std::vector<size_t> > class tensor3 : 
    public tensor<data_t,vec_t,vec_size_t> {

  public:

    /// Create an empty tensor
    tensor3() : tensor<data_t,vec_t,vec_size_t>() {}

    /// Create a rank 3 tensor of size \c (sz,sz2,sz3)
    tensor3(size_t sz, size_t sz2, size_t sz3) : 
      tensor<data_t,vec_t,vec_size_t>() {
      this->rk=3;
      this->size.resize(3);
      this->size[0]=sz;
      this->size[1]=sz2;
      this->size[2]=sz3;
      size_t tot=sz*sz2*sz3;
      this->data.resize(tot);
    }

    /// \name Method to check for valid object
    //@{
    /** \brief Check that the \ref o2scl::tensor3 object is valid
     */
    void is_valid() const {
      tensor<double,vec_t,vec_size_t>::is_valid();
      if (this->rk!=0 && this->rk!=3) {
        O2SCL_ERR2("Rank is neither 0 nor 3 in ",
                   "tensor3::is_valid().",
                   o2scl::exc_esanity);
      
      }
      return;
    }
    //@}
  
    /// \name Specialized get and set functions
    //@{
    /// Get the element indexed by \c (ix1,ix2,ix3)
    data_t &get(size_t ix1, size_t ix2, size_t ix3) { 
      size_t sz[3]={ix1,ix2,ix3};
      return tensor<data_t,vec_t,vec_size_t>::get(sz); 
    }

    /// Get the element indexed by \c (ix1,ix2,ix3)
    const data_t &get(size_t ix1, size_t ix2, size_t ix3) const { 
      size_t sz[3]={ix1,ix2,ix3};
      return tensor<data_t,vec_t,vec_size_t>::get(sz); 
    }

    /// Set the element indexed by \c (ix1,ix2,ix3) to value \c val
    void set(size_t ix1, size_t ix2, size_t ix3, data_t val) {
      size_t sz[3]={ix1,ix2,ix3};
      tensor<data_t,vec_t,vec_size_t>::set(sz,val); 
      return;
    }

    /** \brief Set the element indexed by \c index to value \c val

        (We have to explicitly provide this version since the set()
        function is overloaded in this child of \ref tensor.)
    */
    template<class size_vec_t>
    void set(const size_vec_t &index, data_t val) {
      tensor<data_t,vec_t,vec_size_t>::set(index,val);
      return;
    }
    //@}

  };
  
  /** \brief Rank 4 tensor
   */
  template<class data_t=double, class vec_t=std::vector<data_t>, 
           class vec_size_t=std::vector<size_t> > class tensor4 : 
    public tensor<data_t,vec_t,vec_size_t> {

  public:

    /// Create an empty tensor
    tensor4() : tensor<data_t,vec_t,vec_size_t>() {}

    /// Create a rank 4 tensor of size \c (sz,sz2,sz3,sz4)
    tensor4(size_t sz, size_t sz2, size_t sz3, size_t sz4) : 
      tensor<data_t,vec_t,vec_size_t>() {
      this->rk=4;
      this->size.resize(4);
      this->size[0]=sz;
      this->size[1]=sz2;
      this->size[2]=sz3;
      this->size[3]=sz4;
      size_t tot=sz*sz2*sz3*sz4;
      this->data.resize(tot);
    }
	
    /// \name Method to check for valid object
    //@{
    /** \brief Check that the \ref o2scl::tensor4 object is valid
     */
    void is_valid() const {
      tensor<double,vec_t,vec_size_t>::is_valid();
      if (this->rk!=0 && this->rk!=4) {
        O2SCL_ERR2("Rank is neither 0 nor 4 in ",
                   "tensor4::is_valid().",
                   o2scl::exc_esanity);
      
      }
      return;
    }
    //@}
  
    /// \name Specialized get and set functions
    //@{
    /// Get the element indexed by \c (ix1,ix2,ix3,ix4)
    data_t &get(size_t ix1, size_t ix2, size_t ix3, size_t ix4) { 
      size_t sz[4]={ix1,ix2,ix3,ix4};
      return tensor<data_t,vec_t,vec_size_t>::get(sz); 
    }

    /// Get the element indexed by \c (ix1,ix2,ix3,ix4)
    const data_t &get(size_t ix1, size_t ix2, size_t ix3, 
                      size_t ix4) const { 
      size_t sz[4]={ix1,ix2,ix3,ix4};
      return tensor<data_t,vec_t,vec_size_t>::get(sz); 
    }

    /// Set the element indexed by \c (ix1,ix2,ix3,ix4) to value \c val
    void set(size_t ix1, size_t ix2, size_t ix3, size_t ix4, 
             data_t val) {
      size_t sz[4]={ix1,ix2,ix3,ix4};
      tensor<data_t,vec_t,vec_size_t>::set(sz,val); 
      return;
    }

    /** \brief Set the element indexed by \c index to value \c val

        (We have to explicitly provide this version since the set()
        function is overloaded in this child of \ref tensor.)
    */
    template<class size_vec_t>
    void set(const size_vec_t &index, data_t val) {
      tensor<data_t,vec_t,vec_size_t>::set(index,val);
      return;
    }
    //@}
  };

  /// \name Tensor functions in src/base/tensor.h
  //@{
  /** \brief Output a tensor to a stream
   */
  template<class tensor_t>
  void tensor_out(std::ostream &os, tensor_t &t, bool pretty=true) {

    if (pretty) {
      
      size_t rk=t.get_rank();
      os << "rank: " << rk << " sizes: ";
      for(size_t i=0;i<rk;i++) {
	os << t.get_size(i) << " ";
      }
      os << std::endl;
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
      auto &data=t.get_data();
      for(size_t i=0;i<t.total_size();i++) {
	os << data[i] << " ";
	if (i%10==9) os << std::endl;
      }
      os << std::endl;
      
    }
    
    return;
  }

  /** \brief Compare two tensors for equality
   */
  template<class data_t, class vec_t, class vec_size_t,
           class data2_t, class vec2_t, class vec2_size_t>
  bool operator==(const tensor<data_t,vec_t,vec_size_t> &t1,
                  const tensor<data2_t,vec2_t,vec2_size_t> &t2) {
    if (t1.get_rank()!=t2.get_rank()) return false;
    for(size_t i=0;i<t1.get_rank();i++) {
      if (t1.get_size(i)!=t2.get_size(i)) return false;
    }
    const vec_t &v1=t1.get_data();
    const vec2_t &v2=t2.get_data();
    for(size_t i=0;i<t1.total_size();i++) {
      if (v1[i]!=v2[i]) return false;
    }
    return true;
  }

#ifdef O2SCL_NEVER_DEFINED  
  /** \brief Add two tensors (possibly of different types)
   */
  template<class data_t, class vec_t, class vec_size_t,
           class data2_t, class vec2_t, class vec2_size_t>
  tensor<data_t,vec_t,vec_size_t>
  operator+(const tensor<data_t,vec_t,vec_size_t> &t1,
            const tensor<data2_t,vec2_t,vec2_size_t> &t2) {
    if (t1.get_rank()!=t2.get_rank()) return false;
    vec_size_t ndims(t1.get_rank);
    for(size_t i=0;i<t1.get_rank();i++) {
      if (t1.get_size(i)!=t2.get_size(i)) return false;
      ndims[i]=t1.get_size(i);
    }
    tensor<data_t,vec_t,vec_size_t> t3(t1.get_rank(),ndims);
    const vec_t &v1=t1.get_data();
    const vec_t &v2=t2.get_data();
    const vec_t &v3=t3.get_data();
    for(size_t i=0;i<t1.total_size();i++) {
      v3[i]=v1[i]+v2[i];
    }
    return t3;
  }

  /** \brief Subtract two tensors (possibly of different types)
   */
  template<class data_t, class vec_t, class vec_size_t,
           class data2_t, class vec2_t, class vec2_size_t>
  tensor<data_t,vec_t,vec_size_t>
  operator-(const tensor<data_t,vec_t,vec_size_t> &t1,
            const tensor<data2_t,vec2_t,vec2_size_t> &t2) {
    if (t1.get_rank()!=t2.get_rank()) return false;
    vec_size_t ndims(t1.get_rank);
    for(size_t i=0;i<t1.get_rank();i++) {
      if (t1.get_size(i)!=t2.get_size(i)) return false;
      ndims[i]=t1.get_size(i);
    }
    tensor<data_t,vec_t,vec_size_t> t3(t1.get_rank(),ndims);
    const vec_t &v1=t1.get_data();
    const vec_t &v2=t2.get_data();
    const vec_t &v3=t3.get_data();
    for(size_t i=0;i<t1.total_size();i++) {
      v3[i]=v1[i]-v2[i];
    }
    return t3;
  }
#endif
  //@}
  
}

#endif



