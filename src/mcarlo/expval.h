/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2021, Andrew W. Steiner
  
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
#ifndef EXPVAL_H
#define EXPVAL_H

/** \file expval.h
    \brief File defining 'expectation value' objects
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/vec_stats.h>
#include <o2scl/tensor.h>

// Forward definition of the expval classes for HDF I/O
namespace o2scl {
  class expval_scalar;
  class expval_vector;
  class expval_matrix;
}

// Forward definition of HDF I/O to extend friendship
namespace o2scl_hdf { 
  class hdf_file; 
  void hdf_input_n(hdf_file &hf, o2scl::expval_scalar &t, std::string &name);
  void hdf_output(hdf_file &hf, o2scl::expval_scalar &t, std::string name);
  void hdf_input_n(hdf_file &hf, o2scl::expval_vector &t, std::string &name);
  void hdf_output(hdf_file &hf, o2scl::expval_vector &t, std::string name);
  void hdf_input_n(hdf_file &hf, o2scl::expval_matrix &t, std::string &name);
  void hdf_output(hdf_file &hf, o2scl::expval_matrix &t, std::string name);
}

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Expectation value base class
      
      \verbatim embed:rst
      See the :ref:`Analysis of results from numerical simulations`
      section of the User's guide for basic information about this
      class and its children.
      \endverbatim

      This base class need not be directly instantiated by the casual
      end-user, but provides basic functionality for \ref
      expval_scalar, \ref expval_vector, and \ref expval_matrix.

      \hline

      Internally, neither \ref nblocks nor \ref nperblock should
      ever be zero. This is checked by \ref is_valid() .
  */
  class expval_base {

  protected:

    /// Index denoting the current block number
    size_t iblock;
    
    /// Index for the number of values in the current block
    size_t i;

    /** \brief Total number of blocks (default 1)
	
	This should never be zero.
    */
    size_t nblocks;
    
    /** \brief Number of measurements per block (default 1)

	This should never be zero.
    */
    size_t nperblock;

  public:

    /** \brief Create with \c n_blocks blocks and \c n_per_block points
	per block

	If this is called with a value of zero for either \c n_blocks
	or \c n_per_block, then the error handler is called.
    */
    expval_base(size_t n_blocks=1, size_t n_per_block=1);

    virtual ~expval_base();
    
    /// Copy constructor
    expval_base(const expval_base &ev);

    /// Copy constructor with <tt>operator=()</tt>
    expval_base &operator=(const expval_base &ev);

    /// The name of the expectation value
    std::string name;
    
    /// The shortened name
    std::string short_name;

    /** \brief Reset for \c n_blocks blocks and \c n_per_block points
	per block

	This function resets the currently stored data to zero by
	calling \ref free(). If this is called with a value of zero
	for \c n_blocks, then the value 1 is assumed.
    */
    virtual void set_blocks(size_t n_blocks, size_t n_per_block);

    /** \brief Get the number of blocks and the number of points per
	block
    */
    virtual void get_blocks(size_t &n_blocks, size_t &n_per_block) const;
    
    /** \brief Free allocated data (but do not change the current values
	of \c n_blocks or \c n_per_block)
    */
    virtual void free();

    /** \brief Get the block index and the index within the current block
     */
    virtual void get_block_indices(size_t &i_block, 
				   size_t &i_curr_block) const;

    /** \brief Returns true if all blocks have been stored

	This reports true when exactly \c n_blocks times \c
	n_per_block data points have been added.
    */
    virtual bool finished() const;

    /** \brief Report progress as a fraction between zero to one 
	(inclusive)
	
	When \c n_per_block is nonzero, this reports the total
	progress on all blocks, reporting \c 1.0 only when all \c
	n_blocks times \c n_per_block data points have been added. If
	more data is added after this function reports 1.0, then the
	blocks are rearranged and progress() will report something
	near 0.5 again.
    */
    virtual double progress() const;

    /// Internal consistency check
    void is_valid() const {
      if (nblocks==0 || nperblock==0) {
	O2SCL_ERR2("Either 'nblocks' or 'n_per_block' is zero in ",
		   "expval_base::is_valid().",exc_efailed);
      } 
    }

  };
  
  /** \brief Scalar expectation value

      See \ref expval_base for some general notes on 
      this and related classes. 

      This represents the expectation value of a scalar
      double-precision quantity over several measurements.
  */
  class expval_scalar : public expval_base {
    
  public:

    typedef boost::numeric::ublas::vector<double> ubvector;

  protected:
    
    /** \brief The average for each block
	
	This is a vector of length \c nblocks.
    */
    ubvector vals;

  public:
    
    /// The current rolling average
    double current;

    /** \brief Create with \c n_blocks blocks and \c n_per_block points
	block
    */
    expval_scalar(size_t n_blocks=1, size_t n_per_block=1);

    virtual ~expval_scalar();
    
    /// Copy constructor
    expval_scalar(const expval_scalar &ev);

    /// Copy constructor
    expval_scalar &operator=(const expval_scalar &ev);

    /** \brief Reset for \c n_blocks blocks and \c n_per_block points
	block
    */
    virtual void set_blocks(size_t n_blocks, size_t n_per_block);

    /** \brief Free allocated data (but do not change the current values
	of \c n_blocks or \c n_per_block)
    */
    virtual void free();

    /// Add measurement of value \c val
    virtual void add(double val);

    /// \name Report statistics
    //@{
    /** \brief Report current average, standard deviation, and 
	the error in the average and include block information
    */
    virtual void current_avg_stats(double &avg, double &std_dev, 
				   double &avg_err, size_t &m_block,
				   size_t &m_per_block) const;

    /** \brief Report current average, standard deviation, and 
	the error in the average
    */
    virtual void current_avg(double &avg, double &std_dev, 
			     double &avg_err) const;
    
    /** \brief Report average, standard deviation, and 
	the error in the average assuming a new block size

	\future Use recurrence relation for averages here
	rather than dividing at the end.
    */
    virtual void reblock_avg_stats(size_t new_blocks, double &avg, 
				   double &std_dev, double &avg_err,
				   size_t &m_per_block) const;

    /** \brief Report average, standard deviation, and 
	the error in the average assuming a new block size
    */
    virtual void reblock_avg(size_t new_blocks, double &avg, 
			     double &std_dev, double &avg_err) const;
    //@}

    /// \name Direct manipulation of the stored data 
    //@{
    /// Return the current data for all blocks
    const ubvector &get_data() const;

    /// Return the current data for block with index \c i_block
    const double &operator[](size_t i_block) const;

    /// Return the current data for block with index \c i_block
    double &operator[](size_t i_block);

    /// Set the data for all blocks
    template<class vec_t> void set_data(vec_t &v) {
      for(size_t ib=0;ib<nblocks;ib++) {
	vals[ib]=v[ib];
      }
      return;
    }
    //@}

    /// Internal consistency check
    void is_valid() const;

    friend void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf, 
				      expval_scalar &t, 
				      std::string name);
    
    friend void o2scl_hdf::hdf_input_n(o2scl_hdf::hdf_file &hf, expval_scalar &t, 
				     std::string &name);

  };

  /** \brief Vector expectation value

      See \ref expval_base for some general notes on this and related
      classes.

      This is a similar to \ref expval_scalar, except that it allows
      updating and statistics for a set of scalars en masse. The data
      is stored internally in ublas vector and matrix object,
      but the public member functions operate with template types
      which are compatible with any vector class which provides
      <tt>double &operator[]</tt>. It is assumed that each
      call to \ref add() contains a new measurement for all of
      the vector indices. 
  */
  class expval_vector : public expval_base {
    
  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

  protected:
    
    /** \brief The average for each block

	The first (row) index is the user-specified vector index and
	the second (column) index as the block index.
    */
    ubmatrix vals;
    
    /// The current rolling average
    ubvector current;

    /// The size of the vector
    size_t nvec;

  public:
    
    expval_vector();

    /** \brief Create for a vector of size \c n with \c n_blocks
	blocks and \c n_per_block points block
    */
    expval_vector(size_t n, size_t n_blocks=1, size_t n_per_block=1);

    virtual ~expval_vector();

    /// Copy constructor
    expval_vector(const expval_vector &ev);

    /// Copy constructor
    expval_vector &operator=(const expval_vector &ev);

    /** \brief Set for a vector of size \c n with \c n_blocks blocks
	and \c n_per_block points block

	\comment
	This is named differently from expval_base::set_blocks() 
	because of hiding overloaded virtual function warnings.
	\endcomment
    */
    virtual void set_blocks_vec(size_t n, size_t n_blocks, size_t n_per_block);

    /** \brief Free allocated data (but do not change the current values
	of \c n_blocks or \c n_per_block)
    */
    virtual void free();

    /// Add measurement of value \c val
    template<class vec_t> void add(vec_t &val) {

      // If all blocks are full
      if (iblock==nblocks) {

	for(size_t iv=0;iv<nvec;iv++) {
	  if (current[iv]!=0.0 || i!=0) {
	    O2SCL_ERR2("Current or 'i' nonzero with full blocks in ",
		       "expval_vector::add()",exc_esanity);
	  }
	}
    
	// Double up the data
	for(size_t iv=0;iv<nvec;iv++) {
	  for(size_t j=0;j<nblocks/2;j++) {
	    vals(iv,j)=(vals(iv,2*j)+vals(iv,2*j+1))/2.0;
	  }
	}
	// If the number of blocks is even
	if (nblocks%2==0) {
	  // Just leave current as is and clear out last half of 'vals'
	  i=0;
	  iblock=nblocks/2;
	} else {
	  // Take the odd block from vals and move it to current
	  for(size_t iv=0;iv<nvec;iv++) {
	    current[iv]=vals(iv,nblocks-1);
	  }
	  i=nperblock;
	  iblock=nblocks/2;
	}
	for(size_t iv=0;iv<nvec;iv++) {
	  for(size_t j=nblocks/2;j<nblocks;j++) {
	    vals(iv,j)=0.0;
	  }
	}
	// Double nperblock
	nperblock*=2;

      }

      // Keep track of the rolling average and increment the index
      for(size_t iv=0;iv<nvec;iv++) {
	current[iv]+=(val[iv]-current[iv])/((double)(i+1));
      }
      i++;

      // If the block is full
      if (i==nperblock) {
    
	// Store in vals and clear out current
	for(size_t iv=0;iv<nvec;iv++) {
	  vals(iv,iblock)=current[iv];
	  current[iv]=0.0;
	}
	iblock++;
	i=0;
    
      }

      return;
    }

    /// \name Report statistics
    //@{
    /** \brief Report current average, standard deviation, and 
	the error in the average and include block information

	\future This can't be const because of ubmatrix_row,
	but should be made const later.
    */
    template<class vec_t, class vec2_t, class vec3_t> 
      void current_avg_stats(vec_t &avg, vec2_t &std_dev, vec3_t &avg_err, 
			     size_t &m_block, size_t &m_per_block) {

      for(size_t k=0;k<nvec;k++) {

	// Only one block that is partially full
	if (iblock==0) {
	
	  if (i==0) {
	    O2SCL_ERR("No data in expval_scalar::current_avg_stats().",
		      exc_efailed);
	  } else {
	    m_block=1;
	    m_per_block=i;
	    avg[k]=current[k];
	    std_dev[k]=0.0;
	    avg_err[k]=0.0;
	  }
	
	} else if (iblock==1) {
	  // We're blocking, but have only one full block so far
	  m_block=1;
	  m_per_block=nperblock;
	  avg[k]=vals(k,0);
	  std_dev[k]=0.0;
	  avg_err[k]=0.0;
	
	} else if (iblock<=nblocks) {
	
	  // Report the average and std. dev.
	  // for all the blocks which have been finished
	  m_block=iblock;
	  m_per_block=nperblock;
	  typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
	  ubmatrix_row row(vals,k);
	  avg[k]=vector_mean(iblock,row);
	  std_dev[k]=vector_stddev(iblock,row);
	  avg_err[k]=std_dev[k]/sqrt(((double)iblock));
	
	}

      }

      return;
    }

    /** \brief Report current average, standard deviation, and 
	the error in the average

	\future This can't be const because of ubmatrix_row in
	current_avg_stats(), but should be made const later.
    */
    template<class vec_t, class vec2_t, class vec3_t> 
      void current_avg(vec_t &avg, vec2_t &std_dev, vec3_t &avg_err) {
      size_t m_per_block, m_block;
      return current_avg_stats(avg,std_dev,avg_err,m_block,m_per_block);
    }

    /** \brief Report average, standard deviation, and 
	the error in the average assuming a new block size
    */
    template<class vec_t, class vec2_t, class vec3_t> 
      void reblock_avg_stats(size_t new_blocks, vec_t &avg, vec2_t &std_dev, 
			     vec3_t &avg_err, size_t &m_per_block) const {

      if (new_blocks==0) {
	O2SCL_ERR2("Requested zero blocks in ",
		   "expval_vector::reblock_avg_stats().",exc_einval);
      }
  
      ubmatrix dat(nvec,new_blocks);
      for(size_t ii=0;ii<nvec;ii++) {
	for(size_t j=0;j<new_blocks;j++) {
	  dat(ii,j)=0.0;
	}
      }
  
      // The ratio of the old to new block size
      size_t fact=iblock/new_blocks;
      if (fact==0) {
	O2SCL_ERR2("Not enough data for reblocking ",
		   "in expval_vector::reblock_avg_stats().",exc_einval);
      }
      for(size_t ik=0;ik<nvec;ik++) {
	size_t iblock2=0;
	// Compute the sum
	for(size_t k=0;k<new_blocks;k++) {
	  for(size_t j=0;j<fact;j++) {
	    dat(ik,k)+=vals(ik,iblock2);
	    iblock2++;
	  }
	  // Divide to get averages
	  dat(ik,k)/=((double)fact);
	}
      }
      for(size_t ik=0;ik<nvec;ik++) {
	typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
	ubmatrix_row row(dat,ik);
	// Compute average
	avg[ik]=vector_mean(new_blocks,row);
	// Compute std. dev. and avg. err. if available
	if (new_blocks>1) {
	  std_dev[ik]=vector_stddev(new_blocks,row);
	  avg_err[ik]=std_dev[ik]/sqrt(((double)new_blocks));
	} else {
	  std_dev[ik]=0.0;
	  avg_err[ik]=0.0;
	}
      }
      // Compute m_per_block
      m_per_block=fact*nperblock;
  
      return;
    }
    
    /** \brief Report average, standard deviation, and 
	the error in the average assuming a new block size
    */
    template<class vec_t, class vec2_t, class vec3_t> 
      void reblock_avg(size_t new_blocks, vec_t &avg, 
		       vec2_t &std_dev, vec3_t &avg_err) const {
      size_t m_per_block;
      return reblock_avg_stats(new_blocks,avg,std_dev,avg_err,
			       m_per_block);
    }

    //@}

    /// Return the current data for all blocks
    const ubmatrix &get_data() const;

    friend void o2scl_hdf::hdf_output
      (o2scl_hdf::hdf_file &hf, expval_vector &t, std::string name);
    
    friend void o2scl_hdf::hdf_input_n
      (o2scl_hdf::hdf_file &hf, expval_vector &t, std::string &name);

  };

  /** \brief Matrix expectation value

      See \ref expval_base for some general notes on 
      this and related classes. 

      This is a similar to \ref expval_scalar, except that it allows
      updating and statistics for a set of matrices en masse. The data
      is stored internally in ublas matrix and \ref tensor3 objects,
      but the public member functions operate with template types
      which are compatible with any vector class which provides
      <tt>double &operator[]</tt>. It is assumed that each
      call to \ref add() contains a new measurement for all of
      the matrix entries. 
  */
  class expval_matrix : public expval_base {

  public:    

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::vector_slice<ubvector> ubvector_slice;
    typedef boost::numeric::ublas::slice slice;

  protected:
    
    /** \brief The average for each block
     */
    tensor3<> vals;
    
    /// The current rolling average
    ubmatrix current;

    /// The number of rows (zero for an empty expval_matrix object)
    size_t nr;

    /// The number of columns (zero for an empty expval_matrix object)
    size_t nc;

  public:

    expval_matrix();
    
    /** \brief Create for a vector of size \c n with \c n_blocks
	blocks and \c n_per_block points block
    */
    expval_matrix(size_t rows, size_t cols, size_t n_blocks=1, 
	      size_t n_per_block=1);

    virtual ~expval_matrix();

    /// Copy constructor
    expval_matrix(const expval_matrix &ev);

    /// Copy constructor
    expval_matrix &operator=(const expval_matrix &ev);

    /** \brief Set for a matrix with \c n_blocks blocks and \c
	n_per_block points block

	\comment
	This is named differently from expval_base::set_blocks() 
	because of hiding overloaded virtual function warnings.
	\endcomment
    */
    virtual void set_blocks_mat(size_t rows, size_t cols, 
				size_t n_blocks, size_t n_per_block);

    /** \brief Free allocated data (but do not change the current values
	of \c n_blocks or \c n_per_block)
    */
    virtual void free();

    /// Add measurement of value \c val
    template<class mat_t> void add(mat_t &val) {

      // If all blocks are full
      if (iblock==nblocks) {

	for(size_t iv=0;iv<nr;iv++) {
	  for(size_t jv=0;jv<nc;jv++) {
	    if (current(iv,jv)!=0.0 || i!=0) {
	      O2SCL_ERR2("Current or 'i' nonzero with full blocks in ",
			 "expval_matrix::add()",exc_esanity);
	    }
	  }
	}
    
	// Double up the data
	for(size_t iv=0;iv<nr;iv++) {
	  for(size_t jv=0;jv<nc;jv++) {
	    for(size_t j=0;j<nblocks/2;j++) {
	      vals.get(iv,jv,j)=(vals.get(iv,jv,2*j)+
				 vals.get(iv,jv,2*j+1))/2.0;
	    }
	  }
	}
	// If the number of blocks is even
	if (nblocks%2==0) {
	  // Just leave current as is and clear out last half of 'vals'
	  i=0;
	  iblock=nblocks/2;
	} else {
	  // Take the odd block from vals and move it to current
	  for(size_t iv=0;iv<nr;iv++) {
	    for(size_t jv=0;jv<nc;jv++) {
	      current(iv,jv)=vals.get(iv,jv,nblocks-1);
	    }
	  }
	  i=nperblock;
	  iblock=nblocks/2;
	}
	for(size_t iv=0;iv<nr;iv++) {
	  for(size_t jv=0;jv<nc;jv++) {
	    for(size_t j=nblocks/2;j<nblocks;j++) {
	      vals.get(iv,jv,j)=0.0;
	    }
	  }
	}
	// Double nperblock
	nperblock*=2;

      }

      // Keep track of the rolling average and increment the index
      for(size_t iv=0;iv<nr;iv++) {
	for(size_t jv=0;jv<nc;jv++) {
	  current(iv,jv)+=(val(iv,jv)-current(iv,jv))/((double)(i+1));
	}
      }
      i++;

      // If the block is full
      if (i==nperblock) {
    
	// Store in vals and clear out current
	for(size_t iv=0;iv<nr;iv++) {
	  for(size_t jv=0;jv<nc;jv++) {
	    vals.get(iv,jv,iblock)=current(iv,jv);
	    current(iv,jv)=0.0;
	  }
	}
	iblock++;
	i=0;
    
      }

      return;
    }

    /// \name Report statistics
    //@{
    /** \brief Report current average, standard deviation, and 
	the error in the average and include block information

	\future This should be made const.
	\future Avoid the copy associated with vector_slice().
    */
    template<class mat_t, class mat2_t, class mat3_t> 
      void current_avg_stats(mat_t &avg, mat2_t &std_dev, 
			     mat3_t &avg_err, size_t &m_block,
			     size_t &m_per_block) {

      for(size_t j=0;j<nr;j++) {
	for(size_t k=0;k<nc;k++) {
	  
	  // Only one block that is partially full
	  if (iblock==0) {
	    
	    if (i==0) {
	      O2SCL_ERR("No data in expval_scalar::current_avg_stats().",
			exc_efailed);
	    } else {
	      m_block=1;
	      m_per_block=i;
	      avg(j,k)=current(j,k);
	      std_dev(j,k)=0.0;
	      avg_err(j,k)=0.0;
	    }
	    
	  } else if (iblock==1) {
	    
	    // We're blocking, but have only one full block so far
	    m_block=1;
	    m_per_block=nperblock;
	    avg(j,k)=vals.get(j,k,0);
	    std_dev(j,k)=0.0;
	    avg_err(j,k)=0.0;
	    
	  } else if (iblock<=nblocks) {
	    
	    // Report the average and std. dev.
	    // for all the blocks which have been finished
	    m_block=iblock;
	    m_per_block=nperblock;
	    
	    // Create a vector from vals which leaves the first
	    // index free

	    // The function vector_slice() doesn't quite work this way
	    // with the new tensor class based on std::vector<double>
	    // objects, so we make copy for now as a replacement.

	    //size_t dims[3]={j,k,0};
	    //ubvector_slice col=vals.vector_slice(2,dims);

	    std::vector<double> col(iblock);
	    for (size_t ik=0;ik<iblock;ik++) {
	      col[ik]=vals.get(j,k,ik);
	    }

	    // Obtain stats from that vector
	    avg(j,k)=vector_mean(iblock,col);
	    std_dev(j,k)=vector_stddev(iblock,col);
	    avg_err(j,k)=std_dev(j,k)/sqrt(((double)iblock));

	  }

	}
      }
      
      return;
    }

    /** \brief Report current average, standard deviation, and 
	the error in the average

	\future This should be made const.
    */
    template<class mat_t, class mat2_t, class mat3_t> 
      void current_avg(mat_t &avg, mat2_t &std_dev, mat3_t &avg_err) {
      size_t m_per_block, m_block;
      return current_avg_stats(avg,std_dev,avg_err,m_block,m_per_block);
    }

    /** \brief Report average, standard deviation, and 
	the error in the average assuming a new block size
    */
    template<class mat_t, class mat2_t, class mat3_t> 
      void reblock_avg_stats(size_t new_blocks, mat_t &avg, 
			     mat2_t &std_dev, mat3_t &avg_err,
			     size_t &m_per_block) const {

      if (new_blocks==0) {
	O2SCL_ERR2("Requested zero blocks in ",
		   "expval_vector::reblock_avg_stats().",exc_einval);
      }
  
      tensor3<double> dat(nr,nc,new_blocks);
      dat.set_all(0.0);
  
      // The ratio of the old to new block size
      size_t fact=iblock/new_blocks;
      if (fact==0) {
	O2SCL_ERR2("Not enough data for reblocking ",
		   "in expval_vector::reblock_avg_stats().",exc_einval);
      }
      for(size_t ik=0;ik<nr;ik++) {
	for(size_t jk=0;jk<nc;jk++) {
	  size_t iblock2=0;
	  // Compute the sum
	  for(size_t k=0;k<new_blocks;k++) {
	    for(size_t j=0;j<fact;j++) {
	      dat.get(ik,jk,k)+=vals.get(ik,jk,iblock2);
	      iblock2++;
	    }
	    // Divide to get averages
	    dat.get(ik,jk,k)/=((double)fact);
	  }
	}
      }
      for(size_t ik=0;ik<nr;ik++) {
	for(size_t jk=0;jk<nc;jk++) {

	  //size_t dim[3]={ik,jk,0};
	  //ubvector_slice vec=dat.vector_slice(2,dim);
	  
	  std::vector<double> vec(new_blocks);
	  for (size_t ii=0;ii<new_blocks;ii++) {
	    vec[ii]=dat.get(ik,jk,ii);
	  }

	  // Compute average
	  avg(ik,jk)=vector_mean(new_blocks,vec);
	  // Compute std. dev. and avg. err. if available
	  if (new_blocks>1) {
	    std_dev(ik,jk)=vector_stddev(new_blocks,vec);
	    avg_err(ik,jk)=std_dev(ik,jk)/sqrt(((double)new_blocks));
	  } else {
	    std_dev(ik,jk)=0.0;
	    avg_err(ik,jk)=0.0;
	  }
	}
      }
      // Compute m_per_block
      m_per_block=fact*nperblock;
  
      return;
    }
    
    /** \brief Report average, standard deviation, and 
	the error in the average assuming a new block size
    */
    template<class mat_t, class mat2_t, class mat3_t> 
      void reblock_avg(size_t new_blocks, mat_t &avg, 
		       mat2_t &std_dev, mat3_t &avg_err) const {
      size_t m_per_block;
      return reblock_avg_stats(new_blocks,avg,std_dev,avg_err,
			       m_per_block);
    }
    //@}

    /// Return the current data for all blocks
    const tensor3<double> &get_data() const;

    friend void o2scl_hdf::hdf_output
      (o2scl_hdf::hdf_file &hf, expval_matrix &t, std::string name);
    
    friend void o2scl_hdf::hdf_input_n(o2scl_hdf::hdf_file &hf, expval_matrix &t, 
				     std::string &name);

  };

}

#endif
