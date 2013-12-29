/*
  -------------------------------------------------------------------
  
  Copyright (C) 2011-2013 Andrew W. Steiner
  
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
#ifndef HIST_EV_H
#define HIST_EV_H

#include <o2scl/tensor.h>
#include <o2scl/expval.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/hist.h>
#include <o2scl/hist_2d.h>
#if O2SCL_HDF_SVAR
#include <o2scl/hdf_file.h>
#endif

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A set expectation values for histogram bins

      This class is experimental.

      See \ref expval_base for some general notes on this and related
      classes.

      This class computes expectation values of bins in a histogram.
      It is most useful in cases where one does not know a priori how
      many measurements one is going to get for each bin. The class
      automatically arranges each bin into blocks, so that the
      standard deviation and error and the average can be computed
      even when not all bins have been filled with the same number of
      measurements.

      \todo Test set_grid_blocks() function.

      \todo Create copy constructors as in the scalar_ev class.
  */
  class hist_ev : public expval_base {
    
  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    typedef boost::numeric::ublas::vector<int> ubvector_int;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix<size_t> ubmatrix_size_t;
    typedef boost::numeric::ublas::matrix<int> ubmatrix_int;

  protected:

    /// Running average for each block and each bin
    ubmatrix vals;

    /// The value of \c iblock for each bin
    ubvector_size_t iblock_bins;

    /// The value of \c nperblock for each bin
    ubvector_size_t nperblock_bins;

    /// The value of \c i for each bin
    ubvector_size_t i_bins;

    /** \brief This should always be the same as the size as the
	histogram
    */
    size_t hsize;
    
    /** \brief The associated histogram

	This object is only currently used internally for binning. No
	data is added to it.
    */
    hist h;
    
  public:

    /// Create an empty object
    hist_ev();

    /// Create a histogram expectation value 
    hist_ev(uniform_grid<double> g, size_t n_blocks, size_t n_per_block);
    
    /// Set the histogram grid and the number of blocks
    void set_grid_blocks(uniform_grid<double> g, size_t n_blocks,
			 size_t n_per_block);

    /// Add measurement of value \c val at location \c x
    virtual void add(double x, double val=1.0);
    
    /** \brief Report current average, standard deviation, and 
	the error in the average

	This function deallocates any space already allocated for the
	vector parameters and reallocates space as needed. Information
	previously stored in these vectors will be lost.
    */
    virtual void current_avg_stats(ubvector &reps, ubvector &avg, 
				   ubvector &std_dev, ubvector &avg_err, 
				   ubvector_int &m_block,
				   ubvector_int &m_per_block);

    /** \brief Report current average, standard deviation, and 
	the error in the average

	This function deallocates any space already allocated for the
	vector parameters and reallocates space as needed. Information
	previously stored in these vectors will be lost.
    */
    virtual void current_avg(ubvector &reps, ubvector &avg, ubvector &std_dev, 
			     ubvector &avg_err);

#if O2SCL_HDF_SVAR

    /// \name HDF I/O functions
    //@{
    /// Output to an hdf file using group named \c name
    void hdf_output(o2scl_hdf::hdf_file &hf, std::string name);

    /// Input from an hdf file using group named \c name
    void hdf_input(o2scl_hdf::hdf_file &hf, std::string name="");
    //@}

#endif

  };

  /** \brief A two-dimensional histogram of expectation values

      See \ref expval_base for some general notes on 
      this and related classes. 

      This class is experimental.
   */
  class hist_2d_ev : public expval_base {
    
  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::vector<size_t> ubvector_size_t;
    typedef boost::numeric::ublas::vector<int> ubvector_int;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix<size_t> ubmatrix_size_t;
    typedef boost::numeric::ublas::matrix<int> ubmatrix_int;

  protected:

    /// Running average for each block and each bin
    tensor3 vals;

    /// The value of \c iblock for each bin
    ubmatrix_size_t iblock_bins;

    /// The value of \c nperblock for each bin
    ubmatrix_size_t nperblock_bins;

    /// The value of \c i for each bin
    ubmatrix_size_t i_bins;

    /// Histogram size
    size_t hsize_x;

    /// Histogram size
    size_t hsize_y;

    /// The histogram
    hist_2d h;

  public:
    
    /// Create a histogram expectation value 
    hist_2d_ev(uniform_grid<double> hxg, uniform_grid<double> hyg, 
	       size_t n_blocks, size_t n_per_block);

    /// Add measurement of value \c val at location \c x
    virtual void add(double x, double y, double val);

    /** \brief Report current average, standard deviation, and 
	the error in the average
    */
    void current_avg_stats
      (ubvector &rep_x, ubvector &rep_y, ubmatrix &avg, ubmatrix &std_dev, 
       ubmatrix &avg_err, ubmatrix_int &m_block, ubmatrix_int &m_per_block);

#if O2SCL_HDF_SVAR
#ifdef O2SCL_NEVER_DEFINED

    /// \name HDF I/O functions
    //@{
    /// Output to an hdf file using group named \c name
    void hdf_output(o2scl_hdf::hdf_file &hf, std::string name);

    /// Input from an hdf file using group named \c name
    void hdf_input(o2scl_hdf::hdf_file &hf, std::string name="");
    //@}

#endif
#endif
    
  };

}

#endif
