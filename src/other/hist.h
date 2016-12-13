/*
  -------------------------------------------------------------------
  
  Copyright (C) 2010-2016, Andrew W. Steiner
  
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
#ifndef O2SCL_HIST_H
#define O2SCL_HIST_H

/** \file hist.h
    \brief File defining \ref o2scl::hist
*/
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/convert_units.h>
#include <o2scl/interp.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/table.h>

// Forward definition of the hist class for HDF I/O
namespace o2scl {
  class hist;
}

// Forward definition of HDF I/O to extend friendship in hist
namespace o2scl_hdf { 
  class hdf_file; 
  void hdf_input(hdf_file &hf, o2scl::hist &t, std::string name);
  void hdf_output(hdf_file &hf, o2scl::hist &t, std::string name);
}

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief A one-dimensional histogram class
      
      See discussion in the User's guide in the \ref hist_section 
      section.

      One may set the histogram bins using \ref set_bin_edges() or one
      may manually set the limit of one bin using the reference
      returned by get_bin_low(), get_bin_low_i(), get_bin_high(), or
      get_bin_high_i(). Note that if one attempts to set the bins on a
      histogram where the bins have already been set, one must ensure
      that the new and old bin sets have the same size. This ensures
      that there is no ambiguity in rebinning the data and also 
      prevents accidental data loss. One may set the bin edges
      either with a generic vector, or as a \ref uniform_grid object.

      To save space, representative vectors are not allocated until
      they are used.

      \note In order to ensure the histogram does not employ
      user-specified representative values that are not defined, the
      function \ref set_rep_mode() does not allow one to change the
      mode to \ref rmode_user directly. Instead, use \ref set_reps()
      which automatically sets the mode to \ref rmode_user and allows
      the user to specify the representatives.

      \note If the user changes the bin edges and the histogram is in
      mode \ref rmode_user, the bin weights will not be modified and
      the same representative values will be assumed for the new bin
      edges.

      \hline
      
      \todo Check implementation of <tt>hist::extend_lhs</tt>.
      \todo More testing.

      \future 
      - Add a counter which counts the number of calls to update()?
      - Add conversions back and forth from GSL histograms
      - Create extend_lhs too?
      - Would be nice not to have to create a new \ref
      o2scl::search_vec object in \ref o2scl::hist::get_bin_index()
      (make a search_vec data member?)
      - Consider adding the analogs of the GSL histogram
      sampling functions (separate class?)
      - Add a function which computes the bin sizes?
      - Allow rebinning?
      - Add histograms of float and integer values
      - Allow addition and other operations for two histograms.
      - Make the interpolation functions \c const (this is a bit
      complicated because of \ref o2scl::hist::set_reps_auto() ).

      \hline 

      Internally, none of the vectors should have memory allocated for
      them when hsize is zero, and the vector sizes should match the
      histogram size. These and other checks are performed by \ref
      is_valid() . Also, the function \ref set_reps_auto() should not
      be called when mode is \ref rmode_user.
  */
  class hist {

  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;

  protected:

    /// Bin locations (N+1)
    ubvector ubin;
    
    /// Bin contents (N)
    ubvector uwgt;

    /// Bin representative values (N)
    ubvector urep;
    
    /// User-defined representative values (N)
    ubvector user_rep;
    
    /// Number of bins
    size_t hsize;

    /// Representative mode
    size_t rmode;

    /// Interpolation type
    size_t itype;

    /** \brief Set the representative array according to current 
	rmode (if not in user rep mode)
     */
    void set_reps_auto();

    /// Interpolation typedef
    typedef interp_vec<ubvector> interp_t;

    /** \brief Allocate vectors for a histogram of size \c n

	This function also sets all the weights to zero.
     */
    void allocate(size_t n);

  public:

    /// Create an empty histogram
    hist();

    ~hist();

    /// Copy constructor
    hist(const hist &h);

    /// Copy constructor
    hist &operator=(const hist &h);

    /// Create from a vectors of data
    template<class vec_t> hist(size_t nv, const vec_t &v, size_t n_bins) {
			       
      itype=1;
      rmode=rmode_avg;
      hsize=0;
      extend_lhs=true;
      extend_rhs=true;
      
      double min, max;
      o2scl::vector_minmax_value(nv,v,min,max);
      uniform_grid<double> ug=uniform_grid_end<double>(min,max,n_bins);
      set_bin_edges(ug);

      for(size_t i=0;i<nv;i++) {
	update(v[i]);
      }
      return;
    }
    
    /// Create from a vectors of data
    template<class vec_t, class vec2_t>
      hist(size_t nv, const vec_t &v,
	   const vec2_t &v2, size_t n_bins) {
      
      itype=1;
      rmode=rmode_avg;
      hsize=0;
      extend_lhs=true;
      extend_rhs=true;
      
      double min, max;
      o2scl::vector_minmax_value(nv,v,min,max);
      uniform_grid<double> ug=uniform_grid_end<double>(min,max,n_bins);
      set_bin_edges(ug);

      for(size_t i=0;i<nv;i++) {
	update(v[i],v2[i]);
      }
      return;
    }
    
    /// Create from vectors of data
    template<class vec_t, class vec2_t> hist(const vec_t &v, size_t n_bins) {
      size_t nv=v.size();
      hist(nv,v,n_bins);
      return;
    }

    /// Create from vectors of data
    template<class vec_t, class vec2_t> hist
      (const vec_t &v, const vec2_t &v2, size_t n_bins) {
					     
      size_t nv=v.size();
      hist(nv,v,v2,n_bins);
      return;
    }

    // Create a histogram from a \ref o2scl::table object
    void from_table(o2scl::table<> &t, std::string colx, 
		    size_t n_bins) {
      *this=hist(t.get_nlines(),t.get_column(colx),n_bins);
      return;
    }
    
    // Create a histogram from a \ref o2scl::table object
    void from_table(o2scl::table<> &t, std::string colx, std::string coly,
		    size_t n_bins) {
      *this=hist(t.get_nlines(),t.get_column(colx),t.get_column(coly),
		 n_bins);
      return;
    }
    
    /// The histogram size
    size_t size() const {
      return hsize;
    }

    /** \brief If true, allow abcissae beyond the last bin (default false)

	If this is true, the histogram will allow data with
	corresponding to bins larger than the largest bin 
	(for increasing bin settings) or smaller than the 
	smallest bin (for decreasing bin settings). 
    */
    bool extend_rhs;

    /** \brief If true, allow abcissae before the first bin (default false)
     */
    bool extend_lhs;

    /// \name Initial bin setup
    //@{
    /** \brief Set bins from a \ref uniform_grid object

	If the current histogram is not empty, then the 
	number of bins reported by \ref uniform_grid<>::get_nbins()
	should be equal to the current histogram size so that the
	number of bins is equal and we can use the same weights.

	If either the histogram is empty, or the current
	representative mode is not \ref rmode_user, then the
	representative mode is automatically set to \ref rmode_avg (or
	\ref rmode_gmean if \ref uniform_grid::is_log() returns \c
	true ) .
    */
    void set_bin_edges(uniform_grid<double> g);

    /** \brief Set the bins from a vector

	The parameter \c n is the size of the vector, equal to the
	number of edges (always one more than the number of bins). If
	the current histogram is not empty, then \c n should be equal
	one larger to the size reported by \ref size() so that the
	number of bins is equal and we can use the same weights.
    */
    template<class vec_t> void set_bin_edges(size_t n, const vec_t &v) {
      if (n!=hsize+1) {
	if (hsize!=0) {
	  O2SCL_ERR2("Requested binning change in non-empty ",
		     "histogram in hist::set_bin_edges().",exc_efailed);
	}
	allocate(n-1);
      }
      for(size_t i=0;i<n;i++) ubin[i]=v[i];
      // Reset internal reps
      if (urep.size()>0) urep.clear();
      return;
    }
    //@}

    /// \name Weight functions
    //@{
    /// Increment bin for \c x by value \c val
    void update(double x, double val=1.0);

    /// Increment bin with index \c i by value \c val
    void update_i(size_t i, double val=1.0) {
      uwgt[i]+=val;
      return;
    }

    /// Return contents of bin with index \c i
    const double &get_wgt_i(size_t i) const;

    /// Return contents of bin with index \c i
    double &get_wgt_i(size_t i);

    /// Return contents of bin for \c x
    const double &get_wgt(double x) const {
      return get_wgt_i(get_bin_index(x));
    }

    /// Return contents of bin for \c x
    double &get_wgt(double x) {
      return get_wgt_i(get_bin_index(x));
    }

    /// Set contents of bin with index \c i to value \c val
    void set_wgt_i(size_t i, double val);

    /// Set contents of bin for \c x to value \c val
    void set_wgt(double x, double val) {
      set_wgt_i(get_bin_index(x),val);
      return;
    }
    
    /// Get a reference to the full y vector
    const ubvector &get_wgts() const {
      return uwgt;
    }

    /// Get a reference to the weight for the bin at index \c i
    const double &operator[](size_t i) const {
      return uwgt[i];
    }
    
    /// Get a reference to the weight for the bin at index \c i
    double &operator[](size_t i) {
      return uwgt[i];
    }
    //@}

    /// \name Bin manipulation
    //@{
    /** \brief Get the index of the bin which holds \c x

	Always returns a value between 0 and size() (inclusive)
     */
    size_t get_bin_index(double x) const;

    /// Get the lower edge of bin of index \c i
    double &get_bin_low_i(size_t i);

    /// Get the lower edge of bin of index \c i
    const double &get_bin_low_i(size_t i) const;

    /// Get the upper edge of bin of index \c i
    double &get_bin_high_i(size_t i);
    
    /// Get the upper edge of bin of index \c i
    const double &get_bin_high_i(size_t i) const;

    /// Get the lower edge of bin of index \c i
    double &get_bin_low(double x) {
      return get_bin_low_i(get_bin_index(x));
    }

    /// Get the lower edge of bin of index \c i
    const double &get_bin_low(double x) const {
      return get_bin_low_i(get_bin_index(x));
    }

    /// Get the upper edge of bin of index \c i
    double &get_bin_high(double x) {
      return get_bin_high_i(get_bin_index(x));
    }
    
    /// Get the upper edge of bin of index \c i
    const double &get_bin_high(double x) const {
      return get_bin_high_i(get_bin_index(x));
    }

    /// Get a reference to the full vector of bin specifications
    const ubvector &get_bins() const {
      return ubin;
    }
    //@}

    /// \name Max and min functions
    //@{
    /** \brief Get maximum weight
     */
    double get_max_wgt() const;

    /** \brief Get the bin index of the maximum weight
     */
    size_t get_max_index() const;

    /** \brief Get the representative for the bin with maximum weight
     */
    double get_max_rep();

    /** \brief Get minimum weight
     */
    double get_min_wgt() const;

    /** \brief Get the bin index of the minimum weight
     */
    size_t get_min_index() const;

    /** \brief Get the representative for the bin with minimum weight
     */
    double get_min_rep();
    //@}

    /// \name Delete functions
    //@{
    /// Clear the data, but leave the bins as is
    void clear_wgts();

    /// Clear the entire histogram
    void clear();
    //@}

    /// \name Representative modes (default is rmode_avg)
    // Using \c rmode_avg in documentation doesn't work.
    //@{
    /// Average lower and upper edge
    static const size_t rmode_avg=0;
    /// Use user-specified representative
    static const size_t rmode_user=1;
    /// Use lower edge
    static const size_t rmode_low=2;
    /// Use upper edge
    static const size_t rmode_high=3;
    /// Use the geometric mean of the lower and upper edges
    static const size_t rmode_gmean=4;
    //@}
    
    /// \name Representative functions
    //@{
    /// Set the representative x-values for each bin
    template<class vec_t> void set_reps(size_t n, const vec_t &v) {
      if (n!=hsize) {
	std::string s="Expected a vector of size "+itos(hsize)+
	  " and got a vector of size "+itos(n)+" in hist::set_reps().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
      rmode=rmode_user;
      if (user_rep.size()>0) user_rep.clear();
      user_rep.resize(n);
      for(size_t i=0;i<n;i++) user_rep[i]=v[i];
      return;
    }

    /** \brief Set mode used to compute bin representatives

	Acceptable inputs are \ref rmode_avg, \ref rmode_low,
	\ref rmode_high, and \ref rmode_gmean .
    */
    void set_rep_mode(size_t mode);
    
    /// Get mode used to compute bin representatives
    size_t get_rep_mode() const {
      return rmode;
    }

    /** \brief Return the representative of bin of index \c i

	Note that this function returns a value and not a reference.
	This is because we can't return a reference to the internally
	computed representatives, since they don't always exist.
    */
    double get_rep_i(size_t i);

    /// Return the representative of bin containing \c x
    double get_rep(double x) {
      return get_rep_i(get_bin_index(x));
    }

    /** \brief Create a vector filled with the representatives for 
	each bin
     */
    template<class resize_vec_t> void create_rep_vec(resize_vec_t &v) {
      v.resize(hsize);
      for(size_t i=0;i<hsize;i++) {
	v[i]=get_rep_i(i);
      }
      return;
    }
    //@}

    /// \name Evaluation and interpolation functions
    //@{
    /// Return the value of the function at \c x
    double operator()(double x);

    /// Return the value of the function at \c x
    double interp(double x);

    /// Return the derivative of the function at \c x
    double deriv(double x);

    /// Return the second derivative of the function at \c x
    double deriv2(double x);

    /// Return the integral of the function between \c x and \c y
    double integ(double x, double y);

    /// Set the interpolation type
    void set_interp_type(size_t interp_type);
    //@}

    /// \name Other functions
    //@{
    double bin_sum();
    
    /** \brief Renormalize the weights to fix the integral
	
	This computes the integral using \ref integ() and so the
	action of this function depends on the interpolation type.
	If the histogram is empty, an exception is thrown. 
    */
    void normalize(double new_sum);

    /// Internal consistency check
    void is_valid() const;
    
    /** \brief Copy histogram data to a table

	First, if the table \c t has less rows than the histogram has
	bins, then new rows are added to the table and values in the
	new rows of the current columns are set to zero. Second, this
	creates new columns in the table named \c reps, \c
	lower_edges, \c upper_edges, and \c weights . Finally,
	the histogram data is copied to the four new columns. 
    */
    void copy_to_table(table<> &t, std::string reps, std::string lower_edges, 
		       std::string upper_edges, std::string weights);
    //@}

    // Allow HDF I/O function to access hist data
    friend void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf, o2scl::hist &h, 
				      std::string name);
    
    // Allow HDF I/O function to access hist data
    friend void o2scl_hdf::hdf_input(o2scl_hdf::hdf_file &hf, o2scl::hist &h, 
				     std::string name);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
