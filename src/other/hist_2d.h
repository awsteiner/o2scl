/*
  -------------------------------------------------------------------
  
  Copyright (C) 2010-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_HIST_2D_H
#define O2SCL_HIST_2D_H

/** \file hist_2d.h
    \brief File defining \ref o2scl::hist_2d
*/
#include <iostream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/convert_units.h>
#include <o2scl/interp.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/table3d.h>

// Forward definition of the hist_2d class for HDF I/O
namespace o2scl {
  class hist_2d;
}

// Forward definition of HDF I/O to extend friendship
namespace o2scl_hdf { 
  class hdf_file; 
  void hdf_input(hdf_file &hf, o2scl::hist_2d &t, std::string name);
  void hdf_output(hdf_file &hf, o2scl::hist_2d &t, std::string name);
}

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief A two-dimensional histogram class

      See discussion in the User's guide in the \ref hist_section 
      section.

      Typical usage begins with setting the histogram bins using 
      \ref hist_2d::set_bin_edges(). Note that if one attempts to set
      the bins on a histogram where the bins have already been set,
      one must ensure that the new and old bin settings have the same
      size (in both x and y directions). This ensures that there is no
      ambiguity in rebinning the data and also prevents accidental
      data loss. One may set the bin edges either with generic
      vectors, or with \ref uniform_grid objects.

      \note In order to ensure the histogram does not employ
      user-specified representative values that are not defined, the
      function \ref set_rep_mode() does not allow one to change the
      mode to \ref hist::rmode_user directly. Instead, use \ref
      set_reps() which automatically sets the mode to \ref
      hist::rmode_user and allows the user to specify the
      representatives.

      \comment

      This is commented out for now. The \ref twod_intp
      object stores a reference to xrep and yrep, and thus 
       can't be used since xrep and yrep don't always exist.

      Interpolation for \ref hist_2d objects is performed by creating
      a \ref twod_intp object from the internal histogram data. This
      is done using \ref o2scl::hist_2d::setup_interp() .

      \endcomment

      Internally, either \ref hsize_x and \ref hsize_y should
      both be zero or both be non-zero. 

      \future Write a function to create a 1-d histogram 
      from a 2-d histogram either by selecting one bin
      in one axis or by marginalizing over one direction.

      \future Note that here, there is a conflict between implementing
      operator(size_t,size_t) to give matrix indexing of the histogram
      weights, and operator(double,double) to implement
      two-dimensional interpolation using the weights and the
      representatives. Currently neither is implemented, but maybe
      both should be implemented instead?
   */
  class hist_2d {

  public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    
  protected:

    /// Bin locations (Nx+1)
    ubvector xa;

    /// Bin locations (Ny+1)
    ubvector ya;

    /// Values (Nx,Ny)
    ubmatrix wgt;

    /// "Central" values for x-axis (N)
    ubvector xrep;

    /// "Central" values for y-axis (N)
    ubvector yrep;
    
    /// User-defined central values for x-axis (N)
    ubvector user_xrep;
    
    /// User-defined central values for y-axis (N)
    ubvector user_yrep;
    
    /// Number of x-bins
    size_t hsize_x;

    /// Number of y-bins
    size_t hsize_y;

    /// Rep mode for x
    size_t xrmode;

    /// Rep mode for y
    size_t yrmode;

    /// Interpolation type
    size_t itype;

    /** \brief Allocate for a histogram of size \c nx, \c ny
	
	This function also sets all the weights to zero.
    */
    void allocate(size_t nx, size_t ny);

    /** \brief An internal function to automatically set
	\ref xrep and \ref yrep
    */
    void set_reps_auto();

  public:

    hist_2d();

    virtual ~hist_2d();

    /// Copy constructor
    hist_2d(const hist_2d &h);

    /// Copy from <tt>operator=()</tt>
    hist_2d &operator=(const hist_2d &h);
    
    /** \brief If true, allow abcissa larger than largest bin limit
	to correspond to the highest bin (default false).
    */
    bool extend_rhs;
    
    /** \brief If true, allow abcissa smaller than smallest bin
	limit to correspond to the lowest bin (default false).
    */
    bool extend_lhs;

    /// \name Initial bin setup
    //@{
    /// Set the bins from two \ref uniform_grid objects
    void set_bin_edges(uniform_grid<double> gx, uniform_grid<double> gy);

    /// Set the bins from a vector
    template<class vec_t> void set_bin_edges(size_t nx, vec_t &vx,
					    size_t ny, vec_t &vy) {
      if (nx!=hsize_x+1 || ny!=hsize_y+1) {
	if (hsize_x!=0 || hsize_y!=0) {
	  O2SCL_ERR2("Requested binning change in non-empty ",
			 "histogram in hist_2d::set_bin_edges().",
			 exc_efailed);
	}
	allocate(nx-1,ny-1);
      }
      for(size_t i=0;i<nx;i++) xa[i]=vx[i];
      for(size_t i=0;i<ny;i++) ya[i]=vy[i];
      // Reset internal reps
      if (xrep.size()>0) xrep.resize(0);
      if (yrep.size()>0) yrep.resize(0);
      return;
    }
    //@}

    /// \name Weight functions
    //@{
    /// Increment bin at <tt>(i,j)</tt> by value \c val
    void update_i(size_t i, size_t j, double val=1.0) {
      wgt(i,j)+=val;
      return;
    }

    /// Increment bin for \c x by value \c val
    void update(double x, double y, double val=1.0) {
      size_t i, j;
      get_bin_indices(x,y,i,j);
      update_i(i,j,val);
      return;
    }

    /// Return contents of bin at <tt>(i,j)</tt>
    const double &get_wgt_i(size_t i, size_t j) const;

    /// Return contents of bin for \c x
    const double &get_wgt(double x, double y) const {
      size_t i, j;
      get_bin_indices(x,y,i,j);
      return get_wgt_i(i,j);
    }

    /// Return contents of bin at <tt>(i,j)</tt>
    double &get_wgt_i(size_t i, size_t j);

    /// Return contents of bin for \c x
    double &get_wgt(double x, double y) {
      size_t i, j;
      get_bin_indices(x,y,i,j);
      return get_wgt_i(i,j);
    }

    /// Set contents of bin at <tt>(i,j)</tt> to value \c val
    void set_wgt_i(size_t i, size_t j, double val);
    
    /// Set contents of bin for \c x to value \c val
    void set_wgt(double x, double y, double val) {
      size_t i, j;
      get_bin_indices(x,y,i,j);
      set_wgt_i(i,j,val);
      return;
    }
    
    /// Get a const reference to the full matrix of data
    const ubmatrix &get_wgts() const {
      return wgt;
    }

    /// Get a reference to the full matrix of data
    ubmatrix &get_wgts() {
      return wgt;
    }
    //@}
    
    /// \name Delete functions
    //@{
    /// Clear the data, but leave the bins as is
    void clear_wgts();

    /// Clear the entire histogram
    void clear();
    //@}

    /// \name Bin manipulation
    //@{
    /** \brief Get the index of the bin which holds \c x and 
	the bin which holds \c y
    */
    void get_bin_indices(double x, double y, size_t &i, size_t &j) const;

    /// Get the index of the bin which holds \c x
    size_t get_x_bin_index(double x) const;
    
    /// Get the indey of the bin which holds \c y
    size_t get_y_bin_index(double y) const;

    /// Get the lower edge of bin of index \c i
    double &get_x_low_i(size_t i);

    /// Get the lower edge of bin of index \c i
    const double &get_x_low_i(size_t i) const;

    /// Get the upper edge of bin of index \c i
    double &get_x_high_i(size_t i);
    
    /// Get the upper edge of bin of index \c i
    const double &get_x_high_i(size_t i) const;

    /// Get the lower edge of bin of index \c j
    double &get_y_low_i(size_t j);

    /// Get the lower edge of bin of index \c j
    const double &get_y_low_i(size_t j) const;

    /// Get the upper edge of bin of index \c j
    double &get_y_high_i(size_t j);
    
    /// Get the upper edge of bin of index \c j
    const double &get_y_high_i(size_t j) const;
    //@}

    /// \name Rep modes (default is \c rmode_avg)
    //@{
    static const size_t rmode_avg=0;
    static const size_t rmode_user=1;
    static const size_t rmode_low=2;
    static const size_t rmode_high=3;
    static const size_t rmode_gmean=4;
    //@}
    
    /// \name Representative functions
    //@{
    /// Set the representative x-values for each bin
    template<class vec_t> void set_reps(size_t nx, vec_t &vx,
					size_t ny, vec_t &vy) {
      if (user_xrep.size()!=hsize_x || user_yrep.size()!=hsize_y) {
	std::string s="Expected vectors of size "+itos(hsize_x)+
	  ", "+itos(hsize_y)+" and got a vectors of size "+itos(nx)+
	  ", "+itos(ny)+" in hist_2d::set_reps().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
      xrmode=rmode_user;
      yrmode=rmode_user;
      if (user_xrep.size()>0) user_xrep.clear();
      if (user_yrep.size()>0) user_yrep.clear();
      user_xrep.resize(nx);
      user_yrep.resize(ny);
      for(size_t i=0;i<nx;i++) user_xrep[i]=vx[i];
      for(size_t i=0;i<ny;i++) user_yrep[i]=vy[i];
      return;
    }

    /// Set the representative x-values for each bin
    template<class vec_t> void set_x_reps(size_t nx, vec_t &vx) {
      if (hsize_x!=nx) {
	std::string s="Expected vector of size "+itos(hsize_x)+
	  " and got a vector of size "+itos(nx)+" in hist_2d::set_reps().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
      xrmode=rmode_user;
      if (user_xrep.size()>0) user_xrep.clear();
      user_xrep.resize(nx);
      for(size_t i=0;i<nx;i++) user_xrep[i]=vx[i];
      return;
    }

    /// Set the representative y-values for each bin
    template<class vec_t> void set_y_reps(size_t ny, vec_t &vy) {
      if (hsize_y!=ny) {
	std::string s="Expected vector of size "+itos(hsize_y)+
	  " and got a vector of size "+itos(ny)+" in hist_2d::set_reps().";
	O2SCL_ERR(s.c_str(),exc_einval);
      }
      yrmode=rmode_user;
      if (user_yrep.size()>0) user_yrep.clear();
      user_yrep.resize(ny);
      for(size_t i=0;i<ny;i++) user_yrep[i]=vy[i];
      return;
    }

    /// Set mode used to compute bin reps
    void set_rep_mode(size_t x_mode, size_t y_mode);

    /// Get mode used to compute bin reps
    size_t get_x_rep_mode() const {
      return xrmode;
    }

    /// Get mode used to compute bin reps
    size_t get_y_rep_mode() const {
      return yrmode;
    }

    /// Get a reference to the full vector of bin specifications
    const ubvector &get_x_bins() const {
      return xa;
    }

    /// Get a reference to the full vector of bin specifications
    const ubvector &get_y_bins() const {
      return ya;
    }

    /// Return the histogram size of the x coordinate
    size_t size_x() const {
      return hsize_x;
    }

    /// Return the histogram size of the y coordinate
    size_t size_y() const {
      return hsize_y;
    }

    /** \brief Get a reference to the user-specified reps for x coordinates

	This function will call the error handler if the x-axis
	representative mode is not \ref hist::rmode_user .

	\warning This vector reference is only valid so long as
	the representative mode is unchanged and the function
	clear() is not called. 

	This member function is used by the \o2 HDF I/O functions.
    */
    const ubvector &get_user_reps_x() const {
      if (xrmode!=rmode_user) {
	O2SCL_ERR("Not user mode in hist::get_user_reps_x().",
		  exc_efailed);
      }
      return user_xrep;
    }

    /** \brief Get a reference to the user-specified reps for y coordinates

	This function will call the error handler if the y-axis
	representative mode is not \ref hist::rmode_user .

	\warning This vector reference is only valid so long as
	the representative mode is unchanged and the function
	clear() is not called. 

	This member function is used by the \o2 HDF I/O functions.
    */	
    const ubvector &get_user_reps_y() const {
      if (yrmode!=rmode_user) {
	O2SCL_ERR("Not user mode in hist::get_user_reps_y().",
		  exc_efailed);
      }
      return user_yrep;
    }

    /** \brief Return the rep of bin of index \c i

	Note that this function returns a value and not a reference.
	This is because we can't return a reference to the internally
	computed representatives, since they don't always exist.
     */
    double get_x_rep_i(size_t i);

    /** \brief Return the rep of bin of index \c j

	Note that this function returns a value and not a reference.
	This is because we can't return a reference to the internally
	computed representatives, since they don't always exist.
     */
    double get_y_rep_i(size_t j);
    //@}

    /* \brief Set up a twod_intp object for interpolation
       
       \future This is commented out for now. The \ref twod_intp
       object stores a reference to xrep and yrep, and thus 
       can't be used since xrep and yrep don't always exist.
    */ 
    //void setup_interp(twod_intp &ti, bool x_first=true) {
    //ti.set_data(hsize_x,hsize_y,xrep,yrep,wgt,x_first);
    //return;
    //}

    /// Internal consistency check
    void is_valid() const;
    
    /** \brief Create a table3d object based on the histogram data
     */
    void copy_to_table(table3d &t, std::string xreps_name, 
		       std::string yreps_name, std::string weights);
    
    friend void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf, 
				      o2scl::hist_2d &h, std::string name);
    friend void o2scl_hdf::hdf_input(o2scl_hdf::hdf_file &hf, 
				     o2scl::hist_2d &h, std::string name);

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
