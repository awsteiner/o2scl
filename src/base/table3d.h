/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#ifndef O2SCL_TABLE3D_H
#define O2SCL_TABLE3D_H

/** \file table3d.h
    \brief File defining \ref o2scl::table3d
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>
#include <regex>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/misc.h>
#include <o2scl/err_hnd.h>
#include <o2scl/search_vec.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/interp.h>
#include <o2scl/table_units.h>
#include <o2scl/contour.h>
#include <o2scl/vector.h>
#include <o2scl/hist.h>

// Forward definition of the table3d class for HDF I/O
namespace o2scl {
  class table3d;
  // Define for to_hist_2d() function below
  class hist_2d;
  
}

// Forward definition of HDF I/O to extend friendship
namespace o2scl_hdf { 
  class hdf_file; 
  void hdf_input(hdf_file &hf, o2scl::table3d &t, std::string name);
  void hdf_output(hdf_file &hf, const o2scl::table3d &t, std::string name);
}

namespace o2scl {
  
  /** \brief A data structure containing one or more slices of
      two-dimensional data points defined on a grid

      \verbatim embed:rst

      .. todo:: 

         In class table3d:

         - Future: Improve interpolation and derivative caching, possibly
           through non-const versions of the interpolation functions.
         - Future: Should there be a clear_grid() function separate from
           clear_data() and clear()?
         - Future: Allow the user to more clearly probe 'size_set' vs.
           'xy_set'? (AWS 07/18: This is apparently resolved.)

      \endverbatim
  */
  class table3d {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    // This is used for the interpolation classes
    typedef boost::numeric::ublas::matrix_row<const ubmatrix> ubmatrix_row;
    typedef boost::numeric::ublas::matrix_column<const ubmatrix>
      ubmatrix_column;
  
    /** \brief Create a new 3D \table
     */
    table3d();

    virtual ~table3d();

    /** \brief Create a table3d object from a table, assuming \c scolx
	and \c scoly store the x- and y-grid data, respectively.
    */
    table3d(o2scl::table_units<> &t, std::string colx, std::string coly);
    
    /// Copy constructor
    table3d(const table3d &t);
    
    /// Copy constructor
    table3d &operator=(const table3d &t);

    /// \name Initialization
    //@{
    /** \brief Initialize the x-y grid
	
	This function will not allow you to redefine the grid when
	there is data in the \table if a grid of a different size was
	already set from a previous call to either set_xy() or
	set_size(). However, you may freely redefine the grid after a
	call to clear_data() or clear_table(). You may change
	individual grid points at any time with set_grid_x() and
	set_grid_y().
    */
    template<class vec_t, class vec2_t> 
      void set_xy(std::string x_name, size_t nx, const vec_t &x, 
		  std::string y_name, size_t ny, const vec2_t &y) {
      
      if (has_slice && (size_set || xy_set) && (nx!=numx || ny!=numy)) {
	O2SCL_ERR("Size cannot be reset in table3d::set_xy().",
		  o2scl::exc_einval);
	return;
      }
      if (xy_set) {
	xval.clear();
	yval.clear();
      }
      numx=nx;
      numy=ny;
      xname=x_name;
      yname=y_name;
      xval.resize(nx);
      yval.resize(ny);
      for(size_t i=0;i<nx;i++) (xval)[i]=x[i];
      for(size_t i=0;i<ny;i++) (yval)[i]=y[i];
      size_set=true;
      xy_set=true;
      return;
    }

    /** \brief Initialize the x-y grid with \ref uniform_grid objects
	
	This function will not allow you to redefine the grid when
	there is data in the \table if a grid of a different size was
	already set from a previous call to either set_xy() or
	set_size(). However, you may freely redefine the grid after a
	call to clear_data() or clear_table(). You may change
	individual grid points at any time with set_grid_x() and
	set_grid_y().
    */
    void set_xy(std::string x_name, uniform_grid<double> gx, 
		std::string y_name, uniform_grid<double> gy);

    /** \brief Initialize \table size

	This function will not allow you to resize the \table if it
	already has data or if the size has already been set with the
	set_xy() function, unless you clear the data with clear_data()
	or the \table with clear_table() first.
    */
    void set_size(size_t nx, size_t ny);
    //@}

    // --------------------------------------------------------
    /// \name On-grid get and set methods 
    //@{

    /** \brief Set element in slice \c name at location <tt>ix,iy</tt> 
	to value \c val
    */
    void set(size_t ix, size_t iy, std::string name, double val);

    /** \brief Set element in slice of index \c z at location
	<tt>ix,iy</tt> to value \c val .
    */
    void set(size_t ix, size_t iy, size_t z, double val);

    /** \brief Get element in slice \c name at location
	<tt>ix,iy</tt>
    */
    double &get(size_t ix, size_t iy, std::string name);

    /** \brief Get element in slice \c name at location
	<tt>ix,iy</tt> (const version)
    */
    const double &get(size_t ix, size_t iy, std::string name) const;

    /** \brief Get element in slice of index \c z at location
	<tt>ix,iy</tt>
    */
    double &get(size_t ix, size_t iy, size_t z);

    /** \brief Get element in slice of index \c z at location
	<tt>ix,iy</tt> (const version)
    */
    const double &get(size_t ix, size_t iy, size_t z) const;
    //@}

    // --------------------------------------------------------
    /** \name Off-grid get and set methods 

	These methods return the value of a slice on the grid
	point nearest to a user-specified location. For
	interpolation into a point off the grid, use 
	\ref table3d::interp().
    */
    //@{

    /** \brief Set element in slice \c name at the nearest location to
	<tt>x,y</tt> to value \c val
    */
    void set_val(double x, double y, std::string name, double val);

    /** \brief Set element in slice of index \c z at the nearest
	location to <tt>x,y</tt> to value \c val
    */
    void set_val(double x, double y, size_t z, double val);
    
    /** \brief Get element in slice \c name at location closest to
	<tt>x,y</tt>
    */
    double &get_val(double x, double y, std::string name);
    
    /** \brief Get element in slice \c name at location closest to
	<tt>x,y</tt>
    */
    const double &get_val(double x, double y, std::string name) const;
    
    /** \brief Get element in slice of index \c z at location closest
	to <tt>x,y</tt>
    */
    double &get_val(double x, double y, size_t z);

    /** \brief Get element in slice of index \c z at location closest
	to <tt>x,y</tt>
    */
    const double &get_val(double x, double y, size_t z) const;

    /** \brief Set elements in the first <tt>nv</tt> slices at the
	nearest location to <tt>x,y</tt> to value \c val
    */
    template<class vec_t> 
      void set_slices(double x, double y, size_t nv, vec_t &vals) {
      size_t ix, iy;
      lookup_x(x,ix);
      lookup_y(y,iy);

      for(size_t i=0;i<nv && i<list.size();i++) {
	list[i](ix,iy)=vals[i];
      }
      return;
    }
    
    /** \brief Get the data for every slice at the nearest location to
	<tt>x,y</tt>
    */
    template<class vec_t>
      void get_slices(double x, double y, size_t nv, vec_t &v) {

      size_t ix, iy;

      lookup_x(x,ix);
      lookup_y(y,iy);

      for(size_t i=0;i<nv && i<list.size();i++) {
	v[i]=list[i](ix,iy);
      }
      
      return;
    }
    //@}

    // --------------------------------------------------------
    /// \name Off-grid get and set methods returning nearest point
    //@{
    
    /** \brief Set element in slice \c name at the nearest location to
	<tt>x,y</tt> to value \c val
    */
    void set_val_ret(double &x, double &y, std::string name, double val);

    /** \brief Set element in slice of index \c z at the nearest
	location to <tt>x,y</tt> to value \c val
    */
    void set_val_ret(double &x, double &y, size_t z, double val);
    
    /** \brief Get element in slice \c name at location closest to
	<tt>x,y</tt>, and also return the corresponding values of \c x
	and \c y
    */
    double &get_val_ret(double &x, double &y, std::string name);
    
    /** \brief Get element in slice \c name at location closest to
	<tt>x,y</tt>, and also return the corresponding values of \c x
	and \c y
    */
    const double &get_val_ret(double &x, double &y, std::string name) const;
    
    /** \brief Get element in slice of index \c z at location closest
	to <tt>x,y</tt>, and also return the corresponding values of
	\c x and \c y
    */
    double &get_val_ret(double &x, double &y, size_t z);

    /** \brief Get element in slice of index \c z at location closest
	to <tt>x,y</tt>, and also return the corresponding values of
	\c x and \c y
    */
    const double &get_val_ret(double &x, double &y, size_t z) const;

    /** \brief This function adds a slice from a different table3d
	object, interpolating the results into the current 
	table3d object
    */
    void add_slice_from_table(table3d &source, std::string slice,
			      std::string dest_slice="",
                              int verbose=0);
    
    /** \brief Set elements in the first <tt>nv</tt> slices at the
	nearest location to <tt>x,y</tt> to values \c vals
    */
    template<class vec_t> 
      void set_slices_ret(double &x, double &y, size_t nv, vec_t &vals) {
      size_t ix, iy;
      lookup_x(x,ix);
      lookup_y(y,iy);
      x=xval[ix];
      y=yval[iy];

      for(size_t i=0;i<nv && i<list.size();i++) {
	list[i](ix,iy)=vals[i];
      }
      return;
    }

    /** \brief Get elements in the first <tt>nv</tt> slices at the
	nearest location to <tt>x,y</tt> to value \c val
    */
    template<class vec_t> 
      void get_slices_ret(double &x, double &y, size_t nv, vec_t &vals) {

      size_t ix, iy;
      lookup_x(x,ix);
      lookup_y(y,iy);
      x=xval[ix];
      y=yval[iy];
      
      for(size_t i=0;i<nv && i<list.size();i++) {
	vals[i]=list[i](ix,iy);
      }
      return;
    }
    
    //@}

    // --------------------------------------------------------
    /// \name Grid information get and set methods
    //@{

    /// Set x grid point at index \c ix
    void set_grid_x(size_t ix, double val);
    
    /// Set y grid point at index \c iy
    void set_grid_y(size_t iy, double val);
    
    /// Get x grid point at index \c ix
    double get_grid_x(size_t ix) const;

    /// Get y grid point at index \c iy
    double get_grid_y(size_t iy) const;

    /// Get the name of the x grid variable
    std::string get_x_name() const;
    
    /// Get the name of the y grid variable
    std::string get_y_name() const;

    /// Set the name of the x grid variable
    void set_x_name(std::string name);
    
    /// Set the name of the y grid variable
    void set_y_name(std::string name);

    /// Get a const reference to the full x grid
    const ubvector &get_x_data() const;

    /// Get a const reference to the full y grid
    const ubvector &get_y_data() const;
    //@}

    // --------------------------------------------------------
    /// \name Size get methods
    //@{
    /// Get the size of the slices
    void get_size(size_t &nx, size_t &ny) const;

    /// Get the x size
    size_t get_nx() const;
    
    /// Get the y size
    size_t get_ny() const;

    /// Get the number of slices
    size_t get_nslices() const;

    /// True if the size of the table has been set
    bool is_size_set() const;

    /// True if the grid has been set
    bool is_xy_set() const;
    //@}

    // --------------------------------------------------------
    /// \name Slice manipulation 
    //@{
    
    /// Create a set of new slices specified in the string \c names
    void line_of_names(std::string names);

    /** \brief Returns the name of slice with index \c z
     */
    std::string get_slice_name(size_t z) const;

    /** \brief Add a new slice
     */
    void new_slice(std::string name);

    /** \brief Set all of the values in slice \c name to \c val
     */
    void set_slice_all(std::string name, double val);

    /** \brief Find the index for slice named \c name
     */
    size_t lookup_slice(std::string name) const;
    
    /** \brief Return true if slice is already present
     */
    bool is_slice(std::string name, size_t &ix) const;
    
    /** \brief Rename slice named \c olds to \c news
	
	This is slow since we have to delete the column and re-insert
	it. This process in turn mangles all of the iterators in the
	list.
    */
    void rename_slice(std::string olds, std::string news);
    
    /** \brief Make a new slice named \c dest which is a copy 
	of the slice with name given in \c src. 
    */
    void copy_slice(std::string src, std::string dest);
    
    /** \brief Initialize all values of slice named \c scol to \c val 

	\note This will call the error handler if the value \c val is
	not finite (i.e. either <tt>Inf</tt> or <tt>NaN</tt>).
    */
    void init_slice(std::string slice, double val);

    /// Return a constant reference to a slice
    const ubmatrix &get_slice(std::string slice) const;

    /// Return a constant reference to a slice
    const ubmatrix &get_slice(size_t iz) const;

    /// Return a constant reference to a slice
    ubmatrix &get_slice(std::string slice);
    
    /// Return a constant reference to a slice
    ubmatrix &get_slice(size_t iz);

    /** \brief Return a constant reference to all the slice data

	\comment
	This isn't designated const, i.e. as
	const std::vector<ubmatrix> &get_data() const;
	because it would then have to be 
	const std::vector<const ubmatrix> &get_data() const;
	\endcomment
    */
    const std::vector<ubmatrix> &get_data();

    /** \brief Copy to a slice from a generic matrix object
	
	The type <tt>mat_t</tt> can be any type with an
	<tt>operator(,)</tt> method.
    */
    template<class mat_t> 
      void copy_to_slice(mat_t &m, std::string slice_name) {
      for(size_t i=0;i<numx;i++) {
	for(size_t j=0;j<numy;j++) {
	  this->set(i,j,slice_name,m(i,j));
	}
      }
      return;
    }
    //@}
  
    // --------------------------------------------------------
    /// \name Lookup and search methods
    //@{
    /** \brief Look for a value in the x grid
     */
    void lookup_x(double val, size_t &ix) const;
  
    /** \brief Look for a value in the y grid
     */
    void lookup_y(double val, size_t &iy) const;

    /** \brief Look for a value in a specified slice
     */
    void lookup(double val, std::string slice, size_t &ix,
		size_t &iy) const;
    //@}

    // --------------------------------------------------------
    /// \name Interpolation, differentiation, and integration
    //@{

    /** \brief Specify the interpolation type
     */
    void set_interp_type(size_t interp_type);

    /** \brief Get the interpolation type
     */
    size_t get_interp_type() const;
    
    /** \brief Interpolate \c x and \c y in slice named \c name
     */
    double interp(double x, double y, std::string name) const;

    /** \brief Interpolate the derivative of the data with respect to
	the x grid at point \c x and \c y in slice named \c name
    */
    double deriv_x(double x, double y, std::string name) const;

    /** \brief Interpolate the derivative of the data with respect to
	the y grid at point \c x and \c y in slice named \c name
    */
    double deriv_y(double x, double y, std::string name) const;

    /** \brief Interpolate the mixed second derivative of the data at
	point \c x and \c y in slice named \c name
    */
    double deriv_xy(double x, double y, std::string name) const;

    /** \brief Interpolate the integral of the data 
	respect to the x grid 
    */
    double integ_x(double x1, double x2, double y, std::string name) const;

    /** \brief Interpolate the integral of the data 
	respect to the y grid 
    */
    double integ_y(double x, double y1, double y2, std::string name) const;

    /** \brief Fill a vector of interpolated values from each slice at the
	point <tt>x,y</tt>
    */
    template<class vec_t>
      void interp_slices(double x, double y, size_t nv, vec_t &v) {
      
      for (size_t i=0;i<list.size();i++) {
	std::string name=get_slice_name(i);
	v[i]=interp(x,y,name);
      }

      return;
    }

    /** \brief Create a new slice, named \c fpname, containing the 
	derivative of \c fname with respect to the x coordinate
    */
    void deriv_x(std::string fname, std::string fpname);

    /** \brief Create a new slice, named \c fpname, containing the 
	derivative of \c fname with respect to the y coordinate
    */
    void deriv_y(std::string fname, std::string fpname);
    //@}

    // --------------------------------------------------------
    /// \name Extract 2-dimensional tables
    //@{
    /** \brief Extract a table at a fixed x grid point 
	
	\note All of the information previously stored in \c t will
	be lost.
    */
    void extract_x(double x, table<> &t);
    
    /** \brief Extract a table at a fixed y grid point 

	\note All of the information previously stored in \c t will
	be lost.
    */
    void extract_y(double y, table<> &t);
    //@}

    // --------------------------------------------------------
    /// \name Clear methods 
    //@{
    /** \brief Zero the data entries but keep the slice names
	and grid
    */
    void zero_table();

    /** \brief Clear everything
     */
    void clear();

    /** \brief Remove all of the data by setting the number
	of lines to zero

	This leaves the column names intact and does not remove
	the constants.
    */
    void clear_data();
    //@}

    // --------------------------------------------------------
    /// \name Summary method 
    //@{
    /** \brief Output a summary of the information stored
	
	Outputs the number of constants, the grid information,
	and a list of the slice names
    */
    void summary(std::ostream *out=0, int ncol=79) const;
    //@}

    // ---------
    // Allow HDF I/O functions to access table3d data

    friend void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf,
				      const table3d &t, 
				      std::string name);
    
    friend void o2scl_hdf::hdf_input(o2scl_hdf::hdf_file &hf, table3d &t, 
				     std::string name);
    
    // --------------------------------------------------------
    /// \name Contour lines method
    //@{

    /** \brief Create contour lines from the slice named \c name

	This uses \ref contour to compute contour lines (stored in \c
	clines) from slice \c name given \c nlev contour levels in \c
	levs .
    */
    template<class vec_t> 
    void slice_contours(std::string name, size_t nlev, vec_t &levs,
                        std::vector<contour_line> &clines) {
      
      size_t z=lookup_slice(name);
      
      contour co;
      co.set_data(numx,numy,xval,yval,list[z]);
      co.set_levels(nlev,levs);
      co.calc_contours(clines);

      return;
    }
    //@}

    // --------------------------------------------------------
    /// \name Manipulating constants
    //@{
    /** \brief Add a constant, or if the constant already exists, change 
	its value
    */
    virtual void add_constant(std::string name, double val);

    /// Remove a constant
    virtual void remove_constant(std::string name);
    
    /** \brief Set a constant equal to a value, but don't add it if
	not already present
	
	If \c err_on_notfound is <tt>true</tt> (the default), then
	this function throws an exception if a constant with
	name \c name is not found. If \c err_on_notfound is
	<tt>false</tt>, then if a constant with name \c name
	is not found this function just silently returns
	\ref o2scl::exc_enotfound.
    */
    virtual int set_constant(std::string name, double val,
			     bool err_on_notfound=true);

    /// Test if \c name is a constant
    virtual bool is_constant(std::string name) const;

    /// Get a constant
    virtual double get_constant(std::string name);

    /// Get a constant by index
    virtual void get_constant(size_t ix, std::string &name,
			      double &val) const;

    /// Get the number of constants
    virtual size_t get_nconsts() const {
      return constants.size();
    }
    //@}

    /// \name Miscellaneous methods
    //@{
    /** \brief Read a generic table3d object specified as a 
	text file

	This function reads a set of columns of numerical values,
	presuming that the first column is the x-grid value, the
	second column is the y-grid value, and the remaining columns
	are slices to be added. If the first row appears to be strings
	rather than numerical values, then the first row is used for
	the x name, y name, and slice names. Values in the first two
	columns which differ by less than \c eps are assumed to refer
	to the same grid point. If not all combinations of x and y are
	found, then those entries are left unchanged in all slices.

        \verbatim embed:rst
        
        .. todo:: 
        
           In class table3d::read_gen3_list():
           
           Future: It would be great to add a function which generates
           a text file in this format as well. 

           Future: This function is very similar to from_table() below,
           so it might be possible just to avoid code duplication
           between these two functions.

        \endverbatim
    */
    int read_gen3_list(std::istream &fin, int verbose=0,
		       double eps=1.0e-12);

    /** \brief Set the current table3d object by reading a 
	\ref o2scl::table 

        The function reads the table \c tab and attempts to convert it
        to a table3d object by using column \c xname2 and \c yname2 as
        the values for the x- and y-grids. This function is particularly
        useful for a table which has the following structure:
        \verbatim
        x y z
        1 1 4.5
        1 2 2.0
        1 3 1.6
        2 1 1.5
        2 3 4.7
        3 1 3.9
        3 2 4.5
        3 3 4.8
        \endverbatim
        If \c xname2 or \c yname2 are empty strings, then the first or
        second column is presumed to hold the value for the x- or
        y-grid, respectively. The grids in the table3d object are set
        to include all possible values in the associated column,
        treating any values with a relative deviation of \c eps as
        identical. In the example above, using "x" and "y" as the
        columns for the table3d grid, the grids would both be [1,2,3].

        New slices are created in the table3d for each column in the
        table object which is not part of the grid. Any current data
        stored in the table3d object is destroyed. All slices are
        filled with \c empty_value before being assigned values from
        \c tab. For example, in the table above, the slice named "z"
        in the new table3d object would have a final entry of 0.0 for
        (x,y) at (2,2) because there is no entry in the original table
        for that location.
     */
    template<class vec_t>
      int read_table(const o2scl::table<vec_t> &tab, 
		     std::string xname2="", std::string yname2="",
		     double empty_value=0.0, int verbose=0,
		     bool err_on_fail=true, double eps=1.0e-12) {

      clear();
      
      if (tab.get_ncolumns()<3) {
	if (err_on_fail) {
	  O2SCL_ERR2("Not enough columns of data in ",
		     "table3d::read_table().",o2scl::exc_einval);
	} else {
	  return o2scl::exc_einval;
	}
      }

      // Set up xname and yname if unspecified
      if (xname2.length()==0) {
	if (yname2==tab.get_column_name(0)) {
	  xname2=tab.get_column_name(1);
	} else {
	  xname2=tab.get_column_name(0);
	  yname2=tab.get_column_name(1);
	}
      }
      if (yname2.length()==0) {
	if (xname2==tab.get_column_name(1)) {
	  yname2=tab.get_column_name(0);
	} else {
	  yname2=tab.get_column_name(1);
	}
      }
      
      // Setup x and y grid vectors from data
      const vec_t &xdata=tab[xname2];
      const vec_t &ydata=tab[yname2];
      
      std::vector<double> xgrid, ygrid;
      
      // Note it is important that this loop ends at tab.get_nlines()
      // rather than xdata.size() because the internal vector storage
      // can be larger than the actual table size
      for(size_t i=0;i<tab.get_nlines();i++) {

	// Look for x value in x grid
	bool found=false;
	for(size_t j=0;j<xgrid.size();j++) {
	  if (fabs(xdata[i]-xgrid[j])/fabs(xgrid[j])<eps) {
	    found=true;
	  }
	}

	// If not found, add it
	if (found==false) {
	  xgrid.push_back(xdata[i]);
	}

	// Now look for y value in y grid
	found=false;
	for(size_t j=0;j<ygrid.size();j++) {
	  if (fabs(ydata[i]-ygrid[j])/fabs(ygrid[j])<eps) {
	    found=true;
	  }
	}

	// If not found, add it
	if (found==false) {
	  ygrid.push_back(ydata[i]);
	}
      }
      
      if (verbose>1) {
	std::cout << "X grid: " << std::endl;
	for(size_t k=0;k<xgrid.size();k++) {
	  std::cout << k << " " << xgrid[k] << std::endl;
	}
	std::cout << "Y grid: " << std::endl;
	for(size_t k=0;k<ygrid.size();k++) {
	  std::cout << k << " " << ygrid[k] << std::endl;
	}
      }

      // Sor the grids
      vector_sort_double(xgrid.size(),xgrid);
      vector_sort_double(ygrid.size(),ygrid);
      
      // Set grid from x and y grid vectors
      set_xy(xname2,xgrid.size(),xgrid,yname2,ygrid.size(),ygrid);
      
      // Create new slices
      std::vector<std::string> sl_names;
      for(size_t i=0;i<tab.get_ncolumns();i++) {
	if (tab.get_column_name(i)!=xname2 &&
            tab.get_column_name(i)!=yname2) {
	  std::string sl=tab.get_column_name(i);
	  if (verbose>0) {
	    std::cout << "New slice: " << sl << std::endl;
	  }
	  sl_names.push_back(sl);
	  new_slice(sl);
	  set_slice_all(sl,empty_value);
	}
      }
      
      // Set the data
      for(size_t i=0;i<tab.get_ncolumns();i++) {
	if (tab.get_column_name(i)!=xname2 &&
            tab.get_column_name(i)!=yname2) {
	  std::string sl=tab.get_column_name(i);
	  for(size_t j=0;j<tab.get_nlines();j++) {
	    set_val(xdata[j],ydata[j],sl,tab.get(i,j));
	  }
	}
      }
      
      return 0;
    }
    
    /** \brief Create a table3d object by histogramming a series of
        columns from a \ref o2scl::table_units object

        Create a new table3d object from a histogram a series of
        columns from a \ref table_units object. If \c direction is
        "x", then arrange these histograms "vertically", so that the
        x-coordinate (named \c name) of the ith column is taken from
        the ith entry of \c grid. If \c direction is "y", then arrange
        these histograms "horizontally", so that the y-coordinate
        (named \c name) of the ith column is taken from the ith entry
        of \c grid. The histograms (in either case) are all created
        using the bin edges from \c bin_edges. When \c direction is
        "x" ("y"), \c bin_grid is used for the y-coordinate
        (x-coordinate) of the new table3d object. This coordinate is
        named \c bin_name. The columns are taken from all those
        columns in \c t which match the regular expression in 
        \c pattern. All of the new histogram data is imported
        into a slice named \c slide in a new \ref table3d object. one

        The vector \c bin_grid must have a size which is exactly 1
        smaller than the size of the vector \c bin_edges. The number
        of columns matched from table \c t by the pattern specified in
        \c pattern must be exactly equal to the size of the vector \c
        grid.

        Any data in the current table3d object is destroyed.
    */
    template<class vec_t, class vec2_t, class vec3_t>
    void create_table_hist_set(vec_t &grid, std::string direction,
                               std::string name, vec2_t &bin_edges,
                               vec3_t &bin_grid, std::string bin_name,
                               o2scl::table_units<> &t, std::string pattern,
                               std::string slice, bool use_regex=false,
                               int verbose=0) {

      clear();

      if (bin_grid.size()+1!=bin_edges.size()) {
        O2SCL_ERR2("bin_grid and bin_edges vectors not properly sized ",
                   "in table3d::create_table_hist_set().",
                   o2scl::exc_einval);
      }
      
      // First, get all the columns which match the pattern
      std::vector<std::string> matched;

      if (use_regex) {
        for(size_t j=0;j<t.get_ncolumns();j++) {
          std::regex r(pattern);
          if (std::regex_search(t.get_column_name(j),r)) {
            matched.push_back(t.get_column_name(j));
          }
        }
      } else {
        for(size_t j=0;j<t.get_ncolumns();j++) {
          if (fnmatch(pattern.c_str(),
                      t.get_column_name(j).c_str(),0)==0) {
            matched.push_back(t.get_column_name(j));
          }
        }
      }

      if (verbose>0) {
        std::cout << "In table3d::create_table_hist_set(), "
                  << "columns matched: " << std::endl;
        o2scl::vector_out(std::cout,matched,true);
      }

      if (grid.size()!=matched.size()) {
        std::string errs=o2scl::szttos(matched.size())+
          " columns matched, but the grid had size "+
          o2scl::szttos(grid.size())+
          " in table3d::create_table_hist_set().";
        O2SCL_ERR(errs.c_str(),o2scl::exc_efailed);
      }

      if (direction=="x") {

        numx=grid.size();
        xname=name;
        xval.resize(numx);
        o2scl::vector_copy(numx,grid,xval);

        if (verbose>0) {
          std::cout << "In table3d::create_table_hist_set(), x grid: ";
          std::cout << numx << " " << xname << " ";
          o2scl::vector_out(std::cout,xval,true);
        }
        
        numy=bin_grid.size();
        yname=bin_name;
        yval.resize(numy);
        o2scl::vector_copy(numy,bin_grid,yval);

        if (verbose>0) {
          std::cout << "In table3d::create_table_hist_set(), y grid: ";
          std::cout << numy << " " << yname << " ";
          o2scl::vector_out(std::cout,yval,true);
        }
        
        xy_set=true;
        size_set=true;
        has_slice=false;

        new_slice(slice);

        double hist_min=bin_edges[0];
        double hist_max=bin_edges[bin_edges.size()-1];
        
        // Create the data
        for(size_t i=0;i<numx;i++) {
          // Create the histogram for this x-coordinate
          hist h;
          h.set_bin_edges(bin_edges.size(),bin_edges);
          for(size_t j=0;j<t.get_nlines();j++) {
            double val=t[matched[i]][j];
            if (val>=hist_min && val<=hist_max) {
              h.update(val);
            }
          }
          // Now copy the histogram to the table3d object
          for(size_t j=0;j<h.size();j++) {
            this->set(i,j,slice,h[j]);
          }
        }
        
      } else if (direction=="y") {

        numx=bin_grid.size();
        xname=bin_name;
        xval.resize(numx);
        o2scl::vector_copy(numx,bin_grid,xval);

        if (verbose>0) {
          std::cout << "In table3d::create_table_hist_set(), x grid: ";
          std::cout << numx << " " << xname << " ";
          o2scl::vector_out(std::cout,xval,true);
        }

        numy=grid.size();
        yname=name;
        yval.resize(numy);
        o2scl::vector_copy(numy,grid,yval);
        
        if (verbose>0) {
          std::cout << "In table3d::create_table_hist_set(), y grid: ";
          std::cout << numy << " " << yname << " ";
          o2scl::vector_out(std::cout,yval,true);
        }

        xy_set=true;
        size_set=true;
        has_slice=false;

        new_slice(slice);

        double hist_min=bin_edges[0];
        double hist_max=bin_edges[bin_edges.size()-1];

        // Create the data
        for(size_t i=0;i<numy;i++) {
          // Create the histogram for this x-coordinate
          hist h;
          h.set_bin_edges(bin_edges.size(),bin_edges);
          for(size_t j=0;j<t.get_nlines();j++) {
            double val=t[matched[i]][j];
            if (val>=hist_min && val<=hist_max) {
              h.update(val);
            }
          }
          // Now copy the histogram to the table3d object
          for(size_t j=0;j<h.size();j++) {
            this->set(j,i,slice,h[j]);
          }
        }
        
      } else {
        O2SCL_ERR2("Direction must be \"x\" or \"y\" in ",
                   "table3d::table3d().",o2scl::exc_efailed);
      }

      return;
    }
    
    /** \brief Create a table3d object by histogramming a series of
        columns from a \ref o2scl::table_units object

        This function works very similarly to the more detailed
        function with the same name, but it uses the minimum and
        maximum values of the table columns in order to automatically
        create the histogram bin edges from a set of \c n_bins bins.
        It uses these bin edges to create the \c bin_grid and \c
        bin_edges objects.
    */
    template<class vec_t>
    void create_table_hist_set_minmax
    (vec_t &grid, std::string direction, std::string name, size_t n_bins,
     std::string bin_name, o2scl::table_units<> &t, std::string pattern,
     std::string slice, double factor=1.0e-10,
     bool use_regex=false, int verbose=0) {
      
      // First, get all the columns which match the pattern
      std::vector<std::string> matched;

      if (use_regex) {
        for(size_t j=0;j<t.get_ncolumns();j++) {
          std::regex r(pattern);
          if (std::regex_search(t.get_column_name(j),r)) {
            matched.push_back(t.get_column_name(j));
          }
        }
      } else {
        for(size_t j=0;j<t.get_ncolumns();j++) {
          if (fnmatch(pattern.c_str(),
                      t.get_column_name(j).c_str(),0)==0) {
            matched.push_back(t.get_column_name(j));
          }
        }
      }

      if (matched.size()==0) {
        O2SCL_ERR("No columns matched in create_table_hist_set().",
                  o2scl::exc_einval);
      }
      
      double min=t.get(matched[0],0);
      double max=min;
      for(size_t j=0;j<matched.size();j++) {
        for(size_t i=0;i<t.get_nlines();i++) {
          double val=t.get(matched[j],i);
          if (val<min) min=val;
          else if (val>max) max=val;
        }
      }

      double delta=(max-min)/n_bins;
      
      if (verbose>0) {
        std::cout << "In table3d::create_table_hist_set(), matched: ";
        o2scl::vector_out(std::cout,matched,true);
        std::cout << "  min,max,n_bins,delta: "
                  << min << " " << max << " " << n_bins << " "
                  << delta << std::endl;
      }
        
      uniform_grid_end<double> ug(min-delta*factor,max+delta*factor,n_bins);
      std::vector<double> bin_edges(n_bins+1), bin_grid(n_bins);
      ug.vector(bin_edges);
      for(size_t i=0;i<n_bins;i++) {
        bin_grid[i]=(bin_edges[i]+bin_edges[i+1])/2.0;
      }
      if (verbose>0) {
        std::cout << "In table3d::create_table_hist_set(): " << std::endl;
        std::cout << "  bin_grid: ";
        o2scl::vector_out(std::cout,bin_grid,true);
        std::cout << "  bin_edges: ";
        o2scl::vector_out(std::cout,bin_edges,true);
      }
      
      create_table_hist_set(grid,direction,name,bin_edges,bin_grid,
                            bin_name,t,pattern,slice,use_regex);
      
      return;
    }

    /** \brief Create a table3d object by histogramming a series of
        columns from a \ref o2scl::table_units object

        This function works very similarly to the more detailed
        function with the same name, but it uses \c bin_edges to
        automatically compute the \c bin_grid argument. If \c
        bin_edges appears logarithmic, it uses the geometric mean of
        adjacent edges, otherwise it uses the arithmetic mean.
     */
    template<class vec_t, class vec2_t>
    void create_table_hist_set_edgeonly
    (vec_t &grid, std::string direction, std::string name, vec2_t &bin_edges,
     std::string bin_name, o2scl::table_units<> &t, std::string pattern,
     std::string slice, bool use_regex=false, int verbose=0) {
      
      bool log=false;
      linear_or_log(bin_edges,log);
      std::vector<double> bin_grid(bin_edges.size()-1);
      if (log) {
        for(size_t i=0;i<bin_edges.size()-1;i++) {
          bin_grid[i]=sqrt(bin_edges[i]*bin_edges[i+1]);
        }
      } else {
        for(size_t i=0;i<bin_edges.size()-1;i++) {
          bin_grid[i]=(bin_edges[i]+bin_edges[i+1])/2.0;
        }
      }

      if (verbose>0) {
        std::cout << "In table3d::create_table_hist_set(), log: "
                  << log << std::endl;
        std::cout << "  bin_grid: ";
        o2scl::vector_out(std::cout,bin_grid,true);
        std::cout << "  bin_edges: ";
        o2scl::vector_out(std::cout,bin_edges,true);
      }

      create_table_hist_set(grid,direction,name,bin_edges,bin_grid,
                            bin_name,t,pattern,slice,use_regex,verbose);
      
      return;
    }
    
    /// Return the type, \c "table3d".
    virtual const char *type() { return "table3d"; }
    //@}

    /** \name Parsing mathematical functions specified as strings
	
	\comment
	Note that doxygen doesn't format the documentation propertly
	if the \name specification covers more than one line
	\endcomment
    */
    //@{
    /** \brief Fill a matrix from the function specified in \c function

	\comment
	This function must return an int rather than void because
	of the presence of the 'throw_on_err' mechanism
	\endcomment
    */
    template<class resize_mat_t>
      int function_matrix(std::string function, resize_mat_t &mat,
			  bool throw_on_err=true) {
      
      calc_utf8<> calc;
      std::map<std::string,double> vars;

      std::map<std::string,double>::const_iterator mit;
      for(mit=constants.begin();mit!=constants.end();mit++) {
	vars[mit->first]=mit->second;
      }

      calc.compile(function.c_str(),&vars);

      if (mat.size1()!=numx || mat.size2()!=numy) {
	mat.resize(numx,numy);
      }

      for(size_t i=0;i<numx;i++) {
	for(size_t j=0;j<numy;j++) {
	  vars[xname]=xval[i];
	  vars[yname]=yval[j];
	  
	  for(size_t k=0;k<list.size();k++) {
	    vars[get_slice_name(k)]=list[k](i,j);
	  }
	  
	  mat(i,j)=calc.eval(&vars);
	}
      }
      
      
      return 0;
    }

    /** \brief Make a column from <tt>function</tt> and add it to the table.
	
	If the column already exists, the data already present is 
	overwritten with the result.
    */
    void function_slice(std::string function, std::string col);
    //@}

    /** \brief Copy slice named \c slice to a new \ref o2scl::table3d
	object with a uniform grid using the current interpolation type
    */
    table3d slice_to_uniform_grid(std::string slice, size_t xpts,
				  bool log_x, size_t ypts, bool log_y);
    
    /** \brief Copy entire table to a new \ref o2scl::table3d
	object with a uniform grid using the current interpolation type
    */
    table3d table_to_uniform_grid(size_t xpts, bool log_x, 
				  size_t ypts, bool log_y);

    /** \brief Convert slice named \c slice to a \ref hist_2d object
     */
    hist_2d to_hist_2d(std::string slice, int verbose=1);

#ifdef O2SCL_NEVER_DEFINED
    
    /** \brief Desc
     */
    template<class vec_t> 
      void test_uniform_grid_log(std::string slice, bool &log_x, bool &log_y,
				 vec_t &x, vec_t &y, vec_t &s) {
      vector<double> dev_x;
      for(size_t i=0;i<numx-1;i++) {
	dev_x.push_back(xval[i+1]-xval[i]);
      }
      double avg=vector_mean(dev_x.size(),dev_x);
      double stddev=vector_stddev(dev_x.size(),dev_x);
      bool x_set=false;
      log_x=false;
      if (stddev<1.0e-8*fabs(avg)) {
	x.resize(numx);
	for(size_t i=0;i<numx;i++) {
	  x[i]=xval[i];
	}
	x_set=true;
      } else {
	dev_x.clear();
	for(size_t i=0;i<numx-1;i++) {
	  dev_x.push_back(xval[i+1]/xval[i]);
	}
	avg=vector_mean(dev_x.size(),dev_x);
	stddev=vector_stddev(dev_x.size(),dev_x);
	if (stddev<1.0e-8*fabs(avg)) {
	  x.resize(numx);
	  for(size_t i=0;i<numx;i++) {
	    x[i]=xval[i];
	  }
	  x_set=true;
	  log_x=true;
	}
      }
      if (x_set==false) {
	uniform_grid_end
	  }
      return;
    }
#endif
    
  protected:

#ifndef DOXYGEN_INTERNAL
    
    /// \name Iterator types
    //@{
    typedef std::map<std::string,size_t,
      std::greater<std::string> >::iterator map_iter;
    typedef std::map<std::string,size_t,
      std::greater<std::string> >::const_iterator map_const_iter;
    //@}
  
    /// \name Data storage
    //@{
    /// The list of constants
    std::map<std::string,double> constants;

    /// The size of the x grid
    size_t numx;

    /// The size of the y grid
    size_t numy;

    /// A tree connecting column names to list indexes
    std::map<std::string,size_t,std::greater<std::string> > tree;

    /// The name for the x grid
    std::string xname;

    /// The name for the y grid
    std::string yname;

    /// The pointers to the matrices
    std::vector<ubmatrix> list;

    /// The x grid
    ubvector xval;

    /// The y grid
    ubvector yval;

    /// True if the grid has been set
    bool xy_set;

    /// True if the size of the grid has been set
    bool size_set;

    /// True if the table has at least one slice
    bool has_slice;

    /// The interpolation type
    size_t itype;
    //@}
  
    /// \name Tree iterator boundaries
    //@{
    /// Return the beginning of the slice tree
    map_iter begin() {return tree.begin();};
    /// Return the end of the slice tree
    map_iter end() {return tree.end();};
    /// Return the beginning of the slice tree
    map_const_iter const_begin() const {return tree.begin();};
    /// Return the end of the slice tree
    map_const_iter const_end() const {return tree.end();};
    //@}

#endif

  };

}

#endif
