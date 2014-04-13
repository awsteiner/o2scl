/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <o2scl/misc.h>
#include <o2scl/err_hnd.h>
#include <o2scl/multi_funct.h>
#include <o2scl/search_vec.h>
#include <o2scl/uniform_grid.h>
#include <o2scl/interp.h>
#include <o2scl/table_units.h>
#include <o2scl/contour.h>

// Forward definition of the table3d class for HDF I/O
namespace o2scl {
  class table3d;
}

// Forward definition of HDF I/O to extend friendship
namespace o2scl_hdf { 
  class hdf_file; 
  void hdf_input(hdf_file &hf, o2scl::table3d &t, std::string name);
  void hdf_output(hdf_file &hf, o2scl::table3d &t, std::string name);
}

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief A data structure containing many slices of two-dimensional
      data points defined on a grid

      \future Improve interpolation and derivative caching
      \future Make a 'const' version of the interpolation functions
      \future Should there be a clear_grid() function separate from
      clear_data() and clear_table()?
      \future Allow the user to more clearly probe 'size_set' vs.
      'xy_set'?
  */
  class table3d {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::matrix_row<ubmatrix> ubmatrix_row;
    typedef boost::numeric::ublas::matrix_column<ubmatrix> ubmatrix_column;
  
    /** \brief Create a new 3D \table
     */
    table3d();

    virtual ~table3d();

    /** \brief Create a table3d object from a table, assuming \c scolx
	and \c scoly store the x- and y-grid data, respectively.
     */
    table3d(o2scl::table_units<> &t, std::string colx, std::string coly);

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
		std::string y_name, uniform_grid<double> gy) {

      if (has_slice && (size_set || xy_set) && 
	  (gx.get_npoints()!=numx || gy.get_npoints()!=numy)) {
	O2SCL_ERR("Size cannot be reset in table3d::set_xy().",
		  o2scl::exc_einval);
	return;
      }

      if (xy_set) {
	xval.clear();
	yval.clear();
      }
      numx=gx.get_npoints();
      numy=gy.get_npoints();
      xname=x_name;
      yname=y_name;
      xval.resize(numx);
      yval.resize(numy);
      gx.vector(xval);
      gy.vector(yval);
      size_set=true;
      xy_set=true;
    }

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

    /** \brief Desc
     */
    void add_slice_from_table(table3d &source, std::string slice,
			      std::string dest_slice="") {

      if (dest_slice.length()==0) dest_slice=slice;

      if (xy_set==false) {
	set_xy(source.get_x_name(),source.get_nx(),source.get_x_data(),
	       source.get_y_name(),source.get_ny(),source.get_y_data());
	new_slice(dest_slice);
	for(size_t i=0;i<numx;i++) {
	  for(size_t j=0;j<numx;j++) {
	    set(i,j,dest_slice,source.get(i,j,slice));
	  }
	}
	return;
      }

      size_t szt_tmp;
      if (!is_slice(dest_slice,szt_tmp)) new_slice(dest_slice);
      for(size_t i=0;i<numx;i++) {
	for(size_t j=0;j<numx;j++) {
	  set(i,j,dest_slice,source.interp(get_grid_x(i),get_grid_y(j),
					   slice));
	}
      }
      return;
    }

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
    double get_grid_x(size_t ix);

    /// Get y grid point at index \c iy
    double get_grid_y(size_t iy);

    /// Get the name of the x grid variable
    std::string get_x_name() const {
      return xname;
    }

    /// Get the name of the y grid variable
    std::string get_y_name() const {
      return yname;
    }

    /// Set the name of the x grid variable
    void set_x_name(std::string name) {
      xname=name;
      return;
    }
    
    /// Set the name of the y grid variable
    void set_y_name(std::string name) {
      yname=name;
      return;
    }

    /// Get a const reference to the full x grid
    const ubvector &get_x_data() const { return xval; }

    /// Get a const reference to the full y grid
    const ubvector &get_y_data() const { return yval; }

    //@}

    // --------------------------------------------------------
    /// \name Size get methods
    //@{

    /// Get the size of the slices
    void get_size(size_t &nx, size_t &ny) const;

    /// Get the x size
    size_t get_nx() const {
      return numx;
    }
    
    /// Get the y size
    size_t get_ny() const {
      return numy;
    }

    /// Get the number of slices
    size_t get_nslices() const;
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

    /** \brief Find the index for column named \c name
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
    void init_slice(std::string scol, double val);

    /// Return a constant reference to a slice
    const ubmatrix &get_slice(std::string scol) const;

    /// Return a constant reference to a slice
    const ubmatrix &get_slice(size_t iz) const;

    /// Return a constant reference to a slice
    ubmatrix &get_slice(std::string scol);
    
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
    void lookup(double val, std::string slice, size_t &ix, size_t &iy) const;
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
    double interp(double x, double y, std::string name);

    /** \brief Interpolate the derivative of the data with respect to
	the x grid at point \c x and \c y in slice named \c name
    */
    double deriv_x(double x, double y, std::string name);

    /** \brief Interpolate the derivative of the data with respect to
	the y grid at point \c x and \c y in slice named \c name
    */
    double deriv_y(double x, double y, std::string name);

    /** \brief Interpolate the mixed second derivative of the data at
	point \c x and \c y in slice named \c name
    */
    double deriv_xy(double x, double y, std::string name);

    /** \brief Interpolate the integral of the data 
	respect to the x grid 
    */
    double integ_x(double x1, double x2, double y, std::string name);

    /** \brief Interpolate the integral of the data 
	respect to the y grid 
    */
    double integ_y(double x, double y1, double y2, std::string name);

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

    /** \brief Clear the table and the slice names
     */
    void clear_table();

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
    void summary(std::ostream *out, int ncol=79) const;
    //@}

    /// True if the size of the table has been set
    bool is_size_set() const {
      return size_set;
    }

    /// True if the grid has been set
    bool is_xy_set() const {
      return xy_set;
    }
  
    // ---------
    // Allow HDF I/O functions to access table3d data

    friend void o2scl_hdf::hdf_output(o2scl_hdf::hdf_file &hf, table3d &t, 
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
    /// Add a constant
    virtual void add_constant(std::string name, double val);

    /// Remove a constant
    virtual void remove_constant(std::string name);
    
    /// Add a constant
    virtual void set_constant(std::string name, double val);

    /// Get a constant
    virtual double get_constant(std::string name);

    /// Get a constant by index
    virtual void get_constant(size_t ix, std::string &name, double &val) const;

    /// Get the number of constants
    virtual size_t get_nconsts() const {
      return constants.size();
    }
    //@}

    /// Return the type, \c "table3d".
    virtual const char *type() { return "table3d"; }

    /** \name Parsing mathematical functions specified as strings
	
	\comment
	Note that doxygen doesn't format the documentation propertly
	if the \name specification covers more than one line
	\endcomment
     */
    //@{
    /** \brief Make a column from <tt>formula</tt> and add it to the table.
	
	If the column already exists, the data already present is 
	overwritten with the result.
    */
    void function_slice(std::string formula, std::string col);
    //@}

  protected:

#ifndef DOXYGEN_INTERNAL
    
    /** \brief Internal function to set function parser constants equal to 
	internal constants
    */
    void set_fp_consts(FunctionParser &fp) const {
      std::map<std::string,double>::const_iterator mit;
      for(mit=constants.begin();mit!=constants.end();mit++) {
	fp.AddConstant(mit->first,mit->second);
      }
      return;
    }

    /// \name Interpolation data
    //@{
    /// The array of interp_sm pointers 
    interp_vec<ubvector> **si;

    /// Matrices for interpolation
    ubmatrix_column **aci;
    //@}
    
    /// \name Iterator types
    //@{
    typedef std::map<std::string,size_t, string_comp>::iterator map_iter;
    typedef std::map<std::string,size_t, string_comp>::const_iterator 
      map_const_iter;
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
    std::map<std::string,size_t,string_comp> tree;

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

    size_t itype;

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
