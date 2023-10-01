/*
  ───────────────────────────────────────────────────────────────────
  
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
  
  You should have received a coxpy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  ───────────────────────────────────────────────────────────────────
*/
#ifndef O2SCL_TABLE_H
#define O2SCL_TABLE_H

/** \file table.h
    \brief File defining \ref o2scl::table
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include <map>

#include <o2scl/set_openmp.h>

#ifdef O2SCL_SET_OPENMP
#include <omp.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/misc.h>
#include <o2scl/interp_vec.h>
#include <o2scl/vec_stats.h>
#include <o2scl/lib_settings.h>

namespace o2scl {
  // Forward definition of the table class for HDF5 I/O
  template<class vec_t, class fp_t> class table;
  // Forward definition of matrix_view_table to extend
  // friendship
  template<class vec_t, class fp_t> class matrix_view_table;
  // Forward definition of matrix_view_table_transpose to extend
  // friendship
  template<class vec_t, class fp_t> class matrix_view_table_transpose;
}

// Forward definition of HDF5 I/O to extend friendship in table
namespace o2scl_hdf {
  // Forward definition of hdf_file for function template argument below
  class hdf_file;
  // Forward definition of HDF5 I/O functions to extend friendship
  template<class vec_t, class fp_t>
  void hdf_input(hdf_file &hf, o2scl::table<vec_t,fp_t> &t,
                 std::string name="");
  void hdf_output(hdf_file &hf, 
		  o2scl::table<std::vector<double>,double> &t, 
		  std::string name);
  template<class vec_t, class fp_t>
  void hdf_input_data(hdf_file &hf, o2scl::table<vec_t,fp_t> &t);
  void hdf_output_data(hdf_file &hf, 
		       o2scl::table<std::vector<double>,double> &t);

}

namespace o2scl {

  /** \brief Data \table class

      \b Summary \n 

      A class to contain and manipulate several equally-sized columns
      of data. The purpose of this class is to provide a structure
      which allows one to refer to the columns using a name
      represented by a string. Thus for a table object named \c t with
      3 columns (named "colx", "coly" and "colz") and three rows, one
      could do the following:
      \verbatim
      // Set the 1st row of column "colx" to 1.0
      t.set("colx",0,1.0);
      // Set the 2nd row of column "colz" to 2.0
      t.set("colz",1,2.0);
      // Set the 3rd row of column "coly" to 4.0
      t.set("coly",2,4.0);
      // This will print out 2.0
      cout << t.get("colz",1) << endl;
      \endverbatim
      Note that the rows are numbered starting with 0 instead of 
      starting with 1.
      To output all the rows of entire column, one can use
      \verbatim
      for(size_t i=0;i<t.get_nlines();i++) {
      cout << i << " " << t.get("colx",i) << endl;
      }
      \endverbatim
      To output all the columns of an entire row (in the following
      example it is the second row), labeled by their column name, one
      can use:
      \verbatim
      for(size_t i=0;i<t.get_ncolumns();i++) {
      cout << t.get_column_name(i) << " ";
      }
      cout << endl;
      for(size_t i=0;i<t.get_ncolumns();i++) {
      cout << t.get(i,1) << " ";
      }
      cout << endl;
      \endverbatim

      Methods are provided for interpolating columns, sorting 
      columns, finding data points, and several other manipulations
      of the data. 

      <B> Lookup, differentiation, integration, and 
      interpolation </b> \n

      Lookup, differentiation, integration, and interpolation are
      automatically implemented using splines from the class \ref
      interp_vec . A caching mechanism is implemented so that
      successive interpolations, derivative evaluations or
      integrations over the same two columns are fast.

      <B> Sorting </b>\n

      The columns are automatically sorted by name for speed, the
      results can be accessed from \ref get_sorted_name() . Individual
      columns can be sorted (\ref sort_column() ), or the entire table
      can be sorted by one column (\ref sort_table() ).

      <B> Data representation </b> \n

      Each individual column is just a vector object.
      The columns can be referred to in one of two ways:
      - A numerical index from 0 to C-1 (where C is the number of
      columns). For example, data can be accessed through \ref
      table::get() and \ref table::set(size_t c, size_t r,
      double val), or the overloaded [] operator, <tt>table[c][r]</tt>. 
      - A name of the column which is a string.
      For example, data can be accessed with table::get(string cname,
      int r) and table::set(string cname, int r, double val).
      
      The columns are organized in a both a \<map\> and a \<vector\>
      structure so that finding a column by its index, using either of
      \code
      std::string table::get_column_name(size_t index);
      ubvector &table::operator[](size_t index);
      \endcode
      takes only constant time, and finding a column by its name
      using either of
      \code
      size_t lookup_column(std::string name) const;
      const ubvector &get_column(std::string col) const;      
      \endcode
      is O(log(C)). Insertion of a column ( \ref new_column() ) is
      O(log(C)), but deletion ( \ref delete_column() ) is O(C). Adding
      a row of data can be either O(1) or O(C), but row insertion and
      deletion is slow, since the all of the rows must be shifted
      accordingly.
      
      Because of the structure, this class is not suitable for the
      matrix manipulation.

      <b>Vector types</b> \n

      The type <tt>vec_t</tt> can be any vector type with
      <tt>operator[]</tt>, <tt>size()</tt> and <tt>resize()</tt>
      methods. HDF5 I/O with vector types other than
      <tt>std::vector<double> </tt> requires a copy. 

      \verbatim embed:rst
      See the the discussion in the sections :ref:`Vector and Matrix
      Introduction` and :ref:`I/O and contiguous storage` of the
      User's Guide for more details.
      \endverbatim

      <b>Thread-safety</b> \n

      Generally, the member functions are only thread-safe 
      if they are <tt>const</tt> .  

      \b I/O \b and \b command-line \b manipulation \n

      When data from an object of type \ref table is output to a file
      through the <tt>hdf_output() function</tt> in \ref o2scl_hdf,
      the table can be manipulated on the command-line through the \c
      acol utility 

      \verbatim embed:rst
      See :ref:`The acol Command Line Utility` for more details
      on `acol`.
      \endverbatim

      There is an example for the usage of this class given
      in <tt>examples/ex_table.cpp</tt>.

      \verbatim embed:rst

      .. todo:: 

         In class table:

         - Future: Create a sort_column_names() or a function to 
           arbitrarily rearrange the columns

         - The present structure,
           std::map<std::string,col,string_comp> atree and \c
           std::vector<aiter> alist; could be replaced with \c
           std::vector<col> list and \c std::map<std::string,int> tree
           where the map just stores the index of the the column in
           the list.

      \endverbatim
  */
  template<class vec_t=std::vector<double>, class fp_t=double>
  class table {
    
  public:
  
    /// \name Constructors, destructors
    //@{
    /** \brief Create a new table with space for nlines<=cmaxlines.
     */
    table(size_t cmaxlines=0) {
      nlines=0;
      intp_set=false;
      maxlines=cmaxlines;
      itype=itp_cspline;
    }

    /** \brief Table destructor
     */
    virtual ~table() {
      if (intp_set==true) {
        delete si;
      }
    }

    /// Copy constructor
    table(const table &t) {
  
      // Copy constants
      constants=t.constants;

      // Copy interpolation type
      itype=t.itype;
    
      // Copy the columns and data
      nlines=t.get_nlines();
      maxlines=nlines;

      for(size_t i=0;i<t.get_ncolumns();i++) {

        // Column name
        std::string cname=t.get_column_name(i);

        // Insert column into tree
        col s;
        s.dat.resize(nlines);
        s.index=atree.size();
        atree.insert(make_pair(cname,s));

        // Insert in iterator index
        aiter it=atree.find(cname);
        alist.push_back(it);
    
        // Fill the data
        for(size_t j=0;j<t.get_nlines();j++) {
          it->second.dat[j]=t.get(cname,j);
        }
    
      }

      intp_set=false;

      is_valid();
    
      return;
    }

    /// Copy constructor
    table &operator=(const table &t) {

      if (this!=&t) {

        clear();
      
        // Copy constants
        constants=t.constants;
      
        // Copy the columns and data
        nlines=t.get_nlines();
        maxlines=nlines;
      
        // Copy interpolation type
        itype=t.itype;
      
        for(size_t i=0;i<t.get_ncolumns();i++) {
	
          // Column name
          std::string cname=t.get_column_name(i);
	
          // Insert column into tree
          col s;
          s.dat.resize(nlines);
          s.index=atree.size();
          atree.insert(make_pair(cname,s));
	
          // Insert in iterator index
          aiter it=atree.find(cname);
          alist.push_back(it);
	
          // Fill the data
          for(size_t j=0;j<t.get_nlines();j++) {
            it->second.dat[j]=t.get(cname,j);
          }
	
        }
      
        if (intp_set) {
          intp_set=false;
          delete si;
        }
      
      }

      is_valid();
    
      return *this;
    }

    /** \brief Swap method
     */
    friend void swap(table &t1, table &t2) {

      using std::swap;

      // The data
      swap(t1.maxlines,t2.maxlines);
      swap(t1.nlines,t2.nlines);
      swap(t1.atree,t2.atree);

      // Take care of interpolation
      swap(t1.itype,t2.itype);
      t1.intp_set=false;
      t2.intp_set=false;
      swap(t1.intp_colx,t2.intp_colx);
      swap(t1.intp_coly,t2.intp_coly);

      // Constants
      swap(t1.constants,t2.constants);

      // Recreate iterator lists
      t1.alist.resize(t1.atree.size());
      t1.reset_list();
      t2.alist.resize(t2.atree.size());
      t2.reset_list();

      // Check that it worked
      t1.is_valid();
      t2.is_valid();
    
      return;
    }
    //@}

    // --------------------------------------------------------
    /** \name Basic get and set methods */
    //@{
    /** \brief Set row \c row of column named \c col to value \c val .
        \f$ {\cal O}(\log(C)) \f$

        This function calls the error handler if the row is beyond
        the end of the table or if the specified column is not found.
    */
    void set(std::string scol, size_t row, fp_t val) {

      if (row>=nlines) {
        std::string str=((std::string)"Row ")+o2scl::szttos(row)+
          " beyond end of table (nlines="+o2scl::szttos(nlines)+") in "+
          "table::set(string,size_t,fp_t).";
        O2SCL_ERR(str.c_str(),exc_einval);
        return;
      }

      aiter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   "' not found in table::set(string,size_t,fp_t).").c_str(),
                  exc_enotfound);
        return;
      }

      if ((intp_colx==scol || intp_coly==scol) && intp_set==true) {
        delete si;
        intp_set=false;
      }

#if !O2SCL_NO_RANGE_CHECK
      if (row>=it->second.dat.size()) {
        O2SCL_ERR("Vector size failure in table::set(string,size_t,fp_t).",
                  exc_esanity);
      }
#endif
      it->second.dat[row]=val;
    
      return;
    }

    /** \brief Set row \c row of column number \c icol to value \c val .
        \f$ {\cal O}(1) \f$
    */
    void set(size_t icol, size_t row, fp_t val) {
    
      if (row>=nlines) {
        O2SCL_ERR2("Specified row beyond end of table in ",
                   "table::set(size_t,size_t,fp_t).",exc_einval);
        return;
      }
  
      if (icol>=atree.size()) {
        std::string err=((std::string)"Column index ")+szttos(icol)+
          ">="+szttos(atree.size())+", in table::set(size_t,size_t,fp_t).";
        O2SCL_ERR(err.c_str(),exc_einval);
      }

      std::string scol=get_column_name(icol);
      if ((intp_colx==scol || intp_coly==scol) && intp_set==true) {
        delete si;
        intp_set=false;
      }

#if !O2SCL_NO_RANGE_CHECK
      if (row>=alist[icol]->second.dat.size()) {
        std::string errs=((std::string)"Vector size failure, row ")+
          o2scl::szttos(row)+" >= "+
          o2scl::szttos(alist[icol]->second.dat.size())+
          "in table::set(size_t,size_t,fp_t).";
        O2SCL_ERR(errs.c_str(),exc_esanity);
      }
#endif
      alist[icol]->second.dat[row]=val;
      return;
    }

    /** \brief Set an entire row of data

        This function goes through \c v copying data until it runs out
        of columns in the table or it runs out of entries in \c v,
        whichever comes first.

        The type <tt>size_vec_t</tt> must be a type which has a
        <tt>size()</tt> method. 
    */
    template<class size_vec_t> void set_row(size_t row, size_vec_t &v) {
      if (row>=get_nlines()) {
        std::string err=((std::string)"Row out of range, ")+
          szttos(row)+">="+szttos(get_nlines())+", in table::set_row().";
        O2SCL_ERR(err.c_str(),exc_einval);
      }
      for(size_t i=0;i<get_ncolumns() && i<v.size();i++) {
        alist[i]->second.dat[row]=v[i];
      }
      return;
    }

    /** \brief Get value from row \c row of column named \c col.
        \f$ {\cal O}(\log(C)) \f$
    */
    fp_t get(std::string scol, size_t row) const {
      fp_t tmp;
      aciter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   "' not found in table::get(string,size_t).").c_str(),
                  exc_enotfound);
        return 0.0;
      } else {
        if (row>=nlines) {
          std::string err=((std::string)"Row out of range, ")+
            szttos(row)+">="+szttos(nlines)+
            ", in table::get(string,size_t).";
          O2SCL_ERR(err.c_str(),exc_einval);
        }
#if !O2SCL_NO_RANGE_CHECK
        if (row>=it->second.dat.size()) {
          O2SCL_ERR("Vector size failure in table::get().",
                    exc_esanity);
        }
#endif
        tmp=it->second.dat[row];
      }
      return tmp;
    }
    
    /** \brief Get value from row \c row of column number \c icol.
        \f$ {\cal O}(1) \f$
    */
    fp_t get(size_t icol, size_t row) const {
      if (icol>=atree.size()) {
        std::string err=((std::string)"Column out of range, ")+
          szttos(icol)+">="+szttos(atree.size())+
          ", in table::get(size_t,size_t).";
        O2SCL_ERR(err.c_str(),exc_einval);
      }
      if (row>=nlines) {
        std::string err=((std::string)"Row out of range, ")+
          szttos(row)+">="+szttos(nlines)+", in table::get(size_t,size_t).";
        O2SCL_ERR(err.c_str(),exc_einval);
      }
      return alist[icol]->second.dat[row];
    }

    /** \brief Return the number of columns
     */
    size_t get_ncolumns() const {return atree.size();};
    //@}

    // --------------------------------------------------------
    /** \name Manipulate current and maximum number of rows */
    //@{
    /** \brief Return the number of lines
     */
    size_t get_nlines() const {return nlines;};

    /** \brief Set the number of lines
	
        This function is stingy about increasing the table memory
        space and will only increase it enough to fit \c il lines.
        Using it in succession to slowly increase the number of lines
        in the table is likely to be inefficient compared to \ref
        set_nlines_auto() in this case.
    */
    void set_nlines(size_t il) {
      
      // Try to increase the number of lines
      if (il>maxlines) {
        inc_maxlines(il-maxlines);
      }
      
      // Now that maxlines is large enough, set the number of lines 
      nlines=il;
      
      // Reset the interpolation object for future interpolations
      if (intp_set) {
        intp_set=false;
        delete si;
      }
      
      return;
    }

    /** \brief Return the maximum number of lines before a reallocation
        is required
    */
    size_t get_maxlines() {return maxlines; };

    /** \brief Returns a copy of the row with value \c val in column
        \c col. \f$ {\cal O}(R C) \f$

        This function searches the entire table for the row which has
        the entry in column \c col which is closest to the value \c
        val, and copies that row to the vector \c row.

        If the object \c row previously contains any data, it will be
        lost.

        The type <tt>resize_vec_t</tt> must be a type which has
        <tt>size()</tt> and <tt>resize()</tt> methods.
    */
    template<class resize_vec_t>
    void get_row(std::string scol, fp_t val, resize_vec_t &row) const {
      
      int irow=lookup(scol,val);
      if (irow==exc_enotfound) {
        O2SCL_ERR((((std::string)"Column '")+scol+"' not found in "+
                   "table::get_row(string,fp_t,vec_t) const.").c_str(),
                  exc_enotfound);
        return;
      } 
      get_row(irow,row);
      return; 
    }
    
    /** \brief Returns a copy of row number \c irow. \f$ {\cal O}(C) \f$
	
        This function returns a copy of row with index \c irow,
        where \c irow ranges from 0 to <tt>get_nlines()-1</tt>,
        inclusive.

        If the object \c row previously contains any data, it will be
        lost.

        The type <tt>resize_vec_t</tt> must be a type which has
        <tt>size()</tt> and <tt>resize()</tt> methods.
    */
    template<class resize_vec_t>
    void get_row(size_t irow, resize_vec_t &row) const {
      
      if (irow+1>nlines) {
        O2SCL_ERR((((std::string)"Row '")+ szttos(irow)+
                   "' not found in table::get_row(size_t,vec_t).").c_str(),
                  exc_enotfound);
        return;
      }
      
      row.resize(atree.size());

      for(size_t i=0;i<get_ncolumns();i++) {
        row[i]=(alist[i]->second.dat)[irow];
      }
      
      return;
    }

    /** \brief Set the number of lines, increasing the size more 
        agressively

        This function is like set_nlines(), but doubles the maximum
        column size if an increase in the maximum size is required
        instead of simply making enough room for the current number of
        lines. This function is used internally by \ref set() to
        ensure that the cost of setting lines in sequence is linear
        and not quadratic.
    */
    void set_nlines_auto(size_t il) {
      
      // Try to increase the number of lines
      if (il>maxlines) {
        size_t inc=il-maxlines;
        if (inc<maxlines) inc=maxlines;
        inc_maxlines(inc);
      }
      
      // Now that maxlines is large enough, set the number of lines 
      nlines=il;
      
      // Reset the interpolation object for future interpolations
      if (intp_set) {
        intp_set=false;
        delete si;
      }
      
      return;
    }

    /** \brief Manually increase the maximum number of lines
     */
    void inc_maxlines(size_t llines) {

      vec_t temp_col;
    
      // For the moment, we assume resizes are destructive, so
      // we have to copy the data to a temporary and then
      // copy it back
      for(aiter it=atree.begin();it!=atree.end();it++) {
	
        // Copy data to temporary array
        temp_col.resize(maxlines+llines);
        for(size_t j=0;j<maxlines;j++) {
          temp_col[j]=it->second.dat[j];
        }

        // Resize
        it->second.dat.resize(maxlines+llines);

        // Copy data back to resized array
        for(size_t j=0;j<maxlines;j++) {
          it->second.dat[j]=temp_col[j];
        }
	
      }
  
      maxlines+=llines;

      return;
    }

    /** \brief Manually set the maximum number of lines

        \note This function will call the error handler if
        the argument <tt>llines</tt> is smaller than the
        current number of lines in the table.
    */
    void set_maxlines(size_t llines) {

      if (llines==maxlines) return;
      if (llines<nlines) {
        O2SCL_ERR2("Cannot set maximum number of lines to be smaller ",
                   "than current size in table::set_maxlines().",
                   o2scl::exc_einval);
		 
      }
    
      vec_t temp_col;
    
      // For the moment, we assume resizes are destructive, so
      // we have to copy the data to a temporary and then
      // copy it back
      for(aiter it=atree.begin();it!=atree.end();it++) {
	
        // Copy data to temporary array
        temp_col.resize(llines);
        for(size_t j=0;j<nlines;j++) {
          temp_col[j]=it->second.dat[j];
        }

        // Resize
        it->second.dat.resize(llines);

        // Copy data back to resized array
        for(size_t j=0;j<nlines;j++) {
          it->second.dat[j]=temp_col[j];
        }
	
      }
  
      maxlines=llines;

      return;
    }
    //@}

    // --------------------------------------------------------
    /** \name Column manipulation */
    //@{
    /** \brief Returns a reference to the column named \c col.
        \f$ {\cal O}(\log(C)) \f$

        Note also that the vector object returned by this function may
        be larger than the number of rows in the table, as the table
        resizes these vectors automatically to make room for more data.
        To get the number of valid data entries in the object, 
        use \ref get_nlines() instead of <tt>get_column().size()</tt>.
    */
    const vec_t &get_column(std::string scol) const {
      aciter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   "' not found in table::get_column() const.").c_str(),
                  exc_enotfound);
        return empty_col;
      }
      return it->second.dat;
    }
    
    /** \brief Returns the column of index \c icol (const
        version). \f$ {\cal O}(1) \f$

        Note that several of the methods require reallocation of
        memory and refereces previously returned by this function will
        be incorrect.

        Note also that the vector object returned by this function may
        be larger than the number of rows in the table, as the table
        resizes these vectors automatically to make room for more data.
        To get the number of valid data entries in the object, 
        use \ref get_nlines() instead of <tt>operator[].size()</tt>.

        Unlike set(), this function will not automatically result in
        an increase in the size of the table if the user attempts to
        set an element beyond the current column range.

        This function will throw an exception if \c icol is out
        of range unless <tt>O2SCL_NO_RANGE_CHECK</tt> is defined.
    */
    const vec_t &operator[](size_t icol) const {
#if !O2SCL_NO_RANGE_CHECK
      if (icol>=atree.size()) {
        O2SCL_ERR((((std::string)"Array index ")+szttos(icol)+
                   " out of bounds"+
                   " in table::operator[size_t] const. Size: "+
                   szttos(atree.size())+
                   " (index should be less than size).").c_str(),exc_eindex);
        return empty_col;
      }
#endif
      return (alist[icol]->second.dat);
    }
    
    /** \brief Returns the column named \c scol (const version).
        \f$ {\cal O}(\log(C)) \f$

        Note that several of the methods require reallocation of
        memory and refereces previously returned by this function will
        be incorrect.

        Note also that the vector object returned by this function may
        be larger than the number of rows in the table, as the table
        resizes these vectors automatically to make room for more data.
        To get the number of valid data entries in the object, 
        use \ref get_nlines() instead of <tt>operator[].size()</tt>.

        This function will throw an exception if \c icol is out
        of range unless <tt>O2SCL_NO_RANGE_CHECK</tt> is defined.
    */
    const vec_t &operator[](std::string scol) const {
      aciter it=atree.find(scol);
#if !O2SCL_NO_RANGE_CHECK
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+"' not found in table::"+
                   "operator[string] const.").c_str(),exc_enotfound);
        return empty_col;
      }
#endif
      return (it->second.dat);
    }

    /** \brief Add a new column owned by the \table \f$ {\cal O}(\log(C)) \f$

        \note This function does not set all the column entries to
        zero in the case that a new column is added to a table which
        already contains data.
    */
    void new_column(std::string head) {
      if (head.length()==0) {
        O2SCL_ERR2("Cannot add column with empty name in ",
                   "table::new_column().",exc_einval);
      }
      if (is_column(head)==true) {
        O2SCL_ERR((((std::string)"Column '")+head+
                   "' already present in table::new_column().").c_str(),
                  exc_einval);
      }
      for(int i=0;i<((int)head.size());i++) {
        if (head[i]==' ' || head[i]=='\t' || head[i]=='\n' || head[i]=='\r' 
            || head[i]=='\v' || head[i]=='\f') {
          O2SCL_ERR((((std::string)"Invalid column name '")+head+
                     "' in table::new_column().").c_str(),
                    exc_einval);
        }
      }
      col s;
      s.dat.resize(maxlines);
      s.index=((int)atree.size());
      atree.insert(make_pair(head,s));
      aiter it=atree.find(head);
      alist.push_back(it);
      return;
    }
  
    /** \brief Add a new column by copying data from another vector
      
        This function copies \c sz elements of vector \c v into the
        table in a new column named \c name. If \c sz is larger than
        the current number of lines (as given, e.g. in \ref
        get_nlines() ), then only the first part of the vector \c v is
        copied, up to the current number of lines.
      
        This function calls the error handler if \c sz is zero.

        The type <tt>vec2_t</tt> can be any type with an
        <tt>operator[]</tt> method.
    */
    template<class vec2_t> int new_column(std::string name, 
                                          size_t sz, vec2_t &v) {
    
      if (sz==0) {
        O2SCL_ERR2("Sent column of zero size in ",
                   "table::new_column(string,size_t,vec2_t)",
                   exc_einval);
      }
    
      // Create the space
      int ret=new_column(name);
      if (ret!=0) return ret;

      // Copy the data over
      size_t mxl=sz;
      if (sz>get_nlines()) mxl=get_nlines();
      for(size_t i=0;i<mxl;i++) {
        set(name,i,v[i]);
      }
	
      return 0;
    }

    /** \brief Returns the name of column \c col \f$ {\cal O}(1) \f$

        This will throw if \c icol is larger than or equal to
        the number of columns.
    */
    std::string get_column_name(size_t icol) const {
      if (icol+1>atree.size()) {
        O2SCL_ERR((((std::string)"Index '")+o2scl::szttos(icol)+
                   " larger than number of "+
                   "columns in table::get_column_name().").c_str(),
                  exc_enotfound);
      }
      return alist[icol]->first;
    }
  
    /** \brief Swap the data in column \c scol with that in vector \c v

        This requires that the column \c v must have the correct size,
        that returned by \ref get_maxlines().
      
        This function is useful, in part, because if objects of type
        <tt>vec_t</tt> have <tt>std::move</tt> defined, then the swap
        doesn't require a full copy.
    */
    virtual void swap_column_data(std::string scol, vec_t &v) {
      aiter its=atree.find(scol);
      if (its==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   " not found in table::swap_column_data().").c_str(),
                  exc_enotfound);
        return;
      }
      if (v.size()!=its->second.dat.size()) {
        O2SCL_ERR2("Vector sizes not commensurate in ",
                   "table::swap_column_data().",exc_einval);
      }
      std::swap(its->second.dat,v);
      return;
    }

    /** \brief Rename column named \c src to \c dest
        \f$ {\cal O}(C) \f$
    */
    virtual void rename_column(std::string src, std::string dest) {
      aiter its=atree.find(src);
      if (its==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+src+
                   " not found in table::rename_column().").c_str(),
                  exc_enotfound);
        return;
      }
      new_column(dest);
      aiter itd=atree.find(dest);
      std::swap(its->second.dat,itd->second.dat);
      delete_column(src);
      return;
    }

    /** \brief Delete column named \c scol \f$ {\cal O}(C) \f$
    
        This is slow because the iterators in \ref alist are mangled
        and we have to call \ref reset_list() to get them back.
    */
    virtual void delete_column(std::string scol) {
      
      // Find the tree iterator for the element we want to erase
      aiter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   " not found in table::delete_column().").c_str(),
                  exc_enotfound);
        return;
      }
        
      if (true) {

        int ix_match=it->second.index;
        
        for(aiter it2=atree.begin();it2!=atree.end();it2++) {
          if (it2->second.index>ix_match) {
            it2->second.index=it2->second.index-1;
          }
        }
        
        // Erase the elements from the list and the tree
        atree.erase(it);

        
        // Resize the list to the correct size
        alist.resize(atree.size());
        
      } else {

        // Find the corresponding list iterator
        aviter vit=alist.begin();
        vit+=it->second.index;
        
        // Find the last element in the list and it's corresponding table
        // entry. Change it's index to the index of the column to be
        // deleted.
        alist[alist.size()-1]->second.index=it->second.index;
        
        // Erase the elements from the list and the tree
        atree.erase(it);
        alist.erase(vit);

      }
      
      // Reset the list to reflect the proper iterators
      reset_list();
      
      if ((intp_colx==scol || intp_coly==scol) && intp_set==true) {
        delete si;
        intp_set=false;
      }

      is_valid();

      return;
    }

    /** \brief Returns the name of column \c col in sorted order.
        \f$ {\cal O}(1) \f$
    */
    std::string get_sorted_name(size_t icol) const {
      if (icol+1>atree.size()) {
        return "";
      }
      aciter it=atree.begin();
      for(size_t i=0;i<icol;i++) it++;
      return it->first;
    }

    /** \brief Initialize all values of column named \c scol to \c val 
        \f$ {\cal O}(R \log(C)) \f$

        Note that this does not initialize elements beyond nlines so
        that if the number of rows is increased afterwards, the new
        rows will have uninitialized values.
    */
    void init_column(std::string scol, fp_t val) {
      /*
        if (!std::isfinite(val)) {
        O2SCL_ERR((((std::string)"Value '")+dtos(val)+
        "' not finite for column '"+
        scol+"' in table::init_column()").c_str(),exc_einval);
        return;
        }
      */
      aiter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   "' not found in table::init_column()").c_str(),
                  exc_enotfound);
        return;
      }
      for(size_t i=0;i<nlines;i++) {
        it->second.dat[i]=val;
      }
      if (intp_set && (scol==intp_colx || scol==intp_coly)) {
        intp_set=false;
        delete si;
      }
      return;
    }

    /** \brief Return true if \c scol is a column in the current \table

        This function does not call the error handler if the column is
        not found, but just silently returns false.
    */
    bool is_column(std::string scol) const {
      aciter it=atree.find(scol);
      if (it==atree.end()) return false;
      return true;
    }

    /** \brief Find the index for column named \c name 
        \f$ {\cal O}(\log(C)) \f$

        If the column is not present, this function calls the error
        handler.
    */
    size_t lookup_column(std::string lname) const {
      aciter it=atree.find(lname);
      if (it==atree.end()) {
        O2SCL_ERR("Column not found in table::lookup_column().",
                  exc_enotfound);
      }
      return it->second.index;
    }

    /** \brief Copy data from column named \c src to column 
        named \c dest, creating a new column if necessary
        \f$ {\cal O}(R \log(C)) \f$
    */
    virtual void copy_column(std::string src, std::string dest) {
      if (!is_column(dest)) new_column(dest);
      aiter its=atree.find(src);
      if (its==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+src+
                   " not found in table::copy_column().").c_str(),
                  exc_enotfound);
        return;
      }
      aiter itd=atree.find(dest);
      if (itd==atree.end()) {
        O2SCL_ERR((((std::string)"Destination column '")+dest+
                   " not found in table::copy_column().").c_str(),
                  exc_esanity);
        return;
      }
      for(size_t i=0;i<nlines;i++) {
        itd->second.dat[i]=its->second.dat[i];
      }
      return;
    }

    /** \brief Copy a column to a generic vector object

        \note It is assumed that the vector type is one that can be
        resized with <tt>resize()</tt>.
    */
    template<class resize_vec_t> 
    void column_to_vector(std::string scol, resize_vec_t &v) const {
      v.resize(nlines);
      for(size_t i=0;i<nlines;i++) {
        v[i]=this->get(scol,i);
      }
      return;
    }

    /** \brief Copy to a column from a generic vector object

        The type <tt>vec2_t</tt> can be any type with an
        <tt>operator[]</tt> method.
    */
    template<class vec2_t> 
    void copy_to_column(vec2_t &v, std::string scol) {

      aiter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   "' not found in table::copy_to_column(0.").c_str(),
                  exc_enotfound);
        return;
      }

      if ((intp_colx==scol || intp_coly==scol) && intp_set==true) {
        delete si;
        intp_set=false;
      }

      for(size_t i=0;i<nlines;i++) {
        it->second.dat[i]=v[i];
      }
    
      return;
    }

    /** \brief Insert a column from a separate table, interpolating
        it into a new column
	
        Given a pair of columns ( \c src_index, \c src_col ) in a
        separate table (\c source), this creates a new column in the
        present table named \c src_col which interpolates \c loc_index
        into \c src_index.  The interpolation objects from the \c
        source table will be used. If there is already a column in the
        present table named \c src_col, then this will fail.
    */
    template<class vec2_t>
    void add_col_from_table(table<vec2_t> &source,
                            std::string src_index, std::string src_col,
                            std::string dest_index="",
                            std::string dest_col="") {
    
      if (dest_col=="") dest_col=src_col;
      if (dest_index=="") dest_index=src_index;

      // Add the new column
      if (!is_column(dest_col)) new_column(dest_col);
  
      // Fill the new column
      for(size_t i=0;i<nlines;i++) {
        set(dest_col,i,source.interp(src_index,get(dest_index,i),src_col));
      }
  
      return;
    }
    //@}

    /** \brief Insert columns from a source table into the new
        table by interpolation (or extrapolation)

        This takes all of the columns in \c source, and adds them into
        the current table using interpolation, using the columns \c
        src_index and \c dest_index as the independent variable. The
        column named \c src_index is the column of the independent
        variable in \c source and the column named \c dest_index
        is the column of the independent variable in the current table.
        If \c dest_index is empty (the default) then the names in the
        two tables are taken to be the same.

        If necessary, columns are created in the current table for the
        dependent variable columns in \c source. Columns in the current
        table which do not correspond to dependent variable columns in
        \c source are left unchanged.

        If \c allow_extrap is false, then extrapolation is not allowed,
        and rows in the current table which have values of the independent
        variable which are outside the source table are unmodified. 

        If a column for a dependent variable in \c source has the
        same name as \c dest_index, then it is ignored and not inserted
        into the current table.

        If the column named \c src_index cannot be found in 
        \c source or the column names \c dest_index cannot be found
        in the current table, then the error handler is called.

        If the \c allow_extrap is false and either the minimum or
        maximum values of the column named \c src_index in the \c source
        table are not finite, then the error handler is called.
    */
    template<class vec2_t>
    void insert_table(table<vec2_t> &source, std::string src_index,
                      bool allow_extrap=true, std::string dest_index="") {

      if (dest_index=="") dest_index=src_index;

      if (!source.is_column(src_index)) {
        O2SCL_ERR("Source indep. var. column not found.",o2scl::exc_einval);
      }
      if (!is_column(dest_index)) {
        O2SCL_ERR("Dest. indep. var. column not found.",o2scl::exc_einval);
      }

      // Find limits to avoid extrapolation if necessary
      fp_t min=source.min(src_index);
      fp_t max=source.max(src_index);
      if (allow_extrap==false) {
        if (!std::isfinite(min) || !std::isfinite(max)) {
          O2SCL_ERR2("Minimum or maximum of source index not finite ",
                     "in table::insert_table().",exc_einval);
        }
      }

      // Create list of columns to interpolate
      std::vector<std::string> col_list;
      for(size_t i=0;i<source.get_ncolumns();i++) {
        std::string col_name=source.get_column_name(i);
        if (col_name!=src_index && col_name!=dest_index) {
          col_list.push_back(col_name);
        }
      }

      // Create new columns and perform interpolation
      for(size_t i=0;i<col_list.size();i++) {
        if (!is_column(col_list[i])) new_column(col_list[i]);
        for(size_t j=0;j<get_nlines();j++) {
          fp_t val=get(dest_index,j);
          if (allow_extrap || (val>=min && val<=max)) {
            set(col_list[i],j,source.interp(src_index,val,col_list[i]));
          }
        }
      }
    
      return;
    }

    /** \brief Add a second table to this one, adding all 
        constants, columns and data from the source
    */
    template<class vec2_t>
    void add_table(table<vec2_t> &source) {

      // Add constants
      for(size_t i=0;i<source.get_nconsts();i++) {
        std::string tnam;
        fp_t tval;
        source.get_constant(i,tnam,tval);
        add_constant(tnam,tval);
      }

      // Increase the size of the new table
      size_t n1=get_nlines();
      size_t n2=source.get_nlines();
      set_nlines(n1+n2);

      // Initialize new rows to zero for columns already in the table
      for(size_t i=0;i<get_ncolumns();i++) {
        for(size_t k=n1;k<n1+n2;k++) {
          set(i,k,0.0);
        }
      }
    
      // Go through each column in the source
      for(size_t j=0;j<source.get_ncolumns();j++) {
      
        std::string col_name=source.get_column_name(j);

        // If the column is not present, add it and
        // initialize it to zero
        if (!is_column(col_name)) {
          new_column(col_name);
          for(size_t i=0;i<n1+n2;i++) set(col_name,i,0.0);
        }

        // Copy the source data over for this column
        for(size_t i=0;i<n2;i++) {
          set(col_name,i+n1,source.get(col_name,i));
        }
      }
      
      return;
    }
    
    // --------------------------------------------------------
    /** \name Row maninpulation and data input */
    //@{

    /** \brief Insert a row before row \c n 

        Acceptable values for \c n are between 0 and
        <tt>get_nlines()</tt> inclusive, with the maximum value
        denoting the addition of a row after the last row presently in
        the table.
    */
    void new_row(size_t n) {

      // If we're already at the maximum number of lines,
      // double it so that we can easily add more data later
      if (maxlines==0) inc_maxlines(1);
      else if (nlines>=maxlines) inc_maxlines(maxlines);

      // Increase the nlines parameter
      nlines++;

      // Shift the data if necessary. Note that if n is equal to the
      // original value of nlines, then n==nlines-1, thus i=nlines-2
      // means that i<n and thus the loop is bypassed altogether.
      for(int i=((int)nlines)-2;i>=((int)n);i--) {
        copy_row(i,i+1);
      }

      // Reset the interpolation objects
      if (intp_set) {
        intp_set=false;
        delete si;
      }

      return;
    }

    /** \brief Copy the data in row \c src to row \c dest 

        \comment
        AWS 4/12/18: No need to throw an exception here since get() and
        set() will take care of that already.
        \endcomment
    */
    void copy_row(size_t src, size_t dest) {
      for(int i=0;i<((int)atree.size());i++) {
        set(i,dest,get(i,src));
      }
      if (intp_set) {
        intp_set=false;
        delete si;
      }
      return;
    }

    /** \brief Delete the row with the entry closest to
        the value \c val in column \c scol \f$ {\cal O}(R C) \f$
    */
    void delete_row(std::string scol, fp_t val) {
      // If lookup() fails, it will throw an exception,
      // so there's no need to double-check it here
      size_t irow=lookup(scol,val);
      delete_row(irow);
      return;
    }

    /** \brief Delete the row of index \c irow \f$ {\cal O}(R C) \f$
     */
    void delete_row(size_t irow) {
      if (nlines==0) {
        O2SCL_ERR2("No lines in table in ",
                   "table::delete_row(size_t).",o2scl::exc_einval);
      }
      if (irow>=nlines) {
        std::string str=((std::string)"Cannot delete row ")+
          o2scl::szttos(irow)+" since there are only "+
          o2scl::szttos(nlines)+" lines in the table in "+
          "table::delete_row(size_t).";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      for(aiter it=atree.begin();it!=atree.end();it++) {
        // Can't do size_t because we have to compare to nlines-1
        for(int i=((int)irow);i<((int)nlines)-1;i++) {
          it->second.dat[i]=it->second.dat[i+1];
        }
      }
      nlines--;
      if (intp_set==true) {
        delete si;
        intp_set=false;
      }
      return;
    }

    /** \brief Delete all rows where \c func evaluates to a number greater
        than 0.5 \f$ {\cal O}(R C) \f$

        If no rows match the delete condition, this function silently
        performs no changes to the table.
    */
    void delete_rows_func(std::string func) {
      size_t new_nlines=0;
      for(size_t i=0;i<nlines;i++) {
        fp_t val=row_function(func,i);
        if (val<0.5) {
          // If val<0.5, then the function was evaluated to false and
          // we want to keep the row, but if i==new_nlines, then
          // we don't need to copy because the row is already in
          // the correct place.
          if (i!=new_nlines) {
            for(aiter it=atree.begin();it!=atree.end();it++) {
              it->second.dat[new_nlines]=it->second.dat[i];
            }
          }
          new_nlines++;
        }
      }
      nlines=new_nlines;
      if (intp_set==true) {
        delete si;
        intp_set=false;
      }
      return;
    }

    /** \brief Copy row \c ix from table \c src to the end of the
        current table
    */
    template<class vec2_t>
    void copy_row(table<vec2_t> &src, size_t ix) {
      if (ix>=src.get_nlines()) {
        O2SCL_ERR2("Row not present in ",
                   "table::copy_row().",o2scl::exc_einval);
      }
      size_t n=this->get_nlines();
      // Add one new row at the end of the table
      this->new_row(n);
      for(size_t i=0;i<src.get_ncolumns() && i<this->get_ncolumns();i++) {
        this->set(i,n,src[i][ix]);
      }
      // The interpolation objects need to be reset, but this
      // is already done in new_row() above.
      return;
    }
  
    /** \brief Copy all rows matching a particular condition to
        a new table

        This function begins by ensuring that all columns in the current
        table are present in \c dest, creating new columns in \c dest if
        necessary. It then copies all rows where \c func evaluates to a
        number greater than 0.5 to table \c dest by adding rows at
        the end of the table.
    */
    template<class vec2_t>
    void copy_rows(std::string func, table<vec2_t> &dest) {

      // Set up columns
      for(size_t i=0;i<get_ncolumns();i++) {
        std::string cname=get_column_name(i);
        if (dest.is_column(cname)==false) {
          dest.new_column(cname);
        }
      }

      size_t new_lines=dest.get_nlines();
      for(size_t i=0;i<nlines;i++) {
        fp_t val=row_function(func,i);
        if (val>0.5) {
          dest.set_nlines_auto(new_lines+1);
          for(size_t j=0;j<get_ncolumns();j++) {
            std::string cname=get_column_name(j);
            dest.set(cname,new_lines);
          }
          new_lines++;
        }
      }

      return;
    }

    /** \brief Delete all rows between \c row_start and \c row_end
        \f$ {\cal O}(R C) \f$

        If <tt>row_start</tt> is less or equal to <tt>row_end</tt>,
        then all rows beginnning with <tt>row_start</tt> and
        ending with <tt>row_end</tt> are deleted (inclusive). If
        <tt>row_start</tt> is greater than <tt>row_end</tt>, then
        rows from the start of the table until <tt>row_end</tt> are
        deleted, as well as all rows from <tt>row_start</tt>
        through the end of the table.

        If either <tt>row_start</tt> or <tt>row_end</tt> are beyond the
        end of the table (greater than or equal to the value given by
        \ref get_nlines() ), an exception is thrown.
    */
    void delete_rows_ends(size_t row_start, size_t row_end) {
      if (row_start>=nlines || row_end>=nlines) {
        O2SCL_ERR2("Row specifications beyond end of table in ",
                   "table::delete_rows_ends(size_t,size_t).",exc_einval);
      }
      size_t new_nlines=0;
      for(size_t i=0;i<nlines;i++) {
        if ((row_start<=row_end && (i<row_start || i>row_end)) ||
            (row_start>row_end && (i<row_start && i>row_end))) {
          for(aiter it=atree.begin();it!=atree.end();it++) {
            it->second.dat[new_nlines]=it->second.dat[i];
          }
          new_nlines++;
        }
      }
      nlines=new_nlines;
      if (intp_set==true) {
        delete si;
        intp_set=false;
      }
      return;
    }
  
    /** \brief Delete all rows in a specified list
      
        Given a list of rows in \c row_list, this function deletes all
        of the specified rows. If a row beyond the end of the table is
        in the list, the error handler is called.

        \note This function will proceed normally if a row is specified
        more than once in <tt>row_list</tt>.
    */
    template<class vec_size_t> 
    void delete_rows_list(vec_size_t &row_list) {

      // First, check that they're all valid rows
      for(size_t j=0;j<row_list.size();j++) {
        if (row_list[j]>nlines) {
          O2SCL_ERR("Invalid row in table<>::delete_rows_list(vec_size_t &)",
                    o2scl::exc_einval);
        }
      }

      // Copy the data over, ensuring rows which are to be
      // deleted are skipped
      size_t new_nlines=0;
      for(size_t i=0;i<nlines;i++) {
        bool found=false;
        for(size_t j=0;j<row_list.size();j++) {
          if (row_list[j]==i) found=true;
        }
        if (found==false) {
          for(aiter it=atree.begin();it!=atree.end();it++) {
            it->second.dat[new_nlines]=it->second.dat[i];
          }
          new_nlines++;
        }
      }
    
      // Set the new line number and reset the interpolator
      nlines=new_nlines;
      if (intp_set==true) {
        delete si;
        intp_set=false;
      }
    
      return;
    }

    /** \brief Exaustively search for groups of rows which match within a
        specified tolerance and remove all but one of each group

        For each column, the entries in the two rows do not match if 
        either of their absolute values are greater the absolute
        tolerance and their relative deviation is greater than
        the relative tolerance.

        This function returns the number of rows deleted.
    */
    size_t delete_rows_tolerance(double tol_rel=1.0e-12,
                                 double tol_abs=1.0e-20,
                                 int verbose=0) {
      std::vector<size_t> list;
      for(size_t i=0;i<nlines;i++) {
        for(size_t j=i+1;j<nlines;j++) {
          bool match=true;
          if (i<nlines && j<nlines && j>i) {
            for(aiter it=atree.begin();it!=atree.end() && match==true;it++) {
              if (fabs(it->second.dat[i])>tol_abs ||
                  fabs(it->second.dat[j])>tol_abs) {
                if (fabs(it->second.dat[i]-it->second.dat[j])/
                    fabs(it->second.dat[i]+it->second.dat[j])>tol_rel) {
                  match=false;
                }
              }
            }
          }
          if (match==true) {
            
            list.push_back(j);
            
            // AWS, 1/3/22: The verbose parameter here is from acol,
            // and I don't want this output to occur in acol when
            // verbose=1, so I require verbose to be 2 or larger in
            // order for this output to show up.
            if (verbose>1) {
              std::cout << "Match between rows " << i << " and " << j
                        << std::endl;
              if (verbose>2) {
                for(size_t k=0;k<get_ncolumns();k++) {
                  std::cout << k << " " << get(k,i) << " " << get(k,j)
                            << std::endl;
                }
              }
            }
          }
        }
      }
      if (list.size()>0) {
        delete_rows_list(list);
      } else if (verbose>0) {
        std::cout << "No matches found." << std::endl;
      }
      return list.size();
    }
  
    /** \brief Delete all rows which are identical to
        adjacent rows

        This function does silently does nothing if there are
        less than 2 rows in the table.
    */
    void delete_idadj_rows() {
      if (nlines>=2) {
        // Find duplicate rows
        std::vector<size_t> row_list;
        for(size_t i=0;i<nlines-1;i++) {
          bool match=true;
          for(aiter it=atree.begin();it!=atree.end() && match==true;it++) {
            if (it->second.dat[i+1]!=it->second.dat[i]) match=false;
          }
          if (match) row_list.push_back(i+1);
        }
        // Delete duplicates
        delete_rows_list(row_list);
      }
      return;
    }
  
    /** \brief Read a new set of names from \c newheads

        This function reads a set of white-space delimited column
        names from the string \c newheads, and creates a new column
        for each name which is specified.
	
        For example
        \code
        table t;
        t.line_of_names("position velocity acceleration");
        \endcode
        will create three new columns with the names "position",
        "velocity", and "acceleration".
    */
    void line_of_names(std::string newheads) {
      int ret=0;
      std::string head;

      std::istringstream is(newheads);
      while(is >> head) {
        new_column(head);
      } 
  
      if (ret!=0) {
        O2SCL_ERR2("At least one new column failed in ",
                   "table::line_of_names().",exc_efailed);
      }

      return;
    }

    /** \brief Read a line of data from the first \c nv entries in 
        a vector and store as a new row in the table

        The type <tt>vec2_t</tt> can be any type with an
        <tt>operator[]</tt> method. Note that this function does not
        verify that \c nv is equal to the number of columns, so some of
        the columns may be uninitialized in the new row which is
        created.

        Similar to <tt>std::vector</tt>'s <tt>push_back()</tt> method,
        this function now internally increases the maximum table size
        geometrically to help avoid excessive memory rearrangements.
    */
    template<class vec2_t> void line_of_data(size_t nv, const vec2_t &v) {
      if (maxlines==0) inc_maxlines(1);
      if (nlines>=maxlines) inc_maxlines(maxlines);
    
      if (intp_set) {
        intp_set=false;
        delete si;
      }
      
      if (nlines<maxlines && nv<=(atree.size())) {

        set_nlines_auto(nlines+1);
        for(size_t i=0;i<nv;i++) {
          (*this).set(i,nlines-1,v[i]);
        }
	
        return;
      }

      O2SCL_ERR("Not enough lines or columns in line_of_data().",exc_einval);
      return;
    }

    /** \brief Insert \c nv entries from \c v as a new row in 
        the table at row \c row

        \warning This function potentially requires reallocating
        and copying a portion of the entire table.
    */
    template<class vec2_t> void insert_row(size_t nv, const vec2_t &v,
                                           size_t row) {

      if (row>nlines) {
        O2SCL_ERR("Row greater than nlines in insert_row().",
                  o2scl::exc_einval);
      }
      
      if (maxlines==0) inc_maxlines(1);
      if (nlines>=maxlines) inc_maxlines(maxlines);
    
      if (intp_set) {
        intp_set=false;
        delete si;
      }
      
      if (nlines<maxlines && nv<=(atree.size())) {
        
        set_nlines_auto(nlines+1);
        
        if (nlines>=2 && nlines-2>=row) {
          for(int i=((int)nlines-2);i>=((int)row);i--) {
            for(size_t j=0;j<atree.size();j++) {
              (*this).set(j,i+1,(*this).get(j,i));
            }
          }
        }

        for(size_t i=0;i<nv;i++) {
          (*this).set(i,row,v[i]);
        }
        
        return;
      }

      O2SCL_ERR("Not enough lines or columns in insert_row().",exc_einval);
      return;
    }

    /** \brief Read a line of data and store in a new row of the
        table

        The type <tt>vec2_t</tt> can be any type with an
        <tt>operator[]</tt> method. Note that this function does not
        verify that the vector size is equal to the number of columns,
        so some of the columns may be uninitialized in the new row which
        is created.
    */
    template<class vec2_t> void line_of_data(const vec2_t &v) {
      line_of_data(v.size(),v);
      return;
    }
    //@}

    // --------------------------------------------------------
    /** \name Lookup and search methods */
    //@{

    /** \brief Look for a value in an ordered column 
        \f$ {\cal O}(\log(C) \log(R)) \f$

        This uses the function search_vec::ordered_lookup(), which
        offers caching and assumes the vector is monotonic. If you
        don't have monotonic data, you can still use the
        table::lookup() function, which is more general.
    */
    size_t ordered_lookup(std::string scol, fp_t val) const {
      int ret;
      if (!std::isfinite(val)) {
        O2SCL_ERR((((std::string)"Value '")+dtos(val)+
                   "' not finite for column '"+
                   scol+"' in table::ordered_lookup()").c_str(),exc_einval);
        return exc_einval;
      }
      aciter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   " not found in table::ordered_lookup().").c_str(),
                  exc_enotfound);
        return exc_enotfound;
      }

      search_vec<vec_t,fp_t> se(nlines,it->second.dat);
      ret=se.ordered_lookup(val);
      return ret;
    }

    /** \brief Exhaustively search column \c col for the value \c val
        \f$ {\cal O}(R \log(C)) \f$
    */
    size_t lookup(std::string scol, fp_t val) const {
      if (!std::isfinite(val)) {
        O2SCL_ERR((((std::string)"Value '")+dtos(val)+
                   "' not finite for column '"+
                   scol+"' in table::lookup()").c_str(),exc_einval);
        return exc_einval;
      }
      aciter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+" not found in "+
                   "table::lookup().").c_str(),exc_enotfound);
        return exc_enotfound;
    
      }

      // Note that we cannot use the vector lookup() method here, because
      // the vector size may be larger than the actual table size.

      const vec_t &ov=it->second.dat;
      size_t row=0, i=0;

      // Find first finite row
      while(!std::isfinite(ov[i]) && i<nlines-1) i++;
      if (i==nlines-1) {
        O2SCL_ERR2("Entire array not finite in ",
                   "table::lookup()",exc_einval);
        return 0;
      }

      // Beginning with that row, look for the closest value
      fp_t bdiff=fabs(ov[i]-val);
      for(;i<nlines;i++) {
        if (std::isfinite(ov[i]) && fabs(ov[i]-val)<bdiff) {
          row=i;
          bdiff=fabs(ov[i]-val);
        }
      }
  
      return row;
    }

    /// Search column \c col for the value \c val and return value in \c col2
    fp_t lookup_val(std::string scol, fp_t val, std::string scol2) const {
      int i, indx=0;
      if (!std::isfinite(val)) {
        O2SCL_ERR((((std::string)"Value '")+dtos(val)+
                   "' not finite for column '"+
                   scol+"' in table::lookup_val()").c_str(),exc_einval);
        return exc_einval;
      }
      aciter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+" not found in "+
                   "table::lookup().").c_str(),exc_enotfound);
        return exc_enotfound;
      }
      return get(scol2,this->lookup(scol,val));
    }

    /** \brief Exhaustively search column \c col for the value \c val
        \f$ {\cal O}(R \log(C)) \f$
    */
    size_t lookup(int icol, fp_t val) const {
      return lookup(get_column_name(icol),val);
    }      

    /** \brief Exhaustively search column \c col for many occurences 
        of \c val \f$ {\cal O}(R \log(C)) \f$
    */
    size_t mlookup(std::string scol, fp_t val, std::vector<size_t> &results,
                   fp_t threshold=0.0) const {
      size_t i;
      if (!std::isfinite(val)) {
        O2SCL_ERR((((std::string)"Value '")+dtos(val)+
                   "' not finite for column '"+
                   scol+"' in table::mlookup()").c_str(),exc_einval);
        return exc_einval;
      }
      aciter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+" not found in "+
                   "table::mlookup().").c_str(),exc_enotfound);
        return exc_enotfound;
      }
      if (threshold==0.0) {
        for(i=0;i<nlines;i++) {
          if (it->second.dat[i]==val) {
            results.push_back(i);
          }
        }
      } else {
        for(i=0;i<nlines;i++) {
          if (fabs(it->second.dat[i]-val)<threshold) {
            results.push_back(i);
          }
        }
      }
      return results.size();
    }
    //@}

    // --------------------------------------------------------
    /** \name Interpolation, differentiation, integration, max, min */
    //@{

    /// Set the base interpolation objects
    void set_interp_type(size_t interp_type) {
      itype=interp_type;
      if (intp_set) {
        delete si;
        intp_set=false;
      }
      return;
    }

    /** \brief Get the interpolation type
     */
    size_t get_interp_type() const {
      return itype;
    }

    /** \brief Set the interpolation object for strings \c sx and 
        \c sy to the specified pointer
    */
    int set_interp_obj(std::string sx, std::string sy,
                       interp_base<vec_t,vec_t,fp_t> *si2) {
      if (intp_set==true) {
        delete si;
      } else {
        intp_set=true;
      }
      intp_colx=sx;
      intp_coly=sy;
      si=si2;
      return 0;
    }

    /** \brief Clear the interpolation object
     */
    void clear_interp_obj() {
      if (intp_set) {
        delete si;
        intp_set=false;
      }
      return;
    }
    
    /** \brief Interpolate value \c x0 from column named \c sx 
        into column named \c sy

        This function is \f$ {\cal O}(\log(R) \log(C)) \f$
        but can be as bad as \f$ {\cal O}(C \log(R) \f$ if the
        relevant columns are not well ordered.
    */
    fp_t interp(std::string sx, fp_t x0, std::string sy) {
      fp_t ret;
      aiter itx=atree.find(sx), ity=atree.find(sy);
      if (itx==atree.end() || ity==atree.end()) {
        O2SCL_ERR((((std::string)"Columns '")+sx+"' or '"+sy+
                   "' not found in table::interp().").c_str(),
                  exc_enotfound);
        return 0.0;
      }
      if (!std::isfinite(x0)) {
        O2SCL_ERR("x0 not finite in table::interp().",exc_einval);
        return exc_einval;
      }
      if (intp_set==false || sx!=intp_colx || sy!=intp_coly) {
        if (intp_set==true) {
          delete si;
        } else {
          intp_set=true;
        }
        si=new interp_vec<vec_t>(nlines,itx->second.dat,
                                 ity->second.dat,itype);
        intp_colx=sx;
        intp_coly=sy;
      }
      ret=si->eval(x0);
      return ret;
    }

    /** \brief Interpolate value \c x0 from column named \c sx 
        into column named \c sy (const version)

        This function is \f$ {\cal O}(\log(R) \log(C)) \f$
        but can be as bad as \f$ {\cal O}(C \log(R) \f$ if the
        relevant columns are not well ordered.
    */
    fp_t interp_const(std::string sx, fp_t x0, std::string sy) const {
      fp_t ret;
      aciter itx=atree.find(sx), ity=atree.find(sy);
      if (itx==atree.end() || ity==atree.end()) {
        O2SCL_ERR((((std::string)"Columns '")+sx+"' or '"+sy+
                   "' not found in table::interp_const().").c_str(),
                  exc_enotfound);
        return 0.0;
      }
      if (!std::isfinite(x0)) {
        O2SCL_ERR("x0 not finite in table::interp_const().",exc_einval);
        return exc_einval;
      }
      interp_vec<vec_t> sic(nlines,itx->second.dat,ity->second.dat,itype);

      ret=sic.eval(x0);
      return ret;
    }
    
    /** \brief Interpolate value \c x0 from column with index
        \c ix into column with index \c iy 
        \f$ {\cal O}(\log(R)) \f$
    */
    fp_t interp(size_t ix, fp_t x0, size_t iy) {
      return interp(get_column_name(ix),x0,get_column_name(iy));
    }

    /** \brief Interpolate value \c x0 from column with index \c ix 
        into column with index \c iy \f$ {\cal O}(\log(R)) \f$
    */
    fp_t interp_const(size_t ix, fp_t x0, size_t iy) const {
      return interp_const(get_column_name(ix),x0,get_column_name(iy));
    }

    /** \brief Make a new column named \c yp which is the 
        derivative \f$ y^{\prime}(x) \f$ formed from columns
        named \c x and \c y \f$ {\cal O}(R \log(C)) \f$

        If the column \c yp is not already in the table it is
        automatically created.
    */
    void deriv(std::string x, std::string y, std::string yp) {

      aiter itx, ity, ityp;

      itx=atree.find(x);
      ity=atree.find(y);
      ityp=atree.find(yp);

      if (ityp==atree.end()) {
        new_column(yp);
        ityp=atree.find(yp);
      }
    
      if (itx==atree.end() || ity==atree.end() || ityp==atree.end()) {
        O2SCL_ERR("Column not found in table::deriv(string,string,string).",
                  exc_enotfound);
        return;
      }
  
      size_t ix=lookup_column(x);
      size_t iy=lookup_column(y);
      for(int i=0;i<((int)nlines);i++) {
        ityp->second.dat[i]=deriv(ix,(itx->second.dat)[i],iy);
      }
  
      return;
    }

    /** \brief Compute the first derivative of the function defined
        by x-values stored in column named \c sx and y-values stored
        in column named \c sy at the value \c x0

        This function is O(log(C)*log(R)) but can be as bad as
        O(log(C)*R) if the relevant columns are not well ordered.
    */
    fp_t deriv(std::string sx, fp_t x0, std::string sy) {
      fp_t ret;
      aiter itx=atree.find(sx), ity=atree.find(sy);
      if (itx==atree.end() || ity==atree.end()) {
        O2SCL_ERR((((std::string)"Columns '")+sx+"' or '"+sy+
                   "' not found in table::deriv("+
                   "string,fp_t,string).").c_str(),
                  exc_enotfound);
        return 0.0;
      }
      if (!std::isfinite(x0)) {
        O2SCL_ERR("x0 not finite in table::deriv(string,fp_t,string).",
                  exc_einval);
        return exc_einval;
      }
      if (intp_set==false || sx!=intp_colx || sy!=intp_coly) {
        if (intp_set==true) {
          delete si;
        } else {
          intp_set=true;
        }
        si=new interp_vec<vec_t>
          (nlines,itx->second.dat,ity->second.dat,itype);
			 
        intp_colx=sx;
        intp_coly=sy;
      }
      ret=si->deriv(x0);
      return ret;
    }

    /** \brief Compute the first derivative of the function defined
        by x-values stored in column named \c sx and y-values stored
        in column named \c sy at the value \c x0 (const version)

        O(log(C)*log(R)) but can be as bad as O(log(C)*R) if 
        the relevant columns are not well ordered.
    */
    fp_t deriv_const(std::string sx, fp_t x0, std::string sy) const {
      fp_t ret;
      aciter itx=atree.find(sx), ity=atree.find(sy);
      if (itx==atree.end() || ity==atree.end()) {
        O2SCL_ERR((((std::string)"Columns '")+sx+"' or '"+sy+
                   "' not found in table::deriv_const().").c_str(),
                  exc_enotfound);
        return 0.0;
      }
      if (!std::isfinite(x0)) {
        O2SCL_ERR("x0 not finite in table::deriv_const().",exc_einval);
        return exc_einval;
      }
      interp_vec<vec_t> sic
        (nlines,itx->second.dat,ity->second.dat,itype);
      ret=sic.deriv(x0);
      return ret;
    }
  
    /** \brief Compute the first derivative of the function defined
        by x-values stored in column with index \c ix and y-values stored
        in column with index \c iy at the value \c x0
      
        O(log(R)) but can be as bad as O(R) if the relevant columns
        are not well ordered.
    */
    fp_t deriv(size_t ix, fp_t x0, size_t iy) {
      return deriv(get_column_name(ix),x0,get_column_name(iy));
    }

    /** \brief Compute the first derivative of the function defined
        by x-values stored in column with index \c ix and y-values stored
        in column with index \c iy at the value \c x0 (const version)
      
        O(log(R)) but can be as bad as O(R) if 
        the relevant columns are not well ordered.
    */
    fp_t deriv_const(size_t ix, fp_t x0, size_t iy) const {
      return deriv_const(get_column_name(ix),x0,get_column_name(iy));
    }

    /** \brief Create a new column named \c yp which is 
        equal to the second derivative of the function defined by 
        x-values stored in column named \c x and y-values 
        stored in column named \c y, i.e.
        \f$ y^{\prime \prime}(x) \f$ - O(log(C)*R).

        If the column \c yp is not already in the table it is
        automatically created.
    */
    void deriv2(std::string x, std::string y, std::string yp) {
      aiter itx, ity, ityp;

      itx=atree.find(x);
      ity=atree.find(y);
      ityp=atree.find(yp);
    
      if (ityp==atree.end()) {
        new_column(yp);
        ityp=atree.find(yp);
      }

      if (itx==atree.end() || ity==atree.end() || ityp==atree.end()) {
        O2SCL_ERR("Column not found in table::deriv2(string,string,string).",
                  exc_enotfound);
        return;
      }
  
      size_t ix=lookup_column(x);
      size_t iy=lookup_column(y);
      for(int i=0;i<((int)nlines);i++) {
        ityp->second.dat[i]=deriv2(ix,itx->second.dat[i],iy);
      }
  
      return;
    }

    /** \brief Compute the second derivative of the function defined
        by x-values stored in column named \c sx and y-values stored
        in column named \c sy at the value \c x0

        O(log(C)*log(R)) but can be as bad as O(log(C)*R) if 
        the relevant columns are not well ordered.
    */
    fp_t deriv2(std::string sx, fp_t x0, std::string sy) {
      fp_t ret;
      aiter itx=atree.find(sx), ity=atree.find(sy);
      if (itx==atree.end() || ity==atree.end()) {
        O2SCL_ERR((((std::string)"Columns '")+sx+"' or '"+sy+
                   "' not found in table::deriv2("+
                   "string,fp_t,string).").c_str(),
                  exc_enotfound);
        return 0.0;
      }
      if (!std::isfinite(x0)) {
        O2SCL_ERR("x0 not finite in table::deriv2(string,fp_t,string).",
                  exc_einval);
        return exc_einval;
      }
      if (intp_set==false || sx!=intp_colx || sy!=intp_coly) {
        if (intp_set==true) {
          delete si;
        } else {
          intp_set=true;
        }
        si=new interp_vec<vec_t>
          (nlines,itx->second.dat,ity->second.dat,itype);
			 
        intp_colx=sx;
        intp_coly=sy;
      }
      ret=si->deriv2(x0);
      return ret;
    }

    /** \brief The Compute the second derivative of the function defined
        by x-values stored in column named \c sx and y-values stored
        in column named \c sy at the value \c x0 (const version)

        O(log(C)*log(R)) but can be as bad as O(log(C)*R) if 
        the relevant columns are not well ordered.
    */
    fp_t deriv2_const(std::string sx, fp_t x0, std::string sy) const {
      fp_t ret;
      aciter itx=atree.find(sx), ity=atree.find(sy);
      if (itx==atree.end() || ity==atree.end()) {
        O2SCL_ERR((((std::string)"Columns '")+sx+"' or '"+sy+
                   "' not found in table::deriv2_const().").c_str(),
                  exc_enotfound);
        return 0.0;
      }
      if (!std::isfinite(x0)) {
        O2SCL_ERR("x0 not finite in table::deriv2_const().",exc_einval);
        return exc_einval;
      }
      interp_vec<vec_t> sic
        (nlines,itx->second.dat,ity->second.dat,itype);
      ret=sic.deriv2(x0);
      return ret;
    }

    /** \brief Compute the second derivative of the function defined
        by x-values stored in column with index \c ix and y-values stored
        in column with index \c iy at the value \c x0

        O(log(R)) but can be as bad as O(R) if 
        the relevant columns are not well ordered.
    */
    fp_t deriv2(size_t ix, fp_t x0, size_t iy) {
      return deriv2(get_column_name(ix),x0,get_column_name(iy));
    }

    /** \brief Compute the second derivative of the function defined
        by x-values stored in column with index \c ix and y-values stored
        in column with index \c iy at the value \c x0 (const version)

        O(log(R)) but can be as bad as O(R) if 
        the relevant columns are not well ordered.
    */
    fp_t deriv2_const(size_t ix, fp_t x0, size_t iy) const {
      return deriv2_const(get_column_name(ix),x0,get_column_name(iy));
    }

    /** \brief Compute the integral of the function defined
        by x-values stored in column named \c sx and y-values stored
        in column named \c sy between the values \c x1 and \c x2
      
        O(log(C)*log(R)) but can be as bad as O(log(C)*R) if 
        the relevant columns are not well ordered.
    */
    fp_t integ(std::string sx, fp_t x1, fp_t x2, std::string sy) {
      fp_t ret;
      aiter itx=atree.find(sx), ity=atree.find(sy);
      if (itx==atree.end() || ity==atree.end()) {
        O2SCL_ERR((((std::string)"Columns '")+sx+"' or '"+sy+
                   "' not found in table::integ"+
                   "(string,fp_t,fp_t,string).").c_str(),
                  exc_enotfound);
        return 0.0;
      }
      if (!std::isfinite(x1) || !std::isfinite(x2)) {
        std::string msg=((std::string)"Value x1=")+dtos(x1)+" or x2="+
          dtos(x2)+" not finite in table.integ(string,fp_t,fp_t,string).";
        O2SCL_ERR(msg.c_str(),exc_einval);
      }
      if (intp_set==false || sx!=intp_colx || sy!=intp_coly) {
        if (intp_set==true) {
          delete si;
        } else {
          intp_set=true;
        }
        si=new interp_vec<vec_t>
          (nlines,itx->second.dat,ity->second.dat,itype);
			 
        intp_colx=sx;
        intp_coly=sy;
      }
      ret=si->integ(x1,x2);
      return ret;
    }
  
    /** \brief Compute the integral of the function defined
        by x-values stored in column named \c sx and y-values stored
        in column named \c sy between the values \c x1 and \c x2 
        (const version)
	
        O(log(C)*log(R)) but can be as bad as O(log(C)*R) if 
        the relevant columns are not well ordered.
    */
    fp_t integ_const(std::string sx, fp_t x1, fp_t x2, 
                       std::string sy) const {
      fp_t ret;
      aciter itx=atree.find(sx), ity=atree.find(sy);
      if (itx==atree.end() || ity==atree.end()) {
        O2SCL_ERR((((std::string)"Columns '")+sx+"' or '"+sy+
                   "' not found in table::integ_const().").c_str(),
                  exc_enotfound);
        return 0.0;
      }
      if (!std::isfinite(x1) || !std::isfinite(x2)) {
        O2SCL_ERR("x1 or x2 not finite in table::integ_const().",exc_einval);
        return exc_einval;
      }
      interp_vec<vec_t> sic
        (nlines,itx->second.dat,ity->second.dat,itype);
      ret=sic.integ(x1,x2);
      return ret;
    }
  
    /** \brief Compute the integral of the function defined
        by x-values stored in column with index \c ix and y-values stored
        in column with index \c iy between the values \c x1 and \c x2 

        O(log(R)) but can be as bad as O(R) if 
        the relevant columns are not well ordered.
    */
    fp_t integ(size_t ix, fp_t x1, fp_t x2, size_t iy) {
      return integ(get_column_name(ix),x1,x2,
                   get_column_name(iy));
    }

    /** \brief Compute the integral of the function defined
        by x-values stored in column with index \c ix and y-values stored
        in column with index \c iy between the values \c x1 and \c x2 
        (const version)

        O(log(R)) but can be as bad as O(R) if 
        the relevant columns are not well ordered.
    */
    fp_t integ_const(size_t ix, fp_t x1, fp_t x2, size_t iy) const {
      return integ_const(get_column_name(ix),x1,x2,
                         get_column_name(iy));
    }

    /** \brief Create a new column named \c ynew which is 
        equal to the integral of the function defined by 
        x-values stored in column named \c x and y-values 
        stored in column named \c y

        This function is O(log(R)) but can be as bad as O(R) if the
        relevant columns are not well ordered.

        If the column \c ynew is not already in the table it is
        automatically created.
    */
    void integ(std::string x, std::string y, std::string ynew) {
      aiter itx, ity, itynew;

      itx=atree.find(x);
      ity=atree.find(y);
      itynew=atree.find(ynew);
    
      if (itynew==atree.end()) {
        new_column(ynew);
        itynew=atree.find(ynew);
      }

      if (itx==atree.end() || ity==atree.end() || itynew==atree.end()) {
        O2SCL_ERR("Column not found in table::integ(string,string,string).",
                  exc_enotfound);
        return;
      }
  
      size_t ix=lookup_column(x);
      size_t iy=lookup_column(y);
      for(size_t i=0;i<nlines;i++) {
        itynew->second.dat[i]=integ(ix,(itx->second.dat)[0],
                                    (itx->second.dat)[i],iy);
      }
  
      return;
    }

    /** \brief Return column maximum. Makes no assumptions about 
        ordering, \f$ {\cal O}(R) \f$
    */
    fp_t max(std::string scol) const {
      fp_t ret=0.0;
      int i;
      if (is_column(scol)==false) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   "' not found in table::max().").c_str(),exc_enotfound);
        return 0.0;
      }
      const vec_t &dcol=get_column(scol);
      bool setb=false;
      for(i=0;i<((int)nlines);i++) {
        if (std::isfinite(dcol[i])) {
          if (setb==false) {
            ret=dcol[i];
            setb=true;
          } else if (dcol[i]>ret) {
            ret=dcol[i];
          }
        }
      }
      if (setb==false) {
        O2SCL_ERR((((std::string)"No finite values in column '")+scol+
                   "' in table::max().").c_str(),exc_efailed);
        return 0.0;
      }
      return ret;
    }

    /** \brief Return column minimum. Makes no assumptions about 
        ordering, \f$ {\cal O}(R) \f$ 
    */
    fp_t min(std::string scol) const {
      fp_t ret=0.0;
      int i;
      if (is_column(scol)==false) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   "' not found in table::min().").c_str(),exc_enotfound);
        return 0.0;
      }
      const vec_t &dcol=get_column(scol);
      bool setb=false;
      for(i=0;i<((int)nlines);i++) {
        if (std::isfinite(dcol[i])) {
          if (setb==false) {
            ret=dcol[i];
            setb=true;
          } else if (dcol[i]<ret) {
            ret=dcol[i];
          }
        }
      }
      if (setb==false) {
        O2SCL_ERR((((std::string)"No finite values in column '")+scol+
                   "' in table::min().").c_str(),exc_efailed);
        return 0.0;
      }
      return ret;
    }
    //@}

    // --------------------------------------------------------
    /** \name Subtable method */
    //@{

    /** \brief Make a subtable

        Uses the columns specified in \c list from the row \c top
        to the row of index \c bottom to generate a new table
        which is a copy of part of the original.
    */
    void subtable(std::string list, size_t top, 
                  size_t bottom, table<vec_t,fp_t> &tnew) const {

      tnew.clear_all();
      int sublines, i;
      std::string head;
      aciter it;
  
      if (top>bottom) {
        size_t tmp=bottom;
        bottom=top;
        top=tmp;
      }
      sublines=bottom-top+1;
      if (nlines==0) {
        O2SCL_ERR2("Can't make a subtable of an empty table. ",
                   "Returning 0 in table::subtable().",
                   exc_einval);
        return;
      }
      if (bottom+1>nlines) {
        O2SCL_ERR2("Requested row beyond nlines. Adjusting ",
                   "and continuing in table::subtable().",exc_einval);
        bottom=nlines-1;
      }
  
      std::istringstream is(list);

      tnew.set_nlines(sublines);
      while(is >> head) {
        it=atree.find(head);
        if (it==atree.end()) {
          O2SCL_ERR
            ((((std::string)"Couldn't find column named ")+head+
              " in table::subtable(). Returning 0.").c_str(),
             exc_einval);
        }
        tnew.new_column(head);
        vec_t &dcol=tnew.get_column(head);
        for(i=0;i<sublines;i++) {
          dcol[i]=it->second.dat[i+top];
        }
      } 
      if (tnew.get_ncolumns()==0) {
        O2SCL_ERR("Subtable has no columns in table::subtable().",
                  exc_einval);
      }
      //}
      tnew.nlines=sublines;
  
      return;
    }
    //@}
  
    // --------------------------------------------------------
    /** \name Clear methods */
    //@{

    /** \brief Zero the data entries but keep the column names 
        and nlines fixed
    */
    void zero_table() {
      aiter it;
      for(it=atree.begin();it!=atree.end();it++) {
        for(int j=0;j<((int)nlines);j++) {
          it->second.dat[j]=0.0;
        }
      }

      if (intp_set) {
        intp_set=false;
        delete si;
      }

      return;
    }

    /** \brief Clear everything
     */
    virtual void clear() {
      clear_table();
      clear_constants();
      return;
    }

    /** \brief Clear the table and the column names (but leave constants)
     */
    virtual void clear_table() {
      atree.clear();
      alist.clear();
      nlines=0;
      if (intp_set==true) {
        delete si;
        intp_set=false;
      }
      return;
    }

    /** \brief Remove all of the data by setting the number
        of lines to zero

        This leaves the column names intact and does not remove
        the constants.
    */
    void clear_data() {
      nlines=0;   
      if (intp_set==true) {
        delete si; 
        intp_set=false;
      }
      return;
    }

    /// CLear all constants
    void clear_constants() {
      constants.clear();
      return;
    }
    //@}

    // --------------------------------------------------------
    /** \name Sorting methods */
    //@{

    /** \brief Sort the entire table by the column \c scol

        \note This function works by allocating space for an entirely
        new chunk of memory for the data in the table.
    */
    void sort_table(std::string scol) {

      size_t ncols=get_ncolumns(), nlins=get_nlines();

      // Make a copy of the table
      boost::numeric::ublas::matrix<fp_t> data_copy(ncols,nlins);
      for(size_t i=0;i<ncols;i++) {
        for(size_t j=0;j<nlins;j++) {
          data_copy(i,j)=get(i,j);
        }
      }

      permutation order(nlins);
      aiter it=atree.find(scol);
      vec_t &data=it->second.dat;
      vector_sort_index(nlins,data,order);
      for(size_t i=0;i<ncols;i++) {
        for(size_t j=0;j<nlins;j++) {
          set(i,j,data_copy(i,order[j]));
        }
      }
  
      if (intp_set) {
        intp_set=false;
        delete si;
      }

      return;
    }

    /** \brief Individually sort the column \c scol
     */
    void sort_column(std::string scol) {
      int i;
      aiter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   " not found in table::sort_column().").c_str(),
                  exc_enotfound);
        return;
      }

      vector_sort_double(nlines,it->second.dat);

      if (intp_set && (scol==intp_colx || scol==intp_coly)) {
        intp_set=false;
        delete si;
      }

      return;
    }
    //@}

    // --------------------------------------------------------
    /** \name Summary method */
    //@{
    /** \brief Output a summary of the information stored

        Outputs the number of constants, the number of columns,
        a list of the column names, and the number of lines of
        data.
    */
    virtual void summary(std::ostream *out=0, size_t ncol=79) const {

      if (out==0) {
        out=&std::cout;
      }
      if (constants.size()==1) {
        (*out) << "1 constant:" << std::endl;
      } else {
        (*out) << constants.size() << " constants:" << std::endl;
      }
      typename std::map<std::string,fp_t>::const_iterator mit;
      for(mit=constants.begin();mit!=constants.end();mit++) {
        (*out) << mit->first << " " << mit->second << std::endl;
      }

      // Output number of columns and preprend column numbers
      size_t nh=get_ncolumns(), nh2;

      if (nh==0) {

        (*out) << "No columns." << std::endl;

      } else {
    
        if (nh==1) {
          (*out) << "1 column: " << std::endl;
        } else {
          (*out) << nh << " columns: " << std::endl;
        }
        std::vector<std::string> h(nh);
        for(size_t i=0;i<nh;i++) {
          h[i]=szttos(i)+". "+get_column_name(i);
        }
    
        std::vector<std::string> h2;
        // Convert to string with width 'ncol'
        screenify(nh,h,h2,ncol);
        nh2=h2.size();
    
        // Output column names
        for(size_t i=0;i<nh2;i++) {
          (*out) << h2[i] << std::endl;
        }

      }
  
      if (get_nlines()==0) (*out) << "No lines of data." << std::endl;
      else if (get_nlines()==1) (*out) << "One line of data." << std::endl;
      (*out) << get_nlines() << " lines of data." << std::endl;
  
      return;
    }
    //@}

    /// \name Constant manipulation
    //@{
    /** \brief Add a constant, or if the constant already exists, change 
        its value
    */
    virtual void add_constant(std::string name, fp_t val) {
      if (constants.find(name)!=constants.end()) {
        constants.find(name)->second=val;
        return;
      }
      constants.insert(make_pair(name,val));
      return;
    }

    /** \brief Set a constant equal to a value, but don't add it if
        not already present

        If \c err_on_notfound is <tt>true</tt> (the default), then
        this function throws an exception if a constant with
        name \c name is not found. If \c err_on_notfound is
        <tt>false</tt>, then if a constant with name \c name
        is not found this function just silently returns
        \ref o2scl::exc_enotfound.
    */
    virtual int set_constant(std::string name, fp_t val,
                             bool err_on_notfound=true) {
      if (constants.find(name)!=constants.end()) {
        constants.find(name)->second=val;
        return 0;
      }
      if (err_on_notfound) {
        std::string err=((std::string)"No constant with name '")+name+
          "' in table::set_constant().";
        O2SCL_ERR(err.c_str(),exc_enotfound);
      }
      return exc_enotfound;
    }

    /// Test if \c name is a constant
    virtual bool is_constant(std::string name) const {
      if (constants.find(name)==constants.end()) {
        return false;
      }
      return true;
    }

    /// Get a constant
    virtual fp_t get_constant(std::string name) const {
      if (constants.find(name)==constants.end()) {
        std::string err=((std::string)"No constant with name '")+name+
          "' in table::get_constant(string).";
        O2SCL_ERR(err.c_str(),exc_einval);
      }
      return constants.find(name)->second;
    }

    /// Get the number of constants
    virtual size_t get_nconsts() const {
      return constants.size();
    }

    /// Get a constant by index
    virtual void get_constant(size_t ix, std::string &name, fp_t &val) const {
      if (ix<constants.size()) {
        typename std::map<std::string,fp_t>::const_iterator cit=
          constants.begin();
        for(size_t i=0;i<ix;i++) cit++;
        name=cit->first;
        val=cit->second;
        return;
      }
      O2SCL_ERR("Index too large in table::get_constant(size_t,string,fp_t).",
                exc_eindex);
      return;
    }

    /// Remove a constant
    virtual void remove_constant(std::string name) {
      if (constants.find(name)==constants.end()) {
        O2SCL_ERR2("Constant not present in ",
                   "table::remove_constant().",o2scl::exc_einval);
      }
      constants.erase(name);
      return;
    }
    //@}
    
    /// \name Miscellaneous methods
    //@{
    /** \brief Compute the rolling average of column named \c col
     */
    virtual void average_col_roll(std::string col_name, size_t window) {
      if (window<2) {
        O2SCL_ERR("Window less than 2 in table::average_rows().",
                  o2scl::exc_einval);
      }

      size_t nl=get_nlines();
      size_t k=lookup_column(col_name);

      // A non-const reference to the data vector
      vec_t &v=alist[k]->second.dat;
    
      o2scl::vector_roll_avg<vec_t,fp_t>(nl,v,window);
    
      return;
    }
  
    /** \brief Average nearby rows together over the entire table
     */
    virtual void average_rows(size_t window, bool rolling=false) {
    
      if (window<2) {
        O2SCL_ERR("Window less than 2 in table::average_rows().",
                  o2scl::exc_einval);
      }
      size_t nl=get_nlines();
      size_t nc=get_ncolumns();

      int verbose=0;
    
      if (rolling==false) {
      
        // First pass, combine rows together into the first row
        // of every group
        for(size_t i=0;i<nl;i+=window) {
          for(size_t j=1;j<window;j++) {
            for(size_t k=0;k<nc;k++) {
              if (i+j<nl) {
                set(k,i,get(k,i)+get(k,i+j));
                if (k==0 && verbose>0) {
                  std::cout << "Adding row " << i+j << " to row "
                            << i << std::endl;
                }
              }
            }
          }
        }
        // Variable nw is the number of windows
        size_t nw=nl/window;
        // Correctly handle the last window
        if (nw*window<nl) nw++;
        if (verbose>0) {
          std::cout << "nl: " << nl << " window: " << window << " nw: "
                    << nw << std::endl;
        }

        // Move the first row of every group to the top of the table
        for(int i=((int)nw)-1;i>=0;i--) {
          // The size of the window may be different for the last window
          // so handle that accordingly
          int win=window;
          if (i==((int)nw)-1) {
            win=((int)nl)-(((int)nw)-1)*window;
          } else {
            win=window;
          }
          if (verbose>0) {
            std::cout << "i: " << i << " win: " << win << " i*window: "
                      << i*window << std::endl;
          }
          for(size_t k=0;k<nc;k++) {
            set(k,i,get(k,i*window)/((fp_t)win));
          }
        }

        // Finally, change the number of lines
        set_nlines(nw);
      
      } else {

        // Proceed one column at a time
        for(size_t k=0;k<nc;k++) {

          // A non-const reference to the data vector
          vec_t &v=alist[k]->second.dat;
	
          o2scl::vector_roll_avg<vec_t,fp_t>(nl,v,window);
        }
      
      }
      return;
    }
  
    /// Clear the current table and read from a generic data file
    virtual int read_generic(std::istream &fin, int verbose=0) {

      clear();
      
      fp_t data;
      std::string line;
      std::string cname;

      // Read the first line
      getline(fin,line);
      
      // Determine if there are constants
      std::vector<std::string> vsc;
      split_string_delim(line,vsc,' ');
      if (vsc.size()>1 &&
          (vsc[1]=="constants." || vsc[1]=="constant.")) {
        size_t n_const=o2scl::stoszt(vsc[0]);
        std::string name;
        fp_t val;
        for(size_t i=0;i<n_const;i++) {
          fin >> name >> val;
          set_constant(name,val);
        }
        // Read the remaining carriage return at the end of the
        // constant list
        getline(fin,line);
        // Read the next full line
        getline(fin,line);
      }

      // Determine if the interpolation type was specified
      std::vector<std::string> vsi;
      split_string_delim(line,vsi,' ');
      if (vsi.size()>1 && vsi[0]=="Interpolation:") {
        set_interp_type(o2scl::stoszt(vsi[1]));
        // Read the next full line
        getline(fin,line);
      }
      
      // See if the next line has column names
      std::vector<std::string> onames, nnames;
      std::istringstream is(line);
      while (is >> cname) {
        onames.push_back(cname);
        if (verbose>2) {
          std::cout << "Read possible column name: " << cname << std::endl;
        }
      }

      // Count number of likely numbers in the first row
      size_t n_nums=0;
      for(size_t i=0;i<onames.size();i++) {
        if (is_number(onames[i])) n_nums++;
      }

      int irow=0;

      if (n_nums==onames.size()) {

        if (verbose>0) {
          std::cout << "First row looks like it contains numerical values." 
                    << std::endl;
          std::cout << "Creating generic column names: ";
        }

        for(size_t i=0;i<onames.size();i++) {
          nnames.push_back(((std::string)"c")+szttos(i+1));
          if (verbose>0) std::cout << nnames[i] << " ";
      
        }
        if (verbose>0) std::cout << std::endl;

        // Make columns
        for(size_t i=0;i<nnames.size();i++) {
          new_column(nnames[i]);
        }

        // Add first row of data
        set_nlines_auto(irow+1);
        for(size_t i=0;i<onames.size();i++) {
          set(i,irow,o2scl::stod(onames[i]));
        }
        irow++;

      } else {

        // Ensure good column names
        for(size_t i=0;i<onames.size();i++) {
          std::string temps=onames[i];
          make_fp_varname(temps);
          make_unique_name(temps,nnames);
          nnames.push_back(temps);
          if (temps!=onames[i] && verbose>0) {
            std::cout << "Converted column named '" << onames[i] << "' to '" 
                      << temps << "'." << std::endl;
          }
        }

        // Make columns
        for(size_t i=0;i<nnames.size();i++) {
          new_column(nnames[i]);
        }

      }

      // Read remaining rows
      while ((fin) >> data) {
        set_nlines_auto(irow+1);
        set(0,irow,data);
        for(size_t i=1;i<get_ncolumns();i++) {
          (fin) >> data;
          set(i,irow,data);
        }
        irow++;
      }

      if (intp_set) {
        intp_set=false;
        delete si;
      }

      return 0;
    }

    /** \brief Check if the table object appears to be valid
     */
    void is_valid() const {
      if (maxlines<nlines) {
        O2SCL_ERR2("Value of maxlines smaller than nlines ",
                   "in table::is_valid().",exc_esanity);
      }
      if (atree.size()!=alist.size()) {
        O2SCL_ERR2("Size of table and list do not match in ",
                   "table::is_valid().",exc_esanity);
        return;
      }
      std::vector<int> index_check;
      for(aciter it=atree.begin();it!=atree.end();it++) {
        if (it->second.dat.size()!=maxlines) {
          O2SCL_ERR2("Vector with size different than maxlines ",
                     "in table::is_valid().",exc_esanity);
        }
        if (it->second.index!=alist[it->second.index]->second.index) {
          O2SCL_ERR((((std::string)"Problem with iterator for entry '")+
                     it->first+"' in list in table::is_valid().").c_str(),
                    exc_esanity);
        }
        index_check.push_back(it->second.index);
      }
      vector_sort<std::vector<int>,int>(atree.size(),index_check);
      for(int i=0;i<((int)atree.size());i++) {
        if (i!=index_check[i]) {
          O2SCL_ERR("Column indexes wrong.",exc_esanity);
        }
      }
      for(int i=0;i<((int)atree.size());i++) {
        if (alist[i]->second.index!=i) {
          O2SCL_ERR((((std::string)"Problem with index of entry ")+
                     itos(i)+" in list in table::is_valid().").c_str(),
                    exc_esanity);
          return;
        }
      }
      return;
    }

    /// Return the type, \c "table".
    virtual const char *type() { return "table"; }
    //@}

    /** \name Parsing mathematical functions specified as strings
     */
    //@{
    /** \brief Create new columns or recompute from a list of functions
	
        The list should be a space-delimited list of entries of the
        form <tt>name=function</tt> where <tt>name</tt> is the
        column name and <tt>function</tt> the function specifing the
        values for the column. If a column named <tt>name</tt> is 
        already present, it is overwritten. Otherwise, a new column
        is created. 

        \comment
        The formulas in \c list may depend on any of the column names
        that will be defined later in \c list. For example, for a
        table initially containing two columns, \c x and \c y, the
        calls
        \code
        function_columns("a=2*z z=x+y");
        \endcode
        \code
        function_columns("z=x+y a=2*z");
        \endcode
        both work.
        Circular dependencies do not work, for example
        \code
        function_columns("a=2*z z=a*3");
        \endcode
        will cause the error handler to be thrown. 
        \endcomment
    */
    void functions_columns(std::string list) {

      // Separate the list into names and functions
      std::vector<std::string> funcs, names;
      {
        std::string stemp;
        std::istringstream is(list);
        while(is >> stemp) funcs.push_back(stemp);
        for(size_t i=0;i<(funcs.size());i++) {
          names.push_back(funcs[i].substr(0,funcs[i].find("=")));
          funcs[i]=funcs[i].substr(funcs[i].find("=")+1,
                                   funcs[i].length()-funcs[i].find("=")-1);
          if (names[i].length()==0 || funcs[i].length()==0) {
            O2SCL_ERR2("Name or function blank in ",
                       "table::functions_columns().",exc_einval);
          }
        }
      }

      typename std::map<std::string,fp_t> vars;
      typename std::map<std::string,fp_t>::const_iterator mit;
      for(mit=constants.begin();mit!=constants.end();mit++) {
        vars[mit->first]=mit->second;
      }
    
      std::vector<calc_utf8<> > calcs(funcs.size());
      std::vector<vec_t> newcols(funcs.size());
    
      for(size_t j=0;j<funcs.size();j++) {
        calcs[j].compile(funcs[j].c_str(),&vars);
        newcols[j].resize(maxlines);
      }
    
      // Calculate all of the columns in the newcols list:
      for(size_t i=0;i<nlines;i++) {
      
        // Record the values of the variables for this line:
        for(size_t j=0;j<atree.size();j++) {
          vars[get_column_name(j)]=(*this)[j][i];
        }
      
        // Evaluate the new columns
        for(size_t j=0;j<funcs.size();j++) {
          newcols[j][i]=calcs[j].eval(&vars);
        }
      }

      for(size_t j=0;j<funcs.size();j++) {
        if (!is_column(names[j])) {
          new_column(names[j]);
        }
        swap_column_data(names[j],newcols[j]);
      }
    
      return;
    }

    /** \brief Make a column from the function specified in
        <tt>function</tt> and add it to the table.
	
        If a column named \c scol already exists, the data already
        present is overwritten with the result. Otherwise, a new
        column is created and filled with the result.
    */
    void function_column(std::string function, std::string scol) {
      int ret, i, j;
      std::string vlist;
      aiter it;

      // Create new column if necessary
      if (!is_column(scol)) {
        new_column(scol);
      }

      // Find vector reference
      aiter it2=atree.find(scol);
      vec_t &colp=it2->second.dat;

      // Fill vector with result of function
      function_vector(function,colp);

      return;
    }

    /** \brief Compute a column from a function specified 
        in a string
      
        The type \c resize_vec_t must have <tt>resize()</tt> and
        <tt>size()</tt> methods. If \c vec does not have enough space to
        hold the number of entries given by \ref get_nlines(), it is
        resized.

        \verbatim embed:rst
        .. todo:: 

           In table::function_vector(): FIXME: there may be a problem
           with the OpenMP code if an exception is thrown in the
           calculator class and there is not a unique error handler
           for each thread.

        \endverbatim

        \comment
        This function must return an int rather than void because
        of the presence of the 'throw_on_err' mechanism
        \endcomment
    */
    template<class resize_vec_t>
    int function_vector(std::string function, resize_vec_t &vec,
                        bool throw_on_err=true) {

      int n_threads=1;
      int i_thread=0;
    
      // Resize vector if necessary (outside the parallel region)
      if (vec.size()<nlines) vec.resize(nlines);
      
#ifdef O2SCL_SET_OPENMP
#pragma omp parallel private(i_thread)
#endif
      {

#ifdef O2SCL_SET_OPENMP
        n_threads=omp_get_num_threads();
        i_thread=omp_get_thread_num();
#endif
        // Create a new random number generator for each thread,
        // and ensure a different seed for each thread
        rng<> r;
        unsigned int seed=time(0);
        r.set_seed(seed*(i_thread+1));

        // Parse function, separate calculator for each thread
        calc_utf8<> calc;
        calc.set_rng(r);

        std::map<std::string,fp_t> vars;
        
        if (false) {

          // Old version (slower)
          
          typename std::map<std::string,fp_t>::const_iterator mit;
          for(mit=constants.begin();mit!=constants.end();mit++) {
            vars[mit->first]=mit->second;
          }
          calc.compile(function.c_str(),&vars);

          // Create column from function
          for(int j=i_thread;j<((int)nlines);j+=n_threads) {
            for(aciter it=atree.begin();it!=atree.end();it++) {
              vars[it->first]=it->second.dat[j];
            }
            vec[j]=calc.eval(&vars);
          }
          
        } else {

          // New (hopefully faster) version which uses
          // calc_utf8::get_var_list() to obtain a list of variables
          // needed to compute the user-specified function.
          
          calc.compile(function.c_str(),0);

          // Get the variable list as a list of u32string
          std::vector<std::u32string> cols32=calc.get_var_list();

          // Convert it to a list of utf8 strings
          std::vector<std::string> cols(cols32.size());
          for(size_t ij=0;ij<cols32.size();ij++) {
            char32_to_utf8(cols32[ij],cols[ij]);
          }

          // At this point, the vector \c cols may contain the names
          // of constants which are not columns, so we have to use
          // is_column() below to double check.
          
          typename std::map<std::string,fp_t>::const_iterator mit;
          for(mit=constants.begin();mit!=constants.end();mit++) {
            vars[mit->first]=mit->second;
          }

          // Create column from function
          for(int j=i_thread;j<((int)nlines);j+=n_threads) {
            for(size_t k=0;k<cols.size();k++) {
              // Skip entries that are constants because they're
              // already taken care of above.
              if (this->is_column(cols[k])) {
                vars[cols[k]]=this->get(cols[k],j);
              }
            }
            vec[j]=calc.eval(&vars);
          }
          
        }
      
        // End of parallel region
      }

      return 0;
    }

    /** \brief Compute a value by applying a function to a row
     */
    fp_t row_function(std::string function, size_t row) const {

      // Parse function
      calc_utf8<> calc;
      typename std::map<std::string,fp_t> vars;
      typename std::map<std::string,fp_t>::const_iterator mit;
      for(mit=constants.begin();mit!=constants.end();mit++) {
        vars[mit->first]=mit->second;
      }
      calc.compile(function.c_str(),&vars);

      for(aciter it=atree.begin();it!=atree.end();it++) {
        vars[it->first]=it->second.dat[row];
      }

      fp_t dret=calc.eval(&vars);
      return dret;
    }

    /** \brief Find a row which maximizes a function
     */
    size_t function_find_row(std::string function) const {

      // Parse function
      calc_utf8<> calc;
      typename std::map<std::string,fp_t> vars;
      typename std::map<std::string,fp_t>::const_iterator mit;
      for(mit=constants.begin();mit!=constants.end();mit++) {
        vars[mit->first]=mit->second;
      }
      calc.compile(function.c_str(),&vars);

      fp_t best_val=0.0;
      size_t best_row=0;
      for(size_t row=0;row<nlines-1;row++) {
        for(aciter it=atree.begin();it!=atree.end();it++) {
          vars[it->first]=it->second.dat[row];
        }
        fp_t dtemp=calc.eval(&vars);
        if (row==0) {
          best_val=dtemp;
        } else {
          if (dtemp>best_val) {
            best_val=dtemp;
            best_row=row;
          }
        }
      }
  
      return best_row;
    }
    //@}

    // --------------------------------------------------------
    // Allow HDF5 I/O functions to access table data
  
#ifndef DOXYGEN  
    friend void o2scl_hdf::hdf_output
    (o2scl_hdf::hdf_file &hf, table<> &t, std::string name);
  
    template<class vecf_t, class fpf_t> friend void o2scl_hdf::hdf_input
    (o2scl_hdf::hdf_file &hf, table<vecf_t,fpf_t> &t, std::string name);
  
    friend void o2scl_hdf::hdf_output_data
    (o2scl_hdf::hdf_file &hf, table<> &t);
  
    template<class vecf_t, class fpf_t> friend void o2scl_hdf::hdf_input_data
    (o2scl_hdf::hdf_file &hf, table<vecf_t,fpf_t> &t);
  
    // --------------------------------------------------------
    // Allow matrix_view_table and matrix_view_table_transpose access
  
    template<typename vecf_t, typename fpf_t> friend class
    matrix_view_table;
    template<typename vecf_t, typename fpf_t> friend class
    matrix_view_table_transpose;
#endif

    // --------------------------------------------------------
  
#ifndef DOXYGEN_INTERNAL
  
  protected:
  
    /** \brief Returns a non-const reference to the column named \c col.
        \f$ {\cal O}(\log(C)) \f$
    */
    vec_t &get_column_no_const(std::string scol) {
      aiter it=atree.find(scol);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+scol+
                   "' not found in table::get_column() const.").c_str(),
                  exc_enotfound);
        return empty_col;
      }
      return it->second.dat;
    }
  
    /** \brief Set the elements of alist with the appropriate 
        iterators from atree. \f$ {\cal O}(C) \f$

        Generally, the end-user shouldn't need this method. It is 
        only used in delete_column() to rearrange the list when
        a column is deleted from the tree.
    */
    void reset_list() {
      aiter it;
      for(it=atree.begin();it!=atree.end();it++) {
        alist[it->second.index]=it;
      }
      return;
    }

    /** \brief Ensure a variable name does not match a function or contain 
        non-alphanumeric characters
    */
    void make_fp_varname(std::string &s) {
      if (s=="abs" || s=="acos" || s=="acosh" || s=="asin" ||
          s=="asinh" || s=="atan" || s=="atan2" || s=="atanh" ||
          s=="ceil" || s=="cos" || s=="cosh" || s=="cot" || s=="csc" ||
          s=="eval" || s=="exp" || s=="floor" || s=="if" || s=="int" ||
          s=="log" || s=="log10" || s=="max" || s=="min" || s=="sec" ||
          s=="sin" || s=="sinh" || s=="sqrt" || s=="tan" || s=="tanh") {
        s=((std::string)"v_")+s;
      } else if (s[0]>='0' && s[0]<='9') {
        s=((std::string)"v_")+s;
      }
  
      for(size_t i=0;i<s.length();i++) {
        if (!isalpha(s[i]) && !isdigit(s[i]) && s[i]!='_') s[i]='_';
      }
  
      return;
    }

    /// Make sure a name is unique
    void make_unique_name(std::string &colx, std::vector<std::string> &cnames) {
      bool done;

      do {
        done=true;
        for(size_t i=0;i<cnames.size();i++) {
          if (colx==cnames[i]) {
            done=false;
            i=cnames.size();
          }
        }
        if (done==false) {
          colx+='_';
        }
      } while (done==false);

      return;
    }

    /** \brief The list of constants 
     */
    std::map<std::string,fp_t> constants;

    /** \brief Column structure for \ref table [protected]

        This struct is used internally by \ref table to organize the
        columns and need not be instantiated by the casual end-user.
    */
    class col {
    
    public:
    
      /// Column of data
      vec_t dat;
      /// Column index
      int index;
    
      col() {
      }
    
      /** \brief Copy constructor 
       */
      col(const col &c) {
        dat=c.dat;
        index=c.index;
      }
      /** \brief Copy constructor for assignment operator
       */
      col &operator=(const col &c) {
        if (this!=&c) {
          dat=c.dat;
          index=c.index;
        }
        return *this;
      }
    
      /** \brief Swap method
       */
      friend void swap(col &t1, col &t2) {
        using std::swap;
        swap(t1.dat,t2.dat);
        swap(t1.index,t2.index);
        return;
      }
    };
  
    /// \name Iterator types
    //@{
    /// Map iterator type
    typedef typename std::map<std::string,col,
                              std::greater<std::string> >::iterator aiter;
    /// Const map iterator type
    typedef typename std::map<std::string,col,
                              std::greater<std::string> >::const_iterator 
    aciter;
    /// Vector iterator type
    typedef typename std::vector<aiter>::iterator aviter;
    //@}
  
    /// \name Actual data
    //@{
    /// The size of allocated memory
    size_t maxlines;
    /// The size of presently used memory
    size_t nlines;
    /// The tree of columns
    std::map<std::string,col,std::greater<std::string> > atree;
    /// The list of tree iterators
    std::vector<aiter> alist;
    //@}
  
    /// \name Column manipulation methods
    //@{
    /// Return the iterator for a column
    aiter get_iterator(std::string lname) {
      aiter it=atree.find(lname);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+lname+
                   " not found in table::get_iterator().").c_str(),
                  exc_enotfound);
      }
      return it;
    }
    /// Return the column structure for a column
    col *get_col_struct(std::string lname) {
      aiter it=atree.find(lname);
      if (it==atree.end()) {
        O2SCL_ERR((((std::string)"Column '")+lname+
                   " not found in table::get_col_struct().").c_str(),
                  exc_enotfound);
        return 0;
      }
      return &(it->second);
    }
    /// Return the beginning of the column tree
    aiter begin() { return atree.begin(); }
    /// Return the end of the column tree
    aiter end() { return atree.end(); }
    //@}

    /// An empty vector for get_column()
    vec_t empty_col;

    /// \name Interpolation
    //@{
    /// True if the interpolation object is up-to-date
    bool intp_set;
    
    /// Current interpolation type
    size_t itype;

    /// Interpolation object
    interp_base<vec_t,vec_t,fp_t> *si;

    /// The last x-column interpolated
    std::string intp_colx;

    /// The last y-column interpolated
    std::string intp_coly;
    //@}

#endif

  };

  /** \brief View a o2scl::table object as a matrix

      When connected to a table, the \ref size1() function returns the
      number of lines in the table, and the \ref size2() function
      returns the number of columns which were originally specified in
      either the constructor or the \ref set() function.

      \note This stores a pointer to the table and the user must ensure
      that the pointer is valid with the matrix view is accessed.
  */
  template<class vec_t=std::vector<double>, class fp_t=double> 
  class const_matrix_view_table : public const_matrix_view {
  
  protected:
  
    /// The number of columns
    size_t nc;
    /// The number of lines in the table
    size_t nlines;
    /// Pointers to each column
    std::vector<const vec_t *> col_ptrs;
    
  public:
    
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    const_matrix_view_table() {
      nc=0;
      nlines=0;
    }
    
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    const_matrix_view_table(o2scl::table<vec_t,fp_t> &t,
                            std::vector<std::string> cols) {
      set(t,cols);
    }
  
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    void set(o2scl::table<vec_t,fp_t> &t,
             std::vector<std::string> cols) {
      nc=cols.size();
      col_ptrs.resize(nc);
      for(size_t i=0;i<nc;i++) {
        col_ptrs[i]=&t[cols[i]];
      }
      nlines=t.get_nlines();
    }
  
    /** \brief Return the number of rows
     */
    size_t size1() const {
      return nlines;
    }
  
    /** \brief Return the number of columns
     */
    size_t size2() const {
      if (nlines==0) return 0;
      return nc;
    }
  
    /** \brief Return a reference to the element at row \c row
        and column \c col
    */
    const fp_t &operator()(size_t row, size_t col) const {
      if (row>=nlines) {
        std::string str=((std::string)"Row ")+o2scl::szttos(row)+
          " >= "+o2scl::szttos(nlines)+" in const_matrix_view_table"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      if (col>=nc) {
        std::string str=((std::string)"Column ")+o2scl::szttos(col)+
          " >= "+o2scl::szttos(nc)+" in const_matrix_view_table"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      const vec_t *cp=col_ptrs[col];
      return (*cp)[row];
    }

    /** \brief Swap method
     */
    friend void swap(const_matrix_view_table &t1,
                     const_matrix_view_table &t2) {
      using std::swap;
      swap(t1.nc,t2.nc);
      swap(t1.nlines,t2.nlines);
      swap(t1.col_ptrs,t2.col_ptrs);
      return;
    }
  
  };
  
  /** \brief View a o2scl::table object as a matrix

      When connected to a table, the \ref size1() function returns the
      number of lines in the table, and the \ref size2() function
      returns the number of columns which were originally specified in
      either the constructor or the \ref set() function.

      \note This stores a pointer to the table and the user must ensure
      that the pointer is valid with the matrix view is accessed.
  */
  template<class vec_t=std::vector<double>, class fp_t=double> 
  class matrix_view_table : public matrix_view {
  
  protected:
  
    /// The number of columns
    size_t nc;
    /// The number of lines in the table
    size_t nlines;
    /// Pointers to each column
    std::vector<vec_t *> col_ptrs;
    
  public:
    
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    matrix_view_table() {
      nc=0;
      nlines=0;
    }
    
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    matrix_view_table(o2scl::table<vec_t,fp_t> &t,
                      std::vector<std::string> cols) {
      set(t,cols);
    }
  
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    void set(o2scl::table<vec_t,fp_t> &t,
             std::vector<std::string> cols) {
      nc=cols.size();
      if (nc==0) {
        O2SCL_ERR("No columns in matrix_view_table::set().",
                  o2scl::exc_einval);
      }
      col_ptrs.resize(nc);
      for(size_t i=0;i<nc;i++) {
        col_ptrs[i]=&t.get_column_no_const(cols[i]);
      }
      nlines=t.get_nlines();
    }
  
    /** \brief Return the number of rows
     */
    size_t size1() const {
      return nlines;
    }
  
    /** \brief Return the number of columns
     */
    size_t size2() const {
      if (nlines==0) return 0;
      return nc;
    }
  
    /** \brief Return a reference to the element at row \c row
        and column \c col
    */
    const fp_t &operator()(size_t row, size_t col) const {
      if (row>=nlines) {
        std::string str=((std::string)"Row ")+o2scl::szttos(row)+
          " >= "+o2scl::szttos(nlines)+" in matrix_view_table"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      if (col>=nc) {
        std::string str=((std::string)"Column ")+o2scl::szttos(col)+
          " >= "+o2scl::szttos(nc)+" in matrix_view_table"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      const vec_t *cp=col_ptrs[col];
      return (*cp)[row];
    }

    /** \brief Return a reference to the element at row \c row
        and column \c col
    */
    fp_t &operator()(size_t row, size_t col) {
      if (row>=nlines) {
        std::string str=((std::string)"Row ")+o2scl::szttos(row)+
          " >= "+o2scl::szttos(nlines)+" in matrix_view_table"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      if (col>=nc) {
        std::string str=((std::string)"Column ")+o2scl::szttos(col)+
          " >= "+o2scl::szttos(nc)+" in matrix_view_table"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      vec_t *cp=col_ptrs[col];
      return (*cp)[row];
    }

    template<class vec2_t> friend void std::swap
    (matrix_view_table<vec2_t> &t1,
     matrix_view_table<vec2_t> &t2);
    
  };
  
  /** \brief View a o2scl::table object as a transposed matrix

      When connected to a table, the \ref size1() function returns the
      number of columns which were originally specified in either the
      constructor or the \ref set() function, and the \ref size2()
      function returns number of lines in the table.

      \note This stores a pointer to the table and the user must ensure
      that the pointer is valid with the matrix view is accessed.
  */
  template<class vec_t=std::vector<double>, class fp_t=double> 
  class matrix_view_table_transpose : public matrix_view {
  
  protected:
  
    /// The number of rows
    size_t nr;
    /// The number of lines in the table
    size_t nlines;
    /// Pointers to each column
    std::vector<vec_t *> col_ptrs;
    
  public:
    
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    matrix_view_table_transpose() {
      nr=0;
      nlines=0;
    }
    
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    matrix_view_table_transpose(o2scl::table<vec_t,fp_t> &t,
                                std::vector<std::string> rows) {
      set(t,rows);
    }
  
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    void set(o2scl::table<vec_t,fp_t> &t,
             std::vector<std::string> rows) {
      nr=rows.size();
      col_ptrs.resize(nr);
      for(size_t i=0;i<nr;i++) {
        col_ptrs[i]=&t.get_column_no_const(rows[i]);
      }
      nlines=t.get_nlines();
    }
  
    /** \brief Return the number of rows
     */
    size_t size1() const {
      if (nlines==0) return 0;
      return nr;
    }
  
    /** \brief Return the number of columns
     */
    size_t size2() const {
      return nlines;
    }
  
    /** \brief Return a reference to the element at row \c row
        and column \c col
    */
    const fp_t &operator()(size_t row, size_t col) const {
      if (row>=nr) {
        std::string str=((std::string)"Row ")+o2scl::szttos(row)+
          " >= "+o2scl::szttos(nr)+" in matrix_view_table_transpose"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      if (col>=nlines) {
        std::string str=((std::string)"Column ")+o2scl::szttos(col)+
          " >= "+o2scl::szttos(nlines)+" in matrix_view_table_transpose"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      const vec_t *rp=col_ptrs[row];
      return (*rp)[col];
    }

    /** \brief Return a reference to the element at row \c row
        and column \c col
    */
    fp_t &operator()(size_t row, size_t col) {
      if (row>=nr) {
        std::string str=((std::string)"Row ")+o2scl::szttos(row)+
          " >= "+o2scl::szttos(nr)+" in matrix_view_table_transpose"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      if (col>=nlines) {
        std::string str=((std::string)"Column ")+o2scl::szttos(col)+
          " >= "+o2scl::szttos(nlines)+" in matrix_view_table_transpose"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      vec_t *rp=col_ptrs[row];
      return (*rp)[col];
    }

    template<class vec2_t> friend
    void std::swap(matrix_view_table_transpose<vec2_t> &t1,
                   matrix_view_table_transpose<vec2_t> &t2);
    
  };
  
  /** \brief Swap part of a o2scl::table object into a matrix
   */
  class matrix_swap_table : public matrix_view {
  
  protected:
  
    /// The number of columns
    size_t nc;
    /// The number of lines in the table
    size_t nlines;
    /// Array of columns
    std::vector<std::vector<double> > cols;
    
  public:
    
    /** \brief Create a matrix view object from the specified 
	table and list of columns
    */
    matrix_swap_table() {
      nc=0;
      nlines=0;
    }
    
    /** \brief Create a matrix view object from the specified 
	table and list of columns
    */
    matrix_swap_table(o2scl::table<std::vector<double> > &t,
		      const std::vector<std::string> &col_names) {
      set(t,col_names);
    }
  
    /** \brief Create a matrix view object from the specified 
	table and list of columns
    */
    void set(o2scl::table<std::vector<double> > &t,
	     const std::vector<std::string> &col_names) {
      nc=cols.size();
      nlines=t.get_nlines();
      cols.resize(nc);
      for(size_t i=0;i<nc;i++) {
	// We have to make space for the vector before we use
	// table::swap_column_data()
	cols[i].resize(nlines);
	t.swap_column_data(col_names[i],cols[i]);
      }
    }
  
    /** \brief Return the number of rows
     */
    size_t size1() const {
      return nlines;
    }
  
    /** \brief Return the number of columns
     */
    size_t size2() const {
      if (nlines==0) return 0;
      return nc;
    }
  
    /** \brief Return a reference to the element at row \c row
	and column \c col
    */
    const double &operator()(size_t row, size_t col) const {
      if (row>=nlines) {
	std::string str=((std::string)"Row ")+o2scl::szttos(row)+
	  " >= "+o2scl::szttos(nlines)+" in matrix_swap_table"+
	  "::operator().";
	O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      if (col>=nc) {
	std::string str=((std::string)"Column ")+o2scl::szttos(col)+
	  " >= "+o2scl::szttos(nc)+" in matrix_swap_table"+
	  "::operator().";
	O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      return cols[col][row];
    }
  
    /** \brief Return a reference to the element at row \c row
	and column \c col
    */
    double &operator()(size_t row, size_t col) {
      if (row>=nlines) {
	std::string str=((std::string)"Row ")+o2scl::szttos(row)+
	  " >= "+o2scl::szttos(nlines)+" in matrix_swap_table"+
	  "::operator().";
	O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      if (col>=nc) {
	std::string str=((std::string)"Column ")+o2scl::szttos(col)+
	  " >= "+o2scl::szttos(nc)+" in matrix_swap_table"+
	  "::operator().";
	O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      return cols[col][row];
    }

    /** \brief Swap method
     */
    friend void swap(matrix_swap_table &t1,
		     matrix_swap_table &t2) {
      using std::swap;
      swap(t1.nc,t2.nc);
      swap(t1.nlines,t2.nlines);
      swap(t1.cols,t2.cols);
      return;
    }
  
  };
  
  /** \brief View a o2scl::table object as a matrix

      \note This stores a pointer to the table and the user must
      ensure that the pointer is valid with the matrix view is
      accessed.
  */
  template<class vec_t=std::vector<double>, class fp_t=double> 
  class const_matrix_view_table_transpose :
    public const_matrix_view {
  
  protected:
  
    /** \brief The number of rows in the matrix (equal to the number of 
        pointers to table columns)
    */
    size_t nr;
    /// Pointers to each row
    std::vector<const vec_t *> col_ptrs;
    /// Number of lines in the table
    size_t nlines;
  
  public:

    /** \brief Create a matrix view object from the specified 
        table and list of rows
    */
    const_matrix_view_table_transpose() {
      nr=0;
      nlines=0;
    }

    /** \brief Create a matrix view object from the specified 
        table and list of rows
    */
    const_matrix_view_table_transpose(o2scl::table<vec_t,fp_t> &t,
                                      std::vector<std::string> rows) {
      set(t,rows);
    }
  
    /** \brief Create a matrix view object from the specified 
        table and list of columns
    */
    void set(o2scl::table<vec_t,fp_t> &t,
             std::vector<std::string> rows) {
      nr=rows.size();
      col_ptrs.resize(nr);
      for(size_t i=0;i<nr;i++) {
        col_ptrs[i]=&t[rows[i]];
      }
      nlines=t.get_nlines();
    }
  
    /** \brief Return the number of rows
     */
    size_t size1() {
      if (nlines==0) return 0;
      return nr;
    }
  
    /** \brief Return the number of columns
     */
    size_t size2() {
      return nlines;
    }
  
    /** \brief Swap method
     */
    friend void swap(const_matrix_view_table_transpose &t1,
                     const_matrix_view_table_transpose &t2) {
      using std::swap;
      swap(t1.nr,t2.nr);
      swap(t1.nlines,t2.nlines);
      swap(t1.col_ptrs,t2.col_ptrs);
      return;
    }
  
    /** \brief Return a reference to the element at row \c row
        and column \c col
    */
    const double &operator()(size_t row, size_t col) const {
      if (col>=nlines) {
        std::string str=((std::string)"Column ")+o2scl::szttos(col)+
          " >= "+o2scl::szttos(nlines)+" in const_matrix_view_table_transpose"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      if (row>=nr) {
        std::string str=((std::string)"Row ")+o2scl::szttos(row)+
          " >= "+o2scl::szttos(nr)+" in const_matrix_view_table_transpose"+
          "::operator().";
        O2SCL_ERR(str.c_str(),o2scl::exc_einval);
      }
      const vec_t *rp=col_ptrs[row];
      return (*rp)[col];
    }
    
  };

}

namespace std {
  
  /** \brief Swap method
   */
  template<class vec2_t> void swap
  (o2scl::matrix_view_table<vec2_t> &t1,
   o2scl::matrix_view_table<vec2_t> &t2) {
    
    using std::swap;
    swap(t1.nc,t2.nc);
    swap(t1.nlines,t2.nlines);
    swap(t1.col_ptrs,t2.col_ptrs);
    return;
  }

  /** \brief Swap method
   */
  template<class vec2_t> void swap
  (o2scl::matrix_view_table_transpose<vec2_t> &t1,
   o2scl::matrix_view_table_transpose<vec2_t> &t2) {
    
    using std::swap;
    swap(t1.nr,t2.nr);
    swap(t1.nlines,t2.nlines);
    swap(t1.col_ptrs,t2.col_ptrs);
    return;
  }
  
}

#endif
