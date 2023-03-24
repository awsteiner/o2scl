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
#ifndef O2SCL_COLUMNIFY_H
#define O2SCL_COLUMNIFY_H

/** \file columnify.h
    \brief Class which formats strings into columns
*/

#include <iostream>
#include <string>
#include <vector>

#define BOOST_DISABLE_ASSERTS
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/misc.h>
#include <o2scl/string_conv.h>

namespace o2scl {

  /** \brief Create nicely formatted columns from a table of strings
    
      This is a brute-force approach of order \f$ \mathrm{ncols}
      \times \mathrm{nrows} \f$. The column widths and spacings of are
      computed by exhaustively examining all strings in every column.

      \future Create a single column version of add_spaces().
      \future Create a function which accepts delimited strings
      (e.g. like csv) instead of vector<vector<string>>. 
      \future Move the screenify() functionality from misc.h into 
      this class?
      \future It might be better to allow the string table
      to be specified with iterators.
      \future Consider a function which takes a \ref o2scl::table
      object as input?
  */
  class columnify {

  public:

    columnify() {
      table_lines=0;
    }

    /// Specification for table lines (experimental)
    int table_lines;

    /// Align the left-hand sides
    static const int align_left=1;
    /// Align the right-hand sides
    static const int align_right=2;
    /// Center, slightly to the left if spacing is uneven
    static const int align_lmid=3;
    /// Center, slightly to the right if spacing is uneven
    static const int align_rmid=4;
    /// Align with decimal points
    static const int align_dp=5;
    /** \brief Align negative numbers to the left and use a space for 
	positive numbers
    */
    static const int align_lnum=6;
  
    /** \brief Take \c table and create a new object \c ctable with 
	appropriately formatted columns

	The table of strings should be stored in \c table in
	"column-major" order (<tt>table[ncols][nrows]</tt>), so that
	\c table has the interpretation of a set of columns to be
	aligned. Before calling align(), \c ctable should be allocated
	so that at least the first \c nrows entries can be assigned,
	and \c align_spec should contain \c ncols entries specifying
	the style of alignment for each column.

	The first argument can be any type which is accessible
	using two applications of <tt>operator[]</tt>, such 
	as <tt>string **</tt>, <tt>vector<string>[]</tt>, or
	<tt>vector<vector<string> > </tt>
    */
    template<class mat_string_t, class vec_string_t, class vec_int_t>
      int align(const mat_string_t &table, size_t ncols, size_t nrows, 
		vec_string_t &ctable, vec_int_t &align_spec) {

      terminal ter;
      
      // Make space for the size information
      boost::numeric::ublas::vector<size_t> csizes(ncols);
      boost::numeric::ublas::vector<size_t> csizes2(ncols);
      for(size_t i=0;i<ncols;i++) {
	csizes[i]=0;
	csizes2[i]=0;
      }
  
      // Compute the sizes of all the entries in all of the columns so
      // we know how many spaces to add
      for(size_t i=0;i<ncols;i++) {
	for(size_t j=0;j<nrows;j++) {

	  // If we're aligning with decimal points, we need to compute
	  // the maximum width to the left and the right of the 
	  // decimal point separately

	  if (align_spec[i]==align_dp) {
	    size_t loc=table[i][j].find('.');
	    std::string left, right;
	    if (loc!=std::string::npos) {
	      left=table[i][j].substr(0,loc+1);
	      right=table[i][j].substr(loc+1,
				       table[i][j].length()-loc-1);
	    } else {
	      left=table[i][j]+' ';
	      right="";
	    }
	    if (left.length()>csizes[i]) csizes[i]=left.length();
	    if (right.length()>csizes2[i]) csizes2[i]=right.length();

	  } else {

	    // Otherwise just find the maximum width of each column
	    if (table[i][j].length()>csizes[i]) {
	      csizes[i]=table[i][j].length();
	    }

	  }

	}
      }

      // Go through row by row, adding enough spaces to make one string
      // per row
      for(size_t j=0;j<nrows;j++) {

	std::string tmp="";

	for(size_t i=0;i<ncols;i++) {

	  // Handle each alignment case separately
	  if (align_spec[i]==align_right) {

	    for(size_t k=table[i][j].length();k<csizes[i];k++) {
	      tmp+=' ';
	    }
	    tmp+=table[i][j];

	  } else if (align_spec[i]==align_left) {

	    tmp+=table[i][j];
	    for(size_t k=table[i][j].length();k<csizes[i];k++) {
	      tmp+=' ';
	    }

	  } else if (align_spec[i]==align_lmid) {

	    size_t le=(csizes[i]-table[i][j].length())/2;
	    size_t ri=csizes[i]-table[i][j].length()-le;
	    for(size_t k=0;k<le;k++) tmp+=' ';
	    tmp+=table[i][j];
	    for(size_t k=0;k<ri;k++) tmp+=' ';

	  } else if (align_spec[i]==align_rmid) {

	    size_t ri=(csizes[i]-table[i][j].length())/2;
	    size_t le=csizes[i]-table[i][j].length()-ri;
	    for(size_t k=0;k<le;k++) tmp+=' ';
	    tmp+=table[i][j];
	    for(size_t k=0;k<ri;k++) tmp+=' ';

	  } else if (align_spec[i]==align_dp) {

	    size_t loc=table[i][j].find('.');
	    std::string left, right;
	    if (loc!=std::string::npos) {
	      left=table[i][j].substr(0,loc+1);
	      right=table[i][j].substr(loc+1,
				       table[i][j].length()-loc-1);
	    } else {
	      left=table[i][j]+' ';
	      right="";
	    }

	    for(size_t k=left.length();k<csizes[i];k++) tmp+=' ';
	    tmp+=left;
	    tmp+=right;
	    for(size_t k=right.length();k<csizes2[i];k++) tmp+=' ';

	  } else if (align_spec[i]==align_lnum) {

	    if (table[i][j].length()==csizes[i]) {
	      tmp+=table[i][j];
	    } else {
	      if (table[i][j][0]>='0' && table[i][j][0]<='9') {
		tmp+=' ';
		tmp+=table[i][j];
		for(size_t k=table[i][j].length();k<csizes[i]-1;k++) {
		  tmp+=' ';
		}
	      } else {
		tmp+=table[i][j];
		for(size_t k=table[i][j].length();k<csizes[i];k++) {
		  tmp+=' ';
		}
	      }
	    }
	  }
	  
	  if (i!=ncols-1) {
	    if (table_lines>0) {
	      tmp+=ter.alt_font()+" x "+ter.normal_font();
	    } else {
	      tmp+=' ';
	    }
	  }

	  // Proceed to the next column
	}

	// Add the row to the user-specified array and go to the next row
	ctable[j]=tmp;
      
      }

      return 0;
    }

    /** \brief Add enough spaces to ensure all columns have the
	same width
     */
    template<class mat_string_t, class vec_int_t>
    int add_spaces(const mat_string_t &table_in, size_t ncols, size_t nrows, 
		   vec_int_t &align_spec, mat_string_t &table_out) {

      // Use this class to avoid counting vt100 sequences
      terminal ter;
      
      // Make space for the size information
      boost::numeric::ublas::vector<size_t> csizes(ncols);
      boost::numeric::ublas::vector<size_t> csizes2(ncols);
      for(size_t i=0;i<ncols;i++) {
	csizes[i]=0;
	csizes2[i]=0;
      }
  
      // Compute the sizes of all the entries in all of the columns so
      // we know how many spaces to add
      for(size_t i=0;i<ncols;i++) {
	for(size_t j=0;j<nrows;j++) {

	  // If we're aligning with decimal points, we need to compute
	  // the maximum width to the left and the right of the 
	  // decimal point separately

	  if (align_spec[i]==align_dp) {
	    size_t loc=table_in[i][j].find('.');
	    std::string left, right;
	    if (loc!=std::string::npos) {
	      left=table_in[i][j].substr(0,loc+1);
	      right=table_in[i][j].substr(loc+1,
					  ter.str_len(table_in[i][j])-loc-1);
	    } else {
	      left=table_in[i][j]+' ';
	      right="";
	    }
	    if (ter.str_len(left)>csizes[i]) csizes[i]=ter.str_len(left);
	    if (ter.str_len(right)>csizes2[i]) csizes2[i]=ter.str_len(right);

	  } else {

	    // Otherwise just find the maximum width of each column
	    if (ter.str_len(table_in[i][j])>csizes[i]) {
	      csizes[i]=ter.str_len(table_in[i][j]);
	    }

	  }

	}
      }

      // Go through row by row, adding enough spaces to make one string
      // per row
      for(size_t j=0;j<nrows;j++) {

	for(size_t i=0;i<ncols;i++) {

	  std::string tmp="";
	  
	  // Handle each alignment case separately
	  if (align_spec[i]==align_right) {

	    for(size_t k=ter.str_len(table_in[i][j]);k<csizes[i];k++) {
	      tmp+=' ';
	    }
	    tmp+=table_in[i][j];

	  } else if (align_spec[i]==align_left) {

	    tmp+=table_in[i][j];
	    for(size_t k=ter.str_len(table_in[i][j]);k<csizes[i];k++) {
	      tmp+=' ';
	    }

	  } else if (align_spec[i]==align_lmid) {

	    size_t le=(csizes[i]-ter.str_len(table_in[i][j]))/2;
	    size_t ri=csizes[i]-ter.str_len(table_in[i][j])-le;
	    for(size_t k=0;k<le;k++) tmp+=' ';
	    tmp+=table_in[i][j];
	    for(size_t k=0;k<ri;k++) tmp+=' ';

	  } else if (align_spec[i]==align_rmid) {

	    size_t ri=(csizes[i]-ter.str_len(table_in[i][j]))/2;
	    size_t le=csizes[i]-ter.str_len(table_in[i][j])-ri;
	    for(size_t k=0;k<le;k++) tmp+=' ';
	    tmp+=table_in[i][j];
	    for(size_t k=0;k<ri;k++) tmp+=' ';

	  } else if (align_spec[i]==align_dp) {

	    size_t loc=table_in[i][j].find('.');
	    std::string left, right;
	    if (loc!=std::string::npos) {
	      left=table_in[i][j].substr(0,loc+1);
	      right=table_in[i][j].substr(loc+1,
				       ter.str_len(table_in[i][j])-loc-1);
	    } else {
	      left=table_in[i][j]+' ';
	      right="";
	    }

	    for(size_t k=ter.str_len(left);k<csizes[i];k++) tmp+=' ';
	    tmp+=left;
	    tmp+=right;
	    for(size_t k=ter.str_len(right);k<csizes2[i];k++) tmp+=' ';

	  } else if (align_spec[i]==align_lnum) {

	    if (ter.str_len(table_in[i][j])==csizes[i]) {
	      tmp+=table_in[i][j];
	    } else {
	      if (table_in[i][j][0]>='0' && table_in[i][j][0]<='9') {
		tmp+=' ';
		tmp+=table_in[i][j];
		for(size_t k=ter.str_len(table_in[i][j]);k<csizes[i]-1;k++) {
		  tmp+=' ';
		}
	      } else {
		tmp+=table_in[i][j];
		for(size_t k=ter.str_len(table_in[i][j]);k<csizes[i];k++) {
		  tmp+=' ';
		}
	      }
	    }
	  }
	  
	  table_out[i][j]=tmp;
	  
	  // Proceed to the next column
	}

	// Proceed to the next row
      }

      return 0;
    }

  };
    
  /// \name Matrix output functions from src/base/columnify.h
  //@{
  /** \brief A operator for simple matrix output using \c operator()
      
      The type \c mat_t can be any matrix type which allows 
      individual element access using <tt>operator()(size_t,size_t)</tt>.
    
      This outputs all of the matrix elements using output settings
      specified by \c os. The alignment performed by \ref columnify
      using columnify::align_dp, i.e. the numbers are aligned by
      their decimal points. If the numbers have no decimal points,
      then the decimal point is assumed to be to the right of the
      last character in the string representation of the number.

      This function outputs the matrix assuming the first index is the
      row index and the second index is the column index. For the
      opposite convention, use \ref matrix_trans_out().
  */
  template<class mat_t> void matrix_out(std::ostream &os, size_t nrows, 
                                        size_t ncols, const mat_t &A,
                                        std::string prefix="") {
    
    columnify co;
    std::vector<std::vector<std::string> > stab(ncols);
    std::vector<std::string> ctable(nrows);
    std::vector<int> alig(ncols);
    
    for(size_t j=0;j<ncols;j++) {
      alig[j]=columnify::align_dp;
      for(size_t i=0;i<nrows;i++) {
	stab[j].push_back(dtos(A(i,j),os));
      }
    }
    co.align(stab,ncols,nrows,ctable,alig);
    for(size_t i=0;i<nrows;i++) {
      os << prefix << ctable[i] << std::endl;
    }

    return;
  }

  /** \brief A operator for simple matrix output using \c operator()
      
      The type \c mat_t can be any matrix type which allows 
      individual element access using <tt>operator()(size_t,size_t)</tt>
      and access to the number of columns and rows using
      <tt>A.size1()</tt> and <tt>A.size2()</tt>.
      
      This outputs all of the matrix elements using output settings
      specified by \c os. The alignment performed by \ref columnify
      using columnify::align_dp, i.e. the numbers are aligned by
      their decimal points. If the numbers have no decimal points,
      then the decimal point is assumed to be to the right of the
      last character in the string representation of the number.

      This function outputs the matrix assuming the first index is the
      row index and the second index is the column index. For the
      opposite convention, use \ref matrix_trans_out().
  */
  template<class mat_t> void matrix_out(std::ostream &os, const mat_t &A,
                                        std::string prefix="") {

    size_t nrows=A.size1();
    size_t ncols=A.size2();
    
    columnify co;
    std::vector<std::vector<std::string> > stab(ncols);
    std::vector<std::string> ctable(nrows);
    std::vector<int> alig(ncols);
    
    for(size_t j=0;j<ncols;j++) {
      alig[j]=columnify::align_dp;
      for(size_t i=0;i<nrows;i++) {
	stab[j].push_back(dtos(A(i,j),os));
      }
    }
    co.align(stab,ncols,nrows,ctable,alig);
    for(size_t i=0;i<nrows;i++) {
      os << prefix << ctable[i] << std::endl;
    }

    return;
  }

  /** \brief A operator for simple matrix output using \c operator()
      
      The type \c mat_t can be any matrix type which allows 
      individual element access using <tt>operator()(size_t,size_t)</tt>.
      
      This outputs all of the matrix elements using output settings
      specified by \c os. The alignment performed by \ref columnify
      using columnify::align_dp, i.e. the numbers are aligned by
      their decimal points. If the numbers have no decimal points,
      then the decimal point is assumed to be to the right of the
      last character in the string representation of the number.

      This function outputs the matrix assuming the first index is the
      column index and the second index is the row index. For the
      opposite convention, use \ref matrix_out().
  */
  template<class mat_t> void matrix_trans_out(std::ostream &os, size_t nrows, 
                                              size_t ncols, const mat_t &A,
                                              std::string prefix="") {
   
    columnify co;
    std::vector<std::vector<std::string> > stab(nrows);
    std::vector<std::string> ctable(ncols);
    std::vector<int> alig(nrows);
    
    for(size_t i=0;i<nrows;i++) {
      alig[i]=columnify::align_dp;
      for(size_t j=0;j<ncols;j++) {
	stab[i].push_back(dtos(A(i,j),os));
      }
    }
    co.align(stab,nrows,ncols,ctable,alig);
    for(size_t i=0;i<ncols;i++) {
      os << prefix << ctable[i] << std::endl;
    }

    return;
  }

  /** \brief A operator for simple matrix output using \c operator()
      
      The type \c mat_t can be any matrix type which allows 
      individual element access using <tt>operator()(size_t,size_t)</tt>
      and access to the number of columns and rows using
      <tt>A.size1()</tt> and <tt>A.size2()</tt>.
      
      This outputs all of the matrix elements using output settings
      specified by \c os. The alignment performed by \ref columnify
      using columnify::align_dp, i.e. the numbers are aligned by
      their decimal points. If the numbers have no decimal points,
      then the decimal point is assumed to be to the right of the
      last character in the string representation of the number.

      This function outputs the matrix assuming the first index is the
      column index and the second index is the row index. For the
      opposite convention, use \ref matrix_out().
  */
  template<class mat_t> void matrix_trans_out(std::ostream &os,
                                              const mat_t &A,
                                              std::string prefix="") {

    size_t nrows=A.size1();
    size_t ncols=A.size2();
    
    columnify co;
    std::vector<std::vector<std::string> > stab(nrows);
    std::vector<std::string> ctable(ncols);
    std::vector<int> alig(nrows);
  
    for(size_t i=0;i<nrows;i++) {
      alig[i]=columnify::align_dp;
      for(size_t j=0;j<ncols;j++) {
	stab[i].push_back(dtos(A(i,j),os));
      }
    }
    co.align(stab,nrows,ncols,ctable,alig);
    for(size_t i=0;i<ncols;i++) {
      os << prefix << ctable[i] << std::endl;
    }

    return;
  }

  /** \brief A operator for simple matrix output using \c operator[]
      
      The type \c mat_t can be any 2d-array type which allows 
      individual element access using \c [size_t][size_t]
      
      This outputs all of the matrix elements using output settings
      specified by \c os. The alignment performed by \ref columnify
      using columnify::align_dp, i.e. the numbers are aligned by
      their decimal points. If the numbers have no decimal points,
      then the decimal point is assumed to be to the right of the
      last character in the string representation of the number.

      This function outputs the matrix assuming the first index is the
      row index and the second index is the column index. For the
      opposite convention, use \ref array_2d_trans_out().

      \future If all of the matrix elements are positive integers 
      and scientific mode is not set, then we can avoid printing
      the extra spaces.
  */
  template<class mat_t> void array_2d_out(std::ostream &os, size_t nrows, 
					 size_t ncols, mat_t &A) {

    columnify co;
    std::vector<std::string> *stab=new std::vector<std::string>[ncols];
    std::vector<std::string> ctable(nrows);
    std::vector<int> alig(ncols);
    
    for(size_t j=0;j<ncols;j++) {
      alig[j]=columnify::align_dp;
      for(size_t i=0;i<nrows;i++) {
	stab[j].push_back(dtos(A[i][j],os));
      }
    }
    co.align(stab,ncols,nrows,ctable,alig);
    for(size_t i=0;i<nrows;i++) {
      os << ctable[i] << std::endl;
    }
    
    delete[] stab;

    return;
  }

  /** \brief A operator for simple matrix output using \c operator[]
      
      The type \c mat_t can be any 2d-array type which allows 
      individual element access using \c [size_t][size_t]
      
      This outputs all of the matrix elements using output settings
      specified by \c os. The alignment performed by \ref columnify
      using columnify::align_dp, i.e. the numbers are aligned by
      their decimal points. If the numbers have no decimal points,
      then the decimal point is assumed to be to the right of the
      last character in the string representation of the number.

      \future If all of the matrix elements are positive integers 
      and scientific mode is not set, then we can avoid printing
      the extra spaces.

      This function outputs the matrix assuming the first index is the
      column index and the second index is the row index. For the
      opposite convention, use \ref array_2d_out().
  */
  template<class mat_t> void array_2d_trans_out(std::ostream &os, size_t nrows, 
					       size_t ncols, mat_t &A) {
   
    columnify co;
    std::vector<std::string> *stab=new std::vector<std::string>[nrows];
    std::vector<std::string> ctable(ncols);
    std::vector<int> alig(nrows);
    
    for(size_t i=0;i<nrows;i++) {
      alig[i]=columnify::align_dp;
      for(size_t j=0;j<ncols;j++) {
	stab[i].push_back(dtos(A[i][j],os));
      }
    }
    co.align(stab,nrows,ncols,ctable,alig);
    for(size_t i=0;i<ncols;i++) {
      os << ctable[i] << std::endl;
    }
    
    delete[] stab;

    return;
  }
  //@}

}

#endif
