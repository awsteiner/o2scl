/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
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
#include <o2scl/vector.h>

namespace o2scl {

  /** \brief Create nicely formatted columns from a table of strings
    
      This is a brute-force approach of order \f$ \mathrm{ncols}
      \times \mathrm{nrows} \f$. The column widths and spacings of are
      computed by exhaustively examining all strings in every column.

      \verbatim embed:rst
      .. todo:: 

         - Future: Create a function which accepts delimited strings
           (e.g. like csv) instead of vector<vector<string>>. 
         - Future: Consider a function which takes a \ref o2scl::table
           object as input?
         - Future: Make a LaTeX output?

      \endverbatim
  */
  class columnify {

  public:

    columnify() {
      table_lines=0;
      verbose=0;
    }

    /** \brief Specification for table lines (experimental; default 0)

        If table_lines is 1, then ascii characters are used to create
        table lines. If table_lines is 2, then terminal sequences are
        used to create graphical lines.
     */
    int table_lines;

    /// Align the left-hand sides
    static const int align_left=1;
    /// Align the right-hand sides
    static const int align_right=2;
    /// Center, slightly to the left if spacing is uneven
    static const int align_lmid=3;
    /// Center, slightly to the right if spacing is uneven
    static const int align_rmid=4;
    /** \brief Align entries with decimal points. For header rows,
        align to the left.
    */
    static const int align_dp=5;
    /** \brief Align negative numbers to the left and use a space for 
	positive numbers. Align headers to the left-most digit.
    */
    static const int align_lnum=6;
    /** \brief Pre-process the columns to decide their
        alignment

        This option uses \ref align_left for columns with text, \ref
        align_lnum for columns of all floating point numbers, \ref
        align_right for integers, and \ref align_dp otherwise.
        Header rows are ignored when deciding the alignment.
    */
    static const int align_auto=7;

    /// Verbosity parameter (default 0)
    int verbose;
    
    /** \brief Take \c table and create a new object \c ctable with 
	appropriately formatted columns

        This function takes the input table and aligns all of the
        columns according to the specifications in \c align_spec
        and then places the strings for all of the final rows
        in ctable.
        
	The table of strings should be stored in \c table in
	"column-major" order (<tt>table[ncols][nrows]</tt>), so that
	\c table has the interpretation of a set of columns to be
	aligned. Before calling align(), \c ctable should be allocated
	so that at least the first \c nrows entries can be assigned,
	and \c align_spec should contain \c ncols entries specifying
	the style of alignment for each column.

	The first argument, of type \c mat_string_t, can be any type
	which is accessible using two applications of
	<tt>operator[]</tt>, such as <tt>string **</tt>,
	<tt>vector<string>[]</tt>, or <tt>vector<vector<string>></tt>.
	The fourth type must be a type which can be resized with one
	<tt>resize()</tt> call (of size \c ncols) and then for the
	first index and then an additional \c ncols calls to
	<tt>resize()</tt> of size \c nrows for example
	<tt>vector<string></tt>. The fifth type, \c vec_int_t, may be
	any integer array type.
    */
    template<class mat_string_t, class vec_string_t, class vec_int_t>
    int align(const mat_string_t &table, size_t ncols, size_t nrows, 
              vec_string_t &ctable, vec_int_t &align_spec,
              size_t n_headers=0) {

      // For the terminal table lines
      terminal ter;

      std::vector<std::vector<std::string>> table_out;
      
      int xret=add_spaces(table,ncols,nrows,align_spec,
                          table_out,n_headers);
      if (xret!=0) {
        return xret;
      }
      ctable.clear();

      for(size_t j=0;j<nrows;j++) {
        std::string tmp=table_out[0][j];
        for(size_t i=1;i<ncols;i++) {
          if (table_lines>1) {
            tmp+=ter.alt_font()+" x "+ter.normal_font()+table_out[i][j];
          } else if (table_lines>0) {
            tmp+=" | "+table_out[i][j];
          } else {
            tmp+=' '+table_out[i][j];
          }
        }
        ctable.push_back(tmp);
        if (table_lines>1 && j+1==n_headers) {
          tmp=ter.alt_font();
          for(size_t k=0;k<ter.str_len(table_out[0][j]);k++) {
            tmp+="q";
          }
          for(size_t i=1;i<ncols;i++) {
            tmp+="qnq";
            for(size_t k=0;k<ter.str_len(table_out[i][j]);k++) {
              tmp+="q";
            }
          }
          tmp+=ter.normal_font();
          ctable.push_back(tmp);
        } else if (table_lines>0 && j+1==n_headers) {
          for(size_t k=0;k<ter.str_len(table_out[0][j]);k++) {
            tmp+="-";
          }
          for(size_t i=1;i<ncols;i++) {
            tmp+="-+-";
            for(size_t k=0;k<ter.str_len(table_out[i][j]);k++) {
              tmp+="-";
            }
          }
          ctable.push_back(tmp);
        }
      }
      return 0;
    }

    /** \brief For input \c col_in, add spaces to each entry to
        ensure that each string has the same width, and store
        the result in \c col_out.
    */
    template<class vec_string_t>
    int add_spaces_one(const vec_string_t &col_in, int align_spec,
                       vec_string_t &col_out, size_t n_headers=0) {
      std::vector<std::vector<std::string>> table_in, table_out;
      table_in.push_back(col_in);
      std::vector<int> valign;
      valign.push_back(align_spec);
      int ret=add_spaces(table_in,1,col_in.size(),valign,
                         table_out,n_headers);
      col_out=table_out[0];
      return ret;
    }
      
    /** \brief Add enough spaces to ensure all columns have the
	same width

	The first argument can be any type which is accessible
	using two applications of <tt>operator[]</tt>, such 
	as <tt>string **</tt>, <tt>vector<string>[]</tt>, or
	<tt>vector<vector<string>></tt>. The third type must be
        a type which can be resized with one <tt>resize()</tt>
        call (of size \c ncols) and then for the first index and
        then an additional \c ncols calls to <tt>resize()</tt>
        (of size \c nrows).

        This function makes two passes through the input table.
        The first pass constructs the maximum width for each
        column. The second pass adds spaces sufficient to ensure
        each row has the same width.
    */
    template<class mat_string_t, class mat_string2_t, class vec_int_t>
    int add_spaces(const mat_string_t &table_in, size_t ncols,
                   size_t nrows, vec_int_t &align_spec,
                   mat_string2_t &table_out, size_t n_headers=0) {

      if (n_headers>nrows) {
        O2SCL_ERR2("Number of headers larger than rows in ",
                   "columnify::add_spaces().",o2scl::exc_einval);
      }
      
      // Use this class to avoid counting vt100 sequences
      terminal ter;

      std::vector<int> align2;

      for(size_t i=0;i<ncols;i++) {
        align2[i]=align_spec[i];
        if (align_spec[i]==align_auto) {
          int cnt_int=0, cnt_fp=0, cnt_fp_sci=0, cnt_has_minus=0;
          int cnt_not_num=0;
          for(size_t j=n_headers;j<nrows;j++) {
            bool is_int, is_fp, is_fp_sci, has_minus, not_num;
            guess_type(table_in[i][j],is_int,is_fp,is_fp_sci,
                       has_minus,not_num);
            if (is_int) cnt_int++;
            if (is_fp) cnt_fp++;
            if (is_fp_sci) cnt_fp_sci++;
            if (has_minus) cnt_has_minus++;
            if (not_num) cnt_not_num++;
          }
          if (cnt_not_num>0) {
            align2[i]=align_left;
          } else if (cnt_int==0) {
            align2[i]=align_lnum;
          } else if (cnt_int==((int)(nrows-n_headers))) {
            align2[i]=align_right;
          } else {
            align2[i]=align_dp;
          }
          if (verbose>0) {
            std::cout << "columnify::add_spaces(): "
                      << "Using align " << align2[i]
                      << " for column " << i << "." << std::endl;
          }
        }
      }
      
      // Ensure table_out has enough space to hold the output
      // data
      if (table_out.size()!=ncols) table_out.resize(ncols);
      for(size_t i=0;i<ncols;i++) {
        if (table_out[i].size()!=nrows) {
          table_out[i].resize(nrows);
        }
      }
      
      // Make space for the size information. The value csizes[i] is
      // the maximum number of chararcters for column i. However, if
      // the alignment is of type "dp", then csizes[i] is the maximum
      // size of the space to the left of the decimal point, plus one
      // for the decimal point itself, and csizes2[i] is the maximum
      // size of the space to the right of the decimal point. For
      // rows in the header, when the alignment is of type "dp", the
      // location of the decimal point is ignored. 
      boost::numeric::ublas::vector<size_t> csizes(ncols);
      boost::numeric::ublas::vector<size_t> csizes2(ncols);
      for(size_t i=0;i<ncols;i++) {
	csizes[i]=0;
	csizes2[i]=0;
      }

      // If true, the alignment is of type "align_lnum" and there is a
      // negative number in the first character.
      std::vector<bool> has_negative(ncols);
      o2scl::vector_set_all(has_negative,false);
      
      // First pass: compute the sizes of all the entries in all of
      // the columns so we know how many spaces to add.
      for(size_t i=0;i<ncols;i++) {
        
        // We have to go backwards, because for align_dp, we don't
        // quite know how to align the headers until after we
        // get the sizes for all of the body rows. Similarly,
        // for align_lnum, we don't know if any of the numbers has
        // a minus sign until we go through all of the body rows.
	for(int j=nrows-1;j>=0;j--) {

          if (verbose>1) {
            std::cout << "columnify::add_spaces(): Processing: x"
                      << table_in[i][j] << "x" << std::endl;
          }
          
          // For align_lnum, see if we need to add a space to the
          // left-hand side of the header for the minus sign
	  if (align2[i]==align_lnum) {
            if (has_negative[i]==false && table_in[i][j].length()>0 &&
                table_in[i][j][0]=='-') {
              if (verbose>1) {
                std::cout << "columnify::add_spaces(): "
                          << "Setting has_negative for column "
                          << i << " to true." << std::endl;
              }
              has_negative[i]=true;
            }
          }
          
	  // If we're aligning with decimal points, we need to compute
	  // the maximum width to the left and the right of the 
	  // decimal point separately
	  if (align2[i]==align_dp) {

            //std::cout << "j,nh: " << j << " " << n_headers
            //<< std::endl;
            
            if (j<((int)n_headers)) {
              
              // If we're in a header row, ignore the decimal
              // point, but make more space for the header if
              // necessary,
              if (ter.str_len(table_in[i][j])>csizes[i]+csizes2[i]) {
                csizes[i]=ter.str_len(table_in[i][j])-csizes2[i];
              }
              
            } else {
              
              size_t loc=table_in[i][j].find('.');
              //std::cout << "loc: " << loc << std::endl;
              std::string left, right;
              if (loc!=std::string::npos) {
                left=table_in[i][j].substr(0,loc+1);
                right=table_in[i][j].substr
                  (loc+1,ter.str_len(table_in[i][j])-loc-1);
                   
              } else {
                left=table_in[i][j]+' ';
                right="";
              }
              if (ter.str_len(left)>csizes[i]) {
                csizes[i]=ter.str_len(left);
              }
              if (ter.str_len(right)>csizes2[i]) {
                csizes2[i]=ter.str_len(right);
              }

            }

	  } else {

	    // If the alignment is not of type "align_dp", then just
	    // find the maximum width of each column.
	    if (ter.str_len(table_in[i][j])>csizes[i]) {
	      csizes[i]=ter.str_len(table_in[i][j]);
	    }

            /// Adjust for lnum if the largest width field is
            // a positive number
            if (align2[i]==align_lnum &&
                table_in[i][j].length()>0 &&
                table_in[i][j][0]!='-' &&
                ter.str_len(table_in[i][j])==csizes[i]) {
              csizes[i]++;
            }

	  }

          if (verbose>1) {
            std::cout << "columnify::add_spaces(): "
                      << "i,csizes[i],csizes2[i]: " << i
                      << " " << csizes[i] << " " << csizes2[i]
                      << std::endl;
          }
	}

      }

      // Second pass: go through row by row, adding enough spaces to
      // ensure every row size matches for each column.
      
      for(size_t j=0;j<nrows;j++) {
	for(size_t i=0;i<ncols;i++) {
          
          if (verbose>1) {
            std::cout << "columnify::add_spaces(): ";
            std::cout << "Processing (second round): x" << table_in[i][j]
                      << "x " << csizes[i] << " " << csizes2[i] << std::endl;
          }

          // The final string for this row and column
	  std::string tmp="";
	  
	  // Handle each alignment case separately
	  if (align2[i]==align_right) {

	    for(size_t k=ter.str_len(table_in[i][j]);k<csizes[i];k++) {
	      tmp+=' ';
	    }
	    tmp+=table_in[i][j];

	  } else if (align2[i]==align_left) {

	    tmp+=table_in[i][j];
	    for(size_t k=ter.str_len(table_in[i][j]);k<csizes[i];k++) {
	      tmp+=' ';
	    }

	  } else if (align2[i]==align_lmid) {

	    size_t le=(csizes[i]-ter.str_len(table_in[i][j]))/2;
	    size_t ri=csizes[i]-ter.str_len(table_in[i][j])-le;
	    for(size_t k=0;k<le;k++) tmp+=' ';
	    tmp+=table_in[i][j];
	    for(size_t k=0;k<ri;k++) tmp+=' ';

	  } else if (align2[i]==align_rmid) {

	    size_t ri=(csizes[i]-ter.str_len(table_in[i][j]))/2;
	    size_t le=csizes[i]-ter.str_len(table_in[i][j])-ri;
	    for(size_t k=0;k<le;k++) tmp+=' ';
	    tmp+=table_in[i][j];
	    for(size_t k=0;k<ri;k++) tmp+=' ';

	  } else if (align2[i]==align_dp) {

            if (j<n_headers) {

              tmp=table_in[i][j];
              while (ter.str_len(tmp)<csizes[i]+csizes2[i]) {
                tmp+=' ';
              }
              
            } else {
            
              size_t loc=table_in[i][j].find('.');
              std::string left, right;
              if (loc!=std::string::npos) {
                left=table_in[i][j].substr(0,loc+1);
                right=table_in[i][j].substr
                  (loc+1,ter.str_len(table_in[i][j])-loc-1);
                
              } else {
                left=table_in[i][j]+' ';
                right="";
              }
              
              for(size_t k=ter.str_len(left);k<csizes[i];k++) tmp+=' ';
              tmp+=left;
              tmp+=right;
              for(size_t k=ter.str_len(right);k<csizes2[i];k++) tmp+=' ';
            }

	  } else if (align2[i]==align_lnum) {
            
            if (j<n_headers) {
              if (table_in[i][j].length()>0 && has_negative[i]==true &&
                  table_in[i][j].length()<csizes[i]) {
              }
            } 
            
	    if (ter.str_len(table_in[i][j])==csizes[i]) {
              
              // If the string is already max width, then just
              // append it directly
	      tmp+=table_in[i][j];
              
	    } else {

              // In header rows, add a space at the beginning if there
              // is a negative value. Otherwise, add a space at the
              // beginning if the first character is a digit.
              if (has_negative[i]==true &&
                  (j<n_headers ||
                   (table_in[i][j].length()>0 &&
                    table_in[i][j][0]>='0' && table_in[i][j][0]<='9'))) {
                tmp+=' ';
              }
              
              // Otherwise, we can just add spaces at the right
              tmp+=table_in[i][j];
              for(size_t k=ter.str_len(tmp);k<csizes[i];k++) {
                tmp+=' ';
              }
              
            }
	  }
	  
          if (verbose>1) {
            std::cout << "columnify::add_spaces(): ";
            std::cout << "Storing: x" << tmp << "x" << std::endl;
                      
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

      \verbatim embed:rst
      .. todo:: 
         Future: If all of the matrix elements are positive integers 
         and scientific mode is not set, then we can avoid printing
         the extra spaces.
      \endverbatim
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

      \verbatim embed:rst
      .. todo:: 
         Future: If all of the matrix elements are positive integers 
         and scientific mode is not set, then we can avoid printing
         the extra spaces.
      \endverbatim

      This function outputs the matrix assuming the first index is the
      column index and the second index is the row index. For the
      opposite convention, use \ref array_2d_out().
  */
  template<class mat_t> void array_2d_trans_out
  (std::ostream &os, size_t nrows, size_t ncols, mat_t &A) {
   
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
