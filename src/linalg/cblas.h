/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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
#ifndef O2SCL_CBLAS_H
#define O2SCL_CBLAS_H

/** \file cblas.h
    \brief Header wrapper for \ref cblas_base.h and documentation of
    \ref o2scl_cblas, \ref o2scl_cblas_bracket, \ref o2scl_linalg
    and \ref o2scl_linalg_bracket namespaces
*/

#include <cmath>
#include <o2scl/permutation.h>

/** \brief Namespace for O2scl CBLAS function templates

    These functions are only intended as a fallback for situations
    where an optimized BLAS is not available.

    <b>Level-1 BLAS functions</b>

    Some functionality which would otherwise appear here is 
    already given in \ref vector.h. 
    - The equivalent of <tt>dcopy()</tt> is given in \ref vector_copy()
    except that the ordering is reversed (in \ref vector_copy() the
    source preceeds the destination in the function argument list). 
    - The equivalent of <tt>dswap()</tt> is given in \ref vector_swap().
    - The equivalent of <tt>idamax()</tt> is given in 
    \ref vector_max_index().
    
    <b>Level-2 BLAS functions</b>

    Currently only \ref dgemv(), \ref dtrmv(), and \ref dtrsv() are 
    implemented.

    <b>Level-3 BLAS functions</b>

    Currently only \ref dgemm() and \ref dtrsm() are implemented.

    <b>Helper BLAS functions</b>

    There are several basic BLAS functions which are helpful
    to operate on only a part of a vector or matrix to 
    ensure that the linear algebra routines are flexible
    with the types that they can handle.

    The subvector functions operate only one set of adjacent 
    vector elements. For a vector defined by
    with
    \f[
    {\vec x} = 
    \left(
    \begin{array}{c}
    x_0 \\
    x_1 \\
    . \\
    . \\
    x_{\mathrm{ie}} \\
    x_{\mathrm{ie}+1} \\
    . \\
    . \\
    x_{\mathrm{N}-1} \\
    x_{\mathrm{N}} \\
    x_{\mathrm{N}+1} \\
    . \\
    .
    \end{array}
    \right)
    \f]
    the functions with suffix <tt>subvec</tt> operate only on 
    elements from \f$ x_{\mathrm{ie}} \f$ to \f$ x_{\mathrm{N}-1} \f$
    (inclusive). 

    The subcolumn functions operate only on a part of a column of
    a matrix. For a matrix defined by 
    \f[
    m = \left(						
    \begin{array}{ccccccc}
    m_{0,0} & m_{0,1} & . & . & m_{0,\mathrm{ic}} & . & . \\
    m_{1,0} & m_{1,1} & . & . & m_{1,\mathrm{ic}} & . & . \\
    . & . & . & . & . & . & . \\
    . & . & . & . & m_{\mathrm{ir},\mathrm{ic}} & . & . \\
    . & . & . & . & m_{\mathrm{ir}+1,\mathrm{ic}} & . & . \\
    . & . & . & . & . & . & . \\
    . & . & . & . & . & . & . \\
    . & . & . & . & m_{\mathrm{N}-1,\mathrm{ic}} & . & . \\
    . & . & . & . & m_{\mathrm{N},\mathrm{ic}} & . & . \\
    . & . & . & . & m_{\mathrm{N}+1,\mathrm{ic}} & . & . \\
    . & . & . & . & . & . & . \\
    \end{array}
    \right)
    \f]
    the functions with suffix <tt>subcol</tt> operate only
    on elements in the column from \f$ m_{\mathrm{ir},\mathrm{ic}} \f$
    to \f$ m_{\mathrm{N}-1,\mathrm{ic}} \f$ inclusive.

    The subrow functions operate only on a part of a row of
    a matrix. For a matrix defined by 
    \f[
    m = \left(						
    \begin{array}{ccccccccccc}
    m_{0,0} & m_{0,1} & . & . & . & . & . & . & . & . & . \\
    m_{1,0} & m_{1,1} & . & . & . & . & . & . & . & . & . \\
    . & . & . & . & . & . & . & . & . & . & . \\
    . & . & . & . & . & . & . & . & . & . & . \\
    m_{\mathrm{ir},0} & . & . &
    m_{\mathrm{ir},\mathrm{ic}} &
    m_{\mathrm{ir},\mathrm{ic}+1} &
    . & . & 
    m_{\mathrm{ir},\mathrm{N}-1} &
    m_{\mathrm{ir},\mathrm{N}} &
    m_{\mathrm{ir},\mathrm{N}+1} &
    . \\
    . & . & . & . & . & . & . & . & . & . & . \\
    \end{array}
    \right)
    \f]
    the functions with suffix <tt>subrow</tt> operate only
    on elements in the column from \f$ m_{\mathrm{ir},\mathrm{ic}} \f$
    to \f$ m_{\mathrm{ir},\mathrm{N}-1} \f$ inclusive.

    This namespace is documented inside <tt>src/linalg/cblas.h</tt>.
*/
namespace o2scl_cblas {
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M(i,j)
#include <o2scl/cblas_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2
  
}

/** \brief Namespace for O2scl CBLAS function templates with operator[]

    This namespace contains an identical copy of all the functions given 
    in the \ref o2scl_cblas namespace, but perform matrix indexing
    with <tt>[][]</tt> rather than <tt>(,)</tt>. See \ref o2scl_cblas
    for the function listing and documentation.
 */
namespace o2scl_cblas_bracket {
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M[i][j]
#include <o2scl/cblas_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2

}

/** \brief The namespace for linear algebra classes and functions

    See \ref linalg_section for more complete information
    about linear algebra in \o2. 

    This namespace documentation is in the file 
    <tt>src/base/cblas.h</tt>
*/
namespace o2scl_linalg {
}

/** \brief The namespace for linear algebra classes and functions 
    with operator()

    This namespace contains an identical copy of all the functions given 
    in the \ref o2scl_cblas namespace, but perform matrix indexing
    with <tt>[][]</tt> rather than <tt>(,)</tt>. See \ref o2scl_linalg
    for the function listing and documentation.

    See \ref linalg_section for more complete information
    about linear algebra in \o2. 

    This namespace documentation is in the file 
    <tt>src/base/cblas.h</tt>
*/
namespace o2scl_linalg_bracket {
}


#endif
