/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
/** \file givens.h
    \brief Header wrapper for \ref givens_base.h
*/
#ifndef O2SCL_GIVENS_H
#define O2SCL_GIVENS_H

#include <gsl/gsl_math.h>

#include <o2scl/err_hnd.h>
#include <o2scl/permutation.h>
#include <o2scl/cblas.h>

namespace o2scl_linalg_bracket {
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M[i][j]
#include <o2scl/givens_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2

}

namespace o2scl_linalg {
  
  /** \brief Create a Givens rotation matrix

      Given values \c a and \c b, create entries \c c and 
      \c s of a matrix for which
      \f[
      \left[ \begin{array}{cc} c & -s \\ s & c \end{array} \right]
      \left[ \begin{array}{c} a \\ b \end{array} \right] =
      \left[ \begin{array}{c} r \\ 0 \end{array} \right]
      \f]
      
   */
  void create_givens(const double a, const double b, double &c, 
                     double &s);
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M(i,j)
#include <o2scl/givens_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2
  
}

#endif
