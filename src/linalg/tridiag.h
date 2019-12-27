/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
/** \file tridiag.h
    \brief Header wrapper for \ref tridiag_base.h
*/
#ifndef O2SCL_TRIDIAG_H
#define O2SCL_TRIDIAG_H

#include <o2scl/householder.h>
#include <o2scl/givens.h>

namespace o2scl_linalg {
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M(i,j)
#include <o2scl/tridiag_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2
  
}

namespace o2scl_linalg_bracket {
  
#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M[i][j]
#include <o2scl/tridiag_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2

}

#endif
