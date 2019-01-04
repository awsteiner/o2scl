/*
  -------------------------------------------------------------------
  
  Copyright (C) 2010-2019, Andrew W. Steiner
  
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
/** \file svd.h
    \brief Header wrapper for \ref svd_base.h
*/
#ifndef O2SCL_SVD_H
#define O2SCL_SVD_H

// For gsl_isnan()
#include <gsl/gsl_sys.h>
// For GSL_SQRT_DBL_MAX, etc.
#include <gsl/gsl_machine.h>

#include <o2scl/err_hnd.h>
#include <o2scl/permutation.h>
#include <o2scl/cblas.h>
#include <o2scl/vector.h>
#include <o2scl/bidiag.h>
#include <o2scl/svdstep.h>

namespace o2scl_linalg {

#define O2SCL_CBLAS_NAMESPACE o2scl_cblas

#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M(i,j)
#include <o2scl/svd_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2

#undef O2SCL_CBLAS_NAMESPACE

}

namespace o2scl_linalg_bracket {
  
#define O2SCL_CBLAS_NAMESPACE o2scl_cblas_bracket

#define O2SCL_IX(V,i) V[i]
#define O2SCL_IX2(M,i,j) M[i][j]
#include <o2scl/svd_base.h>  
#undef O2SCL_IX
#undef O2SCL_IX2
  
#undef O2SCL_CBLAS_NAMESPACE

}

#endif
