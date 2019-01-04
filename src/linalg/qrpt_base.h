/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2019, Andrew W. Steiner

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
/* linalg/qrpt.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, 
 * Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
 * 02110-1301, USA.
 */
/** \file qrpt_base.h
    \brief File for QR decomposition and associated solver
*/

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Compute the QR decomposition of matrix \c A
   */
  template<class mat_t, class vec_t, class vec2_t>
    void QRPT_decomp(size_t M, size_t N, mat_t &A, vec_t &tau,
		     o2scl::permutation &p, int &signum, vec2_t &norm) {

    signum=1;
    p.init();

    for(size_t i=0;i<N;i++) {
      norm[i]=O2SCL_CBLAS_NAMESPACE::dnrm2_subcol(A,0,i,M);
    }
    
    size_t imax;
    if (M<N) imax=M;
    else imax=N;
    for(size_t i=0;i<imax;i++) {
      double max_norm=norm[i];
      size_t kmax=i;
      for(size_t j=i+1;j<N;j++) {
	double x=norm[j];
	if (x>max_norm) {
	  max_norm=x;
	  kmax=j;
	}
      }

      if (kmax!=i) {
	o2scl::matrix_swap_cols_double(M,A,i,kmax);
	p.swap(i,kmax);
	o2scl::vector_swap_double(norm,i,kmax);
	signum=-signum;
      }
      
      tau[i]=householder_transform_subcol(A,i,i,M);
      if (i+1<N) {
	householder_hm_subcol(A,i,i+1,M,N,A,i,i,tau[i]);
      }

      if (i+1<M) {
	for(size_t j=i+1;j<N;j++) {
	  double x=norm[j];
	  if (x>0.0) {
	    double y=0.0;
	    double temp=O2SCL_IX2(A,i,j)/x;
	    if (fabs(temp)>=1.0) {
	      y=0.0;
	    } else {
	      y=x*sqrt(1.0-temp*temp);
	    }
	    
	    double sqrt_dbl_eps=sqrt(std::numeric_limits<double>::epsilon());
	    
	    if (fabs(y/x)<sqrt(20.0)*sqrt_dbl_eps) {
	      y=O2SCL_CBLAS_NAMESPACE::dnrm2_subcol(A,i+1,j,M);
	    }
	    norm[j]=y;
	    
	  }
	}
      }

    }

    return;
  }
  
#ifdef DOXYGEN
}
#endif
