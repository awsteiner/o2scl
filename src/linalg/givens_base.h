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
/* linalg/givens.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough
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
/** \file givens_base.h
    \brief File for Givens rotations

    \todo Make sure create_givens() in givens.h is documented.
*/

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Apply a rotation to matrices from the QR decomposition

      This performs \f$ Q \rightarrow Q G \f$ and 
      \f$ R \rightarrow G^{T} R \f$.
   */
  template<class mat1_t, class mat2_t>
    void apply_givens_qr(size_t M, size_t N, mat1_t &Q, mat2_t &R, 
			 size_t i, size_t j, double c, double s) {

    for(size_t k=0;k<M;k++) {
      double qki=O2SCL_IX2(Q,k,i);
      double qkj=O2SCL_IX2(Q,k,j);
      O2SCL_IX2(Q,k,i)=qki*c-qkj*s;
      O2SCL_IX2(Q,k,j)=qki*s+qkj*c;
    }
    size_t kstart;
    if (i<j) kstart=i;
    else kstart=j;
    for(size_t k=kstart;k<N;k++) {
      double rik=O2SCL_IX2(R,i,k);
      double rjk=O2SCL_IX2(R,j,k);
      O2SCL_IX2(R,i,k)=c*rik-s*rjk;
      O2SCL_IX2(R,j,k)=s*rik+c*rjk;
    }

    return;
  }

  /** \brief Apply a rotation to matrices from the LQ decomposition

      This performs \f$ Q \rightarrow Q G \f$ and 
      \f$ L \rightarrow L G^{T} \f$.
   */
  template<class mat1_t, class mat2_t>
    void apply_givens_lq(size_t M, size_t N, mat1_t &Q, mat2_t &L, 
			 size_t i, size_t j, double c, double s) {

    for(size_t k=0;k<M;k++) {
      double qik=O2SCL_IX2(Q,i,k);
      double qjk=O2SCL_IX2(Q,j,k);
      O2SCL_IX2(Q,i,k)=qik*c-qjk*s;
      O2SCL_IX2(Q,j,k)=qik*s+qjk*c;
    }
    size_t kstart;
    if (i<j) kstart=i;
    else kstart=j;
    for(size_t k=kstart;k<N;k++) {
      double lki=O2SCL_IX2(L,k,i);
      double lkj=O2SCL_IX2(L,k,j);
      O2SCL_IX2(L,k,i)=c*lki-s*lkj;
      O2SCL_IX2(L,k,j)=s*lki+c*lkj;
    }

    return;
  }
  
  /// Apply a rotation to a vector, \f$ v \rightarrow G^{T} v \f$
  template<class vec_t>
    void apply_givens_vec(vec_t &v, size_t i, size_t j,
			  double c, double s) {
    double vi=O2SCL_IX(v,i);
    double vj=O2SCL_IX(v,j);
    O2SCL_IX(v,i)=c*vi-s*vj;
    O2SCL_IX(v,j)=s*vi+c*vj;
    return;
  }

#ifdef DOXYGEN
}
#endif
