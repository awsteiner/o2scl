/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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
/* linalg/tridiag.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2002, 2004, 
 * 2007 Gerard Jungman, Brian Gough, David Necas
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
/** \file tridiag_base.h
    \brief File for solving tridiagonal systems
*/

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Allocation object for 2 arrays of equal size

      This class is used in solve_tridiag_nonsym().
   */
  class ubvector_2_mem {
  public:
    typedef boost::numeric::ublas::vector<double> ubvector;
    ubvector v1, v2;
    void resize(size_t n) {
      v1.resize(n);
      v2.resize(n);
    }
  };

  /** \brief Allocation object for 4 arrays of equal size

      This class is used in \ref solve_tridiag_sym() and 
      \ref solve_cyc_tridiag_nonsym().
  */
  class ubvector_4_mem {
  public:
    typedef boost::numeric::ublas::vector<double> ubvector;
    ubvector v1, v2, v3, v4;
    void resize(size_t n) {
      v1.resize(n);
      v2.resize(n);
      v3.resize(n);
      v4.resize(n);
    }
  };

  /** \brief Allocation object for 5 arrays of equal size

      This class is used in solve_cyc_tridiag_sym().
   */
  class ubvector_5_mem {
  public:
    typedef boost::numeric::ublas::vector<double> ubvector;
    ubvector v1, v2, v3, v4, v5;
    void resize(size_t n) {
      v1.resize(n);
      v2.resize(n);
      v3.resize(n);
      v4.resize(n);
      v5.resize(n);
    }
  };

  /** \brief Solve a symmetric tridiagonal linear system with 
      user-specified memory 

      This function solves the system \f$ A x = b \f$ where
      \f$ A \f$ is a matrix of the form
      \verbatim
      *
      *     diag[0]  offdiag[0]             0   .....
      *  offdiag[0]     diag[1]    offdiag[1]   .....
      *           0  offdiag[1]       diag[2]
      *           0           0    offdiag[2]   .....
      \endverbatim
      given the \c N diagonal elements in \c diag, \c N-1 diagonal
      elements in \c offdiag, and the \c N elements \c b from the RHS.

      See \ref EngelnMullges96 .
  */
  template<class vec_t, class vec2_t, class vec3_t, 
    class vec4_t, class mem_t, class mem_vec_t> 
    void solve_tridiag_sym(const vec_t &diag, const vec2_t &offdiag, 
			   const vec3_t &b, vec4_t &x, size_t N, mem_t &m) {
    
    mem_vec_t &gamma=m.v1;
    mem_vec_t &alpha=m.v2;
    mem_vec_t &c=m.v3;
    mem_vec_t &z=m.v4;
    
    size_t i, j;
    
    /* [GSL] Cholesky decomposition
       A = L.D.L^t
       lower_diag(L) = gamma
       diag(D) = alpha
    */
    alpha[0] = O2SCL_IX(diag,0);
    gamma[0] = O2SCL_IX(offdiag,0) / alpha[0];
    
    for (i = 1; i < N - 1; i++) {
      alpha[i] = O2SCL_IX(diag,i) - O2SCL_IX(offdiag,i-1) * gamma[i - 1];
      gamma[i] = O2SCL_IX(offdiag,i) / alpha[i];
    }
    
    if (N > 1) {
      alpha[N-1]=O2SCL_IX(diag,N-1)-O2SCL_IX(offdiag,N-2)*gamma[N-2];
    }
    
    // [GSL] update RHS 
    z[0] = b[0];
    for (i = 1; i < N; i++) {
      z[i] = O2SCL_IX(b,i) - gamma[i - 1] * z[i - 1];
    } for (i = 0; i < N; i++) {
      c[i] = z[i] / alpha[i];
    }
    
    // [GSL] backsubstitution 
    O2SCL_IX(x,N-1)=c[N-1];
    if (N >= 2) {
      for (i = N - 2, j = 0; j <= N - 2; j++, i--) {
	O2SCL_IX(x,i) = c[i] - gamma[i] * O2SCL_IX(x,i+1);
      }
    }
    
    return;
  }

  /** \brief Solve an asymmetric tridiagonal linear system with
      user-specified memory
    
      This function solves the system \f$ A x = b \f$ where
      \f$ A \f$ is a matrix of the form
      \verbatim
      * 
      *       diag[0]  abovediag[0]             0   .....
      *  belowdiag[0]       diag[1]  abovediag[1]   .....
      *             0  belowdiag[1]       diag[2]
      *             0             0  belowdiag[2]   .....
      \endverbatim
      This function uses plain Gauss elimination, not bothering
      with the zeroes.
      
      \future Offer an option to avoid throwing on divide by zero?
  */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t,
    class vec5_t, class mem_t, class mem_vec_t> 
    void solve_tridiag_nonsym(const vec_t &diag, const vec2_t &abovediag, 
			      const vec3_t &belowdiag, const vec4_t &rhs, 
			      vec5_t &x, size_t N, mem_t &m) {
    mem_vec_t &alpha=m.v1;
    mem_vec_t &z=m.v2;

    size_t i, j;

    /* [GSL] Bidiagonalization (eliminating belowdiag)
       & rhs update
       diag' = alpha
       rhs' = z
    */
    alpha[0] = O2SCL_IX(diag,0);
    z[0] = O2SCL_IX(rhs,0);

    for (i = 1; i < N; i++) {
      const double t = O2SCL_IX(belowdiag,i-1)/alpha[i-1];
      alpha[i] = O2SCL_IX(diag,i) - t*O2SCL_IX(abovediag,i-1);
      z[i] = O2SCL_IX(rhs,i) - t*z[i-1];
      if (alpha[i] == 0) {
	O2SCL_ERR2("Divide by zero in ",
		   "solve_tridiag_nonsym().",o2scl::exc_ezerodiv);
	return;
      }
    }
      
    // [GSL] backsubstitution
    O2SCL_IX(x,N-1) = z[N - 1]/alpha[N - 1];
    if (N >= 2) {
      for (i = N - 2, j = 0; j <= N - 2; j++, i--) {
	O2SCL_IX(x,i) = (z[i] - O2SCL_IX(abovediag,i) * 
			 O2SCL_IX(x,i+1))/alpha[i];
      }
    }
      
    return;
  }

  /** \brief Solve a symmetric cyclic tridiagonal linear system with
      user specified memory
    
      This function solves the system \f$ A x = b \f$ where
      \f$ A \f$ is a matrix of the form
      \verbatim
      *
      *      diag[0]  offdiag[0]             0   .....  offdiag[N-1]
      *   offdiag[0]     diag[1]    offdiag[1]   .....
      *            0  offdiag[1]       diag[2]
      *            0           0    offdiag[2]   .....
      *          ...         ...
      * offdiag[N-1]         ...
      \endverbatim

      See \ref EngelnMullges96 .
  */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t,
    class mem_t, class mem_vec_t> 
    void solve_cyc_tridiag_sym(const vec_t &diag, const vec2_t &offdiag, 
			       const vec3_t &b, vec4_t &x, size_t N,
			       mem_t &m) {

      mem_vec_t &delta=m.v1;
      mem_vec_t &gamma=m.v2;
      mem_vec_t &alpha=m.v3;
      mem_vec_t &c=m.v4;
      mem_vec_t &z=m.v5;

      size_t i, j;
      double sum = 0.0;
      
      // [GSL] factor

      if (N == 1)  {
	x[0] = b[0] / O2SCL_IX(diag,0);
	return;
      }

      alpha[0] = O2SCL_IX(diag,0);
      gamma[0] = O2SCL_IX(offdiag,0) / alpha[0];
      delta[0] = O2SCL_IX(offdiag,N-1) / alpha[0];
      
      for (i = 1; i < N - 2; i++) {
	alpha[i] = O2SCL_IX(diag,i) - O2SCL_IX(offdiag,i-1) * gamma[i - 1];
	gamma[i] = O2SCL_IX(offdiag,i) / alpha[i];
	delta[i] = -delta[i - 1] * O2SCL_IX(offdiag,i-1) / alpha[i];
      }

      for (i = 0; i < N - 2; i++) {
	sum += alpha[i] * delta[i] * delta[i];
      }

      alpha[N - 2] = diag[ (N - 2)] - O2SCL_IX(offdiag,N-3) * gamma[N - 3];

      gamma[N - 2] = (offdiag[(N - 2)] - offdiag[(N - 3)] * 
		      delta[N - 3]) / alpha[N - 2];
      
      alpha[N - 1] = diag[(N - 1)] - sum - alpha[(N - 2)] * 
	gamma[N - 2] * gamma[N - 2];
      
      // [GSL] update

      z[0] = b[0];
      for (i = 1; i < N - 1; i++) {
	z[i] = O2SCL_IX(b,i) - z[i - 1] * gamma[i - 1];
      }
      sum = 0.0;
      for (i = 0; i < N - 2; i++) {
	sum += delta[i] * z[i];
      }
      z[N - 1] = b[(N - 1)] - sum - gamma[N - 2] * z[N - 2];
      for (i = 0; i < N; i++) {
	c[i] = z[i] / alpha[i];
      }

      // [GSL] backsubstitution 
      O2SCL_IX(x,N-1) = c[N - 1];
      x[(N - 2)] = c[N - 2] - gamma[N - 2] * O2SCL_IX(x,N-1);
      if (N >= 3) {
	for (i = N - 3, j = 0; j <= N - 3; j++, i--) {
	  O2SCL_IX(x,i) = c[i] - gamma[i] * x[(i + 1)] - 
	    delta[i] * O2SCL_IX(x,N-1);
	}
      }

      return;
    }

  /** \brief Solve an asymmetric cyclic tridiagonal linear system
      with user-specified memory
    
      This function solves the system \f$ A x = b \f$ where
      \f$ A \f$ is a matrix of the form
      \verbatim
      *
      *        diag[0]  abovediag[0]             0   .....  belowdiag[N-1]
      *   belowdiag[0]       diag[1]  abovediag[1]   .....
      *              0  belowdiag[1]       diag[2]
      *              0             0  belowdiag[2]   .....
      *            ...           ...
      * abovediag[N-1]           ...
      \endverbatim

      This function solves the following system without the corner
      elements and then use Sherman-Morrison formula to compensate for
      them.

      \comment
      Note that the three FIXME!!! entries are from the GSL original. 
      \endcomment

      \future Offer an option to avoid throwing on divide by zero?
  */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t, 
    class vec5_t, class mem_t, class mem_vec_t> 
    void solve_cyc_tridiag_nonsym(const vec_t &diag, const vec2_t &abovediag, 
				  const vec3_t &belowdiag, const vec4_t &rhs, 
				  vec5_t &x, size_t N, mem_t &m) {

    double beta;
    mem_vec_t &alpha=m.v1;
    mem_vec_t &zb=m.v2;
    mem_vec_t &zu=m.v3;
    mem_vec_t &w=m.v4;

    /* [GSL] Bidiagonalization (eliminating belowdiag)
       & rhs update
       diag' = alpha
       rhs' = zb
       rhs' for Aq=u is zu
    */
    zb[0] = O2SCL_IX(rhs,0);
    if (O2SCL_IX(diag,0) != 0) {
      beta = -O2SCL_IX(diag,0); 
    } else {
      beta = 1;
    }
    const double q = 1 - O2SCL_IX(abovediag,0)*O2SCL_IX(belowdiag,0)/
      (O2SCL_IX(diag,0)*diag[1]);
    if (fabs(q/beta) > 0.5 && fabs(q/beta) < 2) {
      beta *= (fabs(q/beta) < 1) ? 0.5 : 2;
    }
    zu[0] = beta;
    alpha[0] = O2SCL_IX(diag,0) - beta;

    {
      size_t i;
      for (i = 1; i+1 < N; i++) {
	const double t = O2SCL_IX(belowdiag,i-1)/alpha[i-1];
	alpha[i] = O2SCL_IX(diag,i) - t*O2SCL_IX(abovediag,i-1);
	zb[i] = O2SCL_IX(rhs,i) - t*zb[i-1];
	zu[i] = -t*zu[i-1];
	// [GSL] FIXME!!!
	if (alpha[i] == 0) {
	  O2SCL_ERR2("Divide by zero (1) in ",
		     "solve_cyc_tridiag_nonsym().",o2scl::exc_ezerodiv);
	}
      }
    }

    {
      const size_t i = N-1;
      const double t = O2SCL_IX(belowdiag,i-1)/alpha[i-1];
      alpha[i]=O2SCL_IX(diag,i)-O2SCL_IX(abovediag,i)*
	O2SCL_IX(belowdiag,i)/beta-
	t*O2SCL_IX(abovediag,i-1);
      zb[i] = O2SCL_IX(rhs,i) - t*zb[i-1];
      zu[i] = O2SCL_IX(abovediag,i) - t*zu[i-1];

      // [GSL] FIXME!!!
      if (alpha[i] == 0) {
	O2SCL_ERR2("Divide by zero (2) in ",
		   "solve_cyc_tridiag_nonsym().",o2scl::exc_ezerodiv);
      }
    }

    {
      // [GSL] backsubstitution
      size_t i, j;
      w[N-1] = zu[N-1]/alpha[N-1];
      O2SCL_IX(x,N-1) = zb[N-1]/alpha[N-1];
      for (i = N - 2, j = 0; j <= N - 2; j++, i--) {
	w[i] = (zu[i] - O2SCL_IX(abovediag,i) * w[i+1])/alpha[i];
	O2SCL_IX(x,i) = (zb[i] - O2SCL_IX(abovediag,i) * 
			 O2SCL_IX(x,i + 1))/alpha[i];
      }
    }
      
    // [GSL] Sherman-Morrison
    const double vw = w[0] + O2SCL_IX(belowdiag,N-1)/beta * w[N-1];
    const double vx = O2SCL_IX(x,0) + 
      O2SCL_IX(belowdiag,N-1)/beta * O2SCL_IX(x,N-1);

    // [GSL] FIXME!!!
    if (vw + 1 == 0) {
      O2SCL_ERR2("Divide by zero (3) in ",
		 "solve_cyc_tridiag_nonsym().",o2scl::exc_ezerodiv);
    }

    {
      size_t i;
      for (i = 0; i < N; i++)
	O2SCL_IX(x,i) -= vx/(1 + vw)*w[i];
    }

    return;
  }

  /** \brief Solve a symmetric tridiagonal linear system 
   */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t> 
    void solve_tridiag_sym(const vec_t &diag, const vec2_t &offdiag, 
			   const vec3_t &b, vec4_t &x, size_t N) {
    typedef boost::numeric::ublas::vector<double> ubvector;
    ubvector_4_mem v4m;
    v4m.resize(N);
    solve_tridiag_sym<vec_t,vec2_t,vec3_t,vec4_t,ubvector_4_mem,
      ubvector>(diag,offdiag,b,x,N,v4m);
    return;
  }
  
  /** \brief Solve an asymmetric tridiagonal linear system 
   */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t,
    class vec5_t> 
    void solve_tridiag_nonsym(const vec_t &diag, const vec2_t &abovediag, 
			      const vec3_t &belowdiag, const vec4_t &rhs, 
			      vec5_t &x, size_t N) {
    typedef boost::numeric::ublas::vector<double> ubvector;
      ubvector_2_mem v2m;
      v2m.resize(N);
      solve_tridiag_nonsym<vec_t,vec2_t,vec3_t,vec4_t,vec5_t,ubvector_2_mem,
      ubvector>(diag,abovediag,belowdiag,rhs,x,N,v2m);
      return;
    }
  
  /** \brief Solve a symmetric cyclic tridiagonal linear system 
   */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t> 
    void solve_cyc_tridiag_sym(const vec_t &diag, const vec2_t &offdiag, 
			       const vec3_t &b, vec4_t &x, size_t N) {
    typedef boost::numeric::ublas::vector<double> ubvector;
    ubvector_5_mem v5m;
    v5m.resize(N);
    solve_cyc_tridiag_sym<vec_t,vec2_t,vec3_t,vec4_t,ubvector_5_mem,
      ubvector>(diag,offdiag,b,x,N,v5m);
    return;
  }
  
  /** \brief Solve an asymmetric cyclic tridiagonal linear system
   */
  template<class vec_t, class vec2_t, class vec3_t, class vec4_t,
    class vec5_t> 
    void solve_cyc_tridiag_nonsym(const vec_t &diag, const vec2_t &abovediag, 
				  const vec3_t &belowdiag, const vec4_t &rhs, 
				  vec5_t &x, size_t N) {
    typedef boost::numeric::ublas::vector<double> ubvector;
      ubvector_4_mem v4m;
      v4m.resize(N);
      solve_cyc_tridiag_nonsym<vec_t,vec2_t,vec3_t,vec4_t,vec5_t,
	ubvector_4_mem,ubvector>(diag,abovediag,belowdiag,rhs,x,N,v4m);
      return;
    }
  
#ifdef DOXYGEN
}
#endif
