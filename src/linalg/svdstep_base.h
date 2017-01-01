/*
  -------------------------------------------------------------------

  Copyright (C) 2008-2017, Andrew W. Steiner

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
/* linalg/svdstep.c
 * 
 * Copyright (C) 2007 Brian Gough
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
/** \file svdstep_base.h
    \brief File for SVD decomposition
*/

#ifdef DOXYGEN
namespace o2scl_linalg {
#endif

  /** \brief Zero out small elements in \c f according to the 
      scales set in \c d
      
      The parameter \c N is the size of \c d.

      Used in \ref SV_decomp().
  */
  template<class vec_t, class vec2_t>
    void chop_small_elements(size_t N, vec_t &d, vec2_t &f) {
    
    double d_i=O2SCL_IX(d,0);
    
    for (size_t i=0;i<N-1;i++) {
      
      double f_i=O2SCL_IX(f,i);
      double d_ip1=O2SCL_IX(d,i+1);
      
      double dbl_eps=std::numeric_limits<double>::epsilon();

      if (fabs (f_i)<dbl_eps*(fabs(d_i)+fabs(d_ip1))) {
	O2SCL_IX(f,i)=0.0;
      }
      d_i=d_ip1;
    }
    
    return;
  }
  
  /** \brief Desc
      
      The parameter \c n is the size of the vector \c d.

      Used in \ref qrstep() and \ref qrstep_sub().
  */
  template<class vec_t, class vec2_t>
    double trailing_eigenvalue(size_t n, const vec_t &d, const vec2_t &f) {
    
    double da=O2SCL_IX(d,n-2);
    double db=O2SCL_IX(d,n-1);
    double fa=(n>2) ? O2SCL_IX(f,n-3) : 0.0;
    double fb=O2SCL_IX(f,n-2);
    
    double ta=da*da+fa*fa;
    double tb=db*db+fb*fb;
    double tab=da*fb;
    
    double dt=(ta-tb)/2.0;
    
    double S=ta+tb;
    double da2=da*da,db2=db*db;
    double fa2=fa*fa,fb2=fb*fb;
    double P=(da2*db2)+(fa2*db2)+(fa2*fb2);
    double D=hypot(dt,tab);
    double r1=S/2+D;
    
    double mu;
    if (dt>=0) {
      // [GSL] tb < ta, choose smaller root
      mu=(r1>0) ? P / r1 : 0.0;
    } else {
      // [GSL] tb > ta, choose larger root
      mu=r1;
    }
    return mu;
  }
  
  /** \brief Desc

      Used in \ref svd2() and \ref svd2_sub().
   */
  void create_schur(double d0, double f0, double d1, double &c, 
		    double &s) {
    
    double apq=2.0*d0*f0;
    
    if (d0 == 0 || f0 == 0) {
      c=1.0;
      s=0.0;
      return;
    }
    
    // [GSL] Check if we need to rescale to avoid underflow/overflow 
    if (fabs(d0) < GSL_SQRT_DBL_MIN || fabs(d0) > GSL_SQRT_DBL_MAX
	|| fabs(f0) < GSL_SQRT_DBL_MIN || fabs(f0) > GSL_SQRT_DBL_MAX
	|| fabs(d1) < GSL_SQRT_DBL_MIN || fabs(d1) > GSL_SQRT_DBL_MAX) {
      
      double scale;
      int d0_exp,f0_exp;
      frexp(d0,&d0_exp);
      frexp(f0,&f0_exp);
      // [GSL] Bring |d0*f0| into the range GSL_DBL_MIN to GSL_DBL_MAX
      scale=ldexp(1.0,-(d0_exp+f0_exp)/4);
      d0*=scale;
      f0*=scale;
      d1*=scale;
      apq=2.0*d0*f0;
    }
    
    if (apq != 0.0) {
      double t;
      double tau=(f0*f0+(d1+d0)*(d1-d0))/apq;
      
      if (tau >= 0.0) {
	t=1.0/(tau+hypot(1.0,tau));
      } else {
	t=-1.0/(-tau+hypot(1.0,tau));
      }
      
      c=1.0/hypot(1.0,t);
      s=t*(c);
    } else {
      c=1.0;
      s=0.0;
    }
    return;
  }
  
  /** \brief 2-variable SVD
      
      The parameter \c M is the number of rows in \c U and \c N 
      is the number of rows in \c V. Both U and V assumed to have
      two columns.

      Used in \ref qrstep().
  */
  template<class vec_t, class vec2_t, class mat_t, class mat2_t>
    void svd2(size_t M, size_t N, vec_t &d, vec2_t &f, mat_t &U, mat2_t &V) {
    
    size_t i;
    double c,s,a11,a12,a21,a22;
    
    double d0=O2SCL_IX(d,0);
    double f0=O2SCL_IX(f,0);
    
    double d1=O2SCL_IX(d,1);
    
    if (d0 == 0.0) {
      
      // [GSL] Eliminate off-diagonal element in [0,f0;0,d1] to 
      // make [d,0;0,0]
      o2scl_linalg::create_givens(f0,d1,c,s);
      
      // [GSL] compute B <= G^T B X, where X=[0,1;1,0]

      O2SCL_IX(d,0)=c*f0-s*d1;
      O2SCL_IX(f,0)=s*f0+c*d1;
      O2SCL_IX(d,1)=0.0;
      
      // [GSL] Compute U <= U G
      
      for (size_t i=0;i<M;i++) {
	
	double Uip=O2SCL_IX2(U,i,0);
	double Uiq=O2SCL_IX2(U,i,1);
	O2SCL_IX2(U,i,0)=c*Uip-s*Uiq;
	O2SCL_IX2(U,i,1)=s*Uip+c*Uiq;
      }
      
      // [GSL] Compute V <= V X 

      double temp;
      for(size_t ik=0;ik<N;ik++) {
	temp=O2SCL_IX2(V,ik,0);
	O2SCL_IX2(V,ik,0)=O2SCL_IX2(V,ik,1);
	O2SCL_IX2(V,ik,1)=temp;
      }

      return;

    } else if (d1 == 0.0) {

      // [GSL] Eliminate off-diagonal element in [d0,f0;0,0]

      o2scl_linalg::create_givens(d0,f0,c,s);

      // [GSL] compute B <= B G

      O2SCL_IX(d,0)=d0*c-f0*s;
      O2SCL_IX(f,0)=0.0;

      // [GSL] Compute V <= V G 
      
      for (size_t i=0;i<N;i++) {
	double Vip=O2SCL_IX2(V,i,0);
	double Viq=O2SCL_IX2(V,i,1);
	O2SCL_IX2(V,i,0)=c*Vip-s*Viq;
	O2SCL_IX2(V,i,1)=s*Vip+c*Viq;
      }
      
      return;

    } else {

      // [GSL] Make columns orthogonal, A = [d0, f0; 0, d1] * G 

      create_schur(d0,f0,d1,c,s);

      // [GSL] compute B <= B G 
      
      a11=c*d0-s*f0;
      a21=-s*d1;
      a12=s*d0+c*f0;
      a22=c*d1;
      
      // [GSL] Compute V <= V G 
      
      for (size_t i=0;i<N;i++) {
	
	double Vip=O2SCL_IX2(V,i,0);
	double Viq=O2SCL_IX2(V,i,1);
	O2SCL_IX2(V,i,0)=c*Vip-s*Viq;
	O2SCL_IX2(V,i,1)=s*Vip+c*Viq;
      }
      
      // [GSL] Eliminate off-diagonal elements, bring column with largest
      // norm to first column
      
      if (hypot(a11,a21) < hypot(a12,a22)) {
	
	double t1,t2;
	
	// [GSL] B <= B X
	
	t1=a11; a11=a12; a12=t1;
	t2=a21; a21=a22; a22=t2;
	
	// [GSL] V <= V X
	
	double temp;
	for(size_t ik=0;ik<N;ik++) {
	  temp=O2SCL_IX2(V,ik,0);
	  O2SCL_IX2(V,ik,0)=O2SCL_IX2(V,ik,1);
	  O2SCL_IX2(V,ik,1)=temp;
	}
      } 
      
      o2scl_linalg::create_givens(a11,a21,c,s);
      
      // [GSL] compute B <= G^T B
      
      O2SCL_IX(d,0)=c*a11-s*a21;
      O2SCL_IX(f,0)=c*a12-s*a22;
      O2SCL_IX(d,1)=s*a12+c*a22;
      
      // [GSL] Compute U <= U G
      
      for (size_t i=0;i<M;i++) {
	double Uip=O2SCL_IX2(U,i,0);
	double Uiq=O2SCL_IX2(U,i,1);
	O2SCL_IX2(U,i,0)=c*Uip-s*Uiq;
	O2SCL_IX2(U,i,1)=s*Uip+c*Uiq;
      }
      
      return;
    }
  }

  /** \brief Shifted 2-variable SVD
      
      The parameter \c M is the number of rows in \c U and \c N 
      is the number of rows in \c V. Both U and V assumed to have
      two columns.

      Used in \ref qrstep_sub().
   */
  template<class vec_t, class vec2_t, class mat_t, class mat2_t>
    void svd2_sub(size_t M, size_t N, vec_t &d, vec2_t &f, mat_t &U, 
		  mat2_t &V, size_t a) {
    
    size_t i;
    double c,s,a11,a12,a21,a22;
    
    double d0=O2SCL_IX(d,a);
    double f0=O2SCL_IX(f,a);
    
    double d1=O2SCL_IX(d,a+1);
    
    if (d0 == 0.0) {
      
      // [GSL] Eliminate off-diagonal element in [0,f0;0,d1] to 
      // make [d,0;0,0]
      o2scl_linalg::create_givens(f0,d1,c,s);
      
      // [GSL] compute B <= G^T B X, where X=[0,1;1,0]

      O2SCL_IX(d,a)=c*f0-s*d1;
      O2SCL_IX(f,a)=s*f0+c*d1;
      O2SCL_IX(d,a+1)=0.0;
      
      // [GSL] Compute U <= U G
      
      for (size_t i=0;i<M;i++) {
	
	double Uip=O2SCL_IX2(U,i,a);
	double Uiq=O2SCL_IX2(U,i,a+1);
	O2SCL_IX2(U,i,a)=c*Uip-s*Uiq;
	O2SCL_IX2(U,i,a+1)=s*Uip+c*Uiq;
      }
      
      // [GSL] Compute V <= V X 

      double temp;
      for(size_t ik=0;ik<N;ik++) {
	temp=O2SCL_IX2(V,ik,a);
	O2SCL_IX2(V,ik,a)=O2SCL_IX2(V,ik,a+1);
	O2SCL_IX2(V,ik,a+1)=temp;
      }

      return;

    } else if (d1 == 0.0) {

      // [GSL] Eliminate off-diagonal element in [d0,f0;0,0]

      o2scl_linalg::create_givens(d0,f0,c,s);

      // [GSL] compute B <= B G

      O2SCL_IX(d,a)=d0*c-f0*s;
      O2SCL_IX(f,a)=0.0;

      // [GSL] Compute V <= V G 
      
      for (size_t i=0;i<N;i++) {
	double Vip=O2SCL_IX2(V,i,a);
	double Viq=O2SCL_IX2(V,i,a+1);
	O2SCL_IX2(V,i,a)=c*Vip-s*Viq;
	O2SCL_IX2(V,i,a+1)=s*Vip+c*Viq;
      }
      
      return;

    } else {

      // [GSL] Make columns orthogonal, A = [d0, f0; 0, d1] * G 

      create_schur(d0,f0,d1,c,s);

      // [GSL] compute B <= B G 
      
      a11=c*d0-s*f0;
      a21=-s*d1;
      a12=s*d0+c*f0;
      a22=c*d1;
      
      // [GSL] Compute V <= V G 
      
      for (size_t i=0;i<N;i++) {
	
	double Vip=O2SCL_IX2(V,i,a);
	double Viq=O2SCL_IX2(V,i,a+1);
	O2SCL_IX2(V,i,a)=c*Vip-s*Viq;
	O2SCL_IX2(V,i,a+1)=s*Vip+c*Viq;
      }
      
      // [GSL] Eliminate off-diagonal elements, bring column with largest
      // norm to first column
      
      if (hypot(a11,a21)<hypot(a12,a22)) {
	
	double t1, t2;
	
	// [GSL] B <= B X
	
	t1=a11; a11=a12; a12=t1;
	t2=a21; a21=a22; a22=t2;
	
	// [GSL] V <= V X
	
	double temp;
	for(size_t ik=0;ik<N;ik++) {
	  temp=O2SCL_IX2(V,ik,a);
	  O2SCL_IX2(V,ik,a)=O2SCL_IX2(V,ik,a+1);
	  O2SCL_IX2(V,ik,a+1)=temp;
	}
      } 
      
      o2scl_linalg::create_givens(a11,a21,c,s);
      
      // [GSL] compute B <= G^T B
      
      O2SCL_IX(d,a)=c*a11-s*a21;
      O2SCL_IX(f,a)=c*a12-s*a22;
      O2SCL_IX(d,a+1)=s*a12+c*a22;
      
      // [GSL] Compute U <= U G
      
      for (size_t i=0;i<M;i++) {
	double Uip=O2SCL_IX2(U,i,a);
	double Uiq=O2SCL_IX2(U,i,a+1);
	O2SCL_IX2(U,i,a)=c*Uip-s*Uiq;
	O2SCL_IX2(U,i,a+1)=s*Uip+c*Uiq;
      }
      
      return;
    }
  }
  
  /** \brief Desc
      
      The vector \c d should be of size <tt>n</tt>, the vector \c f
      should be of size <tt>n</tt>, and the matrix U should be of size
      <tt>(M,n)</tt>

      Used in \ref qrstep() and \ref qrstep_sub().
  */
  template<class vec_t, class vec2_t, class mat_t>
    void chase_out_intermediate_zero(size_t M, size_t n, vec_t &d,
				     vec2_t &f, mat_t &U, size_t k0) {
    
    double c, s;
    
    double x=O2SCL_IX(f,k0);
    double y=O2SCL_IX(d,k0+1);
    
    for (size_t k=k0;k<n-1;k++) {
      
      o2scl_linalg::create_givens(y,-x,c,s);
      
      // [GSL] Compute U <= U G 
      for (size_t i=0; i < M; i++) {
	double Uip=O2SCL_IX2(U,i,k0);
	double Uiq=O2SCL_IX2(U,i,k+1);
	//std::cout << "Uip,Uiq: " << Uip << " " << Uiq << std::endl;
	O2SCL_IX2(U,i,k0)=c*Uip-s*Uiq;
	O2SCL_IX2(U,i,k+1)=s*Uip+c*Uiq;
      }
      
      // [GSL] compute B <= G^T B
      
      O2SCL_IX(d,k+1)=s*x+c*y;
      
      if (k == k0) {
        O2SCL_IX(f,k)=c*x-s*y;
      }
      
      if (k<n-2) {
	double z=O2SCL_IX(f,k+1);
	O2SCL_IX(f,k+1)=c*z; 
	
	x=-s*z;
	y=O2SCL_IX(d,k+2); 
      }
    }

    return;
  }

  /** \brief Desc
      
      The vector \c d should be of size <tt>n</tt>, the vector \c f
      should be of size <tt>n</tt>, and the matrix V should be of size
      <tt>(N,n)</tt>

      Used in \ref qrstep().
  */
  template<class vec_t, class vec2_t, class mat_t>
    void chase_out_trailing_zero(size_t N, size_t n, vec_t &d, 
				vec2_t &f, mat_t &V) {
    
    double c, s;

    double x=O2SCL_IX(d,n-2);
    double y=O2SCL_IX(f,n-2);
    
    for (size_t k=n-1;k-- > 0;) {
      
      o2scl_linalg::create_givens(x,y,c,s);
      
      // [GSL] Compute V <= V G where G = [c, s ; -s, c] 
      
      for (size_t i=0;i<N;i++) {
	double Vip=O2SCL_IX2(V,i,k);
	double Viq=O2SCL_IX2(V,i,n-1);
	O2SCL_IX2(V,i,k)=c*Vip-s*Viq;
	O2SCL_IX2(V,i,n-1)=s*Vip+c*Viq;
      }
      
      // [GSL] compute B <= B G 
      
      O2SCL_IX(d,k)=c*x-s*y;
      
      if (k==n-2) {
        O2SCL_IX(f,k)=s*x+c*y;
      }
      
      if (k>0) {
	double z=O2SCL_IX(f,k-1);
	O2SCL_IX(f,k-1)=c*z; 
	x=O2SCL_IX(d,k-1); 
	y=s*z;
      }
    }

    return;
  }

  /** \brief Desc

      The vector \c d should be of size <tt>n</tt>, the vector \c f
      should be of size <tt>n</tt>, and the matrix V should be of size
      <tt>(N,n)</tt>

      Used in \ref qrstep_sub().
  */
  template<class vec_t, class vec2_t, class mat_t>
    void chase_out_trailing_zero_sub(size_t N, size_t n, vec_t &d, 
				     vec2_t &f, mat_t &V, size_t a) {
    
    double c, s;

    double x=O2SCL_IX(d,n-2);
    double y=O2SCL_IX(f,n-2);
    
    for (int k=((int)n)-1;k>=((int)a);k--) {
      
      o2scl_linalg::create_givens(x,y,c,s);
      
      // [GSL] Compute V <= V G where G = [c, s ; -s, c] 
      
      for (size_t i=0;i<N;i++) {
	double Vip=O2SCL_IX2(V,i,k);
	double Viq=O2SCL_IX2(V,i,n-1);
	O2SCL_IX2(V,i,k)=c*Vip-s*Viq;
	O2SCL_IX2(V,i,n-1)=s*Vip+c*Viq;
      }
      
      // [GSL] compute B <= B G 
      
      O2SCL_IX(d,k)=c*x-s*y;
      
      if (k==((int)n)-2) {
        O2SCL_IX(f,k)=s*x+c*y;
      }
      
      if (k>((int)a)) {
	double z=O2SCL_IX(f,k-1);
	O2SCL_IX(f,k-1)=c*z; 
	x=O2SCL_IX(d,k-1); 
	y=s*z;
      }
    }

    return;
  }
  
  /** \brief Desc
      
      The vector \c d should be of size \c n, the vector \c f should
      be of size \c n, the matrix U should be of size <tt>(M,N)</tt>,
      and the matrix \c V should be of size <tt>(N,N)</tt>.
  */
  template<class vec_t, class vec2_t, class mat_t, class mat2_t>
    void qrstep(size_t M, size_t N, size_t n, 
		vec_t &d, vec2_t &f, mat_t &U, mat2_t &V) {
    
    double y, z;
    double ak, bk, zk, ap, bp, aq, bq;
    size_t i, k;
    
    if (n==1) {
      // [GSL] shouldn't happen
      return;
    }
    
    // [GSL] Compute 2x2 svd directly
    
    if (n==2) {
      svd2(M,N,d,f,U,V);
      return;
    }
    
    // [GSL] Chase out any zeroes on the diagonal
    
    for (i=0;i<n-1;i++) {
      double d_i=O2SCL_IX(d,i);
      if (d_i==0.0) {
	chase_out_intermediate_zero(M,n,d,f,U,i);
	return;
      }
    }
    
    // [GSL] Chase out any zero at the end of the diagonal
    double d_nm1=O2SCL_IX(d,n-1);
    
    if (d_nm1==0.0) {
      chase_out_trailing_zero(N,n,d,f,V);
      return;
    }

    // [GSL] Apply QR reduction steps to the diagonal and offdiagonal
    
    double d0=O2SCL_IX(d,0);
    double f0=O2SCL_IX(f,0);
    
    double d1=O2SCL_IX(d,1);
    double f1=O2SCL_IX(f,1);
    
    double mu=trailing_eigenvalue(n,d,f);
    y=d0*d0-mu;
    z=d0*f0;
    
    // [GSL] Set up the recurrence for Givens rotations on a bidiagonal 
    // matrix
    
    ak=0;
    bk=0;
    
    ap=d0;
    bp=f0;
    
    aq=d1;
    bq=f1;
    
    for (k=0; k < n-1; k++) {

      double c, s;
      o2scl_linalg::create_givens(y,z,c,s);

      // [GSL] Compute V <= V G
      
      for (i=0; i < N; i++) {
	double Vip=O2SCL_IX2(V,i,k);
	double Viq=O2SCL_IX2(V,i,k+1);
	O2SCL_IX2(V,i,k)=c*Vip-s*Viq;
	O2SCL_IX2(V,i,k+1)=s*Vip+c*Viq;
      }
      
      // [GSL] compute B <= B G
      
      {
        double bk1=c*bk-s*z;

        double ap1=c*ap-s*bp;
        double bp1=s*ap+c*bp;
        double zp1=-s*aq;

        double aq1=c*aq;

        if (k > 0) O2SCL_IX(f,k-1)=bk1;
	
        ak=ap1;
        bk=bp1;
        zk=zp1;
	
        ap=aq1;
	
        if (k<n-2) bp=O2SCL_IX(f,k+1);
	else bp=0.0;
	
        y=ak;
        z=zk;
      }
      
      o2scl_linalg::create_givens(y,z,c,s);
      
      // [GSL] Compute U <= U G
      
      for (i=0;i<M;i++) {
	double Uip=O2SCL_IX2(U,i,k);
	double Uiq=O2SCL_IX2(U,i,k+1);
	O2SCL_IX2(U,i,k)=c*Uip-s*Uiq;
	O2SCL_IX2(U,i,k+1)=s*Uip+c*Uiq;
      }
      
      // [GSL] compute B <= G^T B
      
      double ak1=c*ak-s*zk;
      double bk1=c*bk-s*ap;
      double zk1=-s*bp;

      double ap1=s*bk+c*ap;
      double bp1=c*bp;

      O2SCL_IX(d,k)=ak1;

      ak=ak1;
      bk=bk1;
      zk=zk1;

      ap=ap1;
      bp=bp1;

      if (k < n-2) {
	aq=O2SCL_IX(d,k+2);
      } else {
	aq=0.0;
      }

      y=bk;
      z=zk;
    }

    O2SCL_IX(f,n-2)=bk;
    O2SCL_IX(d,n-1)=ap;

    return;
  }

  /** \brief A special form of qrstep() for SV_decomp()
      
      The vector \c d should be of size <tt>n</tt>, the vector \c f
      should be of size <tt>n</tt>, the matrix U should be of size
      <tt>(M,n)</tt>, and the matrix \c V should be of size
      <tt>(N,n)</tt>.

      A version of qrstep(), but starting at the a'th column of U, the
      a'th column of V, and the a'th entries of \c d and \c f.

      This function is used by \ref SV_decomp().
  */
  template<class vec_t, class vec2_t, class mat_t, class mat2_t>
    void qrstep_sub(size_t M, size_t N, size_t n, 
		    vec_t &d, vec2_t &f, mat_t &U, mat2_t &V, size_t a) {
    
    double y, z;
    double ak, bk, zk, ap, bp, aq, bq;
    size_t i, k;
    
    //std::cout << "M,N,n: " << M << " " << N << " " << n << std::endl;
    
    if (n-a==1) {
      // [GSL] shouldn't happen
      return;
    }
    
    // [GSL] Compute 2x2 svd directly
    
    if (n-a==2) {
      svd2_sub(M,N,d,f,U,V,a);
      return;
    }
    
    // [GSL] Chase out any zeroes on the diagonal
    
    for (i=a;i<n-1;i++) {
      double d_i=O2SCL_IX(d,i);
      //std::cout << "d_i: " << i << " " << n << " "
      //<< d_i << std::endl;
      if (d_i==0.0) {
	chase_out_intermediate_zero(M,n,d,f,U,i);
	return;
      }
    }
    
    // [GSL] Chase out any zero at the end of the diagonal
    double d_nm1=O2SCL_IX(d,n-1);
    //std::cout << "d_nm1: " << d_nm1 << std::endl;
    if (d_nm1==0.0) {
      chase_out_trailing_zero_sub(N,n,d,f,V,a);
      return;
    }

    // [GSL] Apply QR reduction steps to the diagonal and offdiagonal
    
    double d0=O2SCL_IX(d,a);
    double f0=O2SCL_IX(f,a);
    
    double d1=O2SCL_IX(d,a+1);
    double f1=O2SCL_IX(f,a+1);

    //std::cout << "d0,f0,d1,f1: " << d0 << " " << f0 << " " << d1 << " "
    //<< f1 << std::endl;
    
    double mu=trailing_eigenvalue(n,d,f);
    y=d0*d0-mu;
    z=d0*f0;
    
    // [GSL] Set up the recurrence for Givens rotations on a bidiagonal 
    // matrix
    
    ak=0;
    bk=0;
    
    ap=d0;
    bp=f0;
    
    aq=d1;
    bq=f1;
    
    for (k=a;k<n-1;k++) {
      
      double c, s;
      o2scl_linalg::create_givens(y,z,c,s);

      // [GSL] Compute V <= V G
      
      for (i=0;i<N;i++) {
	double Vip=O2SCL_IX2(V,i,k);
	double Viq=O2SCL_IX2(V,i,k+1);
	//std::cout << "Vip,Viq: " << Vip << " " << Viq << std::endl;
	O2SCL_IX2(V,i,k)=c*Vip-s*Viq;
	O2SCL_IX2(V,i,k+1)=s*Vip+c*Viq;
      }
      
      // [GSL] compute B <= B G

      {
	double bk1=c*bk-s*z;
	
	double ap1=c*ap-s*bp;
	double bp1=s*ap+c*bp;
	double zp1=-s*aq;
	
	double aq1=c*aq;
	
	if (k>a) O2SCL_IX(f,k-1)=bk1;
	
	ak=ap1;
	bk=bp1;
	zk=zp1;
	
	ap=aq1;
	
	if (k<n-2) bp=O2SCL_IX(f,k+1);
	else bp=0.0;
	
	y=ak;
	z=zk;
      }

      o2scl_linalg::create_givens(y,z,c,s);
      
      // [GSL] Compute U <= U G
      
      for (i=0;i<M;i++) {
	double Uip=O2SCL_IX2(U,i,k);
	double Uiq=O2SCL_IX2(U,i,k+1);
	//std::cout << "Uip2,Uiq2: " << Uip << " " << Uiq << std::endl;
	O2SCL_IX2(U,i,k)=c*Uip-s*Uiq;
	O2SCL_IX2(U,i,k+1)=s*Uip+c*Uiq;
      }
      
      // [GSL] compute B <= G^T B
      
      //std::cout << "k,bk,ap2: " << k << " " << bk << " " << ap << std::endl;
      //std::cout << "ak,zk,bp: " << ak << " " << zk << " " 
      //<< bp << std::endl;
      
      {
	//std::cout << "prod1: " << c*ak << " " << s*zk << std::endl;
	//std::cout << "prod2: " << c*bk << " " << s*ap << std::endl;
	//std::cout << "prod3: " << s*bk << " " << c*ap << std::endl;
	double ak1=c*ak-s*zk;
	double bk1=c*bk-s*ap;
	double zk1=-s*bp;
	
	double ap1=s*bk+c*ap;
	double bp1=c*bp;
	
	O2SCL_IX(d,k)=ak1;
	
	ak=ak1;
	bk=bk1;
	zk=zk1;
	
	ap=ap1;
	bp=bp1;
	//std::cout << "c,s: " << c << " " << s << std::endl;
	//std::cout << "k,bk,ap: " << k << " " << bk << " " << ap << std::endl;

	if (k < n-2) {
	  aq=O2SCL_IX(d,k+2);
	} else {
	  aq=0.0;
	}
	
	y=bk;
	z=zk;
      }
    }

    O2SCL_IX(f,n-2)=bk;
    O2SCL_IX(d,n-1)=ap;
    //std::cout << "bk,ap: " << bk << " " << ap << std::endl;

    return;
  }
  
#ifdef DOXYGEN
}
#endif
