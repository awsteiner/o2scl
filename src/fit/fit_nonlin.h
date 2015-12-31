/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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
/*
 * Based on the multifit routines in GSL 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004, 2007 Brian Gough
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
#ifndef O2SCL_FIT_NONLIN_H
#define O2SCL_FIT_NONLIN_H

/** \file fit_nonlin.h
    \brief File defining \ref o2scl::fit_nonlin
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/vector.h>
#include <o2scl/fit_base.h>
#include <o2scl/permutation.h>
#include <o2scl/cblas.h>
#include <o2scl/qr.h>
#include <o2scl/qrpt.h>
#include <o2scl/columnify.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Base routines for the nonlinear fitting classes
   */
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_t=boost::numeric::ublas::matrix<double> > class fit_nonlin_b {
    
  protected:

  /** \brief Desc
   */
  double compute_actual_reduction(double fnorm0, double fnorm1) {
    double actred;
  
    if (0.1*fnorm1<fnorm0) {
      double u=fnorm1/fnorm0;
      actred=1-u*u;
    } else {
      actred=-1;
    }
  
    return actred;
  }

  /** \brief Desc
   */
  size_t count_nsing(const size_t ncols, const mat_t &r2) {
    
    size_t i;
    for (i=0;i<ncols;i++) {
      double rii=r2(i,i);
      if (rii == 0) {
	break;
      }
    }
      
    return i;
  }
  
  /** \brief Desc
   */
  void compute_newton_direction
  (size_t n, const mat_t &r2, const permutation &perm2,
   const vec_t &qtf2, vec_t &x) {
    
    size_t i, j, nsing;
    
    for (i=0;i<n;i++) {
      double qtfi=qtf2[i];
      x[i]=qtfi;
    }
    
    nsing=count_nsing(n,r2);
    
    for (i=nsing;i<n;i++) {
      x[i]=0.0;
    }
    
    if (nsing > 0) {
      for (j=nsing;j>0 && j--;) {
	double rjj=r2(j,j);
	double temp=x[j]/rjj;
	    
	x[j]=temp;
	    
	for (i=0;i<j;i++) {
	  double rij=r2(i,j);
	  double xi=x[i];
	  x[i]=xi-rij*temp;
	}
      }
    }
  
    perm2.apply_inverse(x);
    return;
  }

  /** \brief Desc
   */
  void compute_newton_bound
  (size_t nd, size_t np, const mat_t &r2, const vec_t &x, double dxnorm, 
   const permutation &perm, const vec_t &diag, vec_t &w) {

    size_t i, j;
      
    size_t nsing=count_nsing(np,r2);
      
    if (nsing<np) {
      for(i=0;i<nd;i++) {
	w[i]=0.0;
      }
      return;
    }
      
    for (i=0;i<np;i++) {
      size_t pi=perm[i];
      double dpi=diag[pi];
      double xpi=x[pi];
      w[i]=dpi*(dpi*xpi/dxnorm);
    }
      
    for (j=0;j<np;j++) {
      double sum=0.0;
	
      for (i=0;i<j;i++) {
	sum +=r2(i,j)*w[i];
      }
	
      {
	double rjj=r2(j,j);
	double wj=w[j];
	w[j]=(wj-sum)/rjj;
      }
    }
    return;
  }

  /** \brief Desc
   */
  void compute_gradient_direction
  (size_t n, const mat_t &r, const permutation &p, const vec_t &qtf2, 
   const vec_t &diag, vec_t &g) {
   
    size_t i, j;
      
    for (j=0;j<n;j++) {
      double sum=0;
	
      for (i=0;i <= j;i++) {
	sum +=r(i,j)*qtf2[i];
      }
	
      {
	size_t pj=p[j];
	double dpj=diag[pj];
	  
	g[j]=sum/dpj;
      }
    }
    return;
  }

  /** \brief Desc 
   */
  void update_diag(size_t n, const mat_t &J, vec_t &diag2) {
    
    for (size_t j=0;j<n;j++) {
      double cnorm,diagj,sum=0;
      for (size_t i=0;i<n;i++) {
	double Jij=J(i,j);
	sum+=Jij*Jij;
      }
      if (sum == 0) {
	sum=1.0;
      }
	
      cnorm=sqrt(sum);
      diagj=diag2[j];
	
      if (cnorm>diagj) {
	diag2[j]=cnorm;
      }
    }
    return;
  }
    
  /** \brief Euclidean norm of vector \c f of length \c n, scaled by
      vector \c d
  */
  double scaled_enorm(const vec_t &d, size_t n, const vec_t &f) {
    double e2=0;
    for (size_t i=0;i<n;i++) {
      double di=d[i];
      double u=di*f[i];
      e2+=u*u;
    }
    return sqrt(e2);
  }

  /** \brief Desc 
   */
  double compute_delta(vec_t &diag2, size_t n, const vec_t &x) {
    double Dx=scaled_enorm(diag2,n,x);
    // [GSL] The generally recommended value from MINPACK is 100
    double factor=100; 
    return (Dx>0) ? factor*Dx : factor;
  }
    
  /** \brief Desc 
   */
  void compute_rptdx(const mat_t &r2, const permutation &p,
		     size_t N, vec_t &dx, vec_t &rptdx2) {
    
    for (size_t i=0;i<N;i++) {
      double sum=0.0;
      for (size_t j=i;j<N;j++) {
	sum+=r2(i,j)*dx[p[j]];
      }
      rptdx2[i]=sum;
    }

    return;
  }
    
  /** \brief Compute the solution to a least squares system

      \verbatim
      This function computes the solution to the least squares system
      
      phi=[ A x=b ,lambda D x=0 ]^2
      
      where A is an M by N matrix,D is an N by N diagonal matrix,lambda
      is a scalar parameter and b is a vector of length M.
      
      The function requires the factorization of A into A=Q R P^T,
      where Q is an orthogonal matrix,R is an upper triangular matrix
      with diagonal elements of non-increasing magnitude and P is a
      permuation matrix. The system above is then equivalent to
      
      [ R z=Q^T b,P^T (lambda D) P z=0 ]
      
      where x=P z. If this system does not have full rank then a least
      squares solution is obtained. On output the function also provides
      an upper triangular matrix S such that
      
      P^T (A^T A+lambda^2 D^T D) P=S^T S
      
      Parameters,
      
      r: On input,contains the full upper triangle of R. On output the
      strict lower triangle contains the transpose of the strict upper
      triangle of S,and the diagonal of S is stored in sdiag.  The full
      upper triangle of R is not modified.
      
      p: the encoded form of the permutation matrix P. column j of P is
      column p[j] of the identity matrix.
      
      lambda,diag: contains the scalar lambda and the diagonal elements
      of the matrix D
      
      qtb: contains the product Q^T b
      
      x: on output contains the least squares solution of the system
      
      wa: is a workspace of length N
      \endverbatim
      
  */
  int qrsolv(size_t n, mat_t &r2, const permutation &p,
	     const double lambda, const vec_t &diag2,
	     const vec_t &qtb, vec_t &x, vec_t &sdiag2, vec_t &wa) {

    size_t i, j, k, nsing;
      
    // [GSL] Copy r and qtb to preserve input and initialise s. In 
    // particular, save the diagonal elements of r in x
      
    for (j=0;j<n;j++) {
      double rjj=r2(j,j);
      double qtbj=qtb[j];
	
      for (i=j+1;i<n;i++) {
	r2(i,j)=r2(j,i);
      }

      x[j]=rjj;
      wa[j]=qtbj;
    }
      
    // [GSL] Eliminate the diagonal matrix d using a Givens rotation
      
    for (j=0;j<n;j++) {

      double qtbpj;
	
      size_t pj=p[j];
      double diagpj=lambda*diag2[pj];
      if (diagpj == 0) {
	continue;
      }
	
      sdiag2[j]=diagpj;
	
      for (k=j+1;k<n;k++) {
	sdiag2[k]=0.0;
      }
	
      // [GSL] The transformations to eliminate the row of d modify only a
      // single element of qtb beyond the first n, which is initially
      // zero 
	
      qtbpj=0;
	
      for (k=j;k<n;k++) {

	// [GSL] Determine a Givens rotation which eliminates the
	// appropriate element in the current row of d
	  
	double sine, cosine;
	  
	double wak=wa[k];
	double rkk=r2(k,k);
	double sdiagk=sdiag2[k];
	  
	if (sdiagk == 0) {
	  continue;
	}
	  
	if (fabs(rkk)<fabs(sdiagk)) {
	  double cotangent=rkk/sdiagk;
	  sine=0.5/sqrt(0.25+0.25*cotangent*cotangent);
	  cosine=sine*cotangent;
	} else {
	  double tangent=sdiagk/rkk;
	  cosine=0.5/sqrt(0.25+0.25*tangent*tangent);
	  sine=cosine*tangent;
	}
	  
	// [GSL] Compute the modified diagonal element of r and the
	// modified element of [qtb,0]
	  
	{
	  double new_rkk=cosine*rkk+sine*sdiagk;
	  double new_wak=cosine*wak+sine*qtbpj;
            
	  qtbpj=-sine*wak+cosine*qtbpj;

	  r2(k,k)=new_rkk;
	  wa[k]=new_wak;
	}
	  
	// [GSL] Accumulate the transformation in the row of s
	  
	for (i=k+1;i<n;i++) {
	  double rik=r2(i,k);
	  double sdiagi=sdiag2[i];
            
	  double new_rik=cosine*rik+sine*sdiagi;
	  double new_sdiagi=-sine*rik+cosine*sdiagi;

	  r2(i,k)=new_rik;
	  sdiag2[i]=new_sdiagi;
	}
      }
	
      // [GSL] Store the corresponding diagonal element of s and
      // restore the corresponding diagonal element of r
      {
	double rjj=r2(j,j);
	double xj=x[j];
	  
	sdiag2[j]=rjj;
	r2(j,j)=xj;
      }
	
    }
      
    // [GSL] Solve the triangular system for z. If the system is singular 
    // then obtain a least squares solution 

    nsing=n;
      
    for (j=0;j<n;j++) {
      double sdiagj=sdiag2[j];
	
      if (sdiagj == 0) {
	nsing=j;
	break;
      }
    }
      
    for (j=nsing;j<n;j++) {
      wa[j]=0.0;
    }
      
    for (k=0;k<nsing;k++) {
      double sum=0;
	
      j=(nsing-1)-k;
	
      for (i=j+1;i<nsing;i++) {
	sum+=r2(i,j)*wa[i];
      }
	
      {
	double waj=wa[j];
	double sdiagj=sdiag2[j];
	  
	wa[j]=(waj-sum)/sdiagj;
      }
    }
      
    // [GSL] Permute the components of z back to the components of x 
      
    for (j=0;j<n;j++) {
      size_t pj=p[j];
      double waj=wa[j];
	
      x[pj]=waj;
    }
      
    return success;
  }
    
  /** \brief Desc 
   */
  void compute_newton_correction
  (size_t n, const mat_t &r2, const vec_t &sdiag2, const permutation &p,
   vec_t &x, double dxnorm, const vec_t &diag2, vec_t &w2) {
    
    for (size_t i=0;i<n;i++) {
      size_t pi=p[i];
      double dpi=diag2[pi];
      double xpi=x[pi];
	
      w2[i]=dpi*(dpi*xpi)/dxnorm;
    }
      
    for (size_t j=0;j<n;j++) {
      double sj=sdiag2[j];
      double wj=w2[j];
      double tj=wj/sj;
	
      w2[j]=tj;
	
      for (size_t i=j+1;i<n;i++) {
	double rij=r2(i,j);
	double wi=w2[i];
	  
	w2[i]=wi-rij*tj;
      }
    }
    return;
  }
    
  /** \brief Determine Levenburg-Marquardt parameter
   */
  void lmpar(mat_t &r2, const permutation &perm2, const vec_t &qtf2, 
	     const vec_t &diag2, double delta2, double *par_inout, 
	     vec_t &newton2, vec_t &gradient2, vec_t &sdiag2,
	     vec_t &x, vec_t &w2, size_t nparm, size_t ndata) {
    
    double dxnorm, gnorm, fp, fp_old, par_lower, par_upper, par_c;
      
    double par2=*par_inout;
      
    size_t iter2=0;
      
    this->compute_newton_direction(nparm,r2,perm2,qtf2,newton2);
      
    // [GSL] Evaluate the function at the origin and test for acceptance of
    // the Gauss-Newton direction.
      
    dxnorm=this->scaled_enorm(diag2,nparm,newton2);
      
    fp=dxnorm-delta2;
      
    if (fp<=0.1*delta2) {

      for(size_t i=0;i<nparm;i++) {
	x[i]=newton2[i];
      }
	
      *par_inout=0;
	
      return;
    }
      
    this->compute_newton_bound(ndata,nparm,r2,newton2,dxnorm,perm2,diag2,w2);
      
    {
      double wnorm=o2scl_cblas::dnrm2(ndata,w2);
      double phider=wnorm*wnorm;
	
      // [GSL] w == zero if r rank-deficient,
      // then set lower bound to zero form MINPACK
      // Hans E. Plesser 2002-02-25 (hans.plesser@itf.nlh.no)
      if (wnorm>0) {
	par_lower=fp/(delta2*phider);
      } else {
	par_lower=0.0;
      }
    }
      
    this->compute_gradient_direction(nparm,r2,perm2,qtf2,diag2,gradient2);
      
    gnorm=o2scl_cblas::dnrm2(nparm,gradient2);
      
    par_upper=gnorm/delta2;
      
    if (par_upper == 0) {
      double mint;
      if (delta2<0.1) mint=delta2;
      else mint=0.1;
      double dbl_min=std::numeric_limits<double>::min();
      par_upper=dbl_min/mint;
    }
      
    if (par2>par_upper) {
      par2=par_upper;
    } else if (par2<par_lower) {
      par2=par_lower;
    }
      
    if (par2 == 0) {
      par2=gnorm/dxnorm;
    }
      
    // [GSL] Beginning of iteration 
      
    while (true) {
      
      iter2++;

#ifdef O2SCL_NEVER_DEFINED      
#ifdef BRIANSFIX
      // [GSL] Seems like this is described in the paper but 
      // not in the MINPACK code
      
      if (par2<par_lower || par2>par_upper) {
	par2=GSL_MAX_DBL(0.001*par_upper,sqrt(par_lower*par_upper));
      }
#endif
#endif
      
      // [GSL] Evaluate the function at the current value of par
      
      if (par2 == 0) {
	par2=GSL_MAX_DBL(0.001*par_upper,GSL_DBL_MIN);
      }
      
      // [GSL] Compute the least squares solution of 
      // [ R P x-Q^T f, sqrt(par) D x] for A=Q R P^T 
      
      {
	double sqrt_par=sqrt(par2);
	
	qrsolv(nparm,r2,perm2,sqrt_par,diag2,qtf2,x,sdiag2,w2);
      }
      
      dxnorm=scaled_enorm(diag2,nparm,x);
      
      fp_old=fp;
      
      fp=dxnorm-delta2;
      
      // [GSL] If the function is small enough, accept the current
      // value of par
      
      if (fabs(fp)<=0.1*delta2) {
	*par_inout=par2;
	return;
      }
      if (par_lower == 0 && fp<=fp_old && fp_old<0) {
	*par_inout=par2;
	return;
      }
      
      // [GSL] Check for maximum number of iterations
      
      if (iter2 == 10) {
	*par_inout=par2;
	return;
      }
      
      // [GSL] Compute the Newton correction
      
      this->compute_newton_correction(nparm,r2,sdiag2,perm2,x,
				      dxnorm,diag2,w2);
      {
	double wnorm=o2scl_cblas::dnrm2(ndata,w2);
	par_c=fp/(delta2*wnorm*wnorm);
      }
      
      // [GSL] Depending on the sign of the function,
      // update par_lower or par_upper
      
      if (fp>0) {
	if (par2>par_lower) {
	  par_lower=par2;
	}
      } else if (fp<0) {
	if (par2<par_upper) {
	  par_upper=par2;
	}
      }
      
      // [GSL] Compute an improved estimate for par
      
      par2=GSL_MAX_DBL(par_lower,par2+par_c);
      
    }

    O2SCL_ERR("Sanity check failed in fit_nonlin_b::lmpar().",
	      exc_esanity);
    return;
  }

  /// Compute trial step, \f$ \mathrm{trial}=\mathrm{x}+\mathrm{dx} \f$
  void compute_trial_step(size_t N, vec_t &x, vec_t &dx,
			  vec_t &trial) {
    for (size_t i=0;i<N;i++) {
      trial[i]=x[i]+dx[i];
    }
    return;
  }

  /** \brief Compute the root of the sum of the squares of 
      the columns of \c J

      This computes 
      \f[
      \mathrm{diag\_vec}_j = \sqrt{\sum_{i=0}^{\mathrm{ndata}-1} J_{ij}}
      \f]
      for \f$ 0\leq j \leq \mathrm{nparm}-1 \f$ . If any of the
      columns of \c J is all zero, then the corresponding entry in
      \c diag_vec is set to one instead.
  */
  int compute_diag(size_t nparm, size_t ndata, 
		   const mat_t &J, vec_t &diag_vec) {
      
    for (size_t j=0;j<nparm;j++) {
      double sum=0.0;
	
      for (size_t i=0;i<ndata;i++) {
	sum+=J(i,j)*J(i,j);
      }
      if (sum == 0) {
	sum=1.0;
      }
	
      diag_vec[j]=sqrt(sum);
    }
    return 0;
  }

  /** \brief Compute the covarance matrix \c covar given 
      the Jacobian \c J

      Given a \c m by \c n Jacobian matrix \c J (where \c m must not
      be less than \c n), and a relative tolerance \c epsrel, this
      function computes the entries of the \c n by \c n covariance
      matrix \c covar. The allocation for \c covar must be performed
      beforehand.

      This function is basically the equivalent of the function
      <tt>gsl_multifit_covar()</tt>, but rewritten for generic
      vector and matrix types.

      The workspace \c work1 is used here.

      \comment
      Note that in the GSL example for non-linear fitting, the value
      of \c epsrel is set to zero, so this class sets \ref
      tol_rel_covar to zero by default. We could remove the \c epsrel
      parameter in this function, but it may be one day useful to
      separate this function so we leave \c epsrel as a parameter.
      \endcomment
  */
  int covariance(size_t m, size_t n, const mat_t &J, 
		 mat_t &covar, vec_t &norm, mat_t &r, 
		 vec_t &tau, permutation &perm, double epsrel) {
      
    if (m<n) {
      O2SCL_ERR2("Jacobian must have m>=n in ",
		     "fit_nonlin_b::covariance().",exc_efailed);
    }
      
    double tolr;
      
    size_t kmax=0;
    
    int signum=0;
      
    for(size_t i=0;i<m;i++) {
      for(size_t j=0;j<n;j++) {
	r(i,j)=J(i,j);
      }
    }
    o2scl_linalg::QRPT_decomp(m,n,r,tau,perm,signum,norm);
      
    // [GSL] Form the inverse of R in the full upper triangle of R
      
    tolr=epsrel*fabs(r(0,0));

    for (size_t k=0;k<n;k++) {
      double rkk=r(k,k);
	  
      if (fabs(rkk)<=tolr) {
	break;
      }
	  
      r(k,k)=1.0/rkk;
	  
      for (size_t j=0;j<k;j++) {
	double t=r(j,k)/rkk;
	r(j,k)=0.0;
	    
	for (size_t i=0;i<=j;i++) {
	  double rik=r(i,k);
	  double rij=r(i,j);
	  r(i,k)=rik-t*rij;
	}
      }
      kmax=k;
    }
      
    // [GSL] Form the full upper triangle of the inverse of R^T R in 
    // the full upper triangle of R
      
    for (size_t k=0;k<=kmax;k++) {
      for (size_t j=0;j<k;j++) {
	double rjk=r(j,k);
	  
	for (size_t i=0;i<=j;i++) {
	  double rij=r(i,j);
	  double rik=r(i,k);
	  r(i,j)=rij+rjk*rik;
	}
      }
	
      double t=r(k,k);
	
      for (size_t i=0;i<=k;i++) {
	r(i,k)*=t;
      }
    }
      
    // [GSL] Form the full lower triangle of the covariance matrix in 
    // the strict lower triangle of R and in w
      
    for (size_t j=0;j<n;j++) {
      size_t pj=perm[j];
	
      for (size_t i=0;i<=j;i++) {
	size_t pi=perm[i];
	  
	double rij;
	  
	if (j>kmax) {
	  r(i,j)=0.0;
	  rij=0.0;
	} else {
	  rij=r(i,j);
	}
	  
	if (pi>pj) {
	  r(pi,pj)=rij;
	} else if (pi<pj) {
	  r(pj,pi)=rij;
	}
	  
      }
	
      double rjj=r(j,j);
      covar(pj,pj)=rjj;
    }
      
      
    // [GSL] symmetrize the covariance matrix
      
    for (size_t j=0;j<n;j++) {
      for (size_t i=0;i<j;i++) {
	double rji=r(j,i);
	covar(j,i)=rji;
	covar(i,j)=rji;
      }
    }
    
    return success;
  }
  
  public:

  fit_nonlin_b() {
    tol_rel_covar=0.0;
  }
  
  /** \brief The relative tolerance for the computation of 
      the covariance matrix (default 0)
  */
  double tol_rel_covar;

  /// Test if converged
  int test_delta_f(size_t nparm, vec_t &dx, vec_t &x, 
		   double l_epsabs, double l_epsrel) {
    
    if (l_epsrel<0.0) {
      O2SCL_ERR2("Relative tolerance less than zero ",
		     "in fit_nonlin_b::test_delta_f().",exc_einval);
    } 
    
    for(size_t i=0;i<nparm;i++) {
      double tolerance=l_epsabs+l_epsrel*fabs(x[i]);
      if (fabs(dx[i])>=tolerance) {
	return gsl_continue;
      }
    }
    return success;
  }
  
  /// Test if converged
  int test_gradient_f(size_t nparm, vec_t &g, double l_epsabs) {
    
    double residual=0.0;
    
    if (l_epsabs<0.0) {
      O2SCL_ERR2("Absolute tolerance less than zero ",
		     "in fit_nonlin_b::test_gradient_f().",exc_einval);
    } 
    
    for(size_t i=0;i<nparm;i++) {
      residual+=fabs(g[i]);
    }
    if (residual<l_epsabs) {
      return success;
    } 
    return gsl_continue;
  }
  
  };
  
  /** \brief Non-linear least-squares fitting class (GSL)
      
      The GSL-based fitting class using a Levenberg-Marquardt type
      algorithm. The algorithm stops when
      \f[
      |dx_i| < \mathrm{tol\_abs}+\mathrm{tol\_rel}\times|x_i|
      \f]
      where \f$dx\f$ is the last step and \f$x\f$ is the current
      position. If test_gradient is true, then additionally fit()
      requires that
      \f[
      \sum_i |g_i| < \mathrm{tol\_abs}
      \f]
      where \f$g_i\f$ is the \f$i\f$-th component of the gradient of
      the function \f$\Phi(x)\f$ where
      \f[
      \Phi(x) = || F(x) ||^2
      \f]

      Default template arguments
      - \c func_t - \ref gen_fit_funct\<\>
      - \c vec_t - \ref boost::numeric::ublas::vector \<double \>
      - \c mat_t - \ref boost::numeric::ublas::matrix \<double \>
      
      \todo Allow the user to specify the derivatives
      \todo Fix so that the user can specify automatic
      scaling of the fitting parameters, where the initial
      guess are used for scaling so that the fitting parameters
      are near unity.
      
      \future Some of these member functions (like
      update_diag()) don't depend on member data and could be
      possibly be moved to a parent class?
  */
  template<class func_t=gen_fit_funct<>, 
    class vec_t=boost::numeric::ublas::vector<double>,
    class mat_t=boost::numeric::ublas::matrix<double> >
    class fit_nonlin : public fit_nonlin_b<vec_t,mat_t>,
    public fit_base<func_t,vec_t,mat_t> {
    
  protected:

  /// Function to fit
  func_t *cff;
  
  /// Trial step
  vec_t x_trial;

  /// Trial function value
  vec_t f_trial;

  /// Number of iterations 
  size_t iter;

  /** \brief Desc 
   */
  double xnorm;

  /** \brief Desc */
  double fnorm;

  /** \brief Desc */
  double delta;

  /** \brief Desc */
  double par;

  /** \brief Desc */
  mat_t r;

  /** \brief Desc */
  vec_t tau;

  /** \brief Desc */
  vec_t diag;

  /** \brief Desc */
  vec_t qtf;

  /** \brief Desc */
  vec_t df;

  /** \brief Desc */
  vec_t rptdx;

  /** \brief Desc */
  vec_t newton;

  /** \brief Desc */
  vec_t gradient;

  /** \brief Desc */
  vec_t sdiag;

  /** \brief Desc */
  vec_t w;

  /** \brief Desc */
  vec_t work1;

  /** \brief Desc */
  permutation perm;
    
  /// Number of data points
  size_t ndata;

  /// Number of parameters
  size_t nparm;
    
  /// Desc
  vec_t g_;

  /// Desc
  mat_t J_;

  /// Desc
  vec_t *x_;

  /// Free allocated memory
  void free() {
    if (ndata>0) {
      perm.free();
      f_trial.clear();
      x_trial.clear();
      f_.clear();
      dx_.clear();
      g_.clear();
      J_.clear();
      r.clear();
      rptdx.clear();
      w.clear();
      work1.clear();
      newton.clear();
      gradient.clear();
      sdiag.clear();
      qtf.clear();
      df.clear();
      diag.clear();
      tau.clear();
      ndata=0;
      nparm=0;
    }
  }

  public:
    
  fit_nonlin() {
    use_scaled=true;
    test_gradient=false;
    ndata=0;
    nparm=0;
    err_nonconv=true;
  }

  virtual ~fit_nonlin() {
    free();
  }

  /** \brief If true, call the error handler if fit()
      does not converge (default true)
  */
  bool err_nonconv;
  
  /// Print the progress in the current iteration
  virtual int print_iter(size_t nv, vec_t &x, vec_t &dx, int iter2,
			 double l_epsabs, double l_epsrel) {

    if (this->verbose<=0) return 0;
  
    // Variables 'val' and 'lim' set to zero to avoid 
    // uninitialized variable errors. In reality, they're
    // always set in the loop below.
    double val=0.0, lim=0.0, diff_max=0.0;
    for(size_t i=0;i<nv;i++) {
      double tolerance=l_epsabs+l_epsrel*fabs(x[i]);
      if (fabs(dx[i])>=tolerance) {
	val=fabs(dx[i]);
	lim=tolerance;
	i=nv;
      } else {
	double diff=tolerance-fabs(dx[i]);
	if (diff>diff_max) {
	  diff_max=diff;
	  val=fabs(dx[i]);
	  lim=tolerance;
	}
      }
    }

    std::cout << "Iteration: " << iter2 <<  std::endl;
    std::cout << "x: ";
    for(size_t i=0;i<nv;i++) {
      std::cout << x[i] << " ";
    }
    std::cout << std::endl;
    std::cout << " Val: " << val << " Lim: " << lim << std::endl;
    if (this->verbose>1) {
      std::cout << "Press a key and type enter to continue. ";
      char ch;
      std::cin >> ch;
    }

    return 0;
  }
    
  /// Allocate memory with \c n data points and \c p parameters
  void resize(size_t n, size_t p) {

    free();
      
    ndata=n;
    nparm=p;

    perm.allocate(nparm);
    r.resize(ndata,nparm);
    if (ndata<nparm) {
      tau.resize(ndata);
      for(size_t i=0;i<ndata;i++) tau[i]=0.0;
    } else {
      tau.resize(nparm);
      for(size_t i=0;i<nparm;i++) tau[i]=0.0;
    }
    df.resize(ndata);
    rptdx.resize(ndata);
    qtf.resize(ndata);
    w.resize(ndata);
    diag.resize(nparm);
    work1.resize(nparm);
    newton.resize(nparm);
    gradient.resize(nparm);
    sdiag.resize(nparm);
    for(size_t i=0;i<ndata;i++) {
      for(size_t j=0;j<nparm;j++) {
	r(i,j)=0.0;
      }
      rptdx[i]=0.0;
      w[i]=0.0;
      qtf[i]=0.0;
      df[i]=0.0;
    }
    for(size_t i=0;i<nparm;i++) {
      work1[i]=0.0;
      newton[i]=0.0;
      gradient[i]=0.0;
      sdiag[i]=0.0;
      diag[i]=0.0;
    }

    x_trial.resize(nparm);
    f_trial.resize(ndata);
    f_.resize(ndata);
    dx_.resize(nparm);
    g_.resize(ndata);
    J_.resize(ndata,nparm);
    
    return;
  }

  /// The last step taken in parameter space
  vec_t dx_;

  /// Desc
  vec_t f_;

  /** \brief Set the initial values of the parameters and 
      the fitting function to use for the next call to 
      iterate()
  */
  int set(size_t npar, vec_t &parms, func_t &fitfun) {
    
    cff=&fitfun;

    if (fitfun.get_ndata()==0 || npar==0) {
      O2SCL_ERR2("Either zero data or zero parameters in ",
		 "fit_nonlin::fit().",exc_einval);
    }

    if (ndata!=fitfun.get_ndata() || nparm!=npar) {
      resize(fitfun.get_ndata(),npar);
    }

    x_=&parms;

    int signum;
      
    // Evaluate function and Jacobian at x
    (*cff)(nparm,*x_,ndata,f_);
    cff->jac(nparm,*x_,ndata,f_,J_);

    par=0;
    iter=1;
    fnorm=o2scl_cblas::dnrm2(ndata,f_);
      
    for(size_t i=0;i<nparm;i++) {
      dx_[i]=0.0;
    }
      
    // [GSL] store column norms in diag 
      
    if (use_scaled) {
      this->compute_diag(nparm,ndata,J_,diag);
    } else {
      for(size_t i=0;i<nparm;i++) {
	diag[i]=1.0;
      }
    }
      
    // [GSL] set delta to 100 |D x| or to 100 if |D x| is zero 
      
    xnorm=this->scaled_enorm(diag,nparm,*x_);
    delta=this->compute_delta(diag,nparm,*x_);
      
    // [GSL] Factorize J into QR decomposition 

    for(size_t i=0;i<ndata;i++) {
      for(size_t j=0;j<nparm;j++) {
	r(i,j)=J_(i,j);
      }
    }

    o2scl_linalg::QRPT_decomp(ndata,nparm,r,tau,perm,signum,work1);
      
    for(size_t ii=0;ii<rptdx.size();ii++) rptdx[ii]=0.0;
    for(size_t ii=0;ii<w.size();ii++) w[ii]=0.0;
      
    // [GSL] Zero the trial vector, as in the alloc function 
      
    for(size_t i=0;i<ndata;i++) f_trial[i]=0.0;

    return success;
  }

  /** \brief Perform an iteration
   */
  int iterate() {

    vec_t &f=f_;
    mat_t &J=J_;
    vec_t &dx=dx_;
    vec_t &x=*x_;

    double prered, actred;
    double pnorm, fnorm1, fnorm1p, gnorm;
    double ratio;
    double dirder;
      
    int niter=0;
      
    double p1=0.1, p25=0.25, p5=0.5, p75=0.75, p0001=0.0001;
      
    if (fnorm == 0.0) {
      return success;
    }
      
    // [GSL] Compute qtf=Q^T f 
      
    for(size_t i=0;i<ndata;i++) {
      qtf[i]=f[i];
    }
    o2scl_linalg::QR_QTvec(ndata,nparm,r,tau,qtf);
    
    // [GSL] Compute norm of scaled gradient 

    this->compute_gradient_direction(nparm,r,perm,qtf,diag,gradient);
    size_t iamax=vector_max_index<vec_t,double>(nparm,gradient);
    gnorm=fabs(gradient[iamax]/fnorm);

    // [GSL] Determine the Levenberg-Marquardt parameter

    bool lm_iteration=true;
    while (lm_iteration) {
      
      niter++;
      
      this->lmpar(r,perm,qtf,diag,delta,&par,newton,
		  gradient,sdiag,dx,w,nparm,ndata);
		  
      // [GSL] Take a trial step 

      // [GSL] reverse the step to go downhill 
      for(size_t i=0;i<nparm;i++) dx[i]*=-1.0;
      
      this->compute_trial_step(nparm,x,dx,x_trial);

      pnorm=this->scaled_enorm(diag,nparm,dx);

      if (iter == 1) {
	if (pnorm<delta) {
	  delta=pnorm;
	}
      }
      
      // [GSL] Evaluate function at x+p 
      (*cff)(nparm,x_trial,ndata,f_trial);
      
      fnorm1=o2scl_cblas::dnrm2(ndata,f_trial);
      
      // [GSL] Compute the scaled actual reduction 
      
      actred=this->compute_actual_reduction(fnorm,fnorm1);

      // [GSL] Compute rptdx=R P^T dx,noting that |J dx|=|R P^T dx|
      
      this->compute_rptdx(r,perm,nparm,dx,rptdx);
      
      fnorm1p=o2scl_cblas::dnrm2(ndata,rptdx);

      // [GSL] Compute the scaled predicted reduction,
      // |J dx|^2+2 par |D dx|^2

      double t1=fnorm1p/fnorm;
      double t2=(sqrt(par)*pnorm)/fnorm;
      
      prered=t1*t1+t2*t2/p5;
      dirder=-(t1*t1+t2*t2);

      // [GSL] compute the ratio of the actual to predicted reduction 

      if (prered>0) {
	ratio=actred/prered;
      } else {
	ratio=0;
      }
      
      // [GSL] update the step bound 
      
      if (ratio>p25) {
	if (par == 0 || ratio >= p75) {
	  delta=pnorm/p5;
	  par *= p5;
	}
      } else {
	double temp=(actred >= 0) ? p5 : p5*dirder/
	  (dirder+p5*actred);
	
	if (p1*fnorm1 >= fnorm || temp<p1 ) {
	  temp=p1;
	}
	
	if (delta<pnorm/p1) delta=temp*delta;
	else delta=temp*pnorm/p1;
	
	par /= temp;
      }
      
      
      // Test for successful iteration, termination and 
      // stringent tolerances.

      lm_iteration=false;

      double dbl_eps=std::numeric_limits<double>::epsilon();

      if (ratio >= p0001) {

	vector_copy(nparm,x_trial,x);
	vector_copy(ndata,f_trial,f);
	
	cff->jac(nparm,x_trial,ndata,f_trial,J_);

	// [GSL] wa2_j =diag_j*x_j 
	xnorm=this->scaled_enorm(diag,nparm,x);
	fnorm=fnorm1;
	iter++;
	
	// [GSL] Rescale if necessary 

	if (use_scaled) {
	  this->update_diag(nparm,J,diag);
	}
	  
	{
	  int signum;
	  for(size_t i=0;i<ndata;i++) {
	    for(size_t j=0;j<nparm;j++) {
	      r(i,j)=J(i,j);
	    }
	  }
	  o2scl_linalg::QRPT_decomp(ndata,nparm,r,tau,perm,signum,work1);
	}
	  
	return success;
      } else if (fabs(actred)<=dbl_eps && 
		 prered<=dbl_eps && p5*ratio<=1.0) {
	return exc_etolf;
      } else if (delta<=dbl_eps*xnorm) {
	return exc_etolx;
      } else if (gnorm<=dbl_eps) {
	return exc_etolg;
      } else if (niter<10) {
	lm_iteration=true;
      }

    }
      
    return gsl_continue;
  }
    
  /** \brief Fit the data specified in (xdat,ydat) to
      the function \c fitfun with the parameters in \c par.
	
      The covariance matrix for the parameters is returned in \c covar
      and the value of \f$ \chi^2 \f$ is returned in \c chi2.
  */
  virtual int fit(size_t npar, vec_t &parms, mat_t &covar, 
		  double &chi2, func_t &fitfun) {

    set(npar,parms,fitfun);
    
    int status;
    size_t niter=0;
    do {
      
      niter++;
      status=iterate();

      if (status) {
	break;
      }

      status=this->test_delta_f(nparm,dx_,parms,this->tol_abs,
				this->tol_rel);
      if (test_gradient && status==gsl_continue) {
	o2scl_cblas::dgemv(o2scl_cblas::o2cblas_RowMajor,
			   o2scl_cblas::o2cblas_NoTrans,
			   fitfun.get_ndata(),npar,1.0,J_,f_,0.0,g_);
	status=this->test_gradient_f(nparm,g_,this->tol_abs);
      }
	
      this->print_iter(npar,parms,dx_,niter,this->tol_abs,this->tol_rel);

    } while (status == gsl_continue && niter<this->ntrial);

    if (niter>=this->ntrial) {
      O2SCL_CONV2_RET("Function fit_nonlin::fit() ",
		      "exceeded max. number of iterations.",
		      exc_emaxiter,err_nonconv);
    }

    // Compute covariance matrix
    this->covariance(fitfun.get_ndata(),npar,J_,covar,work1,
		     r,tau,perm,this->tol_rel_covar);

    // Compute chi squared
    chi2=o2scl_cblas::dnrm2(fitfun.get_ndata(),f_);
    chi2*=chi2;

    return 0;
  }
    
  /// If true, test the gradient also (default false)
  bool test_gradient;

  /** \brief Use the scaled routine if true (default true)
   */
  bool use_scaled;

  /// Return string denoting type ("fit_nonlin")
  virtual const char *type() { return "fit_nonlin"; }

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
