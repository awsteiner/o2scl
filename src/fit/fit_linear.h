/*
  -------------------------------------------------------------------
  
  Copyright (C) 2013-2016, Andrew W. Steiner
  
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
/* multifit/multilinear.c
 * 
 * Copyright (C) 2000, 2007, 2010 Brian Gough
 * Copyright (C) 2013, 2015 Patrick Alken
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
#ifndef O2SCL_FIT_LINEAR_H
#define O2SCL_FIT_LINEAR_H

/** \file fit_linear.h
    \brief File defining \ref o2scl::fit_linear
*/

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/vector.h>
#include <o2scl/fit_base.h>
#include <o2scl/permutation.h>
#include <o2scl/cblas.h>
#include <o2scl/svd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Linear least-squares fitting class (GSL)
   */
  template<class vec_t=boost::numeric::ublas::vector<double>,
    class mat_t=boost::numeric::ublas::matrix<double> >
    class fit_linear {
    
  protected:

  /// \name Storage
  //@{
  /** \brief Local copy of <tt>xpred</tt> 
      (and used as workspace by the SV decomposition)
  */
  mat_t A;
  /// The first unitary matrix from the SV decomposition
  mat_t Q;
  /// Workspace for the SV decomposition and storage for \f$ Q S^{-1} \f$
  mat_t QSI;

  /// The singular values from the SV decomposition
  vec_t S; 
  /// SV decomposition workspace and also used to store new predicted values
  vec_t xt; 
  /// Balancing factors for A
  vec_t D; 
  /// Only used for the weighted fit (not yet implemented)
  vec_t t;
  //@}

  /// Number of parameters
  size_t size_par;
  /// Number of data points
  size_t size_dat;

  public:
    
  fit_linear() {
    column_scaling=true;
    tol=2.2204460492503131e-16;
    size_par=0;
    size_dat=0;
  }

  virtual ~fit_linear() {
  }

  /** \brief If true, discard fit components if the associated singular 
      value becomes too small (default true)
  */
  bool column_scaling;

  /// Tolerance (default \f$ \sim 2.22\times 10^{-16} \f$)
  double tol;

  /** \brief The rank of the linear system from the last call to 
      <tt>fit_linear()</tt>
  */
  size_t rank;

  /** \brief Perform the SV decomposition
   */
  virtual void fit_svd(size_t ndat, size_t npar) {
    o2scl_linalg::SV_decomp_mod(ndat,npar,A,QSI,Q,S,xt);
    return;
  }
  
  /** \brief Perform a least-squares fit of a linear system 

      This function performs a least-squares fit of the 
      system
      \f[
      \mathrm{xpred} \cdot \mathrm{parms} = \mathrm{ydat}
      \f]
	
      The variance-covariance matrix for the parameters is returned in
      \c covar and the value of \f$ \chi^2 \f$ is returned in \c chi2.
  */
  virtual void fit(size_t npar, size_t ndat, const vec_t &ydat, 
		   const mat_t &xpred, vec_t &parms, 
		   mat_t &covar, double &chi2) {

    if (npar>ndat) {
      O2SCL_ERR2("Insufficient data points in ",
		 "fit_linear::fit_linear().",exc_einval);
    }

    if (npar!=size_par || ndat!=size_dat) {
      A.resize(ndat,npar);
      Q.resize(npar,npar);
      QSI.resize(npar,npar);
      S.resize(npar);
      t.resize(ndat);
      xt.resize(npar);
      D.resize(npar);
      size_par=npar;
      size_dat=ndat;
    }

    // necessary?
    for(size_t i=0;i<npar;i++) {
      xt[i]=0.0;
      D[i]=0.0;
    }

    for(size_t i=0;i<ndat;i++) {
      for(size_t j=0;j<npar;j++) {
	A(i,j)=xpred(i,j);
      }
    }

    if (column_scaling) {
      o2scl_linalg::balance_columns(ndat,npar,A,D);
    } else {
      for(size_t i=0;i<npar;i++) {
	D[i]=1.0;
      }
    }

    // [GSL] Decompose A into U S Q^T
    fit_svd(ndat,npar);
    
    // [GSL] Solve y = A c for c
    o2scl_cblas::dgemv
      (o2scl_cblas::o2cblas_RowMajor,
       o2scl_cblas::o2cblas_Trans,ndat,npar,1.0,A,ydat,0.0,xt);

    // [GSL] Scale the matrix Q, Q' = Q S^-1
    for(size_t i=0;i<npar;i++) {
      for(size_t j=0;j<npar;j++) {
	QSI(i,j)=Q(i,j);
      }
    }

    double alpha0=S[0];
    size_t p_eff=0;

    for(size_t j=0;j<npar;j++) {
      double alpha=S[j];
      if (alpha<=tol*alpha0) {
	alpha=0.0;
      } else {
	alpha=1.0/alpha;
	p_eff++;
      }
      o2scl_cblas::dscal_subcol(QSI,0,j,npar,alpha);
    }

    rank=p_eff;
    for(size_t i=0;i<npar;i++) {
      parms[i]=0.0;
    }
    
    o2scl_cblas::dgemv
      (o2scl_cblas::o2cblas_RowMajor,
       o2scl_cblas::o2cblas_NoTrans,npar,npar,
       1.0,QSI,xt,0.0,parms);
    
    // [GSL] Unscale the balancing factors

    for(size_t i=0;i<npar;i++) {
      parms[i]/=D[i];
    }

    // [GSL] Compute chi-squared from residual, r = y - X c

    double s2=0.0, r2=0.0;
    for(size_t i=0;i<ndat;i++) {
      double yi=ydat[i];
      double y_est=o2scl_cblas::ddot_subrow(npar,xpred,i,0,parms);
      double ri=yi-y_est;
      r2+=ri*ri;
    }
    s2=r2/(ndat-p_eff);
    chi2=r2;

    // [GSL] Form variance-covariance matrix, cov = s2 * (Q S^-1) (Q S^-1)^T

    for(size_t i=0;i<npar;i++) {
      double d_i=D[i];
      for(size_t j=i;j<npar;j++) {
	double d_j=D[j];
	double s=0.0;
	for(size_t k=0;k<npar;k++) {
	  s+=QSI(i,k)*QSI(j,k);
	}
	covar(i,j)=s*s2/(d_i*d_j);
	covar(j,i)=s*s2/(d_i*d_j);
      }
    }
    
    return;
  }
    
  /// Return string denoting type ("fit_linear")
  virtual const char *type() { return "fit_linear"; }
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
