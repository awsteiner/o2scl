/*
  ───────────────────────────────────────────────────────────────────
  
  Copyright (C) 2006-2025, Andrew W. Steiner
  
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

  ───────────────────────────────────────────────────────────────────
*/
/** \file linear_solver.h
    \brief File for linear solvers
*/
#ifndef O2SCL_LINEAR_SOLVER_H
#define O2SCL_LINEAR_SOLVER_H

#include <gsl/gsl_linalg.h>

#include <o2scl/permutation.h>
#include <o2scl/lu.h>
#include <o2scl/qr.h>
#include <o2scl/hh.h>
#include <o2scl/set_cuda.h>

namespace o2scl_linalg {

  /** \brief A generic solver for the linear system \f$ A x = b \f$
      [abstract base]

      A generic solver for dense linear systems.

      \verbatim embed:rst
      Those writing production level code should consider calling
      LAPACK directly as described in the
      :ref:`Linear Algebra` section of the User's Guide.
      \endverbatim

      \comment
      AWS, 10/28/24, actually this might be a good idea because
      it also has an analytic inverse and it can demonstrate
      multiprecision.
      
      \future The test code uses a Hilbert matrix, which is known
      to be ill-conditioned, especially for the larger sizes. This
      could be changed.
      \endcomment
  */
  template<class vec_t=boost::numeric::ublas::vector<double>, 
           class mat_t=boost::numeric::ublas::matrix<double>> 
  class linear_solver {
    
  public:
    
    virtual ~linear_solver() {}
    
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const mat_t &a, const vec_t &b,
                       vec_t &x)=0;
    
    /** \brief Solve square linear system \f$ A x = b \f$ of size \c n,
        possibly destroying \c a and \c b
    */
    virtual void solve_dest(size_t n, mat_t &a, vec_t &b,
                            vec_t &x) {
      solve(n,a,b,x);
      return;
    }
    
  };

  /** \brief Generic linear solver using LU decomposition
   */
  template <class vec_t=boost::numeric::ublas::vector<double>, 
            class mat_t=boost::numeric::ublas::matrix<double>,
            class fp_t=double> 
  class linear_solver_LU : public linear_solver<vec_t, mat_t> {
    
  public:
    
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const mat_t &A, 
                       const vec_t &b, vec_t &x) {
      mat_t A2=A;
      vec_t b2=b;
      this->solve_dest(n,A2,b2,x);
      return;
    };
    
    /** \brief Solve square linear system \f$ A x = b \f$ of size \c n,
        possibly destroying \c a and \c b
    */
    virtual void solve_dest(size_t n, mat_t &A, 
                            vec_t &b, vec_t &x) {
      int sig;
      o2scl::permutation p(n);
      LU_decomp<mat_t,fp_t>(n,A,p,sig);
      LU_solve<mat_t,vec_t,vec_t>(n,A,p,b,x);
      return;
    };
    
    virtual ~linear_solver_LU() {}
    
  };

  /** \brief Generic linear solver using QR decomposition
   */
  template <class vec_t=boost::numeric::ublas::vector<double>, 
            class mat_t=boost::numeric::ublas::matrix<double>,
            class fp_t=double> 
  class linear_solver_QR : public linear_solver<vec_t, mat_t> {
    
  public:
    
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const mat_t &A,
                       const vec_t &b, vec_t &x) {
      mat_t A2=A;
      vec_t b2=b;
      solve_dest(n,A2,b2,x);
      return;
    };
    
    /** \brief Solve square linear system \f$ A x = b \f$ of size \c n,
        possibly destroying \c a and \c b
    */
    virtual void solve_dest(size_t n, mat_t &A,
                            vec_t &b, vec_t &x) {
      boost::numeric::ublas::vector<fp_t> tau(n);
      QR_decomp<mat_t,vec_t,fp_t>(n,n,A,tau);
      QR_solve<mat_t,vec_t,vec_t,vec_t,fp_t>(n,A,tau,b,x);
      return;
    };
    
    virtual ~linear_solver_QR() {}
    
  };
  
  /** \brief Generic Householder linear solver 
   */
  template <class vec_t=boost::numeric::ublas::vector<double>, 
            class mat_t=boost::numeric::ublas::matrix<double> > 
  class linear_solver_HH : 
    public linear_solver<vec_t, mat_t> {
    
  public:
    
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const mat_t &A, const vec_t &b, vec_t &x) {
      mat_t A2=A;
      vec_t b2=b;
      solve_dest(n,A2,b2,x);
      return;
    };
    
    /** \brief Solve square linear system \f$ A x = b \f$ of size \c n,
        possibly destroying \c a and \c b
    */
    virtual void solve_dest(size_t n, mat_t &A, vec_t &b,
                            vec_t &x) {
      HH_solve(n,A,b,x);
      return;
    };
    
    virtual ~linear_solver_HH() {}
    
  };

  // End of namespace o2scl_linalg
}

#if defined (O2SCL_COND_FLAG) || defined (DOXYGEN)
#if defined (O2SCL_SET_ARMA) || defined (DOXYGEN)
#include <armadillo>
namespace o2scl_linalg {
  /** \brief Armadillo linear solver 

      This class is only defined if Armadillo support was enabled
      during installation
  */
  template<class arma_vec_t, class arma_mat_t> class linear_solver_arma : 
    public linear_solver<arma_vec_t,arma_mat_t> {
  public:
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const arma_mat_t &A,
                       const arma_vec_t &b, arma_vec_t &x) {
      x=arma::solve(A,b);
      return;
    }
  };
}
#endif

#if defined (O2SCL_SET_EIGEN) || defined (DOXYGEN)
#include <eigen3/Eigen/Dense>
namespace o2scl_linalg {
  /** \brief Eigen linear solver using QR decomposition with 
      column pivoting

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
  class linear_solver_eigen_houseQR : 
    public linear_solver<eigen_vec_t,eigen_mat_t> {
  public:
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const eigen_mat_t &A,
                       const eigen_vec_t &b, eigen_vec_t &x) {
      x=A.householderQr().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using QR decomposition with 
      column pivoting

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
  class linear_solver_eigen_colQR : 
    public linear_solver<eigen_vec_t,eigen_mat_t> {
  public:
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const eigen_mat_t &A,
                       const eigen_vec_t &b, eigen_vec_t &x) {
      x=A.colPivHouseholderQr().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using QR decomposition with 
      full pivoting

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
  class linear_solver_eigen_fullQR : 
    public linear_solver<eigen_vec_t,eigen_mat_t> {
  public:
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const eigen_mat_t &A,
                       const eigen_vec_t &b, eigen_vec_t &x) {
      x=A.fullPivHouseholderQr().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using LU decomposition with 
      partial pivoting
      
      This requires the matrix \c A to be invertible.

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
  class linear_solver_eigen_partLU : 
    public linear_solver<eigen_vec_t,eigen_mat_t> {
  public:
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const eigen_mat_t &A,
                       const eigen_vec_t &b, eigen_vec_t &x) {
      x=A.partialPivLu().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using LU decomposition with 
      full pivoting

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
  class linear_solver_eigen_fullLU : 
    public linear_solver<eigen_vec_t,eigen_mat_t> {
  public:
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const eigen_mat_t &A,
                       const eigen_vec_t &b, eigen_vec_t &x) {
      x=A.fullPivLu().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using LLT decomposition with 
      full pivoting
      
      This requires the matrix \c A to be positive definite.

      This class is only defined if Eigen support was enabled during
      installation.

  */
  template<class eigen_vec_t, class eigen_mat_t>
  class linear_solver_eigen_LLT : 
    public linear_solver<eigen_vec_t,eigen_mat_t> {
  public:
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const eigen_mat_t &A,
                       const eigen_vec_t &b, eigen_vec_t &x) {
      x=A.llt().solve(b);
      return;
    }
  };
  
  /** \brief Eigen linear solver using LDLT decomposition with 
      full pivoting
      
      This requires the matrix \c A to be positive or negative
      semidefinite.
      
      This class is only defined if Eigen support was enabled during
      installation.
  */
  template<class eigen_vec_t, class eigen_mat_t>
  class linear_solver_eigen_LDLT : 
    public linear_solver<eigen_vec_t,eigen_mat_t> {
  public:
    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const eigen_mat_t &A,
                       const eigen_vec_t &b, eigen_vec_t &x) {
      x=A.ldlt().solve(b);
      return;
    }
  };

  
}

// end of #if defined (O2SCL_SET_EIGEN) || defined (DOXYGEN)
#endif

#if defined (O2SCL_SET_CUDA) || defined (DOXYGEN)
#include <o2scl/solve_cuda.h>

namespace o2scl_linalg {
  
  /** \brief CUDA linear solver using LU decomposition (partial
      pivoting)

      This class is only defined if CUDA support was enabled during
      installation.

  */
  class linear_solver_LU_cuda :
    public linear_solver_LU_cuda_base,
    linear_solver<std::vector<double>,std::vector<double> > {

  public:

    /// Solve square linear system \f$ A x = b \f$ of size \c n
    virtual void solve(size_t n, const std::vector<double> &A,
                       const std::vector<double> &B,
                       std::vector<double> &x) {
      int ret=this->solve_base((int)n,A,B,x);
      if (ret!=0) {
        std::string errs=((std::string)"Error number ")+
          o2scl::itos(ret)+" in linear_solver_LU_cuda::solve().";
        O2SCL_ERR(errs.c_str(),o2scl::exc_einval);
      }
      return;
    }
    
  };

}

#endif
    
#else

#include <o2scl/linear_special.h>

// end of #if defined (O2SCL_COND_FLAG) || defined (DOXYGEN)
#endif

#endif
