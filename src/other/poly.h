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
#ifndef O2SCL_POLY_H
#define O2SCL_POLY_H

/** \file poly.h
    \brief Classes for solving polynomials

    \warning 
    One must be careful about using pow() in functions using
    complex<double> since pow(((complex<double>)0.0),3.0) returns
    (nan,nan). Instead, we should use pow(((complex<double>)0.0),3)
    which takes an integer for the second argument. The sqrt()
    function, always succeeds i.e. sqrt(((complex<double>)0.0))=0.0

    \future The quartics are tested only for a4=1, which should
    probably be generalized.
*/

#include <iostream>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_poly.h>
#include <o2scl/constants.h>
#include <o2scl/err_hnd.h>
#include <o2scl/misc.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/root_brent_gsl.h>

namespace o2scl {

  /** \brief Base class for a polynomial with real coefficients
   */
  template<class fp_t=double> class poly_real_base {
    
  protected:

    /// Solver for root polishing
    o2scl::root_brent_gsl<std::function<fp_t(fp_t)>,fp_t> rbg;
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    poly_real_base() {
      // Flip this switch because we will manually check
      // the return value of the solver and the accuracy of
      // the root
      rbg.err_nonconv=false;
    }
    
    virtual ~poly_real_base() {}
    
    /** \brief Polish the roots
     */
    template<class vec_t=std::vector<fp_t> >
    int polish_r_arr(int n, const vec_t &co, 
                     int nr, vec_t &ro) {
      
      std::function<fp_t(fp_t)> f=std::bind
        (std::mem_fn<fp_t(fp_t,int,const vec_t &)>
         (&poly_real_base<fp_t>::polish_fun<vec_t>),
         this,std::placeholders::_1,n,std::cref(co));

      int verbose=0;
      //rbg.verbose=2;
      rbg.tol_abs=1.0e-12;
      rbg.tol_rel=1.0e-12;

      // Sort the roots first
      o2scl::vector_sort<vec_t,fp_t>(nr,ro);

      if (verbose>1) {
        std::cout.precision(14);
        std::cout << "Sorted roots: ";
        o2scl::vector_out(std::cout,nr,ro,true);
      }
      
      fp_t low, high;
      for(int i=0;i<nr;i++) {
        fp_t f_ro=polish_fun(ro[i],n,co);
        if (verbose>1) {
          std::cout << "Root: " << i << " "
                    << ro[i] << " value: " << f_ro << std::endl;
        }
        if (fabs(f_ro)>1.0e-12) {
          if (i==0 && i==nr-1) {
            if (verbose>1) {
              std::cout << "One root: " << std::endl;
            }
            // If there's only one root, then polish it
            low=ro[0];
            int ret=rbg.solve(low,f);
            if (verbose>1) {
              std::cout << "ret 1: " << ret << std::endl;
            }
            if (ret==0) {
              fp_t f_low=polish_fun(low,n,co);
              if (verbose>1) {
                std::cout << "before, after: " << f_ro << " " << f_low
                          << std::endl;
              }
              if (fabs(f_low)<fabs(f_ro)) {
                // Only accept the solution if the solver succeeded
                // and the final root is better
                ro[0]=low;
              }
            }
          } else if (i==0 && ro[0]==ro[1]) {
            if (verbose>1) {
              std::cout << "0 == 1" << std::endl;
            }
            // If there is more than one root, and the first
            // root is identical to the second
            int ret=rbg.solve(low,f);
            if (verbose>1) {
              std::cout << "ret 2: " << ret << std::endl;
            }
            if (ret==0) {
              fp_t f_low=polish_fun(low,n,co);
              if (verbose>1) {
                std::cout << "before, after: " << f_ro << " " << f_low
                          << std::endl;
              }
              if (fabs(f_low)<fabs(f_ro)) {
                // Only accept the solution if the solver succeeded
                // and the final root is better
                ro[0]=low;
              }
            }
            rbg.solve(ro[0],f);
          } else if (i==nr-1 && ro[nr-1]==ro[nr-2]) {
            if (verbose>1) {
              std::cout << "nr-1 == nr-2" << std::endl;
            }
            // If the last and second-to-last roots are identical
            low=ro[nr-1];
            int ret=rbg.solve(low,f);
            if (verbose>1) {
              std::cout << "ret 3: " << ret << std::endl;
            }
            if (ret==0) {
              fp_t f_low=polish_fun(low,n,co);
              if (verbose>1) {
                std::cout << "before, after: " << f_ro << " " << f_low
                          << std::endl;
              }
              if (fabs(f_low)<fabs(f_ro)) {
                // Only accept the solution if the solver succeeded
                // and the final root is better
                ro[nr-1]=low;
              }
            }
          } else {
            if (i>0 && ro[i]==ro[i-1]) {
              if (i<nr-1 && ro[i]==ro[i+1]) {
                if (verbose>1) {
                  std::cout << "i-1 == i == i+1" << std::endl;
                }
                // The root is equal to the previous and the next root,
                // so use the solver to automatically bracket
                low=ro[i];
                int ret=rbg.solve(low,f);
                if (verbose>1) {
                  std::cout << "ret 4: " << ret << std::endl;
                }
                if (ret==0) {
                  fp_t f_low=polish_fun(low,n,co);
                  if (verbose>1) {
                    std::cout << "before, after: " << f_ro << " " << f_low
                              << std::endl;
                  }
                  if (fabs(f_low)<fabs(f_ro)) {
                    // Only accept the solution if the solver succeeded
                    // and the final root is better
                    ro[i]=low;
                  }
                }
              } else {
                if (verbose>1) {
                  std::cout << "i-1 == i" << std::endl;
                }
                // The root is equal to the previous one, so try to
                // bracket
                low=ro[i];
                high=(ro[i]+ro[i+1])/2;
                fp_t f_high=polish_fun(high,n,co);
                int ret;
                if (f_ro*f_high<0.0) {
                  ret=rbg.solve_bkt(low,high,f);
                  if (verbose>1) {
                    std::cout << "ret 5: " << ret << std::endl;
                  }
                } else {
                  // If bracketing fails, use the solver to automatically
                  // bracket
                  ret=rbg.solve(low,f);
                  if (verbose>1) {
                    std::cout << "ret 6: " << ret << std::endl;
                  }
                }
                if (ret==0) {
                  fp_t f_low=polish_fun(low,n,co);
                  if (verbose>1) {
                    std::cout << "before, after: " << f_ro << " " << f_low
                              << std::endl;
                  }
                  if (fabs(f_low)<fabs(f_ro)) {
                    // Only accept the solution if the solver succeeded
                    // and the final root is better
                    ro[i]=low;
                  }
                }
              }
            } else if (i<nr-1 && ro[i]==ro[i+1]) {
              if (verbose>1) {
                std::cout << "i == i+1" << std::endl;
              }
              // The root is equal to the next one, so try to
              // bracket
              low=(ro[i]+ro[i-1])/2;
              high=ro[i];
              fp_t f_low=polish_fun(low,n,co);
              int ret;
              if (f_ro*f_low<0.0) {
                ret=rbg.solve_bkt(low,high,f);
                if (verbose>1) {
                  std::cout << "ret 7: " << ret << std::endl;
                }
              } else {
                // If bracketing fails, use the solver to automatically
                // bracket
                ret=rbg.solve(high,f);
                if (verbose>1) {
                  std::cout << "ret 8: " << ret << std::endl;
                }
                low=high;
              }
              if (ret==0) {
                f_low=polish_fun(low,n,co);
                if (verbose>1) {
                  std::cout << "before, after: " << f_ro << " " << f_low
                            << std::endl;
                }
                if (fabs(f_low)<fabs(f_ro)) {
                  // Only accept the solution if the solver succeeded
                  // and the final root is better
                  ro[i]=low;
                }
              }
            } else {
              if (verbose>1) {
                std::cout << "No repeats." << std::endl;
              }
              // Otherwise, if we have multiple but not repeated roots,
              // try to bracket
              if (i==0) {
                low=ro[0]-(ro[1]-ro[0]);
                high=ro[1]-(ro[1]-ro[0])/2.0;
              } else if (i==nr-1) {
                low=ro[nr-2]+(ro[nr-1]-ro[nr-2])/2.0;
                high=ro[nr-1]+(ro[nr-1]-ro[nr-2]);
              } else {
                low=ro[i-1]+(ro[i]-ro[i-1])/2.0;
                high=ro[i+1]-(ro[i+1]-ro[i])/2.0;
              }
              if (verbose>1) {
                std::cout << "low ro high: " << low << " " << ro[i] << " "
                          << high << std::endl;
              }
              fp_t f_low=polish_fun(low,n,co);
              fp_t f_high=polish_fun(high,n,co);
              if (verbose>1) {
                std::cout << "f_low f_ro f_high: " << f_low << " "
                          << f_ro << " "
                          << f_high << std::endl;
              }
              int ret;
              if (f_low*f_ro<0) {
                fp_t tmp=ro[i];
                if (verbose>1) {
                  std::cout << "Case 1: " << low << " " << tmp << std::endl;
                  std::cout << f_low << " " << f_ro << std::endl;
                }
                ret=rbg.solve_bkt(low,tmp,f);
                if (verbose>1) {
                  std::cout << "ret 9: " << ret << std::endl;
                }
              } else if (f_high*f_ro<0) {
                fp_t tmp=ro[i];
                if (verbose>1) {
                  std::cout << "Case 2: " << tmp << " " << high << std::endl;
                }
                ret=rbg.solve_bkt(tmp,high,f);
                if (verbose>1) {
                  std::cout << "ret 10: " << ret << std::endl;
                }
                low=tmp;
              } else {
                // If bracketing fails, use the solver to try to automatically
                // bracket
                low=ro[i];
                if (verbose>1) {
                  std::cout << "Case 3: " << low << std::endl;
                }
                ret=rbg.solve(low,f);
                if (verbose>1) {
                  std::cout << "ret 11: " << ret << std::endl;
                }
              }
              if (ret==0) {
                f_low=polish_fun(low,n,co);
                if (verbose>1) {
                  std::cout << "before, after: " << f_ro << " " << f_low
                            << std::endl;
                }
                if (fabs(f_low)<fabs(f_ro)) {
                  // Only accept the solution if the solver succeeded
                  // and the final root is better
                  ro[i]=low;
                }
              }
            }
          }
          if (verbose>1) {
            f_ro=polish_fun(ro[i],n,co);
            std::cout << "Updated root: " << i << " "
                      << ro[i] << " value: " << f_ro << std::endl;
            char ch;
            std::cin >> ch;
          }
        }
      }
      
      return 0;
    }

    /** \brief Desc
     */
    template<class vec_t=std::vector<fp_t> >
    fp_t polish_fun(fp_t x, int n, const vec_t &co) {
      fp_t ret=co[0];
      for(int i=0;i<n;i++) {
        ret=co[i+1]+x*ret;
      }
      return ret;
    }
  };
  
  /** \brief Solve a quadratic polynomial with real coefficients and 
      real roots [abstract base]
  */
  template<class fp_t=double> class quadratic_real :
    public poly_real_base<fp_t> {
    
  public:

    virtual ~quadratic_real() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 

	If two real solutions exist (i.e. if the discriminant is
        non-zero) then the real roots are placed in \c x1 and \c x2
        and the number 2 is returned. Otherwise, 0 is returned and x1
        and x2 are unmodified.
    */
    virtual int solve_r(const fp_t a2, const fp_t b2, const fp_t c2, 
			fp_t &x1, fp_t &x2)=0;

    /** \brief Compute the quadratic discriminant, \f$ b^2-4ac \f$

        If the discriminant is positive, the quadratic has two real
        roots, if it is zero, it has a double real root, and if it is
        negative then the quadratic has complex conjugate roots.
     */
    virtual fp_t disc2_r(fp_t a2, fp_t b2, fp_t c2) {
      return b2*b2-4.0*a2*c2;
    }
    
    /// Return a string denoting the type ("quadratic_real")
    const char *type() { return "quadratic_real"; }
  };

  /** \brief Solve a quadratic polynomial with real coefficients to
      get the real roots (GSL)

      This is just a wrapper for gsl_poly_complex_solve_quadratic().
  */
  class quadratic_real_gsl : public quadratic_real<> {
    
  public:

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const double a2, const double b2, const double c2, 
			double &x1, double &x2);
    
  };
  
  /** \brief Solve a quadratic polynomial with real coefficients and 
      complex roots [abstract base]

      This class is designed to be identical to the GSL code in
      gsl_poly_complex_solve_quadratic(), but generalized to generic
      floating point types.
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quadratic_real_coeff : public quadratic_real<fp_t> {

  public:

    virtual ~quadratic_real_coeff() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const fp_t a2, const fp_t b2, const fp_t c2, 
			fp_t &x1, fp_t &x2) {
      if (a2==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_real_coeff::solve_r().",
           exc_einval);
      }
      cx_t r1, r2;
      int ret=solve_rc(a2,b2,c2,r1,r2);
      if (ret==2) {
        x1=r1.real();
        x2=r2.real();
      }
      return ret;
    }

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_rc(const fp_t a2, const fp_t b2, const fp_t c2, 
			 cx_t &x1, cx_t &x2)=0;
    
    /** \brief Test \f$ n^3 \f$ quadratics with coefficients near a 
        zero discriminant
     */
    size_t test_coeffs_zero_disc(fp_t alpha, fp_t &s1, fp_t &s2, fp_t &m1,
                                 fp_t &m2, size_t n=40) {

      o2scl::gen_test_number<> ga, gb, gc;
      
      size_t count=0;
      
      // First pick random coefficients
      for(size_t j1=0;j1<n;j1++) {
        fp_t ca=ga.gen();
        gb.reset();
        for(size_t j2=0;j2<n;j2++) {
          fp_t cb=gb.gen()*alpha;
          gc.reset();
          for(size_t j3=0;j3<n;j3++) {
            
            // Ensure that several quadratics near b^2=4*a*c are tested
            fp_t cc=cb*cb/4.0/ca+gc.gen();
            
            // Ensure there is a solution
            if (fabs(ca)!=0.0) {

              cx_t cr1, cr2;
              solve_rc(ca,cb,cc,cr1,cr2);
              
              cx_t cbp=-(cr1+cr2)*ca;
              cx_t ccp=(cr1*cr2)*ca;
              
              cx_t czo1=(ca*cr1+cb)*cr1+cc;
              cx_t czo2=(ca*cr2+cb)*cr2+cc;
              fp_t q1=sqrt(fabs(cb-cbp.real())*fabs(cb-cbp.real())+
                           fabs(cc-ccp.real())*fabs(cc-ccp.real()));
              fp_t q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2));
              
              s1+=q1;
              if (q1>m1) m1=q1;
              s2+=q2;
              if (q2>m2) m2=q2;
              count++;
              
            }
          }
        }
      }
      
      return count;
    }
      
    /** \brief Test \f$ n^2 \f$ quadratics with complex roots
     */
    size_t test_complex_roots(fp_t &s1, fp_t &s2, fp_t &m1, fp_t &m2,
                              size_t n=80) {

      gen_test_number<> gd, ge;

      size_t count=0;
      
      // Next, pick random roots which are complex conjugates
      for(size_t j1=0;j1<n;j1++) {
        cx_t cr1, cr2;
        cr1.real(gd.gen());
        ge.reset();
        for(size_t j2=0;j2<n;j2++) {
          cr1.imag(ge.gen());
          cr2.real(cr1.real());
          cr2.imag(-cr1.imag());
          
          fp_t ca=1.0;
          fp_t cb=(-cr1-cr2).real();
          fp_t cc=(cr1*cr2).real();

          cx_t cr1p, cr2p;
          solve_rc(ca,cb,cc,cr1p,cr2p);
          
          // If the roots are flipped
          if (fabs(cr1.imag()-cr2p.imag())<fabs(cr1.imag()-cr1p.imag())) {
            cx_t temp=cr1p;
            cr1p=cr2p;
            cr2p=temp;
          } 
          
          fp_t q1=sqrt(fabs(cr1.real()-cr1p.real())*
                       fabs(cr1.real()-cr1p.real())+
                       fabs(cr2.real()-cr2p.real())*
                       fabs(cr2.real()-cr2p.real()));
          fp_t q2=sqrt(fabs(cr1.imag()-cr1p.imag())*
                       fabs(cr1.imag()-cr1p.imag())+
                       fabs(cr2.imag()-cr2p.imag())*
                       fabs(cr2.imag()-cr2p.imag()));
          
          s1+=q1;
          if (q1>m1) m1=q1;
          s2+=q2;
          if (q2>m2) m2=q2;
          count++;
          
        }
      }
      
      return count;
    }
      
    /** \brief Test \f$ n^2 \f$ quadratics with real roots
     */
    size_t test_real_roots(fp_t &s1, fp_t &s2, fp_t &m1,
                           fp_t &m2, size_t &wrong_ret, size_t n=80) {
      
      o2scl::gen_test_number<> gd, ge;

      size_t count=0;
      
      // Next, pick random roots which are both real
      for(size_t j1=0;j1<n;j1++) {
        cx_t cr1, cr2;
        cr1.real(gd.gen());
        ge.reset();
        for(size_t j2=0;j2<n;j2++) {
          cr2.real(ge.gen());
          cr1.imag(0.0);
          cr2.imag(0.0);
          
          fp_t ca=1.0;
          fp_t cb=(-cr1-cr2).real();
          fp_t cc=(cr1*cr2).real();

          cx_t cr1p, cr2p;
          int ret=solve_rc(ca,cb,cc,cr1p,cr2p);
          if (ret!=2) {
            wrong_ret+=1;
          }
          
          // If the roots are flipped
          if (fabs(cr1.real()-cr2p.real())<fabs(cr1.real()-cr1p.real())) {
            cx_t temp=cr1p;
            cr1p=cr2p;
            cr2p=temp;
          } 
          
          fp_t q1=sqrt(fabs(cr1.real()-cr1p.real())*
                       fabs(cr1.real()-cr1p.real())+
                       fabs(cr2.real()-cr2p.real())*
                       fabs(cr2.real()-cr2p.real()));
          fp_t q2=sqrt(fabs(cr1.imag()-cr1p.imag())*
                       fabs(cr1.imag()-cr1p.imag())+
                       fabs(cr2.imag()-cr2p.imag())*
                       fabs(cr2.imag()-cr2p.imag()));
          
          s1+=q1;
          if (q1>m1) m1=q1;
          s2+=q2;
          if (q2>m2) m2=q2;
          count++;
          
        }
      }
      
      return count;
    }
    
    /// Return a string denoting the type ("quadratic_real_coeff")
    const char *type() { return "quadratic_real_coeff"; }
  };

  /** \brief Solve a quadratic polynomial with complex coefficients and 
      complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quadratic_complex :  public quadratic_real_coeff<fp_t,cx_t> {
    
  public:

    virtual ~quadratic_complex() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ .
    */
    virtual int solve_r(const fp_t a2, const fp_t b2, const fp_t c2, 
			fp_t &x1, fp_t &x2) {
      if (a2==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_complex::solve_r().",
           exc_einval);
      }

      cx_t r1, r2;
      int ret=solve_c(a2,b2,c2,r1,r2);
      if (this->disc2_r(a2,b2,c2)>=0.0) {
        x1=r1.real();
        x2=r2.real();
        return 2;
      }
      return 0;
    }

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_rc(const fp_t a2, const fp_t b2, const fp_t c2, 
			 cx_t &x1, cx_t &x2) {
      if (a2==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_complex::solve_rc().",
           exc_einval);
      }
      
      int ret=solve_c(a2,b2,c2,x1,x2);
      if (this->disc2_r(a2,b2,c2)>=0.0) {
        return 2;
      }
      return 0;
    }

    /** \brief Solves the complex polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_c(const cx_t a2, const cx_t b2, 
			const cx_t c2, cx_t &x1, cx_t &x2)=0;

    /** \brief Test \f$ n^6 \f$ quadratics with complex roots
     */
    size_t test_complex_coeffs(fp_t &s1, fp_t &s2, fp_t &m1, fp_t &m2,
                              size_t n=10) {

      gen_test_number<> ga, gb, gc, gd, ge, gf;
      
      fp_t rca, rcb, rcc, rcd, rce, rcf;
      cx_t i(0.0,1.0);

      size_t count=0;
      
      // Test random coefficients
      for(size_t j1=0;j1<n;j1++) {
        rca=ga.gen();
        gb.reset();
        for(size_t j2=0;j2<n;j2++) {
          rcb=gb.gen();
          gc.reset();
          for(size_t j3=0;j3<n;j3++) {
            rcc=gc.gen();
            gd.reset();
            for(size_t j4=0;j4<n;j4++) {
              rcd=gd.gen();
              ge.reset();
              for(size_t j5=0;j5<n;j5++) {
                rce=ge.gen();
                gf.reset();
                for(size_t j6=0;j6<n;j6++) {
                  rcf=gf.gen();
                  
                  if (fabs(rca)>0.0 || fabs(rcb)>0.0) {
                    cx_t ca=rca+i*rcb;
                    cx_t cb=rcc+i*rcd;
                    cx_t cc=rce+i*rcf;

                    cx_t cr1, cr2;
                    solve_c(ca,cb,cc,cr1,cr2);
                    
                    cx_t cbp=-(cr1+cr2)*ca;
                    cx_t ccp=(cr1*cr2)*ca;
                    
                    cx_t czo1=(ca*cr1+cb)*cr1+cc;
                    cx_t czo2=(ca*cr2+cb)*cr2+cc;
                    fp_t q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
                                 abs(cc-ccp)*abs(cc-ccp));
                    fp_t q2=sqrt(abs(czo1)*abs(czo1)+
                                 abs(czo2)*abs(czo2));
                    
                    s1+=q1;
                    if (q1>m1) m1=q1;
                    s2+=q2;
                    if (q2>m2) m2=q2;
                    
                    count++;
                    
                  }
                  
                }
              }
            }
          }
        }
      }

      return count;
    }

    /** \brief Test \f$ n^4 \f$ quadratics with complex roots
     */
    size_t test_complex_roots(fp_t &s1, fp_t &s2, fp_t &m1, fp_t &m2,
                              size_t n=20) {
      
      gen_test_number<> gg, gh, gi, gj;
      
      fp_t rca, rcb, rcc, rcd, rce, rcf;
      cx_t i(0.0,1.0);

      size_t count=0;
      
      // Test random roots
      for(size_t j1=0;j1<20;j1++) {
        rca=gg.gen();
        gh.reset();
        for(size_t j2=0;j2<20;j2++) {
          rcb=gh.gen();
          gi.reset();
          for(size_t j3=0;j3<20;j3++) {
            rcc=gi.gen();
            gj.reset();
            for(size_t j4=0;j4<20;j4++) {
              rcd=gj.gen();
              
              if (fabs(rca)>0.0 || fabs(rcb)>0.0) {
                cx_t cr1p=rca+i*rcb;
                cx_t cr2p=rcc+i*rcd;
                
                cx_t ca=1.0;
                cx_t cb=-cr1p-cr2p;
                cx_t cc=cr1p*cr2p;

                cx_t cr1, cr2;
                solve_c(ca,cb,cc,cr1,cr2);
                
                // If the roots are flipped
                if (fabs(cr1.real()-cr2p.real())+
                    fabs(cr1.imag()-cr2p.imag())<
                    fabs(cr1.real()-cr1p.real())+
                    fabs(cr1.imag()-cr1p.imag())) {
                    
                  cx_t temp=cr1;
                  cr1=cr2;
                  cr2=temp;
                } 
                
                fp_t q1=sqrt(fabs(cr1.real()-cr1p.real())*
                             fabs(cr1.real()-cr1p.real())+
                             fabs(cr2.real()-cr2p.real())*
                             fabs(cr2.real()-cr2p.real()));
                fp_t q2=sqrt(fabs(cr1.imag()-cr1p.imag())*
                             fabs(cr1.imag()-cr1p.imag())+
                             fabs(cr2.imag()-cr2p.imag())*
                             fabs(cr2.imag()-cr2p.imag()));
                
                s1+=q1;
                if (q1>m1) m1=q1;
                s2+=q2;
                if (q2>m2) m2=q2;
                count++;
                
              }
              
            }
          }
        }
      }
      
      return count;
    }
    
    /// Return a string denoting the type ("quadratic_complex")
    const char *type() { return "quadratic_complex"; }
  };  

  /** \brief Solve a cubic polynomial with real coefficients and real roots
      [abstract base]
  */
  template<class fp_t=double> class cubic_real :
    public poly_real_base<fp_t> {
  public:

    virtual ~cubic_real() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the three 
	solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_r(const fp_t a3, const fp_t b3, const fp_t c3, 
			const fp_t d3, fp_t &x1, fp_t &x2, 
			fp_t &x3)=0;

    /** \brief Compute the cubic discriminant, 
	\f$ b^2 c^2 - 4 a c^3 - 4 b^3 d - 27 a^2 d^2 + 18 a b c d \f$

        If the discriminant is zero, then all roots qre real and
        at least two are equal (possibly all three are identical).
        If the discriminant is positive, then there are 
        three distinct real roots, and if the discriminant is negative
        then there is one real root and two complex conjugate
        roots. 
     */
    virtual fp_t disc3_r(const fp_t a3, const fp_t b3, const fp_t c3, 
			  const fp_t d3) {
      return b3*b3*c3*c3-4.0*a3*c3*c3*c3-4.0*b3*b3*b3*d3-
	27.0*a3*a3*d3*d3+18.0*a3*b3*c3*d3;
    }

    /// Return a string denoting the type ("cubic_real")
    const char *type() { return "cubic_real"; }
  };

  /** \brief Solve a cubic polynomial with real coefficients and 
      complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_real_coeff : public cubic_real<fp_t> {

  public:

    virtual ~cubic_real_coeff() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the three 
	solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_r(const fp_t a3, const fp_t b3, const fp_t c3, 
			const fp_t d3, fp_t &x1, fp_t &x2, fp_t &x3) {
      
      cx_t r2, r3;
  
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in",
                   "cubic_real_coeff::solve_r().",exc_einval);
      }
      
      int ret=solve_rc(a3,b3,c3,d3,x1,r2,r3);
      if (ret==3) {
        x2=r2.real();
        x3=r3.real();
        return 3;
      }
      return 1;
    }

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the real solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &x1, cx_t &x2,
			 cx_t &x3)=0;

    /** \brief Test \f$ n^4 \f$ cubics with real coefficients
     */
    size_t test_cubic_real_coeffs(fp_t alpha, fp_t &s1, fp_t &s2,
                            fp_t &m1, fp_t &m2, size_t &wrong_ret,
                            size_t n=16) {
      
      size_t count=0;
      
      gen_test_number<> ga, gb, gc, gd;
      
      for(size_t j1=0;j1<n;j1++) {
        fp_t ca=ga.gen()*alpha;
        gb.reset();
        for(size_t j2=0;j2<n;j2++) {
          fp_t cb=gb.gen();
          gc.reset();
          for(size_t j3=0;j3<n;j3++) {
            fp_t cc=gc.gen()*alpha;
            gd.reset();
            for(size_t j4=0;j4<n;j4++) {
              fp_t cd=gd.gen();
              
              if (fabs(ca)>0.0) {

                fp_t cr1;
                cx_t cr2, cr3;
                int ret=solve_rc(ca,cb,cc,cd,cr1,cr2,cr3);
                fp_t disc=this->disc3_r(ca,cb,cc,cd);
                if ((disc>=0.0 && ret==1) ||
                    (disc<0.0 && ret==3) ||
                    (ret!=1 && ret!=3)) {
                  wrong_ret++;
                }
                
                cx_t cbp=-(cr1+cr2+cr3)*ca;
                cx_t ccp=(cr1*cr2+cr1*cr3+cr2*cr3)*ca;
                cx_t cdp=-(cr1*cr2*cr3)*ca;
                
                fp_t czo1=((ca*cr1+cb)*cr1+cc)*cr1+cd;
                cx_t czo2=((ca*cr2+cb)*cr2+cc)*cr2+cd;
                cx_t czo3=((ca*cr3+cb)*cr3+cc)*cr3+cd;
                
                fp_t q1=sqrt(fabs(cb-cbp.real())*fabs(cb-cbp.real())+
                        fabs(cc-ccp.real())*fabs(cc-ccp.real())+
                        fabs(cd-cdp.real())*fabs(cd-cdp.real()));
                fp_t q2=sqrt(fabs(czo1)*fabs(czo1)+abs(czo2)*abs(czo2)+
                             abs(czo3)*abs(czo3));
                
                s1+=q1;
                if (q1>m1) m1=q1;
                s2+=q2;
                if (q2>m2) m2=q2;
                
                count++;
              }
            }
          }
        }
      }
      
      return count;
    }

    /// Return a string denoting the type ("cubic_real_coeff")
    const char *type() { return "cubic_real_coeff"; }
  };

  /** \brief Solve a cubic polynomial with complex coefficients and 
      complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_complex : public cubic_real_coeff<fp_t,cx_t> {

  public:

    virtual ~cubic_complex() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the three 
	solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_r(const fp_t a3, const fp_t b3, const fp_t c3, 
			const fp_t d3, fp_t &x1, fp_t &x2, fp_t &x3) {
      
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in",
                   "cubic_complex::solve_r().",exc_einval);
      }
      
      cx_t r1,r2,r3;
      int ret=solve_c(a3,b3,c3,d3,r1,r2,r3);
      if (this->disc3_r(a3,b3,c3,d3)>=0) {
        x1=r1.real();
        x2=r2.real();
        x3=r3.real();
        return 3;
      }
      fp_t r1i=fabs(r1.imag());
      fp_t r2i=fabs(r2.imag());
      fp_t r3i=fabs(r3.imag());
      if (r1i<=r2i && r1i<=r3i) {
        x1=r1.real();
      } else if (r2i<=r1i && r2i<=r3i) {
        x1=r2.real();
      } else {
        x1=r3.real();
      }
      return 1;
    }

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the real 
	solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &x1, cx_t &x2,
			 cx_t &x3) {
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in",
                   "cubic_complex::solve_rc().",exc_einval);
      }
      
      cx_t r1,r2,r3;
      
      int ret=solve_c(a3,b3,c3,d3,r1,r2,r3);
      
      if (ret<0) {
        O2SCL_ERR("Function solve_c() failed.",o2scl::exc_einval);
      }
      
      // If all roots are real
      
      if (r1.imag()==0.0 && r2.imag()==0.0 && r3.imag()==0.0) {
        if (r1.real()==0.0 && r2.real()==0.0 && r3.real()==0.0) {
          x1=0;
          x2=0;
          x3=0;
        } else {
          x1=r1.real();
          x2=r2.real();
          x3=r3.real();
        }
        
        return 3;
      }

      // Otherwise, determine which of the three roots is real
      // and put the real root in x1
      
      fp_t s1=fabs(r1.imag());
      fp_t s2=fabs(r2.imag());
      fp_t s3=fabs(r3.imag());
      
      if (s1<s2 && s1<s3) {
        x1=r1.real();
        x2=r2;
        x3=r3;
      } else if (s2<s1 && s2<s3) {
        x1=r2.real();
        x2=r1;
        x3=r3;
      } else {
        x1=r3.real();
        x2=r1;
        x3=r2;
      }
      
      return 1;
    }
    
    /** \brief Solves the complex polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the three complex solutions \f$ x=x_1 \f$ , 
	\f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_c(const cx_t a3, const cx_t b3, 
			const cx_t c3, const cx_t d3, 
			cx_t &x1, cx_t &x2, cx_t &x3)=0;

    /** \brief Test \f$ 2 n^4 \f$ cubics with complex coefficients
     */
    size_t test_complex_coeffs(fp_t &s1, fp_t &s2,
                               fp_t &m1, fp_t &m2, size_t n=9) {
  
      cx_t i(0.0,1.0);
      
      gen_test_number<> ga, gb, gc, gd;
      fp_t rca, rcb, rcc, rcd;

      size_t count=0;
      
      for(size_t it=0;it<2;it++) {
        for(size_t j1=0;j1<n;j1++) {
          rca=ga.gen();
          gb.reset();
          for(size_t j2=0;j2<n;j2++) {
            rcb=gb.gen();
            gc.reset();
            for(size_t j3=0;j3<n;j3++) {
              rcc=gc.gen();
              gd.reset();
              for(size_t j4=0;j4<n;j4++) {
                rcd=gd.gen();

                cx_t ca, cb, cc, cd;
                if (it==0) {
                  ca=rca+i;
                  cb=rcb+i;
                  cc=rcc+i;
                  cd=rcd+i;
                } else {
                  fp_t one=1.0;
                  ca=one+i*rca;
                  cb=one+i*rcb;
                  cc=one+i*rcc;
                  cd=one+i*rcd;
                }
                
                if (fabs(ca.real())>0.0 || fabs(ca.imag())>0.0) {

                  cx_t cr1, cr2, cr3;
                  solve_c(ca,cb,cc,cd,cr1,cr2,cr3);
                  
                  cx_t cbp=-(cr1+cr2+cr3)*ca;
                  cx_t ccp=(cr1*cr2+cr1*cr3+cr2*cr3)*ca;
                  cx_t cdp=-(cr1*cr2*cr3)*ca;
                  
                  cx_t czo1=((ca*cr1+cb)*cr1+cc)*cr1+cd;
                  cx_t czo2=((ca*cr2+cb)*cr2+cc)*cr2+cd;
                  cx_t czo3=((ca*cr3+cb)*cr3+cc)*cr3+cd;
                  
                  fp_t q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
                               abs(cc-ccp)*abs(cc-ccp)+
                               abs(cd-cdp)*abs(cd-cdp));
                  fp_t q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
                               abs(czo3)*abs(czo3));
                  
                  s1+=q1;
                  if (q1>m1) m1=q1;
                  s2+=q2;
                  if (q2>m2) m2=q2;

                  count++;
                }
                
              }
            }
          }
        }
      }
      return count;
    }
    
    /// Return a string denoting the type ("cubic_complex")
    const char *type() { return "cubic_complex"; }
  };

  /** \brief Solve a quartic polynomial with real coefficients and 
      real roots [abstract base]
  */
  template<class fp_t=double> class quartic_real :
    public poly_real_base<fp_t> {

  public:

    virtual ~quartic_real() {}
    
    /** \brief Solves the polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const fp_t a4, const fp_t b4, const fp_t c4, 
			const fp_t d4, const fp_t e4, 
			fp_t &x1, fp_t &x2, 
			fp_t &x3, fp_t &x4)=0;

    /** \brief Compute the quartic discriminant

	The discriminant is zero if and only if at least two roots are
	equal. If the discriminant is non-zero, the discriminant is
	negative if there are two real roots and two complex conjugate
	roots, and it is positive if the roots are either all real or
	all non-real.
    */
    virtual fp_t disc4_r(const fp_t a, const fp_t b, const fp_t c, 
			  const fp_t d, const fp_t e) {
      fp_t a2=a*a;
      fp_t b2=b*b;
      fp_t c2=c*c;
      fp_t d2=d*d;
      fp_t e2=e*e;
      
      fp_t a3=a2*a;
      fp_t b3=b2*b;
      fp_t c3=c2*c;
      fp_t d3=d2*d;
      fp_t e3=e2*e;
      
      fp_t b4=b2*b2;
      fp_t c4=c2*c2;
      fp_t d4=d2*d2;
      
      return 256*a3*e3-192*a2*b*d*e2-128*a2*c2*e2+144*a2*c*d2*e-
        27*a2*d4+144*a*b2*c*e2-6*a*b2*d2*e-80*a*b*c2*d*e+
        18*a*b*c*d3+16*a*c4*e-4*a*c3*d2-27*b4*e2+18*b3*c*d*e-
        4*b3*d3-4*b2*c3*e+b2*c2*d2;
    }

    /** \brief Evaluate diagnostic quantities for a quartic
        
        If disc is negative, then the quartic has two distinct real
        roots and two complex conjugate roots. If disc is positive and
        P and D are both negative, then all four roots are real and
        distinct. If disc is positive and either P or D are positive,
        then there are two pairs of complex conjugate roots.
        If disc is zero, then the polynomial has a multiple root,
        and the following cases hold
        - P<0 and D<0 and disc_0 non-zero: one real double root and
        two real roots
        - D>0 or (P>0 and (D nonzero or R nonzero)) real double root
        and two complex conjugate roots
        - disc_0=0 and D non-zero: triple real root and a real root
        - D=0 and P<0: two real double roots
        - D=0 and P>0 and R=0: two complex conjugate double roots
        - D=0 and disc_0=0, all four roots are equal to -b/4/a. 
        These possibilities correspond to different values of
        \c n_real and \c real_type
        - root_type=0 (n_real=0, n_real=2 or n_real=4) : all real roots 
        are distinct
        - root_type=1 (n_real=2 or n_real4): one double root
        - root_type=2 : two double roots
        - root_type=3 (n_real=4): 
        one triple root and one fourth distinct real root
        - root_type=4 (n_real=4): all real roots equal to -b/4/a
        - root_type=5 (n_real=0): two complex conjugate double roots
            
        Following https://en.wikipedia.org/wiki/Quartic_function .
     */
    virtual void diag_r(const fp_t a, const fp_t b, const fp_t c, 
                        const fp_t d, const fp_t e, fp_t &disc,
                        fp_t &P, fp_t &R, fp_t &disc_0, fp_t &D,
                        size_t &n_real, size_t &root_type) {

      fp_t a2=a*a;
      fp_t b2=b*b;
      fp_t c2=c*c;
      fp_t d2=d*d;
      fp_t e2=e*e;
      
      fp_t a3=a2*a;
      fp_t b3=b2*b;
      fp_t c3=c2*c;
      fp_t d3=d2*d;
      fp_t e3=e2*e;
      
      fp_t b4=b2*b2;
      fp_t c4=c2*c2;
      fp_t d4=d2*d2;

      P=8*a*c-3*b2;
      R=b3+8*d*a2-4*a*b*c;
      disc_0=c2-3*b*d+12*a*e;
      D=64*a3*e-16*a2*c2+16*a*b2*c-16*a2*b*d-3*b4;
      disc=256*a3*e3-192*a2*b*d*e2-128*a2*c2*e2+144*a2*c*d2*e-
        27*a2*d4+144*a*b2*c*e2-6*a*b2*d2*e-80*a*b*c2*d*e+
        18*a*b*c*d3+16*a*c4*e-4*a*c3*d2-27*b4*e2+18*b3*c*d*e-
        4*b3*d3-4*b2*c3*e+b2*c2*d2;

      if (disc<0) {
        n_real=2;
        root_type=0;
      } else if (disc>0) {
        if (P<0 && D<0) {
          n_real=4;
          root_type=0;
        } else {
          n_real=0;
          root_type=0;
        }
      } else {
        if (P<0 && D<0 && disc_0!=0) {
          n_real=4;
          root_type=1;
        } else if (D>0 || (P>0 && (D!=0 || R!=0))) {
          n_real=2;
          root_type=1;
        } else if (disc_0==0 && D!=0) {
          n_real=4;
          root_type=3;
        } else if (D==0 && P<0) {
          n_real=4;
          root_type=2;
        } else if (D==0 && P>0 && R==0) {
          n_real=0;
          root_type=5;
        } else {
          n_real=4;
          root_type=4;
        }
      }
      
      return;
    }
    
    /** \brief Test \f$ n^4 \f$ quartics with real roots
     */
    size_t test_real_roots(fp_t alpha, fp_t &s1, fp_t &m1, size_t n=9) {
      
      size_t count=0;
      
      gen_test_number<> ga, gb, gc, gd;
      
      for(size_t j1=0;j1<n;j1++) {
        fp_t r1=ga.gen();
        gb.reset();
        for(size_t j2=0;j2<n;j2++) {
          fp_t r2=-r1+alpha*gb.gen();
          gc.reset();
          for(size_t j3=0;j3<n;j3++) {
            fp_t r3=gc.gen();
            gd.reset();
            for(size_t j4=0;j4<n;j4++) {
              fp_t r4=-r3+alpha*gd.gen();
              
              fp_t ca=1.0;
              fp_t cb=-(r1+r2+r3+r4);
              fp_t cc=(r1*r2+r1*r3+r2*r3+r1*r4+r2*r4+r3*r4);
              fp_t cd=-(r1*r2*r3+r1*r2*r4+r1*r3*r4+r2*r3*r4);
              fp_t ce=r1*r2*r3*r4;

              fp_t cr1, cr2, cr3, cr4;
              solve_r(ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4);

              if (false) {
                std::vector<fp_t> co={ca,cb,cc,cd,ce};
                std::vector<fp_t> ro={cr1,cr2,cr3,cr4};
                this->polish_r_arr(4,co,4,ro);
              }
              
              fp_t zo1=(((ca*cr1+cb)*cr1+cc)*cr1+cd)*cr1+ce;
              fp_t zo2=(((ca*cr2+cb)*cr2+cc)*cr2+cd)*cr2+ce;
              fp_t zo3=(((ca*cr3+cb)*cr3+cc)*cr3+cd)*cr3+ce;
              fp_t zo4=(((ca*cr4+cb)*cr4+cc)*cr4+cd)*cr4+ce;
              
              fp_t q1=sqrt(zo1*zo1+zo2*zo2+zo3*zo3+zo4*zo4);
              //if (m1>1.0e-12) {
              //std::cout << "Xere: " << q1 << " " << m1 << std::endl;
              //}
              
              s1+=q1;
              if (q1>m1) m1=q1;

              count++;
            }
          }
        }
      }

      return count;
    }
    
    /// Return a string denoting the type ("quartic_real")
    const char *type() { return "quartic_real"; }
  };

  /** \brief Solve a quartic polynomial with real coefficients and 
      complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quartic_real_coeff : public quartic_real<fp_t> {

  public:

    virtual ~quartic_real_coeff() {}

    /** \brief Solves the polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const fp_t a4, const fp_t b4, const fp_t c4, 
			const fp_t d4, const fp_t e4, fp_t &x1, 
			fp_t &x2, fp_t &x3, fp_t &x4) {
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_real_coeff::solve_r().",
           exc_einval);
      }
      
      cx_t r1,r2,r3,r4;
      int ret=solve_rc(a4,b4,c4,d4,e4,r1,r2,r3,r4);
      x1=r1.real();
      x2=r2.real();
      x3=r3.real();
      x4=r4.real();
      return ret;
    }

    /**  \brief Solves the polynomial 
	 \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	 giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	 \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_rc(const fp_t a4, const fp_t b4, const fp_t c4, 
			 const fp_t d4, const fp_t e4, 
			 cx_t &x1, cx_t &x2, 
			 cx_t &x3, 
			 cx_t &x4)=0;

    /** \brief Desc
     */
    size_t test_quartic_real_coeffs(fp_t &s1, fp_t &s2, fp_t &m1,
                          fp_t &m2, size_t n=9) {
      
      cx_t i(0.0,1.0);

      size_t count=0;
      
      gen_test_number<> ga, gb, gc, gd, ge;
      for(size_t j1=0;j1<n;j1++) {
        fp_t ca=ga.gen();
        gb.reset();
        for(size_t j2=0;j2<n;j2++) {
          fp_t cb=gb.gen();
          gc.reset();
          for(size_t j3=0;j3<n;j3++) {
            fp_t cc=gc.gen();
            gd.reset();
            for(size_t j4=0;j4<n;j4++) {
              fp_t cd=gd.gen();
              ge.reset();
              for(size_t j5=0;j5<n;j5++) {
                fp_t ce=ge.gen();
                
                if (fabs(ca)>0.0) {

                  cx_t cr1, cr2, cr3, cr4;
                  solve_rc(ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4);
                  
                  cx_t cbp=-(cr1+cr2+cr3+cr4)*ca;
                  cx_t ccp=(cr1*cr2+cr1*cr3+cr2*cr3+
                            cr1*cr4+cr2*cr4+cr3*cr4)*ca;
                  cx_t cdp=-(cr1*cr2*cr3+cr1*cr2*cr4+
                             cr1*cr3*cr4+cr2*cr3*cr4)*ca;
                  cx_t cep=cr1*cr2*cr3*cr4*ca;
                  
                  cx_t czo1=(((ca*cr1+cb)*cr1+cc)*cr1+cd)*cr1+ce;
                  cx_t czo2=(((ca*cr2+cb)*cr2+cc)*cr2+cd)*cr2+ce;
                  cx_t czo3=(((ca*cr3+cb)*cr3+cc)*cr3+cd)*cr3+ce;
                  cx_t czo4=(((ca*cr4+cb)*cr4+cc)*cr4+cd)*cr4+ce;
                  
                  fp_t q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
		      abs(cc-ccp)*abs(cc-ccp)+
		      abs(cd-cdp)*abs(cd-cdp)+
		      abs(ce-cep)*abs(ce-cep));
                  fp_t q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
                               abs(czo3)*abs(czo3)+abs(czo4)*abs(czo4));
                  
                  s1+=q1;
                  if (q1>m1) m1=q1;
                  s2+=q2;
                  if (q2>m2) m2=q2;
                  count++;
                  
                }
              }
            }
          }
        }
      }

      return count;
    }
    
    /// Return a string denoting the type ("quartic_real_coeff")
    const char *type() { return "quartic_real_coeff"; }
  };

  /** \brief Solve a quartic polynomial with complex coefficients and 
      complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quartic_complex : public quartic_real_coeff<fp_t,cx_t> {

  public:

    virtual ~quartic_complex() {}

    /** \brief Solves the polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_r(const fp_t a4, const fp_t b4, const fp_t c4, 
			const fp_t d4, const fp_t e4, fp_t &x1, 
			fp_t &x2, fp_t &x3, fp_t &x4) {
			      
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_complex::solve_r().",
           exc_einval);
      }
      
      cx_t r1,r2,r3,r4;
      int ret=solve_c(a4,b4,c4,d4,e4,r1,r2,r3,r4);
      x1=r1.real();
      x2=r2.real();
      x3=r3.real();
      x4=r4.real();
      return ret;
    }

    /**  \brief Solves the polynomial 
	 \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	 giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	 \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_rc(const fp_t a4, const fp_t b4, const fp_t c4, 
			 const fp_t d4, const fp_t e4, 
			 cx_t &x1, cx_t &x2, cx_t &x3, cx_t &x4) {
			 
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in cx_t &x4) {().",
           exc_einval);
      }
      
      return solve_c(a4,b4,c4,d4,e4,x1,x2,x3,x4);
    }      

    /** \brief Solves the complex polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_c(const cx_t a4, const cx_t b4, 
			const cx_t c4, const cx_t d4, 
			const cx_t e4, cx_t &x1, 
			cx_t &x2, cx_t &x3, cx_t &x4)=0;

    /** \brief Desc
     */
    int test_complex_coeffs(fp_t &s1, fp_t &s2, fp_t &m1,
                             fp_t &m2, size_t n=9) {
  
      cx_t i(0.0,1.0);
      cx_t one(1,0);

      size_t count=0;

      gen_test_number<> ga, gb, gc, gd, ge;
      fp_t rca, rcb, rcc, rcd, rce;
      for(int it=0;it<2;it++) {
        for(size_t j1=0;j1<n;j1++) {
          rca=ga.gen();
          gb.reset();
          for(size_t j2=0;j2<n;j2++) {
            rcb=gb.gen();
            gc.reset();
            for(size_t j3=0;j3<n;j3++) {
              rcc=gc.gen();
              gd.reset();
              for(size_t j4=0;j4<n;j4++) {
                rcd=gd.gen();
                ge.reset();
                for(size_t j5=0;j5<n;j5++) {
                  rce=ge.gen();

                  cx_t ca, cb, cc, cd, ce;
                  if (it==0) {
                    ca=rca+i;
                    cb=rcb+i;
                    cc=rcc+i;
                    cd=rcd+i;
                    ce=rce+i;
                  } else {
                    ca=one+i*rca;
                    cb=one+i*rcb;
                    cc=one+i*rcc;
                    cd=one+i*rcd;
                    ce=one+i*rce;
                  }

                  cx_t cr1, cr2, cr3, cr4;
                  solve_c(ca,cb,cc,cd,ce,cr1,cr2,cr3,cr4);
	    
                  cx_t cbp=-(cr1+cr2+cr3+cr4)*ca;
                  cx_t ccp=(cr1*cr2+cr1*cr3+cr2*cr3+
                            cr1*cr4+cr2*cr4+cr3*cr4)*ca;
                  cx_t cdp=-(cr1*cr2*cr3+cr1*cr2*cr4+
                             cr1*cr3*cr4+cr2*cr3*cr4)*ca;
                  cx_t cep=cr1*cr2*cr3*cr4*ca;
	    
                  cx_t czo1=(((ca*cr1+cb)*cr1+cc)*cr1+cd)*cr1+ce;
                  cx_t czo2=(((ca*cr2+cb)*cr2+cc)*cr2+cd)*cr2+ce;
                  cx_t czo3=(((ca*cr3+cb)*cr3+cc)*cr3+cd)*cr3+ce;
                  cx_t czo4=(((ca*cr4+cb)*cr4+cc)*cr4+cd)*cr4+ce;
                  
                  fp_t q1=sqrt(abs(cb-cbp)*abs(cb-cbp)+
                               abs(cc-ccp)*abs(cc-ccp)+
                               abs(cd-cdp)*abs(cd-cdp)+
                               abs(ce-cep)*abs(ce-cep));
                  fp_t q2=sqrt(abs(czo1)*abs(czo1)+abs(czo2)*abs(czo2)+
                               abs(czo3)*abs(czo3)+abs(czo4)*abs(czo4));
                  
                  //std::cout << q1 << " " << q2 << std::endl;
                  if (!o2isfinite(q1) || !o2isfinite(q2)
                      || q1>1.0e-10 || q2>1.0e-10) {
                    std::cout << "ca,cb,cc,cd,ce: "
                              << ca << " " << cb << " " << cc << " "
                              << cd << " " << ce << std::endl;
                    std::cout << "cr1,cr2,cr3,cr4: "
                              << cr1 << " " << cr2 << " "
                              << cr3 << " " << cr4 << std::endl;
                    exit(-1);
                    O2SCL_ERR("Failure in test_quartic_complex().",
                              exc_esanity);
                  }
                  s1+=q1;
                  if (q1>m1) m1=q1;
                  s2+=q2;
                  if (q2>m2) m2=q2;
                  
                  count++;
                }
                
              }
            }
          }
        }
      }
      
      return count;
    }
    
    /// Return a string denoting the type ("quartic_complex")
    const char *type() { return "quartic_complex"; }
  };

  /** \brief Solve a general polynomial with real
      coefficients and complex roots [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t>,
           class coeff_vec_t=std::vector<fp_t>,
           class root_vec_t=std::vector<cx_t> >
  class poly_real_coeff : public quadratic_real_coeff<fp_t,cx_t>,
                          public cubic_real_coeff<fp_t,cx_t>,
                          public quartic_real_coeff<fp_t,cx_t> {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    
    virtual ~poly_real_coeff() {}
    
    /** \brief Solve the n-th order polynomial
	
        The coefficients are stored in co[], with the leading coefficient
	as co[0] and the constant term as co[n]. The roots are returned
	in ro[0],...,ro[n-1].
    */
    virtual int solve_rc_arr(int n, const coeff_vec_t &co, 
                             root_vec_t &ro)=0;

    /** \brief Polish the roots
     */
    virtual int polish_rc_arr(int n, const coeff_vec_t &co, 
                              root_vec_t &ro) {

      /*
        AWS 4/10/21: This doesn't work right now, because the solver
        almost always fails. I think this could be fixed by an exact
        Jacobian, which doesn't yet work.
       */
      
      o2scl::mroot_hybrids<> mh;
      mm_funct mf=std::bind
        (std::mem_fn<int(size_t,const ubvector &,ubvector &,
                         const coeff_vec_t &,size_t)>
         (&poly_real_coeff<fp_t,cx_t,coeff_vec_t,root_vec_t>::polish_fun),
         this,std::placeholders::_1,std::placeholders::_2,
         std::placeholders::_3,std::cref(co),n);
      ubvector x(2), y(2);
      // Polish all each one of n roots
      for(int j=0;j<n;j++) {
        x[0]=ro[j].real();
        x[1]=ro[j].imag();
        mh.verbose=2;
        mh.err_nonconv=false;
        mh.def_jac.err_nonconv=false;
        mh.tol_rel=1.0e-12;
        mh.tol_abs=1.0e-12;
        mf(2,x,y);
        std::cout << y[0] << " x " << y[1] << std::endl;
        if (fabs(y[0])>1.0e-12 || fabs(y[1])>1.0e-12) {
          int mret=mh.msolve(2,x,mf);
          std::cout << "mret: " << mret << std::endl;
          if (mret==0) {
            ro[j].real(x[0]);
            ro[j].imag(x[1]);
          }
        }
      }
      return 0;
    }

    /** \brief Polish roots
     */
    virtual int polish_fun(size_t nv, const ubvector &x,
                           ubvector &y, const coeff_vec_t &co,
                           size_t n) {
      
      // The value x[0] is the real part and x[1] is the imaginary
      // part of the root. The value y[0] is the real part of the
      // evaluated polynomial expression and y[1] is the imaginary
      // part.

      //std::cout << "pf1: " << x[0] << " " << x[1] << std::endl;
      
      // Using horner's method following, e.g. GSL
      //std::cout << co[0] << std::endl;
      y[0]=co[0];
      y[1]=0.0;
      for(size_t i=0;i<n;i++) {
        //std::cout << co[i+1] << std::endl;
        fp_t tmp=co[i+1]+x[0]*y[0]-x[1]*y[1];
        y[1]=x[1]*y[0]+x[0]*y[1];
        y[0]=tmp;
      }
      
      std::cout << "pf2: " << y[0] << " " << y[1] << std::endl;
      
      return 0;
    }

    /** \brief Desc
     */
    int polish_jac(size_t nx, ubvector &x, size_t ny, ubvector &y,
                   ubmatrix &j) {
      
      // In the original function y[0] is the real part of the
      // polynomial, and y[1] is the imaginary part of the polynomial.
      // Thus j(0,1)? is the derivative of the real part of the
      // polynomial with respect to the imaginary part of the root.

      /*
      j(0,0)=co[0]*n;
      y(1,0)=0.0;
      j(0,1)=0.0;
      y(1,1)=co[0]*n;
      for(size_t i=0;i<n-1;i++) {
        fp_t tmp=co[i+1]*(n-1-i)+x[0]*y[0]-x[1]*y[1];
        j(1,0)=x[1]*y[0]+x[0]*y[1];
        j(0,0)=tmp;
        j(1,1)=tmp;
        j(0,1)=x[1]*y[0]+x[0]*y[1];
      }
      */
      
      return 0;
    }
  
    /// Return a string denoting the type ("poly_real_coeff")
    const char *type() { return "poly_real_coeff"; }
    
  };

  /** \brief Solve a general polynomial with complex
      coefficients [abstract base]
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t>,
           class coeff_vec_t=std::vector<fp_t>,
           class root_vec_t=std::vector<cx_t> >
  class poly_complex : public quadratic_complex<fp_t, cx_t>,
                       public cubic_complex<fp_t, cx_t>,
                       public quartic_complex<fp_t, cx_t> {

  public:

    virtual ~poly_complex() {}
    
    /** \brief Solve the n-th order polynomial
	
        The coefficients are stored in co[], with the leading coefficient
	as co[0] and the constant term as co[n]. The roots are returned
	in ro[0],...,ro[n-1].
    */
    virtual int solve_c_arr(int n, const coeff_vec_t &co,
                            root_vec_t &ro)=0;
    
    /// Polish the roots 
    virtual int polish_c_arr(int n, const coeff_vec_t &co,
                             root_vec_t &ro)=0;

    /// Return a string denoting the type ("poly_complex")
    const char *type() { return "poly_complex"; }
    
  };

  /** \brief Solve a cubic with real coefficients and complex roots 
      (CERNLIB)

      \note The function rrteq3() is based on the CERNLIB routine of
      the same name, but differs slightly. See the documentation of
      that function for details.
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_real_coeff_cern : public cubic_real_coeff<fp_t,cx_t> {

  public:

    cubic_real_coeff_cern() {
      eps=1.0e-6;
      delta=1.0e-15;
      improve_scale=true;
    }

    /// Numerical tolerance (default \f$ 10^{-6} \f$)
    fp_t eps;

    /// Numerical tolerance (default \f$ 10^{-15} \f$)
    fp_t delta;

    /// Improve algorithm for poorly-scaled roots (default true)
    bool improve_scale;

    virtual ~cubic_real_coeff_cern() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ giving the real 
	solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &r1, 
			 cx_t &r2,
                         cx_t &r3) {
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in ",
                   "cubic_real_coeff_cern::solve_rc().",exc_einval);
      }
      
      fp_t x[3],d;
      cx_t i(0.0,1.0);
      
      rrteq3(b3/a3,c3/a3,d3/a3,x,d);
      if (d>0.0) {
        r1=x[0];
        r2=x[1]+i*x[2];
        r3=x[1]-i*x[2];
        return 1;
      } else {
        r1=x[0];
        r2=x[1];
        r3=x[2];
      }
      
      return 3;
    }

    // Generic sign function
    inline fp_t sign(fp_t a, fp_t b) {
      if (b>=0.0) return fabs(a);
      return -fabs(a);
    }

    /** \brief The CERNLIB-like interface

	This function computes the roots of the cubic equation
	\f[
	x^3 + r x^2 + s x + t =0
	\f]
	returning the value of the discriminant in \c d and the roots
	in the array \c x. If the discriminant is negative, then all
	three real roots are stored in \c x. Otherwise, the real root
	is stored in <tt>x[0]</tt> and the real and imaginary parts of
	the complex conjugate roots are stored in <tt>x[1]</tt> and
	<tt>x[2]</tt>, respectively. This differs from the CERNLIB
	routine where the results were stored in <tt>x[1]</tt>,
	<tt>x[2]</tt>, and <tt>x[3]</tt> instead.
	
	Another small change is that the discriminant for the
	resolvent cubic is evaluated slightly differently in order to
	improve the properties in the case where the roots are not all
	of order unity. The default CERNLIB behavior can be restored
	by setting improve_scale to \c false.
    */	
    virtual int rrteq3(fp_t r, fp_t s, fp_t t, fp_t x[],
                       fp_t &d) {
      
      fp_t delta2=delta;
      fp_t r1=2.0;
      r1/=27.0;
      fp_t r2=1.0;
      r2/=2.0;
      fp_t r3=1.0;
      r3/=3.0;
      fp_t w3=3.0;
      w3=sqrt(w3);
      fp_t r4=w3/2.0;
      fp_t q1=2.0;
      q1/=27.0;
      fp_t q2=r2;
      fp_t q3=r3;
      fp_t y[3];
      cx_t z[3], i(0.0,1.0);
      fp_t h2, h3;
      int j,k;

      // AWS 4/1/21: These temporaries fix compiling for long double
      // and complex<long double> types
      fp_t one=1.0;
      fp_t two=2.0;
      fp_t three=3.0;
      
      if (s==0.0 && t==0.0) {
        x[0]=-r;
        x[1]=0.0;
        x[2]=0.0;
        d=0;
        return success;
      }
      fp_t p=s-r3*r*r;
      fp_t q=(r1*r*r-r3*s)*r+t;
      d=r2*r2*q*q+r3*p*r3*p*r3*p;
      if (o2abs(d)<=eps) {
        fp_t pp=s-q3*r*r;
        fp_t qq=(q1*r*r-q3*s)*r+t;
        d=q2*q2*qq*qq+q3*pp*q3*pp*q3*pp;
        p=pp;
        q=qq;
      }
      fp_t h=r3*r;
      fp_t h1=r2*q;
      fp_t u,v,d_new;
      
      // The discriminant in 'd' has units of [x]^6 so it is very
      // sensitive to the absolute magnitude of the roots. We attempt to
      // fix this by using the ratio instead of the sum.
      if (improve_scale) {
        fp_t da=r2*r2*q*q;
        fp_t db=r3*p*r3*p*r3*p;
        if (db==0.0) {
          delta2=0.0;
          d_new=da;
        } else if (db>0.0) {
          d_new=da/db+one;
        } else {
          d_new=-da/db-one;
        }
      } else {
        d_new=d;
      }
      
      if (d_new>delta2) {
        h2=o2sqrt(d);
        fp_t u0=-h1+h2;
        fp_t v0=-h1-h2;
        if (o2abs(u0)==0.0) {
          u=sign(0.0,u0);
        }
        else u=sign(o2pow(o2abs(u0),r3),u0);
        if (o2abs(v0)==0.0) v=sign(0.0,v0);
        else v=sign(o2pow(o2abs(v0),r3),v0);
        x[0]=u+v-h;
        x[1]=-r2*(u+v)-h;
        fp_t tmp=u-v;
        x[2]=r4*o2abs(tmp);
        if (o2abs(u0)<=eps || o2abs(v0)<=eps) {
          y[0]=x[0];
          for(k=0;k<=1;k++) {
            y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/
              ((three*y[k]+two*r)*y[k]+s);
          }
          x[0]=y[2];
          z[0]=x[1]+i*x[2];
          for(k=0;k<=1;k++) {
            cx_t den=three*z[k];
            den-=two*r;
            z[k+1]=z[k]-(((z[k]+r)*z[k]+s)*z[k]+t)/
              (den*z[k]+s);
          }
          x[1]=z[2].real();
          x[2]=z[2].imag();
        }
        
      } else if (o2abs(d_new)<=delta2) {
        
        d=0.0;
        if (o2abs(h1)==0.0) u=sign(0.0,-h1);
        else u=sign(o2pow(o2abs(h1),r3),-h1);
        x[0]=u+u-h;
        x[1]=-u-h;
        x[2]=x[1];
        if (o2abs(h1)<=eps) {
          y[0]=x[0];
          for(k=0;k<=1;k++) {
            h1=(three*y[k]+two*r)*y[k]+s;
            if (o2abs(h1)>delta2) {
              y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/h1;
            } else {
              x[0]=-r3*r;
              x[1]=x[0];
              x[2]=x[0];
              return success;
            }
          }
          x[0]=y[2];
          x[1]=-r2*(r+x[0]);
          x[2]=x[1];
        }
        
      } else {

        fp_t tmp2=r3*p;
        h3=o2abs(tmp2);
        fp_t tmp3=h3*h3*h3;
        h3=o2sqrt(tmp3);
        h2=r3*acos(-h1/h3);
        if (h3==0.0) h1=0.0;
        else h1=o2pow(h3,r3);
        u=h1*cos(h2);
        v=w3*h1*sin(h2);
        x[0]=u+u-h;
        x[1]=-u-v-h;
        x[2]=-u+v-h;
        if (h3<=eps || x[0]<=eps || x[1]<=eps || x[2]<=eps) {
          for(j=0;j<3;j++) {
            y[0]=x[j];
            for(k=0;k<=1;k++) {
              y[k+1]=y[k]-(((y[k]+r)*y[k]+s)*y[k]+t)/
                ((three*y[k]+two*r)*y[k]+s);
            }
            x[j]=y[2];
          }
        }
      }
      
      return success;
    }
    
    /// Return a string denoting the type ("cubic_real_coeff_cern")
    const char *type() { return "cubic_real_coeff_cern"; }
  };

  /** \brief Solve a quartic with real coefficients and complex 
      roots (CERNLIB)
  */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quartic_real_coeff_cern : public quartic_real_coeff<fp_t,cx_t> {
    
  public:

    virtual ~quartic_real_coeff_cern() {}

    /** \brief Solves the polynomial \f$ a_4 x^4 + b_4 x^3 + c_4 x^2 +
	d_4 x + e_4= 0 \f$ giving the four complex solutions \f$ x=x_1
	\f$ , \f$ x=x_2 \f$ , \f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_rc(const fp_t a4, const fp_t b4, const fp_t c4, 
			 const fp_t d4, const fp_t e4, 
			 cx_t &x1, cx_t &x2, cx_t &x3, cx_t &x4) {
      
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_real_coeff_cern::solve_rc().",
           exc_einval);
      }
      
      int mt;
      fp_t dc;
      cx_t x[4];
      
      rrteq4(b4/a4,c4/a4,d4/a4,e4/a4,x,dc,mt);
      x1=x[0];
      x2=x[1];
      x3=x[2];
      x4=x[3];
      
      return success;
    }      

    /** \brief The CERNLIB-like interface

        There are a couple differences with the original routine.
        The arrays z[] and u[] are now zero-indexed.

        \future Make v[] zero-indexed as well.
    */
    virtual int rrteq4(fp_t a, fp_t b, fp_t c, fp_t d, 
		       cx_t z[], fp_t &dc, int &mt) {
      
      cx_t i(0.0,1.0), z0[5];
      cx_t w1(0.0,0.0), w2(0.0,0.0), w3;

      fp_t q2=1;
      q2/=2;
      fp_t r4=q2/2;
      fp_t q4=r4;
      fp_t q8=r4/2;
      fp_t r12=r4/3;
      fp_t q1=q8*3;
      fp_t q3=q1/2;
      fp_t eight=8;
        
      fp_t u[3], v[4], v1, v2;
      int j, k1=0, k2=0;
      
      // degenerate cases
      if (b==0 && c==0) {
        if (d==0) {
          mt=1;
          z[0]=-a;
          z[1]=0;
          z[2]=0;
          z[3]=0;
          dc=0;
          return success;
        } else if (a==0) {
          if (d>0) {
            mt=2;
            z[0]=sqrt(i*sqrt(d));
            z[1]=-z[0];
            z[3]=sqrt(-z[0]*z[0]);
            z[2]=-z[3];
          } else {
            mt=3;
            z[0]=sqrt(sqrt(-d));
            z[1]=-z[0];
            z[2]=sqrt(-z[0]*z[0]);
            z[3]=-z[2];
          }
          dc=-r12*d*r12*d*r12*d;
          return success;
        }
      }
      
      // Solve the resolvent cubic
      fp_t aa=a*a;
      fp_t pp=b-q1*aa;
      fp_t qq=c-q2*a*(b-q4*aa);
      fp_t rr=d-q4*(a*c-q4*aa*(b-q3*aa));
      fp_t rc=q2*pp;
      fp_t sc=q4*(q4*pp*pp-rr);
      fp_t tc=-(q8*qq*q8*qq);
      fp_t zero=0.0;

      cub_obj.rrteq3(rc,sc,tc,u,dc);

      fp_t q=qq;
      fp_t h=r4*a;
      if (dc==0) u[2]=u[1];
      if (dc<=0) {
        mt=2;
        v[1]=o2abs(u[0]);
        v[2]=o2abs(u[1]);
        v[3]=o2abs(u[2]);
        v1=std::max(std::max(v[1],v[2]),v[3]);
        if (v1==v[1]) {
          k1=0;
          v2=std::max(v[2],v[3]);
        } else if (v1==v[2]) {
          k1=1;
          v2=std::max(v[1],v[3]);
        } else {
          k1=2;
          v2=std::max(v[1],v[2]);
        }
        if (v2==v[1]) {
          k2=0;
        } else if (v2==v[2]) {
          k2=1;
        } else {
          k2=2;
        }
        w1=sqrt(((cx_t)(u[k1])));
        w2=sqrt(((cx_t)(u[k2])));
      } else {
        mt=3;
        w1=sqrt(u[1]+i*u[2]);
        w2=sqrt(u[1]-i*u[2]);
      }
      w3=0;
      if (w1*w2!=zero) w3=-q/(eight*w1*w2);
      z0[1]=w1+w2+w3-h;
      z0[2]=-w1-w2+w3-h;
      z0[3]=-w1+w2-w3-h;
      z0[4]=w1-w2-w3-h;
      if (mt==2) {
        if (u[k1]>=0 && u[k2]>=0) {
          mt=1;
          for(j=1;j<=4;j++) {
            z[j-1]=z0[j].real();
          }
        } else if (u[k1]>=0 && u[k2]<0) {
          z[0]=z0[1];
          z[1]=z0[4];
          z[2]=z0[3];
          z[3]=z0[2];
        } else if (u[k1]<0 && u[k2]>=0) {
          z[0]=z0[1];
          z[1]=z0[3];
          z[2]=z0[4];
          z[3]=z0[2];
        } else if (u[k1]<0 && u[k2]<0) {
          z[0]=z0[1];
          z[1]=z0[2];
          z[2]=z0[4];
          z[3]=z0[3];
        }
      } else if (mt==3) {
        for(j=1;j<=2;j++) {
          z[j-1]=z0[j].real();
        }
        z[2]=z0[4];
        z[3]=z0[3];
      }

      return success;
    }
    
    /// Return a string denoting the type ("quartic_real_coeff_cern")
    const char *type() { return "quartic_real_coeff_cern"; }

#ifndef DOXYGEN_INTERNAL

  public:

    /// The object to solve for the associated cubic
    cubic_real_coeff_cern<fp_t,cx_t> cub_obj;

#endif

  };

  /** \brief Solve a quadratic with real coefficients and complex roots (GSL)
   */
  class quadratic_real_coeff_gsl : public quadratic_real_coeff<> {
    
  public:
    
    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_rc(const double a, const double b, const double c,
                         std::complex<double> &x1, std::complex<double> &x2);
    
  };

  /** \brief Solve a quadratic with real coefficients and complex roots 
      (C++ rewrite of GSL algorithm)
      
      This code was based on GSL v2.6's function
      gsl_poly_complex_solve_quadratic().
   */
  template<class fp_t=double, class cx_t=std::complex<double> >
  class quadratic_real_coeff_gsl2 : public quadratic_real_coeff<fp_t,cx_t> {

  public:

    virtual ~quadratic_real_coeff_gsl2() {}

    /** \brief Solves the polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 

        This function returns the number of real roots (either 0
        or 2)
    */
    virtual int solve_rc(const fp_t a, const fp_t b, const fp_t c,
                         cx_t &x1, cx_t &x2) {

      // AWS, 3/31/21: This is nearly equivalent to the GSL code from v2.6,
      // but rewritten in a C++/multiprecision compatible form.
      // It also does not allow a zero value for a (the GSL version
      // does).
      
      if (a==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_real_coeff_gsl2::solve_rc().",
           exc_einval);
      }
      
      fp_t disc=b*b-4*a*c;

      if (disc > 0) {
        if (b == 0) {
          fp_t s=fabs(sqrt(disc)/a/2);
          x1.real(-s);
          x1.imag(0);
          x2.real(s);
          x2.imag(0);
        } else {
          fp_t sgnb=(b > 0 ? 1 : -1);
          fp_t temp=-(b+sgnb*sqrt(disc))/2;
          fp_t r1=temp/a;
          fp_t r2=c/temp;
          
          x1.real(r2);
          x1.imag(0);
          x2.real(r1);
          x2.imag(0);
        }
        return 2;
        
      } else if (disc == 0) {
        
        x1.real(-b/a/2);
        x1.imag(0);
        x2.real(-b/a/2);
        x2.imag(0);
        
        return 2;
        
      }
      
      fp_t s=fabs(sqrt(-disc)/a/2);
      x1.real(-b/a/2);
      x1.imag(-s);
      x2.real(-b/a/2);
      x2.imag(s);
      
      return 0;
    }

    /// Return a string denoting the type ("quadratic_real_coeff_gsl2")
    const char *type() { return "quadratic_real_coeff_gsl2"; }

  };

  /** \brief Solve a cubic with real coefficients and complex roots (GSL)

      This is just a wrapper for gsl_poly_complex_solve_cubic().
   */
  class cubic_real_coeff_gsl : public cubic_real_coeff<double> {

  public:

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the real solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const double a3, const double b3, const double c3, 
			 const double d3, double &x1, 
			 std::complex<double> &x2, std::complex<double> &x3);
    
  };

  /** \brief Solve a cubic with real coefficients and complex roots
      (C++ rewrite of GSL algorithm)

      This class is designed to be identical to the GSL code in
      gsl_poly_complex_solve_cubic(), but generalized to generic
      floating point types.
   */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_real_coeff_gsl2 : public cubic_real_coeff<fp_t,cx_t> {

  protected:
    
    inline void swap(double &a, double &b) {
      double tmp=b;
      b=a;
      a=tmp;
      return;
    }
    
  public:

    virtual ~cubic_real_coeff_gsl2() {}

    /** \brief Solves the polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the real solution \f$ x=x_1 \f$ and two complex solutions 
	\f$ x=x_2 \f$ and \f$ x=x_3 \f$ .
    */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &x1, cx_t &x2, cx_t &x3) {
      
      if (a3==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in cubic_real_coeff_gsl2::solve_rc().",
           exc_einval);
      }
      
      fp_t a=b3/a3;
      fp_t b=c3/a3;
      fp_t c=d3/a3;
      fp_t one=1;
      fp_t third=one/3;
      fp_t half=one/2;
      fp_t three=3;
      fp_t root3=sqrt(three);
      
      fp_t q=(a*a-3*b);
      fp_t r=(2*a*a*a-9*a*b+27*c);
      
      fp_t Q=q/9;
      fp_t R=r/54;
      
      fp_t Q3=Q*Q*Q;
      fp_t R2=R*R;
      
      fp_t CR2=729*r*r;
      fp_t CQ3=2916*q*q*q;

      fp_t pi=boost::math::constants::pi<fp_t>();
      
      if (R == 0 && Q == 0) {
        
        x1=-a/3;
        x2.real(-a/3);
        x2.imag(0);
        x3.real(-a/3);
        x3.imag(0);

        if (!o2isfinite(x1) || !o2isfinite(x2.real()) ||
            !o2isfinite(x2.imag())) {
          std::cout << "1. " << x1 << " " << x2 << " " << x3 << std::endl;
          exit(-1);
        }
        return 3;
        
      } else if (CR2 == CQ3)  {

        /* [GSL] This test is actually r2 == q3, written in a form suitable
           for exact computation with integers */
        
        /* [GSL] Due to finite precision some double roots may be
           missed, and will be considered to be a pair of complex
           roots z=x +/- epsilon i close to the real axis. */
        
        fp_t sqrtQ=sqrt(Q);
        
        if (R > 0) {
          x1=-2*sqrtQ-a/3;
          x2.real(sqrtQ-a/3);
          x2.imag(0);
          x3.real(sqrtQ-a/3);
          x3.imag(0);
        } else {
          x1=-sqrtQ-a/3;
          x2.real(-sqrtQ-a/3);
          x2.imag(0);
          x3.real(2*sqrtQ-a/3);
          x3.imag(0);
        }
        
        if (!o2isfinite(x1) || !o2isfinite(x2.real()) ||
            !o2isfinite(x2.imag())) {
          std::cout << "2. " << x1 << " " << x2 << " " << x3 << std::endl;
          exit(-1);
        }
        return 3;
        
      } else if (CR2 < CQ3)  {

        /* [GSL] equivalent to R2 < Q3 */
        
        fp_t sqrtQ=sqrt(Q);
        fp_t sqrtQ3=sqrtQ*sqrtQ*sqrtQ;
        fp_t theta=acos(R/sqrtQ3);
        
        // [AWS] Modified from the original GSL routine
        // sometimes R/sqrtQ3 is slightly larger than one
        // when the coefficients are arranged just right
        // so theta becomes not finite.
        
        if (R/sqrtQ3>=1.0) theta=0.0;
        if (R/sqrtQ3<=-1.0) theta=pi;
        
        fp_t norm=-2*sqrtQ;
        
        fp_t r0=norm*cos(theta/3)-a/3;
        fp_t r1=norm*cos((theta+2*pi)/3)-a/3;
        fp_t r2=norm*cos((theta-2*pi)/3)-a/3;
        
        x1=r0;
        x2.real(r1);
        x2.imag(0);
        x3.real(r2);
        x3.imag(0);
        
        if (!o2isfinite(x1) || !o2isfinite(x2.real()) ||
            !o2isfinite(x2.imag())) {
          std::cout << "3. " << x1 << " " << x2 << " " << x3 << std::endl;
          std::cout << Q << std::endl;
          std::cout << norm << " " << sqrtQ << " " << theta << std::endl;
          std::cout << r0 << " " << r1 << " " << r2 << std::endl;
          exit(-1);
        }
        return 3;
        
      }
        
      fp_t sgnR=(R>=0?1:-1);
      fp_t A;
      
      // [AWS] Modification from original GSL behavior: Just in case
      // R2=Q3, finite precision can cause the argument of the sqrt()
      // function to be negative. We correct for this here.
      if (R2<=Q3) {
        A=-sgnR*pow(fabs(R),third);
      } else {
        A=-sgnR*pow(fabs(R)+sqrt(R2-Q3),third);
      }
      
      fp_t B=Q/A;
      
      x1=A+B-a/3;
      x2.real(-half*(A+B)-a/3);
      x2.imag(-(root3*half)*fabs(A-B));
      x3.real(-half*(A+B)-a/3);
      x3.imag((root3*half)*fabs(A-B));
      
      if (!o2isfinite(x1) || !o2isfinite(x2.real()) ||
          !o2isfinite(x2.imag())) {
        std::cout << R2-Q3 << std::endl;
        std::cout << fabs(R)+sqrt(R2-Q3) << std::endl;
        std::cout << "4. " << x1 << " " << x2 << " " << x3 << std::endl;
        std::cout << sgnR << " " << R2 << " " << Q3 << " " << R << " "
                  << B << " " << Q << " " << A << " " << a << std::endl;
        exit(-1);
      }
        
      return 1;
    }

    /// Return a string denoting the type ("cubic_real_coeff_gsl2")
    const char *type() { return "cubic_real_coeff_gsl2"; }

  };

  /** \brief Use multiprecision to automatically solve a cubic to 
      a specified level of precision

      This class will fail to evalate a function with the requested
      precision if:
      - The user-specified input and result data type does not have enough
      precision to compute or store the result 
      - The requested precision is near to or smaller than 1.0e-50
      - The function is noisy, non-deterministic, or is not 
      continuous in the local neighborhood

      \note The algorithm attempts not to be wasteful, but is not
      necessarily optimized for speed. 
  */
  class cubic_real_coeff_multip {
    
  protected:
    
    /// \name Typedefs for multiprecision types
    //@{
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<25> >
    cpp_dec_float_25;
  
    typedef
    boost::multiprecision::number<boost::multiprecision::cpp_dec_float<35> >
    cpp_dec_float_35;
  
    typedef boost::multiprecision::cpp_dec_float_50 cpp_dec_float_50;
  
    typedef boost::multiprecision::cpp_dec_float_100 cpp_dec_float_100;
    //@}

    /// \name Solvers
    //@{
    cubic_real_coeff_cern<double,std::complex<double>> q_d;
    cubic_real_coeff_cern<long double,std::complex<long double>> q_ld;
    cubic_real_coeff_cern<cpp_dec_float_25,std::complex<cpp_dec_float_25>>
    q_cdf25;
    cubic_real_coeff_cern<cpp_dec_float_35,std::complex<cpp_dec_float_35>>
    q_cdf35;
    cubic_real_coeff_cern<cpp_dec_float_50,std::complex<cpp_dec_float_50>>
    q_cdf50;
    cubic_real_coeff_cern<cpp_dec_float_100,std::complex<cpp_dec_float_100>>
    q_cdf100;
    //@}
    
  public:

    cubic_real_coeff_multip() {
      verbose=0;
      tol_rel=-1.0;
      err_nonconv=true;
    }

    /** \brief Verbosity parameter
     */
    int verbose;

    /** \brief Relative tolerance
     */
    double tol_rel;

    /** \brief If true, call the error handler if the function
        evaluation fails
     */
    bool err_nonconv;

    /** \brief Evaluate the function and return the error estimate
        with the specified tolerance
     */
    template<class fp_t, class cx_t>
    int solve_rc_tol_err(const fp_t &a, const fp_t &b,
                                 const fp_t &c, const fp_t &d,
                                 fp_t &x1, cx_t &x2, cx_t &x3,
                                 fp_t &err, double tol_loc=-1) {
        
      /// Tolerance choice and verification logic
      
      if (tol_loc<=0.0 && tol_rel<=0.0) {
        tol_loc=pow(10.0,-std::numeric_limits<fp_t>::digits10);
        if (verbose>0) {
          std::cout << "Set tolerance from data type to: "
                    << tol_loc << std::endl;
        }
      } else if (tol_loc<=0.0) {
        if (tol_rel<pow(10.0,-std::numeric_limits<fp_t>::digits10)) {
          std::cerr << "Class data member tol_rel is " << tol_rel
                    << " but data type only stores "
                    << std::numeric_limits<fp_t>::digits10
                    << " digits." << std::endl;
          O2SCL_ERR("Cannot compute to required precision",
                    o2scl::exc_einval);
        }
        tol_loc=tol_rel;
        if (verbose>0) {
          std::cout << "Set tolerance from value of tol_rel to: "
                    << tol_loc << std::endl;
        }
      } else {
        if (tol_loc<pow(10.0,-std::numeric_limits<fp_t>::digits10)) {
          std::cerr << "Caller requested tolerance " << tol_loc
                    << " but data type only stores "
                    << std::numeric_limits<fp_t>::digits10
                    << " digits." << std::endl;
          O2SCL_ERR("Cannot compute to required precision",
                    o2scl::exc_einval);
        } else if (verbose>0) {
          std::cout << "Set tolerance from user-specified value to: "
                    << tol_loc << std::endl;
        }
      }

      if (verbose>0) {
        std::cout << "Cubic: "
                  << dtos(a,0) << " "
                  << dtos(b,0) << " "
                  << dtos(c,0) << " "
                  << dtos(d,0) << std::endl;
      }
          
      /// Degenerate cubics
      if (b==0 && c==0 && d==0) {
        x1=0;
        x2=0;
        x3=0;
        return 3;
      }

      /// First pass, compare double and long double
      
      double a_d=static_cast<double>(a);
      double b_d=static_cast<double>(b);
      double c_d=static_cast<double>(c);
      double d_d=static_cast<double>(d);
      double x1_d;
      std::complex<double> x2_d, x3_d;
      long double a_ld=static_cast<long double>(a);
      long double b_ld=static_cast<long double>(b);
      long double c_ld=static_cast<long double>(c);
      long double d_ld=static_cast<long double>(d);
      long double x1_ld;
      std::complex<long double> x2_ld, x3_ld;
      
      int ret_d=q_d.solve_rc(a_d,b_d,c_d,d_d,x1_d,x2_d,x3_d);
      int ret_ld=q_ld.solve_rc(a_ld,b_ld,c_ld,d_ld,x1_ld,x2_ld,x3_ld);
      
      if (ret_d==ret_ld) {
        err=0;
        if (d==0 && (ret_ld==1 ||
                     (abs(x1_ld)<abs(x2_ld) && abs(x1_ld)<abs(x3_ld)))) {
          x1_ld=0;
        } else {
          err=static_cast<fp_t>(abs(x1_ld-x1_d)/abs(x1_ld));
        }
        if (d==0 && abs(x2_ld)<abs(x1_ld) && abs(x2_ld)<abs(x3_ld)) {
          x2_ld=0;
        } else if (x2_ld.real()!=0 || x2_ld.imag()!=0 ||
                   x2_d.real()!=0 || x2_d.imag()!=0) {
          err+=static_cast<fp_t>(abs(static_cast<std::complex<double>>(x2_ld)-
                                     x2_d)/abs(x2_ld));
        }
        if (d==0 && abs(x3_ld)<abs(x1_ld) && abs(x3_ld)<abs(x2_ld)) {
          x3_ld=0;
        } else if (x3_ld.real()!=0 || x3_ld.imag()!=0 ||
                   x3_d.real()!=0 || x3_d.imag()!=0) {
          err+=static_cast<fp_t>
            (abs(static_cast<std::complex<double>>(x3_ld)-
                 x3_d)/abs(x3_ld));
        }
        if (err<tol_loc) {
          x1=static_cast<fp_t>(x1_ld);
          x2=static_cast<cx_t>(x2_ld);
          x3=static_cast<cx_t>(x3_ld);
          return ret_d;
        }
      }
      
      if (verbose>0) {
        std::cout << "Failed 1: " << ret_d << " " << ret_ld 
                  << "\n  " << dtos(x1_d,0) << " " << dtos(x1_ld,0) << " "
                  << "\n  (" << dtos(x2_d.real(),0) << ","
                  << dtos(x2_d.imag(),0) << ") ("
                  << dtos(x2_ld.real(),0) << ","
                  << dtos(x2_ld.imag(),0) << ")"
                  << "\n  (" << dtos(x3_d.real(),0) << ","
                  << dtos(x3_d.imag(),0) << ") ("
                  << dtos(x3_ld.real(),0) << ","
                  << dtos(x3_ld.imag(),0) << ")"
                  << "\n  " << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      /// Second pass, compare long double and 25-digit precision

      cpp_dec_float_25 a_cdf25=static_cast<cpp_dec_float_25>(a);
      cpp_dec_float_25 b_cdf25=static_cast<cpp_dec_float_25>(b);
      cpp_dec_float_25 c_cdf25=static_cast<cpp_dec_float_25>(c);
      cpp_dec_float_25 d_cdf25=static_cast<cpp_dec_float_25>(d);
      cpp_dec_float_25 x1_cdf25;
      std::complex<cpp_dec_float_25> x2_cdf25, x3_cdf25;
      
      int ret_cdf25=q_cdf25.solve_rc(a_cdf25,b_cdf25,c_cdf25,d_cdf25,
                                     x1_cdf25,x2_cdf25,x3_cdf25);

      if (ret_ld==ret_cdf25) {
        // If d is 0 and the first root is near 0, then avoid dividing
        // by it
        err=0;
        if (d==0 && (ret_cdf25==1 ||
                     (abs(x1_cdf25)<abs(x2_cdf25) &&
                      abs(x1_cdf25)<abs(x3_cdf25)))) {
          x1_cdf25=0;
        } else {
          err=static_cast<fp_t>(abs(x1_cdf25-x1_ld)/abs(x1_cdf25));
        }
        if (d==0 && abs(x2_cdf25)<abs(x1_cdf25) &&
            abs(x2_cdf25)<abs(x3_cdf25)) {
          x2_cdf25=0;
        } else if (x2_cdf25.real()!=0 || x2_cdf25.imag()!=0 ||
                   x2_ld.real()!=0 || x2_ld.imag()!=0) {
          err+=static_cast<fp_t>(abs(abs(x2_cdf25)-abs(x2_ld))/
                                 abs(x2_cdf25));
        }
        if (d==0 && abs(x3_cdf25)<abs(x1_cdf25) &&
            abs(x3_cdf25)<abs(x2_cdf25)) {
          x3_cdf25=0;
        } else if (x3_cdf25.real()!=0 || x3_cdf25.imag()!=0 ||
                   x3_ld.real()!=0 || x3_ld.imag()!=0) {
          err+=static_cast<fp_t>(abs(abs(x3_cdf25)-abs(x3_ld))/
                                 abs(x3_cdf25));
                                     
        }
        if (err<tol_loc) {
          x1=static_cast<fp_t>(x1_cdf25);
          x2.real(static_cast<fp_t>(x2_cdf25.real()));
          x2.imag(static_cast<fp_t>(x2_cdf25.imag()));
          x3.real(static_cast<fp_t>(x3_cdf25.real()));
          x3.imag(static_cast<fp_t>(x3_cdf25.imag()));
          return ret_ld;
        }
      }
      
      if (verbose>0) {
        std::cout << "Failed 2: " << ret_ld << " " << ret_cdf25 
                  << "\n  " << dtos(x1_ld,0) << " "
                  << dtos(x1_cdf25,0) << " "
                  << "\n  (" << dtos(x2_ld.real(),0) << ","
                  << dtos(x2_ld.imag(),0) << ") ("
                  << dtos(x2_cdf25.real(),0) << ","
                  << dtos(x2_cdf25.imag(),0) << ")"
                  << "\n  (" << dtos(x3_ld.real(),0) << ","
                  << dtos(x3_ld.imag(),0) << ") ("
                  << dtos(x3_cdf25.real(),0) << ","
                  << dtos(x3_cdf25.imag(),0) << ")"
                  << "\n  " << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      /// Third pass, compare 25- and 35-digit precision

      cpp_dec_float_35 a_cdf35=static_cast<cpp_dec_float_35>(a);
      cpp_dec_float_35 b_cdf35=static_cast<cpp_dec_float_35>(b);
      cpp_dec_float_35 c_cdf35=static_cast<cpp_dec_float_35>(c);
      cpp_dec_float_35 d_cdf35=static_cast<cpp_dec_float_35>(d);
      cpp_dec_float_35 x1_cdf35;
      std::complex<cpp_dec_float_35> x2_cdf35, x3_cdf35;
      
      int ret_cdf35=q_cdf35.solve_rc(a_cdf35,b_cdf35,c_cdf35,d_cdf35,
                                     x1_cdf35,x2_cdf35,x3_cdf35);

      if (ret_cdf25==ret_cdf35) {
        // If d is 0 and the first root is near 0, then avoid dividing
        // by it
        err=0;
        if (d==0 && (ret_cdf35==1 ||
                     (abs(x1_cdf35)<abs(x2_cdf35) &&
                      abs(x1_cdf35)<abs(x3_cdf35)))) {
          x1_cdf35=0;
        } else {
          err=static_cast<fp_t>(abs(x1_cdf35-x1_cdf25)/abs(x1_cdf35));
        }
        if (d==0 && abs(x2_cdf35)<abs(x1_cdf35) &&
            abs(x2_cdf35)<abs(x3_cdf35)) {
          x2_cdf35=0;
        } else if (x2_cdf35.real()!=0 || x2_cdf35.imag()!=0 ||
                   x2_cdf25.real()!=0 || x2_cdf25.imag()!=0) {
          err+=static_cast<fp_t>(abs(abs(x2_cdf35)-abs(x2_cdf25))/
                                 abs(x2_cdf35));
        }
        if (d==0 && abs(x3_cdf35)<abs(x1_cdf35) &&
            abs(x3_cdf35)<abs(x2_cdf35)) {
          x3_cdf35=0;
        } else if (x3_cdf35.real()!=0 || x3_cdf35.imag()!=0 ||
                   x3_cdf25.real()!=0 || x3_cdf25.imag()!=0) {
          err+=static_cast<fp_t>(abs(abs(x3_cdf35)-abs(x3_cdf25))/
                                 abs(x3_cdf35));
                                     
        }
        if (err<tol_loc) {
          x1=static_cast<fp_t>(x1_cdf35);
          x2.real(static_cast<fp_t>(x2_cdf35.real()));
          x2.imag(static_cast<fp_t>(x2_cdf35.imag()));
          x3.real(static_cast<fp_t>(x3_cdf35.real()));
          x3.imag(static_cast<fp_t>(x3_cdf35.imag()));
          return ret_cdf25;
        }
      }
      
      if (verbose>0) {
        std::cout << "Failed 3: " << ret_cdf25 << " " << ret_cdf35 
                  << "\n  " << dtos(x1_cdf25,0) << " "
                  << dtos(x1_cdf35,0) << " "
                  << "\n  (" << dtos(x2_cdf25.real(),0) << ","
                  << dtos(x2_cdf25.imag(),0) << ") ("
                  << dtos(x2_cdf35.real(),0) << ","
                  << dtos(x2_cdf35.imag(),0) << ")"
                  << "\n  (" << dtos(x3_cdf25.real(),0) << ","
                  << dtos(x3_cdf25.imag(),0) << ") ("
                  << dtos(x3_cdf35.real(),0) << ","
                  << dtos(x3_cdf35.imag(),0) << ")"
                  << "\n  " << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      /// Fourth pass, compare 35- and 50-digit precision
      
      cpp_dec_float_50 a_cdf50=static_cast<cpp_dec_float_50>(a);
      cpp_dec_float_50 b_cdf50=static_cast<cpp_dec_float_50>(b);
      cpp_dec_float_50 c_cdf50=static_cast<cpp_dec_float_50>(c);
      cpp_dec_float_50 d_cdf50=static_cast<cpp_dec_float_50>(d);
      cpp_dec_float_50 x1_cdf50;
      std::complex<cpp_dec_float_50> x2_cdf50, x3_cdf50;
      
      int ret_cdf50=q_cdf50.solve_rc(a_cdf50,b_cdf50,c_cdf50,d_cdf50,
                                     x1_cdf50,x2_cdf50,x3_cdf50);

      if (ret_cdf35==ret_cdf50) {
        // If d is 0 and the first root is near 0, then avoid dividing
        // by it
        err=0;
        if (d==0 && (ret_cdf50==1 ||
                     (abs(x1_cdf50)<abs(x2_cdf50) &&
                      abs(x1_cdf50)<abs(x3_cdf50)))) {
          x1_cdf50=0;
        } else {
          err=static_cast<fp_t>(abs(x1_cdf50-x1_cdf35)/abs(x1_cdf50));
        }
        if (d==0 && abs(x2_cdf50)<abs(x1_cdf50) &&
            abs(x2_cdf50)<abs(x3_cdf50)) {
          x2_cdf50=0;
        } else if (x2_cdf50.real()!=0 || x2_cdf50.imag()!=0 ||
                   x2_cdf35.real()!=0 || x2_cdf35.imag()!=0) {
          err+=static_cast<fp_t>(abs(abs(x2_cdf50)-abs(x2_cdf35))/
                                 abs(x2_cdf50));
        }
        if (d==0 && abs(x3_cdf50)<abs(x1_cdf50) &&
            abs(x3_cdf50)<abs(x2_cdf50)) {
          x3_cdf50=0;
        } else if (x3_cdf50.real()!=0 || x3_cdf50.imag()!=0 ||
                   x3_cdf35.real()!=0 || x3_cdf35.imag()!=0) {
          err+=static_cast<fp_t>(abs(abs(x3_cdf50)-abs(x3_cdf35))/
                                 abs(x3_cdf50));
                                     
        }
        if (err<tol_loc) {
          x1=static_cast<fp_t>(x1_cdf50);
          x2.real(static_cast<fp_t>(x2_cdf50.real()));
          x2.imag(static_cast<fp_t>(x2_cdf50.imag()));
          x3.real(static_cast<fp_t>(x3_cdf50.real()));
          x3.imag(static_cast<fp_t>(x3_cdf50.imag()));
          return ret_cdf35;
        }
      }
      
      if (verbose>0) {
        std::cout << "Failed 4: " << ret_cdf35 << " " << ret_cdf50 
                  << "\n  " << dtos(x1_cdf35,0) << " "
                  << dtos(x1_cdf50,0) << " "
                  << "\n  (" << dtos(x2_cdf35.real(),0) << ","
                  << dtos(x2_cdf35.imag(),0) << ") ("
                  << dtos(x2_cdf50.real(),0) << ","
                  << dtos(x2_cdf50.imag(),0) << ")"
                  << "\n  (" << dtos(x3_cdf35.real(),0) << ","
                  << dtos(x3_cdf35.imag(),0) << ") ("
                  << dtos(x3_cdf50.real(),0) << ","
                  << dtos(x3_cdf50.imag(),0) << ")"
                  << "\n  " << dtos(err,0) << " " << tol_loc << std::endl;
      }
        
      /// Final pass, compare 50- and 100-digit precision
      
      cpp_dec_float_100 a_cdf100=static_cast<cpp_dec_float_100>(a);
      cpp_dec_float_100 b_cdf100=static_cast<cpp_dec_float_100>(b);
      cpp_dec_float_100 c_cdf100=static_cast<cpp_dec_float_100>(c);
      cpp_dec_float_100 d_cdf100=static_cast<cpp_dec_float_100>(d);
      cpp_dec_float_100 x1_cdf100;
      std::complex<cpp_dec_float_100> x2_cdf100, x3_cdf100;
      
      int ret_cdf100=q_cdf100.solve_rc(a_cdf100,b_cdf100,c_cdf100,d_cdf100,
                                       x1_cdf100,x2_cdf100,x3_cdf100);

      if (ret_cdf50==ret_cdf100) {
        // If d is 0 and the first root is near 0, then avoid dividing
        // by it
        err=0;
        if (d==0 && (ret_cdf100==1 ||
                     (abs(x1_cdf100)<abs(x2_cdf100) &&
                      abs(x1_cdf100)<abs(x3_cdf100)))) {
          x1_cdf100=0;
        } else {
          err=static_cast<fp_t>(abs(x1_cdf100-x1_cdf50)/abs(x1_cdf100));
        }
        if (d==0 && abs(x2_cdf100)<abs(x1_cdf100) &&
            abs(x2_cdf100)<abs(x3_cdf100)) {
          x2_cdf100=0;
        } else if (x2_cdf100.real()!=0 || x2_cdf100.imag()!=0 ||
                   x2_cdf50.real()!=0 || x2_cdf50.imag()!=0) {
          err+=static_cast<fp_t>(abs(abs(x2_cdf100)-abs(x2_cdf50))/
                                 abs(x2_cdf100));
        }
        if (d==0 && abs(x3_cdf100)<abs(x1_cdf100) &&
            abs(x3_cdf100)<abs(x2_cdf100)) {
          x3_cdf100=0;
        } else if (x3_cdf100.real()!=0 || x3_cdf100.imag()!=0 ||
                   x3_cdf50.real()!=0 || x3_cdf50.imag()!=0) {
          err+=static_cast<fp_t>(abs(abs(x3_cdf100)-abs(x3_cdf50))/
                                 abs(x3_cdf100));
                                     
        }
        if (err<tol_loc) {
          x1=static_cast<fp_t>(x1_cdf100);
          x2.real(static_cast<fp_t>(x2_cdf100.real()));
          x2.imag(static_cast<fp_t>(x2_cdf100.imag()));
          x3.real(static_cast<fp_t>(x3_cdf100.real()));
          x3.imag(static_cast<fp_t>(x3_cdf100.imag()));
          return ret_cdf50;
        }
      }
      
      if (verbose>0) {
        std::cout << "Failed 5: " << ret_cdf50 << " " << ret_cdf100 
                  << "\n  " << dtos(x1_cdf50,0) << " "
                  << dtos(x1_cdf100,0) << " "
                  << "\n  (" << dtos(x2_cdf50.real(),0) << ","
                  << dtos(x2_cdf50.imag(),0) << ") ("
                  << dtos(x2_cdf100.real(),0) << ","
                  << dtos(x2_cdf100.imag(),0) << ")"
                  << "\n  (" << dtos(x3_cdf50.real(),0) << ","
                  << dtos(x3_cdf50.imag(),0) << ") ("
                  << dtos(x3_cdf100.real(),0) << ","
                  << dtos(x3_cdf100.imag(),0) << ")"
                  << "\n  " << dtos(err,0) << " " << tol_loc << std::endl;
      }
    
      /// Algorithm failed
      
      O2SCL_CONV2("Failed to compute with requested accuracy ",
                  "in cubic_real_coeff_multip::solve_rc_tol_err().",
                  o2scl::exc_efailed,err_nonconv);
      return o2scl::exc_efailed;
    }

    /** \brief Evaluate the function and return the error estimate
        with the default tolerance for the specified type
     */
    template<class fp_t, class cx_t>
    int solve_rc_err(const fp_t &a, const fp_t &b,
                     const fp_t &c, const fp_t &d,
                     fp_t &x1, cx_t &x2, cx_t &x3,
                     fp_t &err) {
      return solve_rc_tol_err(a,b,c,d,x1,x2,x3,err);
    }
  
    /** \brief Evalulate the function without an error estimate
     */
    template<class fp_t, class cx_t>
    int solve_rc(const fp_t &a, const fp_t &b,
                         const fp_t &c, const fp_t &d,
                         fp_t &x1, cx_t &x2, cx_t &x3) {
      fp_t err;
      return solve_rc_err(a,b,c,d,x1,x2,x3,err);
    }

    /** \brief Compute the cubic discriminant, 
	\f$ b^2 c^2 - 4 a c^3 - 4 b^3 d - 27 a^2 d^2 + 18 a b c d \f$

        If the discriminant is zero, then all roots qre real and
        at least two are equal (possibly all three are identical).
        If the discriminant is positive, then there are 
        three distinct real roots, and if the discriminant is negative
        then there is one real root and two complex conjugate
        roots. 
     */
    template<class fp_t>
    fp_t disc3_r(const fp_t a3, const fp_t b3, const fp_t c3, 
                        const fp_t d3) {
      return b3*b3*c3*c3-4.0*a3*c3*c3*c3-4.0*b3*b3*b3*d3-
	27.0*a3*a3*d3*d3+18.0*a3*b3*c3*d3;
    }

    
  };

  /** \brief Solve a general polynomial with real coefficients (GSL)
   */
  template<class fp_t=double, class cx_t=std::complex<fp_t>,
           class coeff_vec_t=std::vector<double>,
           class root_vec_t=std::vector<std::complex<double> > >
  class poly_real_coeff_gsl : public poly_real_coeff<fp_t,cx_t> {

  public:

    /// Check each root and automatically refine
    bool check_refine;
    
    poly_real_coeff_gsl() {
      w2=gsl_poly_complex_workspace_alloc(3);
      w3=gsl_poly_complex_workspace_alloc(4);
      w4=gsl_poly_complex_workspace_alloc(5);
      gen_size=0;
      check_refine=false;
    }

    virtual ~poly_real_coeff_gsl() {
      gsl_poly_complex_workspace_free(w2);
      gsl_poly_complex_workspace_free(w3);
      gsl_poly_complex_workspace_free(w4);
      if (gen_size>0) {
        gsl_poly_complex_workspace_free(wgen);
      }
    }

    /** \brief Solve a generic polynomial given <tt>n+1</tt> coefficients

	\note In order to be consistent with the other solve_rc()
	functions, the ordering of the coefficients is reversed with
	respect to gsl_poly_complex_solve(). The leading coefficient
	is stored in <tt>co[0]</tt> and the constant term is stored in
	<tt>co[n]</tt>.
     */
    virtual int solve_rc_arr(int n, const coeff_vec_t &co,
                             root_vec_t &ro) {
      
      int j;
      // We're sending a pointer to the first element, so
      // this must be a std::vector<double>
      std::vector<double> a(n+1), z(2*n);
      cx_t i(0.0,1.0);
      
      for(j=0;j<n+1;j++) {
        a[j]=co[n-j];
      }  
      if (gen_size!=n) {
        if (gen_size>0) {
          gsl_poly_complex_workspace_free(wgen);
        }
        wgen=gsl_poly_complex_workspace_alloc(n+1);
        gen_size=n;
      }
      if (a[n]==0.0) {
        O2SCL_ERR2("Leading coefficient zero in ",
                   "poly_real_coeff_gsl::solve_rc().",
                   exc_einval);
      }
      gsl_poly_complex_solve(&a[0],n+1,wgen,&z[0]);
      
      for(j=0;j<n;j++) {
        ro[j]=z[2*j]+i*z[2*j+1];
      }

      return success;
    }      

    /** \brief Solve a quadratic polynomial with real coefficients
     */
    virtual int solve_rc(const fp_t a2, const fp_t b2, const fp_t c2, 
			 cx_t &r1, cx_t &r2) {
			 
      if (a2==0.0) {
        O2SCL_ERR2("Leading coefficient zero in ",
		   "poly_real_coeff_gsl::solve_rc().",
		   exc_einval);
      }
      
      fp_t a[3]={c2,b2,a2};
      fp_t z[4];
      cx_t i(0.0,1.0);

      gsl_poly_complex_solve(a,3,w2,z);
      
      int nreal=0;
      if (z[1]==0) nreal++;
      if (z[3]==0) nreal++;
      
      r1=z[0]+i*z[1];
      r2=z[2]+i*z[3];
      
      return nreal;
    }
    
    /** \brief Solve a cubic polynomial with real coefficients
     */
    virtual int solve_rc(const fp_t a3, const fp_t b3, const fp_t c3, 
			 const fp_t d3, fp_t &r1, cx_t &r2, 
			 cx_t &r3) {
      
      if (a3==0.0) {
        O2SCL_ERR2("Leading coefficient zero in ",
		   "poly_real_coeff_gsl::solve_rc().",
		   exc_einval);
      }

      fp_t a[4]={d3,c3,b3,a3};  
      fp_t z[6],s1,s2,s3;
      cx_t i(0.0,1.0);
      
      gsl_poly_complex_solve(a,4,w3,z);
      
      s1=fabs(z[1]);
      s2=fabs(z[3]);
      s3=fabs(z[5]);

      if (s1==0 && s2==0 && s3==0) {
        r1=z[0];
        r2=z[2];
        r3=z[4];
        return 3;
      }
      
      if (s1<s2 && s1<s3) {
        r1=z[0];
        r2=z[2]+i*z[3];
        r3=z[4]+i*z[5];
      } else if (s2<s1 && s2<s3) {
        r1=z[2];
        r2=z[0]+i*z[1];
        r3=z[4]+i*z[5];
      } else {
        r1=z[4];
        r2=z[0]+i*z[1];
        r3=z[2]+i*z[3];
      }

      if (check_refine) {
        coeff_vec_t co(4);
        co[0]=a3;
        co[1]=b3;
        co[2]=c3;
        co[3]=d3;
        root_vec_t ro(3);
        ro[0]=r1;
        ro[1]=r2;
        ro[2]=r3;
        
        this->polish_rc_arr(3,co,ro);
      }
      
      return 1;
    }
    
    /** \brief Solve a quartic polynomial with real coefficients
     */
    virtual int solve_rc(const fp_t a4, const fp_t b4, const fp_t c4, 
			 const fp_t d4, const fp_t e4, 
			 cx_t &r1, cx_t &r2, 
			 cx_t &r3, cx_t &r4) {
      if (a4==0.0) {
        O2SCL_ERR2("Leading coefficient zero in ",
                   "poly_real_coeff_gsl::solve_rc().",
                   exc_einval);
      }
      
      fp_t a[5]={e4,d4,c4,b4,a4};
      fp_t z[8];
      cx_t i(0.0,1.0);
      
      gsl_poly_complex_solve(a,5,w4,z);

      int nreal=0;
      if (z[1]==0) nreal++;
      if (z[3]==0) nreal++;
      if (z[5]==0) nreal++;
      if (z[7]==0) nreal++;
      
      r1=z[0]+i*z[1];
      r2=z[2]+i*z[3];
      r3=z[4]+i*z[5];
      r4=z[6]+i*z[7];
      
      return nreal;
    }

    /// Return a string denoting the type ("poly_real_coeff_gsl")
    const char *type() { return "poly_real_coeff_gsl"; }

  protected:

#ifndef DOXYGEN_INTERNAL

    /// Workspace for quadratic polynomials
    gsl_poly_complex_workspace *w2;

    /// Workspace for cubic polynomials
    gsl_poly_complex_workspace *w3;

    /// Workspace for quartic polynomials
    gsl_poly_complex_workspace *w4;

    /// Workspace for general polynomials
    gsl_poly_complex_workspace *wgen;

    /// The size of the workspace \ref wgen
    int gen_size;

#endif

  };

  /** \brief Solve a quadratic with complex coefficients and complex roots
   */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class quadratic_complex_std : public quadratic_complex<fp_t,cx_t> {

  public:

    virtual ~quadratic_complex_std() {}

    /** \brief Solves the complex polynomial \f$ a_2 x^2 + b_2 x + c_2 = 0 \f$ 
	giving the two complex solutions \f$ x=x_1 \f$ and \f$ x=x_2 \f$ 
    */
    virtual int solve_c(const cx_t a2, const cx_t b2, 
			const cx_t c2, cx_t &x1, cx_t &x2) {
      
      if (a2.real()==0.0 && a2.imag()==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quadratic_complex_std::solve_c().",
           exc_einval);
      }
      // AWS 4/1/21: These temporaries fix compiling for long double
      // and complex<long double> types
      fp_t two=2.0;
      fp_t four=4.0;
      x1=(-b2+sqrt(b2*b2-four*a2*c2))/two/a2;
      x2=(-b2-sqrt(b2*b2-four*a2*c2))/two/a2;
      return success;
    }

    /// Return a string denoting the type ("quadratic_complex_std")
    const char *type() { return "quadratic_complex_std"; }
  };

  /** \brief Solve a cubic with complex coefficients and complex roots
   */
  template<class fp_t=double, class cx_t=std::complex<fp_t> >
  class cubic_complex_std : public cubic_complex<fp_t,cx_t> {

  public:

    virtual ~cubic_complex_std() {}
    
    /** \brief Solves the complex polynomial 
	\f$ a_3 x^3 + b_3 x^2 + c_3 x + d_3= 0 \f$ 
	giving the three complex solutions \f$ x=x_1 \f$ , 
	\f$ x=x_2 \f$ , and \f$ x=x_3 \f$ .
    */
    virtual int solve_c(const cx_t a3, const cx_t b3, 
			const cx_t c3, const cx_t d3, 
			cx_t &x1, cx_t &x2, cx_t &x3) {

      if (a3.real()==0.0 && a3.imag()==0.0) {
        quadratic_complex_std<fp_t,cx_t> qsc;
        qsc.solve_c(b3,c3,d3,x1,x2);
        x3=0.0;
        return success;
      }
      
      cx_t p3, q3, mo(0,1);
      cx_t alpha, beta, cbrta, cbrtb;
      cx_t e2, e4, cacb;
      fp_t test, re_p3;

      fp_t two=2.0;
      fp_t three=3.0;
      fp_t four=4.0;
      fp_t nine=9.0;
      fp_t twoseven=27.0;
      fp_t onethird=nine/twoseven;
      
      fp_t pi=boost::math::constants::pi<fp_t>();

      if (false) {

        // AWS 4/4/21: This code appeared to be slightly more accurate
        // on average than the code below, but it occasionally gave
        // incorrect roots. 
        
        p3=(three*a3*c3-b3*b3)/nine/a3/a3;
        q3=(two*b3*b3*b3-nine*a3*b3*c3+twoseven*a3*a3*d3)/twoseven/a3/a3/a3;
        std::cout << "p,q: " << p3 << " " << q3 << std::endl;
        
        alpha=(-q3+sqrt(q3*q3+four*p3*p3*p3))/two;
        beta=(q3+sqrt(q3*q3+four*p3*p3*p3))/two;
        std::cout << "alpha,beta: " << alpha << " " << beta << std::endl;
        
        if (alpha.real()==0.0) cbrta=0.0;
        else cbrta=pow(alpha,onethird);
        if (beta.real()==0.0) cbrtb=0.0;
        else cbrtb=pow(beta,onethird);
        
        cx_t exp_i_pi_three=exp(mo*pi/three);
        
        // It seems that if the real part of alpha is < 0 and the imaginary
        // part is zero, then cbrta is NaN (esp. w/Cygwin). We fix this
        // here:
        if (!o2isfinite(cbrta.real())) {
          cbrta=pow(-alpha,onethird)*exp_i_pi_three;
        }
        if (!o2isfinite(cbrtb.real())) {
          cbrtb=pow(-beta,onethird)*exp_i_pi_three;
        }
        
        e2=exp_i_pi_three*exp_i_pi_three;
        e4=e2*e2;
        
        std::cout << "alpha,beta: " << cbrta << " " << cbrtb << std::endl;
        
        // This next section is nessary to ensure that the code
        // selects the correct cube roots of alpha and beta.
        // I changed this because I wanted code that had no chance
        // of accidentally falling into an infinite loop.
        re_p3=p3.real();
        cacb=cbrta*cbrtb;
        test=fabs((cacb.real()-re_p3)/re_p3);
        if (fabs(((cacb*e2).real()-re_p3)/re_p3)<test) {
          cbrta*=e2;
          cacb=cbrta*cbrtb;
          test=fabs((cacb.real()-re_p3)/re_p3);
          if (fabs(((cacb*e2).real()-re_p3)/re_p3)<test) {
            cbrta*=e2;
            std::cout << "1" << std::endl;
          }
          std::cout << "2" << std::endl;
        } else {
          if (fabs(((cacb*e2*e2).real()-re_p3)/re_p3)<test) {
            cbrta*=e2*e2;
            std::cout << "3" << std::endl;
          }
        }
        
        x1=cbrta-cbrtb-b3/three/a3;
        x2=cbrta*e2-cbrtb*e4-b3/three/a3;
        x3=cbrta*e4-cbrtb*e2-b3/three/a3;
        std::cout << "x1,x2,x3: " << x1 << " " << x2 << " " << x3
                  << std::endl;
        
      }

      if (true) {

        // Following
        // https://en.wikipedia.org/wiki/Cubic_equation#Vieta's_substitution

        // Write the cubic as a "depressed" cubic: x^3+px+q=0
        cx_t p=(three*a3*c3-b3*b3)/three/a3/a3;
        cx_t q=(two*b3*b3*b3-nine*a3*b3*c3+twoseven*a3*a3*d3)/
          twoseven/a3/a3/a3;

        //std::cout << "p,q: " << p << " " << q << std::endl;
        
        cx_t exp_i_pi_three=exp(mo*pi*two/three);

        // Roots of the depressed cubic
        cx_t r1, r2, r3;

        // AWS 4/11/21: I previously hacked a comparison here between
        // 'one2' and 'one', but I didn't find it made any improvement.
        //
        // Check if p is too small:
        // fp_t one=1.0+abs(p);
        // fp_t one2=1.0;
        //
        // if (one2==one || (p.real()==0.0 && p.imag()==0.0)) {
        
        if (p.real()==0.0 && p.imag()==0.0) {
          
          // If p is zero, then the roots of the depressed
          // cubic are trivial
          
          r1=pow(-q,onethird);
          r2=r1*exp_i_pi_three;
          r3=r2*exp_i_pi_three;
          
        } else {
          
          // Formulate into a quadratic and then W is one of
          // the roots of the quadratic
          cx_t W=-q/two+sqrt(p*p*p/twoseven+q*q/four);
          //std::cout << "W: " << W << std::endl;
          
          // Construct the roots of the depressed cubic from the
          // solution of the quadratic
          cx_t w1=pow(W,onethird);
          cx_t w2=w1*exp_i_pi_three;
          cx_t w3=w2*exp_i_pi_three;
          
          // Vieta's substitution
          if (W.real()==0.0 && W.imag()==0.0) {
            // If W=0, then there are three identical roots of the
            // original cubic, and they are all -b3/a3/3
            r1=0;
            r2=0;
            r3=0;
          } else {
            r1=w1-p/three/w1;
            r2=w2-p/three/w2;
            r3=w3-p/three/w3;
          }
        }

        // Construct the roots of the original cubic from those of the
        // depressed cubic
        x1=r1-b3/three/a3;
        x2=r2-b3/three/a3;
        x3=r3-b3/three/a3;

      }

      if (false) {
        cx_t check1=a3*x1*x1*x1+b3*x1*x1+c3*x1+d3;
        cx_t check2=a3*x2*x2*x2+b3*x2*x2+c3*x2+d3;
        cx_t check3=a3*x3*x3*x3+b3*x3*x3+c3*x3+d3;
        if (abs(check1)>1.0e-8) {
          std::cout << "ccheck: " << check1 << " " << check2 << " "
                    << check3 << std::endl;
          std::cout << x1 << " " << x2 << " " << x3
                    << std::endl;
          std::cout << a3 << " " << b3 << " " << c3 << " "
                    << d3 << std::endl;
          exit(-1);
        }
      }
      
      return success;
    }      
    
    /// Return a string denoting the type ("cubic_complex_std")
    const char *type() { return "cubic_complex_std"; }
  };

  /** \brief Solve a quartic with real coefficients and real roots
   */
  template<class fp_t=double> 
  class quartic_real_std : public quartic_real<fp_t> {
    
  protected:

    cubic_real_coeff_gsl2<fp_t,std::complex<fp_t> > cub2;
    quadratic_real_coeff_gsl2<fp_t,std::complex<fp_t> > quad2;
    
  public:

    quartic_real_std() {
    }

    virtual ~quartic_real_std() {}

    /** \brief Desc
     */
    virtual int solve_r(const fp_t a4, const fp_t b4, const fp_t c4, 
			const fp_t d4, const fp_t e4, fp_t &x1, 
			fp_t &x2, fp_t &x3, fp_t &x4) {
      if (a4==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_real_std::solve_r().",
           exc_einval);
      }
      
      fp_t a34=b4/a4;
      fp_t a24=c4/a4;
      fp_t a14=d4/a4;
      fp_t a04=e4/a4;
      
      //---------------------------------------
      // Solve the resolvent cubic:
      
      fp_t a23=-a24;
      fp_t a13=(a14*a34-4.0*a04);
      fp_t a03=-(a14*a14+a04*a34*a34-4.0*a04*a24);
      
      fp_t u1, u2, u3;
      cub2.solve_r(1,a23,a13,a03,u1,u2,u3);

      // FIXME: if the cubic above has only one real root,
      // then u2 will be uninitialized
      fp_t u4=u2;
      
      //---------------------------------------
      // Now construct the two quadratics:
      
      fp_t t1=u4+a34*a34/4-a24;
      if (t1>0.0) {
        t1=sqrt(t1);
      } else {
        t1=0.0;
      }
      
      fp_t b2a=-t1+a34/2;
      fp_t b2b=t1+a34/2;
      t1=u4*u4/4;
      
      // When numerical errors make t1 slightly smaller than a04.
      if (t1>a04) {
        t1=sqrt(t1-a04);
      } else {
        t1=0;
      }
      
      fp_t c2a=u4/2-t1;
      fp_t c2b=u4/2+t1;
      
      if (fabs((b2a*c2b+c2a*b2b-d4)/d4)>1.0e-4) {
        t1=u4+a34*a34/4-a24;
        if (t1>0) {
          t1=-sqrt(t1);
        } else {
          t1=0;
        }
        
        b2a=-t1+a34/2;
        b2b=t1+a34/2;
        
        t1=u4*u4/4;
        if (fabs((u4*u4/4-a04)/a04)<1.0e-6) {
          t1=0.0;
        } else {
          t1=sqrt(t1-a04);
        }
        c2a=u4/2-t1;
        c2b=u4/2+t1;
      }
      
      //---------------------------------------
      // The solutions to the two quadratics:

      int ir1=quad2.solve_r(1,b2a,c2a,x1,x2);
      int ir2=quad2.solve_r(1,b2b,c2b,x3,x4);

      if (ir1>=2 && (!o2isfinite(x1) || !o2isfinite(x2))) {
        /*
          std::cout << "Zere1: " << std::endl;
          std::cout << u1 << " " << u2 << " " << u3 << std::endl;
          std::cout << t1 << " " << a04 << " " <<  a34 << std::endl;
          std::cout << "b2a,b2b: " << b2a << " " << b2b << std::endl;
          std::cout << "c2a,c2b: " << c2a << " " << c2b << std::endl;
          std::cout << x1 << " " << x2 << std::endl;
          exit(-1);
        */
        O2SCL_ERR("Quartic failure 1 in quartic_real_std.",
                  o2scl::exc_einval);
      }
      if (ir2>=2 && (!o2isfinite(x3) || !o2isfinite(x4))) {
        /*
          std::cout << "Zere1: " << std::endl;
          std::cout << u1 << " " << u2 << " " << u3 << std::endl;
          std::cout << t1 << " " << a04 << " " <<  a34 << std::endl;
          std::cout << "b2a,b2b: " << b2a << " " << b2b << std::endl;
          std::cout << "c2a,c2b: " << c2a << " " << c2b << std::endl;
          std::cout << x3 << " " << x4 << std::endl;
          exit(-1);
        */
        O2SCL_ERR("Quartic failure 2 in quartic_real_std.",
                  o2scl::exc_einval);
      }

      if (ir1==0 && ir2==0) {
        return 0;
      } else if (ir1==0) {
        x1=x3;
        x2=x4;
        return 2;
      }

      return success;
    }

    /// Return a string denoting the type ("quartic_real_std")
    const char *type() { return "quartic_real_std"; }

  };
  
  /** \brief Solve a quartic with complex coefficients and complex roots
   */
  template <class fp_t=double, class cx_t=std::complex<fp_t> >
  class quartic_complex_std : public quartic_complex<fp_t,cx_t> {

  public:

    virtual ~quartic_complex_std() {}

    /** \brief Solves the complex polynomial 
	\f$ a_4 x^4 + b_4 x^3 + c_4 x^2 + d_4 x + e_4 = 0 \f$ 
	giving the four complex solutions \f$ x=x_1 \f$ , \f$ x=x_2 \f$ ,
	\f$ x=x_3 \f$ , and \f$ x=x_4 \f$ .
    */
    virtual int solve_c(const cx_t a4, const cx_t b4, 
			const cx_t c4, const cx_t d4, 
			const cx_t e4, cx_t &x1, cx_t &x2, 
			cx_t &x3, cx_t &x4) {
      
      cx_t p4, q4, r4;
      cx_t a3, b3, c3, d3;
      cx_t b2a, c2a, b2b, c2b;
      cx_t u4, u41, u42;
      
      if (a4.real()==0.0 && a4.imag()==0.0) {
        O2SCL_ERR
          ("Leading coefficient zero in quartic_complex_std::solve_c().",
           exc_einval);
      }

      if (c4.real()==0.0 && c4.imag()==0.0 &&
          d4.real()==0.0 && d4.imag()==0.0 &&
          e4.real()==0.0 && e4.imag()==0.0) {
        x1=-b4/a4;
        x2=0;
        x3=0;
        x4=0;
        return 0;
      }

      fp_t two=2.0;
      fp_t eight=8.0;
      fp_t four=4.0;
      fp_t sixteen=16.0;
      fp_t sixfour=64.0;
      fp_t twofivesix=256.0;
      fp_t three=3.0;

      // Construct the depressed quartic
      p4=(eight*a4*c4-three*b4*b4)/eight/a4/a4;
      q4=(b4*b4*b4-four*a4*b4*c4+eight*a4*a4*d4)/eight/(a4*a4*a4);
      r4=(sixteen*a4*b4*b4*c4+twofivesix*a4*a4*a4*e4-three*b4*b4*b4*b4-
          sixfour*a4*a4*b4*d4)/twofivesix/(a4*a4*a4*a4);

      // If q4 is zero, then the depressed quartic
      // is a biquadratic
      if (q4.real()==0.0 && q4.imag()==0.0) {
        
        cx_t z1, z2;
        quad_obj.solve_c(1.0,p4,r4,z1,z2);
        
        x1=sqrt(z1)-b4/four/a4;
        x2=-sqrt(z1)-b4/four/a4;
        x3=sqrt(z2)-b4/four/a4;
        x4=-sqrt(z2)-b4/four/a4;
        
        cx_t check1=a4*x1*x1*x1*x1+b4*x1*x1*x1+c4*x1*x1+d4*x1+e4;
        cx_t check2=a4*x2*x2*x2*x2+b4*x2*x2*x2+c4*x2*x2+d4*x2+e4;
        cx_t check3=a4*x3*x3*x3*x3+b4*x3*x3*x3+c4*x3*x3+d4*x3+e4;
        cx_t check4=a4*x4*x4*x4*x4+b4*x4*x4*x4+c4*x4*x4+d4*x4+e4;
        if (abs(check1)>1.0e-4 || abs(check2)>1.0e-4 ||
            abs(check3)>1.0e-4 || abs(check4)>1.0e-4) {
          std::cout << "Yere." << std::endl;
          std::cout << check1 << " " << check2 << " "
                    << check3 << " " << check4 << std::endl;
          exit(-1);
        }
        
        return 0;
      }
      
      if (p4.real()==0.0 && p4.imag()==0.0 &&
          q4.real()==0.0 && q4.imag()==0.0 &&
          r4.real()==0.0 && r4.imag()==0.0) {
        x1=-b4/four/a4;
        x2=-b4/four/a4;
        x3=-b4/four/a4;
        x4=-b4/four/a4;
        return 0;
      }
      
      //---------------------------------------
      // Solve the resolvent cubic:
      
      a3=1.0;
      b3=-p4;
      c3=-four*r4;
      d3=four*p4*r4-q4*q4;
      //std::cout << "a3,b3,c3,d3: " << a3 << " " << b3 << " "
      //<< c3 << " " << d3 << std::endl;
      
      cub_obj.solve_c(a3,b3,c3,d3,u4,u41,u42);

      /*
      cx_t check5=a3*u4*u4*u4+b3*u4*u4+c3*u4+d3;
      cx_t check6=a3*u41*u41*u41+b3*u41*u41+c3*u41+d3;
      cx_t check7=a3*u42*u42*u42+b3*u42*u42+c3*u42+d3;
      std::cout << "check5, check6, check7: "
                << check5 << " " << check6 << " " << check7
                << std::endl;
      */
      
      //---------------------------------------
      
      //std::cout << "u4,u41,u42: " << u4 << " "
      //<< u41 << " " << u42 << std::endl;

      // AWS, 4/8/21: I'd prefer not to have a comparison to a number
      // in here, but I find that sometimes u4 and p4 are extremely
      // close but not exactly equal and this causes failures in
      // computing the correct roots.
      if (u4==p4 || abs(u4-p4)<1.0e-15) {
        b2a=sqrt(u41-p4);
        b2b=-sqrt(u41-p4);
        c2a=-sqrt(u41-p4)*q4/two/(u41-p4)+u41/two;
        c2b=sqrt(u41-p4)*q4/two/(u41-p4)+u41/two;
        /*
          b2a=0.0;
          b2b=0.0;
          c2a=u4/two;
          c2b=u4/two;
        */
        //std::cout << "Here." << std::endl;
      } else {
        b2a=sqrt(u4-p4);
        b2b=-sqrt(u4-p4);
        c2a=-sqrt(u4-p4)*q4/two/(u4-p4)+u4/two;
        c2b=sqrt(u4-p4)*q4/two/(u4-p4)+u4/two;
      }
      //std::cout << "u4-p4: " << u4-p4 << " " << abs(u4-p4) << std::endl;
      
      x1=(-b2a+sqrt(b2a*b2a-four*c2a))/two-b4/four/a4;
      x2=(-b2a-sqrt(b2a*b2a-four*c2a))/two-b4/four/a4;
      x3=(-b2b+sqrt(b2b*b2b-four*c2b))/two-b4/four/a4;
      x4=(-b2b-sqrt(b2b*b2b-four*c2b))/two-b4/four/a4;

      //std::cout << "b2a,b2b,c2a,c2b: "
      //<< b2a << " " << b2b << " " << c2a << " "
      //<< c2b << std::endl;

      if (false) {
        cx_t check1=a4*x1*x1*x1*x1+b4*x1*x1*x1+c4*x1*x1+d4*x1+e4;
        cx_t check2=a4*x2*x2*x2*x2+b4*x2*x2*x2+c4*x2*x2+d4*x2+e4;
        cx_t check3=a4*x3*x3*x3*x3+b4*x3*x3*x3+c4*x3*x3+d4*x3+e4;
        cx_t check4=a4*x4*x4*x4*x4+b4*x4*x4*x4+c4*x4*x4+d4*x4+e4;
        if (abs(check1)>1.0e-4 || abs(check2)>1.0e-4 ||
            abs(check3)>1.0e-4 || abs(check4)>1.0e-4) {
          std::cout << "Xere." << std::endl;
          std::cout << "a4,b4,c4,d4,e4: " << a4 << " " << b4 << " "
                    << c4 << " " << d4 << " " << e4 << std::endl;
          std::cout << "a3,b3,c3,d3: " << a3 << " " << b3 << " "
                    << c3 << " " << d3 << std::endl;
          std::cout << "p4,q4,r4: " << p4 << " " << q4 << " "
                    << r4 << std::endl;
          std::cout << "u4,u41,u42: " << u4 << " "
                    << u41 << " " << u42 << std::endl;
          std::cout << "b2a,b2b,c2a,c2b: "
                    << b2a << " " << b2b << " " << c2a << " "
                    << c2b << std::endl;
          std::cout << "x1,x2,x3,x4: " << x1 << " " << x2 << " "
                    << x3 << " " << x4 << std::endl;
          std::cout << "check1,2,3,4: " << check1 << " " << check2 << " "
                    << check3 << " " << check4 << std::endl;
          exit(-1);
        }
      }
      //cx_t check1=a3*u4*u4*u4+b3*u4*u4+c3*u4+d3;
      //std::cout << "check1: " << check1 << std::endl;
      //cx_t check2=a3*u41*u41*u41+b3*u41*u41+c3*u41+d3;
      //std::cout << "check2: " << check2 << std::endl;
      //cx_t check3=a3*u42*u42*u42+b3*u42*u42+c3*u42+d3;
      //std::cout << "check3: " << check3 << std::endl;
      
      if (!o2isfinite(x1.real()) || !o2isfinite(x1.imag()) ||
          !o2isfinite(x2.real()) || !o2isfinite(x2.imag()) ||
          !o2isfinite(x3.real()) || !o2isfinite(x3.imag()) ||
          !o2isfinite(x4.real()) || !o2isfinite(x4.imag())) {
        std::cout << "Pt.1" << std::endl;
        std::cout << "a3,b3,c3,d3: " << a3 << " " << b3 << " "
                  << c3 << " " << d3 << std::endl;
        std::cout << "a4,b4,c4,d4,e4: " << a4 << " " << b4 << " "
                  << c4 << " " << d4 << " "
                  << e4 << std::endl;
        std::cout << "p4,q4,r4: " << p4 << " " << q4 << " "
                  << r4 << std::endl;
        std::cout << "x1,x2,x3,x4: "
                  << x1 << " " << x2 << " " << x3 << " " << x4
                  << std::endl;
        exit(-1);
      }
      
      return success;
    }

    /// Return a string denoting the type ("quartic_complex_std")
    const char *type() { return "quartic_complex_std"; }

  protected:

    /// The object to solve for the associated cubic
    cubic_complex_std<fp_t,cx_t> cub_obj;
    
    /// The object to solve for the associated cubic
    quadratic_complex_std<fp_t,cx_t> quad_obj;
    
  };

}

#endif
