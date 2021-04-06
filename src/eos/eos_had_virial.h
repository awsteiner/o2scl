/*
  -------------------------------------------------------------------
  
  Copyright (C) 2018-2021, Xingfu Du and Andrew W. Steiner
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#ifndef EOS_HAD_VIRIAL_H
#define EOS_HAD_VIRIAL_H

#include <cmath>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <o2scl/test_mgr.h>
#include <o2scl/mm_funct.h>
#include <o2scl/mroot_hybrids.h>
#include <o2scl/mroot_cern.h>
#include <o2scl/linear_solver.h>
#include <o2scl/poly.h>

namespace o2scl {

  typedef boost::numeric::ublas::vector<double> ubvector;
  typedef boost::numeric::ublas::matrix<double> ubmatrix;

  /** \brief Compute the virial EOS
   */
  class eos_had_virial {

  public:

    // Store the number of function and derivative evaluations
    int nf, nd;
    double nn, pn, lambda, T, b_n, b_pn;
    double mfn2_mu_n, mfn2_mu_p, dbndT, zp, zn;
    double dbpndT, dlambdadT, npt, nnt, a, b, c, d, e;
    ubmatrix A;
    ubvector B;

    /// Linear system solver
    o2scl_linalg::linear_solver_LU<ubvector,ubmatrix> lsol;  

    /// Quartic polynomial solver
    o2scl::quartic_real_coeff_cern<> quart;

    // Generic polynomial solver
    o2scl::poly_real_coeff_gsl<> quart2;
  
    /// Storage for the four roots
    std::complex<double> res_zp[4],res_zn[4]; 
    /** \brief Solve for the fugacities given the densities

	This function computes zn and zp from pn and nn
	presuming that lambda, b_n and b_pn have already been specified.
    */
    void solve_fugacity(ubvector &x) {

      npt=pow(lambda,3)/2.0*pn;
      nnt=pow(lambda,3)/2.0*nn;

      if (npt>=nnt) {
      
	// Coefficients for quartic equation of zp in descending order
      
	a=pow(b_n,3)*2.0/b_pn/b_pn-2.0*b_n;
	b=-1.0+b_n*b_n*2.0/b_pn/b_pn-b_n/b_pn;
	c=b_n/(2.0*b_pn*b_pn)-0.5/b_pn-nnt+npt-b_n*b_n*npt*2.0/b_pn/b_pn;
	d=-b_n*npt/b_pn/b_pn+npt/2.0/b_pn;
	e=b_n*npt*npt/2.0/b_pn/b_pn;
      
	quart2.solve_rc(a,b,c,d,e,res_zp[0],res_zp[1],res_zp[2],res_zp[3]);
	//std::cout << "Here1: " << a << " " << b << " " << c << " "
	//<< d << " " << e << std::endl;
	int root_count=0;
	std::complex<double> eval_zp[4];
	ubvector res_zn_real;
	std::complex<double>  eval_zn[4];
	res_zn_real.resize(4);
	/*
	  AWS: Changed on 1/9/18 to select largest fugacity instead
	  of exiting when multiple roots are found
	*/

	for (int i=0;i<4;i++) {

	  // Check that the root is positive and that the imaginary
	  // part is sufficiently small. Note I use fabs() rather
	  // than abs().
	  if(res_zp[i].real()>0 && 
	     fabs(res_zp[i].imag()/res_zp[i].real())<1.0e-6) {

	    // Make sure that the specified root is really a solution
	    // of the original polynomial
	    eval_zp[i]=(a*pow(res_zp[i],4.0)+b*pow(res_zp[i],3.0)+
			c*pow(res_zp[i],2.0)+d*res_zp[i]+e)/e;
	   
	    // Changed from 1e-8 to 1e-6 because at zero temperature no
	    // solutions
	   
	    if (fabs(eval_zp[i].real())<2.0e-6 &&
		fabs(eval_zp[i].imag())<1.0e-8) {
	      double r0, r1;
	      gsl_poly_solve_quadratic(2.0*b_n,
				       2.0*res_zp[i].real()*b_pn+1.0,
				       -nnt,&r0,&r1);
	      std::complex<double> eval_r0=(res_zp[i].real()+2.0
					    *res_zp[i].real()*
					    res_zp[i].real()*b_n
					    +2.0*res_zp[i].real()*
					    r0*b_pn-npt)/npt;
	      std::complex<double> eval_r1=(res_zp[i].real()+2.0
					    *res_zp[i].real()*
					    res_zp[i].real()*b_n
					    +2.0*res_zp[i].real()*
					    r1*b_pn-npt)/npt;
	      if (fabs(eval_r0.real())<2.0e-6 && r0>0.0) {
		res_zn_real[i]=r0;
		eval_zn[i]=eval_r0;
	      }
	      else if (fabs(eval_r1.real())<2.0e-6 && r1>0.0) {
		res_zn_real[i]=r1;
		eval_zn[i]=eval_r1;
	      } else {
		res_zn_real[i]=-1.0;
		//std::cout << "r0,r1:  "<< r0 << " " << r1 << std::endl;
		//std::cout << "eval_r0,r1: " << eval_r0.real() << " " 
		//			  << eval_r1.real() << std::endl;
	      }


            	                       
	      if (res_zn_real[i]>0.0) { 
		root_count++; 
	      }
	    } else {
	      res_zn_real[i]=1.0e10;
	    } 
	  } else {
            res_zn_real[i]=1.0e10;
	  }
	}
      
	if (root_count==0&&false) {
	  std::cout << "Zn/Zp zero roots: " << root_count 
		    << std::endl;
	  std::cout << "nn: " << nn << " pn: " << pn << " T: " 
		    << T << std::endl;
	  std::cout.setf(std::ios::showpos);
	  std::cout << "zn: " <<res_zn_real[0] << " " << " ";
	  std::cout << eval_zn[0].real() << " " << eval_zn[0].imag() 
		    << std::endl;
	  std::cout << "zp: " <<res_zp[0].real() << " " << res_zp[0].imag() 
		    << " ";
	  std::cout << eval_zp[0].real() << " " << eval_zp[0].imag() 
		    << std::endl;
	  std::cout << "zn: " <<res_zn_real[1] << " ";
	  std::cout << eval_zn[1].real() << " " << eval_zn[1].imag() 
		    << std::endl;
	  std::cout << "zp: " <<res_zp[1].real() << " " << res_zp[1].imag() 
		    << " ";
	  std::cout << eval_zp[1].real() << " " << eval_zp[1].imag() 
		    << std::endl;
	  std::cout << "zn: " <<res_zn_real[2] << " ";
	  std::cout << eval_zn[2].real() << " " << eval_zn[2].imag() 
		    << std::endl;
	  std::cout << "zp: " <<res_zp[2].real() << " " << res_zp[2].imag() 
		    << " ";
	  std::cout << eval_zp[2].real() << " " << eval_zp[2].imag() 
		    << std::endl;
	  std::cout << "zn: " <<res_zn_real[3] << " ";
	  std::cout << eval_zn[3].real() << " " << eval_zn[3].imag() 
		    << std::endl;    
	  std::cout << "zp: " <<res_zp[3].real() << " " << res_zp[3].imag() 
		    << " ";
	  std::cout << eval_zp[3].real() << " " << eval_zp[3].imag() 
		    << std::endl;          
	  std::cout.unsetf(std::ios::showpos);
	  O2SCL_ERR("Zero or more than one root in solve_fugacity().",
		    o2scl::exc_efailed);
	}
	int res_index=0;
	double minsq=1.0e100, temp;
	for (int i=0;i<4;i++) {
	  if (res_zn_real[i]>0.0 && res_zp[i].real()>0.0) {
	    temp = res_zn_real[i]*res_zn_real[i] + res_zp[i].real()
	      *res_zp[i].real();
	    if (temp<minsq) {
	      minsq=temp;
	      res_index=i;
	    }
	  }
	}
	zn=res_zn_real[res_index];
	zp=res_zp[res_index].real();  
	x[1]=log(zp)*T;
	x[0]=log(zn)*T; 

      } else {
    
     
	// Coefficients for quartic equation of zp in descending order
      
	a=pow(b_n,3)*2.0/b_pn/b_pn-2.0*b_n;
	b=-1.0+b_n*b_n*2.0/b_pn/b_pn-b_n/b_pn;
	c=b_n/(2.0*b_pn*b_pn)-0.5/b_pn-npt+nnt-b_n*b_n*nnt*2.0/b_pn/b_pn;
	d=-b_n*nnt/b_pn/b_pn+nnt/2.0/b_pn;
	e=b_n*nnt*nnt/2.0/b_pn/b_pn;
      
	quart2.solve_rc(a,b,c,d,e,res_zn[0],res_zn[1],res_zn[2],res_zn[3]);
	/*
	  std::cout << "Here2: " << a << " " << b << " " << c << " "
	  << d << " " << e << std::endl;
	  std::cout << res_zn[0] << " " << res_zn[1] << " " << res_zn[2] << " "
	  << res_zn[3] << std::endl;
	*/
	int root_count=0;
	std::complex<double> eval_zn[4];
	ubvector res_zp_real;
	std::complex<double> eval_zp[4];
	res_zp_real.resize(4);
	/*
	  AWS: Changed on 1/9/18 to select largest fugacity instead
	  of exiting when multiple roots are found
	*/

	for (int i=0;i<4;i++) {

	  // Check that the root is positive and that the imaginary
	  // part is sufficiently small. Note I use fabs() rather
	  // than abs().
	  if(res_zn[i].real()>0 && 
	     fabs(res_zn[i].imag()/res_zn[i].real())<1.0e-6) {

	    // Make sure that the specified root is really a solution
	    // of the original polynomial
	    eval_zn[i]=(a*pow(res_zn[i],4.0)+b*pow(res_zn[i],3.0)+
			c*pow(res_zn[i],2.0)+d*res_zn[i]+e)/e;
	   
	    // Changed from 1e-8 to 1e-6 because at zero temperature no
	    // solutions
	   
	    if (fabs(eval_zn[i].real())<2.0e-6 &&
		fabs(eval_zn[i].imag())<1.0e-8) {

	      double r0, r1;
	      gsl_poly_solve_quadratic(2.0*b_n,
				       2.0*res_zn[i].real()*b_pn+1.0,
				       -npt,&r0,&r1);
	      std::complex<double> eval_r0=(res_zn[i].real()+2.0
					    *res_zn[i].real()*
					    res_zn[i].real()*b_n
					    +2.0*res_zn[i].real()*
					    r0*b_pn-nnt)/nnt;
	      std::complex<double> eval_r1=(res_zn[i].real()+2.0
					    *res_zn[i].real()*
					    res_zn[i].real()*b_n
					    +2.0*res_zn[i].real()*
					    r1*b_pn-nnt)/nnt;
	      if (fabs(eval_r0.real())<2.0e-6 && r0>0.0) {
		res_zp_real[i]=r0;
		eval_zp[i]=eval_r0;
	      } else if (fabs(eval_r1.real())<2.0e-6 && r1>0.0) {
		res_zp_real[i]=r1;
		eval_zp[i]=eval_r1;
	      } else {
		res_zp_real[i]=-1.0;
		//std::cout << "r0,r1:  "<< r0 << " " << r1 << std::endl;
		//std::cout << "eval_r0,r1: " << eval_r0.real() << " "
		//                            << eval_r1.real() << std::endl;
	      }
  
	      if (res_zp_real[i]>0.0) { 
		root_count++; 
	      }
	    } else {
	      res_zp_real[i]=1.0e10;
	    } 
	  } else {
            res_zp_real[i]=1.0e10;
	  }
	}
      
	if (root_count==0&&false) {
	  std::cout << "Zn/Zp zero roots: " << root_count 
		    << std::endl;
	  std::cout << "nn: " << nn << " pn: " << pn << " T: " 
		    << T << std::endl;
	  std::cout.setf(std::ios::showpos);
	  std::cout << "zp: " <<res_zp_real[0] << " ";
	  std::cout << eval_zp[0].real() << " " << eval_zp[0].imag() 
		    << std::endl;
	  std::cout << "zn: " <<res_zn[0].real() << " " << res_zn[0].imag() 
		    << " ";
	  std::cout << eval_zn[0].real() << " " << eval_zn[0].imag() 
		    << std::endl;
	  std::cout << "zp: " <<res_zp_real[1] << " ";
	  std::cout << eval_zp[1].real() << " " << eval_zp[1].imag() 
		    << std::endl;
	  std::cout << "zn: " <<res_zn[1].real() << " " << res_zn[1].imag() 
		    << " ";
	  std::cout << eval_zn[1].real() << " " << eval_zn[1].imag() 
		    << std::endl;
	  std::cout << "zp: " <<res_zp_real[2] << " ";
	  std::cout << eval_zp[2].real() << " " << eval_zp[2].imag() 
		    << std::endl;
	  std::cout << "zn: " <<res_zn[2].real() << " " << res_zn[2].imag() 
		    << " ";
	  std::cout << eval_zn[2].real() << " " << eval_zn[2].imag() 
		    << std::endl;
	  std::cout << "zp: " <<res_zp_real[3] << " ";
	  std::cout << eval_zp[3].real() << " " << eval_zp[3].imag() 
		    << std::endl;       
	  std::cout << "zn: " <<res_zn[3].real() << " " << res_zn[3].imag() 
		    << " ";
	  std::cout << eval_zn[3].real() << " " << eval_zn[3].imag() 
		    << std::endl;          
	  std::cout.unsetf(std::ios::showpos);
	  O2SCL_ERR("Zero or more than one root in solve_fugacity().",
		    o2scl::exc_efailed);
	}
	int res_index=0;
	double minsq=1.0e100, temp;
	for (int i=0;i<4;i++) {
	  if (res_zp_real[i]>0.0 && res_zn[i].real()>0.0) {
	    temp = res_zp_real[i]*res_zp_real[i] + res_zn[i].real()
	      *res_zn[i].real();
	    if (temp<minsq) {
	      minsq=temp;
	      res_index=i;
	    }
	  }
	}
	zp=res_zp_real[res_index];
	zn=res_zn[res_index].real();  
	x[1]=log(zp)*T;
	x[0]=log(zn)*T; 
    
      }

      //std::cout << "zn,zp: " << zn << " " << zp << std::endl;
    
      return;
    }
   
    void mfn_e(ubvector &x) {
      npt=pow(lambda,3)/2*pn;
      //nnt=pow(lambda,3)/2*nn; pn equals nn;
      //mu_p is x[1] while mu_n is x[0]
      a=2*b_n+2*b_pn;
      b=1;
      c=-npt;
      d=(sqrt(b*b+4*a*c)-b)/2;
      x[0]=log(d)*T;
      zn=d;
      zp=d;
      x[1]=x[0];
      /*std::cout<<"------------------"<<std::endl;
	std::cout<<a<<" "<<b<<" "<<c<<" "<<d<<" "<<log(d)<<" "<<T<<std::endl;
	std::cout<<x[1]<<std::endl;
	std::cout<<"------------------"<<std::endl;*/
      return;
    }

    /** \brief Here, a brief description of this function
	derivative with respect to nn of mfn (linear solver)
    */
    double mfn21(ubvector &x2) {
      zn=exp(mfn2_mu_n/T);
      zp=exp(mfn2_mu_p/T);
      A.resize(2,2);
      B.resize(2);
      A(0,0)=2/pow(lambda,3)*(zn/T+4*zn*zn*b_n/T+2*zp*zn*b_pn/T);
      A(0,1)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
      A(1,0)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
      A(1,1)=2/pow(lambda,3)*(zp/T+4*zp*zp*b_n/T+2*zp*zn*b_pn/T);
      /*std::cout<<"mfn21 start: "<<std::endl;
	std::cout<<mfn2_mu_n<<" "<<mfn2_mu_p<<std::endl;
	std::cout<<zn<<" "<<zp<<" "<<T<<std::endl;
	std::cout<<"mfn21 end: "<<std::endl;*/
      /*std::cout<<A(0,0)<<" "<<A(0,1)<<" "<<A(1,0)<<" "
	<<A(1,1)<<std::endl;*/
      B(0)=1;
      B(1)=0;
      if (true) {
	double den=A(0,0)*A(1,1)-A(0,1)*A(1,0);
	x2[0]=A(1,1)/den;
	x2[1]=-A(1,0)/den;
      } else {
	lsol.solve(2,A,B,x2);
      }
      return 0;  
     
    }

    /** \brief Here, a brief description of this function
	derivative with respect to pn of mfn (linear solver)
    */
    double mfn31(ubvector &x3) {
      zn=exp(mfn2_mu_n/T);
      zp=exp(mfn2_mu_p/T);
      A.resize(2,2);
      B.resize(2);
      A(0,0)=2/pow(lambda,3)*(zn/T+4*zn*zn*b_n/T+2*zp*zn*b_pn/T);
      A(0,1)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
      A(1,0)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
      A(1,1)=2/pow(lambda,3)*(zp/T+4*zp*zp*b_n/T+2*zp*zn*b_pn/T);
      B(0)=0;
      B(1)=1;
      if (true) {
	double den=A(0,0)*A(1,1)-A(0,1)*A(1,0);
	x3[0]=-A(0,1)/den;
	x3[1]=A(0,0)/den;
      } else {
	lsol.solve(2,A,B,x3);
      }

      return 0;  
    }
   
    /** \brief Here, a brief description of this function
	derivative with respect to T of mfn
    */
    int mfn41(ubvector &x4) {
      // dmundT=x4[0];
      // dmupdT=x4[1];

      A(0,0)=2/pow(lambda,3)*(zn/T+4*zn*zn*b_n/T+2*zp*zn*b_pn/T);
      A(0,1)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
      A(1,0)=2/pow(lambda,3)*2*zp*zn*b_pn/T;
      A(1,1)=2/pow(lambda,3)*(zp/T+4*zp*zp*b_n/T+2*zp*zn*b_pn/T);
      B(0)=-(2/pow(lambda,4)*(-3)*dlambdadT*
	     (zn+2*zn*zn*b_n+
	      2*zp*zn*b_pn)+2/pow(lambda,3)*
	     (zn*(-mfn2_mu_n/T/T)+2*zn*zn*dbndT-
	      4*zn*zn*b_n*mfn2_mu_n/T/T+
	      2*zp*zn*dbpndT-
	      2*zp*zn*b_pn*mfn2_mu_p/T/T-
	      2*zp*zn*b_pn*mfn2_mu_n/T/T));
      B(1)=-(2/pow(lambda,4)*(-3)*dlambdadT*
	     (zp+2*zp*zp*b_n+
	      2*zp*zn*b_pn)+2/pow(lambda,3)*
	     (zp*(-mfn2_mu_p/T/T)+2*zp*zp*dbndT-
	      4*zp*zp*b_n*mfn2_mu_p/T/T+
	      2*zp*zn*dbpndT-
	      2*zp*zn*b_pn*mfn2_mu_p/T/T-
	      2*zp*zn*b_pn*mfn2_mu_n/T/T));
      nf++;
      if (true) {
	double den=A(0,0)*A(1,1)-A(0,1)*A(1,0);
	x4[0]=(A(1,1)*B(0)-A(0,1)*B(1))/den;
	x4[1]=(A(0,0)*B(1)-A(1,0)*B(0))/den;
      } else {
	lsol.solve(2,A,B,x4);
      }
      return 0;
    }

  };

  /** \brief Virial solver with derivatives
   */
  class eos_had_virial_deriv {
  
  protected:
  
    // Generic polynomial solver
    o2scl::poly_real_coeff_gsl<> quart;
  
    /// Storage for the four roots
    std::complex<double> res[4];

  public:

    /// \name First derivatives of the fugacities
    //@{
    double dzndnn;
    double dzndnp;
    double dzpdnn;
    double dzpdnp;
    double dzndT;
    double dzpdT;
    //@}
  
    /// \name Second derivatives of the fugacities
    //@{
    double d2zndnn2;
    double d2zndnndnp;
    double d2zndnp2;
    double d2zpdnn2;
    double d2zpdnndnp;
    double d2zpdnp2;
    //@}

    /// \name Main functions
    //@{
    /** \brief Solve for the fugacity given the density
     */
    virtual void solve_fugacity(double nn, double np,
				double lam_n, double lam_p,
				double b_n, double b_pn,
				double &zn, double &zp) {

      double npt=pow(lam_n,3)/2.0*np;
      double nnt=pow(lam_p,3)/2.0*nn;

      // At high densities or very low densities, just use the
      // non-interacting result
      //
      // AWS: 9/13/2020: I added the "|| nnt>1.0e5 || npt>1.0e5" option
      // later, and this might not be the best method
      // 
      if (nnt<5.0e-6 || npt<5.0e-6 || nnt>1.0e5 || npt>1.0e5) {
	zn=nnt;
	zp=npt;
	return;
      }

      zn=0.0;
      zp=0.0;
    
      double a=pow(b_n,3)*2.0/b_pn/b_pn-2.0*b_n;
      double b=-1+b_n*b_n*2.0/b_pn/b_pn-b_n/b_pn;
      double c=b_n/(2.0*b_pn*b_pn)-0.5/b_pn-nnt+npt-b_n*b_n*npt*2.0/b_pn/b_pn;
      double d=-b_n*npt/b_pn/b_pn+npt/2.0/b_pn;
      double e=b_n*npt*npt/2.0/b_pn/b_pn;
    
      quart.solve_rc(a,b,c,d,e,res[0],res[1],res[2],res[3]);
    
      std::vector<double> zp_list, zn_list;
      for(size_t k=0;k<4;k++) {
	if (res[k].imag()==0.0 && res[k].real()>0.0 && res[k].real()<500.0) {
	  double r0, r1;
	  gsl_poly_solve_quadratic(2.0*b_n,2.0*res[k].real()*
				   b_pn+1.0,-nnt,&r0,&r1);
	  if (r0>0.0 && r0<500.0) {
	    if (r1>0.0 && r1<500.0) {
	      O2SCL_ERR2("Unexpected pair of roots in ",
			 "virial_solver_deriv::solve_fugacity().",
			 o2scl::exc_einval);
	    }
	    std::cout << res[k] << "," << r0 << " ";
	    zp_list.push_back(res[k].real());
	    zn_list.push_back(r0);
	  }
	  if (r1>0.0 && r1<500.0) {
	    zp_list.push_back(res[k].real());
	    zn_list.push_back(r1);
	  }
	}
      }
      if (zp_list.size()==1) {
	zp=zp_list[0];
	zn=zn_list[0];
      } else if (zp_list.size()==2) {
	double norm_0=zp_list[0]*zp_list[0]+zn_list[0]*zn_list[0];
	double norm_1=zp_list[1]*zp_list[1]+zn_list[1]*zn_list[1];
	if (norm_0<norm_1) {
	  zp=zp_list[0];
	  zn=zn_list[0];
	} else {
	  zp=zp_list[1];
	  zn=zn_list[1];
	}
      } else {
	std::cout << "virial_solver_deriv::solve_fugacity "
		  << "multiplicity problem:\n\t"
		  << "res0,res1: " << res[0] << " " << res[1] << "\n\t"
		  << "res2,res3: " << res[2] << " " << res[3] << std::endl;
	std::cout << "\tnn,np,lam_n,lam_p: " << nn << " " << np << " "
		  << lam_n << " " << lam_p << "\n\t"
		  << "nnt,npt,zp_list.size(): " << nnt << " " << npt << " "
		  << zp_list.size() << std::endl;
	O2SCL_ERR2("Unexpected root multiplicity in ",
		   "virial_solver_deriv::solve_fugacity().",o2scl::exc_einval);
      }
      return;
    }

    /** \brief Compute the derivatives given the densities and the
	fugacities (i.e. after a call to \ref solve_fugacity())
    */
    virtual void calc_deriv(double nn, double np,
			    double lam_n, double lam_p,
			    double b_n, double b_pn,
			    double zn, double zp,
			    double dbndT, double dbpndT,
			    double dlamndT, double dlampdT) {
    
      double npt=pow(lam_n,3)/2.0*np;
      double nnt=pow(lam_p,3)/2.0*nn;
    
      // At high densities or very low densities, just use the
      // non-interacting result

      if (nnt<5.0e-6 || npt<5.0e-6) {
      
	dzndnn=nnt/nn;
	dzpdnp=npt/np;
	dzndnp=0.0;
	dzpdnn=0.0;
      
	d2zndnn2=0.0;
	d2zndnndnp=0.0;
	d2zndnp2=0.0;
	d2zpdnn2=0.0;
	d2zpdnndnp=0.0;
	d2zpdnp2=0.0;
	dzndT=1.5*lam_n*lam_n*nn*dlamndT;
	dzpdT=1.5*lam_p*lam_p*np*dlampdT;
      
	return;
      }

      // Code automatically generated by scipy:
      
      dzndnn=-pow(lam_n,3)*(2*b_n*zp+b_pn*zn+1.0/2.0)/
	      (4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*
	       (4*b_n*zp+2*b_pn*zn+1));
      dzndnp=b_pn*pow(lam_p,3)*zn/
	      (4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*
	       (4*b_n*zp+2*b_pn*zn+1));
      dzpdnn=b_pn*pow(lam_n,3)*zp/
	      (4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*
	       (4*b_n*zp+2*b_pn*zn+1));
      dzpdnp=-pow(lam_p,3)*(2*b_n*zn+b_pn*zp+1.0/2.0)/
	      (4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*
	       (4*b_n*zp+2*b_pn*zn+1));
    
      double dzndnn_dzn=-b_pn*pow(lam_n,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1))-pow(lam_n,3)*(2*b_n*zp+b_pn*zn+1.0/2.0)*(-4*b_n*(-4*b_n*zp-2*b_pn*zn-1)-4*pow(b_pn,2)*zp+2*b_pn*(4*b_n*zn+2*b_pn*zp+1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzndnn_dzp=-2*b_n*pow(lam_n,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1))-pow(lam_n,3)*(2*b_n*zp+b_pn*zn+1.0/2.0)*(4*b_n*(4*b_n*zn+2*b_pn*zp+1)-4*pow(b_pn,2)*zn-2*b_pn*(-4*b_n*zp-2*b_pn*zn-1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzndnp_dzn=b_pn*pow(lam_p,3)*zn*(-4*b_n*(-4*b_n*zp-2*b_pn*zn-1)-4*pow(b_pn,2)*zp+2*b_pn*(4*b_n*zn+2*b_pn*zp+1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2)+b_pn*pow(lam_p,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1));
      double dzndnp_dzp=b_pn*pow(lam_p,3)*zn*(4*b_n*(4*b_n*zn+2*b_pn*zp+1)-4*pow(b_pn,2)*zn-2*b_pn*(-4*b_n*zp-2*b_pn*zn-1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzpdnn_dzn=b_pn*pow(lam_n,3)*zp*(-4*b_n*(-4*b_n*zp-2*b_pn*zn-1)-4*pow(b_pn,2)*zp+2*b_pn*(4*b_n*zn+2*b_pn*zp+1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzpdnn_dzp=b_pn*pow(lam_n,3)*zp*(4*b_n*(4*b_n*zn+2*b_pn*zp+1)-4*pow(b_pn,2)*zn-2*b_pn*(-4*b_n*zp-2*b_pn*zn-1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2)+b_pn*pow(lam_n,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1));
      double dzpdnp_dzn=-2*b_n*pow(lam_p,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1))-pow(lam_p,3)*(2*b_n*zn+b_pn*zp+1.0/2.0)*(-4*b_n*(-4*b_n*zp-2*b_pn*zn-1)-4*pow(b_pn,2)*zp+2*b_pn*(4*b_n*zn+2*b_pn*zp+1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);
      double dzpdnp_dzp=-b_pn*pow(lam_p,3)/(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1))-pow(lam_p,3)*(2*b_n*zn+b_pn*zp+1.0/2.0)*(4*b_n*(4*b_n*zn+2*b_pn*zp+1)-4*pow(b_pn,2)*zn-2*b_pn*(-4*b_n*zp-2*b_pn*zn-1))/pow(4*pow(b_pn,2)*zn*zp-(4*b_n*zn+2*b_pn*zp+1)*(4*b_n*zp+2*b_pn*zn+1),2);

      d2zndnn2=dzndnn_dzn*dzndnn+dzndnn_dzp*dzpdnn;
      d2zndnndnp=dzndnn_dzn*dzndnp+dzndnn_dzp*dzpdnp;
      d2zndnp2=dzndnp_dzn*dzndnp+dzndnp_dzp*dzpdnp;
      d2zpdnn2=dzpdnn_dzn*dzndnn+dzpdnn_dzp*dzpdnn;
      d2zpdnndnp=dzpdnn_dzn*dzndnp+dzpdnn_dzp*dzpdnp;
      d2zpdnp2=dzpdnp_dzn*dzndnp+dzpdnp_dzp*dzpdnp;

      dzndT=(1.0/2.0)*(-16*b_n*dbndT*pow(zn,2)*zp-16*b_n*dbpndT*zn*pow(zp,2)+12*b_n*dlamndT*pow(lam_n,2)*nn*zp-8*b_pn*dbndT*pow(zn,3)+8*b_pn*dbndT*zn*pow(zp,2)+6*b_pn*dlamndT*pow(lam_n,2)*nn*zn-6*b_pn*dlampdT*pow(lam_p,2)*np*zn-4*dbndT*pow(zn,2)-4*dbpndT*zn*zp+3*dlamndT*pow(lam_n,2)*nn)/(16*pow(b_n,2)*zn*zp+8*b_n*b_pn*pow(zn,2)+8*b_n*b_pn*pow(zp,2)+4*b_n*zn+4*b_n*zp+2*b_pn*zn+2*b_pn*zp+1);
    
      dzpdT=(-b_pn*dzndT*zp+(3.0/4.0)*dlamndT*pow(lam_n,2)*nn-1.0/2.0*dzndT-zn*(2*b_n*dzndT+dbndT*zn+dbpndT*zp))/(b_pn*zn);
    
      return;
    }
    //@}

  };

}

#endif
