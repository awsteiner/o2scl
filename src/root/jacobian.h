/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2013, Andrew W. Steiner
  
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
#ifndef O2SCL_JACOBIAN_H
#define O2SCL_JACOBIAN_H

/** \file jacobian.h
    \brief File for Jacobian evaluation and function classes
*/

#include <string>
#include <o2scl/mm_funct.h>
#include <o2scl/deriv_gsl.h>
#include <o2scl/columnify.h>
#include <o2scl/vector.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

#if defined (O2SCL_CPP11) || defined (DOXYGENP)
  
  /// Jacobian function (not necessarily square) (C++11 version)
  typedef std::function<
    int(size_t,boost::numeric::ublas::vector<double> &,
	size_t,boost::numeric::ublas::vector<double> &,
	boost::numeric::ublas::matrix<double> &) > jac_funct11;

#endif

  /** \brief Base for a Jacobian (not necessarily square) where J is
      computed at x given y=f(x) [abstract base]

      Compute
      \f[
      J_{ij} = \frac{\partial f_i}{\partial x_j}
      \f]

      The \c vec_t objects in operator() could have been written to be
      \c const, but they are not \c const so that they can be used as
      temporary workspace. They are typically restored to their
      original values before operator() exits. The Jacobian is stored
      in the order <tt>J[i][j]</tt> where the rows have index \c i
      which runs from 0 to <tt>ny-1</tt> and the columns have index \c
      j with runs from 0 to <tt>nx-1</tt>.

      Default template arguments
      - \c vec_t - boost::numeric::ublas::vector<double>
      - \c mat_t - boost::numeric::ublas::matrix<double>
  */
  template<class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > 
    class jac_funct {

  public:  

  jac_funct() {}

  virtual ~jac_funct() {}
    
  /** \brief The operator()
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, 
			 mat_t &j)=0;
    
#ifndef DOXYGEN_NO_O2NS

  private:

  jac_funct(const jac_funct &);
  jac_funct& operator=(const jac_funct&);

#endif

  };

  /** \brief Function pointer to jacobian

      Default template arguments
      - \c vec_t - boost::numeric::ublas::vector<double>
      - \c mat_t - boost::numeric::ublas::matrix<double>
   */
  template<class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> >
    class jac_funct_fptr : public jac_funct<vec_t,mat_t> {
    
  public:
    
  /** \brief Specify the function pointer
   */
  jac_funct_fptr(int (*fp)(size_t nx, vec_t &x, size_t ny, vec_t &y, 
			   mat_t &j)) {
    fptr=fp;
  }
    
  virtual ~jac_funct_fptr();
    
  /** \brief The operator()
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, mat_t &j) {
    return fptr(nx,x,ny,y,j);
  }

#ifndef DOXYGEN_NO_O2NS

  protected:
  
  jac_funct_fptr() {};
    
  /// Function pointer
  int (*fptr)(size_t nx, vec_t &x, size_t ny, vec_t &y, mat_t &j);
    
  private:

  jac_funct_fptr(const jac_funct_fptr &);
  jac_funct_fptr& operator=(const jac_funct_fptr&);

#endif

  };

  /** \brief Member function pointer to a Jacobian
   */
  template <class tclass, class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class jac_funct_mfptr : 
  public jac_funct<vec_t,mat_t> {

  public:
  
  /** \brief Specify the member function pointer
   */
  jac_funct_mfptr(tclass *tp, int (tclass::*fp)
		  (size_t nx, vec_t &x, size_t ny, vec_t &y, mat_t &j)) {
    tptr=tp;
    fptr=fp;
  }
  
  virtual ~jac_funct_mfptr() {};
  
  /** \brief The operator()
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, mat_t &j) {
    return (*tptr.*fptr)(nx,x,ny,y,j);
  }
  
#ifndef DOXYGEN_INTERNAL

  protected:
  
  /// Member function pointer
  int (tclass::*fptr)(size_t nx, vec_t &x, size_t ny, vec_t &y, mat_t &j);

  /// Class pointer
  tclass *tptr;

#endif

#ifndef DOXYGEN_NO_O2NS
  
  private:

  jac_funct_mfptr(const jac_funct_mfptr &);
  jac_funct_mfptr& operator=(const jac_funct_mfptr&);

#endif

  };

  /** \brief Const member function pointer to a Jacobian
   */
  template <class tclass, class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class jac_funct_cmfptr : 
  public jac_funct<vec_t,mat_t> {

  public:
  
  /** \brief Specify the member function pointer
   */
  jac_funct_cmfptr(tclass *tp, int (tclass::*fp)
		   (size_t nx, vec_t &x, size_t ny, vec_t &y,
		    mat_t &j) const) {
    tptr=tp;
    fptr=fp;
  }
  
  virtual ~jac_funct_cmfptr() {};
  
  /** \brief The operator()
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, mat_t &j) {
    return (*tptr.*fptr)(nx,x,ny,y,j);
  }
  
#ifndef DOXYGEN_INTERNAL

  protected:
  
  /// Member function pointer
  int (tclass::*fptr)(size_t nx, vec_t &x, size_t ny, vec_t &y, 
		      mat_t &j) const;

  /// Class pointer
  tclass *tptr;

#endif

#ifndef DOXYGEN_NO_O2NS
  
  private:

  jac_funct_cmfptr(const jac_funct_cmfptr &);
  jac_funct_cmfptr& operator=(const jac_funct_cmfptr&);

#endif

  };

  /** \brief Base for providing a numerical jacobian [abstract base]
      
      This is provides a Jacobian which is numerically determined
      by differentiating a user-specified function (typically 
      of the form of \ref mm_funct). 

      By convention, the Jacobian is stored in the order
      <tt>J[i][j]</tt> (or <tt>J(i,j)</tt>) where the rows have index
      \c i which runs from 0 to <tt>ny-1</tt> and the columns have
      index \c j with runs from 0 to <tt>nx-1</tt>.

      Default template arguments
      - \c func_t - \ref mm_funct\<\>
      - \c vec_t - boost::numeric::ublas::vector<double>
      - \c mat_t - boost::numeric::ublas::matrix<double>
  */
  template<class func_t=mm_funct<>, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class jacobian : 
  public jac_funct<vec_t,mat_t> {
    
  public:
    
  jacobian() {
  };
    
  virtual ~jacobian() {};
    
  /// Set the function to compute the Jacobian of
  virtual int set_function(func_t &f) {
    func=&f;
    return 0;
  }
    
  /** \brief Evaluate the Jacobian \c j at point \c y(x)
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, 
			 mat_t &j)=0;
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
  /// A pointer to the user-specified function
  func_t *func;
    
  private:
    
  jacobian(const jacobian &);
  jacobian& operator=(const jacobian&);
    
#endif
    
  };
  
  /** \brief Simple automatic Jacobian

      This class computes a numerical Jacobian by finite differencing.
      The stepsize is chosen to be \f$ h_j = \mathrm{epsrel}~x_j \f$ or
      \f$ h_j = \mathrm{epsmin} \f$ if \f$ \mathrm{epsrel}~x_j <
      \mathrm{epsmin} \f$.
      
      This is nearly equivalent to the GSL method for computing
      Jacobians as in \c multiroots/fdjac.c. To obtain the GSL
      behavior, set \ref epsrel to \c GSL_SQRT_DBL_EPSILON and set
      \ref epsmin to zero. The \ref mroot_hybrids class sets \ref
      epsrel to \c GSL_SQRT_DBL_EPSILON in its constructor,
      but does not set \ref epsmin to zero.
      
      This class does not separately check the vector and matrix sizes
      to ensure they are commensurate. 

      Default template arguments
      - \c func_t - \ref mm_funct\<\>
      - \c vec_t - boost::numeric::ublas::vector<double>
      - \c mat_t - boost::numeric::ublas::matrix<double>
  */
  template<class func_t=mm_funct<>, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > 
    class jacobian_gsl :
    public jacobian<func_t,vec_t,mat_t> {
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
  
  /// Function values
  vec_t f;

  /// Function arguments
  vec_t xx;

  /// Size of allocated memory in x
  size_t mem_size_x;

  /// Size of allocated memory in y
  size_t mem_size_y;

#endif

  public:
    
  jacobian_gsl() {
    epsrel=1.0e-4;
    epsmin=1.0e-15;
    err_nonconv=true;
    mem_size_x=0;
    mem_size_y=0;
  }

  virtual ~jacobian_gsl() {
  }
  
  /** \brief The relative stepsize for finite-differencing 
      (default \f$ 10^{-4} \f$ )
  */
  double epsrel;
    
  /// The minimum stepsize (default \f$ 10^{-15} \f$)
  double epsmin;

  /// If true, call the error handler if the routine does not "converge"
  bool err_nonconv;

  /** \brief The operator()
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, 
			 mat_t &jac) {
      
    size_t i,j;
    double h,temp;
    bool success=true;

    if (mem_size_x!=nx || mem_size_y!=ny) {
      f.resize(ny);
      xx.resize(nx);
      mem_size_x=nx;
      mem_size_y=ny;
    }
      
    vector_copy(nx,x,xx);

    for (j=0;j<nx;j++) {
	
      h=epsrel*fabs(x[j]);
      if (fabs(h)<=epsmin) h=epsrel;
      
      xx[j]=x[j]+h;
      (*this->func)(nx,xx,f);
      xx[j]=x[j];

      // This is the equivalent of GSL's test of
      // gsl_vector_isnull(&col.vector)

      bool nonzero=false;
      for (i=0;i<ny;i++) {
	temp=(f[i]-y[i])/h;
	if (temp!=0.0) nonzero=true;
	jac(i,j)=temp;
      }
      if (nonzero==false) success=false;


    }
      
    if (success==false) {
      O2SCL_CONV2_RET("At least one row of the Jacobian is zero ",
		      "in jacobian_gsl::operator().",exc_esing,
		      this->err_nonconv);
    }
    return 0;
  }
  };
  
  /** \brief A direct calculation of the jacobian using a \ref
      deriv_base object
      
      Note that it is most often wasteful to use this Jacobian in a
      root-finding routine and using more approximate Jacobians is
      more efficient. This class is mostly useful for demonstration
      and testing purposes.

      By default, the stepsize, \ref deriv_gsl::h is set to \f$
      10^{-4} \f$ in the \ref exact_jacobian constructor.

      Default template arguments
      - \c func_t - \ref mm_funct\<\>
      - \c vec_t - boost::numeric::ublas::vector<double>
      - \c mat_t - boost::numeric::ublas::matrix<double>
  */
  template<class func_t=mm_funct<>, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class exact_jacobian : 
  public jacobian<func_t,vec_t,mat_t> {
    
  public:
    
  exact_jacobian() {
    def_deriv.h=1.0e-4;
    dptr=&def_deriv;
  }
    
  /** \brief Parameter structure for passing information
	
      This class is primarily useful for specifying derivatives
      for using the jacobian::set_deriv() function.

      \comment
      This type needs to be publicly available so that the
      user can properly specify a base 1-dimensional derivative
      object. 
      \endcomment
  */
  typedef struct {
    /// The number of variables
    size_t nx;
    /// The number of variables
    size_t ny;
    /// The current x value
    size_t xj;
    /// The current y value
    size_t yi;
    /// The x vector
    vec_t *x;
    /// The y vector
    vec_t *y;
  } ej_parms;

  /// The default derivative object
  deriv_gsl<funct> def_deriv;
  
  /// Set the derivative object
  int set_deriv(deriv_base<funct> &de) {
    dptr=&de;
    return 0;
  }
    
  /** \brief The operator()
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, 
			 mat_t &jac) {

    double h,temp;

    ej_parms ejp;
    ejp.nx=nx;
    ejp.ny=ny;
    ejp.x=&x;
    ejp.y=&y;
    
    funct_mfptr_param<exact_jacobian,ej_parms> 
    dfnp(this,&exact_jacobian::dfn,ejp);
      
    for (size_t j=0;j<nx;j++) {
      ejp.xj=j;
      for (size_t i=0;i<ny;i++) {
	ejp.yi=i;
	double tmp=(*ejp.x)[j];
	jac(i,j)=dptr->deriv(tmp,dfnp);
	(*ejp.x)[j]=tmp;
      }
    }

    return 0;
  }

#ifndef DOXYGEN_INTERNAL

  protected:

  /// Pointer to the derivative object
  deriv_base<funct> *dptr;
    
  /// Function for the derivative object
  double dfn(double x, ej_parms &ejp) {
    (*ejp.x)[ejp.xj]=x;
    (*this->func)(ejp.nx,*ejp.x,*ejp.y);
    return (*ejp.y)[ejp.yi];
  }

#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
