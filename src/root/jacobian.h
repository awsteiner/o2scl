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
  
  /// Jacobian function (not necessarily square) in src/root/jacobian.h
  typedef std::function<
    int(size_t,boost::numeric::ublas::vector<double> &,
	size_t,boost::numeric::ublas::vector<double> &,
	boost::numeric::ublas::matrix<double> &) > jac_funct;

  /** \brief Base for providing a numerical jacobian [abstract base]
      
      This is provides a Jacobian which is numerically determined
      by differentiating a user-specified function (typically 
      of the form of \ref mm_funct). 

      By convention, the Jacobian is stored in the order
      <tt>J[i][j]</tt> (or <tt>J(i,j)</tt>) where the rows have index
      \c i which runs from 0 to <tt>ny-1</tt> and the columns have
      index \c j with runs from 0 to <tt>nx-1</tt>.

      Default template arguments
      - \c func_t - \ref mm_funct
      - \c vec_t - boost::numeric::ublas::vector<double>
      - \c mat_t - boost::numeric::ublas::matrix<double>
  */
  template<class func_t=mm_funct, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class jacobian {
    
  public:
    
  jacobian() {
    err_nonconv=true;
  };
    
  virtual ~jacobian() {};
    
  /// If true, call the error handler if the routine does not converge
  bool err_nonconv;

  /// Set the function to compute the Jacobian of
  virtual int set_function(func_t &f) {
    func=f;
    return 0;
  }
    
  /** \brief Evaluate the Jacobian \c j at point \c y(x)
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, 
			 mat_t &j)=0;
    
#ifndef DOXYGEN_INTERNAL
    
  protected:
    
  /// A pointer to the user-specified function
  func_t func;
    
  private:
    
  jacobian(const jacobian &);
  jacobian& operator=(const jacobian&);
    
#endif
    
  };
  
  /** \brief Simple automatic Jacobian

      This class computes a numerical Jacobian by finite differencing.
      The stepsize is initially chosen to be
      \f$ h_j = \mathrm{max}(\mathrm{epsrel}~|x_j|,\mathrm{epsmin}) \f$.
      Then if \f$ h_j = 0 \f$, the value of \f$ h_j \f$ is set to
      \f$ \mathrm{epsrel}) \f$ .

      Values of \c epsmin which are non-zero are useful, for example,
      in \ref mroot_hybrids when one of the variables is either
      very small or zero, so that the step size doesn't become too
      small. 

      If the function evaluation leads to a non-zero return value,
      then the step size is alternately flipped in sign or decreased
      by a fixed factor (default \f$ 10^2 \f$, set in \ref
      set_shrink_fact() ) in order to obtain a valid result. This
      process is repeated a fixed number of times (default 10, set in
      \ref set_max_shrink_iters() ).
      
      This is equivalent to the GSL method for computing Jacobians as
      in \c multiroots/fdjac.c if one calls \ref
      set_max_shrink_iters() with a parameter value of zero.

      If one row of the Jacobian is all zero, or if there was no
      step-size found which would give a zero return value from
      the user-specified function, then the error handler is called
      depending on the value of \ref err_nonconv.
      
      This class does not separately check the vector and matrix sizes
      to ensure they are commensurate. 

      Default template arguments
      - \c func_t - \ref mm_funct
      - \c vec_t - boost::numeric::ublas::vector<double>
      - \c mat_t - boost::numeric::ublas::matrix<double>
  */
  template<class func_t=mm_funct, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > 
    class jacobian_gsl : public jacobian<func_t,vec_t,mat_t> {
    
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

  /** \brief The relative stepsize for finite-differencing 
  */
  double epsrel;
    
  /// The minimum stepsize 
  double epsmin;

  /// Maximum number of times to shrink the step size
  size_t max_shrink_iters;

  /// Factor to shrink stepsize by
  double shrink_fact;

#endif

  public:
    
  jacobian_gsl() {
    epsrel=sqrt(std::numeric_limits<double>::epsilon());
    epsmin=0.0;
    mem_size_x=0;
    mem_size_y=0;
    max_shrink_iters=10;
    shrink_fact=1.0e2;
  }

  virtual ~jacobian_gsl() {
  }

  /** \brief Get the relative stepsize (default \f$ 10^{-4} \f$ )
   */
  double get_epsrel() { return epsrel; }

  /** \brief Get the minimum stepsize (default \f$ 10^{-15} \f$)
   */
  double get_epsmin() { return epsmin; }

  /** \brief Set the relative stepsize (must be \f$ > 0 \f$)
   */
  void set_epsrel(double l_epsrel) {
    if (l_epsrel<=0.0) {
      O2SCL_ERR2("Negative or zero value specified in ",
		 "jacobian_gsl::set_epsrel().",exc_einval);
    }
    epsrel=l_epsrel;
    return;
  }

  /** \brief Set the minimum stepsize (must be \f$ \geq 0 \f$)
   */
  void set_epsmin(double l_epsmin) {
    if (l_epsmin<0.0) {
      O2SCL_ERR("Negative value specified in jacobian_gsl::set_epsmin().",
		exc_einval);
    }
    epsmin=l_epsmin;
    return;
  }
  
  /** \brief Set shrink factor for decreasing step size
   */
  void set_shrink_fact(double l_shrink_fact) {
    if (l_shrink_fact<0.0) {
      O2SCL_ERR("Negative value specified in jacobian_gsl::set_shrink_fact().",
		exc_einval);
    }
    shrink_fact=l_shrink_fact;
    return;
  }

  /** \brief Set number of times to decrease step size
   */
  void set_max_shrink_iters(size_t it) {
    max_shrink_iters=it;
    return;
  }
  
  /** \brief The operator()
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, 
			 mat_t &jac) {
      
    size_t i,j;
    double h,temp;

    if (mem_size_x!=nx || mem_size_y!=ny) {
      f.resize(ny);
      xx.resize(nx);
      mem_size_x=nx;
      mem_size_y=ny;
    }
      
    vector_copy(nx,x,xx);

    for (j=0;j<nx;j++) {

      // Thanks to suggestion from Conrad Curry.
      h=epsrel*fabs(x[j]);
      if (h<epsmin) h=epsmin;
      if (h==0.0) h=epsrel;

      xx[j]=x[j]+h;
      int ret=(this->func)(nx,xx,f);
      xx[j]=x[j];
      
      // The function returned a non-zero value, so try a different step
      size_t it=0;
      while (ret!=0 && h>=epsmin && it<max_shrink_iters) {

	// First try flipping the sign
	h=-h;
	xx[j]=x[j]+h;
	ret=(this->func)(nx,xx,f);
	xx[j]=x[j];

	if (ret!=0) {

	  // If that didn't work, flip to positive and try a smaller
	  // stepsize
	  h/=-shrink_fact;
	  if (h>=epsmin) {
	    xx[j]=x[j]+h;
	    ret=(this->func)(nx,xx,f);
	    xx[j]=x[j];
	  }
	  
	}

	it++;
      }

      if (ret!=0) {
	O2SCL_CONV2_RET("Jacobian failed to find valid step in ",
			"jacobian_gsl::operator().",exc_ebadfunc,
			this->err_nonconv);
      }

      // This is the equivalent of GSL's test of
      // gsl_vector_isnull(&col.vector)

      bool nonzero=false;
      for (i=0;i<ny;i++) {
	temp=(f[i]-y[i])/h;
	if (temp!=0.0) nonzero=true;
	jac(i,j)=temp;
      }
      if (nonzero==false) {
	O2SCL_CONV_RET((((std::string)"Row ")+o2scl::szttos(j)+
			" of the Jacobian is zero "+
			"in jacobian_gsl::operator().").c_str(),exc_esing,
		       this->err_nonconv);
      }

    }
    
    return 0;
  }

  };
  
  /** \brief A direct calculation of the jacobian using a \ref
      deriv_base object
      
      By default, the stepsize, \ref deriv_gsl::h is set to \f$
      10^{-4} \f$ in the \ref jacobian_exact constructor.

      Note that it is most often wasteful to use this Jacobian in a
      root-finding routine and using more approximate Jacobians is
      more efficient. 

      Default template arguments
      - \c func_t - \ref mm_funct
      - \c vec_t - boost::numeric::ublas::vector<double>
      - \c mat_t - boost::numeric::ublas::matrix<double>
  */
  template<class func_t=mm_funct, 
    class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class jacobian_exact : 
  public jacobian<func_t,vec_t,mat_t> {
    
  public:
    
  jacobian_exact() {
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
  deriv_gsl<> def_deriv;
  
  /// Set the derivative object
  int set_deriv(deriv_base<> &de) {
    dptr=&de;
    return 0;
  }
    
  /** \brief The operator()
   */
  virtual int operator()(size_t nx, vec_t &x, size_t ny, vec_t &y, 
			 mat_t &jac) {

    // The function return value
    int ret=0;

    // Temporary storage fo the derivative uncertainty
    double err;

    ej_parms ejp;
    ejp.nx=nx;
    ejp.ny=ny;
    ejp.x=&x;
    ejp.y=&y;
    
    funct dfnp=std::bind(std::mem_fn<double(double,ej_parms &)>
			   (&jacobian_exact::dfn),
			   this,std::placeholders::_1,std::ref(ejp));

    for (size_t j=0;j<nx;j++) {
      ejp.xj=j;
      for (size_t i=0;i<ny;i++) {
	ejp.yi=i;
	double tmp=(*ejp.x)[j];
	int dret=dptr->deriv_err(tmp,dfnp,jac(i,j),err);
	if (dret!=0) {
	  if (this->err_nonconv==true) {
	    O2SCL_ERR2("Derivative object tailed in jacobian_exact::",
		      "operator().",o2scl::exc_efailed);
	  }
	  ret=1;
	}
	(*ejp.x)[j]=tmp;
      }
    }
    
    return ret;
  }

  /** \brief Compute the Jacobian and its uncertainty
      from the numerical differentiation
   */
  virtual int jac_err(size_t nx, vec_t &x, size_t ny, vec_t &y, 
		      mat_t &jac, mat_t &err) {

    // The function return value
    int ret=0;

    ej_parms ejp;
    ejp.nx=nx;
    ejp.ny=ny;
    ejp.x=&x;
    ejp.y=&y;
    
    funct dfnp=std::bind(std::mem_fn<double(double,ej_parms &)>
			   (&jacobian_exact::dfn),
			   this,std::placeholders::_1,std::ref(ejp));

    for (size_t j=0;j<nx;j++) {
      ejp.xj=j;
      for (size_t i=0;i<ny;i++) {
	ejp.yi=i;
	double tmp=(*ejp.x)[j];
	int dret=dptr->deriv_err(tmp,dfnp,jac(i,j),err(i,j));
	if (dret!=0) {
	  if (this->err_nonconv==true) {
	    O2SCL_ERR2("Derivative object tailed in jacobian_exact::",
		       "jac_err().",o2scl::exc_efailed);
	  }
	  ret=1;
	}
	(*ejp.x)[j]=tmp;
      }
    }
    
    return ret;
  }
  
#ifndef DOXYGEN_INTERNAL

  protected:

  /// Pointer to the derivative object
  deriv_base<> *dptr;
    
  /// Function for the derivative object
  double dfn(double x, ej_parms &ejp) {
    (*ejp.x)[ejp.xj]=x;
    (this->func)(ejp.nx,*ejp.x,*ejp.y);
    return (*ejp.y)[ejp.yi];
  }

#endif

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
