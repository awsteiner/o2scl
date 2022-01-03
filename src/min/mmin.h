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
#ifndef O2SCL_MULTI_MIN_H
#define O2SCL_MULTI_MIN_H

/** \file mmin.h
    \brief File defining \ref o2scl::mmin_base, \ref o2scl::grad_funct,
    \ref o2scl::gradient, and \ref o2scl::gradient_gsl
*/

#include <o2scl/multi_funct.h>
#include <o2scl/mm_funct.h>
#include <o2scl/string_conv.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// Array of multi-dimensional functions typedef in src/min/mmin.h
  typedef std::function<int(size_t,boost::numeric::ublas::vector<double> &,
			    boost::numeric::ublas::vector<double> &)> 
    grad_funct;
  
  /** \brief Class for automatically computing gradients [abstract base]

      Default template arguments
      - \c func_t - (no default)
      - \c vec_t - \ref boost::numeric::ublas::vector \< double \> 

      \future Consider making an exact_grad class for computing exact
      gradients.
  */
  template<class func_t, class vec_t=boost::numeric::ublas::vector<double> > 
    class gradient {
    
  public:

  // (Need to have empty default constructor since we
  // have private copy constructor)
  gradient() {}

  virtual ~gradient() {}

  /// Set the function to compute the gradient of
  virtual int set_function(func_t &f) {
    func=&f;
    return 0;
  }

  /** \brief Compute the gradient \c g at the point \c x
   */
  virtual int operator()(size_t nv, vec_t &x, vec_t &g)=0;

#ifndef DOXYGEN_INTERNAL

  protected:

  /// A pointer to the user-specified function
  func_t *func;
  
  private:
  
  gradient(const gradient &);
  gradient& operator=(const gradient&);
  
#endif

  };

  /** \brief Simple automatic computation of gradient by finite 
      differencing

      \comment
      \endcomment
  */
  template<class func_t, class vec_t> class gradient_gsl :
  public gradient<func_t,vec_t> {
    
  public:
    
    gradient_gsl() {
      epsrel=1.0e-6;
      epsmin=1.0e-15;
    }
    
    virtual ~gradient_gsl() {}

    /** \brief The relative stepsize for finite-differencing
        (default \f$ 10^{-6} \f$ )
    */
    double epsrel;

    /// The minimum stepsize (default \f$ 10^{-15} \f$)
    double epsmin;

    /** \brief Compute the gradient \c g at the point \c x
     */
    virtual int operator()(size_t nv, vec_t &x, vec_t &g) {
      double fv1, fv2, h;

      fv1=(*this->func)(nv,x);
      
      for(size_t i=0;i<nv;i++) {
	
	h=epsrel*fabs(x[i]);
	if (fabs(h)<=epsmin) h=epsrel;
	
	x[i]+=h;
	fv2=(*this->func)(nv,x);
	x[i]-=h;
	g[i]=(fv2-fv1)/h;
	
      }
      
      return 0;
    }

  };
    
  /** \brief Multidimensional minimization [abstract base]

      <b>The template parameters:</b>
      The template parameter \c func_t specifies the function to 
      min and should be a class containing a definition 
      \code
      func_t::operator()(size_t nv, const vec_t &x, double &f);
      \endcode
      where \c f is the value of the function at \c x ,
      where \c x is a array-like class defining \c operator[] of size \c nv.
      The parameter \c dfunc_t (if used) should provide the gradient with
      \code
      func_t::operator()(size_t nv, vec_t &x, vec_t &g);
      \endcode
      where \c g is the gradient of the function at \c x. 

      Verbose I/O is sent through \c std::cout and \c std::cin by
      default, but this can be modified using \ref
      set_verbose_stream(). Note that this function
      stores pointers to the user-specified output streams,
      and these pointers are not copied in child copy
      constructors. 
  */
  template<class func_t=multi_funct, class dfunc_t=func_t,
    class vec_t=boost::numeric::ublas::vector<double> > class mmin_base {

#ifndef DOXYGEN_INTERNAL

  protected:

  /// Stream for verbose output
  std::ostream *outs;
  
  /// Stream for verbose input
  std::istream *ins;
    
#endif

  public:
    
  mmin_base() {
    verbose=0;
    ntrial=100;
    tol_rel=1.0e-4;
    tol_abs=1.0e-4;
    last_ntrial=0;
    err_nonconv=true;
    outs=&std::cout;
    ins=&std::cin;
  }

  virtual ~mmin_base() {}
      
  /// Output control
  int verbose;
      
  /// Maximum number of iterations
  int ntrial;
      
  /// Function value tolerance
  double tol_rel;
      
  /// The independent variable tolerance
  double tol_abs;

  /// The number of iterations for in the most recent minimization
  int last_ntrial;
      
  /// If true, call the error handler if the routine does not "converge"
  bool err_nonconv;
      
  /** \brief Set streams for verbose I/O
      
      Note that this function stores pointers to the user-specified
      output streams, and these pointers are not copied in child copy
      constructors.
  */
  int set_verbose_stream(std::ostream &out, std::istream &in) {
    outs=&out;
    ins=&in;
    return 0;
  }
    
  /** \brief Calculate the minimum \c min of \c func w.r.t. the
      array \c x of size \c nvar.
  */
  virtual int mmin(size_t nvar, vec_t &x, double &fmin, 
		   func_t &func)=0;
      
  /** \brief Calculate the minimum \c min of \c func
      w.r.t. the array \c x of size \c nvar with gradient
      \c dfunc
  */
  virtual int mmin_de(size_t nvar, vec_t &x, double &fmin, 
		      func_t &func, dfunc_t &dfunc)
  {
    return mmin(nvar,x,fmin,func);
  }
      
  /** \brief Print out iteration information.
	  
      Depending on the value of the variable verbose, this prints out
      the iteration information. If verbose=0, then no information is
      printed, while if verbose>1, then after each iteration, the
      present values of x and y are output to std::cout along with the
      iteration number. If verbose>=2 then each iteration waits for a
      character.
  */
  template<class vec2_t> 
  int print_iter(size_t nv, vec2_t &x, double y, int iter,
		 double value, double limit, std::string comment) 
  {

    if (verbose<=0) return 0;
    
    int i;
    char ch;

    (*outs) << comment << " Iteration: " << iter << std::endl;
    {
      (*outs) << "x: " << std::endl;
      for(i=0;i<((int)nv);i++) (*outs) << x[i] << " ";
      (*outs) << std::endl;
    }
    (*outs) << "y: " << y << " Val: " << value << " Lim: " 
    << limit << std::endl;
    if (verbose>1) {
      (*outs) << "Press a key and type enter to continue. ";
      (*ins) >> ch;
    }
	
    return 0;
  }
      
  /// Return string denoting type ("mmin_base")
  const char *type() { return "mmin_base"; }

  /** \brief Copy constructor
   */
  mmin_base<func_t,dfunc_t,vec_t>
  (const mmin_base<func_t,dfunc_t,vec_t> &mb) {
    this->verbose=mb.verbose;
    this->ntrial=mb.ntrial;
    this->tol_rel=mb.tol_rel;
    this->tol_abs=mb.tol_abs;
    this->last_ntrial=mb.last_ntrial;
    this->err_nonconv=mb.err_nonconv;
  }
  
  /** \brief Copy constructor from operator=
   */
  mmin_base<func_t,dfunc_t,vec_t>& operator=
  (const mmin_base<func_t,dfunc_t,vec_t> &mb) {

    if (this != &mb) {
      this->verbose=mb.verbose;
      this->ntrial=mb.ntrial;
      this->tol_rel=mb.tol_rel;
      this->tol_abs=mb.tol_abs;
      this->last_ntrial=mb.last_ntrial;
      this->err_nonconv=mb.err_nonconv;
    }

    return *this;
  }
      
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

