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
#ifndef O2SCL_GEN_INTE_H
#define O2SCL_GEN_INTE_H

/** \file inte_gen.h
    \brief File defining \ref o2scl::inte_gen
*/

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Generalized multi-dimensional integration [abstract base]
 
      Perform the generalized multi-dimensional integral:
      \f[ 
      \int_{x_0=a_0}^{x_0=b_0} f(x_0) \int_{x_1=a_1(x_0)}^{x_1=b_1(x_0)} 
      f(x_0, x_1) ...
      \int_{x_{\mathrm{n}-1}=a_{\mathrm{n}-1}(x_0,x_1,..,x_{\mathrm{n}-2})}^
      {x_{\mathrm{n}-1}=b_{\mathrm{n}-1}(x_0,x_1,..,x_{\mathrm{n}-2})} 
      f(x_0,x_1,...,x_{\mathrm{n-1}})~d x_{\mathrm{n}-1}~...~d x_1~d x_0
      \f]

      The functions \f$ a_i \f$ and \f$ b_i \f$ are specified 
      in the arguments <tt>a</tt> and <tt>b</tt> to the function
      \ref ginteg() or \ref ginteg_err() .

      In order to allow the user to specify only three functions (for
      the integrand, the lower limits, and the upper limits) the first
      argument to the limit and integrand functions is used to
      distinguish among the limits for each separate integral. So
      first argument to <tt>a</tt> for \f$ a_0() \f$ is 0, and the
      first argument to <tt>a</tt> for \f$ a_1() \f$ is 1, etc., and 
      similarly for the upper limits specified in <tt>b</tt> and 
      the integrands specified in <tt>func</tt>.

      \future It might be interesting to construct a child class of
      \ref o2scl::inte_gen which automatically transforms variables to
      a hypercube and then applies a child of \ref o2scl::inte_multi
      to do the integration.
  */
  template<class func_t, class lfunc_t, class ufunc_t,
    class vec_t=boost::numeric::ublas::vector<double> > class inte_gen {

    public:
  
    inte_gen() {
      tol_rel=1.0e-8;
      verbose=0;
      interror=0.0;
      err_nonconv=true;
    }

    virtual ~inte_gen() {}
  
    /** \brief Verbosity
     */
    int verbose;
      
    /** \brief The maximum "uncertainty" in the value of the integral
     */
    double tol_rel;
  
    /// If true, call the error handler if the routine does not "converge"
    bool err_nonconv;

    /** \brief Integrate function \c func from
	\f$ x_i=a_i(x_i) \f$ to \f$ x_i=b_i(x_i) \f$ for
	\f$ 0<i<\mathrm{ndim}-1 \f$.
    */
    virtual double ginteg(func_t &func, size_t ndim, lfunc_t &a, ufunc_t &b)=0;
  
    /** \brief Integrate function \c func from
	\f$ x_i=a_i(x_i) \f$ to \f$ x_i=b_i(x_i) \f$ for
	\f$ 0<i<\mathrm{ndim}-1 \f$.
    */
    virtual int ginteg_err(func_t &func, size_t ndim, lfunc_t &a, 
			   ufunc_t &b, double &res, double &err) {
      res=ginteg(func,ndim,a,b);
      return 0;
    }
  
    /** \brief Return the error in the result from the last call to 
	ginteg() or ginteg_err()
      
	This will quietly return zero if no integrations have been
	performed.
    */
    double get_error() { return interror; }
  
    /// Return string denoting type ("inte_gen")
    const char *type() { return "inte_gen"; }
  
#ifndef DOXYGEN_INTERNAL
  
    protected:
  
    /// The uncertainty for the last integration computation
    double interror;
  
#endif
  
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif





  
