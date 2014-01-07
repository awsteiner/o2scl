/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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
#ifndef O2SCL_ODE_JAC_FUNCT_H
#define O2SCL_ODE_JAC_FUNCT_H

/** \file ode_jac_funct.h
    \brief File defining ODE Jacobian function objects
*/

#include <string>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

#if defined (O2SCL_CPP11) || defined (DOXYGENP)
  
  /// Array of multi-dimensional functions for C++11
  typedef std::function<
    int(double,size_t,const boost::numeric::ublas::vector<double> &,
	boost::numeric::ublas::matrix<double> &,
	boost::numeric::ublas::vector<double> &) > ode_jac_funct11;
  
#endif

  /** \brief Ordinary differential equation Jacobian [abstract base]
      
      This base class provides the basic format for specifying the
      Jacobian for ordinary differential equations to integrate with
      the \o2 ODE solvers. Select the appropriate child of this class
      according to the kind of functions which are to be given to the
      solver.
  */
  template <class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > class ode_jac_funct {

    public:

    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;

    ode_jac_funct() {}

    virtual ~ode_jac_funct() {}

    /** \brief Compute the derivatives \c dfdx and the Jacobian matrix
	\c dfdy given \c y at the point \c x.
    */
    virtual int operator()(double x, size_t nv, const vec_t &y, 
			   mat_t &dfdy, vec_t &dfdx)=0;
  
#ifndef DOXYGEN_INTERNAL

    private:

    ode_jac_funct(const ode_jac_funct &);
    ode_jac_funct& operator=(const ode_jac_funct&);

#endif

  };

  /** \brief Provide ODE Jacobian in the form of function pointers
  */
    template <class vec_t=boost::numeric::ublas::vector<double> , 
      class mat_t=boost::numeric::ublas::matrix<double> >
    class ode_jac_funct_fptr : public ode_jac_funct<vec_t> 
    {

    public:
    
    /** \brief Create an object given a function pointer
     */
    ode_jac_funct_fptr(int (*fp)(double, size_t, const vec_t &, 
				 mat_t &, vec_t &)) {
      fptr=fp;
    }
    
    virtual ~ode_jac_funct_fptr() {}
    
    /** \brief Compute the derivatives \c dfdx and the Jacobian matrix
	\c dfdy given \c y at the point \c x.
    */
    virtual int operator()(double x, size_t nv, const vec_t &y, 
			   mat_t &dfdy, vec_t &dfdx) {
      return fptr(x,nv,y,dfdy,dfdx);
    }
    
#ifndef DOXYGEN_INTERNAL
    
    protected:

    ode_jac_funct_fptr() {};

    /// The function pointer
    int (*fptr)(double x, size_t nv, const vec_t &y, 
		mat_t &dfdy, vec_t &dfdx);

    private:

    ode_jac_funct_fptr(const ode_jac_funct_fptr &);
    ode_jac_funct_fptr& operator=(const ode_jac_funct_fptr&);

#endif

  };

  /** \brief Provide ODE Jacobian in the form of member 
      function pointers
  */
  template <class tclass, class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > 
    class ode_jac_funct_mfptr : public ode_jac_funct<vec_t> {
    
    public:
    
    /** \brief Create an object given a class and member function pointer
     */
    ode_jac_funct_mfptr
    (tclass *tp, int (tclass::*fp)(double x, size_t nv, const vec_t &y, 
				   mat_t &dfdy, vec_t &dfdx)) {
      tptr=tp;
      fptr=fp;
    }
    
    virtual ~ode_jac_funct_mfptr() {};
    
    /** \brief Compute the derivatives \c dfdx and the Jacobian matrix
	\c dfdy given \c y at the point \c x.
    */
    virtual int operator()(double x, size_t nv, const vec_t &y, 
			   mat_t &dfdy, vec_t &dfdx) {
      return (*tptr.*fptr)(x,nv,y,dfdy,dfdx);
    }
  
#ifndef DOXYGEN_INTERNAL

    protected:

    /// The pointer to the member function
    int (tclass::*fptr)(double x, size_t nv, const vec_t &y, 
			mat_t &dfdy, vec_t &dfdx);

    /// The pointer to the class
    tclass *tptr;

    private:

    ode_jac_funct_mfptr(const ode_jac_funct_mfptr &);
    ode_jac_funct_mfptr& operator=(const ode_jac_funct_mfptr&);

#endif

  };

  /** \brief Provide ODE Jacobian in the form of const member 
      function pointers
  */
  template <class tclass, class vec_t=boost::numeric::ublas::vector<double>, 
    class mat_t=boost::numeric::ublas::matrix<double> > 
    class ode_jac_funct_cmfptr : public ode_jac_funct<vec_t> {

    public:
    
    /** \brief Create an object given a class and member function pointer
     */
    ode_jac_funct_cmfptr
    (tclass *tp, int (tclass::*fp)(double x, size_t nv, const vec_t &y, 
				   mat_t &dfdy, vec_t &dfdx) const) {
      tptr=tp;
      fptr=fp;
    }
  
    virtual ~ode_jac_funct_cmfptr() {};
  
    /** \brief Compute the derivatives \c dfdx and the Jacobian matrix
	\c dfdy given \c y at the point \c x.
    */
    virtual int operator()(double x, size_t nv, const vec_t &y, 
			   mat_t &dfdy, vec_t &dfdx) {
      return (*tptr.*fptr)(x,nv,y,dfdy,dfdx);
    }
  
#ifndef DOXYGEN_INTERNAL

    protected:

    /// The pointer to the member function
    int (tclass::*fptr)(double x, size_t nv, const vec_t &y, 
			mat_t &dfdy, vec_t &dfdx) const;

    /// The pointer to the class
    tclass *tptr;

    private:

    ode_jac_funct_cmfptr(const ode_jac_funct_cmfptr &);
    ode_jac_funct_cmfptr& operator=(const ode_jac_funct_cmfptr&);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
