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
#ifndef O2SCL_FUNCT_H
#define O2SCL_FUNCT_H

/** \file funct.h
    \brief Function object classes for one-dimensional functions
*/

#include <string>
#include <functional>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/fparser.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// One-dimensional function typedef
  typedef std::function<double(double)> funct11;

#ifdef O2SCL_NEVER_DEFINED
  
  /** \brief One-dimensional function [abstract base]
      
      This class generalizes a one-dimensional function y(x).

      This class is one of a large number of function object classes
      in \o2 designed to provide a mechanism for the user to 
      supply functions to solvers, minimizers, integrators, etc.
      See \ref funct_section for a general description.

      \todo I would like to make the operator() const, but this causes
      problems with the function parser class in funct_string.
  */
  class funct {
    
  public:  
    
    /*
      This (undocumented) empty constructor is required for the 
      child funct classes
    */
    funct() {}

    virtual ~funct() {}
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const=0;

  };
  
  /** \brief Function pointer to a function
   */
  class funct_fptr : public funct {

  public:
    
    /** \brief Specify the function pointer
     */
    funct_fptr(double (*fp)(double)) {
      fptr=fp;
    }
    
    virtual ~funct_fptr() {}

    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const {
      return fptr(x);
    }

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Function pointer
    double (*fptr)(double x);
    
    /// Copy constructor
  funct_fptr(const funct_fptr &f) : funct() {
      fptr=f.fptr;
    }
    
    /// Copy constructor
    funct_fptr &operator=(const funct_fptr &f) {
      fptr=f.fptr;
      return *this;
    }

#endif

  };

  /** \brief Function pointer to a function with a parameter
   */
  template<class param_t> class funct_fptr_param : public funct {
    
  public:
    
    /** \brief Specify the function pointer
     */
    funct_fptr_param(double (*fp)(double, param_t &), param_t &pa) {
      fptr=fp;
      pp=&pa;
    }
    
    virtual ~funct_fptr_param() {}

    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const {
      return fptr(x,*pp);
    }
    
#ifndef DOXYGEN_INTERNAL

  protected:

    /// The function pointer
    double (*fptr)(double x, param_t &pa);

    /// The parameter
    param_t *pp;

    /// Copy constructor
    funct_fptr_param(const funct_fptr_param &f) {
      fptr=f.fptr;
    }
    
    /// Copy constructor
    funct_fptr_param &operator=(const funct_fptr_param &f) {
      fptr=f.fptr;
      return *this;
    }

#endif

  };

  /** \brief Member function pointer to a one-dimensional function
   */
  template <class tclass> class funct_mfptr : public funct {

  public:
  
    /** \brief Specify the member function pointer
     */
    funct_mfptr(tclass *tp, double (tclass::*fp)(double x)) {
      tptr=tp;
      fptr=fp;
    }
  
    virtual ~funct_mfptr() {};
  
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const {
      return (*tptr.*fptr)(x);
    }

#ifndef DOXYGEN_INTERNAL

  protected:
  
    /// Storage for the member function pointer
    double (tclass::*fptr)(double x);

    /// Store the pointer to the class instance
    tclass *tptr;

    /// Copy constructor
    funct_mfptr(const funct_mfptr &f) {
      fptr=f.fptr;
      tptr=f.tptr;
    }
    
    /// Copy constructor
    funct_mfptr &operator=(const funct_mfptr &f) {
      fptr=f.fptr;
      tptr=f.tptr;
      return *this;
    }

#endif

  };

  /** \brief Member function pointer to a one-dimensional function
      with a parameter
  */
  template <class tclass, class param_t> class funct_mfptr_param : 
  public funct {

  public:
  
    /** \brief Specify the member function pointer
     */
    funct_mfptr_param(tclass *tp, 
		      double (tclass::*fp)(double x, param_t &pa), 
		      param_t &pa) {
      tptr=tp;
      fptr=fp;
      pp=&pa;
    }
  
    virtual ~funct_mfptr_param() {};
  
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const {
      return (*tptr.*fptr)(x,*pp);
    }

#ifndef DOXYGEN_INTERNAL

  protected:
  
    /// Storage for the member function pointer
    double (tclass::*fptr)(double x, param_t &pa);

    /// Store the pointer to the class instance
    tclass *tptr;

    /// The parameter
    param_t *pp;

    /// Copy constructor
    funct_mfptr_param(const funct_mfptr_param &f) {
      fptr=f.fptr;
      tptr=f.tptr;
      pp=f.pp;
    }
    
    /// Copy constructor
    funct_mfptr_param &operator=(const funct_mfptr_param &f) {
      fptr=f.fptr;
      tptr=f.tptr;
      pp=f.pp;
      return *this;
    }

#endif

  };

  /** \brief Const member function pointer to a one-dimensional function

      \note While this is designed to accept a pointer to a const
      member function, the choice of whether the class pointer
      given in the template type <tt>tclass</tt> is const or 
      not is up to the user.
  */
  template <class tclass> class funct_cmfptr : public funct {

  public:
  
    /** \brief Specify the member function pointer
     */
    funct_cmfptr(tclass *tp, double (tclass::*fp)(double x) const) {
      tptr=tp;
      fptr=fp;
    }
  
    virtual ~funct_cmfptr() {};
  
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const {
      return (*tptr.*fptr)(x);
    }

#ifndef DOXYGEN_INTERNAL

  protected:
  
    /// Storage for the const member function pointer
    double (tclass::*fptr)(double x) const;

    /// Store the pointer to the class instance
    tclass *tptr;

    /// Copy constructor
    funct_cmfptr(const funct_cmfptr &f) {
      fptr=f.fptr;
      tptr=f.tptr;
    }
    
    /// Copy constructor
    funct_cmfptr &operator=(const funct_cmfptr &f) {
      fptr=f.fptr;
      tptr=f.tptr;
      return *this;
    }

#endif

  };

  /** \brief Const member function pointer to a one-dimensional function 
      with a parameter

      \note While this is designed to accept a pointer to a const
      member function, the choice of whether the class pointer
      given in the template type <tt>tclass</tt> is const or 
      not is up to the user.
  */
  template <class tclass, class param_t> class funct_cmfptr_param : 
  public funct {

  public:
  
    /** \brief Specify the member function pointer
     */
    funct_cmfptr_param(tclass *tp, 
		       double (tclass::*fp)(double x, param_t &pa) const,
		       param_t &pa) {
      tptr=tp;
      fptr=fp;
      pp=&pa;
    }
  
    virtual ~funct_cmfptr_param() {};
  
    /** \brief Compute the function at point \c x with parameter \c pa 
	and return the result
     */
    virtual double operator()(double x, param_t &pa) const {
      return (*tptr.*fptr)(x,*pp);
    }

#ifndef DOXYGEN_INTERNAL

  protected:
  
    /// Storage for the const member function pointer
    double (tclass::*fptr)(double x, param_t &pa) const;

    /// Store the pointer to the class instance
    tclass *tptr;

    /// The parameter
    param_t *pp;

    /// Copy constructor
    funct_cmfptr_param(const funct_cmfptr_param &f) {
      fptr=f.fptr;
      tptr=f.tptr;
      pp=f.pp;
    }
    
    /// Copy constructor
    funct_cmfptr_param &operator=(const funct_cmfptr_param &f) {
      fptr=f.fptr;
      tptr=f.tptr;
      pp=f.pp;
      return *this;
    }

#endif

  };

  /** \brief One-dimensional function from a string
      
      For example,
      \code
      funct_string f("pi*r^2","r",1,"pi");
      ubvector par(1);
      par[0]=o2scl_const::pi;
      f.set_parms(par);
      for(double r=1.0;r<=2.0;r+=0.1) {
        cout << f(r) << endl;
      }
      \endcode
      will print out the area of circles having radii between 1 and 2.
  */
  class funct_string : public funct {
    
  public:
    
    /** \brief Specify the string and the parameters
     */
    funct_string(std::string formula, std::string var, 
		 int np=0, std::string parms="") {
      if(np<1) {
	fpw.Parse(formula,var);
	st_np=0;
	st_parms="";
      } else {
	std::string all=var+","+parms;
	fpw.Parse(formula,all);
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_form=formula;
      st_var=var;
    }

    virtual ~funct_string() {
      if (st_np>0) {
	delete[] arr;
      }
    };

  
    /** \brief Specify the string and the parameters
     */
    int set_function(std::string formula, std::string var, 
		     int np=0, std::string parms="") {
      if(np<1) {
	fpw.Parse(formula,var);
	st_np=0;
	st_parms="";
      } else {
	std::string all=var+","+parms;
	fpw.Parse(formula,all);
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_form=formula;
      st_var=var;
      return 0;
    }

    /** \brief Set the values of the auxilliary parameters that were
	specified in \c parms in the constructor
    */
    template<class vec_t> int set_parms(const vec_t &p) {
      for(int i=0;i<st_np;i++) {
	arr[i]=p[i];
      }
      return 0;
    }
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const {
      int i;
      double y;
      if(st_np<1) {
	y=fpw.Eval(&x);
      } else {
	double *all=new double[st_np+1];
	all[0]=x;
	for(i=1;i<=st_np;i++) all[i]=arr[i-1];
	y=fpw.Eval(all);
	delete[] all;
      }
      return y;
    }

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The object for evaluating strings
    mutable FunctionParser fpw;

    /// The number of parameters
    int st_np;

    /// Storage for the \ref fpw call.
    double *arr;

    /// The formula
    std::string st_form;
    /// The variables
    std::string st_var;
    /// The parameters
    std::string st_parms;

    funct_string() {};

#endif
#ifndef DOXYGEN_NO_O2NS

  private:

    funct_string(const funct_string &);
    funct_string& operator=(const funct_string&);

#endif

  };

#endif

  /** \brief One-dimensional function from a string
      
      For example,
      \code
      funct11_string f("pi*r^2","r",1,"pi");
      ubvector par(1);
      par[0]=o2scl_const::pi;
      f.set_parms(par);
      for(double r=1.0;r<=2.0;r+=0.1) {
        cout << f(x) << endl;
      }
      \endcode
      will print out the area of circles having radii between 1 and 2.
  */
  class funct11_string {
    
  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;

    /** \brief Specify the string and the parameters
     */
    funct11_string(std::string formula, std::string var, 
		 int np=0, std::string parms="") {
      if(np<1) {
	fpw.Parse(formula,var);
	st_np=0;
	st_parms="";
      } else {
	std::string all=var+","+parms;
	fpw.Parse(formula,all);
	st_np=np;
	st_parms=parms;
	arr.resize(np);
      }
      st_form=formula;
      st_var=var;
    }

    virtual ~funct11_string() {
    };

  
    /** \brief Specify the string and the parameters
     */
    int set_function(std::string formula, std::string var, 
		     int np=0, std::string parms="") {
      if(np<1) {
	fpw.Parse(formula,var);
	st_np=0;
	st_parms="";
      } else {
	std::string all=var+","+parms;
	fpw.Parse(formula,all);
	st_np=np;
	st_parms=parms;
	arr.resize(np);
      }
      st_form=formula;
      st_var=var;
      return 0;
    }

    /** \brief Set the values of the auxilliary parameters that were
	specified in \c parms in the constructor
    */
    template<class vec_t> int set_parms(const vec_t &p) {
      for(int i=0;i<st_np;i++) {
	arr[i]=p[i];
      }
      return 0;
    }
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const {
      int i;
      double y;
      if(st_np<1) {
	y=fpw.Eval(&x);
      } else {
	double *all=new double[st_np+1];
	all[0]=x;
	for(i=1;i<=st_np;i++) all[i]=arr[i-1];
	y=fpw.Eval(all);
	delete[] all;
      }
      return y;
    }

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The object for evaluating strings
    mutable FunctionParser fpw;

    /// The number of parameters
    size_t st_np;

    /// Storage for the \ref fpw call.
    ubvector arr;

    /// The formula
    std::string st_form;
    /// The variables
    std::string st_var;
    /// The parameters
    std::string st_parms;

    funct11_string() {};

#endif
#ifndef DOXYGEN_NO_O2NS

  private:

    funct11_string(const funct11_string &);
    funct11_string& operator=(const funct11_string&);

#endif

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
