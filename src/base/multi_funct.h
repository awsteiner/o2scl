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
#ifndef O2SCL_MULTI_FUNCT_H
#define O2SCL_MULTI_FUNCT_H

/** \file multi_funct.h
    \brief Function object classes for a multi-dimensional function
*/

#include <string>
#include <functional>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/err_hnd.h>
#include <o2scl/fparser.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// Multi-dimensional function typedef (C++11 version)
  typedef std::function<
    double(size_t,const boost::numeric::ublas::vector<double> &)>
    multi_funct11;
  
  /** \brief Multi-dimensional function [abstract base]

      This class generalizes one function of several variables,
      i.e. \f$ y(x_0,x_1,...,x_{nv-1}) \f$ where \c nv is the number
      of variables in the function y. 

      This class is one of a large number of function object classes
      in \o2 designed to provide a mechanism for the user to 
      supply functions to solvers, minimizers, integrators, etc.
      See \ref funct_section for a general description.
      
  */
  template<class vec_t=boost::numeric::ublas::vector<double> > 
    class multi_funct {
    
    public:  

    multi_funct() {}

    virtual ~multi_funct() {}

    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
    virtual double operator()(size_t nv, const vec_t &x)=0;

#ifndef DOXYGEN_NO_O2NS

    private:

    multi_funct(const multi_funct &);
    multi_funct& operator=(const multi_funct&);

#endif

  };

  /** \brief Function pointer to a multi-dimensional function
   */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class multi_funct_fptr : public multi_funct<vec_t> {
    
    public:

    /** \brief Specify the function pointer
     */
    multi_funct_fptr(double (*fp)(size_t nv, const vec_t &x)) {
      fptr=fp;
    }

    virtual ~multi_funct_fptr() {};

    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
    virtual double operator()(size_t nv, const vec_t &x) {
      return fptr(nv,x);
    }

#ifndef DOXYGEN_INTERNAL

    protected:
    
    /// Store the function pointer
    double (*fptr)(size_t nv, const vec_t &x);

    multi_funct_fptr() {}

#ifndef DOXYGEN_NO_O2NS
#endif

    private:

    multi_funct_fptr(const multi_funct_fptr &);
    multi_funct_fptr& operator=(const multi_funct_fptr&);

#endif

  };

  /** \brief Function pointer to a multi-dimensional function
   */
  template<class param_t, class vec_t=boost::numeric::ublas::vector<double> >
    class multi_funct_fptr_param : public multi_funct<vec_t> {
    
    public:

    /** \brief Specify the function pointer
     */
    multi_funct_fptr_param(double (*fp)(size_t nv, const vec_t &x,
					param_t &), param_t &pa) {
      fptr=fp;
      pp=&pa;
    }

    virtual ~multi_funct_fptr_param() {};

    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
    virtual double operator()(size_t nv, const vec_t &x) {
      return fptr(nv,x,*pp);
    }

#ifndef DOXYGEN_INTERNAL

    protected:
    
    /// Store the function pointer
    double (*fptr)(size_t nv, const vec_t &x, param_t &);

    /// Parameter
    param_t *pp;

    multi_funct_fptr_param() {}

#ifndef DOXYGEN_NO_O2NS
#endif

    private:

    multi_funct_fptr_param(const multi_funct_fptr_param &);
    multi_funct_fptr_param& operator=(const multi_funct_fptr_param&);

#endif

  };

  /** \brief Member function pointer to a multi-dimensional function
   */
  template<class tclass, class vec_t=boost::numeric::ublas::vector<double> >
    class multi_funct_mfptr : public multi_funct<vec_t> {
    public:
  
    /** \brief Specify the member function pointer
     */
    multi_funct_mfptr(tclass *tp, double (tclass::*fp)
		      (size_t nv, const vec_t &x)) {
      tptr=tp;
      fptr=fp;
    }
  
    virtual ~multi_funct_mfptr() {}
  
    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
    virtual double operator()(size_t nv, const vec_t &x) {
      return (*tptr.*fptr)(nv,x);
    }

#ifndef DOXYGEN_INTERNAL
  
    protected:
  
    /// Store the function pointer
    double (tclass::*fptr)(size_t nv, const vec_t &x);

    /// Store a pointer to the class instance
    tclass *tptr;
  
#ifndef DOXYGEN_NO_O2NS
#endif

    private:

    multi_funct_mfptr(const multi_funct_mfptr &);
    multi_funct_mfptr& operator=(const multi_funct_mfptr&);

#endif

  };
  
  /** \brief Member function pointer to a multi-dimensional function
   */
  template<class tclass, class param_t, class vec_t=boost::numeric::ublas::vector<double> >
    class multi_funct_mfptr_param : public multi_funct<vec_t> {
    public:
  
    /** \brief Specify the member function pointer
     */
    multi_funct_mfptr_param(tclass *tp, double (tclass::*fp)
			    (size_t nv, const vec_t &x, param_t &),
			    param_t &pa) {
      tptr=tp;
      fptr=fp;
      pp=&pa;
    }
  
    virtual ~multi_funct_mfptr_param() {}
  
    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
    virtual double operator()(size_t nv, const vec_t &x) {
      return (*tptr.*fptr)(nv,x,*pp);
    }

#ifndef DOXYGEN_INTERNAL
  
    protected:
  
    /// Store the function pointer
    double (tclass::*fptr)(size_t nv, const vec_t &x, param_t &);

    /// Parameter
    param_t *pp;

    /// Store a pointer to the class instance
    tclass *tptr;
  
#ifndef DOXYGEN_NO_O2NS
#endif

    private:

    multi_funct_mfptr_param(const multi_funct_mfptr_param &);
    multi_funct_mfptr_param& operator=(const multi_funct_mfptr_param&);

#endif

  };
  
  /** \brief Const member function pointer to a multi-dimensional function
   */
  template<class tclass, class vec_t=boost::numeric::ublas::vector<double> >
    class multi_funct_cmfptr : public multi_funct<vec_t> {

    public:
  
    /** \brief Specify the member function pointer
     */
    multi_funct_cmfptr(tclass *tp, double (tclass::*fp)
		       (size_t nv, const vec_t &x) const) {
			
      tptr=tp;
      fptr=fp;
    }
  
    virtual ~multi_funct_cmfptr() {}
  
    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
    virtual double operator()(size_t nv, const vec_t &x) {
      return (*tptr.*fptr)(nv,x);
    }

#ifndef DOXYGEN_INTERNAL
  
    protected:
  
    /// Store the function pointer
    double (tclass::*fptr)(size_t nv, const vec_t &x) const;

    /// Store a pointer to the class instance
    tclass *tptr;
  
#ifndef DOXYGEN_NO_O2NS
#endif

    private:

    multi_funct_cmfptr(const multi_funct_cmfptr &);
    multi_funct_cmfptr& operator=(const multi_funct_cmfptr&);

#endif

  };

  /** \brief A multi-dimensional function from a string
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class multi_funct_strings : public multi_funct<vec_t> {

  public:
  
  /** \brief Specify the string and the parameters
   */
  multi_funct_strings(std::string formula, int nv, std::string vars, 
		      int np=0, std::string parms="") {
    
    if(np<1) {
      fpw.Parse(formula,vars);
      st_np=0;
      st_parms="";
    } else {
      std::string all=vars+","+parms;
      fpw.Parse(formula,all);
      st_np=np;
      st_parms=parms;
      arr=new double[np];
    }
    st_form=formula;
    st_vars=vars;
    st_nv=nv;
  }
  
  /** \brief Specify the string and the parameters
   */
  int set_function(std::string formula, int nv, std::string vars, 
		   int np=0, std::string parms="") {

      if(np<1) {
	fpw.Parse(formula,vars);
	st_np=0;
	st_parms="";
      } else {
	std::string all=vars+","+parms;
	fpw.Parse(formula,all);
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_form=formula;
      st_vars=vars;
      st_nv=nv;
      return 0;
    }

    virtual ~multi_funct_strings() {
      if (st_np>0) {
	delete[] arr;
      }
    };
  
    /** \brief Set the values of the auxilliary parameters that were
	specified in \c parms  in the constructor
    */
    int set_parms(const vec_t &p) {
      for(int i=0;i<st_np;i++) {
	arr[i]=p[i];
      }
      return 0;
    }

    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
  virtual double operator()(size_t nv, const vec_t &x) {
      int i;
      double y;
      if(st_np<1) {
	y=fpw.Eval(x);
      } else {
	double *all=new double[st_np+st_nv];
	for(i=0;i<st_nv;i++) {
	  all[i]=x[i];
	}
	for(i=st_nv;i<st_np+st_nv;i++) {
	  all[i]=arr[i-st_nv];
	}
	y=fpw.Eval(all);
	delete[] all;
      }
      return y;
    }

#ifndef DOXYGEN_INTERNAL

    protected:

    /// The object for evaluating strings
    FunctionParser fpw;

    /// The number of parameters
    int st_np;
    /// The number of variables
    int st_nv;
    /// Storage for the parameters for \ref fpw
    double *arr;
    /// The formula string
    std::string st_form;
    /// The variable string
    std::string st_vars;
    /// The parameter string
    std::string st_parms;
  
    multi_funct_strings() {}
  
#ifndef DOXYGEN_NO_O2NS
#endif

    private:

    multi_funct_strings(const multi_funct_strings &);
    multi_funct_strings& operator=(const multi_funct_strings&);

#endif

  };

#if !defined (O2SCL_NO_CPP11) || defined (DOXYGENP)

  /** \brief A multi-dimensional function from a string
  */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class multi_funct11_strings {
    
  public:
  
  /** \brief Specify the string and the parameters
   */
  multi_funct11_strings(std::string formula, int nv, std::string vars, 
		      int np=0, std::string parms="") {
    
    if(np<1) {
      fpw.Parse(formula,vars);
      st_np=0;
      st_parms="";
    } else {
      std::string all=vars+","+parms;
      fpw.Parse(formula,all);
      st_np=np;
      st_parms=parms;
      arr=new double[np];
    }
    st_form=formula;
    st_vars=vars;
    st_nv=nv;
  }
  
  /** \brief Specify the string and the parameters
   */
  int set_function(std::string formula, int nv, std::string vars, 
		   int np=0, std::string parms="") {

      if(np<1) {
	fpw.Parse(formula,vars);
	st_np=0;
	st_parms="";
      } else {
	std::string all=vars+","+parms;
	fpw.Parse(formula,all);
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_form=formula;
      st_vars=vars;
      st_nv=nv;
      return 0;
    }

    virtual ~multi_funct11_strings() {
      if (st_np>0) {
	delete[] arr;
      }
    };
  
    /** \brief Set the values of the auxilliary parameters that were
	specified in \c parms  in the constructor
    */
    int set_parms(const vec_t &p) {
      for(int i=0;i<st_np;i++) {
	arr[i]=p[i];
      }
      return 0;
    }

    /** \brief Compute a function \c y of \c nv variables stored in \c x
	with parameter \c pa.
    */
  virtual double operator()(size_t nv, const vec_t &x) {
      int i;
      double y;
      if(st_np<1) {
	y=fpw.Eval(x);
      } else {
	double *all=new double[st_np+st_nv];
	for(i=0;i<st_nv;i++) {
	  all[i]=x[i];
	}
	for(i=st_nv;i<st_np+st_nv;i++) {
	  all[i]=arr[i-st_nv];
	}
	y=fpw.Eval(all);
	delete[] all;
      }
      return y;
    }

#ifndef DOXYGEN_INTERNAL

    protected:

    /// The object for evaluating strings
    FunctionParser fpw;

    /// The number of parameters
    int st_np;
    /// The number of variables
    int st_nv;
    /// Storage for the parameters for \ref fpw
    double *arr;
    /// The formula string
    std::string st_form;
    /// The variable string
    std::string st_vars;
    /// The parameter string
    std::string st_parms;
  
    multi_funct11_strings() {}
  
#ifndef DOXYGEN_NO_O2NS
#endif

    private:

    multi_funct11_strings(const multi_funct11_strings &);
    multi_funct11_strings& operator=(const multi_funct11_strings&);

#endif

  };

#endif

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
