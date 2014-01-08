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
#ifndef O2SCL_ODE_FUNCT_H
#define O2SCL_ODE_FUNCT_H

/** \file ode_funct.h
    \brief Function object classes for ODE functions
*/

#include <string>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/fparser.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

#if !defined (O2SCL_NO_CPP11) || defined (DOXYGENP)

  //template<class vec_y_t, class vec_dydt_t> using ode_funct11 = 
  //std::function<int(double,size_t,const vec_y_t &,vec_dydt_t &)>;

  /** \brief Ordinary differential equation function (C++11 version)
   */
  typedef std::function<int(double,size_t,
			    const boost::numeric::ublas::vector<double> &,
			    boost::numeric::ublas::vector<double> &)> 
    ode_funct11;
  
#endif

  /** \brief Ordinary differential equation function [abstract base]

      This base class provides the basic format for specifying
      ordinary differential equations to integrate with the \o2 ODE
      solvers. Select the appropriate child of this class according
      to the kind of functions which are to be given to the solver.

      \future One could, in principle, allow the function values
      and derivatives to have different vector types, i.e.
      <tt>const vec_t &y, vec2_t &dydx</tt>. 
  */
  template <class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t> class ode_funct {

  public:

  ode_funct() {}

  virtual ~ode_funct() {}

  /** \brief Compute the \c nv derivatives as a function of the \c nv 
      functions specified in \c y at the point \c x.
  */
  virtual int operator()(double x, size_t nv, const vec_y_t &y, 
			 vec_dydx_t &dydx)=0;
  
#ifndef DOXYGEN_INTERNAL

  private:

  ode_funct(const ode_funct &);
  ode_funct& operator=(const ode_funct&);

#endif

  };

  /** \brief Provide ODE functions in the form of function pointers
  */
  template <class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t> class ode_funct_fptr : 
    public ode_funct<vec_y_t,vec_dydx_t> {

  public:

  /** \brief Create an object given a function pointer
   */
  ode_funct_fptr(int (*fp)(double, size_t, const vec_y_t &, 
			   vec_dydx_t &)) {
    fptr=fp;
  }

  virtual ~ode_funct_fptr() {}

  /** \brief Compute the \c nv derivatives as a function of the \c nv 
      functions specified in \c y at the point \c x.
  */
  virtual int operator()(double x, size_t nv, const vec_y_t &y, 
			 vec_dydx_t &dydx) {
    return fptr(x,nv,y,dydx);
  }
  
#ifndef DOXYGEN_INTERNAL

  protected:

  ode_funct_fptr() {};
  
  /// The function pointer
  int (*fptr)(double x, size_t nv, const vec_y_t &y, 
	      vec_dydx_t &dydx);
  
  private:
  
  ode_funct_fptr(const ode_funct_fptr &);
  ode_funct_fptr& operator=(const ode_funct_fptr&);

#endif

  };

  /** \brief Provide ODE functions in the form of function pointers
   */
  template <class param_t, 
    class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t> class ode_funct_fptr_param : 
    public ode_funct<vec_y_t,vec_dydx_t> {

  public:

  /** \brief Create an object given a function pointer
   */
  ode_funct_fptr_param(int (*fp)(double, size_t, const vec_y_t &, 
			   vec_dydx_t &, param_t &), param_t &pa) {
    fptr=fp;
    pp=&pa;
  }

  virtual ~ode_funct_fptr_param() {}

  /** \brief Compute the \c nv derivatives as a function of the \c nv 
      functions specified in \c y at the point \c x.
  */
  virtual int operator()(double x, size_t nv, const vec_y_t &y, 
			 vec_dydx_t &dydx) {
    return fptr(x,nv,y,dydx,*pp);
  }
  
#ifndef DOXYGEN_INTERNAL

  protected:

  ode_funct_fptr_param() {};

  /// The function pointer
  int (*fptr)(double x, size_t nv, const vec_y_t &y, 
	      vec_dydx_t &dydx, param_t &);

  /// The parameter
  param_t *pp;
  
  private:
  
  ode_funct_fptr_param(const ode_funct_fptr_param &);
  ode_funct_fptr_param& operator=(const ode_funct_fptr_param&);

#endif

  };

  /** \brief Provide ODE functions in the form of member 
      function pointers
  */
  template <class tclass, class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t> class ode_funct_mfptr : 
    public ode_funct<vec_y_t,vec_dydx_t> {
    
  public:
    
  /** \brief Create an object given a class and member function pointer
   */
  ode_funct_mfptr
  (tclass *tp, int (tclass::*fp)(double x, size_t nv, 
				 const vec_y_t &y, vec_dydx_t &dydx)) {
    tptr=tp;
    fptr=fp;
  }
  
  virtual ~ode_funct_mfptr() {};
  
  /** \brief Compute the \c nv derivatives as a function of the \c nv 
      functions specified in \c y at the point \c x.
  */
  virtual int operator()(double x, size_t nv, const vec_y_t &y, 
			 vec_dydx_t &dydx) {
    return (*tptr.*fptr)(x,nv,y,dydx);
  }
  
#ifndef DOXYGEN_INTERNAL

  protected:

  /// The pointer to the member function
  int (tclass::*fptr)(double x, size_t nv, const vec_y_t &y, 
		      vec_dydx_t &dydx);

  /// The pointer to the class
  tclass *tptr;

  private:

  ode_funct_mfptr(const ode_funct_mfptr &);
  ode_funct_mfptr& operator=(const ode_funct_mfptr&);

#endif

  };

  /** \brief Provide ODE functions in the form of member 
      function pointers
  */
  template <class tclass, class param_t, 
    class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t> class ode_funct_mfptr_param : 
    public ode_funct<vec_y_t,vec_dydx_t> {
    
  public:
    
  /** \brief Create an object given a class and member function pointer
   */
  ode_funct_mfptr_param
  (tclass *tp, 
   int (tclass::*fp)(double x, size_t nv, 
		     const vec_y_t &y, vec_dydx_t &dydx, 
		     param_t &), param_t &pa) {
    tptr=tp;
    fptr=fp;
    pp=&pa;
  }
  
  virtual ~ode_funct_mfptr_param() {};
  
  /** \brief Compute the \c nv derivatives as a function of the \c nv 
      functions specified in \c y at the point \c x.
  */
  virtual int operator()(double x, size_t nv, const vec_y_t &y, 
			 vec_dydx_t &dydx) {
    return (*tptr.*fptr)(x,nv,y,dydx,*pp);
  }
  
#ifndef DOXYGEN_INTERNAL

  protected:

  /// The pointer to the member function
  int (tclass::*fptr)(double x, size_t nv, const vec_y_t &y, 
		      vec_dydx_t &dydx, param_t &);

  /// The pointer to the class
  tclass *tptr;

  /// The parameter
  param_t *pp;

  private:

  ode_funct_mfptr_param(const ode_funct_mfptr_param &);
  ode_funct_mfptr_param& operator=(const ode_funct_mfptr_param&);

#endif

  };

  /** \brief Provide ODE functions in the form of const member 
      function pointers
  */
  template <class tclass, class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t> class ode_funct_cmfptr : 
    public ode_funct<vec_y_t,vec_dydx_t> {
    
  public:
    
  /** \brief Create an object given a class and member function pointer
   */
  ode_funct_cmfptr
  (tclass *tp, int (tclass::*fp)(double x, size_t nv, 
				 const vec_y_t &y, 
				 vec_dydx_t &dydx) const) {
    tptr=tp;
    fptr=fp;
  }
  
  virtual ~ode_funct_cmfptr() {};
  
  /** \brief Compute the \c nv derivatives as a function of the \c nv 
      functions specified in \c y at the point \c x.
  */
  virtual int operator()(double x, size_t nv, const vec_y_t &y, 
			 vec_dydx_t &dydx) {
    return (*tptr.*fptr)(x,nv,y,dydx);
  }
  
#ifndef DOXYGEN_INTERNAL

  protected:

  /// The pointer to the member function
  int (tclass::*fptr)(double x, size_t nv, const vec_y_t &y, 
		      vec_dydx_t &dydx) const;

  /// The pointer to the class
  tclass *tptr;

  private:

  ode_funct_cmfptr(const ode_funct_cmfptr &);
  ode_funct_cmfptr& operator=(const ode_funct_cmfptr&);

#endif

  };

  /** \brief One-dimensional function from strings
      \nothing
  */
  template <class vec_y_t=boost::numeric::ublas::vector<double>,
    class vec_dydx_t=vec_y_t>
    class ode_funct_strings : public ode_funct<vec_y_t,vec_dydx_t> {
    public:

    /** \brief Specify the string and the parameters
     */
    ode_funct_strings(size_t nv, std::string *formulas, std::string funcs, 
		      std::string var, size_t np=0, std::string parms="") {
      size_t i;
      fpw=new FunctionParser[nv];
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],funcs+","+var);
	}
	st_np=0;
	st_parms="";
      } else {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],funcs+","+var+","+parms);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_forms=formulas;
      st_var=var;
      st_funcs=funcs;
      st_nv=nv;
    }
  
    virtual ~ode_funct_strings() {
      if (st_np>0) {
	delete[] arr;
      }
      delete[] fpw;
    }
  
    /** \brief Specify the string and the parameters
     */
    int set_function(size_t nv, std::string *formulas, std::string funcs, 
		     std::string var, size_t np=0, std::string parms="") {
      size_t i;
      if (nv!=st_nv) {
	delete[] fpw;
	fpw=new FunctionParser[nv];
      }
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],funcs+","+var);
	}
	st_np=0;
	st_parms="";
      } else {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],funcs+","+var+","+parms);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_forms=formulas;
      st_var=var;
      st_funcs=funcs;
      st_nv=nv;
      return 0;
    }

    /** \brief Set the values of the auxilliary parameters that were
	specified in 'parms' in the constructor
    */
    template<class vec_p_t> void set_parms(vec_p_t &p) {
      for(size_t i=0;i<st_np;i++) {
	arr[i]=p[i];
      }
      return;
    }
  
    virtual int operator()(double x, size_t nv, const vec_y_t &y, 
			   vec_dydx_t &dydx) {
      size_t i;
      if(st_np<1) {
	double *all=new double[nv+1];
	for(i=0;i<nv;i++) all[i]=y[i];
	all[nv]=x;
	for(i=0;i<st_nv;i++) {
	  dydx[i]=fpw[i].Eval(all);
	}
	delete[] all;
      } else {
	double *all=new double[st_np+st_nv+1];
	for(i=0;i<st_nv;i++) all[i]=y[i];
	all[st_nv]=x;
	for(i=st_nv+1;i<st_np+st_nv+1;i++) all[i]=arr[i-st_nv-1];
	for(i=0;i<st_nv;i++) {
	  dydx[i]=fpw[i].Eval(all);
	}
	delete[] all;
      }
      return 0;
    }

#ifndef DOXYGEN_INTERNAL

    protected:

    /// The formula parser
    FunctionParser *fpw;
    /// The number of parameters
    size_t st_np;
    /// The number of variables
    size_t st_nv;
    /// The parameters
    double *arr;
    /// The formulas
    std::string *st_forms;
    /// The variables
    std::string st_var;
    /// The parameters
    std::string st_parms;
    /// The function names
    std::string st_funcs;

    ode_funct_strings() {};

    private:

    ode_funct_strings(const ode_funct_strings &);
    ode_funct_strings& operator=(const ode_funct_strings&);

#endif

  };


#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
