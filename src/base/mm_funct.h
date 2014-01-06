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
#ifndef O2SCL_MM_FUNCT_H
#define O2SCL_MM_FUNCT_H

/** \file mm_funct.h
    \brief Function object classes for multi-dimensional functions
*/

#include <string>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/fparser.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

#if defined (O2SCL_CPP11) || defined (DOXYGENP)
  
  /// Array of multi-dimensional functions typedef (C++11 version)
  typedef std::function<
    int(size_t,const boost::numeric::ublas::vector<double> &,
	boost::numeric::ublas::vector<double> &) > mm_funct11;

#endif

  /** \brief Array of multi-dimensional functions [abstract base]
      
      This class generalizes \c nv functions of \c nv variables, i.e.
      \f$ y_j(x_0,x_1,...,x_{nv-1}) \f$ for \f$ 0\leq j \leq 
      nv-1 \f$ .

      This class is one of a large number of function object classes
      in \o2 designed to provide a mechanism for the user to 
      supply functions to solvers, minimizers, integrators, etc.
      See \ref funct_section for a general description.
  */
  template<class vec_t=boost::numeric::ublas::vector<double> > 
    class mm_funct {

    public:  

    mm_funct() {}
    
    virtual ~mm_funct() {}
    
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    virtual int operator()(size_t nv, const vec_t &x, vec_t &y)=0;
    
#ifndef DOXYGEN_NO_O2NS

    private:

    mm_funct(const mm_funct &);
    mm_funct& operator=(const mm_funct&);

#endif

  };
  
  /** \brief Function pointer to array of multi-dimensional functions
   */
  template<class vec_t=boost::numeric::ublas::vector<double> > 
    class mm_funct_fptr : public mm_funct<vec_t> {
    
    public:
    
    mm_funct_fptr() {}

    virtual ~mm_funct_fptr() {}

    /** \brief Specify the function pointer
     */
    mm_funct_fptr(int (*fp)(size_t nv, const vec_t &x, vec_t &y)) {
      fptr=fp;
    }
    
    /** \brief Specify the function pointer
     */
    int set_function(int (*fp)(size_t nv, const vec_t &x, vec_t &y)) {
      fptr=fp;
      return 0;
    }
    
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    virtual int operator()(size_t nv, const vec_t &x, vec_t &y) {
      return fptr(nv,x,y);
    }
    
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// The function pointer to the user-supplied function
    int (*fptr)(size_t nv, const vec_t &x, vec_t &y);
    
    private:
    
    mm_funct_fptr(const mm_funct_fptr &);
    mm_funct_fptr& operator=(const mm_funct_fptr&);
    
#endif
    
  };
  
  /** \brief Function pointer to array of multi-dimensional functions
   */
  template<class param_t, class vec_t=boost::numeric::ublas::vector<double> > 
    class mm_funct_fptr_param : public mm_funct<vec_t> {
    
    public:
    
    mm_funct_fptr_param() {}

    virtual ~mm_funct_fptr_param() {}

    /** \brief Specify the function pointer
     */
    mm_funct_fptr_param(int (*fp)
			(size_t nv, const vec_t &x, vec_t &y, param_t &),
			param_t &pa) {
      fptr=fp;
      pp=&pa;
    }
    
    /** \brief Specify the function pointer
     */
    int set_function(int (*fp)(size_t nv, const vec_t &x, vec_t &y, param_t &),
		     param_t &pa) {
      fptr=fp;
      pp=&pa;
      return 0;
    }
    
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    virtual int operator()(size_t nv, const vec_t &x, vec_t &y) {
      return fptr(nv,x,y,*pp);
    }
    
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// The function pointer to the user-supplied function
    int (*fptr)(size_t nv, const vec_t &x, vec_t &y, param_t &);
    
    /// The parameter
    param_t *pp;

    private:
    
    mm_funct_fptr_param(const mm_funct_fptr_param &);
    mm_funct_fptr_param& operator=(const mm_funct_fptr_param&);
    
#endif
    
  };
  
  /** \brief Member function pointer to an array of 
      multi-dimensional functions
  */
  template<class tclass, class vec_t=boost::numeric::ublas::vector<double> >
    class mm_funct_mfptr : public mm_funct<vec_t> {
    public:
    
    /** \brief Empty constructor
     */
    mm_funct_mfptr() {
    }

    /** \brief Specify the member function pointer
     */
    mm_funct_mfptr(tclass *tp, int (tclass::*fp)
		   (size_t nv, const vec_t &x, vec_t &y)) {
      tptr=tp;
      fptr=fp;
    }

    /** \brief Specify the member function pointer
     */
    int set_function(tclass *tp, int (tclass::*fp)
		     (size_t nv, const vec_t &x, vec_t &y)) 
    {
      tptr=tp;
      fptr=fp;
      return 0;
    }
    
    virtual ~mm_funct_mfptr() {};
    
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    virtual int operator()(size_t nv, const vec_t &x, vec_t &y) {
      return (*tptr.*fptr)(nv,x,y);
    }
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// The member function pointer
    int (tclass::*fptr)(size_t nv, const vec_t &x, vec_t &y);

    /// The class pointer
    tclass *tptr;
    
    private:
    
    mm_funct_mfptr(const mm_funct_mfptr &);
    mm_funct_mfptr& operator=(const mm_funct_mfptr&);
    
#endif
    
  };

  /** \brief Member function pointer to an array of 
      multi-dimensional functions
  */
  template<class tclass, class param_t, class vec_t=boost::numeric::ublas::vector<double> >
    class mm_funct_mfptr_param : public mm_funct<vec_t> {
    public:
    
    /** \brief Empty constructor
     */
    mm_funct_mfptr_param() {
    }

    /** \brief Specify the member function pointer
     */
    mm_funct_mfptr_param(tclass *tp, int (tclass::*fp)
			 (size_t nv, const vec_t &x, vec_t &y, param_t &), 
			 param_t &pa) {
      tptr=tp;
      fptr=fp;
      pp=&pa;
    }

    /** \brief Specify the member function pointer
     */
    int set_function(tclass *tp, int (tclass::*fp)
		     (size_t nv, const vec_t &x, vec_t &y, param_t &), 
		     param_t &pa) {
      tptr=tp;
      fptr=fp;
      pp=&pa;
      return 0;
    }
    
    virtual ~mm_funct_mfptr_param() {};
    
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    virtual int operator()(size_t nv, const vec_t &x, vec_t &y) {
      return (*tptr.*fptr)(nv,x,y,*pp);
    }
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// The member function pointer
    int (tclass::*fptr)(size_t nv, const vec_t &x, vec_t &y, param_t &pa);

    /// The class pointer
    tclass *tptr;

    /// Parameter
    param_t *pp;
    
    private:
    
    mm_funct_mfptr_param(const mm_funct_mfptr_param &);
    mm_funct_mfptr_param& operator=(const mm_funct_mfptr_param&);
    
#endif
    
  };

  /** \brief Const member function pointer to an array of 
      multi-dimensional functions
  */
  template<class tclass, class vec_t=boost::numeric::ublas::vector<double> >
    class mm_funct_cmfptr : public mm_funct<vec_t> {

    public:
    
    /** \brief Empty constructor
     */
    mm_funct_cmfptr() {
    }

    /** \brief Specify the member function pointer
     */
    mm_funct_cmfptr(tclass *tp, int (tclass::*fp)
		   (size_t nv, const vec_t &x, vec_t &y) const) {
      tptr=tp;
      fptr=fp;
    }

    /** \brief Specify the member function pointer
     */
    int set_function(tclass *tp, int (tclass::*fp)
		     (size_t nv, const vec_t &x, vec_t &y) const) 
    {
      tptr=tp;
      fptr=fp;
      return 0;
    }
    
    virtual ~mm_funct_cmfptr() {};
    
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    virtual int operator()(size_t nv, const vec_t &x, vec_t &y) {
      return (*tptr.*fptr)(nv,x,y);
    }
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// The member function pointer
    int (tclass::*fptr)(size_t nv, const vec_t &x, vec_t &y) const;

    /// The class pointer
    tclass *tptr;
    
    private:
    
    mm_funct_cmfptr(const mm_funct_cmfptr &);
    mm_funct_cmfptr& operator=(const mm_funct_cmfptr&);
    
#endif
    
  };

  /** \brief Array of multi-dimensional functions in an array of strings
   */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class mm_funct_strings : public mm_funct<vec_t> {
    public:
    
    /** \brief Specify the strings
     */
    mm_funct_strings(int nv, std::string *formulas, std::string vars, 
		     int np=0, std::string parms="") {
      int i;
      fpw=new FunctionParser[nv];
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],vars);
	}
	st_np=0;
	st_parms="";
      } else {
	std::string all=vars+","+parms;
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],all);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_forms=formulas;
      st_vars=vars;
      st_nv=nv;
    }
    
    virtual ~mm_funct_strings() {
      if (st_np>0) {
	delete[] arr;
      }
      delete[] fpw;
    };
    
    /** \brief Set the values of the auxilliary parameters that were
	specified in 'parms' in the constructor
    */
    int set_parms(const vec_t &p) {
      for(int i=0;i<st_np;i++) {
	arr[i]=p[i];
      }
      return 0;
    }
    
    
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    virtual int operator()(size_t nv, const vec_t &x, vec_t &y) {
      int i;
      if(st_np<1) {
	for(i=0;i<st_nv;i++) {
	  y[i]=fpw[i].Eval(x);
	}
      } else {
	double *all=new double[st_np+st_nv];
	for(i=0;i<st_nv;i++) all[i]=x[i];
	for(i=st_nv;i<st_np+st_nv;i++) {
	  all[i]=arr[i-st_nv];
	}
	for(i=0;i<st_nv;i++) {
	  y[i]=fpw[i].Eval(all);
	}
	delete[] all;
      }
      return 0;
    }

    /// Set the functions
    int set_function(int nv, std::string *formulas, std::string vars, 
		     int np=0, std::string parms="") {
      int i;
      if (nv!=st_nv) {
	delete[] fpw;
	fpw=new FunctionParser[nv];
      }
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],vars);
	}
	st_np=0;
	st_parms="";
      } else {
	std::string all=vars+","+parms;
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],all);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np+nv];
      }
      st_forms=formulas;
      st_vars=vars;
      st_nv=nv;
      return 0;
    }
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// The function parser
    FunctionParser *fpw;

    /// The number of parameters
    int st_np;

    /// The number of variables
    int st_nv;

    /// The arguments to the function parser
    double *arr;

    /// The formulas
    std::string *st_forms;

    /// The variables
    std::string st_vars;

    /// The parameters
    std::string st_parms;
    
    mm_funct_strings() {};

    private:
    
    mm_funct_strings(const mm_funct_strings &);
    mm_funct_strings& operator=(const mm_funct_strings&);
    
#endif
    
  };
  
#if defined (O2SCL_CPP11) || defined (DOXYGENP)

  /** \brief Array of multi-dimensional functions in an array of strings
   */
  template<class vec_t=boost::numeric::ublas::vector<double> >
    class mm_funct11_strings {
    public:
    
    /** \brief Specify the strings
     */
    mm_funct11_strings(int nv, std::string *formulas, std::string vars, 
		     int np=0, std::string parms="") {
      int i;
      fpw=new FunctionParser[nv];
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],vars);
	}
	st_np=0;
	st_parms="";
      } else {
	std::string all=vars+","+parms;
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],all);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np];
      }
      st_forms=formulas;
      st_vars=vars;
      st_nv=nv;
    }
    
    virtual ~mm_funct11_strings() {
      if (st_np>0) {
	delete[] arr;
      }
      delete[] fpw;
    };
    
    /** \brief Set the values of the auxilliary parameters that were
	specified in 'parms' in the constructor
    */
    int set_parms(const vec_t &p) {
      for(int i=0;i<st_np;i++) {
	arr[i]=p[i];
      }
      return 0;
    }
    
    
    /** \brief Compute \c nv functions, \c y, of \c nv variables
	stored in \c x with parameter \c pa.
    */
    virtual int operator()(size_t nv, const vec_t &x, vec_t &y) {
      int i;
      if(st_np<1) {
	for(i=0;i<st_nv;i++) {
	  y[i]=fpw[i].Eval(x);
	}
      } else {
	double *all=new double[st_np+st_nv];
	for(i=0;i<st_nv;i++) all[i]=x[i];
	for(i=st_nv;i<st_np+st_nv;i++) {
	  all[i]=arr[i-st_nv];
	}
	for(i=0;i<st_nv;i++) {
	  y[i]=fpw[i].Eval(all);
	}
	delete[] all;
      }
      return 0;
    }

    /// Set the functions
    int set_function(int nv, std::string *formulas, std::string vars, 
		     int np=0, std::string parms="") {
      int i;
      if (nv!=st_nv) {
	delete[] fpw;
	fpw=new FunctionParser[nv];
      }
      if(np<1) {
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],vars);
	}
	st_np=0;
	st_parms="";
      } else {
	std::string all=vars+","+parms;
	for(i=0;i<nv;i++) {
	  fpw[i].Parse(formulas[i],all);
	}
	st_np=np;
	st_parms=parms;
	arr=new double[np+nv];
      }
      st_forms=formulas;
      st_vars=vars;
      st_nv=nv;
      return 0;
    }
    
#ifndef DOXYGEN_INTERNAL
    
    protected:
    
    /// The function parser
    FunctionParser *fpw;

    /// The number of parameters
    int st_np;

    /// The number of variables
    int st_nv;

    /// The arguments to the function parser
    double *arr;

    /// The formulas
    std::string *st_forms;

    /// The variables
    std::string st_vars;

    /// The parameters
    std::string st_parms;
    
    mm_funct11_strings() {};

    private:
    
    mm_funct11_strings(const mm_funct11_strings &);
    mm_funct11_strings& operator=(const mm_funct11_strings&);
    
#endif
    
  };
#endif

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
