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
      double y;
      if(st_np<1) {
	y=fpw.Eval(&x);
      } else {
	double *all=new double[st_np+1];
	all[0]=x;
	for(size_t i=1;i<=st_np;i++) all[i]=arr[i-1];
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
