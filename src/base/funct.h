/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2016, Andrew W. Steiner
  
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

#include <gsl/gsl_math.h>

#include <boost/numeric/ublas/vector.hpp>

#include <o2scl/shunting_yard.h>

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
    funct11_string(std::string expr, std::string var) {
      calc.compile(expr.c_str(),&vars);
      st_form=expr;
      st_var=var;
    }

    virtual ~funct11_string() {
    };

  
    /** \brief Specify the string and the parameters
     */
    int set_function(std::string expr, std::string var) {
      calc.compile(expr.c_str(),&vars);
      st_form=expr;
      st_var=var;
      return 0;
    }

    /** \brief Set the values of the auxilliary parameters that were
	specified in \c parms in the constructor
    */
    int set_parm(std::string name, double val) {
      vars[name]=val;
      return 0;
    }
    
    /** \brief Compute the function at point \c x and return the result
     */
    virtual double operator()(double x) const {
      vars[st_var]=x;
      return calc.eval(&vars);
    }

#ifndef DOXYGEN_INTERNAL

  protected:

    /// The object for evaluating strings
    mutable calculator calc;

    /// Parameter map
    mutable std::map<std::string,double> vars;
    
    /// The expr
    std::string st_form;
    /// The variable
    std::string st_var;

    funct11_string() {};

#endif
#ifndef DOXYGEN_NO_O2NS

  private:

    funct11_string(const funct11_string &);
    funct11_string& operator=(const funct11_string&);

#endif

  };

  /** \brief A wrapper to specify \ref o2scl::funct11 objects 
      to GSL
   */
  class funct_gsl : public gsl_function {

  protected:
    
    /// The function wrapper
    static double funct_wrap(double x, void *params) {
      funct11 *fp=(funct11 *)params;
      return (*fp)(x);
    }

  public:

    /// Create an object based on the specified function, \c f
    funct_gsl(funct11 &f) {
      function=&funct_wrap;
      params=&f;
    }
    
  };


#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
