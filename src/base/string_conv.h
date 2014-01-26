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
#ifndef O2SCL_STRING_CONV_H
#define O2SCL_STRING_CONV_H
/** \file string_conv.h
    \brief Various string conversion functions
*/

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// \name Functions in string_conv.h
  //@{
  /** \brief Convert a pointer to a string 

      This uses an \c ostringstream to convert a pointer to a string
      and is architecture-dependent.
  */
  std::string ptos(void *p);

  /** \brief Convert an integer to a string 
   */
  std::string itos(int x);

  /** \brief Convert a size_t to a string 
   */
  std::string szttos(size_t x);

  /** \brief Convert an integer to a string (exception-free version)
   */
  std::string itos_nothrow(int x) throw();

  /** \brief Convert a boolean value to a string 

      This returns \c "1" for true and \c "0" for false.
   */
  std::string btos(bool b);

  /** \brief Convert a double to a string 

      If \c auto_prec is false, then the number is converted to 
      a string in the <tt>ios::scientific</tt> mode, otherwise,
      neither the scientific or fixed mode flags are set and the
      number is converted to a string in "automatic" mode.
   */
  std::string dtos(double x, int prec=6, bool auto_prec=false);

  /** \brief Returns the number of characters required to display the 
      exponent of \c x in scientific mode
      
      This returns 2 or 3, depending on whether or not the absolute
      magnitude of the exponent is greater than or equal to 100. It
      uses <tt>stringstream</tt> to convert the number to a string and
      counts the number of characters directly.
   */
  size_t size_of_exponent(double x);

  /** \brief Convert a double to a string using a specified format
   */
  std::string dtos(double x, std::ostream &format);

  /** \brief Convert a string to an integer 
      
      If \c err_on_fail is true and the conversion fails, this
      function calls the error handler, otherwise this function just
      returns zero.

      \warning Note that this is currently different than
      <tt>std::stoi()</tt> so you may have to explicitly provide the
      namespace, i.e. <tt>o2scl::stoi()</tt> in your code.
  */
  int stoi(std::string s, bool err_on_fail=true);

  /** \brief Convert a string to a size_t
      
      If \c err_on_fail is true and the conversion fails, this
      function calls the error handler, otherwise this function just
      returns zero.
  */
  size_t stoui(std::string s, bool err_on_fail=true);

  /** \brief Convert a string to a boolean value
      
      This returns true if only if the string has at least one
      character and the first non-whitespace character is either \c t,
      \c T, or one of the numbers 1 through 9.
      
      If \c err_on_fail is true and the conversion fails, this
      function calls the error handler, otherwise this function just
      returns false.
  */
  bool stob(std::string s, bool err_on_fail=true);

  /** \brief Convert a string to a double 

      If \c err_on_fail is true and the conversion fails, this
      function calls the error handler, otherwise this function just
      returns 0.0.

      \warning Note that this is currently different than
      <tt>std::stod()</tt> so you may have to explicitly provide the
      namespace, i.e. <tt>o2scl::stod()</tt> in your code.
  */
  double stod(std::string s, bool err_on_fail=true);

  /** \brief Find out if the number pointed to by \c x has a minus sign
      
      This function returns true if the number pointed to by \c x has
      a minus sign using the GSL IEEE functions. It is useful, for
      example, in distinguishing "-0.0" from "+0.0".
  */
  bool has_minus_sign(double *x);

  /** \brief Return true if the string \c s is likely a integral or
      floating point number
      
      \note The test employed is not exhaustive and this function may
      return \c true for some numbers and may return \c false for some
      non-numbers.
   */
  bool is_number(std::string s);

  /** \brief Convert a formula to a double 
      
      This uses \c FunctionParser to convert strings like "-1.0e-3",
      "1.0/3.0" and "exp(cos(-1.0e-2))" to floating point numbers.
  */
  double function_to_double(std::string s, bool err_on_fail=true);

  /** \brief Split a string into words using whitespace for delimiters and 
      (partially) respecting quotes

      \todo 
      - More documentation
      - Add user-specified delimiters?
      - Add version which ignores quotes
      - Use this function in acol
  */
  void split_string(std::string str, std::vector<std::string> &sv);
  //@}

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

