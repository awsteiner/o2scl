/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2018, Andrew W. Steiner
  
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

// For screenify()
#include <o2scl/misc.h>

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
      
      This function is now just a wrapper for <tt>std::stoi</tt>.
  */
  int stoi(std::string s);

  /** \brief Convert a string to an integer without throwing an 
      exception
   */
  int stoi_nothrow(std::string s, int &result);

  /** \brief Convert a string to a size_t
      
      If \c err_on_fail is true and the conversion fails, this
      function calls the error handler, otherwise this function just
      returns zero.
  */
  size_t stoszt(std::string s);

  /** \brief Convert a string to a size_t without throwing an
      exception
      
      If \c err_on_fail is true and the conversion fails, this
      function calls the error handler, otherwise this function just
      returns zero.
  */
  int stoszt_nothrow(std::string s, size_t &result);

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

      This function is now just a wrapper for <tt>std::stod</tt>. 

      \warning Because of the presence of <tt>std::stod()</tt> in 
      C++11, you may have to explicitly provide the
      namespace, i.e. <tt>o2scl::stod()</tt> in your code.
  */
  double stod(std::string s);

  /** \brief Convert a string to a double returning non-zero
      value for failure
  */
  int stod_nothrow(std::string s, double &result);

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
      
      This function removes all quotes and apostrophes from the string
      and then uses \ref o2scl::calculator to convert strings like
      "-1.0e-3", "1.0/3.0" and "exp(cos(-1.0e-2))" to floating point
      numbers.
  */
  double function_to_double(std::string s);

  /** \brief Split a string into words using whitespace for delimiters
      and (partially) respecting quotes

      This function separates a string into words, and handles words
      that begin with a <tt>"</tt> by adding more words until finding
      one which ends with another <tt>"</tt>. Strings like
      \code
      this is a test
      \endcode
      get parsed as "this", "is", "a", "test" and strings like
      \code
      "this is" a test
      \endcode
      get parsed as "this is", "a", "test".

      This is used to reformat command descriptions and help text for
      the screen width in cli::comm_option_help(), to process lines
      read from a file in cli::comm_option_run(), and to process input
      in cli::run_interactive().
      
      \future Add user-specified delimiters?
  */
  void split_string(std::string str, std::vector<std::string> &sv);

  /** \brief Rewrap a string into a single column, avoiding
      strings less than a particular number of characters
   */
  void rewrap(std::string str, std::vector<std::string> &sv,
	      size_t ncol=79);
  //@}

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

