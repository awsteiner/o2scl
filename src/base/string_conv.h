/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2023, Andrew W. Steiner
  
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
#include <codecvt>

#include <boost/multiprecision/cpp_dec_float.hpp>
#ifdef O2SCL_MPFR
#include <boost/multiprecision/mpfr.hpp>
#endif

// For numeric_limits for dtos()
#include <limits>

// For screenify()
#include <o2scl/misc.h>
#include <o2scl/rng.h>

namespace o2scl {

  /// \name Functions in src/base/string_conv.h
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

  /** \brief Convert a floating-point number to a string 

      This uses a \c ostringstream object to convert the
      floating-point number to a string. The value of
      <tt>std::numeric_limits::max_digits10</tt> is used to determine
      the maximum precision. If \c prec is greater than this maximum
      value, then the maximum is used. If \c prec is 0, then
      <tt>std::numeric_limits::digits10</tt> is used. The 
      default value of \c prec is 6. 
      
      If \c auto_prec is false (the default), then the number is
      converted to a string in the <tt>ios::scientific</tt> mode,
      otherwise, neither the scientific or fixed mode flags are set
      and the number is converted to a string in "automatic" mode.
  */
  template<class fp_t>
  std::string dtos(const fp_t &x, int prec=6, bool auto_prec=false) {
    
    std::ostringstream strout;
    
    size_t max=std::numeric_limits<fp_t>::max_digits10;
    size_t dig=std::numeric_limits<fp_t>::digits10;

    if (prec>((int)max)) prec=((int)max);
    if (prec==-1) prec=max;
    if (prec==0) prec=dig;
    
    if (!auto_prec) strout.setf(std::ios::scientific);
    strout.precision(prec);

    if (strout << x) {
      return strout.str();
    }
    
    O2SCL_ERR2("Conversion from floating point value to string failed in ",
	       "dtos(fp_t,int,bool).",exc_einval);
    return "";
  }

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

  /** \brief Given a floating-point number, extract the exponent
      and mantissa separately

      This function ensures that the mantissa is always 
      greater than or equal to 1 and less than 10.
   */
  template<class fp_t>
  void float_expo_mant(fp_t x, int &expo, fp_t &mant) {
    
    // Compute exponent and mantissa separately
    expo=((int)log10(x));
    mant=x/pow(10,expo);
    
    // Occasionally, finite precision errors compute the 
    // mantissa incorrectly. Fix this here.
    if (mant<1) {
      mant*=10;
      expo-=1;
    }
    if (mant>=10) {
      mant/=10;
      expo+=1;
    }
    return;
  }
  
  /** \brief Convert a value and an uncertainty to a string, 
      e.g. "1.72634(34)e-12"
  */
  std::string unc_to_string(double val, double err, int verbose=0);
  
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

  /** \brief Convert a string to a double returning non-zero
      value for failure
  */
  int s32tod_nothrow(std::u32string s, double &result);

  /** \brief Convert a string to a long double returning non-zero
      value for failure
  */
  int s32tod_nothrow(std::u32string s, long double &result);

  /** \brief Convert a string to a multiprecision boost number
   */
  template<class fp_t> int s32tod_nothrow
  (std::u32string s, fp_t &result) {
    
    std::string s2;
    bool done=false;
    for (size_t i=0;i<s.length() && done==false;i++) {
      if (s[i]<128) {
        s2+=s[i];
      } else {
        done=true;
      }
    }
    
    fp_t ret(s2);
    
    result=ret;
    
    return 0;
  }
  
  /** \brief Find out if the number pointed to by \c x has a minus sign
      
      This function returns true if the number pointed to by \c x has
      a minus sign, determining this by converting the number to a
      string. It is useful, for example, in distinguishing "-0.0" from
      "+0.0".
  */
  bool has_minus_sign(double *x);

  /** \brief Return true if the string \c s is likely a integral or
      floating point number
      
      \note The test employed is not exhaustive and this function may
      return \c true for some numbers and may return \c false for some
      non-numbers.
  */
  bool is_number(std::string s);

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

      \note The string vector is not emptied before processing,
      and entries from \c str are added to the end of \c sv. 

      \note The parsing algorithm here is simple-minded and can
      produce unexpected results in several ways. For example, it may
      not properly handle nested quotes, like <tt>"""test" test2"
      test3"</tt>.

      \verbatim embed:rst

      .. todo:: 

         In function split_string(), the rules surrounding spaces and
         quotes are not well documented.

         - Future: Replace with a better algorithm. Should quotes be
           escaped?

      \endverbatim
  */
  void split_string(std::string str, std::vector<std::string> &sv);

  /** \brief Split a string into parts using a delimiter
   */
  int split_string_delim(std::string str, std::vector<std::string> &list,
			 char delim);
  
  /** \brief Rewrap a string into a single column, avoiding
      strings less than a particular number of characters
  */
  void rewrap(std::string str, std::vector<std::string> &sv,
	      size_t ncol=79);

  /** \brief Rewrap string \c str splitting at spaces and place in \c
      sv, but ignore vt100 characters which do not occupy space in the
      terminal
  */
  void rewrap_ignore_vt100(std::string str,
			   std::vector<std::string> &sv,
			   size_t ncol=79);

  /** \brief In string \c s, replace all occurrences of \c s1
      with string \c s2, and return the number of replacements

      \note If the string \c s1 can be found inside \c s2, then
      this would lead to an infinite loop, so the error handler
      is called.
  */
  size_t string_replace(std::string &s, const std::string &s1,
                        const std::string &s2);
  
  /** \brief Convert from UTF-8 to 32-bit integers

      \warning This depends on C++ extensions that will 
      eventually be deprecated, but apparently haven't been
      replaced in C++20 yet?
  */
  void utf8_to_char32(const std::string &in,
                      std::u32string &out);

  /** \brief Convert from 32-bit integers to UTF-8

      \warning This depends on C++ extensions that will 
      eventually be deprecated, but apparently haven't been
      replaced in C++20 yet?
  */
  void char32_to_utf8(const std::u32string &in,
                      std::string &out);
  
  /** \brief Rewrap a string into a single column, avoiding
      strings less than a particular number of characters

      This function is used to format the help output 
      in \ref o2scl::cli .

      Note that this treats whitespace other than ' ' and 
      '\n' as it does normal characters.
  */
  void rewrap_keep_endlines(std::string str, std::vector<std::string> &sv,
			    size_t ncol=79, int verbose=0,
			    bool ignore_vt100=true);

  /** \brief Convert a string-based list of unsigned integers
      to a list
  */
  template<class size_vec_t>
  int string_to_uint_list(const std::string &x,
                          size_vec_t &list) {
    
    list.clear();
    std::vector<std::string> ranges;
    size_t k=0;
    while (k<x.length()) {
      size_t loc=x.find(',',k);
      if (loc!=std::string::npos) {
	std::string stemp=x.substr(k,loc-k);
	ranges.push_back(stemp);
	k+=stemp.length()+1;
      } else {
	if (k<x.length()) {
	  ranges.push_back(x.substr(k,x.length()-k));
	}
	k=x.length();
      }
    }
    size_t uitmp, uitmp2;
    for(size_t j=0;j<ranges.size();j++) {
      if (ranges[j].find('-')==std::string::npos) {
	int ret=stoszt_nothrow(ranges[j],uitmp);
	if (ret!=0) return ret;
	list.push_back(uitmp);
      } else {
	size_t loc=ranges[j].find('-');
	std::string sstart=ranges[j].substr(0,loc);
	std::string send=ranges[j].substr(loc+1,ranges[j].size()-loc);
	int ret=stoszt_nothrow(sstart,uitmp);
	if (ret!=0) return ret;
	ret=stoszt_nothrow(send,uitmp2);
	if (ret!=0) return ret;
	for(size_t jk=uitmp;jk<=uitmp2;jk++) {
	  list.push_back(jk);
	}
      }
    }
    
    return 0;
  }
  //@}

  /** \brief Parse \c line into \c entries using the FORTRAN
      format string \c format
  */
  void parse_fortran_format(std::string line, std::string format,
                            std::vector<std::string> &entries);

  /// Copy string \c s to character array \c x of length \c len
  void string_to_char_array(std::string s, char *x, int len);
  
}

#endif

