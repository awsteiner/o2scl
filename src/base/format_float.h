/*
  -------------------------------------------------------------------

  Copyright (C) 2006-2017, Andrew W. Steiner

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
#ifndef O2SCL_FORMAT_FLOAT_H
#define O2SCL_FORMAT_FLOAT_H

/** \file format_float.h
    \brief File defining \ref o2scl::format_float
*/
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <o2scl/err_hnd.h>
#include <o2scl/misc.h>
#include <o2scl/string_conv.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief Format a floating point number into a Latex or HTML string
      
      This class formats floating point strings into something 
      useful for HTML or Latex documents. For basic use,
      simply call either \ref html_mode() or \ref latex_mode() and 
      then use \ref convert().

      The base-10 logarithm of the smallest and largest numbers to be
      represented without a string akin to "times 10 to the nth power"
      can be specified in \ref set_exp_limits(). The number
      of significant figures can be specified with set_sig_figs() (the
      default is 5).

      To force \ref convert() to adds zeros to the right side of the
      mantissa to guarantee that the requested number of significant
      digits is given, call \ref set_pad_zeros() with a <tt>true</tt>
      argument.

      To force scientific notation for all numbers, set the maximum
      exponent to be smaller than the minimum exponent. 

      \note This function does not warn the user if the number
      of significant figures requested is larger than the machine 
      precision.
      \note If the absolute magnitude for either the minimum or
      maximum exponent is larger than or equal to the number of 
      significant figures, then rounding will automatically 
      occur. 

      The format for a normal number is
      \verbatim
      prefix (sign-string) number suffix
      \endverbatim
      and in scientific notation is
      \verbatim 
      sci-prefix (sci-sign-string) number times-string 
      exp-prefix (exp-sign-string) exponent exp-suffix sci-suffix
      \endverbatim
      
      \b Examples

      The code
      \code 
      format_float fd;
      fd.latex_mode();
      cout << fd.convert(-sqrt(2.0)*1.0e-5) << endl;
      cout << fd.convert(sqrt(2.0)*1.0e-2) << endl;
      cout << fd.convert(-sqrt(2.0)*1.0e-1) << endl;
      \endcode
      outputs
      \verbatim
      $-$1.4142 $\\times 10^{-5}$
      0.014142
      $-$0.14142
      \endverbatim
      and the code
      \code
      format_float fd;
      fd.html_mode();
      fd.set_sig_figs(7);
      fd.set_pad_zeros(true);
      cout << fd.convert(1.414e-5) << endl;
      cout << fd.convert(1.414e-2) << endl;
      \endcode
      outputs
      \verbatim
      1.414000 &times; 10<sup>-5</sup>
      0.01414000
      \endverbatim

      \future Handle inf's and nan's correctly.
      \future Allow change of string for the "+" sign for the exponent
  */
  class format_float {
    
  protected:

    /// \name Base text settings
    //@{
    /// Prefix (default "")
    std::string prefx;
    /// Suffix (default "")
    std::string suffx;
    /// Sign string (default "-")
    std::string sgn;
    /// Sign string in scientific mode (default "-")
    std::string sci_sgn;
    /// Sign string for exponent in scientific mode (default "-")
    std::string exp_sgn;
    /// Prefix in scientific mode (default "")
    std::string sci_prefx;
    /// Suffix in scientific mode (default "")
    std::string sci_suffx;
    /// Exponent prefix (default "")
    std::string exp_prefx;
    /// Exponent suffix (default "")
    std::string exp_suffx;
    /// Times symbol for scientific mode (default " x ")
    std::string tmes;
    /// String for numbers which are not finite (default "Nan")
    std::string not_finte;
    /// String for zeros (default "0")
    std::string zeros;
    //@}
    
    /// \name Other settings
    //@{
    /// Number of significant figures (default 5)
    size_t sig_fgs;
    /// Number of digits in exponent (default 0 which prints the minimum)
    size_t exp_dgs;
    /// Lower limit for automatic mode (default -2)
    int ex_mn;
    /// Upper limit for automatic mode (default 3)
    int ex_mx;
    /// If true, pad with zeros (default false)
    bool pad_zeros;
    /** \brief If true, show the sign of the exponent when 
	it's positive (default false)
    */
    bool show_exp_sgn;
    /// The decimal point (default <tt>'.'</tt>)
    std::string dpt;
    //@}
    
    /** \brief Remove extra zeros and decimal point from mantisaa
     */
    int remove_zeros_dpt(std::string &s);

  public:

    format_float();

    /// \name Basic usage
    //@{
    /** \brief Set HTML mode
	
	This function is equivalent to the settings:
	\code
	set_prefix("");
	set_sign("-");
	set_suffix("");
	set_sci_prefix("");
	set_times(" &times; ");
	set_exp_prefix("10<sup>");
	set_exp_sign("-");
	set_sci_sign("-");
	set_exp_suffix("</sup>");
	set_sci_suffix("");
	set_not_finite("Nan");
	set_zero("0");
	set_exp_digits(0);
	set_show_exp_sign(false);
	\endcode
    */
    void html_mode();

    /** \brief Set Latex mode

	This function is equivalent to the settings:
	\code
	set_prefix("");
	set_sign("$-$");
	set_suffix("");
	set_sci_prefix("");
	set_times(" $\\times ");
	set_exp_prefix("10^{");
	set_exp_sign("-");
	set_sci_sign("$-$");
	set_exp_suffix("}");
	set_sci_suffix("");
	set_not_finite("Nan");
	set_zero("0");
	set_exp_digits(0);
	set_show_exp_sign(false);
	\endcode
	
	\note This setting assumes that the user is not 
	in LaTeX's "math mode" already.
    */
    void latex_mode();
    
    /** \brief C-like mode

	This reproduces the default settings of \c cout in automatic
	mode. Obviously it is faster to use iostreams than to format
	numbers with this class. Nevertheless, this mode is very
	useful for testing to ensure that this class processes the
	numbers correctly at the requested precision.

	This function is equivalent to the settings:
	\code
	set_prefix("");
	set_sign("-");
	set_suffix("");
	set_sci_prefix("");
	set_times("e");
	set_exp_prefix("");
	set_exp_sign("-");
	set_sci_sign("-");
	set_exp_suffix("");
	set_sci_suffix("");
	set_not_finite("NaN");
	set_zero("0");
	set_exp_limits(-4,5);
	set_exp_digits(2);
	set_show_exp_sign(true);
	\endcode
    */
    void c_mode();

    /// Convert a floating point number to a string
    std::string convert(double x, bool debug=false);
    //@}
    
    /** \name Set text settings

	These are modified by the functions html_mode() and latex_mode()
    */
    //@{
    /// set prefix
    void set_prefix(std::string prefix) {
      prefx=prefix;
      return;
    }

    /// Set suffix
    void set_suffix(std::string suffix) {
      suffx=suffix;
      return;
    }

    /// Set prefix for scientific notation
    void set_sci_prefix(std::string sci_prefix) {
      sci_prefix=sci_prefx;
      return;
    }

    /// Set suffix for scientific notation
    void set_sci_suffix(std::string sci_suffix) {
      sci_suffx=sci_suffix;
      return;
    }

    /// Set prefix for exponent
    void set_exp_prefix(std::string exp_prefix) {
      exp_prefx=exp_prefix;
      return;
    }

    /// Set suffix for exponent
    void set_exp_suffix(std::string exp_suffix) {
      exp_suffx=exp_suffix;
      return;
    }

    /// Set sign
    void set_sign(std::string sign) {
      sgn=sign;
      return;
    }

    /// Set sign for exponent
    void set_exp_sign(std::string exp_sign) {
      exp_sgn=exp_sign;
      return;
    }

    /// Set policy for showing positive exponent sign
    void set_show_exp_sign(bool b) {
      show_exp_sgn=b;
      return;
    }

    /// Set sign for scientific notation
    void set_sci_sign(std::string sci_sign) {
      sci_sgn=sci_sign;
      return;
    }

    /// Set times
    void set_times(std::string times) {
      tmes=times;
      return;
    }

    /// Set zero
    void set_zero(std::string zero) {
      zeros=zero;
      return;
    }

    /// Set string for numbers which are not finite
    void set_not_finite(std::string not_finite) {
      not_finte=not_finite;
      return;
    }
    //@}

    /** \name Set other settings

	These are not modified by the functions html_mode() and latex_mode()
    */
    //@{
    /// Set the exponent limits
    void set_exp_limits(int min, int max) {
      ex_mn=min;
      ex_mx=max;
      return;
    }
    
    /** \brief Set the number of significant figures (argument has
	maximum of 15 and cannot be zero)
     */
    void set_sig_figs(size_t sig_figs) {
      if (sig_fgs==0 || sig_figs>15) {
	O2SCL_ERR2("Argument must be less than or equal to 15",
		   "in format_float::set_sig_figs().",exc_einval);
      }
      sig_fgs=sig_figs;
      return;
    }

    /// Set pad zeros
    void set_pad_zeros(bool pad) {
      pad_zeros=pad;
      return;
    }

    /// Set decimal point
    void set_dec_point(std::string dec_point) {
      dpt=dec_point;
      return;
    }

    /// Set minimum number of digits in the exponent
    void set_exp_digits(size_t d) {
      exp_dgs=d;
      return;
    }
    //@}

    /** \name Get text settings

	These are modified by the functions html_mode() and latex_mode()
    */
    //@{
    /// Get prefix
    std::string get_prefix() {
      return prefx;
    }
    
    /// Get suffix
    std::string get_suffix() {
      return suffx;
    }
    
    /// Get prefix for scientific notation
    std::string get_sci_prefix() {
      return sci_prefx;
    }
    
    /// Get suffix for scientific notation
    std::string get_sci_suffix() {
      return sci_suffx;
    }
    
    /// Get prefix for exponent
    std::string get_exp_prefix() {
      return exp_prefx;
    }
    
    /// Get suffix for exponent
    std::string get_exp_suffix() {
      return exp_suffx;
    }
    
    /// Get sign
    std::string get_sign() {
      return sgn;
    }
    
    /// Get sign for exponent
    std::string get_exp_sign() {
      return exp_sgn;
    }
    
    /// Get sign for scientific notation
    std::string get_sci_sign() {
      return sci_sgn;
    }
    
    /// Get times
    std::string get_times() {
      return tmes;
    }
    
    /// Get zero
    std::string get_zero() {
      return zeros;
    }
    
    /// Get string for numbers which are not finite
    std::string get_not_finite() {
      return not_finte;
    }
    //@}

    /** \name Get other settings

	These are not modified by the functions html_mode() and latex_mode()
    */
    //@{
    /// Get minimum exponent
    int get_exp_min() {
      return ex_mn;
    }

    /// Get maximum exponent
    int get_exp_max() {
      return ex_mx;
    }
    
    /// Get sig_figs
    size_t get_sig_figs() {
      return sig_fgs;
    }
    
    /// Get pad_zeros
    bool get_pad_zeros() {
      return pad_zeros;
    }

    /// Get decimal point
    std::string get_dec_point() {
      return dpt;
    }
    //@}

  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
