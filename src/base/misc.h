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
#ifndef O2SCL_MISC_H
#define O2SCL_MISC_H
/** \file misc.h
    \brief Miscellaneous functions
*/

#include <cstdlib>
#include <iostream>
#include <string>
// For stringstream for count_words()
#include <sstream>
#include <vector>
// For std::isinf and std::isnan in C++11
#include <cmath>
// For the vec_index class below
#include <map>
#include <initializer_list>

#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /// \name Functions from misc.h
  //@{
  /** \brief Calculate a Fermi-Dirac distribution function safely
      
      \f$ \left[1+\exp\left(E/T-\mu/T\right)\right]^{-1} \f$ 
      
      This calculates a Fermi-Dirac distribution function guaranteeing
      that numbers larger than \f$ \exp(\mathrm{limit}) \f$ and
      smaller than \f$ \exp(-\mathrm{limit}) \f$ will be avoided. The
      default value of <tt>limit=40</tt> ensures accuracy to within 1
      part in \f$ 10^{17} \f$ compared to the maximum of the
      distribution (which is unity).
      
      Note that this function may return Inf or NAN if \c limit is too 
      large, depending on the machine precision.
  */
  double fermi_function(double E, double mu, double T, double limit=40.0);

  /** \brief Store the first line from the output of the shell
      command \c cmd up to \c nmax characters in \c result
      
      \note This function uses popen() and may not work properly on
      some systems. If HAVE_POPEN was not defined during O2scl's
      compilation, then this function will throw an exception (or if
      \c err_on_fail is false, it will return a nonzero value).
  */
  int pipe_cmd_string(std::string cmd, std::string &result,
		      bool err_on_fail=true, int nmax=80);

  /** \brief Return the first line from the output of the shell
      command \c cmd up to \c nmax characters

      This function always throws exceptions if it fails.
  */
  std::string pipe_cmd_string(std::string cmd, int nmax=80);
  
  /** \brief Reformat the columns for output of width \c size 

      Given a string array \c in_cols of size \c nin, screenify()
      reformats the array into columns creating a new string array \c
      out_cols.
      
      For example, for an array of 10 strings 
      \verbatim
      test1
      test_of_string2
      test_of_string3
      test_of_string4
      test5
      test_of_string6
      test_of_string7
      test_of_string8
      test_of_string9
      test_of_string10
      \endverbatim
      screenify() will create an array of 3 new strings:
      \verbatim
      test1            test_of_string4  test_of_string7  test_of_string10
      test_of_string2  test5            test_of_string8
      test_of_string3  test_of_string6  test_of_string9
      \endverbatim
      
      If the value of \c max_size is less than the length of the
      longest input string (plus one for a space character), then the
      output strings may have a larger length than \c max_size.
  */
  template<class string_arr_t>
    void screenify(size_t nin, const string_arr_t &in_cols, 
		   std::vector<std::string> &out_cols,
		   size_t max_size=80) {

    if (nin==0) {
      O2SCL_ERR("No strings specified in screenify().",exc_efailed);
    }
        
    size_t i,j,lmax,itemp;
    std::string *in_spaces=new std::string[nin];
    
    // Determine size of largest string
    lmax=0;
    for(i=0;i<nin;i++) {
      if (lmax<in_cols[i].size()) lmax=in_cols[i].size();
    }

    // Pad with spaces
    for(i=0;i<nin;i++) {
      itemp=in_cols[i].size();
      in_spaces[i]=in_cols[i];
      for(j=0;j<lmax+1-itemp;j++) {
	in_spaces[i]+=' ';
      }
    }

    // Determine number of rows and columns
    size_t row, col;
    col=max_size/(lmax+1);
    if (col==0) col=1;
    if (nin/col*col==nin) row=nin/col;
    else row=nin/col+1;

    // Create outc
    out_cols.reserve(row);
    for(i=0;i<row;i++) {
      out_cols.push_back("");
      for(j=0;j<col;j++) {
	if (i+j*row<nin) {
	  out_cols[i]+=in_spaces[i+j*row];
	}
      }
    }

    delete[] in_spaces;

    return;
  }

  /** \brief Count the number of words in the string \c str 
   
      Words are defined as groups of characters separated by
      whitespace, where whitespace is any combination of adjacent
      spaces, tabs, carriage returns, etc. On most systems, whitespace
      is usually defined as any character corresponding to the
      integers 9 (horizontal tab), 10 (line feed), 11 (vertical tab),
      12 (form feed), 13 (carriage return), and 32 (space bar). The
      test program \c misc_ts enumerates the characters between 0 and
      255 (inclusive) that count as whitespace for this purpose.

      \future Make consistent with split_string().
  */
  size_t count_words(std::string str);
  
  /** \brief Remove all whitespace from the string \c s

      This function removes all characters in \c s which correspond to
      the integer values 9, 10, 11, 12, 13, or 32.
  */
  void remove_whitespace(std::string &s);

  /** \brief Take a string of binary quads and compress them to 
      hexadecimal digits

      This function proceeds from left to right, ignoring parts of the
      string that do not consist of squences of four '1's or '0's.
  */
  std::string binary_to_hex(std::string s);

  /** \brief Convert RGB to HSV color
      
      Taken from Nathan Schaller's webpage at
      http://www.cs.rit.edu/~ncs/color/t_convert.html 
      
      The inputs should be in the ranges \f$ h \in [0,360] \f$, \f$ s
      \in [0,1] \f$, and \f$ v \in [0,1] \f$. The output values \c r,
      \c g, and \c b are \f$ \in [0,1] \f$.
      
      If s == 0, then h = -1 (undefined)
  */
  void RGBtoHSV(double r, double g, double b, 
		double &h, double &s, double &v);
  
  /** \brief Convert RGB to HSV color
      
      Taken from Nathan Schaller's webpage at
      http://www.cs.rit.edu/~ncs/color/t_convert.html 
      
      The inputs should be in the ranges \f$ h \in [0,360] \f$, \f$ s
      \in [0,1] \f$, and \f$ v \in [0,1] \f$. The output values \c r,
      \c g, and \c b are \f$ \in [0,1] \f$.
      
      If s == 0, then h = -1 (undefined)
  */
  void HSVtoRGB(double h, double s, double v, 
		double &r, double &g, double &b);
  //@}
  
  /** \brief Generate number sequence for testing
      
      A class which generates \c tot numbers from -1 to 1, making sure
      to include -1, 1, 0, and numbers near -1, 0 and 1 (so long as \c
      tot is sufficiently large). If gen() is called more than \c tot
      times, it just recycles through the list again.  
      
      This class is used to generate combinations of coefficients for
      testing the polynomial solvers.

      For example, the first 15 numbers generated by
      an object of type gen_test_number<10> are:
      \verbatim
      0  -1.000000e+00
      1  -9.975274e-01
      2  -8.807971e-01
      3  -1.192029e-01
      4  -2.472623e-03
      5  +0.000000e+00
      6  +2.472623e-03
      7  +1.192029e-01
      8  +8.807971e-01
      9  +1.000000e+00
      10 -1.000000e+00
      11 -9.975274e-01
      12 -8.807971e-01
      13 -1.192029e-01
      14 -2.472623e-03
      \endverbatim

      This function is used in <tt>src/other/poly_ts.cpp</tt> which
      tests the polynomial solvers.

      \future Document what happens if \c tot is pathologically small.
  */
  template<size_t tot> class gen_test_number {

#ifndef DOXYGEN_INTERNAL

  protected:

    /// Count number of numbers generated
    int n;

    /** \brief A constant factor for the argument to 
	<tt>tanh()</tt>, equal to \c tot divided by 20.
    */
    double fact;

#endif

  public:

    gen_test_number() {
      n=0;
      fact=((double)tot)/20.0;
    }

    /// Return the next number in the sequence
    double gen() {
      double x, dtot=((double)tot), dn=((double)n);
      if (n==0) {
	x=-1.0;
      } else if (n==tot/2) {
	x=0.0;
      } else if (n==tot-1) {
	x=1.0;
      } else if (n==tot) {
	// Start the sequence over
	x=-1.0;
	n=0;
      } else if (n<((int)tot)/2) {
	// Since we're in the o2scl namespace, we explicitly
	// specify std::tanh() here
	x=(std::tanh((dn-dtot/4.0)/fact)-1.0)/2.0;
      } else {
	x=(std::tanh((dn-0.75*dtot)/fact)+1.0)/2.0;
      }
      n++;
      return x;
    }
  };

  /** \brief Return the x value of the extremum of a quadratic defined by 
      three \f$ (x,y) \f$ pairs

      This function should work for any floating-point data type,
      but will suffer from problems due to lack of precision in
      some cases.
  */
  template<class data_t>
    data_t quadratic_extremum_x(const data_t x1, const data_t x2, 
				const data_t x3, const data_t y1, 
				const data_t y2, const data_t y3) {

    if (x1==x2 || x2==x3 || x1==x3) {
      O2SCL_ERR2("Two abscissae cannot be equal in function ",
		 "quadratic_extremum_x().",exc_einval);
    }
    
    /*
      Start with:
      y1=a x1^2 + b x1 + c
      y2=a x2^2 + b x2 + c
      y3=a x3^2 + b x3 + c
      
      Eliminate 'c':
      (y1-y2)=a(x1^2-x2^2)+b(x1-x2)
      (y3-y2)=a(x3^2-x2^2)+b(x3-x2)
      
      Eliminate 'b':
      (x3-x2)*(y1-y2)=a*(x1^2-x2^2)*(x3-x2)+b*(x1-x2)*(x3-x2)
      (x1-x2)*(y3-y2)=a*(x3^2-x2^2)*(x1-x2)+b*(x3-x2)*(x1-x2)
      
      Alternatively, eliminate 'c' with:
      (y2-y1)=a(x2^2-x1^2)+b(x2-x1)
      (y3-y1)=a(x3^2-x1^2)+b(x3-x1)
      
      Eliminate 'b':
      (x3-x1)*(y2-y1)=a(x2^2-x1^2)*(x3-x1)+b(x2-x1)*(x3-x1)
      (x2-x1)*(y3-y1)=a(x3^2-x1^2)*(x2-x1)+b(x3-x1)*(x2-x1)
    */
    
    data_t a,b,c,den=(x1*x1-x2*x2)*(x3-x2)-(x3*x3-x2*x2)*(x1-x2);
    if (den==0.0) {
      den=(x2*x2-x1*x1)*(x3-x1)-(x3*x3-x1*x1)*(x2-x1);
      a=((x3-x1)*(y2-y1)-(x2-x2)*(y3-y1))/den;
    } else {
      a=((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))/den;
    }
    b=(y1-y2-a*(x1*x1-x2*x2))/(x1-x2);
    c=y2-a*x2*x2-b*x2;
    
    return -b/2.0/a;
  }

  /** \brief Return values related to a quadratic defined by 
      three \f$ (x,y) \f$ pairs
      
      This function provides the three coefficients of the 
      quadratic, \c a, \c b, and \c c, and the denominator
      \c den.

      This function should work for any floating-point data type,
      but will suffer from problems due to lack of precision in
      some cases.
  */
  template<class data_t>
    void quadratic_extremum_y_full(const data_t x1, const data_t x2, 
				     const data_t x3, const data_t y1, 
				     const data_t y2, const data_t y3,
				     const data_t &xmin, const data_t &ymin,
				     const data_t &a, const data_t &b,
				     const data_t &c, const data_t &den) {

    if (x1==x2 || x2==x3 || x1==x3) {
      O2SCL_ERR2("Two abscissae cannot be equal in function ",
		 "quadratic_extremum_y().",exc_einval);
    }
    
    den=(x1*x1-x2*x2)*(x3-x2)-(x3*x3-x2*x2)*(x1-x2);
    if (den==0.0) {
      den=(x2*x2-x1*x1)*(x3-x1)-(x3*x3-x1*x1)*(x2-x1);
      a=((x3-x1)*(y2-y1)-(x2-x2)*(y3-y1))/den;
    } else {
      a=((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))/den;
    }
    b=(y1-y2-a*(x1*x1-x2*x2))/(x1-x2);
    c=y2-a*x2*x2-b*x2;
    xmin=-b/2.0/a;
    ymin=c-b*b/4.0/a;
    return;
  }

  /** \brief Return the y value of the extremum of a quadratic defined by 
      three \f$ (x,y) \f$ pairs

      This function should work for any floating-point data type,
      but will suffer from problems due to lack of precision in
      some cases.
  */
  template<class data_t>
    data_t quadratic_extremum_y(const data_t x1, const data_t x2, 
				const data_t x3, const data_t y1, 
				const data_t y2, const data_t y3) {
    
    if (x1==x2 || x2==x3 || x1==x3) {
      O2SCL_ERR2("Two abscissae cannot be equal in function ",
		 "quadratic_extremum_y().",exc_einval);
    }
    
    data_t a,b,c,den=(x1*x1-x2*x2)*(x3-x2)-(x3*x3-x2*x2)*(x1-x2);
    if (den==0.0) {
      den=(x2*x2-x1*x1)*(x3-x1)-(x3*x3-x1*x1)*(x2-x1);
      a=((x3-x1)*(y2-y1)-(x2-x2)*(y3-y1))/den;
    } else {
      a=((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))/den;
    }
    b=(y1-y2-a*(x1*x1-x2*x2))/(x1-x2);
    c=y2-a*x2*x2-b*x2;
    
    return c-b*b/4/a;    
  }

  /** \brief Return the (x,y) for the extremum of a quadratic defined by 
      three \f$ (x,y) \f$ pairs

      This function should work for any floating-point data type,
      but will suffer from problems due to lack of precision in
      some cases.
  */
  template<class data_t>
    void quadratic_extremum_xy(const data_t x1, const data_t x2, 
			       const data_t x3, const data_t y1, 
			       const data_t y2, const data_t y3,
			       data_t &x, data_t &y) {
    
    if (x1==x2 || x2==x3 || x1==x3) {
      O2SCL_ERR2("Two abscissae cannot be equal in function ",
		 "quadratic_extremum_xy().",exc_einval);
    }
    
    data_t a,b,c,den=(x1*x1-x2*x2)*(x3-x2)-(x3*x3-x2*x2)*(x1-x2);
    if (den==0.0) {
      den=(x2*x2-x1*x1)*(x3-x1)-(x3*x3-x1*x1)*(x2-x1);
      a=((x3-x1)*(y2-y1)-(x2-x2)*(y3-y1))/den;
    } else {
      a=((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))/den;
    }
    b=(y1-y2-a*(x1*x1-x2*x2))/(x1-x2);
    c=y2-a*x2*x2-b*x2;
    
    x=-b/2/a;
    y=c-b*b/4/a;    

    return;
  }

  /** \brief Return the (x,y) for the extremum of a quadratic defined by 
      three \f$ (x,y) \f$ pairs

      This function should work for any floating-point data type,
      but will suffer from problems due to lack of precision in
      some cases.
  */
  template<class data_t>
    void quadratic_extremum_coeffs(const data_t x1, const data_t x2, 
				   const data_t x3, const data_t y1, 
				   const data_t y2, const data_t y3,
				   data_t &a, data_t &b, data_t &c) {
    
    if (x1==x2 || x2==x3 || x1==x3) {
      O2SCL_ERR2("Two abscissae cannot be equal in function ",
		 "quadratic_extremum_coeffs().",exc_einval);
    }
    
    data_t den=(x1*x1-x2*x2)*(x3-x2)-(x3*x3-x2*x2)*(x1-x2);
    if (den==0.0) {
      den=(x2*x2-x1*x1)*(x3-x1)-(x3*x3-x1*x1)*(x2-x1);
      a=((x3-x1)*(y2-y1)-(x2-x2)*(y3-y1))/den;
    } else {
      a=((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))/den;
    }
    b=(y1-y2-a*(x1*x1-x2*x2))/(x1-x2);
    c=y2-a*x2*x2-b*x2;

    return;
  }

#ifndef O2SCL_OLDER_COMPILER
  
  /** \brief A class to assign string labels to array indices
   */
  class vec_index {

  protected:

    /// The map version for string lookup
    std::map<std::string,size_t,std::greater<std::string> > tmap;
    /// The vector version for size_t lookup
    std::vector<std::string> tvec;

  public:

    /// Create an empty assignment
    vec_index() {}
    
    /// Create an assignment based on the strings in \c list
    vec_index(std::vector<std::string> &list) {
      for(size_t i=0;i<list.size();i++) {
	tmap.insert(std::make_pair(list[i],i));
	tvec.push_back(list[i]);
      }
    }
    
    /// Create an assignment based on the strings in \c list
    vec_index(std::initializer_list<std::string> list) {
      size_t ix=0;
      for(std::initializer_list<std::string>::iterator it=list.begin();
	  it!=list.end();it++) {
	tmap.insert(std::make_pair(*it,ix));
	tvec.push_back(*it);
	ix++;
      }
    }
    
    /// Return the string of index \c i
    std::string operator()(size_t i) {
      return tvec[i];
    }

    /// Return the index of string \c s
    size_t operator()(std::string s) {
      std::map<std::string,size_t,std::greater<std::string> >::iterator it;
      it=tmap.find(s);
      if (it==tmap.end()) {
	std::string str=((std::string)"Failed to find '")+s+
	  "' in vec_index::operator().";
	O2SCL_ERR(str.c_str(),o2scl::exc_efailed);
      }
      return it->second;
    }

    /// Return the string of index \c i
    std::string operator[](size_t i) {
      return tvec[i];
    }

    /// Return the index of string \c s
    size_t operator[](std::string s) {
      std::map<std::string,size_t,std::greater<std::string> >::iterator it;
      it=tmap.find(s);
      if (it==tmap.end()) {
	std::string str=((std::string)"Failed to find '")+s+
	  "' in vec_index::operator[].";
	O2SCL_ERR(str.c_str(),o2scl::exc_efailed);
      }
      return it->second;
    }

    /// Add string \c s and assign it the next index
    void append(std::string s) {
      tmap.insert(std::make_pair(s,tvec.size()));
      tvec.push_back(s);
      return;
    }
    
    /// Add a list of strings
    void append(std::initializer_list<std::string> list) {
      size_t ix=tvec.size();
      for(std::initializer_list<std::string>::iterator it=list.begin();
	  it!=list.end();it++) {
	tmap.insert(std::make_pair(*it,ix));
	tvec.push_back(*it);
	ix++;
      }
    }
  };

#endif

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif

