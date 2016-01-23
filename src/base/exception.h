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
#ifndef O2SCL_EXCEPTION_H
#define O2SCL_EXCEPTION_H

/** \file exception.h
    \brief Error handler class \ref o2scl::err_hnd_cpp and 
    the \o2 exception objects

    See also \ref err_hnd.h .
*/

#include <stdexcept>
#include <iostream>

#include <o2scl/err_hnd.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif
  
  /** \brief Generic exception

      This class derives from <tt>std::exception</tt>.

      The errors which are handled with this exception type are
      - \ref o2scl::gsl_failure (-1) Failure \n
      - \ref o2scl::exc_efailed (5) Generic failure \n
      - \ref o2scl::exc_esanity (7) Sanity check failed \n
      - \ref o2scl::exc_eunsup (23) Requested feature is not 
      supported by the hardware \n
      - \ref o2scl::exc_eunimpl (24) Requested feature not 
      (yet) implemented  
    
  */
  class exc_exception : public std::exception {

  public:
    
    /// Create an exception 
  exc_exception() : std::exception() {
    }

    virtual ~exc_exception() throw() {}

    /// Return the error string
    virtual const char* what() const throw()
    {
      return err_hnd->get_str();
    }
  
  };

  /** \brief Logic error exception

      This class derives from <tt>std::logic_error</tt>.

      The error which is handled with this exception type is
      - \ref o2scl::exc_ememtype (34) Incorrect type for memory object

  */
  class exc_logic_error : public std::logic_error {
    
  public:
    
    /// Create an exception with description provided in \c s
  exc_logic_error(const std::string &s) : std::logic_error(s) {
    }
    
    virtual ~exc_logic_error() throw() {}
    
    /// Return the error string
    virtual const char* what() const throw()
    {
      return err_hnd->get_str();
    }
  
  };

  /** \brief Invalid argument exception

      This class derives from <tt>std::invalid_argument</tt>.

      The errors which are handled with this exception type are
      - \ref o2scl::exc_einval (4) invalid argument supplied by user \n
      - \ref o2scl::exc_ebadtol (13) user specified an invalid tolerance \n
      - \ref o2scl::exc_ebadlen (19) matrix, vector lengths are not conformant \n
      - \ref o2scl::exc_enotsqr (20) matrix not square \n
      - \ref o2scl::exc_eindex (36) Invalid index for array or matrix 
  */
  class exc_invalid_argument : public std::invalid_argument {
    
  public:
    
    /// Create an exception with description provided in \c s
  exc_invalid_argument(const std::string &s) : std::invalid_argument(s) {
    }
    
    virtual ~exc_invalid_argument() throw() {}
    
    /// Return the error string
    virtual const char* what() const throw()
    {
      return err_hnd->get_str();
    }
  
  };

  /** \brief Generic runtime error exception

      This class derives from <tt>std::runtime_error</tt>.

      The errors which are handled with this exception type are
      - \ref o2scl::exc_efault (3) invalid pointer \n
      - \ref o2scl::exc_efactor (6) factorization failed \n
      - \ref o2scl::exc_enomem (8) malloc failed \n
      - \ref o2scl::exc_ebadfunc (9) problem with user-supplied function \n
      - \ref o2scl::exc_erunaway (10) iterative process is out of control \n
      - \ref o2scl::exc_emaxiter (11) exceeded max number of iterations \n
      - \ref o2scl::exc_etol (14) failed to reach the specified tolerance \n
      - \ref o2scl::exc_eloss (17) loss of accuracy \n
      - \ref o2scl::exc_eround (18) failed because of roundoff error \n
      - \ref o2scl::exc_esing (21) apparent singularity detected \n
      - \ref o2scl::exc_ediverge (22) integral or series is divergent \n
      - \ref o2scl::exc_ecache (25) cache limit exceeded \n
      - \ref o2scl::exc_etable (26) table limit exceeded \n
      - \ref o2scl::exc_enoprog (27) iteration is not making progress toward solution \n
      - \ref o2scl::exc_enoprogj (28) jacobian evaluations are not
      improving the solution \n
      - \ref o2scl::exc_etolf (29) cannot reach the specified tolerance in f \n
      - \ref o2scl::exc_etolx (30) cannot reach the specified tolerance in x \n
      - \ref o2scl::exc_etolg (31) cannot reach the specified tolerance in gradient \n
      - \ref o2scl::exc_enotfound (33) Generic "not found" result \n
      - exc_outsidecons (37) Outside constraint region

  */
  class exc_runtime_error : public std::runtime_error {

  public:
    
  exc_runtime_error(const std::string &s) : std::runtime_error(s) {
    }
    
    /// Create an exception with description provided in \c s
    virtual ~exc_runtime_error() throw() {}
    
    /// Return the error string
    virtual const char* what() const throw()
    {
      return err_hnd->get_str();
    }
  
  };
  
  /** \brief Range error runtime exception

      This class derives from <tt>std::range_error</tt>.

      The errors which are handled with this exception type are
      - \ref o2scl::exc_edom (1) input domain error, e.g sqrt(-1) \n
      - \ref o2scl::exc_erange (2) output range error, e.g. exp(1e100) \n
      - \ref o2scl::exc_eundrflw (15) underflow   
  */
  class exc_range_error : public std::range_error {

  public:
    
    /// Create an exception with description provided in \c s
  exc_range_error(const std::string &s) : std::range_error(s) {
    }
    
    virtual ~exc_range_error() throw() {}
    
    /// Return the error string
    virtual const char* what() const throw()
    {
      return err_hnd->get_str();
    }
  
  };
  
  /** \brief Overflow error runtime exception

      This class derives from <tt>std::overflow_error</tt>.

      The errors which are handled with this exception type are
      - \ref o2scl::exc_ezerodiv (12) tried to divide by zero \n
      - \ref o2scl::exc_eovrflw (16) overflow 
  */
  class exc_overflow_error : public std::overflow_error {

  public:
    
    /// Create an exception with description provided in \c s
  exc_overflow_error(const std::string &s) : std::overflow_error(s) {
    }
    
    virtual ~exc_overflow_error() throw() {}
    
    /// Return the error string
    virtual const char* what() const throw()
    {
      return err_hnd->get_str();
    }
  
  };
  
  /** \brief I/O failure error exception

      This class derives from <tt>std::ios::failure</tt>.

      The errors which are handled with this exception type are
      - \ref o2scl::exc_eof=32 end of file \n
      - \ref o2scl::exc_efilenotfound=35 File not found
  */
  class exc_ios_failure : public std::ios::failure {
    
  public:

    /// Create an exception with description provided in \c s
  exc_ios_failure(const std::string &s) : std::ios::failure(s) {
    }
    
    virtual ~exc_ios_failure() throw() {}
    
    /// Return the error string
    virtual const char* what() const throw()
    {
      return err_hnd->get_str();
    }
  
  };
  
  /** \brief Error handler to throw C++ exceptions

      The default error handler, \ref def_err_hnd, is of this type.
   */
  class err_hnd_cpp : public err_hnd_gsl {
    
  public:  

    err_hnd_cpp();

    virtual ~err_hnd_cpp() throw() {}
    
    /// Set an error 
    virtual void set(const char *reason, const char *file, 
		     int line, int lerrno);

    /// Return type ("err_hnd_cpp")
    virtual const char *type() const { return "err_hnd_cpp"; }

  };

  /** \brief The default error handler
   */      
  extern err_hnd_cpp def_err_hnd;

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
