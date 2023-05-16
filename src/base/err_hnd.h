/*
  ───────────────────────────────────────────────────────────────────
  
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

  ───────────────────────────────────────────────────────────────────
*/
#ifndef O2SCL_ERR_HND_H
#define O2SCL_ERR_HND_H

/** \file err_hnd.h
    \brief Error handling classes, \ref o2scl::err_hnd_type and 
    \ref o2scl::err_hnd_gsl

    See also \ref exception.h .
*/

#include <iostream>
#include <string>

namespace o2scl {

  /** \brief The integer error definitions in src/base/err_hnd.h
      
      The errors associated with the integers between -2 and 32
      are based on GSL, the rest are specific to \o2 . 
  */
  enum {
    /// Success
    success=0,
    /// Failure
    gsl_failure=-1,
    /// iteration has not converged
    gsl_continue=-2,  
    /// input domain error, e.g sqrt(-1)
    exc_edom=1,   
    /// output range error, e.g. exp(1e100)
    exc_erange=2,   
    /// invalid pointer
    exc_efault=3,   
    /// invalid argument supplied by user
    exc_einval=4,   
    /// generic failure
    exc_efailed=5,   
    /// factorization failed
    exc_efactor=6,   
    /// sanity check failed - shouldn't happen 
    exc_esanity=7,   
    /// malloc failed
    exc_enomem=8,   
    /// problem with user-supplied function
    exc_ebadfunc=9,   
    /// iterative process is out of control
    exc_erunaway=10,  
    /// exceeded max number of iterations
    exc_emaxiter=11,  
    /// tried to divide by zero
    exc_ezerodiv=12,  
    /// user specified an invalid tolerance 
    exc_ebadtol=13,  
    /// failed to reach the specified tolerance 
    exc_etol=14,  
    /// underflow
    exc_eundrflw=15,  
    /// overflow
    exc_eovrflw=16,  
    /// loss of accuracy
    exc_eloss=17,  
    /// failed because of roundoff error 
    exc_eround=18,  
    /// matrix, vector lengths are not conformant 
    exc_ebadlen=19,  
    /// matrix not square
    exc_enotsqr=20,  
    /// apparent singularity detected
    exc_esing=21,  
    /// integral or series is divergent 
    exc_ediverge=22,  
    /// requested feature is not supported by the hardware 
    exc_eunsup=23,  
    /// requested feature not (yet) implemented 
    exc_eunimpl=24,  
    /// cache limit exceeded 
    exc_ecache=25,  
    /// \table limit exceeded 
    exc_etable=26,  
    /// iteration is not making progress toward solution 
    exc_enoprog=27,  
    /// \jacobian evaluations are not improving the solution 
    exc_enoprogj=28,  
    /// cannot reach the specified tolerance in f 
    exc_etolf=29,  
    /// cannot reach the specified tolerance in x 
    exc_etolx=30,  
    /// cannot reach the specified tolerance in \gradient 
    exc_etolg=31,  
    /// end of file 
    exc_eof=32,  
    /// Generic "not found" result
    exc_enotfound=33,
    /// Incorrect type for memory object
    exc_ememtype=34,
    /// File not found
    exc_efilenotfound=35,
    /// Invalid index for array or matrix
    exc_eindex=36,
    /// Outside constraint region
    exc_outsidecons=37
  };

  // Forward declaration
  class err_hnd_type;

  /** \brief The global error handler pointer

      This is set by the \ref def_err_hnd constructor to 
      point to that object.
   */      
  extern err_hnd_type *err_hnd;

  /** \brief Class defining an error handler [abstract base]

      A global object of this type is defined, \ref err_hnd .

      \future There may be an issue associated with the string
      manipulations causing errors in the error handler.
   */
  class err_hnd_type {
    
  public:

    err_hnd_type() {}

    virtual ~err_hnd_type() {}

    /** \brief Set an error 
    
	This is separate from set(), since the gsl error handler
	needs to be a static function.
    */
    static void gsl_hnd(const char *reason, const char *file, 
			int line, int lerrno) {
      err_hnd->set(reason,file,line,lerrno);
    }

    /// Set an error 
    virtual void set(const char *reason, const char *file, 
		     int line, int lerrno)=0;
    
    /// Get the last error
    virtual void get(const char *&reason, const char *&file,
		     int &line, int &lerrno)=0;

    /// Return the last error number
    virtual int get_errno() const=0;

    /// Return the line number of the last error
    virtual int get_line() const=0;

    /// Return the reason for the last error
    virtual const char *get_reason() const=0;

    /// Return the file name of the last error
    virtual const char *get_file() const=0;

    /// Return a string summarizing the last error
    virtual const char *get_str()=0;

    /// Remove last error information
    virtual void reset()=0;

    /// Return type
    virtual const char *type() const=0;

  };
  
  /** \brief The error handler 
      
      An error handler for use in \o2 which replaces the GSL error handler
      
      Note that the string arguments to set() can refer to temporary
      storage, since they are copied when the function is called and
      an error is set.
  */
  class err_hnd_gsl : public err_hnd_type {

  public:  

    err_hnd_gsl();

    virtual ~err_hnd_gsl() {}

    /// Set an error 
    virtual void set(const char *reason, const char *file, 
		     int line, int lerrno);

    /// Get the last error
    virtual void get(const char *&reason, const char *&file,
		     int &line, int &lerrno);

    /// Return the last error number
    virtual int get_errno() const;

    /// Return the line number of the last error
    virtual int get_line() const;

    /// Return the reason for the last error
    virtual const char *get_reason() const;

    /// Return the file name of the last error
    virtual const char *get_file() const;

    /// Return a string summarizing the last error
    virtual const char *get_str();

    /// Remove last error information
    virtual void reset();

    /// Number of characters from filename to print (default 28)
    size_t fname_size;

    /// Return type ("err_hnd_gsl")
    virtual const char *type() const { return "err_hnd_gsl"; }

  protected:

    /// Convert an error number to a string
    std::string errno_to_string(int errnox);

    /// The maximum size of error explanations
    static const int rsize=300;
    /// The maximum size of error explanations with the line and file info
    static const int fsize=400;

    /// The error number
    int a_errno;
    /// The line number
    int a_line;

    /// The filename
    char *a_file;
    /// The error explanation
    char a_reason[rsize];
    /// A full string with explanation and line and file info
    char fullstr[fsize];

  };

  /** \brief An alternate GSL-like error handler
   */      
  extern err_hnd_gsl alt_err_hnd;

  /** \brief Set an error with message \c d and code \c n
   */
#define O2SCL_ERR(d,n) o2scl::set_err_fn(d,__FILE__,__LINE__,n);
  
  /** \brief Set a "convergence" error
   */
#define O2SCL_CONV(d,n,b) {if (b) o2scl::set_err_fn(d,__FILE__,__LINE__,n);}
  
  /** \brief Set an error, two-string version
   */
#define O2SCL_ERR2(d,d2,n) o2scl::set_err_fn((std::string(d)+d2).c_str(), \
					     __FILE__,__LINE__,n);
  
  /** \brief Set an error, three-string version
   */
#define O2SCL_ERR3(d,d2,d3,n) o2scl::set_err_fn(\
  (std::string(d)+d2+d3).c_str(),__FILE__,__LINE__,n);
  
  /** \brief Set a "convergence" error, two-string version
   */
#define O2SCL_CONV2(d,d2,n,b) {if (b)					\
      o2scl::set_err_fn((std::string(d)+d2).c_str(),			\
			__FILE__,__LINE__,n);}
  
  /** \brief Set a "convergence" error and return the error value
   */
#define O2SCL_CONV_RET(d,n,b)						\
  do { if (!b) { return n; } else {					\
      o2scl::set_err_fn(d,__FILE__,__LINE__,n); return n; } } while (0)
  
  /** \brief Set an error and return the error value, two-string version
   */
#define O2SCL_CONV2_RET(d,d2,n,b)					\
  do { if (!b) { return n; } else {					\
      o2scl::set_err_fn((std::string(d)+d2).c_str(),			\
			__FILE__,__LINE__,n); return n; } } while (0)
  
  /// \name The error handler function in src/base/err_hnd.h
  //@{
  /** \brief Call the error handler
   */
  inline void set_err_fn(const char *desc, const char *file, int line,
			 int errnum) {
    err_hnd->set(desc,file,line,errnum);
    return;
  }
  //@}
  
  /// \name The error update function in src/base/err_hnd.h
  //@{
  /** \brief Update an error value \c err with the value in \c ret

      If \c ret is zero, this sets \c ret to the value \c err, and 
      if \c ret is nonzero this function does nothing.
   */
  inline void error_update(int &ret, int err) { if (ret==0) ret=err; }
  //@}

#ifdef O2SCL_NEVER_DEFINED
  
  /** \brief A version of \c assert, i.e. exit if the error value is
      non-zero and do nothing otherwise
      
      \future Make this consistent with assert() using NDEBUG?
  */
#define O2SCL_ASSERT(ev)						\
  do { if (ev!=0) { std::cout << "O2scl: Macro err_assert() causing exit" \
			      << " from error " << ev << " at "		\
			      << __LINE__ << " in file:\n "		\
			      << __FILE__ << std::endl;			\
      std::cout << "Error handler string:\n " << err_hnd->get_str()	\
		<< std::endl; exit(ev); } } while (0)
  
  /** \brief A version of \c assert for bool types. Exit if the argument
      is false
  */
#define O2SCL_BOOL_ASSERT(ev,str)					\
  do { if (ev==false) {							\
      std::cout << "O2scl: Macro bool_assert() causing exit at line "	\
		<< __LINE__ << " in file:\n "				\
		<< __FILE__ << std::endl;				\
      std::cout << "Given string: " << str				\
		<< std::endl; exit(-1); } } while (0)

#endif
  
}

extern "C" {

  void o2scl_python_prep();
  
}

#endif
