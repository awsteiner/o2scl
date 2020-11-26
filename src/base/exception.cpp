/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2020, Andrew W. Steiner
  
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// For strncpy(), etc.
#include <cstring>

// For gsl_set_error_handler()
#include <gsl/gsl_errno.h>

#include <o2scl/err_hnd.h>
#include <o2scl/exception.h>

using namespace std;
using namespace o2scl;

err_hnd_cpp o2scl::def_err_hnd;

err_hnd_cpp::err_hnd_cpp() {

#ifdef O2SCL_USE_GSL_HANDLER
  err_hnd=&alt_err_hnd;
#else
  err_hnd=this;
#endif
  
  gsl_set_error_handler(err_hnd->gsl_hnd);
}

void err_hnd_cpp::set(const char *reason, const char *file, 
		      int line, int lerrno) {

#ifdef O2SCL_NO_EXCEPTIONS

  err_hnd_gsl::set(reason,file,line,lerrno);

#else

  a_errno=lerrno;
  a_file=(char *)file;
  a_line=line;
  strncpy(a_reason,reason,200);
      
  if (lerrno==exc_ememtype) {
    throw exc_logic_error(a_reason);
  } else if (lerrno==exc_einval || lerrno==exc_ebadtol ||
	     lerrno==exc_ebadlen || lerrno==exc_enotsqr ||
	     lerrno==exc_eindex) {
    throw exc_invalid_argument(a_reason);
  } else if (lerrno==gsl_failure || lerrno==exc_efailed ||
	     lerrno==exc_esanity || lerrno==exc_eunsup ||
	     lerrno==exc_eunimpl) {
    throw exc_exception();
  } else if (lerrno==exc_edom || lerrno==exc_eundrflw ||
	     lerrno==exc_erange) {
    throw exc_range_error(a_reason);
  } else if (lerrno==exc_ezerodiv || lerrno==exc_eovrflw) {
    throw exc_overflow_error(a_reason);
  } else if (lerrno==exc_eof || lerrno==exc_efilenotfound) {
    throw exc_ios_failure(a_reason);
  } else {
    throw exc_runtime_error(a_reason);
  }

#endif
      
  return;
}

