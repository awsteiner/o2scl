/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2021, Andrew W. Steiner
  
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

#include <o2scl/err_hnd.h>
#include <o2scl/string_conv.h>
// For exit()
#include <cstdlib>
// For strcpy(), etc.
#include <cstring>

using namespace std;
using namespace o2scl;

err_hnd_type *o2scl::err_hnd;
err_hnd_gsl o2scl::alt_err_hnd;

err_hnd_gsl::err_hnd_gsl() {
  
  a_errno=0;
  a_line=0;
  a_file=0;
  a_reason[0]='\0';

  fname_size=28;
}

void err_hnd_gsl::reset() {

  a_errno=0;
  a_line=0;
  a_file=0;
  a_reason[0]='\0';

  return;
}

void err_hnd_gsl::set(const char *reason, const char *file, 
		      int line, int lerrno) {
  /*
    a_errno=lerrno;
    a_file=(char *)file;
    a_line=line;
    strncpy(a_reason,reason,200);
  */
  std::cerr << "O2scl error: " << reason << "\n  at line "
	    << line << " in file " << file << "." << std::endl;
  std::cerr << "Exiting with code " << lerrno << std::endl;
  exit(lerrno);
  return;
}

void err_hnd_gsl::get(const char *&reason, const char *&file, 
		      int &line, int &lerrno) {

  lerrno=a_errno;
  file=a_file;
  reason=a_reason;
  line=a_line;

  return;
}

int err_hnd_gsl::get_errno() const {
  return a_errno;
}

int err_hnd_gsl::get_line() const {
  return a_line;
}

const char *err_hnd_gsl::get_reason() const {
  return a_reason;
}

const char *err_hnd_gsl::get_file() const {
  return a_file;
}

const char *err_hnd_gsl::get_str() {

  if (a_errno==0) {

    sprintf(fullstr,"No error occured.");

  } else {
    
    std::string temp="Error ";
    temp+=errno_to_string(a_errno);
    temp+=" in file ";

    // Print out only the last part of the filename to keep it on 
    // only 1 line
    std::string ftemp=a_file;
    if (ftemp.length()>fname_size) {
      ftemp=ftemp.substr(ftemp.length()-fname_size,ftemp.length());
      ftemp=((std::string)"...")+ftemp;
    }
    temp+=ftemp;

    temp+=" at line ";
    temp+=itos(a_line);
    temp+=".\n  ";
    temp+=a_reason;
    size_t i;
    // Use fsize-1 here because we want to save
    // space for the string termination character
    for(i=0;i<temp.size() && i<fsize-1;i++) {
      fullstr[i]=temp[i];
    }
    fullstr[i]='\0';
  }
  return fullstr;
}

string err_hnd_gsl::errno_to_string(int errnox) {
  
  if (errnox==0) return "success";
  if (errnox==-1) return "failure";
  if (errnox==-2) return "continue";
  if (errnox==1) return "edom";
  if (errnox==2) return "erange";
  if (errnox==3) return "efault";
  if (errnox==4) return "einval";
  if (errnox==5) return "efailed";
  if (errnox==6) return "efactor";
  if (errnox==7) return "esanity";
  if (errnox==8) return "enomem";
  if (errnox==9) return "ebadfunc";
  if (errnox==10) return "erunaway";
  if (errnox==11) return "emaxiter";
  if (errnox==12) return "ezerodiv";
  if (errnox==13) return "ebadtol";
  if (errnox==14) return "etol";
  if (errnox==15) return "eundrflw";
  if (errnox==16) return "eovrflw";
  if (errnox==17) return "eloss";
  if (errnox==18) return "eround";
  if (errnox==19) return "ebadlen";
  if (errnox==20) return "enotsqr";
  if (errnox==21) return "esing";
  if (errnox==22) return "ediverge";
  if (errnox==23) return "eunsup";
  if (errnox==24) return "eunimpl";
  if (errnox==25) return "ecache";
  if (errnox==26) return "etable";
  if (errnox==27) return "enoprog";
  if (errnox==28) return "enoprogj";
  if (errnox==29) return "etolf";
  if (errnox==30) return "etolx";
  if (errnox==31) return "etolg";
  if (errnox==32) return "eof";
  if (errnox==33) return "enotfound";
  if (errnox==34) return "ememtype";
  if (errnox==35) return "efilenotfound";
  if (errnox==36) return "eindex";
  if (errnox==37) return "outsidecons";
  return "unknown";
}

void o2scl_python_prep() {
  err_hnd=&alt_err_hnd;
  cout.setf(ios::scientific);
  return;
}
