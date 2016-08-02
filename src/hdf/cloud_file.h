/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016, Andrew W. Steiner
  
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
/** \file cloud_file.h
    \brief File for definition of \ref o2scl::cloud_file
*/
#ifndef O2SCL_CLOUD_FILE_H
#define O2SCL_CLOUD_FILE_H

#include <cmath>
#include <o2scl/err_hnd.h>
#include <o2scl/hdf_file.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl_hdf {
#endif

  /** \brief Desc
  */
  class cloud_file {

  public:
    
    /** \brief Desc
     */
    bool allow_wget;
    /** \brief Desc
     */
    bool allow_curl;
    /** \brief Desc
     */
    bool force_subdir;
    /** \brief Desc
     */
    int verbose;
    /** \brief Desc
     */
    bool throw_on_fail;
    /** \brief Desc
     */
    std::string env_var;
    
    cloud_file() {
      alloc_wget=true;
      allow_curl=true;
      force_subdir=true;
      verbose=1;
      throw_on_fail=true;
      env_var="O2SCL_DATA";
    }
    
    /** \brief Desc
     */
    int hdf5_open(hdf_file &hf, std::string file, std::string subdir,
		  std::string url, std::string dir="") {
      get_file(file,subdir,url,dir);
      hf.open(file);
      return 0;
    }
    
    /** \brief Desc
     */
    int hdf5_open_or_create(hdf_file &hf, std::string file, std::string subdir,
			    std::string url, std::string dir="") {
      get_file(file,subdir,url,dir);
      hf.open_or_create(file);
      return 0;
    }
    
    /** \brief Desc
     */
    int get_file(std::string file, std::string subdir, std::string url,
		 std::string dir="") {
      
      if (dir=="") {
	char *dir_ptr=getenv(env_var.c_str());
	if (dir_ptr!=0) {
	  dir=dir_ptr;
	}
      }

      struct stat sb;
      if (dir.length()==0 || stat(dir.c_str(),&sb)!=0 ||
	  !S_ISDIR(sb.st_mode)) {
	std::cout << "No directory specified. Please enter directory."
		  << std::endl;
	std::cin >> dir;
	if (dir.length()==0 || stat(dir.c_str(),&sb)!=0 ||
	    !S_ISDIR(sb.st_mode)) {
	  O2SCL_ERR("Could not find correct directory.",o2scl::exc_efailed);
	}
      }

      std::string full_dir=dir+"/"+subdir;
      bool full_dir_present=true;
      
      if (full_dir.length()==0 || stat(full_dir.c_str(),&sb)!=0 ||
	  !S_ISDIR(sb.st_mode)) {
	if (verbose>0) {
	  std::cout << "Directory did not exist. Trying mkidr."
		    << std::endl;
	}
	std::string cmd=((std::string)"mkdir -p ")+full_dir;
	int ret=system(cmd.c_str());
	if (ret!=0) {
	  O2SCL_ERR("Failed to create directory.",o2scl::exc_efailed);
	}
	full_dir_present=false;
      }
      
      std::string fname=full_dir+"/"+file;
      if (!full_dir_present || stat(fname.c_str(),&sb)!=0 ||
	  !S_ISREG(sb.st_mode)) {
	int ret=1;
	if (allow_curl) {
	  if (verbose>0) {
	    std::cout << "Directory did not exist. Trying curl."
		      << std::endl;
	  }
	  std::string cmd=((std::string)"cd ")+full_dir+"; curl "+url;
	  ret=system(cmd.c_str());
	}
	if (allow_wget && ret!=0) {
	  if (verbose>0) {
	    std::cout << "Directory did not exist. Trying wget."
		      << std::endl;
	  }
	  std::string cmd=((std::string)"cd ")+full_dir+"; wget "+url;
	  ret=system(cmd.c_str());
	}
	if (ret!=0) {
	  O2SCL_ERR("Failed to download file.",o2scl::exc_efailed);
	}
      }
      
      return o2scl::success;
    }
    
  };

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
