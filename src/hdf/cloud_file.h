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
    \brief File for definition of \ref o2scl_hdf::cloud_file
*/
#ifndef O2SCL_CLOUD_FILE_H
#define O2SCL_CLOUD_FILE_H

#include <iostream>
// For getenv() 
#include <cstdlib>
// For struct stat and associated functions
#include <sys/stat.h>
#include <o2scl/err_hnd.h>
#include <o2scl/hdf_file.h>
#ifndef O2SCL_LEGACY_IO
#include <boost/filesystem.hpp>
#endif

#ifndef DOXYGEN_NO_O2NS
namespace o2scl_hdf {
#endif

  /** \brief Read a file and download from a URL if necessary
      
      \note This function requires POSIX I/O calls and a system call
      which uses <tt>mkdir -p</tt>, thus will probably only work on
      unix-like systems.

      \todo Convert to use boost::filesystem .

      \warning This class has several potential security issues 
      and should not be used without due care.
  */
  class cloud_file {
    
  public:
  
    /** \brief If true, allow the use of \c wget to download the file
     */
    bool allow_wget;
    /** \brief If true, allow the use of \c curl to download the file
     */
    bool allow_curl;
    /** \brief Verbosity parameter
     */
    int verbose;
    /** \brief If true, throw an exception on failure
     */
    bool throw_on_fail;
    /** \brief The environment variable which stores the directory
     */
    std::string env_var;
  
    cloud_file();

    /** \brief Open an HDF file named \c file in directory \c dir
	in subdirectory \c subdir, downloading from URL \c url if
	necessary
    */
    int hdf5_open(hdf_file &hf, std::string file, std::string subdir,
		  std::string url, std::string dir="");
    
    /** \brief Get file named \c file in directory \c dir 
	in subdirectory \c subdir from url \c url
	
	This function attempts to find a file named \c file in
	subdirectory \c subdir of the data directory \c dir. If \c dir
	is empty, it attempts to set it equal to the value of the
	environment variable \ref env_var. If that environment
	variable is not present, the user is prompted for the correct
	data directory. If the file is not found, then this function
	uses curl (or wget if curl was unsuccessful) to download the
	file from \c url. If this process was successful at finding or
	downloading the file, then the full filename is returned.
	Otherwise, an exception is thrown.
    */
    int get_file(std::string file, std::string subdir, std::string url,
		 std::string &fname, std::string dir="");
  
  };

  // -------------------------------------------------------------------
  // Class cloud_file function definitions
  
  cloud_file::cloud_file() {
    allow_wget=true;
    allow_curl=true;
    verbose=1;
    throw_on_fail=true;
    env_var="";
  }

  int cloud_file::hdf5_open(hdf_file &hf, std::string file,
			    std::string subdir, std::string url,
			    std::string dir) {
    std::string fname;
    get_file(file,subdir,url,fname,dir);
    hf.open(fname);
    return 0;
  }
    
  int cloud_file::get_file(std::string file, std::string subdir,
			   std::string url, std::string &fname,
			   std::string dir) {

    if (dir=="" && env_var.length()>0) {
      char *dir_ptr=getenv(env_var.c_str());
      if (dir_ptr!=0) {
	dir=dir_ptr;
      }
      if (verbose>1) {
	std::cout << "Obtained directory from environment variable '"
		  << env_var << "':\n\t" << dir << std::endl;
      }
    }

#ifdef O2SCL_LEGACY_IO

    // File status object
    struct stat sb;
    // Return value of stat()
    int sret=1;

    // -------------------------------------------------------
    // Main directory section

    // Look for directory 
    bool dir_present=false;
    if (dir.length()>0) {
      sret=stat(dir.c_str(),&sb);
      if (sret==0) {
	dir_present=S_ISDIR(sb.st_mode);
      }
      if (dir_present==false) {
	if (verbose>1) {
	  std::cout << "Directory specified but not present in filesystem."
		    << std::endl;
	  std::cout << "Trying to create with 'mkdir'." << std::endl;
	}
	// If not found, try to make it with 'mkdir'
	std::string cmd=((std::string)"mkdir -p ")+dir;
	int mret=system(cmd.c_str());
	if (mret!=0) {
	  if (verbose>1) {
	    std::cout << "Command to make directory '" << cmd
		      << "' failed." << std::endl;
	  }
	} else {
	  dir_present=true;
	}
      }
    } 
    if (dir.length()==0 || sret!=0 || dir_present==false) {
      // If not found, prompt user for it
      std::cout << "No directory specified or could not create "
		<< "directory. Please enter directory name."
		<< std::endl;
      std::cin >> dir;
      // Check again
      dir_present=false;
      if (dir.length()>0) {
	sret=stat(dir.c_str(),&sb);
	if (sret==0) {
	  dir_present=S_ISDIR(sb.st_mode);
	}
      }
      // If that failed, then give up
      if (dir.length()==0 || sret!=0 || dir_present==false) {
	if (throw_on_fail) {
	  O2SCL_ERR("Could not find correct directory.",
		    o2scl::exc_efilenotfound);
	} else {
	  return o2scl::exc_efilenotfound;
	}
      }
    } 

    if (verbose>1) {
      std::cout << "Directory " << dir << " found." << std::endl;
    }

    // End of main directory section
    // -------------------------------------------------------

    // The full local directory and subdirectory
    std::string full_dir=dir;

    // -------------------------------------------------------
    // Subdirectory section
    
    if (subdir.length()>0) {
      // Subdirectory was specified, so look for it on the filesystem
      std::string full_dir=dir+"/"+subdir;
      if (verbose>1) {
	std::cout << "Set full_dir to: " << full_dir << std::endl;
      }
      bool full_dir_present=false;
      sret=stat(full_dir.c_str(),&sb);
      if (sret==0) {
	full_dir_present=S_ISDIR(sb.st_mode);
      }
      
      if (full_dir.length()==0 || sret!=0 || full_dir_present==false) {
	// If not found, try to make it with 'mkdir'
	std::string cmd=((std::string)"mkdir -p ")+full_dir;
	if (verbose>0) {
	  std::cout << "Directory did not exist. Trying mkdir with "
		    << "command:\n\t" << cmd << std::endl;
	}
	int ret=system(cmd.c_str());
	if (ret!=0) {
	  if (throw_on_fail) {
	    O2SCL_ERR("Failed to create directory.",o2scl::exc_efilenotfound);
	  } else {
	    return o2scl::exc_efilenotfound;
	  }
	}
	// Check again
	full_dir_present=false;
	sret=stat(full_dir.c_str(),&sb);
	if (sret==0) {
	  full_dir_present=S_ISDIR(sb.st_mode);
	}
	// If that failed, then give up
	if (full_dir.length()==0 || sret!=0 || full_dir_present==false) {
	  if (throw_on_fail) {
	    O2SCL_ERR("Could not create full directory.",
		      o2scl::exc_efilenotfound);
	  } else {
	    return o2scl::exc_efilenotfound;
	  }
	}
      } else if (verbose>1) {
	std::cout << "Full directory " << full_dir << " found." << std::endl;
      }
    }

    // End of subdirectory section
    // -------------------------------------------------------
    // Start of file section
    
    // Now look for the full data file
    fname=full_dir+"/"+file;
    bool file_present=false;
    sret=stat(fname.c_str(),&sb);
    if (sret==0) {
      file_present=S_ISREG(sb.st_mode);
    }
      
    if (sret!=0 || file_present==false) {
      // If it couldn't be found, try to download it
      int ret=1;
      if (allow_curl) {
	std::string cmd=((std::string)"cd ")+full_dir+"; curl -o "+file++url;
	if (verbose>0) {
	  std::cout << "File did not exist. Trying curl command:\n\t"
		    << cmd << std::endl;
	}
	ret=system(cmd.c_str());
      }
      if (allow_wget && ret!=0) {
	std::string cmd=((std::string)"cd ")+full_dir+"; wget -O "+file+
	  url;
	if (verbose>0) {
	  std::cout << "File did not exist. Trying wget command:\n\t"
		    << cmd << std::endl;
	}
	ret=system(cmd.c_str());
      }
      // Check to see if download command succeeded
      if (ret!=0) {
	if (throw_on_fail) {
	  O2SCL_ERR("Failed to download file.",o2scl::exc_efilenotfound);
	} else {
	  return o2scl::exc_efilenotfound;
	}
      }
      // Check to see if the file is there
      sret=stat(fname.c_str(),&sb);
      if (sret==0) {
	file_present=S_ISREG(sb.st_mode);
      }
      // If it's still not there, then give up
      if (sret!=0 || file_present==false) {
	if (throw_on_fail) {
	  O2SCL_ERR("Could not find or download file.",
		    o2scl::exc_efilenotfound);
	} else {
	  return o2scl::exc_efilenotfound;
	}
      }
    }

    // Output full filename
    if (verbose>1) {
      std::cout << "Success with file named " << fname << std::endl;
    }

#else

    std::cout << "Here6 " << dir << std::endl;
    boost::filesystem::path p(dir.c_str());
    std::cout << "Here7." << std::endl;
    std::cout << boost::filesystem::exists(p) << std::endl;
    std::cout << "Here7b." << std::endl;
    std::cout << boost::filesystem::is_directory(p) << std::endl;
    std::cout << "Here8." << std::endl;
    std::cout << boost::filesystem::is_regular_file(p) << std::endl;
    std::cout << "Here9." << std::endl;
    exit(-1);
      
#endif      
      
    return o2scl::success;
  }

  // End of cloud_file function definitions
  // -------------------------------------------------------------------
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
