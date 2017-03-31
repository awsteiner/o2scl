/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016-2017, Andrew W. Steiner
  
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

#ifdef O2SCL_USE_BOOST_FILESYSTEM
#include <boost/filesystem.hpp>
#endif

#include <o2scl/err_hnd.h>
#include <o2scl/hdf_file.h>

#ifndef DOXYGEN_NO_O2NS
namespace o2scl_hdf {
#endif

  /** \brief Read a file and download from a URL if necessary
      
      \note This class requires POSIX I/O calls and a system call
      which uses <tt>mkdir -p</tt>, thus will probably only work on
      unix-like systems.

      \note This class uses system calls to <tt>curl</tt> or
      <tt>wget</tt> which must be installed separatley.

      \future Convert to use boost::filesystem .

      \warning This class has several potential security issues 
      and should not be used without due care. 
  */
  class cloud_file {
    
  public:
  
    /** \brief If true, allow the use of \c wget to download the file
	(default true)
     */
    bool allow_wget;
    /** \brief If true, allow the use of \c curl to download the file
	(default true)
     */
    bool allow_curl;
    /** \brief Verbosity parameter (default 1)
     */
    int verbose;
    /** \brief If true, throw an exception on failure (default true)
     */
    bool throw_on_fail;
    /** \brief The environment variable which stores the directory
	(default "")
     */
    std::string env_var;
    /// \name Specify hash type
    //@{
    /// Current hash type (default sha256)
    int hash_type;
    static const int sha256=0;
    static const int md5=1;
    static const int md5sum=2;
    //@}
  
    cloud_file();

    /** \brief Open an HDF file named \c file in directory \c dir
	downloading from URL \c url if necessary
    */
    int hdf5_open(hdf_file &hf, std::string file, 
		  std::string url, std::string dir="");
    
    /** \brief Open an HDF file named \c file in directory \c dir
	with hash \c hash, downloading from URL \c url if
	necessary
    */
    int hdf5_open_hash(hdf_file &hf, std::string file, std::string hash,
		      std::string url, std::string dir="");

    /** \brief Open an HDF file named \c file in directory \c dir
	in subdirectory \c subdir, downloading from URL \c url if
	necessary
    */
    int hdf5_open_subdir(hdf_file &hf, std::string file, std::string subdir,
			 std::string url, std::string dir="");

    /** \brief Open an HDF file named \c file in directory \c dir
	in subdirectory \c subdir with hash \c hash, 
	downloading from URL \c url if necessary
    */
    int hdf5_open_hash_subdir(hdf_file &hf, std::string file, std::string hash,
			     std::string subdir, std::string url,
			     std::string dir="");
			     

    /** \brief Get file named \c file in directory \c dir 
	in subdirectory \c subdir from url \c url
    */
    int get_file(std::string file, std::string url,
		 std::string &fname, std::string dir="");
    
    /** \brief Get file named \c file in directory \c dir 
	in subdirectory \c subdir from url \c url
     */
    int get_file_hash(std::string file, std::string hash, std::string url,
		 std::string &fname, std::string dir="");
    
    /** \brief Get file named \c file in directory \c dir 
	in subdirectory \c subdir from url \c url
     */
    int get_file_subdir(std::string file, std::string subdir, std::string url,
		 std::string &fname, std::string dir="");
    
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
    int get_file_hash_subdir(std::string file, std::string hash,
			    std::string subdir, std::string url,
			    std::string &fname, std::string dir="");

  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
