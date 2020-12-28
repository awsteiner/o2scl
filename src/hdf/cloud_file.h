/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016-2021, Andrew W. Steiner
  
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

      \future Convert to use boost::filesystem?

      \future Automatically handle compressed files? This turns out to
      be complicated, because does the user specify the compressed
      hash, the decompressed hash, or both? I suppose it's best that
      the user specifies both, and you only download if you need to,
      or uncompress if you need to, but this will require nearly a
      full reworking of the get_file_hash function (or a completely
      new function). For the moment, compression has to be handled
      by the user. 

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
    int hdf5_open_hash(hdf_file &hf, std::string file, std::string url,
		      std::string hash, std::string dir="");

    /** \brief Get file named \c file in directory \c dir 
	from url \c url
    */
    int get_file(std::string file, std::string url,
		 std::string dir="");
    
    /** \brief Get file named \c file in directory \c dir 
	in subdirectory \c subdir from url \c url
	
	This function begins with the directory \c dir. If \c dir is
	not present and cannot be created, the user is prompted for
	the correct data directory. This function then searches for
	file \c file in the directory. If it is found, it is compared
	with the specified hash. If the hash matches, then the full
	filename is returned. If the hash does not match or if the
	file is not found, then this function uses curl (or wget if
	curl was unsuccessful) to download the file from \c url. The
	file is then compared with the hash again, and the full
	filename is returned if the hash matches. Otherwise the error
	handler is called.
    */
    int get_file_hash(std::string file, std::string url, std::string hash="",
		      std::string dir="");
    
  };
  
#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
