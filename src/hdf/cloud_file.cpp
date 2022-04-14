/*
  -------------------------------------------------------------------
  
  Copyright (C) 2016-2022, Andrew W. Steiner
  
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
#include <o2scl/cloud_file.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

cloud_file::cloud_file() {
  allow_wget=true;
  allow_curl=true;
  verbose=1;
  throw_on_fail=true;
  hash_type=sha256;
}

int cloud_file::hdf5_open_hash
(o2scl_hdf::hdf_file &hf, std::string file, 
 std::string url, std::string hash, std::string dir) {
  get_file_hash(file,url,hash,dir);
  if (dir.length()>0) {
    hf.open(dir+"/"+file);
  } else {
    hf.open(file);
  }
  return 0;
}
    
int cloud_file::hdf5_open(o2scl_hdf::hdf_file &hf, std::string file,
			  std::string url, std::string dir) {
  return hdf5_open_hash(hf,file,url,"",dir);
}
  
int cloud_file::get_file(std::string file, std::string url, 
			 std::string dir) {
  return get_file_hash(file,url,"",dir);
}
    
int cloud_file::get_file_hash
(std::string file, std::string url,
 std::string hash, std::string dir) {

#ifndef O2SCL_USE_BOOST_FILESYSTEM

  if (file=="_") {
    if (url.length()<5) {
      if (throw_on_fail) {
	O2SCL_ERR("Could not extract filename.",
		  o2scl::exc_einval);
      } else {
	return o2scl::exc_einval;
      }
    }
    size_t loc=url.length()-1;
    size_t it=0;
    file="";
    while (url[loc]!='/' && it<100 && loc>1) {
      file=url[loc]+file;
      it++;
      loc--;
    }
    if (it>=100 || loc==1) {
      if (throw_on_fail) {
	O2SCL_ERR("Could not extract filename.",
		  o2scl::exc_einval);
      } else {
	return o2scl::exc_einval;
      }
    }
    if (verbose>0) {
      cout << "Extracted filename " << file << " from URL\n  "
	   << url << endl;
    }
  }
  
  // File status object
  struct stat sb;
  // Return value of stat()
  int sret=0;

  // -------------------------------------------------------
  // Main directory section

  // First use wordexp to do tilde expansion
  string dir_old=dir;
  if (dir.length()>0) {
    std::vector<std::string> matches;
    int wret=wordexp_wrapper(dir,matches);
    if (matches.size()>1 || matches.size()==0 || wret!=0) {
      if (verbose>1) {
        cout << "Function cloud_file::get_file_hash() failed to perform "
             << "tilde expansion for directory " << dir << "." << endl;
      }
      return 10;
    }
    dir=matches[0];
    if (verbose>1) {
      cout << "Function wordexp() converted "
           << dir_old << " to " << dir << endl;
    }
  }
  
  // Look for directory 
  bool dir_present=false;
  if (dir.length()>0) {
    if (verbose>1) {
      std::cout << "Using directory " << dir << std::endl;
    }
    sret=stat(dir.c_str(),&sb);
    if (sret==0) {
      dir_present=S_ISDIR(sb.st_mode);
    }
    if (dir_present==false) {
      if (verbose>0) {
	std::cout << "Directory '" << dir
		  << "' not present. Trying to create it." << std::endl;
      }
      // If not found, try to make it with 'mkdir'
      std::string cmd=((std::string)"mkdir -p ")+dir;
      if (verbose>1) {
	std::cout << "Directory specified but not present in filesystem."
		  << std::endl;
	std::cout << "Trying to create with command:\n  "
		  << cmd << std::endl;
      }
      int mret=system(cmd.c_str());
      if (mret!=0) {
	if (verbose>1) {
	  std::cout << "Command to make directory '" << cmd
		    << "' failed." << std::endl;
	}
      } else {
	sret=stat(dir.c_str(),&sb);
	if (sret==0) {
	  dir_present=S_ISDIR(sb.st_mode);
	}
      }
    }
    if (sret!=0 || dir_present==false) {
      // If not found, prompt user for it
      std::cout << "Could not find or create directory '" << dir 
		<< ". Please enter new directory name."
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
  }
  
  // End of directory section
  // -------------------------------------------------------

  // -------------------------------------------------------
  // Start of file section
    
  // Now look for the full data file
  std::string fname;
  if (dir.length()>0) {
    fname=dir+"/"+file;
  } else {
    fname=file;
  }

  // First use wordexp to do tilde expansion
  string fname_old=fname;
  std::vector<std::string> matches;
  int wret=wordexp_wrapper(fname,matches);
  if (matches.size()>1 || matches.size()==0 || wret!=0) {
    return 10;
  }
  fname=matches[0];
  if (verbose>1) {
    cout << "Function wordexp() converted "
	 << fname_old << " to " << fname << endl;
  }
  
  bool file_present=false;
  bool valid_hash=false;
  sret=stat(fname.c_str(),&sb);
  if (sret==0) {
    file_present=S_ISREG(sb.st_mode);
  }

  // If there's no hash, assume it's valid
  if (hash.length()==0) valid_hash=true;

  if (verbose>1) {
    std::cout << "Function cloud_file::get_file_hash(): file_present: "
              << file_present << endl;
  }
  
  // If the file was found, check the hash
  if (sret==0 && file_present && hash.length()>0) {
    std::string cmd;
    if (hash_type==sha256) {
      cmd=((std::string)"openssl dgst -sha256 ")+fname+
	" | awk '{print $2}'";
    } else if (hash_type==md5sum) {
      cmd=((std::string)"md5sum ")+fname+
	" | awk '{print $1}'";
    } else {
      cmd=((std::string)"md5 ")+fname+
	" | awk '{print $4}'";
    }
    if (verbose>1) {
      std::cout << "Checking hash with command:\n  " << cmd
		<< std::endl;
    }
    std::string hash2;
    int pret=o2scl::pipe_cmd_string(cmd,hash2,false);
    if (pret==0) {
      o2scl::remove_whitespace(hash2);
      if (hash2==hash) {
	valid_hash=true;
	if (verbose>1) {
	  std::cout << "Hash valid." << std::endl;
	}
      }
    } else {
      if (verbose>0) {
	std::cout << "Function pipe_cmd_string() failed." << std::endl;
      }
    }
    if (valid_hash==false) {
      if (verbose>1) {
	std::cout << "File hash " << hash2 << " does not match "
		  << hash << "." << std::endl;
      }
    }
  }
      
  if (verbose>1) {
    std::cout << "Function cloud_file::get_file_hash(): valid_hash: "
              << valid_hash << endl;
  }
  
  if (sret!=0 || file_present==false || valid_hash==false) {
    // If it couldn't be found, try to download it
    int ret=1;
    if (allow_curl) {
      std::string cmd;
      if (dir.length()>0) {
        cmd=((std::string)"cd ")+dir+"; curl -o "+
          file+" "+url;
      } else {
        cmd=((std::string)"curl -o ")+file+" "+url;
      }
      if (verbose>0) {
	std::cout << "Trying curl command:\n  "
		  << cmd << std::endl;
      }
      ret=system(cmd.c_str());
    }
    if (allow_wget && ret!=0) {
      std::string cmd;
      if (dir.length()>0) {
        cmd=((std::string)"cd ")+dir+"; wget -O "+file+" "+url;
      } else {
        cmd=((std::string)"wget -O ")+file+" "+url;
      }
      if (verbose>0) {
	std::cout << "File did not exist or read failed or invalid hash."
		  << std::endl;
	std::cout << "Trying wget command:\n  "
		  << cmd << std::endl;
      }
      ret=system(cmd.c_str());
    }
    
    // Check to see if download command succeeded
    if (ret!=0) {
      if (verbose>1) {
        std::cout << "Function cloud_file::get_file_hash() received "
                  << "non-zero return value from system command." << endl;
      }
      if (throw_on_fail) {
	O2SCL_ERR("Failed to download file.",o2scl::exc_efilenotfound);
      } else {
	return o2scl::exc_efilenotfound;
      }
    }
    
    // Check to see if the file is there
    if (verbose>1) {
      std::cout << "Function cloud_file::get_file_hash() checking "
                << "that file exists." << endl;
    }
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

  // -------------------------------------------------------
  // Check hash if specified by the caller and it was originally
  // invalid
    
  if (hash.length()>0 && valid_hash==false) {
    std::string cmd;
    if (hash_type==sha256) {
      cmd=((std::string)"openssl dgst -sha256 ")+fname+
	" | awk '{print $2}'";
    } else if (hash_type==md5sum) {
      cmd=((std::string)"md5sum ")+fname+
	" | awk '{print $1}'";
    } else {
      cmd=((std::string)"md5 ")+fname+
	" | awk '{print $4}'";
    }
    if (verbose>1) {
      std::cout << "Checking hash with command:\n  " << cmd
		<< std::endl;
    }
    std::string hash2;
    int pret=o2scl::pipe_cmd_string(cmd,hash2,false);
    if (pret==0) {
      o2scl::remove_whitespace(hash2);
      if (hash2!=hash) {
	if (throw_on_fail) {
	  O2SCL_ERR("Invalid hash after download in cloud_file. Wrong URL?",
		    o2scl::exc_efailed);
	} else {
	  return o2scl::exc_efilenotfound;
	}
      }
    } else {
      if (throw_on_fail) {
	O2SCL_ERR("Function pipe_cmd_string() failed.",
		  o2scl::exc_efailed);
      } else {
	return o2scl::exc_efailed;
      }
    }
  }
  
  // -------------------------------------------------------
  // Output full filename
    
  if (verbose>0) {
    std::cout << "Function cloud_file::get_file_hash() succeeded "
              << "to obtain file named '" << fname
	      << "'." << std::endl;
  }

#else

  // AWS: I'm debugging the boost filesystem functions here.
  std::cout << "Here0 " << dir << std::endl;
  boost::filesystem::path p(dir.c_str());
  std::cout << "Here1." << std::endl;
  std::cout << boost::filesystem::exists(p) << std::endl;
  std::cout << "Here2." << std::endl;
  std::cout << boost::filesystem::is_directory(p) << std::endl;
  std::cout << "Here3." << std::endl;
  std::cout << boost::filesystem::is_regular_file(p) << std::endl;
  std::cout << "Here4." << std::endl;
  exit(-1);
      
#endif      
      
  return o2scl::success;
}

