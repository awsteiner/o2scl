/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2024, Andrew W. Steiner

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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <regex>
#include <fnmatch.h>

#include <o2scl/err_hnd.h>
#include <o2scl/hdf_file.h>
#include <o2scl/table.h>
#include <o2scl/hdf_io.h>

using namespace std;
using namespace o2scl;
using namespace o2scl_hdf;

typedef unsigned long long o2_u64_t;

hdf_file::hdf_file() {
  file=0;
  current=0;
  file_open=false;
#ifdef O2SCL_HDF5_COMP
  compr_type=1;
#else
  compr_type=0;
#endif
  write_access=false;
  min_compr_size=40;
}

hdf_file::~hdf_file() {
  if (file_open) {
    H5Fclose(file);
  }
}

int hdf_file::open(std::string fname, bool allow_write, bool err_on_fail) {
      
  H5E_BEGIN_TRY
    {
      if (allow_write) {
	file=H5Fopen(fname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
      } else {
	file=H5Fopen(fname.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
      }
    } 
  H5E_END_TRY
    if (file<0) {
      if (err_on_fail) {
	O2SCL_ERR((((string)"Open file named '")+fname+
		   "' failed in hdf_file::open().").c_str(),exc_efilenotfound);
      }
      return exc_efilenotfound;
    }
  write_access=allow_write;
  file_open=true;
  current=file;
  return success;
}

void hdf_file::open_or_create(std::string fname) {
      
  H5E_BEGIN_TRY
    {
      file=H5Fopen(fname.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    } 
  H5E_END_TRY 
    if(file < 0) {
#ifdef O2SCL_NEVER_DEFINED
      if (parallel) {
	int mpi_size, mpi_rank;
	MPI_Comm comm=MPI_COMM_WORLD;
	MPI_Info info=MPI_INFO_NULL;
	MPI_Comm_size(comm,&mpi_size);
	MPI_Comm_rank(comm,&mpi_rank);  
	hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id,comm,info);
	file=H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,plist_id);
	H5Pclose(plist_id);
      } else {
	file=H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
      }
#else
      file=H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
#endif
    }
  if (file<0) {
    O2SCL_ERR((((string)"Open or create file named '")+fname+
	       "' failed in hdf_file::open_or_create().").c_str(),
	      exc_efilenotfound);
  }
  write_access=true;
  file_open=true;
  current=file;
  return;
}

void hdf_file::close() {
  if (file_open) {
    herr_t status=H5Fclose(file);
    file_open=false;
    current=0;
  } else {
    O2SCL_ERR("No file to close in hdf_file::close().",exc_einval);
  }
  return;
}
    
hid_t hdf_file::get_file_id() {
  if (!file_open) {
    O2SCL_ERR("No file open in hdf_file::get_file_id().",-1);
  }
  return file;
}

hid_t hdf_file::get_current_id() {
  return current;
}

void hdf_file::set_current_id(hid_t cur) {
  current=cur;
  return;
}

void hdf_file::setc(std::string name, char c) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::setc().",exc_efailed);
  }
  
  hid_t dset, space=0;
  bool space_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=1;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);
    dset=H5Dcreate(current,name.c_str(),H5T_STD_I8LE,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    if (dset<0) {
      H5Sclose(space);
      O2SCL_ERR2("Failed to create dataspace in ",
		 "hdf_file::setc().",exc_einval);
    }
    space_alloc=true;

  }

  // Write the data 
  int status=H5Dwrite(dset,H5T_NATIVE_CHAR,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,&c);
      
  int status2=H5Dclose(dset);
  if (space_alloc) status2=H5Sclose(space);

  // Throw exception if write failed
  if (status<0) {
    O2SCL_ERR2("Failed to write data in ",
	       "hdf_file::setc().",exc_einval);
  }
      
  return;
}
    
void hdf_file::setd(std::string name, double d) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::setd().",exc_efailed);
  }

  hid_t dset, space=0;
  bool space_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=1;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);
    dset=H5Dcreate(current,name.c_str(),H5T_IEEE_F64LE,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    if (dset<0) {
      O2SCL_ERR2("Failed to create dataspace in ",
		 "hdf_file::setd().",exc_einval);
    }
    space_alloc=true;

  }

  // Write the data 
  int status=H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,&d);
      
  status=H5Dclose(dset);
  if (space_alloc) status=H5Sclose(space);
      
  return;
}
    
void hdf_file::setf(std::string name, float f) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::setf().",exc_efailed);
  }

  hid_t dset, space=0;
  bool space_alloc=false;

  H5E_BEGIN_TRY
    {
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=1;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);
    dset=H5Dcreate(current,name.c_str(),H5T_IEEE_F32LE,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    if (dset<0) {
      O2SCL_ERR2("Failed to create dataspace in ",
		 "hdf_file::setf().",exc_einval);
    }
    space_alloc=true;

  }

  // Write the data 
  int status=H5Dwrite(dset,H5T_NATIVE_FLOAT,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,&f);
      
  status=H5Dclose(dset);
  if (space_alloc) status=H5Sclose(space);
      
  return;
}
    
void hdf_file::seti(std::string name, int i) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::seti().",exc_efailed);
  }

  hid_t dset, space=0;
  bool space_alloc=false;

  H5E_BEGIN_TRY
    {
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=1;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);
    dset=H5Dcreate(current,name.c_str(),H5T_STD_I32LE,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    if (dset<0) {
      O2SCL_ERR2("Failed to create dataspace in ",
		 "hdf_file::seti().",exc_einval);
    }
    space_alloc=true;

  }

  // Write the data 
  int status=H5Dwrite(dset,H5T_NATIVE_INT,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,&i);
      
  status=H5Dclose(dset);
  if (space_alloc) status=H5Sclose(space);
      
  return;
}

void hdf_file::set_szt(std::string name, size_t u) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::set_szt().",exc_efailed);
  }

  hid_t dset, space=0;
  bool space_alloc=false;

  H5E_BEGIN_TRY
    {
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=1;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);
    dset=H5Dcreate(current,name.c_str(),H5T_STD_U64LE,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    if (dset<0) {
      O2SCL_ERR2("Failed to create dataspace in ",
		 "hdf_file::set_szt().",exc_einval);
    }
    space_alloc=true;

  }

  // Write the data 
  int status;
  if (std::numeric_limits<size_t>::digits==64) {
    status=H5Dwrite(dset,H5T_NATIVE_ULLONG,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,&u);
  } else {
    status=H5Dwrite(dset,H5T_NATIVE_ULONG,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,&u);
  }
  
  status=H5Dclose(dset);
  if (space_alloc) status=H5Sclose(space);
      
  return;
}
    
void hdf_file::sets(std::string name, std::string s) { 
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::sets().",exc_efailed);
  }
  setc_arr(name,s.length(),s.c_str());
  return;
}
    
void hdf_file::sets_fixed(std::string name, std::string s) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::sets_fixed().",exc_efailed);
  }

  hid_t dset, space=0, filetype=0;
  hsize_t str_size=0;
  
  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=1;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);
    
    // Create the file type
    filetype=H5Tcopy(H5T_C_S1);
    herr_t statusx=H5Tset_size(filetype,s.length()+1);
    
    dset=H5Dcreate(current,name.c_str(),filetype,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    if (dset<0) {
      O2SCL_ERR2("Failed to create dataspace in ",
		 "hdf_file::sets_fixed().",exc_einval);
    }
    str_size=s.length()+1;
    
  } else {

    filetype=H5Dget_type(dset);
    str_size=H5Tget_size(filetype);
    
    hsize_t dims[1];
    space=H5Dget_space(dset);
    int ndims=H5Sget_simple_extent_dims(space,dims,0);
    if (ndims!=1 || dims[0]!=1) {
      std::cerr << "Incorrect dimensions in hdf_file::sets_fixed():\n  "
                << ndims << " " << dims[0] << std::endl;
      O2SCL_ERR2("Incorrect dimensions in hdf_file::sets_fixed().",
		 "",exc_einval);
    }
    if (str_size<s.length()+1) {
      O2SCL_ERR2("Not enough space in hdf_file::sets_fixed().",
		 "",exc_einval);
    }

  }

  char *c=new char[str_size];
  for(size_t i=0;i<str_size-1;i++) {
    c[i]=s[i];
  }
  c[str_size-1]='\0';

  // Create the memory type
  hid_t memtype=H5Tcopy(H5T_C_S1);
  herr_t status=H5Tset_size(memtype,str_size);

  // Write the data 
  status=H5Dwrite(dset,memtype,H5S_ALL,
		  H5S_ALL,H5P_DEFAULT,c);
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  status=H5Tclose(memtype);
  status=H5Tclose(filetype);
      
  delete[] c;

  return;
}
    
int hdf_file::gets_fixed(std::string name, std::string &s) {
  
  hid_t dset, space=0, filetype=0;
  
  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
    O2SCL_ERR((((string)"Dataspace named '")+name+
	       "' not found in hdf_file::gets_fixed().").c_str(),
	      exc_einval);
  }

  filetype=H5Dget_type(dset);
  size_t str_size=H5Tget_size(filetype);

  hsize_t dims[1];
  space=H5Dget_space(dset);
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  if (ndims!=1 || dims[0]!=1) {
    O2SCL_ERR2("Incorrect dimensions in hdf_file::gets_fixed().",
	       "",exc_einval);
  }
  
  char *c=new char[str_size];
  
  // Create the memory type
  hid_t memtype=H5Tcopy(H5T_C_S1);
  herr_t status=H5Tset_size(memtype,str_size);
  
  // Write the data 
  herr_t status2=H5Dread(dset,memtype,H5S_ALL,
			 H5S_ALL,H5P_DEFAULT,c);
  
  // Make sure it's terminated properly
  c[str_size-1]='\0';
  s=c;
  
  herr_t status3=H5Dclose(dset);
  status=H5Sclose(space);
  status=H5Tclose(memtype);
  status=H5Tclose(filetype);

  delete[] c;
      
  return 0;
}
    
int hdf_file::gets_def_fixed(std::string name, std::string def,
			     std::string &s) { 

  hid_t dset, space=0, filetype=0;
  
  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, return the default string
  if (dset<0) {
    s=def;
    return 0;
  }

  filetype=H5Dget_type(dset);
  size_t str_size=H5Tget_size(filetype);

  hsize_t dims[1];
  space=H5Dget_space(dset);
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  if (ndims!=1 || dims[0]!=1) {
    O2SCL_ERR2("Incorrect dimensions in gets_def_fixed().",
	       "",exc_einval);
  }
  
  char *c=new char[str_size];
  
  // Create the memory type
  hid_t memtype=H5Tcopy(H5T_C_S1);
  herr_t status=H5Tset_size(memtype,str_size);
  
  // Write the data 
  herr_t status2=H5Dread(dset,memtype,H5S_ALL,
			 H5S_ALL,H5P_DEFAULT,c);
  
  // Make sure it's terminated properly
  c[str_size-1]='\0';
  s=c;
  
  herr_t status3=H5Dclose(dset);
  status=H5Sclose(space);
  status=H5Tclose(memtype);
  status=H5Tclose(filetype);

  delete[] c;
      
  return 0;
}
    
int hdf_file::getc(std::string name, char &c) {
  
  // Open the data space
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  
  if (dset<0) {
    O2SCL_ERR((((string)"Dataspace named '")+name+
	       "' not found in hdf_file::getc().").c_str(),exc_einval);
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,&c);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::getc().",
	      exc_einval);
  }

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::getd(std::string name, double &d) {
      
  // Open the data space
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  if (dset<0) {
    O2SCL_ERR((((string)"Dataspace named '")+name+
	       "' not found in hdf_file::getd().").c_str(),exc_einval);
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,&d);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::getd().",
	      exc_einval);
  }

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::getf(std::string name, float &f) {
      
  // Open the data space
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  if (dset<0) {
    O2SCL_ERR((((string)"Dataspace named '")+name+
	       "' not found in hdf_file::getf().").c_str(),exc_einval);
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,&f);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::getf().",
	      exc_einval);
  }

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::geti(std::string name, int &i) {
      
  // Open the data space
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  if (dset<0) {
    O2SCL_ERR((((string)"Dataspace named '")+name+
	       "' not found in hdf_file::geti().").c_str(),exc_einval);
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,&i);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::geti().",
	      exc_einval);
  }

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::get_szt(std::string name, size_t &u) {
      
  // Open the data space
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  if (dset<0) {
    O2SCL_ERR((((string)"Dataspace named '")+name+
	       "' not found in hdf_file::get_szt().").c_str(),exc_einval);
  }

  // Read the data
  int status;
  if (std::numeric_limits<size_t>::digits==64) {
    status=H5Dread(dset,H5T_NATIVE_HSIZE,H5S_ALL,H5S_ALL,
		   H5P_DEFAULT,&u);
  } else {
    o2_u64_t u2;
    status=H5Dread(dset,H5T_NATIVE_ULLONG,H5S_ALL,H5S_ALL,
		   H5P_DEFAULT,&u2);
    if (u2>((o2_u64_t)std::numeric_limits<size_t>::max())) {
      O2SCL_ERR2("Not enough space in native size_t type in ",
		 "hdf_file::get_szt().",exc_efailed);
    }
    u=(size_t)u2;
  }

  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::get_szt().",
	      exc_einval);
  }

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::gets(std::string name, std::string &s) {

  // Open the data space
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  if (dset<0) {
    O2SCL_ERR((((string)"Dataspace named '")+name+
	       "' not found in hdf_file::gets().").c_str(),exc_einval);
  }

  if (true) {
    // Determine if this is a fixed-length string, and if so, use
    // gets_fixed() instead.
    hid_t filetype=H5Dget_type(dset);
    size_t str_size=H5Tget_size(filetype);
    
    hsize_t dims[3], maxdims[3];
    hid_t space=H5Dget_space(dset);
    int ndims=H5Sget_simple_extent_dims(space,dims,maxdims);
    hid_t memtype=-1;
    if (ndims==1 && maxdims[0]!=H5S_UNLIMITED) {
      memtype=H5Tcopy(H5T_C_S1);
    }
    if (memtype>0) {
      int status=H5Tclose(memtype);
      status=H5Sclose(space);
      status=H5Tclose(filetype);
      status=H5Dclose(dset);
      gets_fixed(name,s);
      return 1;
    } else {
      int status=H5Sclose(space);
      status=H5Tclose(filetype);
    }
  }
  
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  if (ndims!=1) {
    O2SCL_ERR2("Dataspace has incorrect number of dimensions ",
	       "in hdf_file::gets().",exc_einval);
  }

  // Allocate memory
  char *c=new char[dims[0]];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,c);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::gets().",
	      exc_einval);
  }

  // Close the dataset
  status=H5Dclose(dset);

  // Copy to the string object
  s="";
  for(size_t i=0;i<dims[0];i++) s+=c[i];

  // Delete char memory
  delete[] c;

  return 0;
}

int hdf_file::gets_var(std::string name, std::string &s) {
      
  // Open the data space
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  if (dset<0) {
    O2SCL_ERR((((string)"Dataspace named '")+name+
	       "' not found in hdf_file::gets().").c_str(),exc_einval);
  }

  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  if (ndims!=1) {
    O2SCL_ERR2("Dataspace has incorrect number of dimensions ",
	       "in hdf_file::gets().",exc_einval);
  }

  // Allocate memory
  char *c=new char[dims[0]];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,c);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::gets().",
	      exc_einval);
  }

  // Close the dataset
  status=H5Dclose(dset);

  // Copy to the string object
  s="";
  for(size_t i=0;i<dims[0];i++) s+=c[i];

  // Delete char memory
  delete[] c;

  return 0;
}

int hdf_file::getc_def(std::string name, char def, char &c) {
   
  hid_t dset=0;
   
  H5E_BEGIN_TRY {
    // See if the dataspace already exists first
    dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  }
  H5E_END_TRY
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
  
  // Not found, return default
  if (dset<0) {
    c=def;
    return success;
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,&c);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::getc_def().",
	      exc_einval);
  }

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::getd_def(std::string name, double def, double &d) {
      
  hid_t dset=0;
   
  H5E_BEGIN_TRY {
    // See if the dataspace already exists first
    dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  }
  H5E_END_TRY
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
  
  // Not found, return default
  if (dset<0) {
    d=def;
    return success;
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,&d);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::getd_def().",
	      exc_einval);
  }

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::getf_def(std::string name, float def, float &f) {
      
  hid_t dset=0;
   
  H5E_BEGIN_TRY {
    // See if the dataspace already exists first
    dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  }
  H5E_END_TRY
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
  
  // Not found, return default
  if (dset<0) {
    f=def;
    return success;
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,&f);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::getf_def().",
	      exc_einval);
  }

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::get_szt_def(std::string name, size_t def, size_t &u) {
      
  hid_t dset=0;
   
  H5E_BEGIN_TRY {
    // See if the dataspace already exists first
    dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  }
  H5E_END_TRY
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
  
  // Not found, return default
  if (dset<0) {
    u=def;
    return success;
  }
  
  // Read the data
  int status;
  if (std::numeric_limits<size_t>::digits==64) {
    status=H5Dread(dset,H5T_NATIVE_HSIZE,H5S_ALL,H5S_ALL,
		   H5P_DEFAULT,&u);
  } else {
    o2_u64_t u2;
    status=H5Dread(dset,H5T_NATIVE_ULLONG,H5S_ALL,H5S_ALL,
		   H5P_DEFAULT,&u2);
    if (u2>((o2_u64_t)std::numeric_limits<size_t>::max())) {
      O2SCL_ERR2("Not enough space in native size_t type in ",
		 "hdf_file::get_szt_def().",exc_efailed);
    }
    u=(size_t)u2;
  }
  
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::get_szt_def().",
	      exc_einval);
  }

  status=H5Dclose(dset);
  
  return 0;
}

int hdf_file::geti_def(std::string name, int def, int &i) {
      
  hid_t dset=0;
   
  H5E_BEGIN_TRY {
    // See if the dataspace already exists first
    dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  }
  H5E_END_TRY
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
  
  // Not found, return default
  if (dset<0) {
    i=def;
    return success;
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,&i);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::geti_def().",
	      exc_einval);
  }

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::gets_def(std::string name, std::string def, std::string &s) {
      
  // Open the data space
  hid_t dset=0;

  H5E_BEGIN_TRY {
    // See if the dataspace already exists first
    dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  }
  H5E_END_TRY
#ifdef O2SCL_NEVER_DEFINED
    }{
#endif

  // Not found, return default
  if (dset<0) {
    s=def;
    return success;
  }

  if (false) {
    /*
      There's a problem here because this appears to 
      improperly detect a 1 character variable length string as a 
      fixed-length string.
    */

    // Determine if this is a fixed-length string, and if so, use
    // gets_fixed() instead.
    hid_t filetype=H5Dget_type(dset);
    size_t str_size=H5Tget_size(filetype);
    
    hsize_t dims[3];
    hid_t space=H5Dget_space(dset);
    int ndims=H5Sget_simple_extent_dims(space,dims,0);
    hid_t memtype=-1;
    if (ndims==1 && dims[0]==1) {
      memtype=H5Tcopy(H5T_C_S1);
    }
    if (memtype>0) {
      int status=H5Tclose(memtype);
      status=H5Sclose(space);
      status=H5Tclose(filetype);
      status=H5Dclose(dset);
      gets_fixed(name,s);
      return 1;
    } else {
      int status=H5Sclose(space);
      status=H5Tclose(filetype);
    }
  }

  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  if (ndims!=1) {
    O2SCL_ERR2("Dataspace has incorrect number of dimensions ",
	       "in hdf_file::gets_def().",exc_einval);
  }

  // Allocate memory
  char *c=new char[dims[0]];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,c);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::gets_def().",
	      exc_einval);
  }

  // Close the dataset
  status=H5Dclose(dset);

  // Copy to the string object
  s="";
  for(size_t i=0;i<dims[0];i++) s+=c[i];

  // Delete char memory
  delete[] c;

  return 0;
}

int hdf_file::gets_var_def(std::string name, std::string def, std::string &s) {
      
  // Open the data space
  hid_t dset=0;

  H5E_BEGIN_TRY {
    // See if the dataspace already exists first
    dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  }
  H5E_END_TRY
#ifdef O2SCL_NEVER_DEFINED
    }{
#endif

  // Not found, return default
  if (dset<0) {
    s=def;
    return success;
  }

  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  if (ndims!=1) {
    O2SCL_ERR2("Dataspace has incorrect number of dimensions ",
	       "in hdf_file::gets_def().",exc_einval);
  }

  // Allocate memory
  char *c=new char[dims[0]];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,c);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::gets_def().",
	      exc_einval);
  }

  // Close the dataset
  status=H5Dclose(dset);

  // Copy to the string object
  s="";
  for(size_t i=0;i<dims[0];i++) s+=c[i];

  // Delete char memory
  delete[] c;

  return 0;
}

int hdf_file::setc_arr_fixed(std::string name, size_t n, const char *c) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::setc_arr_fixed().",exc_efailed);
  }

  if (n==0) {
    O2SCL_ERR("Tried to call setc_arr_fixed() with n=0.",exc_einval);
  }

  hid_t dset, space;
  bool space_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=n;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);

    dset=H5Dcreate(current,name.c_str(),H5T_STD_I8LE,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    space_alloc=true;

  } else {
    
    // Get space requirements, to make sure they coincide
    // with the size specified by the user
    space=H5Dget_space(dset);  
    hsize_t dims[1];
    int ndims=H5Sget_simple_extent_dims(space,dims,0);
    if (ndims!=1 || dims[0]!=n) {
      O2SCL_ERR2("Incompatible dataspace size or dimensions ",
		 "in hdf_file::setc_arr_fixed().",exc_einval);
    }
    
  }

  // Write the data 
  int status=H5Dwrite(dset,H5T_NATIVE_CHAR,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,c);
      
  status=H5Dclose(dset);
  if (space_alloc) status=H5Sclose(space);
      
  return 0;
}

int hdf_file::setd_arr_fixed(std::string name, size_t n, const double *d) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::setd_arr_fixed().",exc_efailed);
  }

  if (n==0) {
    O2SCL_ERR("Tried to call setd_arr_fixed() with n=0.",exc_einval);
  }
  
  hid_t dset, space;
  bool space_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=n;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);

    dset=H5Dcreate(current,name.c_str(),H5T_IEEE_F64LE,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    space_alloc=true;

  } else {
    
    // Get space requirements, to make sure they coincide
    // with the size specified by the user
    space=H5Dget_space(dset);  
    hsize_t dims[1];
    int ndims=H5Sget_simple_extent_dims(space,dims,0);
    if (ndims!=1 || dims[0]!=n) {
      O2SCL_ERR2("Incompatible dataspace size or dimensions ",
		 "in hdf_file::setd_arr_fixed().",exc_einval);
    }

  }

  // Write the data 
  int status=H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,d);
      
  status=H5Dclose(dset);
  if (space_alloc) status=H5Sclose(space);
      
  return 0;
}

int hdf_file::setf_arr_fixed(std::string name, size_t n, const float *f) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::setf_arr_fixed().",exc_efailed);
  }

  if (n==0) {
    O2SCL_ERR("Tried to call setf_arr_fixed() with n=0.",exc_einval);
  }

  hid_t dset, space;
  bool space_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=n;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);

    dset=H5Dcreate(current,name.c_str(),H5T_IEEE_F32LE,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    space_alloc=true;

  } else {
    
    // Get space requirements, to make sure they coincide
    // with the size specified by the user
    space=H5Dget_space(dset);  
    hsize_t dims[1];
    int ndims=H5Sget_simple_extent_dims(space,dims,0);
    if (ndims!=1 || dims[0]!=n) {
      O2SCL_ERR2("Incompatible dataspace size or dimensions ",
		 "in hdf_file::setf_arr_fixed().",exc_einval);
    }
    
  }

  // Write the data 
  int status=H5Dwrite(dset,H5T_NATIVE_FLOAT,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,f);
      
  status=H5Dclose(dset);
  if (space_alloc) status=H5Sclose(space);
      
  return 0;
}

int hdf_file::seti_arr_fixed(std::string name, size_t n, const int *i) { 

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::seti_arr_fixed().",exc_efailed);
  }

  if (n==0) {
    O2SCL_ERR("Tried to call seti_arr_fixed() with n=0.",exc_einval);
  }

  hid_t dset, space;
  bool space_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
	
    // Create the dataspace
    hsize_t dims=n;
    // Arguments to H5Screate_simple are (rank,current_dims,max_dims)
    // and if max_dims is 0 then max_dims=current_dims
    space=H5Screate_simple(1,&dims,0);

    dset=H5Dcreate(current,name.c_str(),H5T_STD_I32LE,space,H5P_DEFAULT,
		   H5P_DEFAULT,H5P_DEFAULT);
    space_alloc=true;

  } else {
    
    // Get space requirements, to make sure they coincide
    // with the size specified by the user
    space=H5Dget_space(dset);  
    hsize_t dims[1];
    int ndims=H5Sget_simple_extent_dims(space,dims,0);
    if (ndims!=1 || dims[0]!=n) {
      O2SCL_ERR2("Incompatible dataspace size or dimensions ",
		 "in hdf_file::seti_arr_fixed().",exc_einval);
    }
    
  }

  // Write the data 
  int status=H5Dwrite(dset,H5T_NATIVE_INT,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,i);
      
  status=H5Dclose(dset);
  if (space_alloc) status=H5Sclose(space);
      
  return 0;
}
    
int hdf_file::setc_arr(std::string name, size_t n, const char *c) { 
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::setc_arr().",exc_efailed);
  }

  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
    
    // Create the dataspace
    hsize_t dims=n;
    hsize_t max=H5S_UNLIMITED;
    space=H5Screate_simple(1,&dims,&max);

    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk=def_chunk(n);
    int status2=H5Pset_chunk(dcpl,1,&chunk);

#ifdef O2SCL_HDF5_COMP
    if (n>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::setc_arr().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_STD_I8LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims;
    int ndims=H5Sget_simple_extent_dims(space,&dims,0);

    // Set error if this dataset is more than 1-dimensional
    if (ndims!=1) {
      O2SCL_ERR2("Tried to set a multidimensional dataset with an ",
		 "array in hdf_file::setc_arr().",exc_einval);
    }

    // If necessary, extend the dataset
    if (n!=dims) {
      hsize_t new_dims=n;
      int status3=H5Dset_extent(dset,&new_dims);
    }
    
  }

  // Write the data 
  int status;
  if (n==0) {
    // *FIXME* Does this character pointer really need to be non-null?
    // Either remove this or comment this more clearly.
    char c2[1]={' '};
    status=H5Dwrite(dset,H5T_NATIVE_CHAR,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,c2);
  } else {
    status=H5Dwrite(dset,H5T_NATIVE_CHAR,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,c);
  }
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }
      
  return 0;
}

int hdf_file::setd_arr(std::string name, size_t n, const double *d) {
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::setd_arr().",exc_efailed);
  }

  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
    
#ifdef O2SCL_NEVER_DEFINED
    int mpi_size, mpi_rank;
    MPI_Comm comm=MPI_COMM_WORLD;
    MPI_Info info=MPI_INFO_NULL;
    MPI_Comm_size(comm,&mpi_size);
    MPI_Comm_rank(comm,&mpi_rank);
#endif
    
    // Create the dataspace
    hsize_t dims=n;
    hsize_t max=H5S_UNLIMITED;
    space=H5Screate_simple(1,&dims,&max);

    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk=def_chunk(n);
    int status2=H5Pset_chunk(dcpl,1,&chunk);

#ifdef O2SCL_HDF5_COMP    
    if (n>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::setd_arr().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_IEEE_F64LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims;
    int ndims=H5Sget_simple_extent_dims(space,&dims,0);

    // Set error if this dataset is more than 1-dimensional
    if (ndims!=1) {
      O2SCL_ERR2("Tried to set a multidimensional dataset with an ",
		 "array in hdf_file::setd_arr().",exc_einval);
    }

    // If necessary, extend the dataset
    if (n!=dims) {
      hsize_t new_dims=n;
      int status3=H5Dset_extent(dset,&new_dims);
    }
    
  }

  // Write the data 
  int status;
  if (n==0) {
    double d2[1]={0.0};
    status=H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,d2);
  } else {
    status=H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,d);
  }
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }
      
  return 0;
}

int hdf_file::setf_arr(std::string name, size_t n, const float *f) { 
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::setf_arr().",exc_efailed);
  }

  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
    
    // Create the dataspace
    hsize_t dims=n;
    hsize_t max=H5S_UNLIMITED;
    space=H5Screate_simple(1,&dims,&max);

    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk=def_chunk(n);
    int status2=H5Pset_chunk(dcpl,1,&chunk);

#ifdef O2SCL_HDF5_COMP    
    if (n>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::setf_arr().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_IEEE_F32LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims;
    int ndims=H5Sget_simple_extent_dims(space,&dims,0);

    // Set error if this dataset is more than 1-dimensional
    if (ndims!=1) {
      O2SCL_ERR2("Tried to set a multidimensional dataset with an ",
		 "array in hdf_file::setf_arr().",exc_einval);
    }

    // If necessary, extend the dataset
    if (n!=dims) {
      hsize_t new_dims=n;
      int status3=H5Dset_extent(dset,&new_dims);
    }
    
  }

  // Write the data 
  int status;
  if (n==0) {
    float f2[1]={0.0};
    status=H5Dwrite(dset,H5T_NATIVE_FLOAT,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,f2);
  } else {
    status=H5Dwrite(dset,H5T_NATIVE_FLOAT,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,f);
  }
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }
      
  return 0;
}

int hdf_file::seti_arr(std::string name, size_t n, const int *i) { 
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::seti_arr().",exc_efailed);
  }

  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
    
    // Create the dataspace
    hsize_t dims=n;
    hsize_t max=H5S_UNLIMITED;
    space=H5Screate_simple(1,&dims,&max);

    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk=def_chunk(n);
    int status2=H5Pset_chunk(dcpl,1,&chunk);

#ifdef O2SCL_HDF5_COMP    
    if (n>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::seti_arr().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_STD_I32LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims;
    int ndims=H5Sget_simple_extent_dims(space,&dims,0);

    // Set error if this dataset is more than 1-dimensional
    if (ndims!=1) {
      O2SCL_ERR2("Tried to set a multidimensional dataset with an ",
		 "array in hdf_file::seti_arr().",exc_einval);
    }

    // If necessary, extend the dataset
    if (n!=dims) {
      hsize_t new_dims=n;
      int status3=H5Dset_extent(dset,&new_dims);
    }
    
  }

  // Write the data 
  int status;
  if (n==0) {
    int i2[1]={0};
    status=H5Dwrite(dset,H5T_NATIVE_INT,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,i2);
  } else {
    status=H5Dwrite(dset,H5T_NATIVE_INT,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,i);
  }
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }
      
  return 0;
}

int hdf_file::set_szt_arr(std::string name, size_t n, const size_t *u) { 
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access in ",
	       "hdf_file::set_szt_arr().",exc_efailed);
  }

  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
    
    // Create the dataspace
    hsize_t dims=n;
    hsize_t max=H5S_UNLIMITED;
    space=H5Screate_simple(1,&dims,&max);

    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk=def_chunk(n);
    int status2=H5Pset_chunk(dcpl,1,&chunk);

#ifdef O2SCL_HDF5_COMP    
    if (n>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::set_szt_arr().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_STD_U64LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims;
    int ndims=H5Sget_simple_extent_dims(space,&dims,0);

    // Set error if this dataset is more than 1-dimensional
    if (ndims!=1) {
      O2SCL_ERR2("Tried to set a multidimensional dataset with an ",
		 "array in hdf_file::set_szt_arr().",exc_einval);
    }

    // If necessary, extend the dataset
    if (n!=dims) {
      hsize_t new_dims=n;
      int status3=H5Dset_extent(dset,&new_dims);
    }
    
  }

  // Write the data 
  int status;
  if (std::numeric_limits<size_t>::digits==64) {
    if (n==0) {
      size_t u2[1]={0};
      status=H5Dwrite(dset,H5T_NATIVE_ULLONG,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,u2);
    } else {
      status=H5Dwrite(dset,H5T_NATIVE_ULLONG,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,u);
    }
  } else {
    if (n==0) {
      size_t u2[1]={0};
      status=H5Dwrite(dset,H5T_NATIVE_ULONG,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,u2);
    } else {
      status=H5Dwrite(dset,H5T_NATIVE_ULONG,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,u);
    }
  }
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }
      
  return 0;
}

int hdf_file::getc_arr(std::string name, size_t n, char *c) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);

  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  if (dims[0]!=n) {
    O2SCL_ERR("Incompatible size in hdf_file::getc_arr().",
	      exc_einval);
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,c);

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::getd_arr(std::string name, size_t n, double *d) {
  int compr;
  return getd_arr_compr(name,n,d,compr);
}

int hdf_file::getd_arr_compr(std::string name, size_t n, double *d,
			     int &compr) {

  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);

  // Get filter information
  compr=0;
  hid_t plist_id=H5Dget_create_plist(dset);
  int num_filters=H5Pget_nfilters(plist_id);
  for(int i=0;i<num_filters;i++) {
    size_t n_elements=0;
    unsigned flags, filter_info;
    H5Z_filter_t filter_type=H5Pget_filter2
      (plist_id,0,&flags,&n_elements,0,0,0,&filter_info);
    if (filter_type==H5Z_FILTER_DEFLATE) {
      compr=1;
    } else if (filter_type==H5Z_FILTER_SZIP) {
      compr=2;
    }
  }

  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  if (dims[0]!=n) {
    string str="Asked for size "+itos(n)+" but file has size "+
      itos(dims[0])+" in hdf_file::getd_arr_compr().";
    O2SCL_ERR(str.c_str(),exc_einval);
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,d);

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::getf_arr(std::string name, size_t n, float *f) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);

  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  if (dims[0]!=n) {
    O2SCL_ERR("Incompatible size in hdf_file::getf_arr().",
	      exc_einval);
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,f);

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::geti_arr(std::string name, size_t n, int *i) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);

  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  if (dims[0]!=n) {
    O2SCL_ERR("Incompatible size in hdf_file::geti_arr().",
	      exc_einval);
  }

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,i);

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::getc_arr_alloc(std::string name, size_t &n, char *c) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);

  // Get space requirements, to allocate memory
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  n=dims[0];
  c=new char[n];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,c);

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::getd_arr_alloc(std::string name, size_t &n, double *d) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);

  // Get space requirements, to allocate memory
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  n=dims[0];
  d=new double[n];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,d);

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::getf_arr_alloc(std::string name, size_t &n, float *f) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);

  // Get space requirements, to allocate memory
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  n=dims[0];
  f=new float[n];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,f);

  status=H5Dclose(dset);

  return 0;
}

int hdf_file::geti_arr_alloc(std::string name, size_t &n, int *i) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);

  // Get space requirements, to allocate memory
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  n=dims[0];
  i=new int[n];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,i);

  status=H5Dclose(dset);

  return 0;
}

hid_t hdf_file::open_group(std::string path) {
  if (!file_open) {
    O2SCL_ERR("File not opened in hdf_file::open_group().",
	      exc_einval);
  }
  if (current<=0) {
    O2SCL_ERR("Invalid current HDF5 id in hdf_file::open_group().",
	      exc_einval);
  }
  hid_t group;
  H5E_BEGIN_TRY
    {
      group=H5Gopen(current,path.c_str(),H5P_DEFAULT);
    }
  H5E_END_TRY 
    if (group<0) {
      if (write_access==false) {
	std::string str=((std::string)"File not opened with write access ")+
	  "and group with path '"+path+"' not found in hdf_file::"+
	  "open_group(std::string).";
	O2SCL_ERR(str.c_str(),exc_efailed);
      }
      group=H5Gcreate(current,path.c_str(),H5P_DEFAULT,
		      H5P_DEFAULT,H5P_DEFAULT);
    }
  if (group<0) {
    O2SCL_ERR2("Failed to open or create group in ",
	       "hdf_file::open_group().",exc_einval);
  }
  return group;
}

hid_t hdf_file::open_group(hid_t init_id, std::string path) {
  hid_t group;
  H5E_BEGIN_TRY
    {
      group=H5Gopen(init_id,path.c_str(),H5P_DEFAULT);
    }
  H5E_END_TRY 
    if (group<0) {
      if (write_access==false) {
	O2SCL_ERR2("File not opened with write access and group not found",
		   "in hdf_file::open_group(hid_t,std::string).",exc_efailed);
      }
      group=H5Gcreate(current,path.c_str(),H5P_DEFAULT,
		      H5P_DEFAULT,H5P_DEFAULT);
    }
  if (group<0) {
    O2SCL_ERR2("Failed to open or create group in ",
	       "hdf_file::open_group().",exc_einval);
  }
  return group;
}

int hdf_file::setd_vec(std::string name, const std::vector<double> &v) {
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::setd_vec().",exc_efailed);
  }
  if (v.size()==0) {
    return setd_arr(name,0,0);
  }
  return setd_arr(name,v.size(),&v[0]);
}

int hdf_file::seti_vec(std::string name, const std::vector<int> &v) {
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::seti_vec().",exc_efailed);
  }
  if (v.size()==0) {
    return seti_arr(name,0,0);
  }
  return seti_arr(name,v.size(),&v[0]);
}

int hdf_file::set_szt_vec(std::string name, const std::vector<size_t> &v) {
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::set_szt_vec().",exc_efailed);
  }
  if (v.size()==0) {
    return set_szt_arr(name,0,0);
  }
  return set_szt_arr(name,v.size(),&v[0]);
}

int hdf_file::getd_vec(std::string name, std::vector<double> &v) {
      
  // Get the dataset with the specified name
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  
  // Get the dimension of the dataset
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  herr_t status;

  if (dims[0]>0) {
    
    // Resize the vector
    v.resize(dims[0]);
    
    // Read the data directly into the pointer
    status=H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
		   H5P_DEFAULT,&v[0]);
    if (status<0) {
      O2SCL_ERR("Could not read dataspace in hdf_file::getd_vec().",
		exc_einval);
    }
    
  }
  
  status=H5Dclose(dset);
  
  return 0;
}

int hdf_file::geti_vec(std::string name, std::vector<int> &v) {
      
  // Get the dataset with the specified name
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  
  // Get the dimension of the dataset
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  herr_t status;

  if (dims[0]>0) {
    
    // Resize the vector
    v.resize(dims[0]);
    
    // Read the data directly into the pointer
    status=H5Dread(dset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
		   H5P_DEFAULT,&v[0]);
    if (status<0) {
      O2SCL_ERR("Could not read dataspace in hdf_file::geti_vec().",
		exc_einval);
    }
    
  }
  
  status=H5Dclose(dset);
  
  return 0;
}

int hdf_file::get_szt_vec(std::string name, std::vector<size_t> &v) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  
  herr_t status;

  if (dims[0]>0) {
    
    v.resize(dims[0]);
 
    // Read the data
    if (std::numeric_limits<size_t>::digits==64) {
      status=H5Dread(dset,H5T_NATIVE_HSIZE,H5S_ALL,H5S_ALL,
		     H5P_DEFAULT,&v[0]);
      if (status<0) {
	O2SCL_ERR("Could not read dataspace in hdf_file::get_szt_vec().",
		  exc_einval);
      }
    } else {
      o2_u64_t *tmp_arr=new o2_u64_t[dims[0]];
      status=H5Dread(dset,H5T_NATIVE_HSIZE,H5S_ALL,H5S_ALL,
		     H5P_DEFAULT,&tmp_arr[0]);
      if (status<0) {
	O2SCL_ERR("Could not read dataspace in hdf_file::get_szt_vec().",
		  exc_einval);
      }
      for(size_t i=0;i<dims[0];i++) {
	if (tmp_arr[i]>((o2_u64_t)std::numeric_limits<size_t>::max())) {
	  O2SCL_ERR2("Not enough space in native size_t type in ",
		     "hdf_file::get_szt_vec().",exc_efailed);
	}
	v[i]=(size_t)tmp_arr[i];
      }
      delete[] tmp_arr;
    }
       
  }
  
  status=H5Dclose(dset);
  
  return 0;
}

int hdf_file::gets_vec_copy(std::string name, std::vector<std::string> &s) {
		  
  int nc, nw;

  // Open the group
  hid_t top=get_current_id();
  hid_t group=open_group(name);
  set_current_id(group);
  
  string o2t;
  gets_fixed("o2scl_type",o2t);
  if (o2t!="string[]") {
    set_current_id(top);
    O2SCL_ERR2("The specified name does not refer to data which ",
	       "can be read by O2scl in hdf_file::gets_vec().",
	       exc_efailed);
  }

  // Get number of words
  geti("nw",nw);

  if (nw>0) {
    
    // Get number of characters
    geti("nc",nc);
    
    if (nc>0) {

      // Allocate space for ip and cp
      vector<int> ip(nw);
      std::string cp;
      
      // Get counter and data
      geti_vec("counter",ip);
      gets("data",cp);
      
      // Copy the data over
      size_t ix=0;
      for(int i=0;i<nw;i++) {
	string tmp;
	for(int j=0;j<ip[i];j++) {
	  tmp+=cp[ix];
	  ix++;
	}
	s.push_back(tmp);
      }

    } else {

      // If nw>0 and nc==0 then there are a number of strings,
      // all of which are empty.
      for(int i=0;i<nw;i++) {
        s.push_back("");
      }
      
    }
    
  }

  close_group(group);

  // Return file location
  set_current_id(top);

  return 0;
}

int hdf_file::sets_vec_copy(std::string name,
                            const std::vector<std::string> &s) {

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::sets_vec().",exc_efailed);
  }

  // Create the group
  hid_t top=current;
  hid_t group=open_group(name);
  set_current_id(group);
  
  sets_fixed("o2scl_type","string[]");

  seti("nw",s.size());
  
  // Compute total length
  int nc=0;
  for(size_t i=0;i<s.size();i++) nc+=s[i].length();
  
  seti("nc",nc);
  
  // count lengths

  if (nc>0) {
    vector<int> ip;
    std::string cp;
    for(size_t i=0;i<s.size();i++) {
      ip.push_back(s[i].length());
      for(size_t j=0;j<s[i].size();j++) {
        cp+=s[i][j];
      }
    }
  
    // Set data
    seti_vec("counter",ip);
    sets("data",cp);
  }

  // Close the group
  close_group(group);
  
  // Return file location
  set_current_id(top);

  return 0;
}

int hdf_file::gets_ten_copy(std::string name, 
                            tensor_base<std::string,
                            std::vector<std::string>,
                            std::vector<size_t>> &t) {
		  
  // Open the group
  hid_t top=get_current_id();
  hid_t group=open_group(name);
  set_current_id(group);

  // Obtain the type and verify that it's "tensor_string"
  string o2t;
  gets_fixed("o2scl_type",o2t);
  if (o2t!="tensor_string") {
    set_current_id(top);
    O2SCL_ERR2("The specified name does not refer to data which ",
	       "can be read by O2scl in hdf_file::gets_vec().",
	       exc_efailed);
  }

  // Get rank
  int rank;
  geti("rank",rank);

  // Get sizes
  vector<size_t> size;
  get_szt_arr("size",size);

  // Resize
  t.resize(rank,size);

  // Read length tensor
  tensor<size_t> lens;
  get_szt_arr("lengths",lens);

  // Read character array as a string object
  std::string s;
  gets("data",s);

  // Set references for data objects
  const std::vector<std::string> &sdat=t.get_data();
  const std::vector<size_t> &ldat=tlens.get_data();

  // Set strings from character array
  size_t ctr=0;
  for(size_t i=0;i<t.total_size();i++) {
    sdat[i]="";
    for(size_t j=0;j<ldat[i];j++) {
      sdat[i]+=s[ctr];
      ctr++;
    }
  }

  // Close the associated group
  close_group(group);

  // Return file location
  set_current_id(top);

  return 0;
}

int hdf_file::sets_ten_copy(std::string name,
                            const tensor_base<std::string,
                            std::vector<std::string>,
                            std::vector<size_t>> &t) {

  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::sets_ten().",exc_efailed);
  }

  // Create the group
  hid_t top=current;
  hid_t group=open_group(name);
  set_current_id(group);
  
  sets_fixed("o2scl_type","tensor_string");

  seti("rank",t.get_rank());

  if (t.get_rank>0) {
    
    const std::vector<size_t> &s=t.get_size_arr();
    set_szt_arr("size",s);

    // Create a tensor to store string lengths
    tensor<size_t> tlens;
    tlens.resize(t.get_rank(),s);

    // Obtain data references
    const std::vector<std::string> &sdat=t.get_data();
    const std::vector<size_t> &ldat=tlens.get_data();

    // Count total number of characters and set string length tensor
    size_t nchars=0;
    for(size_t i=0;i<sdat.size();i++) {
      ldat[i]=sdat[i].length();
      chars+=ldat[i];
    }

    // Store string length tensor
    set_szt_ten("lengths",tlens);

    // Allocate space for all characters
    char *cp=new char[nchars];

    // Fill array with characters
    size_t ctr=0;
    for(size_t i=0;i<sdat.size();i++) {
      for(size_t j=0;j<sdat[i].length();j++) {
        cp[ctr]=sdat[i][j];
        ctr++;
      }
    }

    // Store character array
    setc_arr("data",ctr,cp);

    // Free memory associated with character array
    delete[] cp;
    
  }
  

  // Close the group
  close_group(group);
  
  // Return file location
  set_current_id(top);

  return 0;
}

int hdf_file::gets_vec_vec_copy
(std::string name, std::vector<std::vector<std::string>> &s) {
  
  // Open the group
  hid_t top=get_current_id();
  hid_t group=open_group(name);
  set_current_id(group);

  string o2t;
  gets_fixed("o2scl_type",o2t);
  if (o2t!="vec_vec_string") {
    set_current_id(top);
    O2SCL_ERR2("The specified name does not refer to data which ",
	       "can be read by O₂scl in hdf_file::gets_vec().",
	       exc_efailed);
  }

  size_t n;
  get_szt("n",n);
  
  if (n>0) {
    
    s.resize(n);
    
    vector<int> sizes;
    geti_vec("sizes",sizes);
    vector<int> lengths;
    geti_vec("lengths",lengths);

    std::string cp;
    gets("data",cp);

    size_t ic=0;
    size_t ik=0;
    for(size_t i=0;i<n;i++) {
      s[i].resize(sizes[i]);
      for(size_t j=0;j<s[i].size();j++) {
        int length=lengths[ik];
        
        for(int k=0;k<length;k++) {
          s[i][j]+=cp[ic];
          ic++;
        }
        
        ik++;
      }
    }
  }

  close_group(group);

  // Return file location
  set_current_id(top);

  return 0;
}

int hdf_file::getd_vec_vec_copy
(std::string name, std::vector<std::vector<double>> &vvd) {
  
  // Open the group
  hid_t top=get_current_id();
  hid_t group=open_group(name);
  set_current_id(group);

  string o2t;
  gets_fixed("o2scl_type",o2t);
  if (o2t!="vec_vec_double") {
    set_current_id(top);
    O2SCL_ERR2("The specified name does not refer to data which ",
	       "can be read by O₂scl in hdf_file::gets_vec().",
	       exc_efailed);
  }

  size_t n;
  get_szt("n",n);
  
  if (n>0) {
    
    vvd.resize(n);
    
    vector<int> sizes;
    geti_vec("sizes",sizes);

    vector<double> vd;
    getd_vec("data",vd);

    size_t ic=0;
    for(size_t i=0;i<n;i++) {
      vvd[i].resize(sizes[i]);
      for(size_t j=0;j<vvd[i].size();j++) {
        vvd[i][j]=vd[ic];
        ic++;
      }
    }
  }

  close_group(group);

  // Return file location
  set_current_id(top);

  return 0;
}

int hdf_file::sets_vec_vec_copy
(std::string name, const std::vector<std::vector<std::string>> &s) {
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::sets_vec_vec().",exc_efailed);
  }

  // Create the group
  hid_t top=current;
  hid_t group=open_group(name);
  set_current_id(group);
  
  sets_fixed("o2scl_type","vec_vec_string");
  
  set_szt("n",s.size());

  if (s.size()>0) {
    
    vector<int> sizes;
    vector<int> lengths;
    for(size_t i=0;i<s.size();i++) {
      sizes.push_back(s[i].size());
      for(size_t j=0;j<s[i].size();j++) {
        lengths.push_back(s[i][j].length());
      }
    }
    seti_vec("sizes",sizes);
    seti_vec("lengths",lengths);
    
    std::string cp;
    for(size_t i=0;i<s.size();i++) {
      for(size_t j=0;j<s[i].size();j++) {
        for(size_t k=0;k<s[i][j].length();k++) {
          cp+=s[i][j][k];
        }
      }
    }
    
    sets("data",cp);
  }
  
  // Close the group
  close_group(group);
  
  // Return file location
  set_current_id(top);

  return 0;
}

int hdf_file::setd_vec_vec_copy
(std::string name, const std::vector<std::vector<double>> &vvd) {
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::sets_vec_vec().",exc_efailed);
  }

  // Create the group
  hid_t top=current;
  hid_t group=open_group(name);
  set_current_id(group);
  
  sets_fixed("o2scl_type","vec_vec_double");
  
  set_szt("n",vvd.size());

  if (vvd.size()>0) {
    
    vector<int> sizes;
    
    for(size_t i=0;i<vvd.size();i++) {
      sizes.push_back(vvd[i].size());
    }
    seti_vec("sizes",sizes);

    vector<double> vd;
    for(size_t i=0;i<vvd.size();i++) {
      for(size_t j=0;j<vvd[i].size();j++) {
        vd.push_back(vvd[i][j]);
      }
    }
    
    setd_vec("data",vd);
  }
  
  // Close the group
  close_group(group);
  
  // Return file location
  set_current_id(top);

  return 0;
}

int hdf_file::setd_mat_copy(std::string name, const ubmatrix &m) {
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::setd_mat_copy().",exc_efailed);
  }

  // Copy to a C-style array
  double *d=new double[m.size1()*m.size2()];
  for(size_t i=0;i<m.size1();i++) {
    for(size_t j=0;j<m.size2();j++) {
      d[i*m.size2()+j]=m(i,j);
    }
  }
  
  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
    
    // Create the dataspace
    hsize_t dims[2]={m.size1(),m.size2()};
    hsize_t max[2]={H5S_UNLIMITED,H5S_UNLIMITED};
    space=H5Screate_simple(2,dims,max);

    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[2]={def_chunk(m.size1()),def_chunk(m.size2())};
    int status2=H5Pset_chunk(dcpl,2,chunk);

#ifdef O2SCL_HDF5_COMP    
    if (m.size1()*m.size2()>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::setd_mat_copy().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_IEEE_F64LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims[2];
    int ndims=H5Sget_simple_extent_dims(space,dims,0);

    // Set error if this dataset is more than 1-dimensional
    if (ndims!=2) {
      O2SCL_ERR2("Tried to set a non-matrix dataset with a ",
		 "matrix in hdf_file::setd_mat_copy().",exc_einval);
    }

    // If necessary, extend the dataset
    if (m.size1()!=dims[0] || m.size2()!=dims[1]) {
      hsize_t new_dims[2]={m.size1(),m.size2()};
      int status3=H5Dset_extent(dset,new_dims);
    }
    
  }

  // Write the data 
  int status;
  status=H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,
		  H5S_ALL,H5P_DEFAULT,d);
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }

  // Free the C-style matrix
  delete[] d;
      
  return 0;
}

int hdf_file::getd_mat_copy(std::string name, ubmatrix &m) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[2];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);

  if (ndims!=2) {
    O2SCL_ERR2("Dimensions of dataspace do not match ",
	       "a matrix in hdf_file::getd_mat().",exc_efailed);
  }
  
  m.resize(dims[0],dims[1]);
  double *d=new double[m.size1()*m.size2()];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,d);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::getd_mat_copy().",
	      exc_einval);
  }

  for(size_t i=0;i<m.size1();i++) {
    for(size_t j=0;j<m.size2();j++) {
      m(i,j)=d[i*m.size2()+j];
    }
  }
  delete[] d;
  
  status=H5Dclose(dset);

  return 0;
}

int hdf_file::seti_mat_copy(std::string name, const ubmatrix_int &m) {
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::seti_mat_copy().",exc_efailed);
  }

  int *d=new int[m.size1()*m.size2()];
  for(size_t i=0;i<m.size1();i++) {
    for(size_t j=0;j<m.size2();j++) {
      d[i*m.size2()+j]=m(i,j);
    }
  }
  
  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
      
  // If it doesn't exist, create it
  if (dset<0) {
    
    // Create the dataspace
    hsize_t dims[2]={m.size1(),m.size2()};
    hsize_t max[2]={H5S_UNLIMITED,H5S_UNLIMITED};
    space=H5Screate_simple(2,dims,max);

    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[2]={def_chunk(m.size1()),def_chunk(m.size2())};
    int status2=H5Pset_chunk(dcpl,2,chunk);

#ifdef O2SCL_HDF5_COMP    
    if (m.size1()*m.size2()>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::seti_mat_copy().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_STD_I32LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims[2];
    int ndims=H5Sget_simple_extent_dims(space,dims,0);

    // Set error if this dataset is more than 1-dimensional
    if (ndims!=2) {
      O2SCL_ERR2("Tried to set a non-matrix dataset with a ",
		 "matrix in hdf_file::seti_mat_copy().",exc_einval);
    }

    // If necessary, extend the dataset
    if (m.size1()!=dims[0] || m.size2()!=dims[1]) {
      hsize_t new_dims[2]={m.size1(),m.size2()};
      int status3=H5Dset_extent(dset,new_dims);
    }
    
  }

  // Write the data 
  int status;
  status=H5Dwrite(dset,H5T_NATIVE_INT,H5S_ALL,
		  H5S_ALL,H5P_DEFAULT,d);
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }

  delete[] d;
      
  return 0;
}

int hdf_file::geti_mat_copy(std::string name, ubmatrix_int &m) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[2];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);

  if (ndims!=2) {
    O2SCL_ERR2("Dimensions of dataspace do not match ",
	       "a matrix in hdf_file::getd_mat().",exc_efailed);
  }
  
  m.resize(dims[0],dims[1]); 
  int *d=new int[m.size1()*m.size2()];

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,d);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::geti_mat_copy().",
	      exc_einval);
  }

  for(size_t i=0;i<m.size1();i++) {
    for(size_t j=0;j<m.size2();j++) {
      m(i,j)=d[i*m.size2()+j];
    }
  }
  delete[] d;
  
  status=H5Dclose(dset);

  return 0;
}

int hdf_file::setd_ten(std::string name, 
		       const o2scl::tensor<double,std::vector<double>,
		       std::vector<size_t> > &t) {

  // Check if the tensor object is valid
  t.is_valid();
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::setd_ten().",exc_efailed);
  }

  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
  
  size_t ndims=t.get_rank();
  hsize_t *dims=new hsize_t[ndims];
  for(size_t k=0;k<ndims;k++) {
    dims[k]=t.get_size(k);
  }
  
  // If it doesn't exist, create it
  if (dset<0) {

    // Create dims, chunk, and max arrays
    hsize_t *chunk=new hsize_t[ndims];
    hsize_t *max=new hsize_t[ndims];
    for(size_t k=0;k<ndims;k++) {
      max[k]=H5S_UNLIMITED;
      chunk[k]=def_chunk(t.get_size(k));
    }
    // Create the dataspace
    space=H5Screate_simple(ndims,dims,max);
    
    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    int status2=H5Pset_chunk(dcpl,ndims,chunk);
    
#ifdef O2SCL_HDF5_COMP    
    if (t.total_size()>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::setd_ten().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_IEEE_F64LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

    delete[] chunk;
    delete[] max;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims2[100];
    int ndims2=H5Sget_simple_extent_dims(space,dims2,0);
    if (ndims2!=((int)ndims)) {
      O2SCL_ERR2("Tried to set a tensor on top of a tensor of ",
		 "different rank in hdf_file::setd_ten().",exc_efailed);
    }

    // If necessary, extend the dataset
    int status3=H5Dset_extent(dset,dims);
    
  }

  // Write the data 
  vector<size_t> zero(ndims);
  for(size_t k=0;k<ndims;k++) zero[k]=0;
  const double *ptr=&t.get(zero);
  int status;
  status=H5Dwrite(dset,H5T_NATIVE_DOUBLE,H5S_ALL,
		  H5S_ALL,H5P_DEFAULT,ptr);
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }
  delete[] dims;
      
  return 0;
}

int hdf_file::seti_ten(std::string name, 
		       const o2scl::tensor<int,std::vector<int>,
		       std::vector<size_t> > &t) {
  
  // Check if the tensor object is valid
  t.is_valid();
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::setd_ten().",exc_efailed);
  }

  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
  
  size_t ndims=t.get_rank();
  hsize_t *dims=new hsize_t[ndims];
  for(size_t k=0;k<ndims;k++) {
    dims[k]=t.get_size(k);
  }
  
  // If it doesn't exist, create it
  if (dset<0) {

    // Create dims, chunk, and max arrays
    hsize_t *chunk=new hsize_t[ndims];
    hsize_t *max=new hsize_t[ndims];
    for(size_t k=0;k<ndims;k++) {
      max[k]=H5S_UNLIMITED;
      chunk[k]=def_chunk(t.get_size(k));
    }
    // Create the dataspace
    space=H5Screate_simple(ndims,dims,max);
    
    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    int status2=H5Pset_chunk(dcpl,ndims,chunk);
    
#ifdef O2SCL_HDF5_COMP    
    if (t.total_size()>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::setd_ten().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_STD_I32LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

    delete[] chunk;
    delete[] max;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims2[100];
    int ndims2=H5Sget_simple_extent_dims(space,dims2,0);
    if (ndims2!=((int)ndims)) {
      O2SCL_ERR2("Tried to set a tensor on top of a tensor of ",
		 "different rank in hdf_file::setd_ten().",exc_efailed);
    }

    // If necessary, extend the dataset
    int status3=H5Dset_extent(dset,dims);
    
  }

  // Write the data 
  vector<size_t> zero(ndims);
  for(size_t k=0;k<ndims;k++) zero[k]=0;
  const int *ptr=&t.get(zero);
  int status;
  status=H5Dwrite(dset,H5T_NATIVE_INT,H5S_ALL,
		  H5S_ALL,H5P_DEFAULT,ptr);
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }
  delete[] dims;
      
  return 0;
}

int hdf_file::set_szt_ten(std::string name, 
			  const o2scl::tensor<size_t,std::vector<size_t>,
			  std::vector<size_t> > &t) {
  
  // Check if the tensor object is valid
  t.is_valid();
  
  if (write_access==false) {
    O2SCL_ERR2("File not opened with write access ",
	       "in hdf_file::setd_ten().",exc_efailed);
  }

  hid_t dset, space, dcpl=0;
  bool chunk_alloc=false;

  H5E_BEGIN_TRY
    {
      // See if the dataspace already exists first
      dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
    } 
  H5E_END_TRY 
#ifdef O2SCL_NEVER_DEFINED
  {
  }
#endif
  
  size_t ndims=t.get_rank();
  hsize_t *dims=new hsize_t[ndims];
  for(size_t k=0;k<ndims;k++) {
    dims[k]=t.get_size(k);
  }
  
  // If it doesn't exist, create it
  if (dset<0) {

    // Create dims, chunk, and max arrays
    hsize_t *chunk=new hsize_t[ndims];
    hsize_t *max=new hsize_t[ndims];
    for(size_t k=0;k<ndims;k++) {
      max[k]=H5S_UNLIMITED;
      chunk[k]=def_chunk(t.get_size(k));
    }
    // Create the dataspace
    space=H5Screate_simple(ndims,dims,max);
    
    // Set chunk with size determined by def_chunk()
    dcpl=H5Pcreate(H5P_DATASET_CREATE);
    int status2=H5Pset_chunk(dcpl,ndims,chunk);
    
#ifdef O2SCL_HDF5_COMP    
    if (t.total_size()>=min_compr_size) {
      // Compression part
      if (compr_type==1) {
	int status3=H5Pset_deflate(dcpl,6);
      } else if (compr_type==2) {
	int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
      } else if (compr_type!=0) {
	O2SCL_ERR2("Invalid compression type in ",
		   "hdf_file::setd_ten().",exc_einval);
      }
    }
#endif

    // Create the dataset
    dset=H5Dcreate(current,name.c_str(),H5T_STD_U64LE,space,H5P_DEFAULT,
		   dcpl,H5P_DEFAULT);
    chunk_alloc=true;

    delete[] chunk;
    delete[] max;

  } else {
    
    // Get current dimensions
    space=H5Dget_space(dset);  
    hsize_t dims2[100];
    int ndims2=H5Sget_simple_extent_dims(space,dims2,0);
    if (ndims2!=((int)ndims)) {
      O2SCL_ERR2("Tried to set a tensor on top of a tensor of ",
		 "different rank in hdf_file::setd_ten().",exc_efailed);
    }

    // If necessary, extend the dataset
    int status3=H5Dset_extent(dset,dims);
    
  }

  // Write the data 
  vector<size_t> zero(ndims);
  for(size_t k=0;k<ndims;k++) zero[k]=0;
  const size_t *ptr=&t.get(zero);
  int status;
  if (std::numeric_limits<size_t>::digits==64) {
    status=H5Dwrite(dset,H5T_NATIVE_ULLONG,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,ptr);
  } else {
    status=H5Dwrite(dset,H5T_NATIVE_ULONG,H5S_ALL,
		    H5S_ALL,H5P_DEFAULT,ptr);
  }
  
  status=H5Dclose(dset);
  status=H5Sclose(space);
  if (chunk_alloc) {
    status=H5Pclose(dcpl);
  }
  delete[] dims;
      
  return 0;
}

int hdf_file::getd_ten(std::string name, 
		       o2scl::tensor<double,std::vector<double>,
		       std::vector<size_t> > &t) {
  
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[100];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  if (ndims<1 || ndims>100) {
    O2SCL_ERR2("Dimensions less than 1 or greater than 100 ",
	       "in hdf_file::getd_ten().",exc_efailed);
  }
  
  // Allocate new data
  t.resize(ndims,dims);
  
  // Get pointer to first element
  vector<size_t> zero(ndims);
  for(int k=0;k<ndims;k++) zero[k]=0;
  double *start=&t.get(zero);

  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,start);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::getd_ten().",
	      exc_einval);
  }
  
  status=H5Dclose(dset);

  // Check if the tensor object is valid
  t.is_valid();

  return 0;
}

int hdf_file::geti_ten(std::string name, 
		       o2scl::tensor<int,std::vector<int>,
		       std::vector<size_t> > &t) {
  
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[100];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  if (ndims<1 || ndims>100) {
    O2SCL_ERR2("Dimensions less than 1 or greater than 100 ",
	       "in hdf_file::getd_ten().",exc_efailed);
  }
  
  // Allocate new data
  t.resize(ndims,dims);
  
  // Get pointer to first element
  vector<size_t> zero(ndims);
  for(int k=0;k<ndims;k++) zero[k]=0;
  int *start=&t.get(zero);
  
  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,start);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::geti_ten().",
	      exc_einval);
  }
  
  status=H5Dclose(dset);

  // Check if the tensor object is valid
  t.is_valid();

  return 0;
}

int hdf_file::get_szt_ten(std::string name, 
			  o2scl::tensor<size_t,std::vector<size_t>,
			  std::vector<size_t> > &t) {
  
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
  
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[100];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);
  if (ndims<1 || ndims>100) {
    O2SCL_ERR2("Dimensions less than 1 or greater than 100 ",
	       "in hdf_file::getd_ten().",exc_efailed);
  }
  
  // Allocate new data
  t.resize(ndims,dims);
  
  // Get pointer to first element
  vector<size_t> zero(ndims);
  for(int k=0;k<ndims;k++) zero[k]=0;
  size_t *start=&t.get(zero);
  
  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_HSIZE,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,start);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::get_szt_ten().",
	      exc_einval);
  }
  
  status=H5Dclose(dset);

  // Check if the tensor object is valid
  t.is_valid();

  return 0;
}

int hdf_file::getd_vec_prealloc(std::string name, size_t n, double *d) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
      
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);

  // Make sure sizes match
  if (dims[0]!=n) {
    O2SCL_ERR2("Vector size in file doesn't match allocated ",
	       "size in hdf_file::getd_vec_prealloc().",
	       exc_einval);
  }
      
  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,d);
  if (status<0) {
    O2SCL_ERR("Could not read dataspace in hdf_file::getd_vec_prealloc().",
	      exc_einval);
  }

  status=H5Dclose(dset);
      
  return 0;
}

int hdf_file::geti_vec_prealloc(std::string name, size_t n, int *i) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
      
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[1];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);

  // Make sure sizes match
  if (dims[0]!=n) {
    O2SCL_ERR2("Vector size in file doesn't match allocated ",
	       "size in hdf_file::geti_vec_prealloc().",
	       exc_einval);
  }
      
  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,i);
      
  status=H5Dclose(dset);
      
  return 0;
}

int hdf_file::getd_mat_prealloc(std::string name, size_t n, 
				size_t m, double *d) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
      
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[2];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);

  // Make sure sizes match
  if (dims[0]!=n || dims[1]!=m) {
    O2SCL_ERR2("Matrix size in file doesn't match allocated ",
	       "size in hdf_file::getd_mat_prealloc().",
	       exc_einval);
  }
      
  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_DOUBLE,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,d);

  status=H5Dclose(dset);
      
  return 0;
}

int hdf_file::geti_mat_prealloc(std::string name, size_t n, 
				size_t m, int *i) {
      
  // See if the dataspace already exists first
  hid_t dset=H5Dopen(current,name.c_str(),H5P_DEFAULT);
      
  // Get space requirements, to make sure they coincide
  // with the size specified by the user
  hid_t space=H5Dget_space(dset);  
  hsize_t dims[2];
  int ndims=H5Sget_simple_extent_dims(space,dims,0);

  // Make sure sizes match
  if (dims[0]!=n || dims[1]!=m) {
    O2SCL_ERR2("Matrix size in file doesn't match allocated ",
	       "size in hdf_file::getd_mat_prealloc().",
	       exc_einval);
  }
      
  // Read the data
  herr_t status=H5Dread(dset,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,
			H5P_DEFAULT,i);
      
  status=H5Dclose(dset);
      
  return 0;
}

int hdf_file::find_object_by_type(std::string type, std::string &name,
                                  bool use_regex, int verbose) {
  std::vector<std::string> vs;
  iterate_parms ip={"",this,false,type,verbose,ip_name_from_type,
    use_regex,vs};
  H5Literate(get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
             0,iterate_func,&ip);
  if (ip.found) {
    name=ip.tname;
    return success;
  }
  return exc_enotfound;
}

int hdf_file::list_objects_by_type(std::string type, 
                                   std::vector<std::string> &vs,
                                   bool use_regex, int verbose) {
  iterate_parms ip={"",this,false,type,verbose,ip_name_list_from_type,
    use_regex,vs};
  H5Literate(get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
             0,iterate_func,&ip);
  if (ip.name_list.size()>0) {
    vs=ip.name_list;
    return success;
  }
  return exc_enotfound;
}

int hdf_file::find_object_by_pattern(std::string pattern, std::string &type,
                                     bool use_regex, int verbose) {
				     
  std::vector<std::string> vs;
  iterate_parms ip={pattern,this,false,"",verbose,ip_type_from_pattern,
    use_regex,vs};
  H5Literate(get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
             0,iterate_func,&ip);
  if (ip.found) {
    type=ip.type;
    return success;
  }
  return exc_enotfound;
}

int hdf_file::find_object_by_name(std::string name, std::string &type,
                                  bool use_regex, int verbose) {
  std::vector<std::string> vs;
  iterate_parms ip={name,this,false,"",verbose,ip_type_from_name,
    use_regex,vs};
  H5Literate(get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
             0,iterate_func,&ip);
  if (ip.found) {
    type=ip.type;
    return success;
  }
  return exc_enotfound;
}

void hdf_file::type_process(iterate_parms &ip, int mode, size_t ndims, 
			    hsize_t dims[100], hsize_t max_dims[100],
			    std::string base_type, std::string name) {
  
  if (ndims==1) {

    if (dims[0]==0) {

      // --------------------------------------------------------------
      // Empty case

#ifdef O2SCL_REGEX      
      if (mode==ip_filelist) {
	cout << base_type;
      } else {
        regex r(name);
	if (mode==ip_type_from_name && name==ip.tname) {
	  ip.type=base_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_type_from_pattern &&
                   regex_search(ip.tname,r)) {
	  ip.type=base_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_name_from_type && base_type==ip.type) {
	  ip.tname=name;
	  ip.found=true;
	  return;
	} else if (mode==ip_name_list_from_type && base_type==ip.type) {
          ip.name_list.push_back(name);
	  ip.found=true;
	  return;
	}
      }
#else
      if (mode==ip_filelist) {
	cout << base_type;
      } else {
	if (mode==ip_type_from_name && name==ip.tname) {
	  ip.type=base_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_type_from_pattern &&
		   fnmatch(name.c_str(),ip.tname.c_str(),0)==0) {
	  ip.type=base_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_name_from_type && base_type==ip.type) {
	  ip.tname=name;
	  ip.found=true;
	  return;
	} else if (mode==ip_name_list_from_type && base_type==ip.type) {
          ip.name_list.push_back(name);
	  ip.found=true;
	  return;
	}
      }
#endif
      
    } else if (max_dims[0]==H5S_UNLIMITED) {

      // --------------------------------------------------------------
      // Unlimited array case

#ifdef O2SCL_REGEX      
      if (mode==ip_filelist) {
	cout << base_type << "[" << dims[0] << "/inf]";
      } else {
	string curr_type;
	if (base_type=="char") {
	  curr_type="string";
	} else {
	  curr_type=base_type+"[]";
	}
        regex r(name);
	if (mode==ip_type_from_name && name==ip.tname) {
	  ip.type=curr_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_type_from_pattern &&
                   regex_search(ip.tname,r)) {
	  ip.type=curr_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_name_from_type && curr_type==ip.type) {
	  ip.tname=name;
	  ip.found=true;
          return;
	} else if (mode==ip_name_list_from_type && curr_type==ip.type) {
          ip.name_list.push_back(name);
	  ip.found=true;
	  return;
	}
      }
#else
      if (mode==ip_filelist) {
	cout << base_type << "[" << dims[0] << "/inf]";
      } else {
	string curr_type;
	if (base_type=="char") {
	  curr_type="string";
	} else {
	  curr_type=base_type+"[]";
	}
	if (mode==ip_type_from_name && name==ip.tname) {
	  ip.type=curr_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_type_from_pattern &&
		   fnmatch(name.c_str(),ip.tname.c_str(),0)==0) {
	  ip.type=curr_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_name_from_type && curr_type==ip.type) {
	  ip.tname=name;
	  ip.found=true;
	  return;
	} else if (mode==ip_name_list_from_type && curr_type==ip.type) {
          ip.name_list.push_back(name);
	  ip.found=true;
	  return;
	}
      }
#endif

    } else {

      // --------------------------------------------------------------
      // Fixed-length array case

#ifdef O2SCL_REGEX      
      if (mode==ip_filelist) {
	// Determine if it is a fixed-size vector or single value
	if (dims[0]==1 && max_dims[0]==1) {
	  cout << base_type;
	} else {
	  cout << base_type << "[" << dims[0] << "/" << max_dims[0]
	       << "]";
	}
      } else {
	string curr_type;
	if (dims[0]==1 && max_dims[0]==1) {
	  curr_type=base_type;
	} else {
	  if (base_type=="char") {
	    curr_type="string";
	  } else {
	    curr_type=base_type+"[]";
	  }
	}
        regex r(name);
	if (mode==ip_type_from_name && name==ip.tname) {
	  ip.type=curr_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_type_from_pattern &&
                   regex_search(ip.tname,r)) {
	  ip.type=curr_type;
	  ip.found=true;
	  return;
	} else if (mode==ip_name_from_type && curr_type==ip.type) {
	  ip.tname=name;
	  ip.found=true;
	  return;
	} else if (mode==ip_name_list_from_type && curr_type==ip.type) {
          ip.name_list.push_back(name);
	  ip.found=true;
	  return;
	}
      }
    }
#else
    if (mode==ip_filelist) {
      // Determine if it is a fixed-size vector or single value
      if (dims[0]==1 && max_dims[0]==1) {
        cout << base_type;
      } else {
        cout << base_type << "[" << dims[0] << "/" << max_dims[0]
             << "]";
      }
    } else {
      string curr_type;
      if (dims[0]==1 && max_dims[0]==1) {
        curr_type=base_type;
      } else {
        if (base_type=="char") {
          curr_type="string";
        } else {
          curr_type=base_type+"[]";
        }
      }
      if (mode==ip_type_from_name && name==ip.tname) {
        ip.type=curr_type;
        ip.found=true;
        return;
      } else if (mode==ip_type_from_pattern &&
                 fnmatch(name.c_str(),ip.tname.c_str(),0)==0) {
        ip.type=curr_type;
        ip.found=true;
        return;
      } else if (mode==ip_name_from_type && curr_type==ip.type) {
        ip.tname=name;
        ip.found=true;
        return;
      } else if (mode==ip_name_list_from_type && curr_type==ip.type) {
        ip.name_list.push_back(name);
        ip.found=true;
        return;
      }
    }
  }
#endif
    
} else {
    
  // --------------------------------------------------------------
  // Otherwise, if it has a rank of 2 or larger, then treat it as a
  // tensor-like object

#ifdef O2SCL_REGEX  
  if (mode==ip_filelist) {
    cout << base_type << "[";
    for(int i=0;i<((int)ndims)-1;i++) {
      if (max_dims[i]==H5S_UNLIMITED) {
        cout << dims[i] << "/inf,";
      } else {
        cout << dims[i] << "/" << max_dims[i] << ",";
      }
    }
    if (max_dims[ndims-1]==H5S_UNLIMITED) {
      cout << dims[ndims-1] << "/inf].";
    } else {
      cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "]";
    }
  } else {
    string curr_type=base_type;
    for(int i=0;i<((int)ndims);i++) {
      curr_type+="[]";
    }
    regex r(name);
    if (mode==ip_type_from_name && name==ip.tname) {
      ip.type=curr_type;
      ip.found=true;
      return;
    } else if (mode==ip_type_from_pattern &&
               regex_search(ip.tname,r)) {
      ip.type=curr_type;
      ip.found=true;
      return;
    } else if (mode==ip_name_from_type && curr_type==ip.type) {
      ip.tname=name;
      ip.found=true;
      return;
    } else if (mode==ip_name_list_from_type && curr_type==ip.type) {
      ip.name_list.push_back(name);
      ip.found=true;
      return;
    }
  }
#else
  if (mode==ip_filelist) {
    cout << base_type << "[";
    for(int i=0;i<((int)ndims)-1;i++) {
      if (max_dims[i]==H5S_UNLIMITED) {
        cout << dims[i] << "/inf,";
      } else {
        cout << dims[i] << "/" << max_dims[i] << ",";
      }
    }
    if (max_dims[ndims-1]==H5S_UNLIMITED) {
      cout << dims[ndims-1] << "/inf].";
    } else {
      cout << dims[ndims-1] << "/" << max_dims[ndims-1] << "]";
    }
  } else {
    string curr_type=base_type;
    for(int i=0;i<((int)ndims);i++) {
      curr_type+="[]";
    }
    if (mode==ip_type_from_name && name==ip.tname) {
      ip.type=curr_type;
      ip.found=true;
      return;
    } else if (mode==ip_type_from_pattern &&
               fnmatch(name.c_str(),ip.tname.c_str(),0)==0) {
      ip.type=curr_type;
      ip.found=true;
      return;
    } else if (mode==ip_name_from_type && curr_type==ip.type) {
      ip.tname=name;
      ip.found=true;
      return;
    } else if (mode==ip_name_list_from_type && curr_type==ip.type) {
      ip.name_list.push_back(name);
      ip.found=true;
      return;
    }
  }
#endif    
    
  // End of 'else' for 'if (ndims==1)'
 }
return;
}

void hdf_file::file_list(bool use_regex, int verbose) {
  
  iterate_parms ip={"",this,false,"",verbose,
    ip_filelist,use_regex};
  
  H5Literate(get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
	     0,iterate_func,&ip);
  return;
}

void hdf_file::copy(int verbose, hdf_file &hf2) {
  
  iterate_copy_parms ip={this,&hf2,verbose};
  
  H5Literate(get_current_id(),H5_INDEX_NAME,H5_ITER_NATIVE,
	     0,iterate_copy_func,&ip);
  return;
}

herr_t hdf_file::iterate_copy_func(hid_t loc, const char *name, 
				   const H5L_info_t *inf, void *op_data) {

  // Arrange parameters
  iterate_copy_parms *ip=(iterate_copy_parms *)op_data;
  hdf_file &hf=*(ip->hf);
  hdf_file &hf2=*(ip->hf2);
  int loc_verbose=ip->verbose;

  hid_t top=hf.get_current_id();

  H5O_info_t infobuf;

  // AWS, 3/12/2020: modified to support hdf5 v1.12 which deprecates
  // the 4 parameter form.
#ifdef O2SCL_HDF5_PRE_1_12
  herr_t status=H5Oget_info_by_name(loc,name,&infobuf,H5P_DEFAULT);
#else
  herr_t status=H5Oget_info_by_name(loc,name,&infobuf,H5O_INFO_ALL,
				    H5P_DEFAULT);
#endif

  // If it's a group
  if (infobuf.type==H5O_TYPE_GROUP) {

    // Open the group and see if it's an O2scl object
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);
    string otype;
    hf.gets_def_fixed("o2scl_type","",otype);
    hf.close_group(group);
    hf.set_current_id(top);

    std::string name2=name;
    
    // An o2scl group
    if (otype==((string)"string[]")) {
      vector<string> s;
      hf.gets_vec_copy(name,s);
      hf2.sets_vec_copy(name,s);
    } else if (otype==((string)"table_units")) {
      table_units<> t;
      hdf_input(hf,t,name2);
      hdf_output(hf2,t,name);
    } else if (otype==((string)"table")) {
      table<> t;
      hdf_input(hf,t,name2);
      hdf_output(hf2,t,name);
    } else if (otype==((string)"uniform_grid<double>")) {
      uniform_grid<double> ug;
      hdf_input(hf,ug,name2);
      hdf_output(hf2,ug,name);
    } else if (otype==((string)"tensor")) {
      tensor<> t;
      hf.getd_ten(name,t);
      hf2.setd_ten(name,t);
    } else if (otype==((string)"tensor<int>")) {
      tensor<int> t;
      hf.geti_ten(name,t);
      hf2.seti_ten(name,t);
    } else if (otype==((string)"tensor<size_t>")) {
      tensor<size_t> t;
      hf.get_szt_ten(name,t);
      hf2.set_szt_ten(name,t);
    } else if (otype==((string)"vector<contour_line>")) {
      vector<contour_line> vcl;
      hdf_input(hf,vcl,name2);
      hdf_output(hf2,vcl,name);
    } else if (otype==((string)"hist")) {
      hist h;
      hdf_input(hf,h,name2);
      hdf_output(hf2,h,name);
    } else if (otype==((string)"hist_2d")) {
      hist_2d h;
      hdf_input(hf,h,name2);
      hdf_output(hf2,((const hist_2d &)h),name);
    } else if (otype==((string)"tensor_grid")) {
      tensor_grid<> tg;
      hdf_input(hf,tg,name2);
      hdf_output(hf2,tg,name);
    } else {
      cout << "Non O2scl group \"" << name << "\"." << endl;
    }
    
  } else if (infobuf.type==H5O_TYPE_DATASET) {
    
    // Open data set
    hid_t dset=H5Dopen(loc,name,H5P_DEFAULT);

    // Get type information
    hid_t type_id=H5Dget_type(dset);
    hid_t nat_id=H5Tget_native_type(type_id,H5T_DIR_ASCEND);

    // Get filter information
    int compr=0;
    hid_t plist_id=H5Dget_create_plist(dset);
    int num_filters=H5Pget_nfilters(plist_id);
    for(int i=0;i<num_filters;i++) {
      size_t n_elements=0;
      unsigned flags, filter_info;
      H5Z_filter_t filter_type=H5Pget_filter2
	(plist_id,0,&flags,&n_elements,0,0,0,&filter_info);
      if (filter_type==H5Z_FILTER_DEFLATE) {
	compr=1;
      } else if (filter_type==H5Z_FILTER_SZIP) {
	compr=2;
      }
    }

    // Get dataspace information
    hid_t space_id = H5Dget_space(dset);
    hsize_t dims[100];
    hsize_t max_dims[100];
    int ndims=H5Sget_simple_extent_dims(space_id,dims,max_dims);

    if (H5Tequal(nat_id,H5T_NATIVE_CHAR)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      if (ndims==1) {
	if (dims[0]==0) {
	  //cout << "Name " << name << " copy char." << endl;
	  char ch;
	  hf.getc(name,ch);
	  hf2.setc(name,ch);
	} else if (max_dims[0]==H5S_UNLIMITED) {
	  //cout << "Name " << name << " copy string (unlim)." << endl;
	  std::string str;
	  hf.gets(name,str);
	  hf2.sets(name,str);
	} else { 
	  if (dims[0]==1 && max_dims[0]==1) {
	    //cout << "Name " << name << " copy char (2)." << endl;
	    char ch;
	    hf.getc(name,ch);
	    hf2.setc(name,ch);
	  } else {
	    //cout << "Name " << name << " copy string." << endl;
	    std::string str;
	    hf.gets(name,str);
	    hf2.sets(name,str);
	  }
	}
      } else if (ndims==2) {
	vector<string> vs;
	hf.gets_vec_copy(name,vs);
	hf2.sets_vec_copy(name,vs);
      } else {
	O2SCL_ERR("HDF5 copy char[][][]... failed.",o2scl::exc_einval);
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_SHORT)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy short failed.",o2scl::exc_einval);
    } else if (H5Tequal(nat_id,H5T_NATIVE_INT)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      if (ndims==1) {
	if (dims[0]==0) {
	  //cout << "Name " << name << " copy int." << endl;
	  int i;
	  hf.geti(name,i);
	  hf2.seti(name,i);
	} else if (max_dims[0]==H5S_UNLIMITED) {
	  //cout << "Name " << name << " copy int[] (unlim)." << endl;
	  vector<int> vi;
	  hf.geti_vec(name,vi);
	  hf2.seti_vec(name,vi);
	} else { 
	  if (dims[0]==1 && max_dims[0]==1) {
	    //cout << "Name " << name << " copy int (2)." << endl;
	    int i;
	    hf.geti(name,i);
	    hf2.seti(name,i);
	  } else {
	    //cout << "Name " << name << " copy int[]." << endl;
	    vector<int> vi;
	    hf.geti_vec(name,vi);
	    hf2.seti_vec(name,vi);
	  }
	}
      } else {
	tensor<int> ti;
	hf.geti_ten(name,ti);
	hf2.seti_ten(name,ti);
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_LONG)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy long failed.",o2scl::exc_einval);
    } else if (H5Tequal(nat_id,H5T_NATIVE_LLONG)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy llong failed.",o2scl::exc_einval);
    } else if (H5Tequal(nat_id,H5T_NATIVE_UCHAR)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy uchar failed.",o2scl::exc_einval);
    } else if (H5Tequal(nat_id,H5T_NATIVE_USHORT)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy ushort failed.",o2scl::exc_einval);
    } else if (H5Tequal(nat_id,H5T_NATIVE_UINT)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy uint failed.",o2scl::exc_einval);
    } else if (H5Tequal(nat_id,H5T_NATIVE_ULONG)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      if (ndims==1) {
	if (dims[0]==0) {
	  //cout << "Name " << name << " copy size_t." << endl;
	  size_t szt;
	  hf.get_szt(name,szt);
	  hf2.set_szt(name,szt);
	} else if (max_dims[0]==H5S_UNLIMITED) {
	  //cout << "Name " << name << " copy size_t[] (unlim)." << endl;
	  vector<size_t> vszt;
	  hf.get_szt_vec(name,vszt);
	  hf2.set_szt_vec(name,vszt);
	} else { 
	  if (dims[0]==1 && max_dims[0]==1) {
	    //cout << "Name " << name << " copy size_t (2)." << endl;
	    size_t szt;
	    hf.get_szt(name,szt);
	    hf2.set_szt(name,szt);
	  } else {
	    //cout << "Name " << name << " copy size_t[]." << endl;
	    vector<size_t> vszt;
	    hf.get_szt_vec(name,vszt);
	    hf2.set_szt_vec(name,vszt);
	  }
	}
      } else {
	tensor<size_t> tszt;
	hf.get_szt_ten(name,tszt);
	hf2.set_szt_ten(name,tszt);
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_ULLONG)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy ullong failed.",o2scl::exc_einval);
    } else if (H5Tequal(nat_id,H5T_NATIVE_FLOAT)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy float failed.",o2scl::exc_einval);
    } else if (H5Tequal(nat_id,H5T_NATIVE_DOUBLE)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      if (ndims==1) {
	if (dims[0]==0) {
	  //cout << "Name " << name << " copy double." << endl;
	  double d;
	  hf.getd(name,d);
	  hf2.setd(name,d);
	} else if (max_dims[0]==H5S_UNLIMITED) {
	  //cout << "Name " << name << " copy size_t[] (unlim)." << endl;
	  vector<double> vd;
	  hf.getd_vec(name,vd);
	  hf2.setd_vec(name,vd);
	} else { 
	  if (dims[0]==1 && max_dims[0]==1) {
	    //cout << "Name " << name << " copy double (2)." << endl;
	    double d;
	    hf.getd(name,d);
	    hf2.setd(name,d);
	  } else {
	    //cout << "Name " << name << " copy double[]." << endl;
	    vector<double> vd;
	    hf.getd_vec(name,vd);
	    hf2.setd_vec(name,vd);
	  }
	}
      } else {
	tensor<double> td;
	hf.getd_ten(name,td);
	hf2.setd_ten(name,td);
      }
    } else if (H5Tequal(nat_id,H5T_NATIVE_LDOUBLE)) {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy ldouble failed.",o2scl::exc_einval);
    } else {
      H5Sclose(space_id);
      H5Tclose(nat_id);
      H5Tclose(type_id);
      H5Dclose(dset);
      O2SCL_ERR("HDF5 copy other dataset failed.",o2scl::exc_einval);
    }
    
  } else if (infobuf.type==H5O_TYPE_NAMED_DATATYPE) {
    O2SCL_ERR("HDF5 copy named datatype failed.",o2scl::exc_einval);
  } else {
    O2SCL_ERR("HDF5 copy unexpected type failed.",o2scl::exc_einval);
  }

  return 0;
}

herr_t hdf_file::iterate_func(hid_t loc, const char *name, 
			      const H5L_info_t *inf, void *op_data) {

  // Arrange parameters
  iterate_parms *ip=(iterate_parms *)op_data;
  hdf_file &hf=*(ip->hf);
  int loc_verbose=ip->verbose;
  int mode=ip->mode;
  bool use_regex=ip->use_regex;

  hid_t top=hf.get_current_id();

  H5O_info_t infobuf;

  // AWS, 3/12/2020: modified to support hdf5 v1.12 which deprecates
  // the 4 parameter form.
#ifdef O2SCL_HDF5_PRE_1_12
  herr_t status=H5Oget_info_by_name(loc,name,&infobuf,H5P_DEFAULT);
#else
  herr_t status=H5Oget_info_by_name(loc,name,&infobuf,H5O_INFO_ALL,
				    H5P_DEFAULT);
#endif

  ip->found=false;

  // If it's a group
  if (infobuf.type==H5O_TYPE_GROUP) {

    if (loc_verbose>1) {
      std::cout << "iterate_func(), group:" << std::endl;
    }

    // Open the group and see if it's an O2scl object
    hid_t group=hf.open_group(name);
    hf.set_current_id(group);
    string otype;
    hf.gets_def_fixed("o2scl_type","",otype);
    hf.close_group(group);
    hf.set_current_id(top);

    if (use_regex) {
      if (otype.length()!=0) {
        string sname=name;
        regex r(sname);
        if (mode==ip_filelist) {
          cout << "O2scl object \"" << name << "\" of type " 
               << otype << "." << endl;
        } else if (mode==ip_type_from_name && name==ip->tname) {
          ip->type=otype;
          ip->found=true;
          return 1;
        } else if (mode==ip_type_from_pattern &&
                   regex_search(ip->tname,r)) {
          ip->type=otype;
          ip->found=true;
          return 1;
        } else if (mode==ip_name_from_type && otype==ip->type) {
          ip->tname=name;
          ip->found=true;
          return 1;
        } else if (mode==ip_name_list_from_type && otype==ip->type) {
          ip->name_list.push_back(name);
          ip->found=true;
          return 0;
        }
      } else {
        // Non-O2scl group
        if (mode==ip_filelist) {
          cout << "Group \"" << name << "\"." << endl;
        }
      }
    } else {
      if (otype.length()!=0) {
        if (mode==ip_filelist) {
          if (otype==((string)"string[]")) {
            cout << "O2scl object \"" << name << "\" of type " 
                 << otype << "." << endl;
          } else {
            cout << "O2scl object \"" << name << "\" of type " 
                 << otype << "." << endl;
          }
        } else if (mode==ip_type_from_name && name==ip->tname) {
          ip->type=otype;
          ip->found=true;
          return 1;
        } else if (mode==ip_type_from_pattern &&
                   fnmatch(name,ip->tname.c_str(),0)==0) {
          ip->type=otype;
          ip->found=true;
          return 1;
        } else if (mode==ip_name_from_type && otype==ip->type) {
          ip->tname=name;
          ip->found=true;
          return 1;
        } else if (mode==ip_name_list_from_type && otype==ip->type) {
          ip->name_list.push_back(name);
          ip->found=true;
          return 0;
        }
      } else {
        // Non-O2scl group
        if (mode==ip_filelist) {
          cout << "Group \"" << name << "\"." << endl;
        }
      }
    }
    
    if (loc_verbose>1) {
      std::cout << "iterate_func(), group end." << std::endl;
    }
    
  } else if (infobuf.type==H5O_TYPE_DATASET) {
    
    if (loc_verbose>1) {
      std::cout << "iterate_func(), dataset:" << std::endl;
    }
    
    if (mode==ip_filelist) {
      cout << "Dataset \"" << name << "\" of type ";
    }

    if (loc_verbose>1) {
      cout << "Reading dataset named " << name << endl;
    }

    // Open data set
    hid_t dset=H5Dopen(loc,name,H5P_DEFAULT);

    // Get type information
    hid_t type_id=H5Dget_type(dset);
    hid_t nat_id=H5Tget_native_type(type_id,H5T_DIR_ASCEND);

    // Get filter information
    int compr=0;
    hid_t plist_id=H5Dget_create_plist(dset);
    int num_filters=H5Pget_nfilters(plist_id);
    for(int i=0;i<num_filters;i++) {
      size_t n_elements=0;
      unsigned flags, filter_info;
      H5Z_filter_t filter_type=H5Pget_filter2
	(plist_id,0,&flags,&n_elements,0,0,0,&filter_info);
      if (filter_type==H5Z_FILTER_DEFLATE) {
	compr=1;
      } else if (filter_type==H5Z_FILTER_SZIP) {
	compr=2;
      }
    }

    // Output compression information
    if (loc_verbose>1) {
      if (compr==1) {
	cout << "Compressed with deflate." << endl;
      } else if (compr==2) {
	cout << "Compressed with szip." << endl;
      } else {
	cout << "Uncompressed." << endl;
      }
    }

    if (mode==ip_filelist) {
      if (compr==1) {
	cout << "(1) ";
      } else if (compr==2) {
	cout << "(2) ";
      }
    }
    
    //if (loc_verbose>1) {
    //std::cout << "iterate_func(), 1:" << std::endl;
    //}
    
    // Get dataspace information
    hid_t space_id = H5Dget_space(dset);
    hsize_t dims[100];
    hsize_t max_dims[100];
    int ndims=H5Sget_simple_extent_dims(space_id,dims,max_dims);
    
    //if (loc_verbose>1) {
    //std::cout << "iterate_func(), 2:" << std::endl;
    //}
    
    if (H5Tequal(nat_id,H5T_NATIVE_CHAR)) {

      //if (loc_verbose>1) {
      //std::cout << "iterate_func(), 3:" << std::endl;
      //}
      
      type_process(*ip,mode,ndims,dims,max_dims,"char",name);
      
      //if (loc_verbose>1) {
      //std::cout << "iterate_func(), 3b:" << std::endl;
      //}
      
      if (ndims==1 && mode==ip_filelist) {

        //if (loc_verbose>1) {
        //std::cout << "iterate_func(), 3c: " << ndims << " "
        //<< dims[0] << std::endl;
        //}

        // Allocate memory
        char *c=new char[dims[0]];
        
        // Read the data
        herr_t statusx=H5Dread(dset,H5T_NATIVE_CHAR,H5S_ALL,H5S_ALL,
                              H5P_DEFAULT,c);
        if (statusx<0) {
          O2SCL_ERR("Could not read dataspace in hdf_file::iterate_func().",
                    exc_einval);
        }
        
        // Copy to the string object
        std::string s="";
        for(size_t i=0;i<dims[0];i++) s+=c[i];
        
        // Delete char memory
        delete[] c;
        
        //if (loc_verbose>1) {
        //std::cout << "iterate_func(), 3c2:" << std::endl;
        //}

	//std::string s;
	//hf.gets(name,s);

        //if (loc_verbose>1) {
        //std::cout << "iterate_func(), 3d:" << std::endl;
        //}
        
	if (dims[0]==0) {
	  cout << " is empty.";
	} else if (dims[0]<20) {
	  cout << " with value \"" << s << "\".";
	} else {
	  cout << " with value \""
	       << s[0] << s[1] << s[2] << s[3] << s[4] << s[5]
	       << " ... "
	       << s[dims[0]-6] << s[dims[0]-5] << s[dims[0]-4]
	       << s[dims[0]-3] << s[dims[0]-2] << s[dims[0]-1]
	       << "\".";

        }
      }

      //if (loc_verbose>1) {
      //std::cout << "iterate_func(), 4: " << ip->found << std::endl;
      //}

      if (ip->found==true) return 1;
      
    } else if (H5Tequal(nat_id,H5T_NATIVE_SHORT)) {
      type_process(*ip,mode,ndims,dims,max_dims,"short",name);
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_INT)) {
      type_process(*ip,mode,ndims,dims,max_dims,"int",name);
      if (ndims==1 && mode==ip_filelist) {
	std::vector<int> iarr;
	hf.geti_vec(name,iarr);
	if (dims[0]==0) {
	  cout << " is empty";
	} else if (dims[0]==1) {
	  cout << " with value " << iarr[0];
	} else if (dims[0]==2) {
	  cout << " with value " << iarr[0] << ", " << iarr[1];
	} else if (dims[0]==3) {
	  cout << " with value " << iarr[0] << ", " << iarr[1] << ", "
	       << iarr[2];
	} else {
	  cout << " with value " << iarr[0] << ", " << iarr[1] << ", ..., "
	       << iarr[dims[0]-1];
	}
	cout << ".";
      }
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_LONG)) {
      type_process(*ip,mode,ndims,dims,max_dims,"long",name);
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_LLONG)) {
      type_process(*ip,mode,ndims,dims,max_dims,"llong",name);
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_UCHAR)) {
      type_process(*ip,mode,ndims,dims,max_dims,"uchar",name);
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_USHORT)) {
      type_process(*ip,mode,ndims,dims,max_dims,"ushort",name);
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_UINT)) {
      type_process(*ip,mode,ndims,dims,max_dims,"uint",name);
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_ULONG)) {
      type_process(*ip,mode,ndims,dims,max_dims,"size_t",name);
      if (ndims==1 && mode==ip_filelist) {
	std::vector<size_t> sarr;
	hf.get_szt_vec(name,sarr);
	if (dims[0]==0) {
	  cout << " is empty";
	} else if (dims[0]==1) {
	  cout << " with value " << sarr[0];
	} else if (dims[0]==2) {
	  cout << " with value " << sarr[0] << ", " << sarr[1];
	} else if (dims[0]==3) {
	  cout << " with value " << sarr[0] << ", " << sarr[1] << ", "
	       << sarr[2];
	} else {
	  cout << " with value " << sarr[0] << ", " << sarr[1] << ", ..., "
	       << sarr[dims[0]-1];
	}
	cout << ".";
      }
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_ULLONG)) {
      type_process(*ip,mode,ndims,dims,max_dims,"ullong",name);
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_FLOAT)) {
      type_process(*ip,mode,ndims,dims,max_dims,"float",name);
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_DOUBLE)) {
      type_process(*ip,mode,ndims,dims,max_dims,"double",name);
      if (ndims==1 && mode==ip_filelist) {
	std::vector<double> darr;
	hf.getd_vec(name,darr);
	if (dims[0]==0) {
	  cout << " is empty";
	} else if (dims[0]==1) {
	  cout << " with value " << darr[0];
	} else if (dims[0]==2) {
	  cout << " with value " << darr[0] << ", " << darr[1];
	} else {
	  cout << " with value \n\t" << darr[0] << ", ..., "
	       << darr[dims[0]-1];
	}
	cout << ".";
      }
      if (ip->found==true) return 1;
    } else if (H5Tequal(nat_id,H5T_NATIVE_LDOUBLE)) {
      type_process(*ip,mode,ndims,dims,max_dims,"ldouble",name);
      if (ip->found==true) return 1;
    } else {
      hid_t filetype=H5Dget_type(dset);
      size_t str_size=H5Tget_size(filetype);
      
      hsize_t dimsx[3];
      hid_t space=H5Dget_space(dset);
      int ndimsx=H5Sget_simple_extent_dims(space,dimsx,0);
      hid_t memtype=-1;
      if (ndimsx==1 && dimsx[0]==1) {
	memtype=H5Tcopy(H5T_C_S1);
      }
      if (memtype>0) {
	status=H5Tclose(memtype);
	status=H5Sclose(space);
	status=H5Tclose(filetype);
        if (use_regex) {
          regex r(name);
          if (mode==ip_type_from_name && name==ip->tname) {
            ip->type="char[fixed]";
            ip->found=true;
            return 1;
          } else if (mode==ip_type_from_pattern &&
                     regex_search(ip->tname,r)) {
            ip->type="char[fixed]";
            ip->found=true;
          }
        } else {
          if (mode==ip_type_from_name && name==ip->tname) {
            ip->type="char[fixed]";
            ip->found=true;
            return 1;
          } else if (mode==ip_type_from_pattern &&
                     fnmatch(name,ip->tname.c_str(),0)==0) {
            ip->type="char[fixed]";
            ip->found=true;
          }
        }
	if (mode==ip_filelist) {
	  if (str_size==0) {
	    cout << "Error. Fixed length-string with no space for null.";
	  } else {
	    cout << "fixed-length (" << str_size-1
		 << ") string with value \"";
	    std::string s;
	    hf.gets_fixed(name,s);
	    cout << s << "\".";
	  }
	}
      } else {
	if (mode==ip_filelist) {
	  cout << "<unknown>.";
	}
	status=H5Sclose(space);
	status=H5Tclose(filetype);
      }
    }
    if (mode==ip_filelist) {
      cout << endl;
    }

    if ((mode==ip_type_from_name || mode==ip_type_from_pattern) &&
	loc_verbose>1) {
      cout << "Value of found is " << ip->found << " and type "
	   << ip->type << endl;
    }
    
    if (loc_verbose>1) {
      std::cout << "iterate_func(), 5: " << endl;
    }
    
    H5Sclose(space_id);
    H5Tclose(nat_id);
    H5Tclose(type_id);
    H5Dclose(dset);
    
    if (loc_verbose>1) {
      std::cout << "iterate_func(), dataset end" << std::endl;
    }
    
  } else if (infobuf.type==H5O_TYPE_NAMED_DATATYPE) {
    if (mode==ip_filelist) {
      cout << "Named type \"" << name << "\"." << endl;
    }
  } else {
    if (mode==ip_filelist) {
      cout << "Unexpected HDF type. " << endl;
    }
  }
  
  return 0;
}

