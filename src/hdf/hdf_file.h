/*
  ───────────────────────────────────────────────────────────────────

  Copyright (C) 2006-2025, Andrew W. Steiner

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
#ifndef O2SCL_HDF_FILE_H
#define O2SCL_HDF_FILE_H

/** \file hdf_file.h
    \brief File defining \ref o2scl_hdf::hdf_file
*/
#include <limits>

#ifdef O2SCL_PLAIN_HDF5_HEADER
#include <hdf5.h>
#else
#ifdef O2SCL_LINUX
#include <hdf5/serial/hdf5.h>
#else
#include <hdf5.h>
#endif
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <o2scl/vector.h>
#include <o2scl/tensor.h>
#include <o2scl/format_float.h>

/** \brief The \o2 namespace for I/O with HDF
 */
namespace o2scl_hdf {

  /** \brief Store data in an \o2 compatible HDF5 file

      \verbatim embed:rst
      See also the :ref:`File I/O with HDF5` section of the \o2
      User's guide.
      \endverbatim

      The member functions which write or get data from an HDF file
      begin with either <tt>get</tt> or <tt>set</tt>.  Where
      appropriate, the next character is either \c c for character, \c
      d for double, \c f for float, or \c i for int. 
      
      By default, vectors and matrices are written to HDF files in a
      chunked format, so their length can be changed later as
      necessary. The chunk size is chosen in \ref def_chunk() to be
      the closest power of 10 to the current vector size. 

      All files not closed by the user are closed in the destructor,
      but the destructor does not automatically close groups.

      \note Currently, HDF I/O functions write data to HDF files
      assuming that \c int and \c float have 4 bytes, while \c size_t
      and \c double are 8 bytes. All output is done in little endian
      format. While <tt>get</tt> functions can read data with
      different sizes or in big endian format, the <tt>set</tt>
      functions cannot currently write data this way.

      \note It does make sense to write a zero-length vector to an HDF
      file if the vector does not have a fixed size in order to create
      a placeholder for future output. Thus the <tt>set_vec()</tt> and
      allow zero-length vectors and the <tt>set_arr()</tt> functions
      allow the <tt>size_t</tt> parameter to be zero, in which case
      the pointer parameter is ignored. The <tt>set_vec_fixed()</tt>
      and <tt>set_arr_fixed()</tt> functions do not allow this, and
      will throw an exception if sent a zero-length vector.
      
      \future This class opens all files in R/W mode, which may
      cause I/O problems in file systems. This needs to be
      fixed by allowing the user to open a read-only file. 
      (AWS: 3/16/18 I think this is fixed now.)
      \future The \o2 HDF functions do not always consistently choose
      between throwing \o2 exceptions and throwing HDF5 exceptions.
      Check and/or fix this.
      \future Automatically close groups, e.g. by storing hid_t's in a
      stack?
      \future Rewrite the _arr_alloc() functions so that they 
      return a shared_ptr?
      \future Move the code from the 'filelist' acol command here
      into hdf_file.
  */
  class hdf_file {

  public:
    
    typedef boost::numeric::ublas::vector<double> ubvector;
    typedef boost::numeric::ublas::matrix<double> ubmatrix;
    typedef boost::numeric::ublas::vector<int> ubvector_int;
    typedef boost::numeric::ublas::matrix<int> ubmatrix_int;

  protected:

    /// Properties (protected)
    //@{
    /// File ID
    hid_t file;

    /// True if a file has been opened
    bool file_open;
    
    /// Current file or group location
    hid_t current;
    
    /// If true, then the file has read and write access 
    bool write_access;
    //@}

    /// \name Other functions (protected)
    //@{
    /** \brief Default chunk size

	Choose the closest power of 10 which is greater than or equal
	to 10 and less than or equal to \f$ 10^6 \f$.
    */
    virtual hsize_t def_chunk(size_t n) {
      size_t ch=(size_t)((1.0+1.0e-12)*
			 pow(10.0,floor(log10(((double)n))+0.5)));
      if (ch<10) ch=10;
      if (ch>1000000) ch=1000000;
      return ch;
    }
    //@}
    
  public:

    /// \name Constructor and destructor
    //@{
    hdf_file();

    virtual ~hdf_file();
    //@}

    /// \name Other functions
    //@{
    /// If true, then the file has read and write access 
    bool has_write_access() {
      return write_access;
    }

    /** \brief List datasets and \o2 objects in the top-level
	of the file 
    */
    void file_list(bool use_regex=false, int verbose=0);

    /** \brief Create a copy of the current HDF5 file and place
        the copy in \c hf2
     */
    void copy(int verbose, hdf_file &hf2);
    //@}

    /// \name Compression properties
    //@{
    /// Compression type (support experimental)
    int compr_type;

    /// Minimum size to compress by default
    size_t min_compr_size;
    //@}

    /// \name Open and close files
    //@{
    /** \brief Open a file named \c fname

	If \c err_on_fail is \c true, this calls the error handler if
	opening the file fails (e.g. because the file does not exist).
	If \c err_on_fail is \c false and opening the file fails,
	nothing is done and the function returns the value \ref
	o2scl::exc_efilenotfound. If the open succeeds, this function
	returns \ref o2scl::success.
    */
    int open(std::string fname, bool write_access=false,
	     bool err_on_fail=true);
    
    /// Open a file named \c fname or create if it doesn't already exist
    void open_or_create(std::string fname);

    /// Close the file
    void close();
    //@}
    
    /// \name Manipulate ids
    //@{
    /// Get the current file id
    hid_t get_file_id();

    /// Set the current working id
    void set_current_id(hid_t cur);

    /// Retrieve the current working id
    hid_t get_current_id();
    //@}

    /** \name Simple get functions

	If the specified object is not found, the \o2 error handler
	will be called.
    */
    //@{
    /// Get a character named \c name
    int getc(std::string name, char &c);
    
    /// Get a double named \c name
    int getd(std::string name, double &d);
    
    /// Get a float named \c name
    int getf(std::string name, float &f);

    /// Get a integer named \c name
    int geti(std::string name, int &i);

    /// Get an unsigned integer named \c name
    int get_szt(std::string name, size_t &u);

    /** \brief Get a (either variable or fixed-length)
        string named \c name

	\note Strings are stored as character arrays and thus
	retrieving a string from a file requires loading the
	information from the file into a character array, and then
	copying it to the string. This will be slow for very long
	strings.
    */
    int gets(std::string name, std::string &s);

    /** \brief Get a variable length string named \c name
     */
    int gets_var(std::string name, std::string &s);
    
    /** \brief Get a fixed-length string named \c name
     */
    int gets_fixed(std::string name, std::string &s);

    /** \brief Get a fixed-length string named \c name with default 
	value \c s
    */
    int gets_def_fixed(std::string name, std::string def, std::string &s);
    //@}

    /// \name Simple set functions
    //@{
    /// Set a character named \c name to value \c c
    void setc(std::string name, char c);

    /// Set a double named \c name to value \c d
    void setd(std::string name, double d);

    /// Set a float named \c name to value \c f
    void setf(std::string name, float f);
    
    /// Set an integer named \c name to value \c i
    void seti(std::string name, int i);

    /// Set an unsigned integer named \c name to value \c u
    void set_szt(std::string name, size_t u);

    /** \brief Set a string named \c name to value \c s
	
	The string is stored in the HDF file as an extensible
	character array rather than a string.
    */
    void sets(std::string name, std::string s);

    /** \brief Set a fixed-length string named \c name to value \c s
	
	This function stores <tt>s</tt> as a fixed-length string
	in the HDF file. If a dataset named \c name is already
	present, then \c s must not be longer than the string
	length already specified in the HDF file.
    */
    void sets_fixed(std::string name, std::string s);
    //@}

    /// \name Generic floating point I/O
    //@{
    /** \brief Set a generic floating point named \c name to value \c f
     */
    template<class fp_t> int setfp_copy(std::string name, fp_t &f) {
      std::string s=o2scl::dtos(f,-1);
      sets(name,s);
      return 0;
    }
    
    /** \brief Set a generic floating point vector named \c name to
        value \c f
     */
    template<class vec_fp_t> int setfp_vec_copy(std::string name,
                                                vec_fp_t &f) {
      size_t n=f.size();
      std::vector<std::string> vs;
      for(size_t i=0;i<n;i++) {
        vs.push_back(o2scl::dtos(f[i],-1));
      }
      sets_vec_copy(name,vs);
      return 0;
    }
    
    /** \brief Set a generic floating point tensor
        named \c name to value \c f
     */
    template<class ten_fp_t> int setfp_ten_copy(std::string name,
                                                ten_fp_t &tf) {
      size_t n=tf.total_size();
      o2scl::tensor_string ts;
      std::vector<size_t> size(tf.get_rank());
      tf.get_size_arr(size);
      ts.resize(tf.get_rank(),size);
      for(size_t i=0;i<n;i++) {
        tf.unpack_index(i,size);
        ts.get(size)=o2scl::dtos(tf.get(size),-1);
      }
      sets_ten_copy(name,ts);
      return 0;
    }
    
    /** \brief Get a generic floating point named \c name

        \warning No checks are made to ensure that the stored
        precision matches the precision of the floating point which
        is used.
    */
    template<class fp_t> int getfp_copy(std::string name, fp_t &f) {
      std::string s;
      gets(name,s);
      f=std::stod(s);
      return 0;
    }

    /** \brief Get a long double named \c name

        \warning No checks are made to ensure that the stored
        precision matches the precision of the floating point which is
        used. Note that the precision of the long double type is also
        not platform-independent.
    */
    int getfp_copy(std::string name, long double &f) {
      std::string s;
      gets(name,s);
      f=stold(s);
      return 0;
    }

    /** \brief Get a boost multiprecision floating point named \c name
        (specialization for Boost multiprecision numbers)

        \warning No checks are made to ensure that the stored
        precision matches the precision of the floating point which
        is used.
     */
    template<size_t N> int getfp_copy(std::string name,
                                      boost::multiprecision::number<
                                      boost::multiprecision::cpp_dec_float<
                                      N> > &f) {
      std::string s;
      gets(name,s);
      f=f(s);
      return 0;
    }
    
    /** \brief Get a generic floating point vector named \c name

        \warning No checks are made to ensure that the stored
        precision matches the precision of the floating point which
        is used.
    */
    template<class vec_fp_t> int getfp_vec_copy(std::string name,
                                                vec_fp_t &f) {
      std::vector<std::string> vs;
      gets_vec_copy(name,vs);
      f.resize(vs.size());
      for(size_t i=0;i<vs.size();i++) {
        f[i]=std::stod(vs[i]);
      }
      return 0;
    }

    /** \brief Get a generic floating point vector named \c name
        (specialization for Boost multiprecision numbers)

        \warning No checks are made to ensure that the stored
        precision matches the precision of the floating point which
        is used.
    */
    template<size_t N> int getfp_vec_copy
    (std::string name,
     std::vector<boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<
     N> > > &f) {
      
      std::vector<std::string> vs;
      gets_vec_copy(name,vs);
      f.resize(vs.size());
      for(size_t i=0;i<vs.size();i++) {
        f=f(vs[i]);
      }
      return 0;
    }

    /** \brief Get a generic floating point tensor named \c name
        (specialization for Boost multiprecision numbers)

        \warning No checks are made to ensure that the stored
        precision matches the precision of the floating point which
        is used.
    */
    template<size_t N> int getfp_ten_copy
    (std::string name,
     o2scl::tensor<boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<
     N> >, std::vector<boost::multiprecision::number<
     boost::multiprecision::cpp_dec_float<
     N> > >,
     std::vector<size_t> > &tf) {

      o2scl::tensor_string ts;
      gets_ten_copy(name,ts);
      size_t n=ts.total_size();
      std::vector<size_t> size(ts.get_rank());
      const std::vector<size_t> &csize=ts.get_size_arr();
      size=csize;
      tf.resize(ts.get_rank(),size);
      for(size_t i=0;i<n;i++) {
        ts.unpack_index(i,size);
        tf.get(size)=ts.get(size);
      }
      return 0;
    }
    //@}
    
    /// \name Group manipulation
    //@{
    /** \brief Open a group relative to the location specified in 
	\c init_id

	\note In order to ensure that future objects are written to the
	newly-created group, the user must use set_current_id()
	using the newly-created group ID for the argument. 
    */
    hid_t open_group(hid_t init_id, std::string path);
    
    /** \brief Open a group relative to the current location

	\note In order to ensure that future objects are written to the
	newly-created group, the user must use set_current_id()
	using the newly-created group ID for the argument. 
    */
    hid_t open_group(std::string path); 

    /// Close a previously created group
    int close_group(hid_t group) {
      return H5Gclose(group);
    }
    //@}

    /** \name Vector get functions

	These functions automatically free any previously allocated
	memory in <tt>v</tt> and then allocate the proper space
	required to read the information from the HDF file.
    */
    //@{
    /// Get vector dataset and place data in \c v
    int getd_vec(std::string name, std::vector<double> &v);

    /** \brief Get vector dataset and place data in \c v

	This works with any vector class which has a 
	<tt>resize()</tt> method.
    */
    template<class vec_t>
      int getd_vec_copy(std::string name, vec_t &v) {
      std::vector<double> v2;
      int ret=getd_vec(name,v2);
      v.resize(v2.size());
      o2scl::vector_copy(v2.size(),v2,v);
      return ret;
    }

    /// Get vector dataset and place data in \c v
    int geti_vec(std::string name, std::vector<int> &v);
    
    /** \brief Get vector dataset and place data in \c v
    */
    template<class vec_int_t>
      int geti_vec_copy(std::string name, vec_int_t &v) {
      std::vector<int> v2;
      int ret=geti_vec(name,v2);
      v.resize(v2.size());
      for(size_t i=0;i<v2.size();i++) v[i]=v2[i];
      return ret;
    }
    
    /// Get vector dataset and place data in \c v
    int get_szt_vec(std::string name, std::vector<size_t> &v);
    
    /** \brief Get vector dataset and place data in \c v
    */
    template<class vec_size_t> 
      int get_szt_vec_copy(std::string name, vec_size_t &v) {
      std::vector<int> v2;
      int ret=geti_vec(name,v2);
      v.resize(v2.size());
      for(size_t i=0;i<v2.size();i++) v[i]=v2[i];
      return ret;
    }
    
    /** \brief Get a vector of strings named \c name and store it in \c s
     */
    int gets_vec_copy(std::string name, std::vector<std::string> &s);

    /** \brief Get a vector of a vector of strings named \c name
        and store it in \c s
     */
    int gets_vec_vec_copy(std::string name,
                          std::vector<std::vector<std::string>> &s);
    /** \brief Get a vector of a vector of strings named \c name
        and store it in \c s
     */
    int getd_vec_vec_copy(std::string name,
                          std::vector<std::vector<double>> &s);
    //@}

    /** \name Vector set functions

	These functions automatically write all of the vector elements
	to the HDF file, if necessary extending the data that is
	already present.
    */
    //@{
    /// Set vector dataset named \c name with \c v
    int setd_vec(std::string name, const std::vector<double> &v);

    /** \brief Set vector dataset named \c name with \c v

	This requires a copy before the vector is written to
	the file.
    */
    template<class vec_t>
      int setd_vec_copy(std::string name, const vec_t &v) {
      if (v.size()==0) {
	return setd_arr(name,0,0);
      }
      
      // We have to copy to an std::vector first
      std::vector<double> v2(v.size());
      o2scl::vector_copy(v.size(),v,v2);
      
      return setd_arr(name,v2.size(),&v2[0]);
    }

    /// Set vector dataset named \c name with \c v
    int seti_vec(std::string name, const std::vector<int> &v);
    
    /** \brief Set vector dataset named \c name with \c v
	
	This requires a copy before the vector is written to
	the file.
    */
    template<class vec_int_t> 
      int seti_vec_copy(std::string name, vec_int_t &v) {
      if (v.size()==0) {
	return seti_arr(name,0,0);
      }
      
      // We have to copy to an std::vector first
      std::vector<int> v2(v.size());
      vector_copy(v.size(),v,v2);
      
      return seti_arr(name,v2.size(),&v2[0]);
    }
    /// Set vector dataset named \c name with \c v
    int set_szt_vec(std::string name, const std::vector<size_t> &v);
    /** \brief Set vector dataset named \c name with \c v
	
	This requires a copy before the vector is written to
	the file.
    */
    template<class vec_size_t> 
      int set_szt_vec_copy(std::string name, const vec_size_t &v) {
      if (v.size()==0) {
	return set_szt_arr(name,0,0);
      }
      
      // We have to copy to an std::vector first
      std::vector<size_t> v2(v.size());
      vector_copy(v.size(),v,v2);
      
      return set_szt_arr(name,v2.size(),&v2[0]);
    }

    /** \brief Set a vector of strings named \c name

        \warning This function copies the data in the vector of strings
        to a new string before writing the data to the HDF5 file and
        thus may be less useful for larger vectors or vectors which
        contain longer strings.

        \devnote 
	String vectors are reformatted as a single character array, in
	order to allow each string to have different length and to
	make each string extensible. The size of the vector \c s is
	stored as an integer named <tt>nw</tt>.
    */
    int sets_vec_copy(std::string name, const std::vector<std::string> &s);
    
    /** \brief Set a vector of vectors of strings named \c name

        \warning This function copies the data in the vector of strings
        to a new string before writing the data to the HDF5 file and
        thus may be less useful for larger vectors or vectors which
        contain longer strings.

        \devnote 
	String vectors are reformatted as a single character array, in
	order to allow each string to have different length and to
	make each string extensible. The size of the vector \c s is
	stored as an integer named <tt>nw</tt>.
     */
    int sets_vec_vec_copy(std::string name,
                          const std::vector<std::vector<std::string>> &s);
    
    /** \brief Set a vector of vectors named \c name
     */
    int setd_vec_vec_copy(std::string name,
                          const std::vector<std::vector<double>> &vvd);
    //@}

    /** \name Matrix get functions

	These functions automatically free any previously allocated
	memory in <tt>m</tt> and then allocate the proper space
	required to read the information from the HDF file.
    */
    //@{
    /// Get matrix dataset and place data in \c m
    int getd_mat_copy(std::string name, ubmatrix &m);
    /// Get matrix dataset and place data in \c m
    int geti_mat_copy(std::string name, ubmatrix_int &m);
    //@}

    /** \name Matrix set functions

	These functions automatically write all of the vector elements
	to the HDF file, if necessary extending the data that is
	already present.
    */
    //@{
    /** \brief Set matrix dataset named \c name with \c m
     */
    int setd_mat_copy(std::string name, const ubmatrix &m);

    /** \brief Set matrix dataset named \c name with \c m
     */
    int seti_mat_copy(std::string name, const ubmatrix_int &m);

    /** \brief Set a two-dimensional array dataset named \c name with \c m
     */
    template<class arr2d_t>
      int setd_arr2d_copy(std::string name, size_t r,
			  size_t c, const arr2d_t &a2d) {
      
      if (write_access==false) {
	O2SCL_ERR2("File not opened with write access ",
		   "in hdf_file::setd_arr2d_copy().",o2scl::exc_efailed);
      }
      
      // Copy to a C-style array
      double *d=new double[r*c];
      for(size_t i=0;i<r;i++) {
	for(size_t j=0;j<c;j++) {
	  d[i*c+j]=a2d[i][j];
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
	hsize_t dims[2]={r,c};
	hsize_t max[2]={H5S_UNLIMITED,H5S_UNLIMITED};
	space=H5Screate_simple(2,dims,max);
	
	// Set chunk with size determined by def_chunk()
	dcpl=H5Pcreate(H5P_DATASET_CREATE);
	hsize_t chunk[2]={def_chunk(r),def_chunk(c)};
	int status2=H5Pset_chunk(dcpl,2,chunk);
	
#ifdef O2SCL_HDF5_COMP    
	// Compression part
	if (compr_type==1) {
	  int status3=H5Pset_deflate(dcpl,6);
	} else if (compr_type==2) {
	  int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
	} else if (compr_type!=0) {
	  O2SCL_ERR2("Invalid compression type in ",
		     "hdf_file::setd_arr2d_copy().",o2scl::exc_einval);
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
		     "matrix in hdf_file::setd_arr2d_copy().",
		     o2scl::exc_einval);
	}
	
	// If necessary, extend the dataset
	if (r!=dims[0] || c!=dims[1]) {
	  hsize_t new_dims[2]={r,c};
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
      
      // Free the C-style array
      delete[] d;
      
      return 0;
    }

    /** \brief Set a two-dimensional array dataset named \c name with \c m
     */
    template<class arr2d_t>
      int seti_arr2d_copy(std::string name, size_t r,
			  size_t c, const arr2d_t &a2d) {
      
      if (write_access==false) {
	O2SCL_ERR2("File not opened with write access ",
		   "in hdf_file::seti_arr2d_copy().",o2scl::exc_efailed);
      }
      
      // Copy to a C-style array
      int *d=new int[r*c];
      for(size_t i=0;i<r;i++) {
	for(size_t j=0;j<c;j++) {
	  d[i*c+j]=a2d[i][j];
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
	hsize_t dims[2]={r,c};
	hsize_t max[2]={H5S_UNLIMITED,H5S_UNLIMITED};
	space=H5Screate_simple(2,dims,max);
	
	// Set chunk with size determined by def_chunk()
	dcpl=H5Pcreate(H5P_DATASET_CREATE);
	hsize_t chunk[2]={def_chunk(r),def_chunk(c)};
	int status2=H5Pset_chunk(dcpl,2,chunk);
	
#ifdef O2SCL_HDF5_COMP    
	// Compression part
	if (compr_type==1) {
	  int status3=H5Pset_deflate(dcpl,6);
	} else if (compr_type==2) {
	  int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
	} else if (compr_type!=0) {
	  O2SCL_ERR2("Invalid compression type in ",
		     "hdf_file::seti_arr2d_copy().",o2scl::exc_einval);
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
		     "matrix in hdf_file::seti_arr2d_copy().",
		     o2scl::exc_einval);
	}
	
	// If necessary, extend the dataset
	if (r!=dims[0] || c!=dims[1]) {
	  hsize_t new_dims[2]={r,c};
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
      
      // Free the C-style array
      delete[] d;
      
      return 0;
    }
    
    /** \brief Set a two-dimensional array dataset named \c name with \c m
     */
    template<class arr2d_t>
      int set_szt_arr2d_copy(std::string name, size_t r,
			     size_t c, const arr2d_t &a2d) {
      
      if (write_access==false) {
	O2SCL_ERR2("File not opened with write access ",
		   "in hdf_file::set_szt_arr2d_copy().",o2scl::exc_efailed);
      }
      
      // Copy to a C-style array
      size_t *d=new size_t[r*c];
      for(size_t i=0;i<r;i++) {
	for(size_t j=0;j<c;j++) {
	  d[i*c+j]=a2d[i][j];
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
	hsize_t dims[2]={r,c};
	hsize_t max[2]={H5S_UNLIMITED,H5S_UNLIMITED};
	space=H5Screate_simple(2,dims,max);
	
	// Set chunk with size determined by def_chunk()
	dcpl=H5Pcreate(H5P_DATASET_CREATE);
	hsize_t chunk[2]={def_chunk(r),def_chunk(c)};
	int status2=H5Pset_chunk(dcpl,2,chunk);
	
#ifdef O2SCL_HDF5_COMP    
	// Compression part
	if (compr_type==1) {
	  int status3=H5Pset_deflate(dcpl,6);
	} else if (compr_type==2) {
	  int status3=H5Pset_szip(dcpl,H5_SZIP_NN_OPTION_MASK,16);
	} else if (compr_type!=0) {
	  O2SCL_ERR2("Invalid compression type in ",
		     "hdf_file::set_szt_arr2d_copy().",o2scl::exc_einval);
	}
#endif
	
	// Create the dataset
	dset=H5Dcreate(current,name.c_str(),H5T_STD_U64LE,space,H5P_DEFAULT,
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
		     "matrix in hdf_file::set_szt_arr2d_copy().",
		     o2scl::exc_einval);
	}
	
	// If necessary, extend the dataset
	if (r!=dims[0] || c!=dims[1]) {
	  hsize_t new_dims[2]={r,c};
	  int status3=H5Dset_extent(dset,new_dims);
	}
	
      }
      
      // Write the data 
      int status;
      status=H5Dwrite(dset,H5T_NATIVE_HSIZE,H5S_ALL,
		      H5S_ALL,H5P_DEFAULT,d);
      
      status=H5Dclose(dset);
      status=H5Sclose(space);
      if (chunk_alloc) {
	status=H5Pclose(dcpl);
      }
      
      // Free the C-style array
      delete[] d;
      
      return 0;
    }
    //@}
    
    /// \name Tensor I/O functions
    //@{
    /** \brief Get a tensor of double-precision numbers from an HDF file

	This version does not require a full copy of the tensor.
    */
    int getd_ten(std::string name, 
		 o2scl::tensor<double,std::vector<double>,
		 std::vector<size_t> > &t);

    /** \brief Get a tensor of integers from an HDF file

	This version does not require a full copy of the tensor.
    */
    int geti_ten(std::string name, 
		 o2scl::tensor<int,std::vector<int>,
		 std::vector<size_t> > &t);

    /** \brief Get a tensor of size_t from an HDF file

	This version does not require a full copy of the tensor.
    */
    int get_szt_ten(std::string name, 
		 o2scl::tensor<size_t,std::vector<size_t>,
		 std::vector<size_t> > &t);

    /** \brief Get a tensor of double-precision numbers from an HDF file

	This version requires a full copy of the tensor from the
	HDF5 file into the \ref o2scl::tensor object.
    */
    template<class vec_t, class vec_size_t>
      int getd_ten_copy(std::string name,
			o2scl::tensor<double,vec_t,vec_size_t> &t) {
      o2scl::tensor<double,std::vector<double>,std::vector<size_t> > t2;
      int ret=getd_ten(name,t2);
      t=t2;
      return ret;
    }
		 
    /** \brief Get a tensor of integers from an HDF file

	This version requires a full copy of the tensor from the
	HDF5 file into the \ref o2scl::tensor object.
    */
    template<class vec_t, class vec_size_t>
      int geti_ten_copy(std::string name,
			o2scl::tensor<int,vec_t,vec_size_t> &t) {
      o2scl::tensor<int,std::vector<int>,std::vector<size_t> > t2;
      int ret=geti_ten(name,t2);
      t=t2;
      return ret;
    }
		 
    /** \brief Get a tensor of strings from an HDF file
        
	This version requires a full copy of the tensor from the
	HDF5 file into the \ref o2scl::tensor_base object.
    */
    int gets_ten_copy(std::string name, 
                      o2scl::tensor_base<std::string,
                      std::vector<std::string>,
                      std::vector<size_t>> &t);
    
    /** \brief Write a tensor of double-precision numbers to an HDF file

	You may overwrite a tensor already present in the
	HDF file only if it has the same rank. This version
	does not require a full copy of the tensor.
    */
    int setd_ten(std::string name, 
		 const o2scl::tensor<double,std::vector<double>,
		 std::vector<size_t> > &t);

    /** \brief Write a tensor of integers to an HDF file

	You may overwrite a tensor already present in the
	HDF file only if it has the same rank. This version
	does not require a full copy of the tensor.
    */
    int seti_ten(std::string name, 
		 const o2scl::tensor<int,std::vector<int>,
		 std::vector<size_t> > &t);

    /** \brief Write a tensor of size_t values to an HDF file

	You may overwrite a tensor already present in the
	HDF file only if it has the same rank. This version
	does not require a full copy of the tensor.
    */
    int set_szt_ten(std::string name, 
		    const o2scl::tensor<size_t,std::vector<size_t>,
		    std::vector<size_t> > &t);

    /** \brief Write a tensor of strings to an HDF file

        This stores all of the tensor strings in a monolithic
        character array. The string lengths are saved in a size_t
        tensor object.
    */
    int sets_ten_copy(std::string name,
                      const o2scl::tensor_base<std::string,
                      std::vector<std::string>,
                      std::vector<size_t>> &t);

    /** \brief Write a tensor of double-precision numbers to an HDF file

	You may overwrite a tensor already present in the
	HDF file only if it has the same rank. This version
	requires a full copy of the tensor from the \ref o2scl::tensor
	object into the HDF5 file.
    */
    template<class vec_t, class vec_size_t>
      int setd_ten_copy(std::string name, 
			const o2scl::tensor<double,std::vector<double>,
			std::vector<size_t> > &t) {
      o2scl::tensor<double,std::vector<double>,std::vector<size_t> > t2;
      t2=t;
      int ret=setd_ten(name,t2);
      return ret;
    }

    /** \brief Write a tensor of integers to an HDF file

	You may overwrite a tensor already present in the
	HDF file only if it has the same rank. This version
	requires a full copy of the tensor from the \ref o2scl::tensor
	object into the HDF5 file.
    */
    template<class vec_t, class vec_size_t>
      int seti_ten_copy(std::string name, 
			const o2scl::tensor<int,std::vector<int>,
			std::vector<size_t> > &t) {
      o2scl::tensor<int,std::vector<int>,std::vector<size_t> > t2;
      t2=t;
      int ret=seti_ten(name,t2);
      return ret;
    }
    //@}
    
    /** \name Array get functions

	All of these functions assume that the
	pointer allocated beforehand, and matches the size of the
	array in the HDF file. If the specified object is not found,
	the \o2 error handler will be called.
    */
    //@{
    /** \brief Get a character array named \c name of size \c n

	\note The pointer \c c must be allocated beforehand to 
	hold \c n entries, and \c n must match the size of the
	array in the HDF file.
    */
    int getc_arr(std::string name, size_t n, char *c);

    /** \brief Get a double array named \c name of size \c n 

	\note The pointer \c d must be allocated beforehand to 
	hold \c n entries, and \c n must match the size of the
	array in the HDF file.
    */
    int getd_arr(std::string name, size_t n, double *d);

    /** \brief Get a double array named \c name of size \c n 
	and put the compression type in \c compr

	\note The pointer \c d must be allocated beforehand to 
	hold \c n entries, and \c n must match the size of the
	array in the HDF file.
    */
    int getd_arr_compr(std::string name, size_t n, double *d,
		       int &compr);
    
    /** \brief Get a float array named \c name of size \c n 

	\note The pointer \c f must be allocated beforehand to 
	hold \c n entries, and \c n must match the size of the
	array in the HDF file.
    */
    int getf_arr(std::string name, size_t n, float *f);
    
    /** \brief Get an integer array named \c name of size \c n

	\note The pointer \c i must be allocated beforehand to 
	hold \c n entries, and \c n must match the size of the
	array in the HDF file.
    */
    int geti_arr(std::string name, size_t n, int *i);
    //@}
    
    /** \name Array get functions with memory allocation

        These functions allocate memory with \c new, which
	should be freed by the user with \c delete .
    */
    //@{
    /// Get a character array named \c name of size \c n
    int getc_arr_alloc(std::string name, size_t &n, char *c);

    /// Get a double array named \c name of size \c n 
    int getd_arr_alloc(std::string name, size_t &n, double *d);
    
    /// Get a float array named \c name of size \c n 
    int getf_arr_alloc(std::string name, size_t &n, float *f);
    
    /// Get an integer array named \c name of size \c n
    int geti_arr_alloc(std::string name, size_t &n, int *i);
    //@}
    
    /** \name Array set functions
     */
    //@{
    /// Set a character array named \c name of size \c n to value \c c
    int setc_arr(std::string name, size_t n, const char *c);
    
    /// Set a double array named \c name of size \c n to value \c d
    int setd_arr(std::string name, size_t n, const double *d);

    /// Set a float array named \c name of size \c n to value \c f
    int setf_arr(std::string name, size_t n, const float *f);

    /// Set a integer array named \c name of size \c n to value \c i
    int seti_arr(std::string name, size_t n, const int *i);

    /// Set a integer array named \c name of size \c n to value \c i
    int set_szt_arr(std::string name, size_t n, const size_t *u);
    //@}

    /** \name Fixed-length array set functions
	
	If a dataset named \c name is already present, then the
	user-specified array must not be longer than the array already
	present in the HDF file.
    */
    //@{
    /// Set a character array named \c name of size \c n to value \c c
    int setc_arr_fixed(std::string name, size_t n, const char *c);

    /// Set a double array named \c name of size \c n to value \c d
    int setd_arr_fixed(std::string name, size_t n, const double *c);

    /// Set a float array named \c name of size \c n to value \c f
    int setf_arr_fixed(std::string name, size_t n, const float *f);
    
    /// Set an integer array named \c name of size \c n to value \c i
    int seti_arr_fixed(std::string name, size_t n, const int *i);
    //@}
        
    /** \name Get functions with default values

	If the requested dataset is not found in the HDF file,
	the object is set to the specified default value
	and the error handler is not called.
    */
    //@{
    /// Get a character named \c name
    int getc_def(std::string name, char def, char &c);
    
    /// Get a double named \c name
    int getd_def(std::string name, double def, double &d);
    
    /// Get a float named \c name
    int getf_def(std::string name, float def, float &f);

    /// Get a integer named \c name
    int geti_def(std::string name, int def, int &i);

    /// Get a size_t named \c name
    int get_szt_def(std::string name, size_t def, size_t &i);

    /// Get a string named \c name
    int gets_def(std::string name, std::string def, std::string &s);

    /// Get a variable length string named \c name
    int gets_var_def(std::string name, std::string def, std::string &s);
    //@}

    /** \name Get functions with pre-allocated pointer
     */
    //@{
    /// Get a double array \c d pre-allocated to have size \c n
    int getd_vec_prealloc(std::string name, size_t n, double *d);

    /// Get an integer array \c i pre-allocated to have size \c n
    int geti_vec_prealloc(std::string name, size_t n, int *i);

    /// Get a double matrix \c d pre-allocated to have size <tt>(n,m)</tt>
    int getd_mat_prealloc(std::string name, size_t n, size_t m, double *d);

    /// Get an integer matrix \c i pre-allocated to have size <tt>(n,m)</tt>
    int geti_mat_prealloc(std::string name, size_t n, size_t m, int *i);
    //@}
    
    /// \name Find an object by type or name
    //@{
    /** \brief Look in hdf_file \c hf for an \o2 object of type \c
	type and if found, set \c name to the associated object
	name

	This function returns 0 if an object of type \c type is found
	and \ref o2scl::exc_enoprog if it fails.
    */
    int find_object_by_type(std::string type, std::string &name,
                            bool use_regex=false, int verbose=0);

    /** \brief Find all objects in hdf_file \c hf of type \c
	type and store the names in \c vs
        
	This function returns 0 if an object of type \c type is found
	and \ref o2scl::exc_enoprog if it fails.
    */
    int list_objects_by_type(std::string type,
                             std::vector<std::string> &vs,
                             bool use_regex=false, int verbose=0);

    /** \brief Look in hdf_file \c hf for an \o2 object with name 
	\c name and if found, set \c type to the associated type
	
	This function returns 0 if an object with name \c name is
	found and \ref o2scl::exc_enoprog if it fails.
    */
    int find_object_by_name(std::string name, std::string &type,
                            bool use_regex=false, int verbose=0);
			    
    
    /** \brief Look in hdf_file \c hf for an \o2 object with name 
	which matches a regular expression

	If an object is found, \c type is set to the associated type.
	This function returns 0 if an object with name \c name is
	found and \ref o2scl::exc_enoprog if it fails.
     */
    int find_object_by_pattern(std::string name, std::string &type,
                               bool use_regex=false, int verbose=0);
			       
    //@}

    /// \name Mode values for iteration functions
    //@{
    /// Used for <tt>acol -filelist</tt>
    static const int ip_filelist=1;
    /// Get name of HDF5 object from type
    static const int ip_name_from_type=2;
    /// Get type of HDF5 object with name 
    static const int ip_type_from_name=3;
    /// Get types of HDF5 object from names matching pattern
    static const int ip_type_from_pattern=4;
    /// Get list of names given type
    static const int ip_name_list_from_type=5;
    //@}

    /// \name Functions and structures for iteration over HDF5 objects
    //@{
    /// Parameters for iterate_func()
    typedef struct {
      /// Object name
      std::string tname;
      /// Pointer to HDF5 file
      o2scl_hdf::hdf_file *hf;
      /// True if found
      bool found;
      /// Object type
      std::string type;
      /// Verbose parameter
      int verbose;
      /** \brief Iteration mode, either \ref hdf_file::ip_filelist,
          \ref hdf_file::ip_name_from_type, \ref
          hdf_file::ip_type_from_name or \ref
          hdf_file::ip_type_from_pattern
       */
      int mode;
      /// If true, then use regex to match names
      bool use_regex;
      /// The list of names, used by \ref hdf_file::list_objects_by_type()
      std::vector<std::string> name_list;
    } iterate_parms;

    /// Parameters for iterate_copy_func()
    typedef struct {
      /// Pointer to source HDF5 file
      o2scl_hdf::hdf_file *hf;
      /// Pointer to destination HDF5 file
      o2scl_hdf::hdf_file *hf2;
      /// Verbosity parameter
      int verbose;
    } iterate_copy_parms;
    
    /// Process a type for \ref iterate_func() 
    static void type_process(iterate_parms &ip, int mode, size_t ndims, 
			     hsize_t dims[100], hsize_t max_dims[100],
			     std::string base_type, std::string name);
    
    /// HDF object iteration function
    static herr_t iterate_func(hid_t loc, const char *name, 
			       const H5L_info_t *inf, void *op_data);

    /// HDF5 object iteration function when copying
    static herr_t iterate_copy_func(hid_t loc, const char *name, 
				    const H5L_info_t *inf, void *op_data);
    //@}
    
  private:

    /// Copy constructor (private and empty)
    hdf_file(const hdf_file &);

    /// Copy constructor (private and empty)
    hdf_file& operator=(const hdf_file&);

  };

}

#endif
