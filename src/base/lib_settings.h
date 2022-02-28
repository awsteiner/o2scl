/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2022, Andrew W. Steiner
  
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
#ifndef O2SCL_LIB_SETTINGS_H
#define O2SCL_LIB_SETTINGS_H
#include <iostream>
#include <string>

#include <o2scl/convert_units.h>
#include <o2scl/find_constants.h>

/** \file lib_settings.h
    \brief Library settings class and global settings object
*/

/** \brief The main \o2 namespace
    
    By default, all \o2 classes and functions which are not listed as
    being in one of \o2's smaller specialized namespaces are in this
    namespace. 

    \comment
    \htmlonly
    For a full list of all the 
    O<span style='position: relative; top: 0.3em; font-size: 0.8em'>2</span>scl 
    classes, see
    <a href="annotated.html">Data Structures</a>.
    \endhtmlonly
    \endcomment
    
    This namespace documentation is in the file 
    <tt>src/base/lib_settings.h</tt>
*/
/*
    \comment
    For a full list of all the
    O<span style='position: relative; top: 0.3em; font-size: 0.8em'>2</span>scl 
    global objects which are not
    classes, see <a href="globals.html">Files::Globals</a>.
    \endcomment
*/
namespace o2scl {
}

namespace o2scl {

  /** \brief A class to manage global library settings

      This class reports global \o2 settings such as the current
      version, whether or not sub-libraries were installed and what
      the current parent directory for \o2 data files is.

      A global object of this type is defined in
      <tt>lib_settings.h</tt> named \ref o2scl_settings .

      This class should not typically be instantiated by the end-user.
  */
  class lib_settings_class {

  public:

    lib_settings_class();

    ~lib_settings_class();

    /// \name Library settings
    //@{
    /** \brief Return the data directory */
    std::string get_data_dir() {
      return data_dir;
    }

    /** \brief Set the data directory */
    int set_data_dir(std::string dir) {
      data_dir=dir;
      return 0;
    }
    
    /** \brief Return the doc directory */
    std::string get_doc_dir() {
      return doc_dir;
    }

    /** \brief Set the doc directory */
    int set_doc_dir(std::string dir) {
      doc_dir=dir;
      return 0;
    }
    
    /// Return true if \o2 was installed with OpenMP support
    bool openmp_support();

    /// Return true if \o2 was installed with readline support
    bool readline_support();

    /// Return true if \o2 was installed with mpfr support
    bool mpfr_support();

    /// Return true if \o2 was installed with ncurses support
    bool ncurses_support();

    /// Return true if \o2 was installed with support for GSL V2.0+
    bool gsl2_support();

    /// Return true if \o2 was installed with Armadillo support
    bool armadillo_support();

    /// Return true if \o2 was installed with Eigen support
    bool eigen_support();

    /// Return true if \o2 was installed with FFTW support
    bool fftw_support();

    /// Return true if \o2 was installed with CUBATURE support
    bool cubature_support();

    /// Return true if \o2 was installed with POLYLOGARITHM support
    bool polylogarithm_support();

    // Return true if \o2 was installed with Python support
    bool python_support();

    /// Return true if \o2 was installed with HDF5 compression support
    bool hdf5_compression_support();

    /** \brief Return system type determined by autoconf

	Returns either "OSX", "Linux" or "unknown".
     */
    std::string system_type();
    
    /** \brief Return true if range checking was turned on during 
	installation (default true)
    */
    bool range_check();

    /// Return the time \o2 was compiled
    std::string time_compiled();

    /// Return the date \o2 was compiled
    std::string date_compiled();

    /// Return the library version
    std::string o2scl_version();

    /// Report some of the settings from config.h
    void config_h_report();

    /// Obtain HDF5 version from the library compiled with O2scl
    void hdf5_lib_version(unsigned &maj, unsigned &min, unsigned &rel);
    
    /// Obtain HDF5 version from the current headers
    void hdf5_header_version(unsigned &maj, unsigned &min, unsigned &rel);
    //@}
    
    /// \name Miscellaneous config.h string properties
    //@{
    std::string o2scl_name();
    std::string o2scl_package();
    std::string o2scl_bugreport();
    std::string o2scl_string();
    std::string o2scl_tarname();
    //@}

    /// \name Unit conversion and constant objects
    //@{
    /// Default convert_units object
    convert_units<double> def_cu;

    /// Get the global convert_units object
    convert_units<double> &get_convert_units() {
      return *cup;
    }

    /// Get the global convert_units object
    find_constants &get_find_constants() {
      return *fcp;
    }

    /// Set the global convert_units and find_constants objects
    void set_convert_units(convert_units<double> &cu) {
      cup=&cu;
      fcp=&cu.fc;
      return;
    }
    //@}

    // AWS: 2/22/21: I was originally thinking of using this to
    // control OpenMP threads, but I think for now the best is just to
    // use export OMP_NUM_THREADS to control this
    //size_t omp_num_threads;

    /// True if Python has been initialized (default false)
    bool py_initialized;

    /// Initialize the python interface
    int py_init_nothrow(int verbose=0);

    /// Initialize the python interface
    void py_init(int verbose=0);
    
    /// Finalize the python interface
    int py_final_nothrow(int verbose=0);

    /// Finalize the python interface
    void py_final(int verbose=0);

    /// Add path \c path to the python system search path
    void add_python_path(std::string path, int verbose=0);
    
  protected:

#ifndef DOXYGEN_INTERNAL

    /// \name Internal data set in the constructor [protected]
    //@{
    /// The present data directory
    std::string data_dir;

    /// The present documentation directory
    std::string doc_dir;

    /// Pointer to current \ref convert_units object
    convert_units<double> *cup;

    /// Pointer to current \ref find_constants object
    find_constants *fcp;
    //@}

#endif
  
  };

  /** \brief The global library settings object

      This global object is used by \o2 classes to store the global
      constant database, the global unit conversion database, and the
      directory for \o2p and \o2e data files. It may also be used by
      the end-user to probe details of the \o2 installation.
   */
  extern lib_settings_class o2scl_settings;

  /** \brief

      Note that this template is in lib_settings.h because of the
      hierarchy of header files which must be included
   */
  template<class fp_t> int function_to_double_nothrow
  (std::string s,
   fp_t &result, int verbose=0) {
    
    std::string s2;
    
    // Remove quotes and apostrophes
    for(size_t i=0;i<s.length();i++) {
      if (s[i]!='\"' && s[i]!='\'') {
        s2+=s[i];
      }
    }
    
    calc_utf8<fp_t> calc;
    
    int ret=calc.compile_nothrow(s2.c_str(),0);
    if (ret!=0) return ret;
    
    std::vector<std::u32string> vs=calc.get_var_list();
    
    // The o2scl class find_constants doesn't work for 
    // multiprecision, so we return a non-zero value instead
    if (vs.size()!=0) {

      // There are undefined constants
      return 1;
      
    } else {
      
      // No variables, so just evaluate
      int ret2=calc.eval_nothrow(0,result);
      if (ret2!=0) return ret2;
    }
    
    return 0;
  }
  
}

extern "C" {
  void *o2scl_get_o2scl_settings();
}

#endif
