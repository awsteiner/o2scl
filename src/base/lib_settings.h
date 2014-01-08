/*
  -------------------------------------------------------------------
  
  Copyright (C) 2006-2014, Andrew W. Steiner
  
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

/** \brief The namespace for linear algebra classes and functions

    See \ref linalg_section for more complete information
    about linear algebra in \o2. 

    This namespace documentation is in the file 
    <tt>src/base/lib_settings.h</tt>

    \future Move this documentation to the linalg directory.
*/
namespace o2scl_linalg {
}

/** \brief The namespace for linear algebra classes and functions 
    with operator()

    This namespace contains an identical copy of all the functions given 
    in the \ref o2scl_cblas namespace, but perform matrix indexing
    with <tt>[][]</tt> rather than <tt>(,)</tt>. See \ref o2scl_linalg
    for the function listing and documentation.

    See \ref linalg_section for more complete information
    about linear algebra in \o2. 

    This namespace documentation is in the file 
    <tt>src/base/lib_settings.h</tt>
*/
namespace o2scl_linalg_bracket {
}

#ifndef DOXYGEN_NO_O2NS
namespace o2scl {
#endif

  /** \brief A class to manage global library settings

      This class reports global \o2 settings such as the current
      version, whether or not sub-libraries were installed and what
      the current parent directory for \o2 data files is.

      A global object of this type is defined in
      <tt>lib_settings.h</tt> named \ref o2scl_settings .
  */
  class lib_settings_class {

  public:

    lib_settings_class();

    ~lib_settings_class();

    /** \brief Return the data directory */
    std::string get_data_dir() {
      return data_dir;
    }

    /** \brief Set the data directory */
    int set_data_dir(std::string dir) {
      data_dir=dir;
      return 0;
    }
    
    /// Return true if the EOS library was installed
    bool eos_installed();

    /// Return true if the particle library was installed
    bool part_installed();

    /// Return true if \o2 was installed with HDF support
    bool hdf_support();

    /// Return true if \o2 was installed with OpenMP support
    bool openmp_support();

    /// Return true if \o2 was installed with Armadillo support
    bool armadillo_support();

    /// Return true if \o2 was installed with Eigen support
    bool eigen_support();

    /// Return true if \o2 was installed with C++11 support
    bool cpp11_support();
    
    /** \brief Return true if range checking was turned on during 
	installation (default true)
    */
    bool range_check();

    /// Return the library version
    std::string o2scl_version();

    /// Report some of the settings from config.h
    void config_h_report();

    /// \name Miscellaneous config.h string properties
    //@{
    std::string o2scl_name();
    std::string o2scl_package();
    std::string o2scl_bugreport();
    std::string o2scl_string();
    std::string o2scl_tarname();
    //@}

    /// Default convert units object
    convert_units def_cu;
    
    /// Get the global convert_units object
    convert_units &get_convert_units() {
      return *cup;
    }

    /// Set the global convert_units object
    void set_convert_units(convert_units &cu) {
      cup=&cu;
    }

  protected:

#ifndef DOXYGEN_INTERNAL

    /// The present data directory
    std::string data_dir;

    /// Pointer to current \ref convert_units object
    convert_units *cup;

#endif
  
  };

  /** \brief The global library settings object

      This global object is used by \ref polylog and some of the \o2e
      classes to find data files. It may also be used by the end-user
      to probe details of the \o2 installation.
   */
  extern lib_settings_class o2scl_settings;

#ifndef DOXYGEN_NO_O2NS
}
#endif

#endif
