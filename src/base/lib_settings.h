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
#ifndef O2SCL_LIB_SETTINGS_H
#define O2SCL_LIB_SETTINGS_H
#include <iostream>
#include <string>

/** \file lib_settings.h
    \brief Library settings class and global settings object
*/

#include <o2scl/convert_units.h>
#include <o2scl/find_constants.h>
#include <o2scl/rng.h>
#include <o2scl/set_python.h>

#ifdef O2SCL_MPI
#include <mpi.h>
#endif

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

    /// Return true if \o2 was installed with module support
    bool module_support();

    /// Return true if \o2 was installed with mpfr support
    bool mpfr_support();

    /// Return true if \o2 was installed with ncurses support
    bool ncurses_support();

    /// Return true if \o2 was installed with Armadillo support
    bool armadillo_support();

    /// Return true if \o2 was installed with Eigen support
    bool eigen_support();

    /// Return true if \o2 was installed with FFTW support
    bool fftw_support();

    /// Return true if \o2 was installed with CUBATURE support
    bool cubature_support();

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

    /// \name Unit conversion objects
    //@{
    /// Default convert_units object
    convert_units<double> def_cu;
    
    /// Get the global convert_units object
    convert_units<double> &get_convert_units() {
      return *cup;
    }
    //@}

    /** \brief Seed for thread-safe random number generation 
        (see \ref rng_set_seed() and \ref rng_set_seed_mpi())
     */
    unsigned int seed;
    
    // AWS: 2/22/21: I was originally thinking of using this to
    // control OpenMP threads, but I think for now the best is just to
    // use export OMP_NUM_THREADS to control this
    //size_t omp_num_threads;

    /// True if Python has been initialized (default false)
    bool py_initialized;

    /** \brief Initialize the python interface

        This function uses the environment variable O2SCL_PYTHON_EXT
        to set the Python name via the function Py_SetProgramName()
        before the function Py_Initialize() is called.
     */
    int py_init_nothrow(int verbose=0);

    /** \brief Initialize the python interface

        This function uses the environment variable O2SCL_PYTHON_EXT
        to set the Python name via the function Py_SetProgramName()
        before the function Py_Initialize() is called.
     */
    void py_init(int verbose=0);
    
    /// Finalize the python interface
    int py_final_nothrow(int verbose=0);

    /// Finalize the python interface
    void py_final(int verbose=0);

    /// String containing python version
    std::string py_version();
    
    /// Import arrays for numpy C API
    void *py_import_array();
    
    /// Get the path for the Python module named \c module
    std::string py_get_module_path(std::string module);
    
    /// Add path \c path to the python system search path
    void add_python_path(std::string path, int verbose=0);

    /// Get the python path and place the path strings in \c vs.
    void get_python_path(std::vector<std::string> &vs, int verbose=0);

#ifdef O2SCL_SET_PYTHON
    /** \brief Import module named \c module

        This function outputs some debugging information if the
        module import fails. The function py_init() must be 
        called first.

	This function is only defined if Python support is enabled.
     */
    PyObject *py_import_module(std::string module, int verbose=0);
#endif
    
  protected:

    /// \name Internal data set in the constructor [protected]
    //@{
    /// The present data directory
    std::string data_dir;

    /// The present documentation directory
    std::string doc_dir;

    /// Pointer to current \ref convert_units object
    convert_units<double> *cup;

    /// Pointer to current \ref find_constants object
    find_constants<> *fcp;
    //@}

  };

  /** \brief The global library settings object

      This global object is used by \o2 classes to store the global
      constant database, the global unit conversion database, and the
      directory for \o2p and \o2e data files. It may also be used by
      the end-user to probe details of the \o2 installation.
  */
  extern lib_settings_class o2scl_settings;
  
  /** \brief Find constant named \c name with unit \c unit and
      return the associated value
  */
  template<class fp_t=double>
  fp_t find_constant(std::string name, std::string unit) {
    o2scl::convert_units<fp_t> &cu=o2scl_settings.get_convert_units();
    return cu.find_unique(name,unit);
  }

  /** \brief Set the RNG seed using an OpenMP critical block

      This function sets the random number seed for the given random
      number generator. If OpenMP support is included when O2scl is
      compiled, then the seed generation is enclosed in an OpenMP
      critical block to ensure that multiple threads do not get the
      same seed. This function uses the variable \ref
      lib_settings_class::seed to store a global library seed.

      When OpenMP is enabled, it is possible that this function could 
      create a condition where no OpenMP thread can proceed and thus
      prevent subsequent code from running as normal, however, this 
      is expected to be rare in practice.
  */
  void rng_set_seed(rng<> &r, int mpi_size=1, int mpi_rank=0,
                    int verbose=1);
  
#if defined (O2SCL_MPI) || defined (DOXYGEN)
  
  /** \brief MPI version of rng_set_seed()
      
      \note This function requires that -DO2SCL_MPI is defined and
      that MPI_Init() has been called previously. It is defined as a
      template, but currently double is the only floating point type
      supported.

      This function is a MPI wrapper around \ref rng_set_seed().
   */
  template<class fp_t=double>
  void rng_set_seed_mpi(rng<fp_t> &r, int verbose=1) {
    int mpi_rank=0, mpi_size=1;
    // Get MPI rank, etc.
    MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
    return rng_set_seed(r,mpi_size,mpi_rank,verbose);
  }
#endif
  
}

extern "C" {
  void *o2scl_get_o2scl_settings();
}

#endif
