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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <o2scl/lib_settings.h>
#include <o2scl/prev_commit.h>

#ifdef O2SCL_SET_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#endif

#ifdef O2SCL_PLAIN_HDF5_HEADER
#include <hdf5.h>
#else
#ifdef O2SCL_LINUX
#include <hdf5/serial/hdf5.h>
#else
#include <hdf5.h>
#endif
#endif

#include <o2scl/set_python.h>
#include <o2scl/set_openmp.h>
#include <o2scl/set_readline.h>
#include <o2scl/set_fftw.h>
#include <o2scl/set_mpfr.h>
#include <o2scl/vector_special.h>

#ifdef O2SCL_SET_OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace o2scl;

lib_settings_class o2scl::o2scl_settings;

lib_settings_class::lib_settings_class() {
  data_dir=((std::string)(O2SCL_DATA_DIR));
  doc_dir=((std::string)(O2SCL_DOC_DIR));
  cup=&def_cu;
  fcp=&def_cu.fc;
  
  // Default conversions are given here. Obviously GNU units is better
  // at handling these things, but it's nice to have some of the easy
  // conversions in by default rather than worrying about opening a
  // pipe, etc.
  def_cu.default_conversions();

  // AWS: 2/22/21: I was originally thinking of using this to
  // control OpenMP threads, but I think for now the best is just to
  // use export OMP_NUM_THREADS to control this
  //
  // When this is zero, O2scl functions just use the value returned by
  // omp_get_num_threads()
  //omp_num_threads=0;

  py_initialized=false;
  seed=0;
}

lib_settings_class::~lib_settings_class() {
}

int lib_settings_class::py_init_nothrow(int verbose) {
#ifdef O2SCL_SET_PYTHON
  
  /*
    4/30/24: Removed this, as setprogramname is deprecated.
    
    char *pe=getenv("O2SCL_PYTHON_EXE");
    std::string python_exe;
    if (pe) {
    python_exe=pe;
    // This conversion method works only for normal ASCII strings.
    std::wstring pe2=std::wstring(python_exe.begin(),python_exe.end());
    Py_SetProgramName(pe2.c_str());
    }
  */

  if (py_initialized) {
    return 0;
  }
  
  if (verbose>0) {
    cout << "Running Py_Initialize()." << endl;
  }
  Py_Initialize();
  if (verbose>0) {
    cout << "Finished Py_Initialize()." << endl;
  }
  if (Py_IsInitialized()) {
    py_initialized=true;
    return 0;
  }
#else
  return 2;
#endif
  return 1;
}

void lib_settings_class::py_init(int verbose) {

  int ret=py_init_nothrow(verbose);
  if (ret!=0) {
    O2SCL_ERR("Python initialization failed.",o2scl::exc_efailed);
  }
  return;
}

std::string lib_settings_class::py_get_module_path(std::string module) {
#ifdef O2SCL_SET_PYTHON
  PyObject *p_module;
  p_module=PyImport_ImportModule("o2sclpy");
  //cout << "module: " << p_module << endl;
  PyObject *p_path=PyObject_GetAttrString(p_module,"__path__");
  //cout << "path: " << p_path << endl;
  //cout << PyList_Check(p_path) << endl;
  PyObject *p_item=PyList_GetItem(p_path,0);
  const char *result=PyUnicode_AsUTF8(p_item);
  std::string s_res=result;
  return s_res;
#else
  return "";
#endif
}

void *lib_settings_class::py_import_array() {

#ifdef O2SCL_SET_PYTHON
  // AWS, 2/12/23: The import_array() function is a macro (?!) which
  // returns NULL if it fails, but then throws a python, not a C++
  // exception. This forces us to use "void *" for the return type to
  // this function. This function uses PyErr_Occurred() to throw a C++
  // exception if an error occurs, so this function either fails, or
  // returns 0.

  std::cout << "Running import_array()." << std::endl;
  import_array();
  if (PyErr_Occurred()) {
    O2SCL_ERR("Failed to import numpy Python module(s).",
              o2scl::exc_efailed);
  }
#endif
  
  return 0;
}

int lib_settings_class::py_final_nothrow(int verbose) {
#ifdef O2SCL_SET_PYTHON
  if (verbose>0) {
    cout << "Running Py_Finalize()." << endl;
  }
  Py_Finalize();
  if (verbose>0) {
    cout << "Finished Py_Finalize()." << endl;
  }
  py_initialized=false;
#else
  return 2;
#endif
  return 0;
}

void lib_settings_class::py_final(int verbose) {
  int ret=py_final_nothrow(verbose);
  if (ret!=0) {
    O2SCL_ERR("Python finalization failed.",o2scl::exc_efailed);
  }
  return;
}

bool lib_settings_class::range_check() {
#if !O2SCL_NO_RANGE_CHECK
  return true;
#else
  return false;
#endif
}

#ifdef O2SCL_SET_PYTHON
PyObject *lib_settings_class::py_import_module(std::string module,
                                               int verbose) {

  if (o2scl_settings.py_initialized==false) {
    py_init(verbose);
  }
    
  PyObject *pModule, *pName;
  
  // Get the Unicode name of the user-specified module
  if (verbose>0) {
    cout << "Getting unicode for module name()." << endl;
  }
  pName=PyUnicode_FromString(module.c_str());
  if (pName==0) {
    O2SCL_ERR2("Create module name failed in ",
              "funct_python::set_function().",o2scl::exc_efailed);
  }

  // Import the user-specified module
  if (verbose>0) {
    cout << "Importing module." << endl;
  }
  pModule=PyImport_Import(pName);
  if (pModule==0) {

    std::cerr << "───────────────────────────────────────────"
              << "─────────────────────────────" << std::endl;
    std::cerr << "funct_python::set_function(): Load module "
              << module << " failed.\n" << std::endl;

    if (o2scl_settings.py_initialized==false) {
      std::cerr << "It does not appear that Python was initialized using "
                << "the O2scl function\n  o2scl::o2scl_settings.py_init(), "
                << "but maybe it was done "
                << "somewhere else?\n" << std::endl;
    }

    std::cerr << "Sometimes this error occurs because the desired module is "
              << "in a virtual\n  environment that O2scl is not using. "
              << std::endl;
    // 5/29/24: this section had to be removed because of the
    // deprecation of the Python() setprogramname function (see above).
    /*
              << "This can be fixed by setting the\n  environment variable "
              << "O2SCL_PYTHON_EXE before calling\n  "
              << "o2scl::o2scl_settings.py_init()." << std::endl;
    char *pe=getenv("O2SCL_PYTHON_EXE");
    if (pe) {
      std::string pe2=pe;
      std::cerr << "  The value of O2SCL_PYTHON_EXE is " << pe2
                << " .\n" << std::endl;
    } else {
      std::cerr << "  The environment value O2SCL_PYTHON_EXE is not "
                << "currently set.\n" << std::endl;
    }
    */
    
    std::cerr << "Alternatively, sometimes this error occurs "
              << "because the module is not in\n  the "
              << "Python path, which is currently set to:\n  ";
    /// Get the Python path
    std::vector<std::string> vs;
    o2scl_settings.get_python_path(vs,0);
    for(size_t j=0;j<vs.size();j++) {
      if (j==0) {
        std::cerr << "  [" << vs[j] << "," << std::endl;
      } else {
        std::cerr << "   " << vs[j] << "," << std::endl;
      }
    }
    std::cerr << "  ]" << std::endl;
    std::cerr << "  You can modify the Python path using "
              << "the function\n  o2scl::o2scl_settings.add_python_path().\n"
              << std::endl;
    
    std::cerr << "───────────────────────────────────────────"
              << "─────────────────────────────" << std::endl;
    
    O2SCL_ERR2((((std::string)"Load module named ")+module+
                " failed in ").c_str(),
               "funct_python::set_function().",o2scl::exc_efailed);
  }

  Py_DECREF(pName);
  
  return pModule;
}
#endif

void lib_settings_class::get_python_path(std::vector<std::string> &vs,
                                         int verbose) {

#ifdef O2SCL_SET_PYTHON

  vs.clear();
  
  // Import the system module so that we can modify the search path
  // to import the module where the function is located
  if (verbose>0) {
    cout << "Importing sys." << endl;
  }
  PyObject *sys_mod=PyImport_ImportModule("sys");
  if (sys_mod==0) {
    O2SCL_ERR2("Import system module failed in",
               "funct_python::set_function().",
               o2scl::exc_efailed);
  }
  
  // Obtain the path and the number of elements 
  if (verbose>0) {
    cout << "Getting sys.path" << endl;
  }
  PyObject *sys_path=PyObject_GetAttrString(sys_mod,"path");
  if (sys_path==0) {
    O2SCL_ERR2("Obtain sys.path failed in",
               "funct_python::set_function().",
               o2scl::exc_efailed);
  }
  if (verbose>0) {
    cout << "Getting len(sys.path)" << endl;
  }
  Py_ssize_t path_size=PySequence_Size(sys_path);
  if (path_size==-1) {
    O2SCL_ERR2("Getting sys.path sequence size failed in",
               "funct_python::set_function().",
               o2scl::exc_efailed);
  }

  // Iterate through each element in the path, and compare with
  // the user-specified path
  
  bool path_found=false;
  
  for(int j=0;j<path_size;j++) {
    if (verbose>0) {
      cout << "Getting item." << endl;
    }
    PyObject *item=PySequence_GetItem(sys_path,j);
    if (item==0) {
      O2SCL_ERR2("Getting sequence item failed in",
                 "funct_python::set_function().",
                 o2scl::exc_efailed);
    }
    if (verbose>0) {
      cout << "Getting representation." << endl;
    }
    PyObject *item2=PyObject_Repr(item);
    if (item2==0) {
      O2SCL_ERR2("Getting sequence item unicode failed in",
                 "funct_python::set_function().",
                 o2scl::exc_efailed);
    }
    if (verbose>0) {
      cout << "Getting string." << endl;
    }
    PyObject *str=PyUnicode_AsEncodedString(item2,"utf-8","Error");
    if (str==0) {
      O2SCL_ERR2("Getting encoded sequence item failed in",
                 "funct_python::set_function().",
                 o2scl::exc_efailed);
    }
    if (verbose>0) {
      cout << "Getting C string." << endl;
    }
    const char *cstr=PyBytes_AS_STRING(str);
    if (cstr==0) {
      O2SCL_ERR2("Getting C string from sequence item failed in",
                 "funct_python::set_function().",
                 o2scl::exc_efailed);
    }
    
    if (verbose>0) {
      cout << "Converting " << cstr << " to std::string." << endl;
    }
    string cppstr=cstr;
    
    vs.push_back(cppstr);
    
    if (verbose>0) {
      cout << "Decref() for next item" << endl;
    }
    Py_DECREF(str);
    Py_DECREF(item2);
    Py_DECREF(item);
  }

#endif
  
  return;
}

void lib_settings_class::add_python_path(std::string path, int verbose) {
  
#ifdef O2SCL_SET_PYTHON
  
  // Import the system module so that we can modify the search path
  // to import the module where the function is located
  if (verbose>0) {
    cout << "Importing sys." << endl;
  }
  PyObject *sys_mod=PyImport_ImportModule("sys");
  if (sys_mod==0) {
    O2SCL_ERR2("Import system module failed in",
               "funct_python::set_function().",
               o2scl::exc_efailed);
  }
  
  // Obtain the path and the number of elements 
  if (verbose>0) {
    cout << "Getting sys.path" << endl;
  }
  PyObject *sys_path=PyObject_GetAttrString(sys_mod,"path");
  if (sys_path==0) {
    O2SCL_ERR2("Obtain sys.path failed in",
               "funct_python::set_function().",
               o2scl::exc_efailed);
  }
  if (verbose>0) {
    cout << "Getting len(sys.path)" << endl;
  }
  Py_ssize_t path_size=PySequence_Size(sys_path);
  if (path_size==-1) {
    O2SCL_ERR2("Getting sys.path sequence size failed in",
               "funct_python::set_function().",
               o2scl::exc_efailed);
  }

  // Iterate through each element in the path, and compare with
  // the user-specified path
  
  bool path_found=false;
  
  for(int j=0;j<path_size;j++) {
    if (verbose>0) {
      cout << "Getting item." << endl;
    }
    PyObject *item=PySequence_GetItem(sys_path,j);
    if (item==0) {
      O2SCL_ERR2("Getting sequence item failed in",
                 "funct_python::set_function().",
                 o2scl::exc_efailed);
    }
    if (verbose>0) {
      cout << "Getting representation." << endl;
    }
    PyObject *item2=PyObject_Repr(item);
    if (item2==0) {
      O2SCL_ERR2("Getting sequence item unicode failed in",
                 "funct_python::set_function().",
                 o2scl::exc_efailed);
    }
    if (verbose>0) {
      cout << "Getting string." << endl;
    }
    PyObject *str=PyUnicode_AsEncodedString(item2,"utf-8","Error");
    if (str==0) {
      O2SCL_ERR2("Getting encoded sequence item failed in",
                 "funct_python::set_function().",
                 o2scl::exc_efailed);
    }
    if (verbose>0) {
      cout << "Getting C string." << endl;
    }
    const char *cstr=PyBytes_AS_STRING(str);
    if (cstr==0) {
      O2SCL_ERR2("Getting C string from sequence item failed in",
                 "funct_python::set_function().",
                 o2scl::exc_efailed);
    }
    
    if (verbose>0) {
      cout << "Converting " << cstr << " to std::string." << endl;
    }
    string cppstr=cstr;
    if (cppstr==path) {
      path_found=true;
      if (verbose>0) {
        cout << "Path found." << endl;
      }
    }
    
    if (verbose>0) {
      cout << "Decref() for next item" << endl;
    }
    Py_DECREF(str);
    Py_DECREF(item2);
    Py_DECREF(item);
  }

  if (verbose>0) {
    cout << "import sys" << endl;
  }
  // Is this necessary?
  PyRun_SimpleString("import sys");

  // If necessary, add the user-specified path to sys.path
  if (path_found==false) {
    if (verbose>0) {
      cout << "Adding path" << endl;
    }
    std::string pys=((std::string)"sys.path.append('")+path+"')";
    PyRun_SimpleString(pys.c_str());
  } else {
    if (verbose>0) {
      cout << "Path found" << endl;
    }
  }

  // Free the memory associated with the system path and module
  if (verbose>0) {
    cout << "Decref" << endl;
  }
  Py_DECREF(sys_path);
  Py_DECREF(sys_mod);

  if (verbose>0) {
    cout << "Done in add_python_path()." << endl;
  }
  
#else

  O2SCL_ERR2("Python not supported in ",
             "lib_settings_class::add_python_path().",o2scl::exc_efailed);

#endif  
  
  return;
}

std::string lib_settings_class::o2scl_version() {
  return PACKAGE_VERSION;
}

std::string lib_settings_class::o2scl_name() {
  return PACKAGE_NAME;
}

std::string lib_settings_class::o2scl_package() {
  return PACKAGE;
}

std::string lib_settings_class::o2scl_bugreport() {
  return PACKAGE_BUGREPORT;
}

std::string lib_settings_class::o2scl_string() {
  return PACKAGE_STRING;
}

std::string lib_settings_class::o2scl_tarname() {
  return PACKAGE_TARNAME;
}

std::string lib_settings_class::date_compiled() {
#ifdef O2SCL_UBUNTU_PKG
  return "<package>";
#else
  return __DATE__;
#endif
}

std::string lib_settings_class::time_compiled() {
#ifdef O2SCL_UBUNTU_PKG
  return "<package>";
#else
  return __TIME__;
#endif
}

bool lib_settings_class::hdf5_compression_support() {
#ifdef O2SCL_HDF5_COMP
  return true;
#else
  return false;
#endif
}

void lib_settings_class::hdf5_header_version(unsigned &maj,
					     unsigned &min, unsigned &rel) {
  maj=H5_VERS_MAJOR;
  min=H5_VERS_MINOR;
  rel=H5_VERS_RELEASE;
  return;
}

void lib_settings_class::hdf5_lib_version(unsigned &maj,
					  unsigned &min, unsigned &rel) {
  H5get_libversion(&maj,&min,&rel);
  return;
}

bool lib_settings_class::armadillo_support() {
#ifdef O2SCL_SET_ARMA
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::eigen_support() {
#ifdef O2SCL_SET_EIGEN
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::fftw_support() {
#ifdef O2SCL_SET_FFTW
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::cubature_support() {
#ifdef O2SCL_CUBATURE
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::multiprecision_support() {
#ifdef O2SCL_SET_MULTIP
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::python_support() {
#ifdef O2SCL_SET_PYTHON
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::openmp_support() {
#ifdef O2SCL_SET_OPENMP
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::readline_support() {
#ifdef O2SCL_SET_READLINE
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::cuda_support() {
#ifdef O2SCL_CUDA
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::mpfr_support() {
#ifdef O2SCL_SET_MPFR
  return true;
#else
  return false;
#endif
}

bool lib_settings_class::ncurses_support() {
#ifdef O2SCL_NCURSES
  return true;
#else
  return false;
#endif
}

std::string lib_settings_class::system_type() {
#ifdef O2SCL_LINUX
  return ((string)"Linux");
#else
#ifdef O2SCL_OSX
  return ((string)"OSX");
#else
  return ((string)"unknown");
#endif
#endif
}

void lib_settings_class::config_h_report() {
  //cout << "Previous commit hash: " << O2SCL_PREV_COMMIT_HASH << endl;
  //cout << "Previous commit date: " << O2SCL_PREV_COMMIT_DATE << endl;
  cout << "Commit branch: " << O2SCL_BRANCH << endl;
#ifdef HAVE_ACOSH
  cout << "HAVE_ACOSH: " << HAVE_ACOSH << endl;
#else
  cout << "HAVE_ACOSH: <not defined>" << endl;
#endif
#ifdef HAVE_ASINH
  cout << "HAVE_ASINH: " << HAVE_ASINH << endl;
#else
  cout << "HAVE_ASINH: <not defined>" << endl;
#endif
#ifdef HAVE_ATANH
  cout << "HAVE_ATANH: " << HAVE_ATANH << endl;
#else
  cout << "HAVE_ATANH: <not defined>" << endl;
#endif
#ifdef HAVE_DLFCN_H
  cout << "HAVE_DLFCN_H: " << HAVE_DLFCN_H << endl;
#else
  cout << "HAVE_DLFCN_H: <not defined>" << endl;
#endif
#ifdef HAVE_FINITE
  cout << "HAVE_FINITE: " << HAVE_FINITE << endl;
#else
  cout << "HAVE_FINITE: <not defined>" << endl;
#endif
#ifdef HAVE_FLOOR
  cout << "HAVE_FLOOR: " << HAVE_FLOOR << endl;
#else
  cout << "HAVE_FLOOR: <not defined>" << endl;
#endif
#ifdef HAVE_INTTYPES_H
  cout << "HAVE_INTTYPES_H: " << HAVE_INTTYPES_H << endl;
#else
  cout << "HAVE_INTTYPES_H: <not defined>" << endl;
#endif
#ifdef HAVE_LIBCBLAS
  cout << "HAVE_LIBCBLAS: " << HAVE_LIBCBLAS << endl;
#else
  cout << "HAVE_LIBCBLAS: <not defined>" << endl;
#endif
#ifdef HAVE_LIBGSL
  cout << "HAVE_LIBGSL: " << HAVE_LIBGSL << endl;
#else
  cout << "HAVE_LIBGSL: <not defined>" << endl;
#endif
#ifdef HAVE_LIBGSLCBLAS
  cout << "HAVE_LIBGSLCBLAS: " << HAVE_LIBGSLCBLAS << endl;
#else
  cout << "HAVE_LIBGSLCBLAS: <not defined>" << endl;
#endif
#ifdef HAVE_LIBHDF5
  cout << "HAVE_LIBHDF5: " << HAVE_LIBHDF5 << endl;
#else
  cout << "HAVE_LIBHDF5: <not defined>" << endl;
#endif
#ifdef HAVE_LIBHDF5_HL
  cout << "HAVE_LIBHDF5_HL: " << HAVE_LIBHDF5_HL << endl;
#else
  cout << "HAVE_LIBHDF5_HL: <not defined>" << endl;
#endif
#ifdef HAVE_MALLOC
  cout << "HAVE_MALLOC: " << HAVE_MALLOC << endl;
#else
  cout << "HAVE_MALLOC: <not defined>" << endl;
#endif
#ifdef HAVE_MALLOC_H
  cout << "HAVE_MALLOC_H: " << HAVE_MALLOC_H << endl;
#else
  cout << "HAVE_MALLOC_H: <not defined>" << endl;
#endif
#ifdef HAVE_MEMORY_H
  cout << "HAVE_MEMORY_H: " << HAVE_MEMORY_H << endl;
#else
  cout << "HAVE_MEMORY_H: <not defined>" << endl;
#endif
#ifdef HAVE_POPEN
  cout << "HAVE_POPEN: " << HAVE_POPEN << endl;
#else
  cout << "HAVE_POPEN: <not defined>" << endl;
#endif
#ifdef HAVE_POW
  cout << "HAVE_POW: " << HAVE_POW << endl;
#else
  cout << "HAVE_POW: <not defined>" << endl;
#endif
#ifdef HAVE_SELECT
  cout << "HAVE_SELECT: " << HAVE_SELECT << endl;
#else
  cout << "HAVE_SELECT: <not defined>" << endl;
#endif
#ifdef HAVE_SQRT
  cout << "HAVE_SQRT: " << HAVE_SQRT << endl;
#else
  cout << "HAVE_SQRT: <not defined>" << endl;
#endif
#ifdef HAVE_STDBOOL_H
  cout << "HAVE_STDBOOL_H: " << HAVE_STDBOOL_H << endl;
#else
  cout << "HAVE_STDBOOL_H: <not defined>" << endl;
#endif
#ifdef HAVE_STDINT_H
  cout << "HAVE_STDINT_H: " << HAVE_STDINT_H << endl;
#else
  cout << "HAVE_STDINT_H: <not defined>" << endl;
#endif
#ifdef HAVE_STDLIB_H
  cout << "HAVE_STDLIB_H: " << HAVE_STDLIB_H << endl;
#else
  cout << "HAVE_STDLIB_H: <not defined>" << endl;
#endif
#ifdef HAVE_STRCHR
  cout << "HAVE_STRCHR: " << HAVE_STRCHR << endl;
#else
  cout << "HAVE_STRCHR: <not defined>" << endl;
#endif
#ifdef HAVE_STRINGS_H
  cout << "HAVE_STRINGS_H: " << HAVE_STRINGS_H << endl;
#else
  cout << "HAVE_STRINGS_H: <not defined>" << endl;
#endif
#ifdef HAVE_SYS_SELECT_H
  cout << "HAVE_SYS_SELECT_H: " << HAVE_SYS_SELECT_H << endl;
#else
  cout << "HAVE_SYS_SELECT_H: <not defined>" << endl;
#endif
#ifdef HAVE_SYS_SOCKET_H
  cout << "HAVE_SYS_SOCKET_H: " << HAVE_SYS_SOCKET_H << endl;
#else
  cout << "HAVE_SYS_SOCKET_H: <not defined>" << endl;
#endif
#ifdef HAVE_SYS_STAT_H
  cout << "HAVE_SYS_STAT_H: " << HAVE_SYS_STAT_H << endl;
#else
  cout << "HAVE_SYS_STAT_H: <not defined>" << endl;
#endif
#ifdef HAVE_SYS_TYPES_H
  cout << "HAVE_SYS_TYPES_H: " << HAVE_SYS_TYPES_H << endl;
#else
  cout << "HAVE_SYS_TYPES_H: <not defined>" << endl;
#endif
#ifdef HAVE_UNISTD_H
  cout << "HAVE_UNISTD_H: " << HAVE_UNISTD_H << endl;
#else
  cout << "HAVE_UNISTD_H: <not defined>" << endl;
#endif
#ifdef HAVE__BOOL
  cout << "HAVE__BOOL: " << HAVE__BOOL << endl;
#else
  cout << "HAVE__BOOL: <not defined>" << endl;
#endif
  return;
}

void *o2scl_get_o2scl_settings() {
  return &o2scl::o2scl_settings;
}

std::string lib_settings_class::py_version() {
#ifdef O2SCL_SET_PYTHON
  ostringstream oss;
  oss << PY_MAJOR_VERSION << " "
      << PY_MINOR_VERSION << " " << PY_MICRO_VERSION << " "
      << PY_RELEASE_LEVEL << " " << PY_RELEASE_SERIAL << " "
      << PY_VERSION;
  return oss.str();
#else
  return "";
#endif
}

void o2scl::rng_set_seed(rng<> &r, int mpi_size, int mpi_rank,
                         int verbose) {

  int i_thread=0;
  
#ifdef O2SCL_SET_OPENMP
  i_thread=omp_get_thread_num();
#endif
  
#ifdef O2SCL_SET_OPENMP
#pragma omp critical (o2scl_make_rng_thread_safe)
#endif
  {
    int shift=1;
    if (mpi_size>0) {
      // This works as long as there are not more than 1e6 threads
      shift+=1000001*mpi_rank;
    }
    
    if (o2scl_settings.seed==0) {
      // If the seed has never been used before, set it equal to the time
      o2scl_settings.seed=time(0);
    } else if (o2scl_settings.seed+(unsigned int)shift<
               (unsigned int)mpi_size) {
      // If the next seed flipped around to zero, adjust accordingly
      o2scl_settings.seed=shift;
    } else {
      // Proceed to next seed
      o2scl_settings.seed+=shift;
    }
    if (verbose>0) {
      std::cout << "New RNG for thread " << i_thread << " and rank "
                << mpi_rank << " with seed: "
                << o2scl_settings.seed << std::endl;
    }
    // Set the seed for this rank
    r.set_seed(o2scl_settings.seed);
  }
  return;
}

