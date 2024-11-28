# Interface file for o2scl base classes
#
# Todo:
# 1. Finish o2sclpy/test/test_tensor.py
# 2. Fix interp_vec and interp_krige_optim
# 3. Fix all of the functions at the end
# 4. Check links in python docs to o2scl website
# 5. global functions like abs() and norm() for complex numbers
#
namespace o2scl
py_class_doc |
| Python interface for C++ class ``%name%``.
dll_name o2scl
rst_header |
| .. _base:
|
| Base classes from O2scl
| ==================================
|
| :ref:`O2sclpy <o2sclpy>`
|
| This python interface is not intended to provide the full 
| functionality of the corresponding C++ class.
# 
# Include statements for C++ header file
# 
h_include <o2scl/table_units.h>
h_include <o2scl/table3d.h>
h_include <o2scl/tensor.h>
h_include <o2scl/tensor_grid.h>
h_include <o2scl/find_constants.h>
h_include <o2scl/convert_units.h>
h_include <o2scl/lib_settings.h>
h_include <boost/numeric/ublas/vector.hpp>
h_include <boost/numeric/ublas/matrix.hpp>
h_include <o2scl/format_float.h>
h_include <o2scl/interp_krige.h>
h_include <o2scl/cli.h>
h_include <o2scl/funct.h>
h_include <o2scl/cursesw.h>
h_include <o2scl/funct_to_fp.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/base_python.h>
# 
# Namespace to use in C++ source code
# 
cpp_using std
cpp_using o2scl
#
# Additional python headers
#
py_header from o2sclpy.string import *
py_header 
#
# The interpolation types
#
py_header itp_linear=1
py_header itp_cspline=2
py_header itp_cspline_peri=3
py_header itp_akima=4
py_header itp_akima_peri=5
py_header itp_monotonic=6
py_header itp_steffen=7
py_header itp_nearest_neigh=8
py_header 
#
# Define the force_bytes_string function
#
py_header def force_bytes_string(obj):
py_header     """
py_header     This function returns the bytes object corresponding to ``obj``
py_header     in case it is a string using UTF-8. 
py_header     """
py_header     if (isinstance(obj,numpy.bytes_)==True or
py_header         isinstance(obj,bytes)==True):
py_header         return obj
py_header     if isinstance(obj,o2sclpy.base.std_string):
py_header         return obj.to_bytes()
py_header     return bytes(obj,'utf-8')
#
# ------------------------------------------------------
#
# Class vector<double>
#                              
# Create a python interface to std::vector<double> for vector
# arguments to O2scl functions                             
#
class std::vector<double>
- py_name std_vector
- std_cc                             
- function resize
  - void
  - size_t n                             
- function size
  - size_t
- function operator[]
  - double &
  - size_t n
- function push_back
  - void
  - double x
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: a Python int
|     """
|     return self.size()
|
| def to_numpy(self):
|     """
|     Copy the vector to a numpy array
|
|     Returns: a one-dimensional ``numpy`` array
|     """
|     ret=numpy.zeros((self.size()))
|     for i in range(0,self.size()):
|         ret[i]=self.__getitem__(i)
|     return ret
|
| def append(self,value):
|     """
|     Add an element to the end of the vector
|     """
|     self.push_back(value)
|     return
|
| def from_list(self,lst):
|     """
|     Set the vector with a python list
|     """
|     self.resize(len(lst))
|     for i in range(0,len(lst)):
|         self[i]=lst[i]
|     return
|              
| def erase(self,index):
|     """
|     Erase item at specified index
|     """
|     n=self.size()
|     v2=self.deepcopy()
|     self.resize(n-1)
|     for i in range(0,n-1):
|         if i<index:
|             self[i]=v2[i]
|         else:
|             self[i]=v2[i+1]
|     return
class std::vector<int>
- py_name std_vector_int
- std_cc                             
- function resize
  - void
  - size_t n                             
- function size
  - size_t
- function operator[]
  - int &
  - size_t n
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: a Python int
|     """
|     return self.size()
|
| def to_numpy(self):
|     """
|     Copy the vector to a numpy array
|
|     Returns: a one-dimensional ``numpy`` array with dtype ``int32``
|     """
|     ret=numpy.zeros((self.size()),dtype=numpy.int32)
|     for i in range(0,self.size()):
|         ret[i]=self.__getitem__(i)
|     return ret
class std::vector<size_t>
- py_name std_vector_size_t
- std_cc                             
- function resize
  - void
  - size_t n                             
- function size
  - size_t
- function operator[]
  - size_t &
  - size_t n
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: a Python int
|     """
|     return self.size()
|
| def to_numpy(self):
|     """
|     Copy the vector to a numpy array
|
|     Returns: a one-dimensional ``numpy`` array with dtype ``uint64``
|     """
|     ret=numpy.zeros((self.size()),dtype=numpy.uint64)
|     for i in range(0,self.size()):
|         ret[i]=self.__getitem__(i)
|     return ret
| 
| def init_py(self,v):
|     """
|     Initialize the vector from a python array
|     """
|     self.resize(len(v))
|     for i in range(0,len(v)):
|         self.__setitem__(i,v[i])
|     return
| def __str__(self):
|     """
|     Desc
|     """
|     s='('
|     for i in range(0,len(self)):
|         if i!=len(self)-1:
|             s=s+str(self[i])+','
|         else:
|             s=s+str(self[i])
|     s=s+')'
|     return s
|     
class std::vector<std::string>
- py_name std_vector_string
- std_cc                             
- function resize
  - void
  - size_t n                             
- function size
  - size_t
- function push_back
  - void
  - std::string x
- function operator[]
  - std::string &
  - size_t n
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: a Python int
|     """
|     return self.size()
|
| def append(self,value):
|     """
|     Add an element to the end of the vector
|     """
|     self.push_back(value)
|     return
|
| def set_list(self,ls):
|     """
|     Set the object from a python list
|     """
|     self.resize(len(ls))
|     for i in range(0,len(ls)):
|         self[i]=force_bytes(ls[i])
|     return
|
| def to_list(self):
|     """
|     Set the object from a python list
|     """
|     ret=[]
|     for i in range(0,self.size()):
|         ret.append(self[i])
|     return ret
#
# ------------------------------------------------------
#
# Class ublas_vector
# 
class boost::numeric::ublas::vector<double>
- py_name ublas_vector
- std_cc                             
- function size
  - size_t
- function resize
  - void
  - size_t n
- function operator[]
  - double &
  - size_t i
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: a Python int
|     """
|     return self.size()
|
| def to_numpy(self):
|     """
|     Copy the vector to a numpy array
|
|     Returns: a one-dimensional ``numpy`` array
|     """
|     ret=numpy.zeros((self.size()))
|     for i in range(0,self.size()):
|         ret[i]=self.__getitem__(i)
|     return ret
#
# ------------------------------------------------------
#
# Class ublas_vector_int
# 
class boost::numeric::ublas::vector<int>
- py_name ublas_vector_int
- std_cc                             
- function size
  - size_t
- function resize
  - void
  - size_t n
- function operator[]
  - int &
  - size_t i
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: a Python int
|     """
|     return self.size()
|
| def to_numpy(self):
|     """
|     Copy the vector to a numpy array
|
|     Returns: a one-dimensional ``numpy`` array
|     """
|     ret=numpy.zeros((self.size()),dtype=numpy.intc)
|     for i in range(0,self.size()):
|         ret[i]=self.__getitem__(i)
|     return ret
#
# ------------------------------------------------------
#
# Class ublas_matrix
# 
class boost::numeric::ublas::matrix<double>
- py_name ublas_matrix
- std_cc                             
- function size1
  - size_t
- function size2
  - size_t
- function resize
  - void
  - size_t m
  - size_t n
- function operator()
  - double &
  - size_t m
  - size_t n    
- extra_py |
| def to_numpy(self):
|     """
|     Copy the vector to a numpy matrix
|
|     Returns: a two-dimensional ``numpy`` array, with dimension
|     ``size1(),size2()``.
|     """
|     ret=numpy.zeros((self.size1(),self.size2()))
|     for i in range(0,self.size1()):
|         for j in range(0,self.size2()):
|             ret[i,j]=self.__getitem__((i,j))
|     return ret
#
# ------------------------------------------------------
#
# Class ublas_matrix_int
# 
class boost::numeric::ublas::matrix<int>
- py_name ublas_matrix_int
- std_cc                             
- function size1
  - size_t
- function size2
  - size_t
- function resize
  - void
  - size_t m
  - size_t n
- function operator()
  - int &
  - size_t m
  - size_t n    
- extra_py |
| def to_numpy(self):
|     """
|     Copy the ublas matrix to a numpy matrix
|
|     Returns: a two-dimensional ``numpy`` array, with dimension
|     ``size1(),size2()``.
|     """
|     ret=numpy.zeros(self.size1(),self.size2(),dtype=numpy.intc)
|     for i in range(0,self.size1()):
|         for j in range(0,self.size2()):
|             ret[i,j]=self.__getitem__((i,j))
|     return ret
#
# ------------------------------------------------------
#
# Class vector<vector<double>>
#                              
# Create a python interface to std::vector<std::vector<double>> for
# vector of vector arguments to O2scl functions                             
#
class std::vector<std::vector<double>>
- py_name std_vector_vector
- std_cc                             
- function resize
  - void
  - size_t n                             
- function size
  - size_t
- function operator[]
  - std::vector<double> &
  - size_t n
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: a Python int
|     """
|     return self.size()
#
# ------------------------------------------------------
#
# Class vector<vector<string>>
#                              
# Create a python interface to std::vector<std::vector<std::string>> for
# vector of vector arguments to O2scl functions                             
#
class std::vector<std::vector<std::string>>
- py_name vec_vec_string
- std_cc                             
- function resize
  - void
  - size_t n                             
- function size
  - size_t
- function operator[]
  - std::vector<std::string> &
  - size_t n
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: a Python int
|     """
|     return self.size()
#
# ------------------------------------------------------
#
# Class complex<double>
#                              
# Create a python interface to std::vector<double> for vector
# arguments to O₂scl functions                             
#
class std::complex<double>
- py_class_doc |
| Note that python complex numbers are immutable, but this class is
| not, so the real and imaginary parts can be changed with real_set()
| and imag_set(). 
|     
- py_name std_complex
- cons init
  - double re
  - double im
- function real
  - double
  - py_name real
- function real
  - void
  - py_name real_set
  - double value
- function imag
  - double
  - py_name imag
- function imag
  - void
  - py_name imag_set
  - double value
- extra_py |
| def to_python(self):
|     """
|     Convert to a python complex number
|
|     Returns: a python complex number
|     """
|     ret=self.real()+self.imag()*1j
|     return ret
# -------------------------------------------------------------------
#
# Set the python class documentation for the following classes
# 
py_class_doc |
| Python interface for O₂scl class ``%name%``,
| see
| https://awsteiner.org/code/o2scl/html/class/%name%.html .
|
#
# ------------------------------------------------------
# 
# Class lib_settings_class
#
class lib_settings_class
- function get_data_dir
  - std::string
- function set_data_dir
  - int
  - std::string dir
- function get_doc_dir
  - std::string
- function set_doc_dir
  - int
  - std::string dir
- function openmp_support
  - bool
- function readline_support
  - bool
- function ncurses_support
  - bool
- function armadillo_support
  - bool
- function eigen_support
  - bool
- function fftw_support
  - bool
- function hdf5_compression_support
  - bool
- function system_type
  - std::string
- function range_check
  - bool
- function time_compiled
  - std::string
- function date_compiled
  - std::string
- function o2scl_version
  - std::string
- function config_h_report
  - void
- function get_convert_units
  - convert_units<> &
#
# ------------------------------------------------------
# 
# Class table
#
class table<>
- py_class_doc |
| Python interface for O2scl class ``table``,
| see
| https://awsteiner.org/code/o2scl/html/class/table.html .
- py_name table
- std_cc                             
- function operator[]
  - const std::vector<double> &
  - std::string col
- function set
  - void
  - std::string col
  - size_t row
  - double val
- function get
  - double
  - std::string col
  - size_t row
- function get_ncolumns
  - size_t
- function get_nlines
  - size_t
- function set_nlines
  - void
  - size_t lines
- function get_maxlines
  - size_t
- function set_maxlines
  - void
  - size_t llines
- function set_nlines_auto
  - void
  - size_t il
- function inc_maxlines
  - void
  - size_t llines
- function new_column
  - void
  - std::string col
- function get_column_name
  - std::string
  - size_t icol
- function rename_column
  - void
  - std::string src
  - std::string dest
- function delete_column
  - void
  - std::string col
- function get_sorted_name
  - std::string
  - size_t icol
- function init_column
  - void
  - std::string scol
  - double val                             
- function is_column
  - bool
  - std::string scol
- function lookup_column
  - size_t
  - std::string scol
- function copy_column
  - void
  - std::string src
  - std::string dest
- function add_col_from_table
  - void
  - io table<> &source
  - std::string src_index
  - std::string src_col
  - std::string dest_index
  - std::string dest_col
- function insert_table
  - void
  - io table<> &source
  - std::string src_index
  - bool allow_extrap
  - std::string dest_index
- function add_table
  - void
  - io table<> &source
- function new_row
  - void
  - size_t n
- function copy_row
  - void
  - size_t src
  - size_t dest
- function delete_row
  - void
  - std::string scol
  - double val
- function delete_rows_func
  - void
  - std::string func
- function delete_rows_ends
  - void
  - size_t row_start
  - size_t row_end
- function line_of_names
  - void
  - std::string names
#
# This version requires a std_vector as an argument. An alternate
# version below takes any python array but requires a copy.
# 
- function line_of_data
  - void
  - py_name line_of_data_vector
  - io std::vector<double> &data
- function insert_row
  - void
  - size_t nv
  - io std::vector<double> &data
  - size_t row
- function ordered_lookup
  - size_t
  - std::string scol
  - double val
- function lookup
  - size_t
  - std::string scol
  - double val                             
- function lookup_val
  - size_t
  - std::string scol
  - double val                             
  - std::string scol2
- function set_interp_type
  - void
  - size_t interp_type
- function get_interp_type
  - size_t
- function interp
  - double
  - py_name interp
  - std::string sx
  - double x0
  - std::string sy
- function interp
  - double
  - py_name interp_index
  - size_t ix
  - double x0
  - size_t iy
- function deriv
  - void
  - py_name deriv_col    
  - std::string x    
  - std::string y
  - std::string yp    
- function deriv
  - double
  - py_name deriv
  - std::string sx
  - double x0
  - std::string sy   
- function deriv
  - double
  - py_name deriv_index
  - size_t ix
  - double x0
  - size_t iy
- function deriv2
  - void
  - py_name deriv2_col    
  - std::string x    
  - std::string y
  - std::string yp    
- function deriv2
  - double
  - py_name deriv2
  - std::string sx
  - double x0
  - std::string sy   
- function deriv2
  - double
  - py_name deriv2_index
  - size_t ix
  - double x0
  - size_t iy
- function integ
  - double
  - py_name integ
  - std::string sx
  - double x1
  - double x2                             
  - std::string sy   
- function integ
  - double
  - py_name integ_index
  - size_t ix
  - double x1
  - double x2                             
  - size_t iy
- function integ
  - void
  - py_name integ_col    
  - std::string x    
  - std::string y
  - std::string yi
- function max
  - double
  - std::string max    
- function min
  - double
  - std::string min    
- function zero_table
  - void
- function clear
  - void
- function clear_data
  - void
- function clear_table
  - void
- function clear_constants
  - void
- function sort_table
  - void
  - std::string scol    
- function sort_column
  - void
  - std::string scol
- function average_col_roll
  - void
  - std::string col_name
  - size_t window
- function average_rows
  - void
  - size_t window
  - bool rolling
- function is_valid
  - void
- function functions_columns
  - void
  - std::string list    
- function function_column
  - void
  - std::string function
  - std::string scol
- function row_function
  - double
  - std::string scol
  - size_t row
- function function_find_row
  - size_t
  - std::string function
- function summary
  - void                             
- extra_py |
| def line_of_data(self,v):
|     """
|     Copy ``v`` to an :class:`std_vector` object and add the line of
|     data to the table
|     """
|     # Create a std_vector object and copy the data over
|     vec=std_vector()
|     vec.resize(len(v))
|     for i in range(0,len(v)):
|         vec[i]=v[i]
|     self.line_of_data_vector(vec)
|     return
|                              
| def row_to_dict(self,row):
|     """
|     Convert the specified row to a python dictionary
|     """
|     dct={}
|     for i in range(0,self.get_ncolumns()):
|         dct[self.get_column_name(i)]=self.get(self.get_column_name(i),
|                                               row)
|     return dct
#
# ------------------------------------------------------
# 
# Class table_units<>
#
class table_units<>
- py_class_doc |
| Python interface for O2scl class ``table_units``,
| see
| https://awsteiner.org/code/o2scl/html/class/table_units.html .
- parent table<>
- std_cc                             
- py_name table_units
- function set_unit
  - void
  - std::string col
  - std::string unit
- function get_unit
  - std::string
  - std::string col
- function line_of_units
  - void
  - std::string unit_line    
- function remove_unit
  - void
  - std::string col
- function convert_to_unit
  - int
  - std::string col
  - std::string unit
  - bool err_on_fail [True]
# 
# Shared pointer for table_units<> objects
#
shared_ptr table_units<>
- py_name table_units
#
# ------------------------------------------------------
# 
# Class uniform_grid
#
class uniform_grid<>
- py_class_doc |
| Python interface for O2scl class ``uniform_grid``,
| see
| https://awsteiner.org/code/o2scl/html/class/uniform_grid.html .
- py_name uniform_grid
- function get_nbins
  - size_t                             
- function get_npoints
  - size_t                             
- function is_log
  - bool                             
- function get_start
  - double                             
- function get_end
  - double                             
- function get_width
  - double
- function operator[]
  - double
  - size_t n
- function vector
  - void
  - out std::vector<double> &v
- extra_py |
| def to_numpy(self):
|     """
|     Copy the ``uniform_grid`` object to a numpy array
|
|     Returns: a one-dimensional ``numpy`` array
|     """
|     v=std_vector()
|     self.vector(v)
|     ret=numpy.zeros((self.get_npoints()))
|     for i in range(0,self.get_npoints()):
|         ret[i]=v[i]
|     return ret
#
# ------------------------------------------------------
# 
# Class uniform_grid_end
#
class uniform_grid_end<>
- py_class_doc |
| Python interface for O2scl class ``uniform_grid_end``,
| see
| https://awsteiner.org/code/o2scl/html/class/uniform_grid_end.html .
- py_name uniform_grid_end                             
- parent uniform_grid<>
- cons init
  - double start
  - double end
  - size_t n_bins                             
#
# ------------------------------------------------------
# 
# Class uniform_grid_width
#
class uniform_grid_width<>
- py_class_doc |
| Python interface for O2scl class ``uniform_grid_width``,
| see
| https://awsteiner.org/code/o2scl/html/class/uniform_grid_width.html .
- py_name uniform_grid_width
- parent uniform_grid<>
- cons init
  - double start
  - double width
  - size_t n_bins                             
#
# ------------------------------------------------------
# 
# Class uniform_grid_end_width
#
class uniform_grid_end_width<>
- py_class_doc |
| Python interface for O2scl class ``uniform_grid_end_width``,
| see
| https://awsteiner.org/code/o2scl/html/class/uniform_grid_end_width.html .
- py_name uniform_grid_end_width                             
- parent uniform_grid<>
- cons init
  - double start
  - double end
  - double width
#
# ------------------------------------------------------
# 
# Class uniform_grid_log_end
#
class uniform_grid_log_end<>
- py_class_doc |
| Python interface for O2scl class ``uniform_grid_log_end``,
| see
| https://awsteiner.org/code/o2scl/html/class/uniform_grid_log_end.html .
- py_name uniform_grid_log_end
- parent uniform_grid<>
- cons init
  - double start
  - double end
  - size_t n_bins                             
#
# ------------------------------------------------------
# 
# Class uniform_grid_log_width
#
class uniform_grid_log_width<>
- py_class_doc |
| Python interface for O2scl class ``uniform_grid_log_width``,
| see
| https://awsteiner.org/code/o2scl/html/class/uniform_grid_log_width.html .
- py_name uniform_grid_log_width
- parent uniform_grid<>
- cons init
  - double start
  - double width
  - size_t n_bins                             
#
# ------------------------------------------------------
# 
# Class uniform_grid_log_end_width
#
class uniform_grid_log_end_width<>
- py_class_doc |
| Python interface for O2scl class ``uniform_grid_log_end_width``,
| see
| https://awsteiner.org/code/o2scl/html/class/uniform_grid_log_end_width.html .
- py_name uniform_grid_log_end_width
- parent uniform_grid<>
- cons init
  - double start
  - double end
  - double width
#
# ------------------------------------------------------
# 
# Class table3d
#
class table3d
- std_cc
- function set_size
  - void
  - size_t nx
  - size_t ny
- function set_xy
  - void
  - py_name set_xy
  - std::string x_name
  - size_t nx
  - io std_vector &x
  - std::string y_name
  - size_t ny
  - io std_vector &y
- function set_xy
  - void
  - py_name set_xy_grid
  - std::string x_name
  - io uniform_grid<double> &x_grid
  - std::string y_name
  - io uniform_grid<double> &y_grid
- function set
  - void
  - py_name set
  - size_t ix
  - size_t iy
  - std::string name
  - double val
- function get
  - double
  - py_name get
  - size_t ix
  - size_t iy
  - std::string name
- function get
  - double
  - py_name get_i
  - size_t ix
  - size_t iy
  - size_t iz
- function set
  - void
  - py_name set_i
  - size_t ix
  - size_t iy
  - size_t iz
  - double val
- function set_val
  - void
  - double x
  - double y
  - std::string name
  - double val
- function get_val
  - double
  - double x
  - double y
  - std::string name
- function set_grid_x
  - void
  - size_t ix
  - double val    
- function set_grid_y
  - void
  - size_t iy
  - double val    
- function get_grid_x
  - double
  - size_t ix
- function get_grid_y
  - double
  - size_t iy
- function get_x_name
  - std::string
- function get_y_name
  - std::string
- function set_x_name
  - void
  - std::string name
- function set_y_name
  - void
  - std::string name
- function get_size
  - void
  - out size_t &nx
  - out size_t &ny
- function get_nx
  - size_t
- function get_ny
  - size_t
- function get_nslices
  - size_t
- function is_size_set
  - bool
- function is_xy_set
  - bool
- function get_slice_name
  - std::string
  - size_t i
- function new_slice
  - void
  - std::string slice
- function set_slice_all
  - void
  - std::string name
  - double val
- function lookup_slice
  - size_t
  - std::string name
- function is_slice
  - bool
  - std::string name
  - out size_t &ix
- function rename_slice
  - void
  - std::string name1
  - std::string name2
- function copy_slice
  - void
  - std::string name1
  - std::string name2
- function get_slice
  - boost::numeric::ublas::matrix<double> &
  - py_name get_slice
  - std::string slice
- function get_slice
  - boost::numeric::ublas::matrix<double> &
  - py_name get_slice_i    
  - std::string slice
- function lookup_x
  - void
  - double val
  - size_t ix
- function lookup_y
  - void
  - double val
  - size_t iy
- function interp
  - double
  - double x
  - double y
  - std::string name
- function deriv_x
  - double
  - double x
  - double y
  - std::string name
- function deriv_y
  - double
  - double x
  - double y
  - std::string name
- function deriv_xy
  - double
  - double x
  - double y
  - std::string name
- function integ_x
  - double
  - double x1
  - double x2    
  - double y
  - std::string name
- function integ_y
  - double
  - double x
  - double y1
  - double y2    
  - std::string name
- function zero_table
  - void
- function clear
  - void
- function function_matrix 
  - int
  - std::string function
  - out boost::numeric::ublas::matrix<double> &mat
  - bool throw_on_err
- function function_slice
  - void
  - std::string function
  - std::string slice
- function summary
  - void
#
# ------------------------------------------------------
# 
# Class index_spec
#
class index_spec
- size_t type
- size_t ix1
- size_t ix2
- size_t ix3
- double val1
- double val2
- double val3
#
# ------------------------------------------------------
# 
# Class ix_index
#
class ix_index
- no_def_cons
- cons init
  - size_t ix
#
# ------------------------------------------------------
# 
# Class ix_fixed
#
class ix_fixed
- no_def_cons
- cons init
  - size_t ix
  - size_t ix2
#
# ------------------------------------------------------
# 
# Class ix_sum
#
class ix_sum
- no_def_cons
- cons init
  - size_t ix
#
# ------------------------------------------------------
# 
# Class ix_trace
#
class ix_trace
- no_def_cons
- cons init
  - size_t ix
  - size_t ix2
#
# ------------------------------------------------------
# 
# Class ix_reverse
#
class ix_reverse
- no_def_cons
- cons init
  - size_t ix
#
# ------------------------------------------------------
# 
# Class ix_range
#
class ix_range
- no_def_cons
- cons init
  - size_t ix
  - size_t start
  - size_t end
#
# ------------------------------------------------------
# 
# Class ix_interp
#
class ix_interp
- no_def_cons
- cons init
  - size_t ix
  - double v
#
# ------------------------------------------------------
# 
# Class ix_grid
#
class ix_grid
- no_def_cons
- cons init
  - size_t ix
  - double start
  - double end
  - size_t n_bins
  - bool log [false]
#
# ------------------------------------------------------
# 
# Class ix_gridw
#
class ix_gridw
- no_def_cons
- cons init
  - size_t ix
  - double start
  - double end
  - double width
  - bool log [false]
#
# ------------------------------------------------------
# 
# Class tensor
#
class tensor<>
- py_class_doc |
| Python interface for O2scl class ``tensor``,
| see
| https://awsteiner.org/code/o2scl/html/class/tensor.html .
- std_cc                             
- py_name tensor
- cons create_size
  - py_name create_size_vector
  - size_t rank
  - io std::vector<size_t> &sizes
- function is_valid
  - void
- function clear
  - void
- function set
  - void
  - py_name set_vector
  - vector<size_t> &index
  - double val
- function set_all
  - void
  - double x
- function swap_data
  - void
  - std::vector<double> &data
- function get
  - double
  - py_name get_vector
  - io vector<size_t> &index
- function resize
  - void
  - py_name resize_vector
  - size_t n
  - io vector<size_t> &index
- function get_rank
  - size_t
- function get_size
  - size_t
  - size_t i
- function get_size_arr
  - const std::vector<size_t> &
- function get_data
  - const vector<double> &
- function total_size
  - size_t
- function pack_indices
  - size_t
  - io std::vector<size_t> &index
- function unpack_index
  - void
  - size_t ix
  - out std::vector<size_t> &index
- function min_value
  - double
- function min_index
  - void
  - out std::vector<size_t> &index
- function min
  - void
  - out std::vector<size_t> &ix
  - out double &value    
- function max_value
  - double
- function max_index
  - void
  - out std::vector<size_t> &index
- function max
  - void
  - out std::vector<size_t> &ix
  - out double &value    
- function minmax_value
  - void
  - out double &min
  - out double &max
- function minmax_index
  - void
  - out std::vector<size_t> &min
  - out std::vector<size_t> &max
- function minmax
  - void
  - out std::vector<size_t> &min_ix
  - out double &min_value    
  - out std::vector<size_t> &max_ix
  - out double &max_value    
- function total_sum
  - double
- function copy_table3d_sum
  - void
  - size_t ix_x
  - size_t ix_y
  - io table3d &tab
  - std::string x_name ["x"]
  - std::string y_name ["y"]
  - std::string slice_name ["z"]
- function copy_table3d
  - void
  - size_t ix_x
  - size_t ix_y
  - io table3d &tab
  - std::string x_name ["x"]
  - std::string y_name ["y"]
  - std::string slice_name ["z"]
- extra_py |
| def create_size(self,v):
|     """
|     Copy ``v`` to an :class:`std_vector_size_t` object and add the line of
|     data to the table
|     """
|     # Create a std_vector object and copy the data over
|     vec=std_vector_size_t()
|     vec.resize(len(v))
|     for i in range(0,len(v)):
|         vec[i]=v[i]
|     self.create_size_vector(vec)
|     return
| 
| def set(self,index,val):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object and add the 
|     data to the table
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     self.set_vector(svst,val)
|     return
| 
| def get(self,index):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object and get the 
|     data from the table
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     return self.get_vector(svst)
|
| def resize(self,index):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object 
|     and resize
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     self.resize_vector(len(svst),svst)
|     return
#
# Function rearrange_and_copy() from tensor.h
#
function rearrange_and_copy<tensor<>,double>
- tensor<>
- py_name rearrange_and_copy
- tensor<> &t
- std::string spec
- int verbose [0]
- bool err_on_fail [true]
function rearrange_and_copy<tensor<int>,int>
- tensor<int>
- py_name rearrange_and_copy_int
- tensor<int> &t
- std::string spec
- int verbose [0]
- bool err_on_fail [true]
function rearrange_and_copy<tensor<size_t>,size_t>
- tensor<size_t>
- py_name rearrange_and_copy_size_t
- tensor<size_t> &t
- std::string spec
- int verbose [0]
- bool err_on_fail [true]
#
# ------------------------------------------------------
# 
# Class tensor_int
#
class tensor<int,std::vector<int>>
- py_class_doc |
| Python interface for O2scl class ``tensor``,
| see
| https://awsteiner.org/code/o2scl/html/class/tensor.html .
- std_cc                             
- py_name tensor_int
- cons create_size
  - py_name create_size_vector
  - size_t rank
  - io std::vector<size_t> &sizes
- function is_valid
  - void
- function clear
  - void
- function set
  - void
  - py_name set_vector
  - vector<size_t> &index
  - int val
- function set_all
  - void
  - int x
- function get
  - int
  - py_name get_vector
  - io vector<size_t> &index
- function resize
  - void
  - py_name resize_vector
  - size_t n
  - io vector<size_t> &index
- function get_rank
  - size_t
- function get_size
  - size_t
  - size_t i
- function get_data
  - const std::vector<int> &
- function total_size
  - size_t
- function pack_indices
  - size_t
  - io const std::vector<size_t> &index
- function unpack_index
  - void
  - size_t ix
  - io const std::vector<size_t> &index
- function min_value
  - int
- function min_index
  - void
  - out std::vector<size_t> &index
- function min
  - void
  - out std::vector<size_t> &index
  - out int &val
- function max_value
  - int
- function max_index
  - void
  - out std::vector<size_t> &index
- function max
  - void
  - out std::vector<size_t> &index
  - out int &val
- function minmax_value
  - void
  - out int &min
  - out int &max
- function minmax_index
  - void
  - out std::vector<size_t> &index_min
  - out std::vector<size_t> &index_max
- function minmax
  - void
  - out std::vector<size_t> &index_min
  - out int &min
  - out std::vector<size_t> &index_max
  - out int &max
- function total_sum
  - int
- function copy_table3d_sum
  - void
  - size_t ix_x
  - size_t ix_y
  - io table3d &tab
  - std::string x_name ["x"]
  - std::string y_name ["y"]
  - std::string slice_name ["z"]
- function copy_table3d
  - void
  - size_t ix_x
  - size_t ix_y
  - io table3d &tab
  - std::string x_name ["x"]
  - std::string y_name ["y"]
  - std::string slice_name ["z"]
- extra_py |
| def create_size(self,v):
|     """
|     Copy ``v`` to an :class:`std_vector_size_t` object and add the line of
|     data to the table
|     """
|     # Create a std_vector object and copy the data over
|     vec=std_vector_size_t()
|     vec.resize(len(v))
|     for i in range(0,len(v)):
|         vec[i]=v[i]
|     self.create_size_vector(vec)
|     return
| 
| def set(self,index,val):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object and add the 
|     data to the table
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     self.set_vector(svst,val)
|     return
| 
| def get(self,index):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object and get the 
|     data from the table
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     return self.get_vector(svst)
|
| def resize(self,index):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object 
|     and resize
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     self.resize_vector(svst)
|     return
#
# ------------------------------------------------------
# 
# Class tensor_size_t
#
class tensor<size_t,std::vector<size_t>>
- py_class_doc |
| Python interface for O2scl class ``tensor``,
| see
| https://awsteiner.org/code/o2scl/html/class/tensor.html .
- std_cc                             
- py_name tensor_size_t
- cons create_size
  - py_name create_size_vector
  - size_t rank
  - io std::vector<size_t> &sizes
- function is_valid
  - void
- function clear
  - void
- function set
  - void
  - py_name set_vector
  - vector<size_t> &index
  - size_t val
- function set_all
  - void
  - size_t x
- function get
  - int
  - py_name get_vector
  - io vector<size_t> &index
- function resize
  - void
  - py_name resize_vector
  - size_t n
  - io vector<size_t> &index
- function get_rank
  - size_t
- function get_size
  - size_t
  - size_t i
- function get_data
  - const std::vector<size_t> &
- function total_size
  - size_t
- function min_value
  - size_t
- function min_index
  - void
  - out std::vector<size_t> &index
- function min
  - void
  - out std::vector<size_t> &index
  - out size_t &val
- function max_value
  - size_t
- function max_index
  - void
  - out std::vector<size_t> &index
- function max
  - void
  - out std::vector<size_t> &index
  - out size_t &val
- function minmax_value
  - void
  - out size_t &min
  - out size_t &max
- function minmax_index
  - void
  - out std::vector<size_t> &index_min
  - out std::vector<size_t> &index_max
- function minmax
  - void
  - out std::vector<size_t> &index_min
  - out size_t &min
  - out std::vector<size_t> &index_max
  - out size_t &max
- function total_sum
  - size_t
- function copy_table3d_sum
  - void
  - size_t ix_x
  - size_t ix_y
  - io table3d &tab
  - std::string x_name ["x"]
  - std::string y_name ["y"]
  - std::string slice_name ["z"]
- function copy_table3d
  - void
  - size_t ix_x
  - size_t ix_y
  - io table3d &tab
  - std::string x_name ["x"]
  - std::string y_name ["y"]
  - std::string slice_name ["z"]
- extra_py |
| def create_size(self,v):
|     """
|     Copy ``v`` to an :class:`std_vector_size_t` object and add the line of
|     data to the table
|     """
|     # Create a std_vector object and copy the data over
|     vec=std_vector_size_t()
|     vec.resize(len(v))
|     for i in range(0,len(v)):
|         vec[i]=v[i]
|     self.create_size_vector(vec)
|     return
|   
| def set(self,index,val):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object and add the 
|     data to the table
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     self.set_vector(svst,val)
|     return
| 
| def get(self,index):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object and get the 
|     data from the table
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     return self.get_vector(svst)
|
| def resize(self,index):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object 
|     and resize
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     self.resize_vector(svst)
|     return
#
# ------------------------------------------------------
#
# Class tensor_grid
#
class tensor_grid<>
- py_class_doc |
| Python interface for O2scl class ``tensor_grid``,
| see
| https://awsteiner.org/code/o2scl/html/class/tensor_grid.html .
- parent tensor<>
- std_cc                             
- py_name tensor_grid
- function is_valid
  - void
- function set_val
  - void
  - py_name set_val_vector
  - io vector<double> &grid_point
  - double val
- function get_val
  - double
  - py_name get_val_vector
  - io vector<double> &grid_point
- function resize
  - void
  - py_name resize_vector
  - size_t rank
  - io vector<size_t> &dim  
- function is_grid_set
  - bool
- function set_grid_packed
  - void
  - io vector<double> &grid
- function set_grid
  - void
  - py_name set_grid_vec_vec
  - vector<vector<double>> &grid_vecs
- function default_grid
  - void
- function set_grid_i_vec    
  - void
  - size_t i
  - io vector<double> &grid
- function set_grid_i_func
  - void
  - size_t ix
  - std::string func
- function get_grid
  - double
  - py_name get_grid
  - size_t i
  - size_t j
- function get_grid
  - const vector<double> &
  - py_name get_grid_packed
- function set_grid
  - void
  - py_name set_grid
  - size_t i
  - size_t j
  - double val
#- function lookup_grid_val
#  - void
#  - io const double &val
#  - out double &val2
- function lookup_grid
  - size_t
  - size_t i
  - double val
- function copy_slice_interp
  - tensor_grid<>
  - io std::vector<size_t> &ifix
  - oi std::vector<double> &vals
#- function convert_table3d_sum
#  - void
#  - size_t ix_x
#  - size_t ix_y
#  - out table3d &tab
#  - std::string x_name ["x"]
#  - std::string y_name ["y"]
#  - std::string slice_name ["z"]
- function copy_table3d_align
  - void
  - size_t ix_x
  - size_t ix_y
  - io vector<size_t> &index
  - io table3d &tab
  - std::string z_name ["z"]
- function copy_table3d_align_setxy
  - void
  - size_t ix_x
  - size_t ix_y
  - io vector<size_t> &index
  - io table3d &tab
  - std::string x_name ["x"]
  - std::string y_name ["y"]
  - std::string z_name ["z"]
- function copy_table3d_interp
  - void
  - size_t ix_x
  - size_t ix_y
  - io std::vector<size_t> &index
  - out table3d &tab
  - std::string slice_name ["z"]
- function copy_table3d_interp_values
  - void
  - size_t ix_x
  - size_t ix_y
  - io std::vector<double> &values
  - out table3d &tab
  - std::string slice_name ["z"]
  - int verbose [0]
- function copy_table3d_interp_values_setxy
  - void
  - size_t ix_x
  - size_t ix_y
  - io std::vector<double> &values
  - out table3d &tab
  - std::string x_name ["x"]
  - std::string y_name ["y"]
  - std::string slice_name ["z"]
- function clear
  - void
- function set_interp_type
  - void
  - size_t interp_type
- function interp_linear_partial
  - double
  - io const std::vector<size_t> &ix_to_interp
  - io std::vector<size_t> &ix
  - io const std::vector<double> &val
- function interp_linear
  - double
  - io vector<double> &v
- function from_table3d_fermi
  - void
  - const table3d &t3d
  - std::string slice
  - size_t n_points
  - double low [0.0]
  - double high [0.0]
  - double width [0.0]
- extra_py |
| def resize(self,index):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object 
|     and resize
|     """
|     svst=std_vector_size_t()
|     svst.init_py(index)
|     self.resize_vector(len(svst),svst)
|     return
#
# Function grid_rearrange_and_copy() from tensor_grid.h
#
function grid_rearrange_and_copy<tensor_grid<>,double>
- tensor_grid<>
- py_name grid_rearrange_and_copy
- tensor_grid<> &t
- std::string spec
- int verbose [0]
- bool err_on_fail [true]
#
# ------------------------------------------------------
# 
# Class find_constants::const_entry
#
class find_constants<>::const_entry
- py_name find_constants_const_entry
- std::vector<std::string> names
- std::string unit
- int unit_flag
- double val
- std::string source
- int m
- int k
- int s
- int K
- int A
- int mol
- int cd
#
# ------------------------------------------------------
# 
# Class find_constants
#
class find_constants<>
- py_name find_constants
- function output_list_cout
  - void
- function add_constant
  - void
  - io find_constants<>::const_entry &f
  - int verbose [0]    
- function del_constant
  - void
  - io std::string &name
  - int verbose [0]
#
# ------------------------------------------------------
#
# Class convert_units_der_unit  
#
class convert_units<>::der_unit
- py_name convert_units_der_unit
- extra_py |
| def set(self,label,val,name='',m=0,k=0,s=0,K=0,A=0,mol=0,cd=0):
|     """
|     Set the properties of a derived unit
|     FIXME: better docs here
|     """
|     label2=std_string()
|     label2.init_bytes(force_bytes(label))
|     self.set_label(label2)
|     self.val=val
|     name2=std_string()
|     name2.init_bytes(force_bytes(name))
|     self.set_name(name2)
|     self.m=m
|     self.k=k
|     self.s=s
|     self.K=K
|     self.A=A
|     self.mol=mol
|     self.cd=cd
|     return
- std::string label
- int m
- int k
- int s
- int K
- int A
- int mol
- int cd
- double val
- std::string name  
#
# ------------------------------------------------------
#  
# Class convert_units<>
#
# Note that 'from' is a reserved keyword in python so we
# rename it to 'frm' instead
#
class convert_units<>
- py_class_doc |
| Python interface for O2scl class ``convert_units``,
| see
| https://awsteiner.org/code/o2scl/html/class/convert_units.html .
- py_name convert_units
- extra_py |
| def add_unit(self,label,val,name='',m=0,k=0,s=0,K=0,A=0,mol=0,cd=0):
|     """
|     Add a unit
|     """
|     du=convert_units_der_unit()
|     du.set(label,val,name,m,k,s,K,A,mol,cd)
|     self.add_unit_internal(du)
|     return
- function convert
  - double
  - std::string frm
  - std::string to
  - double val
- function convert_ret
  - int
  - std::string frm
  - std::string to
  - double val
  - double converted
- function del_unit
  - void
  - io std::string &name
- function add_unit
  - void
  - py_name add_unit_internal
  - const io convert_units<>::der_unit &d
- function set_natural_units
  - void
  - bool c_is_one [true]
  - bool hbar_is_one [true]
  - bool kb_is_one [true]
- function is_in_cache
  - int
  - std::string frm
  - std::string to     
- function remove_cache
  - int
  - std::string frm
  - std::string to
- function clear_cache
  - void
- function test_unique
  - void
- function print_cache
  - void
- function print_units_cout
  - void
- function find_print
  - void
  - std::string name
  - std::string unit
  - size_t prec
  - bool use_regex
- function find_unique
  - double
  - std::string name
  - std::string unit
  - bool use_regex [false]
- int verbose
- bool err_on_fail
- bool combine_two_conv
#
# ------------------------------------------------------
#
# Class columnify
#
class columnify
- static const int align_left
- static const int align_right
- static const int align_lmid
- static const int align_rmid
- static const int align_dp
- static const int align_lnum
# - function add_spaces<std::vector<std::vector<std::string>>,std::vector<int>>
#   - int
#   - 
#     template<class mat_string_t, class vec_int_t>
#     int add_spaces(const mat_string_t &table_in, size_t ncols, size_t nrows, 
# 		   vec_int_t &align_spec, mat_string_t &table_out) {
#
# ------------------------------------------------------
#
# Class format_float
#
class format_float
- function set_sig_figs
  - void
  - size_t sig_figs
- function set_exp_limits
  - void
  - int min
  - int max
- function set_pad_zeros
  - void
  - bool pad
- function set_dec_point
  - void
  - std::string dec_point
- function set_exp_digits
  - void
  - int d
- function html_mode
  - void
- function latex_mode
  - void
- function c_mode
  - void
- function convert
  - std::string
  - double x
  - bool debug [false]
#
# ------------------------------------------------------
#
# Class interp_vec
#
class interp_vec<std::vector<double>>
- py_name interp_vec  
#- cons interp_vec_pre
#  - size_t n
#  - std::vector<double> &x
#  - std::vector<double> &y
#  - int interp_type
- function set
  - void
  - size_t n
  - io std::vector<double> &x
  - io std::vector<double> &y
  - int interp_type
- function clear
  - void
- function eval
  - double
  - double x0
- function deriv
  - double
  - double x0
- function deriv2
  - double
  - double x0
- function integ
  - double
  - double x1
  - double x2
#
# ------------------------------------------------------
#
# Class interp_krige_optim
#
class interp_krige_optim<std::vector<double>,std::vector<double>,covar_funct_rbf_noise>
- py_name interp_krige_optim_rbf_noise
- static const size_t mode_loo_cv
- static const size_t mode_loo_cv_bf
- static const size_t mode_max_lml
- int verbose
- size_t mode
#- bool full_min
- function set
  - void
  - size_t size
  - io const std::vector<double> &x
  - io const std::vector<double> &y
- function eval
  - double
  - double x0
- function deriv
  - double
  - double x0
- function deriv2
  - double
  - double x0
- function sigma
  - double
  - double x0
- function sample
  - double
  - double x0
- function sample_vec
  - void
  - io std::vector<double> &x
  - io std::vector<double> &y
#
# ------------------------------------------------------
#
# Functions and classes from misc.h
#
class terminal
- function is_redirected
  - bool
- function str_len
  - size_t
  - std::string str
- function hrule
  - std::string
  - size_t n [78]
- function cyan_fg
  - std::string
- function magenta_fg
  - std::string
- function yellow_fg
  - std::string
- function red_fg
  - std::string
- function green_fg
  - std::string
- function blue_fg
  - std::string
- function cyan_bg
  - std::string
- function magenta_bg
  - std::string
- function yellow_bg
  - std::string
- function red_bg
  - std::string
- function green_bg
  - std::string
- function blue_bg
  - std::string
- function default_fgbg
  - std::string
- function bold
  - std::string
- function eight_bit_fg
  - std::string
  - short col
- function eight_bit_bg
  - std::string
  - short col
#- function three_byte_fg
#  - std::string
#  - short red
#  - short green
#  - short blue
#- function three_byte_bg
#  - std::string
#  - short red
#  - short green
#  - short blue
- function lowint
  - std::string
- function underline
  - std::string
- function reverse
  - std::string
- function alt_font
  - std::string
- function normal_font
  - std::string
- function eight_bit_summ
  - std::string
- function three_byte_summ
  - std::string
- function three_byte_summ_long
  - std::string
#
# ------------------------------------------------------
#  
function fermi_function
- double
- double x
#
# ------------------------------------------------------
#  
function bose_function
- double
- double x
#
# ------------------------------------------------------
#  
function quadratic_extremum_x<double>
- double
- py_name quadratic_extremum_x
- double x1
- double x2
- double x3
- double y1
- double y2
- double y3
#
# ------------------------------------------------------
#  
function quadratic_extremum_y<double>
- double
- py_name quadratic_extremum_y
- double x1
- double x2
- double x3
- double y1
- double y2
- double y3
#
# ------------------------------------------------------
#  
function screenify<vector<std::string>>
- void  
- py_name screenify
- size_t nin
- io vector<std::string> &in_cols
- out vector<std::string> &out_cols
- size_t max_size [80]
#
# ------------------------------------------------------
#  
function file_exists
- bool
- std::string fname
#
# ------------------------------------------------------
#  
function RGBtoHSV
- void
- double r
- double g
- double b
- out double &h
- out double &s
- out double &v
#
# ------------------------------------------------------
#  
function HSVtoRGB
- void
- double h
- double s
- double v
- out double &r
- out double &g
- out double &b
#
# ------------------------------------------------------
#  
function wordexp_single_file
- void
- io std::string &fname  
#
# ------------------------------------------------------
#  
function wordexp_wrapper
- void
- std::string word
- out std::vector<std::string> &matches
#
# ------------------------------------------------------
#  
class gen_test_number<double>
- py_name gen_test_number
- function reset
  - void
- function set_radix
  - void
  - double r
- function gen
  - double
#
# ------------------------------------------------------
#
# Functions from lib_settings.h
#
function function_to_double
- double
- std::string s
- int verbose [0]
#
# ------------------------------------------------------
#  
function function_to_double_nothrow
- int
- std::string s
- out double &result
- int verbose [0]
#
# ------------------------------------------------------
#  
function find_constant
- double
- std::string name
- std::string unit  
#
# ------------------------------------------------------
#
# Functions from string_conv.h
#  
function string_to_uint_list<vector<size_t>>
- int
- py_name string_to_uint_list
- io const std::string &x
- out vector<size_t> &list
#  
#function string_to_char_array
#- void
#- std::string s
#- char *x
#- int len
#
# ------------------------------------------------------
#
function rewrap_keep_endlines
- void
- std::string str
- out std::vector<std::string> &sv
- size_t ncol [79]
- int verbose [0]
- bool ignore_vt100 [True]  
#
# ------------------------------------------------------
#
# Functions from interp.h
# 
function vector_level_count<std::vector<double>,std::vector<double>>
- size_t
- py_name vector_level_count  
- double level
- size_t n
- std::vector<double> &x
- std::vector<double> &y
#
# ------------------------------------------------------
#
function vector_deriv_interp<std::vector<double>,std::vector<double>>
- void
- py_name vector_deriv_interp
- size_t n  
- io std::vector<double> &v
- out std::vector<double> &dv
- size_t interp_type [2]
#
# ------------------------------------------------------
#
function vector_deriv2_interp<std::vector<double>,std::vector<double>>
- void
- py_name vector_deriv2_interp
- size_t n  
- io std::vector<double> &v
- out std::vector<double> &dv
- size_t interp_type [2]
#
# ------------------------------------------------------
#
function vector_deriv_xy_interp<std::vector<double>,std::vector<double>,std::vector<double>>
- void
- py_name vector_deriv_xy_interp  
- size_t n  
- io std::vector<double> &vx
- io std::vector<double> &vy
- out std::vector<double> &dv
- size_t interp_type [2]
#
# ------------------------------------------------------
#
function vector_deriv2_xy_interp<std::vector<double>,std::vector<double>,std::vector<double>>
- void
- py_name vector_deriv2_xy_interp  
- size_t n  
- io std::vector<double> &vx
- io std::vector<double> &vy
- out std::vector<double> &dv
- size_t interp_type [2]
#
# ------------------------------------------------------
#
function vector_integ_interp<std::vector<double>>
- double
- py_name vector_integ_interp  
- size_t n  
- io std::vector<double> &vx
- size_t interp_type [2]
#
# ------------------------------------------------------
#
function vector_integ_xy_interp<std::vector<double>,std::vector<double>>
- double
- py_name vector_integ_xy_interp  
- size_t n  
- io std::vector<double> &vx
- io std::vector<double> &vy
- size_t interp_type [2]
#
# ------------------------------------------------------
#
function vector_integ_ul_interp<std::vector<double>>
- double
- py_name vector_integ_ul_interp  
- size_t n  
- double x2
- io std::vector<double> &v
- size_t interp_type [2]
#
# ------------------------------------------------------
#
function vector_integ_ul_xy_interp<std::vector<double>,std::vector<double>>
- double
- py_name vector_integ_ul_xy_interp  
- size_t n  
- double x2  
- io std::vector<double> &vx
- io std::vector<double> &vy
- size_t interp_type [2]
#
# ------------------------------------------------------
#
function vector_find_level<std::vector<double>,std::vector<double>>
- void
- py_name vector_find_level
- double level
- size_t n
- io std::vector<double> &x
- io std::vector<double> &y
- out std::vector<double> &locs
#
# ------------------------------------------------------
#
function vector_invert_enclosed_sum<std::vector<double>,std::vector<double>>
- void
- py_name vector_invert_enclosed_sum
- double sum
- size_t n
- io std::vector<double> &x
- io std::vector<double> &y
- out double &lev
- int boundaries [0]
- int verbose [0]
- bool err_on_fail [true]
#
# ------------------------------------------------------
#
function vector_region_int<std::vector<double>,std::vector<double>>
- int
- py_name vector_region_int  
- size_t n
- io std::vector<double> &x
- io std::vector<double> &y
- double intl
- out std::vector<double> &locs
- int boundaries [0]
- int verbose [0]
- bool err_on_fail [true]
#
# ------------------------------------------------------
#
function vector_region_fracint<std::vector<double>,std::vector<double>>
- int
- py_name vector_region_fracint  
- size_t n
- io std::vector<double> &x
- io std::vector<double> &y
- double intl
- out std::vector<double> &locs
- int boundaries [0]
- int verbose [0]
- bool err_on_fail [true]
#
# ------------------------------------------------------
#
function vector_bound_fracint<std::vector<double>,std::vector<double>>
- int
- py_name vector_bound_fracint  
- size_t n
- io std::vector<double> &x
- io std::vector<double> &y
- double frac
- out double &low  
- out double &high
- int boundaries [0]
- int verbose [0]
- bool err_on_fail [true]
#
# ------------------------------------------------------
#
function vector_bound_int<std::vector<double>,std::vector<double>>
- int
- py_name vector_bound_int  
- size_t n
- io std::vector<double> &x
- io std::vector<double> &y
- double frac
- out double &low  
- out double &high
- int boundaries [0]
- int verbose [0]
- bool err_on_fail [true]
#
# ------------------------------------------------------
#
function rebin_xy<std::vector<double>,std::vector<double>,std::vector<double>,std::vector<double>>
- void
- py_name rebin_xy
- io std::vector<double> &x
- io std::vector<double> &y
- out std::vector<double> &x_out
- out std::vector<double> &y_out
- size_t n_pts
- size_t interp_type
#
# ------------------------------------------------------
#
function linear_or_log_chi2<std::vector<double>,std::vector<double>>
- double
- py_name linear_or_log_chi2
- io std::vector<double> &x
- io std::vector<double> &y
#
# ------------------------------------------------------
#
function linear_or_log<std::vector<double>,std::vector<double>>
- void
- py_name linear_or_log_pair
- io std::vector<double> &x
- io std::vector<double> &y
- out bool &log_x
- out bool &log_y
#
# ------------------------------------------------------
#
function vector_refine<std::vector<double>,std::vector<double>,double>
- void
- py_name vector_refine
- size_t n
- std::vector<double> &index
- std::vector<double> &data
- size_t factor
- size_t interp_type [2]
#
# ------------------------------------------------------
#
function linear_or_log<std::vector<double>>
- void
- py_name linear_or_log
- io std::vector<double> &x
- out bool &log_x
#
# ------------------------------------------------------
#
class funct_string<double>
- py_name funct_string
- no_def_cons
- cons init
  - std::string expr
  - std::string var
- function set_parm
  - int
  - std::string name
  - double val
- function operator()
  - double
  - double x
#
# ------------------------------------------------------
#
class comm_option_s
- char shrt
- std::string lng
- std::string desc
- int min_parms
- int max_parms
- std::string parm_desc
- std::string help
- int type
#
# ------------------------------------------------------
#
class cmd_line_arg
- std::string arg
- bool is_option
- bool is_valid
- std::vector<std::string> parms
#- comm_option_s *cop
#
# ------------------------------------------------------
#
class cli
- bool sync_verbose
- bool gnu_intro
- std::string desc
- std::string cmd_name
- std::string addl_help_cmd
- std::string addl_help_cli
- function set_verbose
  - int
  - int v
- function parse_for_aliases
  - int
  - io std::vector<std::string> &sv
  - bool allow_undashed
  - bool debug [false]
- function apply_aliases
  - int
  - io std::vector<std::string> &sv
  - size_t istart
  - bool debug [false]
#- function get_parameter_list
#  - std::vector<std::string>
- function get_option_list
  - std::vector<std::string>
- function parameter_desc
  - std::string
  - std::string name
- function option_short_desc
  - std::string
  - std::string name
#
# ------------------------------------------------------
#
# Function from cursesw.h
#
function get_screen_size_ioctl
- int
- out int &row
- out int &col
