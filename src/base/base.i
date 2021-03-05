# Interface file for o2scl base classes
# 
namespace o2scl
py_class_doc |
| Python interface for C++ class ``%name%``.
dll_name o2scl
rst_header |
| .. _base:
|
| Base classes from O\ :sub:`2`\ scl
| ==================================
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
# (none)
#
# Class string
# 
class std::string
- py_name std_string
- function length
  - size_t
- function operator[]
  - char &
  - size_t n
- function resize
  - void
  - size_t n
- extra_py |
| def __len__(self):
|     return length()
| 
| def init_bytes(self,s):
|     """
|     Initialize the string from a Python bytes object
|     """
|     self.resize(len(s))
|     for i in range(0,len(s)):
|         self.__setitem__(i,s[i])
|     return
|
| def to_bytes(self):
|     """
|     Copy the string to a Python bytes object
|     """
|     ret=b''
|     for i in range(0,self.length()):
|         ret=ret+self.__getitem__(i)
|     return ret
#
# Class vector<double>
#                              
# Create a python interface to std::vector<double> for vector
# arguments to O2scl functions                             
#
class std::vector<double>
- py_name std_vector
- function resize
  - void
  - size_t n                             
- function size
  - size_t
- function operator[]
  - double &
  - size_t n
- extra_py |
| def __len__(self):
|     return size()
|
| def to_numpy(self):
|     """
|     Copy the vector to a numpy array
|     """
|     ret=numpy.zeros((self.size()))
|     for i in range(0,self.size()):
|         ret[i]=self.__getitem__(i)
|     return ret
class vector<int>
- py_name std_vector_int
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
|     return size()
|
| def to_numpy(self):
|     """
|     Copy the vector to a numpy array
|     """
|     ret=numpy.zeros((self.size()),dtype=numpy.int32_t)
|     for i in range(0,self.size()):
|         ret[i]=self.__getitem__(i)
|     return ret
class vector<size_t>
- py_name std_vector_size_t
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
|     Test
|     """
|     return size()
|
| def to_numpy(self):
|     """
|     Copy the vector to a numpy array
|     """
|     ret=numpy.zeros((self.size()),dtype=numpy.uint64_t)
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
#class vector<string>
#- py_name std_vector_string
#- function resize
#  - void
#  - size_t n                             
#- function size
#  - size_t
#- function operator[]
#  - string
#  - size_t n
#- extra_py |
#| def __len__(self):
#|     return size()
#
# Class ublas_vector
# 
class boost::numeric::ublas::vector<double>
- py_name ublas_vector
- function size
  - size_t
- function resize
  - void
  - size_t n
- function operator[]
  - double &
  - size_t i
#
# Class ublas_matrix
# 
class boost::numeric::ublas::matrix<double>
- py_name ublas_matrix
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
#
# -------------------------------------------------------------------
#
# Set the python class documentation for the following classes
# 
py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``%name%``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/%name%.html .
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
- function eos_installed
  - bool
- function part_installed
  - bool
- function hdf_support
  - bool
- function openmp_support
  - bool
- function readline_support
  - bool
- function ncurses_support
  - bool
- function gsl2_support
  - bool
- function armadillo_support
  - bool
- function eigen_support
  - bool
- function fftw_support
  - bool
- function python_support
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
# Class table
#
class table<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``table``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/table.html .
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
- function get_column
  - vector<double> &
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
  - table<> &source
  - std::string src_index
  - std::string src_col
  - std::string dest_index
  - std::string dest_col
- function insert_table
  - void
  - table<> &source
  - std::string src_index
  - bool allow_extrap
  - std::string dest_index
- function add_table
  - void
  - table<> &source
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
  - std_vector &data
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
|     vec=std_vector(self._link)
|     vec.resize(len(v))
|     for i in range(0,len(v)):
|         vec[i]=v[i]
|     self.line_of_data_vector(vec)
|     return
# 
# Class table_units<>
#
class table_units<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``table_units``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/table_units.html .
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
  - bool err_on_fail
# 
# Shared pointer for table_units<> objects
#
shared_ptr table_units<>
- py_name table_units
# 
# Class uniform_grid
#
class uniform_grid<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``uniform_grid``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/uniform_grid.html .
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
  - std::vector<double> &v
# 
# Class uniform_grid_end
#
class uniform_grid_end<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``uniform_grid_end``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/uniform_grid_end.html .
- py_name uniform_grid_end                             
- parent uniform_grid<>
- cons init
  - double start
  - double end
  - size_t n_bins                             
# 
# Class uniform_grid_width
#
class uniform_grid_width<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``uniform_grid_width``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/uniform_grid_width.html .
- py_name uniform_grid_width
- parent uniform_grid<>
- cons init
  - double start
  - double width
  - size_t n_bins                             
# 
# Class uniform_grid_end_width
#
class uniform_grid_end_width<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``uniform_grid_end_width``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/uniform_grid_end_width.html .
- py_name uniform_grid_end_width                             
- parent uniform_grid<>
- cons init
  - double start
  - double end
  - double width
# 
# Class uniform_grid_log_end
#
class uniform_grid_log_end<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``uniform_grid_log_end``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/uniform_grid_log_end.html .
- py_name uniform_grid_log_end
- parent uniform_grid<>
- cons init
  - double start
  - double end
  - size_t n_bins                             
# 
# Class uniform_grid_log_width
#
class uniform_grid_log_width<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``uniform_grid_log_width``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/uniform_grid_log_width.html .
- py_name uniform_grid_log_width
- parent uniform_grid<>
- cons init
  - double start
  - double width
  - size_t n_bins                             
# 
# Class uniform_grid_log_end_width
#
class uniform_grid_log_end_width<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``uniform_grid_log_end_width``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/uniform_grid_log_end_width.html .
- py_name uniform_grid_log_end_width
- parent uniform_grid<>
- cons init
  - double start
  - double end
  - double width
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
  - std_vector &x
  - std::string y_name
  - size_t ny
  - std_vector &y
- function set_xy
  - void
  - py_name set_xy_grid
  - std::string x_name
  - uniform_grid<double> &x_grid
  - std::string y_name
  - uniform_grid<double> &y_grid
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
- function get_size
  - void
  - size_t &nx
  - size_t &ny
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
  - size_t &ix
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
  - boost::numeric::ublas::matrix<double> &mat
  - bool throw_on_err
- function function_slice
  - void
  - std::string function
  - std::string slice
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
#- function ix_index
#  - index_spec
#  - size_t ix
# 
# Class tensor
#
class tensor<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``tensor``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/tensor.html .
- std_cc                             
- py_name tensor
#- cons create_size
#  - size_t rank
#  - vector<size_t> &sizes
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
- function get
  - double
  - py_name get_vector
  - vector<size_t> &index
#- function resize
#  - void
#  - py_name resize_vector
#  - vector<size_t> &index
- function get_rank
  - size_t
- function get_size
  - size_t
  - size_t i
#- function get_data
#  - const vector<double> &
- function total_size
  - size_t
- function min_value
  - double
- function max_value
  - double
- function total_sum
  - double
- extra_py |
| def set(self,index,val):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object and add the 
|     data to the table
|     """
|     svst=o2sclpy.std_vector_size_t(self._link)
|     syst.init_py(index)
|     self.set_vector(syst,val)
|     return
| 
| def get(self,index):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object and get the 
|     data from the table
|     """
|     svst=o2sclpy.std_vector_size_t(self._link)
|     syst.init_py(index)
|     return self.get_vector(syst)
|
| def resize(self,index):
|     """
|     Copy ``index`` to an :class:`std_vector_size_t` object 
|     and resize
|     """
|     svst=o2sclpy.std_vector_size_t(self._link)
|     syst.init_py(index)
|     self.resize_vector(syst)
|     return
# 
# Class tensor_grid
#
class tensor_grid<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``tensor_grid``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/tensor_grid.html .
- std_cc                             
- py_name tensor_grid
#- cons create_size
#  - size_t rank
#  - vector<size_t> &sizes
- function is_valid
  - void
- function set_val
  - void
  - py_name set_val_vector
  - vector<double> &grid_point
  - double val
- function get_val
  - double
  - py_name get_val_vector
  - vector<double> &grid_point
- function is_grid_set
  - bool
- function set_grid_packed
  - void
  - vector<double> &grid
- function default_grid
  - void
- function set_grid_i_vec    
  - void
  - size_t i
  - vector<double> &grid
- function get_grid
  - double
  - size_t i
  - size_t j
- function set_grid
  - void
  - size_t i
  - size_t j
  - double val
# 
# Class find_constants
#
class find_constants
- function find_print
  - void
  - std::string name
  - std::string unit
  - size_t prec
  - int verbose
- function find_unique
  - double
  - std::string name
  - std::string unit
# 
# Class convert_units<>
#
# Note that 'from' is a reserved keyword in python so we
# rename it to 'frm' instead
#
class convert_units<>
- py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``convert_units``,
| see
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/convert_units.html .
- py_name convert_units
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
- int verbose
- bool use_gnu_units
- bool err_on_fail
- bool combine_two_conv
- std::string units_cmd_string
- function print_cache
  - void

