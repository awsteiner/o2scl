# Interface file for o2scl base classes
# 
namespace o2scl
py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``%name%``,
| See
| https://neutronstars.utk.edu/code/o2scl-dev/html/class/%name%.html .
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
# py_header from o2sclpy.part import *
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
- py_name table
- std_cc                             
- function operator[]
  - vector<double> &
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
# - function get_column
#  - vector<double> &
# - std::string col
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
# - function init_column
#  - void
#  - std::string scol
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
# 
# Class table_units<>
#
class table_units<>
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
# Class table3d
#
class table3d
- std_cc
#- function set_xy
#  - void
#  - std::string x_name
#  - size_t nx
#  - vector &x
#  - std::string y_name
#  - size_t ny
#  - vector &y
- function set
  - void
  - size_t ix
  - size_t iy
  - std::string name
  - double val
- function get
  - double
  - size_t ix
  - size_t iy
  - std::string name
- function new_slice
  - void
  - std::string slice
- function get_nx
  - size_t
- function get_ny
  - size_t
- function get_nslices
  - size_t
# 
# Class tensor
#
class tensor<>
- std_cc                             
- py_name tensor
- function clear
  - void
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

