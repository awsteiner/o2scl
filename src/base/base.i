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
- function clear
  - void
- function clear_data
  - void
- function clear_table
  - void
- function clear_constants
  - void
- function operator[]
  - vector<double> &
  - std::string col
# 
# Class table_units<>
#
class table_units<>
- parent table<>
- py_name table_units
- function get_unit
  - std::string
  - std::string col
- function set_unit
  - void
  - std::string col
  - std::string unit
- function convert_to_unit
  - int
  - std::string col
  - std::string unit
  - bool err_on_fail
- function clear_table
  - void
# 
# Shared pointer for table_units<> objects
#
shared_ptr table_units<>
- py_name table_units
# 
# Class table3d
#
class table3d
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

