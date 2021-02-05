# Interface file for o2scl base classes
# 
namespace o2scl
py_class_doc_pattern "Python interface for O\ :sub:`2`\ scl class :ref:`%name% <o2scl:%name%>`."
dll_name o2scl
rst_name Base
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
- function eos_installed
  - bool
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
- function new_column
  - void
  - std::string col
- function get_column_name
  - std::string
  - size_t icol
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

