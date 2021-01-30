# Interface file for o2scl base classes
# 
namespace o2scl
py_class_doc_pattern "Python interface for class :ref:`%name% <o2scl:%name%>`."
# 
# Include statements for C++ header file
# 
h_include <o2scl/table_units.h>
h_include <o2scl/table3d.h>
h_include <o2scl/tensor.h>
h_include <o2scl/tensor_grid.h>
h_include <o2scl/find_constants.h>
h_include <o2scl/convert_units.h>
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

