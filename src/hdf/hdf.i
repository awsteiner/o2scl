# Interface file for o2scl HDF5 classes
# 
namespace o2scl_hdf
py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``%name%``,
| See
| https://neutronstars.utk.edu/code/o2scl/html/class/%name%.html .
dll_name o2scl_hdf
rst_header |
| .. _hdf:
|
| HDF5 classes from O\ :sub:`2`\ scl
| ==================================
# 
# Include statements for C++ header file
# 
h_include <o2scl/table.h>
h_include <o2scl/hdf_file.h>
h_include <o2scl/hdf_io.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/hdf_python.h>
# 
# Namespace to use in C++ source code
# 
cpp_using std
cpp_using o2scl
cpp_using o2scl_hdf
#
# Additional python headers
#
py_header from o2sclpy.base import *
#
class hdf_file
- function open
  - void
  - std::string fname
  - bool write_access [false]
  - bool err_on_fail [true]
- function open_or_create
  - void
  - std::string fname
- function close
  - void
function hdf_input
- void                             
- py_name hdf_input_table
- hdf_file &hf
- table<> &t
- std::string name  
function hdf_input_n
- void                             
- py_name hdf_input_n_table
- hdf_file &hf
- table<> &t
- std::string &name  
function hdf_output
- void                             
- py_name hdf_output_table
- hdf_file &hf
- table<> &t
- std::string name  
function hdf_input
- void                             
- py_name hdf_input_table_units
- hdf_file &hf
- table_units<> &t
- std::string name  
function hdf_input_n
- void                             
- py_name hdf_input_n_table_units
- hdf_file &hf
- table_units<> &t
- std::string &name  
function hdf_output
- void                             
- py_name hdf_output_table_units
- hdf_file &hf
- table_units<> &t
- std::string name  
function hdf_input
- void                             
- py_name hdf_input_table3d
- hdf_file &hf
- table3d &t
- std::string name  
function hdf_input_n
- void                             
- py_name hdf_input_n_table3d
- hdf_file &hf
- table3d &t
- std::string &name  
function hdf_output
- void                             
- py_name hdf_output_table3d
- hdf_file &hf
- const table3d &t
- std::string name  
function hdf_input
- void                             
- py_name hdf_input_uniform_grid
- hdf_file &hf
- uniform_grid<> &t
- std::string name  
function hdf_input_n
- void                             
- py_name hdf_input_n_uniform_grid
- hdf_file &hf
- uniform_grid<> &t
- std::string &name  
function hdf_output
- void                             
- py_name hdf_output_uniform_grid
- hdf_file &hf
- uniform_grid<> &t
- std::string name  
