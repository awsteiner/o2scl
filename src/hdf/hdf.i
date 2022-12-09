# Interface file for o2scl HDF5 classes
# 
namespace o2scl_hdf
py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``%name%``,
| See
| https://neutronstars.utk.edu/code/o2scl/html/class/%name%.html .
dll_name o2scl
rst_header |
| .. _hdf:
|
| HDF5 classes from O\ :sub:`2`\ scl
| ==================================
|
| :ref:`O2sclpy <o2sclpy>`
# 
# Include statements for C++ header file
# 
h_include <o2scl/table.h>
h_include <o2scl/hdf_file.h>
h_include <o2scl/hdf_io.h>
h_include <o2scl/acolm.h>
h_include <o2scl/cloud_file.h>
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
cpp_using o2scl_acol
#
# Additional python headers
#
py_header from o2sclpy.base import *
py_header from o2sclpy.other import *
#
#
# Class hdf_file
#
class hdf_file
- int compr_type
- size_t min_compr_size  
- function has_write_access
  - bool                             
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
- function getc
  - int
  - std::string name
  - out char &c
- function getd
  - int
  - std::string name
  - out double &d
- function geti
  - int
  - std::string name
  - out int &i
- function get_szt
  - int
  - std::string name
  - out size_t &u
- function gets
  - int
  - std::string name
  - out std::string &s
- function gets_var
  - int
  - std::string name
  - out std::string &s
- function gets_fixed
  - int
  - std::string name
  - out std::string &s
# Note we rename the second parameter "deft" instead of "def" since
# def is a python keyword and "default" is a C++ keyword.
- function gets_def_fixed
  - int
  - std::string name
  - std::string deft
  - out std::string &s
- function setc
  - void
  - std::string name
  - char c
- function setd
  - void
  - std::string name
  - double d
- function seti
  - void
  - std::string name
  - int i
- function set_szt
  - void
  - std::string name
  - size_t u
- function sets
  - void
  - std::string name
  - std::string s
- function sets_fixed
  - void
  - std::string name
  - std::string s
- function getd_vec
  - int
  - std::string name
  - out std::vector<double> &v
- function geti_vec
  - int
  - std::string name
  - out std::vector<int> &v
- function get_szt_vec
  - int
  - std::string name
  - out std::vector<size_t> &v
- function gets_vec_copy
  - int
  - std::string name
  - out std::vector<std::string> &s
- function setd_vec
  - int
  - std::string name
  - io std::vector<double> &v
- function seti_vec
  - int
  - std::string name
  - io std::vector<int> &v
- function set_szt_vec
  - int
  - std::string name
  - io std::vector<size_t> &v
- function sets_vec_copy
  - int
  - std::string name
  - io std::vector<std::string> &s
- function getd_mat_copy
  - int
  - std::string name
  - out boost::numeric::ublas::matrix<double> &m
- function geti_mat_copy
  - int
  - std::string name
  - out boost::numeric::ublas::matrix<int> &m
- function setd_mat_copy
  - int
  - std::string name
  - io boost::numeric::ublas::matrix<double> &m
- function seti_mat_copy
  - int
  - std::string name
  - io boost::numeric::ublas::matrix<int> &m
- function getd_ten
  - int
  - std::string name
  - out tensor<> &t
- function geti_ten
  - int
  - std::string name
  - out tensor<int> &t
- function get_szt_ten
  - int
  - std::string name
  - out tensor<size_t> &t
- function setd_ten
  - int
  - std::string name
  - io tensor<> &t
- function seti_ten
  - int
  - std::string name
  - io tensor<int> &t
- function set_szt_ten
  - int
  - std::string name
  - io tensor<size_t> &t
- function getc_def
  - int
  - std::string name
  - char deft
  - out char &c
- function getd_def
  - int
  - std::string name
  - double deft
  - out double &d
- function geti_def
  - int
  - std::string name
  - int deft
  - out int &i
- function get_szt_def
  - int
  - std::string name
  - size_t deft
  - out size_t &u
- function gets_def
  - int
  - std::string name
  - std::string deft
  - out std::string &s
- function find_object_by_type
  - int
  - std::string type
  - out std::string &name
  - bool use_regex [false]                             
  - int verbose [0]
- function find_object_by_name
  - int
  - std::string name
  - out std::string &type
  - bool use_regex [false]                             
  - int verbose [0]
- function find_object_by_pattern
  - int
  - std::string pattern
  - out std::string &type
  - bool use_regex [false]                             
  - int verbose [0]
- function file_list
  - void
  - int verbose
- function copy
  - void
  - int verbose
  - out hdf_file &hf2
#
# Functions from hdf_io.h
#
function hdf_input
- void                             
- py_name hdf_input_table
- hdf_file &hf
- table<> &t
- std::string name [""]  
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
- std::string name [""]
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
- std::string name [""]
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
- std::string name [""]
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
function hdf_input
- void                             
- py_name hdf_input_tensor_grid
- hdf_file &hf
- tensor_grid<> &t
- std::string name [""]
function hdf_input_n
- void                             
- py_name hdf_input_n_tensor_grid
- hdf_file &hf
- tensor_grid<> &t
- std::string &name  
function hdf_output
- void                             
- py_name hdf_output_tensor_grid
- hdf_file &hf
- tensor_grid<> &t
- std::string name
function hdf_input
- void                             
- py_name hdf_input_vector_contour_line
- hdf_file &hf
- std::vector<contour_line> &v
- std::string name [""]
function hdf_input_n
- void                             
- py_name hdf_input_n_vector_contour_line
- hdf_file &hf
- std::vector<contour_line> &v
- std::string &name  
function hdf_output
- void                             
- py_name hdf_output_vector_contour_line
- hdf_file &hf
- std::vector<contour_line> &v
- std::string name
function value_spec
- int
- std::string spec
- out double &d
- int verbose [0]
- bool err_on_fail [true]
function vector_spec<std::vector<double>>
- int
- py_name vector_spec  
- std::string spec
- out std::vector<double> &v
- int verbose [0]
- bool err_on_fail [true]
function strings_spec<std::vector<std::string>>
- int
- py_name strings_spec  
- std::string spec
- out std::vector<std::string> &v
- int verbose [0]
- bool err_on_fail [true]
function vector_spec
- std::vector<double>
- py_name vector_spec_v
- std::string spec  
function mult_vector_spec<std::vector<double>>
- int
- py_name mult_vector_spec  
- std::string spec
- out std::vector<std::vector<double>> &v
- bool use_regex [false]                             
- int verbose [0]
- bool err_on_fail [true]
#
# Class acol_manager
#
class acol_manager
- std::string env_var_name
- cli *cl
- int verbose
- std::string def_args
- std::string type
- table_units<> table_obj
- table3d table3d_obj
- hist hist_obj
- hist_2d hist_2d_obj
- int int_obj
- char char_obj
- double double_obj
- size_t size_t_obj
- std::string string_obj
- std::vector<contour_line> cont_obj
- uniform_grid<double> ug_obj
- std::vector<int> intv_obj
- std::vector<double> doublev_obj
- std::vector<size_t> size_tv_obj
- std::vector<std::string> stringv_obj
- std::vector<std::vector<std::string>> vvstring_obj
- std::vector<std::vector<double>> vvdouble_obj
- tensor<> tensor_obj
- tensor<int> tensor_int_obj
- tensor<size_t> tensor_size_t_obj
- tensor_grid<> tensor_grid_obj
- prob_dens_mdim_amr<> pdma_obj  
- prob_dens_mdim_gaussian<> pdmg_obj  
- std::string command_color
- std::string type_color
- std::string param_color
- std::string help_color
- std::string exec_color
- std::string url_color
- std::string default_color
- std::string color_spec
- function help_found
  - bool
  - std::string arg1
  - std::string arg2
- function run_empty
  - int
- function validate_interp_type
  - int
- function parse_vec_string
  - void
  - io std::vector<std::string> &args
- function command_add
  - void
  - std::string new_type    
- function command_del
  - void
  - std::string ltype    
#- function run
#  - int argv [0]
#  - char *argc[] [0]
#  - bool full_process [true]
#
# Class cloud_file
#
class cloud_file
- int hash_type
- int verbose
- bool throw_on_fail
- bool allow_wget
- bool allow_curl
- function get_file
  - int
  - std::string file
  - std::string url
  - std::string dir [""]
- function get_file_hash
  - int
  - std::string file
  - std::string url
  - std::string hash [""]
  - std::string dir [""]
- function hdf5_open    
  - int
  - io hdf_file &hf
  - std::string file
  - std::string url
  - std::string dir [""]
- function hdf5_open_hash
  - int
  - io hdf_file &hf
  - std::string file
  - std::string url
  - std::string hash [""]
  - std::string dir [""]
