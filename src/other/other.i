# Interface file for o2scl classes and functions in src/other
# 
namespace o2scl
py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``%name%``,
| See
| https://neutronstars.utk.edu/code/o2scl/html/class/%name%.html .
dll_name o2scl
rst_header |
| .. _other:
|
| Other classes from O\ :sub:`2`\ scl
| ===================================
|
| :ref:`O2sclpy <o2sclpy>`
# 
# Include statements for C++ header file
# 
h_include <o2scl/slack_messenger.h>
h_include <o2scl/poly.h>
h_include <o2scl/polylog.h>
h_include <o2scl/hist.h>
h_include <o2scl/contour.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/other_python.h>
# 
# Namespace to use in C++ source code
# 
cpp_using std
cpp_using o2scl
#
# Additional python headers
#
py_header from o2sclpy.base import *
#
class slack_messenger
- int verbose
- std::string url                             
- std::string channel
- std::string icon
- std::string username
- double min_time_between
- cons init
  - std::string p_channel
  - std::string p_username
  - std::string p_url
  - bool p_mpi_time
- function set_url_from_env
  - bool
  - std::string env_var  
- function set_channel_from_env
  - bool
  - std::string env_var  
- function set_username_from_env
  - bool
  - std::string env_var
- function send
  - int
  - std::string message
  - bool err_on_fail [true]  
class quadratic_real_coeff_gsl
- function solve_r
  - int
  - const double a2
  - const double b2
  - const double c2
  - out double &r1
  - out double &r2  
- function solve_rc
  - int
  - const double a2
  - const double b2
  - const double c2
  - out std::complex<double> &r1  
  - out std::complex<double> &r2
class quadratic_real_coeff_gsl2<>
- py_name quadratic_real_coeff_gsl2                             
- function solve_r
  - int
  - const double a2
  - const double b2
  - const double c2
  - out double &r1
  - out double &r2  
- function solve_rc
  - int
  - const double a2
  - const double b2
  - const double c2
  - out std::complex<double> &r1  
  - out std::complex<double> &r2
class cubic_real_coeff_cern<>
- py_name cubic_real_coeff_cern                             
- function solve_r
  - int
  - const double a3
  - const double b3
  - const double c3
  - const double d3
  - out double &r1
  - out double &r2  
  - out double &r3
- function solve_rc
  - int
  - const double a3
  - const double b3
  - const double c3
  - const double d3
  - out double &r1
  - out std::complex<double> &r2  
  - out std::complex<double> &r3
class cubic_real_coeff_gsl
- function solve_r
  - int
  - const double a3
  - const double b3
  - const double c3
  - const double d3
  - out double &r1
  - out double &r2  
  - out double &r3
- function solve_rc
  - int
  - const double a3
  - const double b3
  - const double c3
  - const double d3
  - out double &r1
  - out std::complex<double> &r2  
  - out std::complex<double> &r3
class quartic_real_coeff_cern<>
- py_name quartic_real_coeff_cern                             
- function solve_r
  - int
  - const double a4
  - const double b4
  - const double c4
  - const double d4
  - const double e4
  - out double &r1
  - out double &r2  
  - out double &r3
  - out double &r4
- function solve_rc
  - int
  - const double a4
  - const double b4
  - const double c4
  - const double d4
  - const double e4
  - out std::complex<double> &r1
  - out std::complex<double> &r2  
  - out std::complex<double> &r3
  - out std::complex<double> &r4
class fermi_dirac_integ_gsl
- function calc_m1o2
  - double
  - double x
- function calc_1o2
  - double
  - double x
- function calc_3o2
  - double
  - double x
- function calc_2
  - double
  - double x
- function calc_3
  - double
  - double x
class bessel_K_exp_integ_gsl
- function K1exp
  - double
  - double x
- function K2exp
  - double
  - double x
- function K3exp
  - double
  - double x
# class polylog
# function set_tol
# - void
# - double tol
# function calc
# - double
# - double s
# - double y
class hist
- std_cc
#- cons init
#  - size_t nv
#  - io std::vector<double> &v
#  - size_t n_bins
#- function get_wgts
#  - const ublas_vector &
- function create_rep_vec
  - void
  - io std_vector &v                             
- function from_table
  - void                             
  - py_name from_table
  - io table<> &t
  - std::string colx
  - size_t n_bins
- function from_table
  - void                             
  - py_name from_table_twocol
  - io table<> &t
  - std::string colx
  - std::string coly
  - size_t n_bins
- function size
  - size_t
- bool extend_rhs
- bool extend_lhs
- function set_bin_edges
  - void
  - py_name set_bin_edges_grid
  - io uniform_grid<double> &g
- function set_bin_edges
  - void
  - py_name set_bin_edges_vec
  - size_t n
  - io vector<double> &v
- function update
  - void
  - double x
  - double val [1.0]
- function update_i
  - void
  - size_t i
  - double val [1.0]
- function get_wgt_i
  - double
  - size_t i
- function get_wgt
  - double
  - double x
- function set_wgt_i
  - void 
  - size_t i
  - double x
- function set_wgt
  - void
  - double x
  - double val
- function operator[]
  - const double &
  - size_t i
- function get_bin_index
  - size_t
  - double x
- function function
  - int
  - std::string func
- function clear
  - void
class contour_line
- double level
- std::vector<double> x  
- std::vector<double> y
- std_cc
class std::vector<contour_line>
- py_name vector_contour_line
- function operator[]
  - contour_line &
  - size_t n
- function resize
  - void
  - size_t n
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: an int
|     """
|     return self.length()
| 
class contour
- int verbose
- double lev_adjust
- bool debug_next_point
- function set_data
  - void
  - io const uniform_grid<double> &ugx
  - io const uniform_grid<double> &ugy
  - io const boost::numeric::ublas::matrix<double> &udata
- function set_levels
  - void
  - size_t n_levels
  - io vector<size_t> &levels
- function calc_contours
  - void
  - io vector<contour_line> &clines
