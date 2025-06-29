# Interface file for o2scl classes and functions in src/other
# 
namespace o2scl
py_class_doc |
| Python interface for O2scl class ``%name%``,
| See
| https://awsteiner.org/code/o2scl/html/class/%name%.html .
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
h_include <o2scl/hist_2d.h>
h_include <o2scl/contour.h>
h_include <o2scl/vec_stats.h>
h_include <o2scl/prob_dens_func.h>
h_include <o2scl/prob_dens_mdim_amr.h>
h_include <o2scl/root_brent_gsl.h>
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
#
# ------------------------------------------------------
#
# Class slack_messenger
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
#
# ------------------------------------------------------
#
# Class root_brent_gsl
# 
# class root_brent_gsl<funct_string<double>,double>
# - py_name root_brent_gsl
# - double tol_rel
# - double tol_abs
# - int verbose
# - int ntrial
# - bool err_nonconv                               
# - function solve_bkt
#   - int
#   - io double &x1
#   - double x2                               
#   - io funct_string &fp                               
#
# ------------------------------------------------------
#
# Class quadratic_real_coeff_gsl
# 
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
#
# ------------------------------------------------------
#
# Class quadratic_real_coeff_gsl_gsl2
# 
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
#
# ------------------------------------------------------
#
# Class cubic_real_coeff_cern
# 
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
#
# ------------------------------------------------------
#
# Class cubic_real_coeff_gsl
# 
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
#
# ------------------------------------------------------
#
# Class quartic_real_coeff_cern
# 
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
#
# ------------------------------------------------------
#
# Class fermi_dirac_integ_gsl
# 
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
#
# ------------------------------------------------------
#
# Class bessel_K_exp_integ_gsl
# 
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
#                               
# class polylog
# function set_tol
# - void
# - double tol
# function calc
# - double
# - double s
# - double y
#
#
# ------------------------------------------------------
#
# Class hist
# 
class hist
- std_cc
#- cons init
#  - size_t nv
#  - io std::vector<double> &v
#  - size_t n_bins
- function create_rep_vec
  - void
  - io std_vector &v
- function get_wgts
  - const boost::numeric::ublas::vector<double> &
- function get_bins
  - const boost::numeric::ublas::vector<double> &
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
#
# ------------------------------------------------------
#
# Class hist_2d
# 
class hist_2d
- std_cc
- function create_x_rep_vec
  - void
  - io std_vector &v
- function create_y_rep_vec
  - void
  - io std_vector &v
#- function get_wgts
#  - const boost::numeric::ublas::matrix<double> &
- function from_table
  - void                             
  - py_name from_table
  - io table<> &t
  - std::string colx
  - std::string coly
  - size_t n_bins_x
  - size_t n_bins_y
- function from_table
  - void                             
  - py_name from_table_wgt
  - io table<> &t
  - std::string colx
  - std::string coly
  - std::string colz
  - size_t n_bins_x
  - size_t n_bins_y
- function size_x
  - size_t
- function size_y
  - size_t
- bool extend_rhs
- bool extend_lhs
- function set_bin_edges
  - void
  - py_name set_bin_edges_grid
  - io uniform_grid<double> &gx
  - io uniform_grid<double> &gy
- function set_bin_edges
  - void
  - py_name set_bin_edges_vec
  - size_t nx
  - io vector<double> &vx
  - size_t ny
  - io vector<double> &vy
- function update
  - void
  - double x
  - double y
  - double val [1.0]
- function update_i
  - void
  - size_t i
  - size_t j
  - double val [1.0]
- function get_wgt_i
  - double
  - size_t i
  - size_t j
- function get_wgt
  - double
  - double x
  - double y
- function set_wgt_i
  - void 
  - size_t i
  - size_t j
  - double val
- function get_x_low_i
  - double
  - size_t i
- function get_x_high_i
  - double
  - size_t i
- function get_y_low_i
  - double
  - size_t i
- function get_y_high_i
  - double
  - size_t i
- function set_wgt
  - void
  - double x
  - double y
  - double val
#- function get_bin_indices
#  - size_t
#  - double x
#  - double y
#  - out size_t &ix
#  - out size_t &iy
- function get_wgts
  - boost::numeric::ublas::matrix<double> &
- function clear
  - void
- function clear_wgts
  - void
#
# ------------------------------------------------------
#
# Class contour_line
# 
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
- function size
  - size_t                        
- extra_py |
| def __len__(self):
|     """
|     Return the length of the vector
|
|     Returns: an int
|     """
|     return self.length()
| 
#
# ------------------------------------------------------
#
# Class contour
# 
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
#
# ------------------------------------------------------
#
# Class prob_dens_func
# 
class prob_dens_func
- function pdf
  - double
  - double x                               
- function log_pdf
  - double
  - double x                               
- function cdf
  - double
  - double x                               
- function invert_cdf
  - double
  - double x                               
- function entropy
  - double
- function operator()
  - double
- function sample
  - double
#
# ------------------------------------------------------
#
# Class prob_dens_gaussian
# 
class prob_dens_gaussian
- parent prob_dens_func
- function set_center
  - void
  - double cent
- function set_sigma
  - void
  - double sigma
#
# ------------------------------------------------------
#
# Class prob_dens_hist
# 
class prob_dens_hist
- parent prob_dens_func
- function init
  - void
  - io hist &h
#
# ------------------------------------------------------
#
# Class prob_dens_mdim
# 
class prob_dens_mdim<std::vector<double>>
- py_name prob_dens_mdim
- function pdf
  - double
  - io std_vector &x
- function log_pdf
  - double
  - io std_vector &x
- function dim
  - size_t
# The operator() class is parsed by yanic as a __getitem__() magic
# method, which is ok because this implies a similar syntax compared
# to the C++ usage.
- function operator()
  - void
  - io std_vector &x
#
# ------------------------------------------------------
#
# Class prob_dens_mdim_biv_gaussian
# 
class prob_dens_mdim_biv_gaussian<std::vector<double>>
- py_name prob_dens_mdim_biv_gaussian
- parent prob_dens_mdim<std::vector<double>>
- function set
  - void
  - double x_cent
  - double y_cent
  - double x_std
  - double y_std
  - double covar
- function get
  - void
  - out double &x_cent
  - out double &y_cent
  - out double &x_std
  - out double &y_std
  - out double &covar
- function level_fixed_integral
  - double
  - double integral
class prob_dens_mdim_gaussian<>
- py_name prob_dens_mdim_gaussian
- parent prob_dens_mdim<std::vector<double>>
- function make_biv
  - prob_dens_mdim_biv_gaussian<>
#    
#class std::vector<cprob_dens_mdim_gaussian<>>
#- py_name vector_prob_dens_mdim_gaussian
#- function operator[]
#  - prob_dens_mdim_gaussian &
#  - size_t n
#- function resize
#  - void
#  - size_t n
#- function size
#  - size_t                        
#- extra_py |
#| def __len__(self):
#|     """
#|     Return the length of the vector
#|
#|     Returns: an int
#|     """
#|     return self.length()
#|
#
#
# ------------------------------------------------------
#
# Class hypercube
# 
class prob_dens_mdim_amr<>::hypercube
- py_name hypercube
- size_t n_dim
- std::vector<double> low
- std::vector<double> high
- std::vector<size_t> inside
- double frac_vol
- double weight  
#
# ------------------------------------------------------
#
# Class std_vector_hypercube
# 
class std::vector<prob_dens_mdim_amr<>::hypercube>
- py_name std_vector_hypercube
- function resize
  - void
  - size_t n                             
- function size
  - size_t
- function operator[]
  - prob_dens_mdim_amr<>::hypercube &
  - size_t n
#
# ------------------------------------------------------
#
# Class prob_dens_mdim_amr
# 
class prob_dens_mdim_amr<>
- py_name prob_dens_mdim_amr
- parent prob_dens_mdim<std::vector<double>>
- std::vector<prob_dens_mdim_amr<>::hypercube> mesh
- size_t n_dim
- std::vector<double> low
- std::vector<double> high
- bool allow_resampling
- std::vector<double> scale
- int verbose
- function clear
  - void
- function clear_mesh
  - void
- function total_volume
  - double
#
# ------------------------------------------------------
#
# Functions from vec_stats.h
# 
function vector_mean<std::vector<double>,double>
- double
- py_name vector_mean
- size_t n  
- io const vector<double> &v
function vector_stddev<std::vector<double>>
- double
- py_name vector_stddev
- size_t n  
- io const vector<double> &v
function vector_lagk_autocorr<std::vector<double>>
- double
- py_name vector_lagk_autocorr
- size_t n  
- io const vector<double> &v
- size_t k
function vector_autocorr_vector<std::vector<double>,std::vector<double>>
- void
- py_name vector_autocorr_vector
- size_t n  
- io const vector<double> &v
- out vector<double> &ac
function vector_autocorr_vector_fftw<std::vector<double>,std::vector<double>>
- void
- py_name vector_autocorr_vector_fftw
- io const vector<double> &v
- out vector<double> &ac
- double mean
- double stddev
function vector_autocorr_tau<std::vector<double>,std::vector<double>>
- size_t
- py_name vector_autocorr_tau
- io const vector<double> &ac
- out vector<double> &ftom
function vector_acor<std::vector<double>>
- void
- py_name vector_acor
- size_t n  
- io const vector<double> &v
- double mean
- double sigma
- io double &tau
