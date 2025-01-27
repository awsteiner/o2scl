# Interface file for o2scl nuclei classes
# 
namespace o2scl
py_class_doc |
| Python interface for O2scl class ``%name%``,
| See
| https://awsteiner.org/code/o2scl/html/class/%name%.html .
dll_name o2scl
rst_header |
| .. _nuclei:
|
| Nuclei and nuclear mass classes
| ===============================
|
| :ref:`O2sclpy <o2sclpy>`
# 
# Include statements for C++ header file
# 
h_include <o2scl/nucleus.h>
h_include <o2scl/nucmass.h>
h_include <o2scl/nucmass_ame.h>
h_include <o2scl/nucmass_dz.h>
h_include <o2scl/nucmass_frdm.h>
h_include <o2scl/nucmass_dglg.h>
h_include <o2scl/nucmass_sdnp.h>
h_include <o2scl/nucmass_wlw.h>
h_include <o2scl/nucmass_hfb.h>
h_include <o2scl/nucmass_ktuy.h>
h_include <o2scl/nucmass_fit.h>
h_include <o2scl/nucmass_gen.h>
h_include <o2scl/nucdist.h>
h_include <o2scl/hdf_nucmass_io.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/nuclei_python.h>
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
py_header from o2sclpy.part import *
#
# ------------------------------------------------------
# 
# Class nucleus
#
class nucleus
- parent part
- int Z
- int N
- int A
- double mex
- double be
#
# ------------------------------------------------------
# 
# Class nucmass_info
#
class nucmass_info
- function parse_elstring
  - int
  - std::string ela
  - out int &Z
  - out int &N
  - out int &A
- function eltoZ
  - int
  - std::string el
- function Ztoel
  - std::string
  - size_t Z
- function Ztoname
  - std::string
  - size_t Z
- function tostring
  - std::string
  - size_t Z
  - size_t N
- function int_to_spinp
  - std::string
  - int g
- function spinp_to_int
  - int
  - std::string s
#
# ------------------------------------------------------
# 
# Class nucmass
#
class nucmass abstract
- double m_prot
- double m_neut
- double m_elec
- double m_amu
- function is_included
  - bool
  - int Z
  - int N
- function get_nucleus
  - int
  - int Z
  - int N
  - nucleus &n
- function mass_excess
  - double
  - int Z
  - int N
- function mass_excess_d
  - double
  - double Z
  - double N
- function electron_binding
  - double
  - double Z
- function binding_energy
  - double
  - int Z
  - int N
- function binding_energy_d
  - double
  - double Z
  - double N
- function total_mass
  - double
  - int Z
  - int N
- function total_mass_d
  - double
  - double Z
  - double N
- function neutron_sep
  - double
  - int Z
  - int N
- function two_neutron_sep
  - double
  - int Z
  - int N
- function proton_sep
  - double
  - int Z
  - int N
- function two_proton_sep
  - double
  - int Z
  - int N
- function atomic_mass
  - double
  - int Z
  - int N
- function atomic_mass_d
  - double
  - double Z
  - double N
#
# ------------------------------------------------------
# 
# Class nucmass_table
#
class nucmass_table abstract
- parent nucmass
- size_t n
- std::string reference
- function is_loaded
  - bool
- function get_nentries
  - size_t
#
# ------------------------------------------------------
# 
# Class nucmass_fit_base
#
class nucmass_fit_base abstract
- parent nucmass
- size_t nfit
- function fit_fun
  - int
  - size_t nv
  - const io boost::numeric::ublas::vector<double> &x
- function guess_fun
  - int
  - size_t nv
  - out boost::numeric::ublas::vector<double> &x
#
# ------------------------------------------------------
# 
# Class nucmass_semi_empirical
#
class nucmass_semi_empirical
- parent nucmass_fit_base
- double B
- double Sv
- double Ss
- double Ec
- double Epair
- function mass_excess
  - double
  - int Z
  - int N
- function mass_excess_d
  - double
  - double Z
  - double N
#
# ------------------------------------------------------
# 
# Class nucmass_ame
#
class nucmass_ame
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_ame2
#
class nucmass_ame2
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_dz_table
#
class nucmass_dz_table
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_dz_fit
#
class nucmass_dz_fit
- parent nucmass_fit_base
#
# ------------------------------------------------------
# 
# Class nucmass_dz_fit_33
#
class nucmass_dz_fit_33
- parent nucmass_fit_base
#
# ------------------------------------------------------
# 
# Class nucmass_densmat
#
class nucmass_densmat abstract
- parent nucmass_fit_base
#
# ------------------------------------------------------
# 
# Class nucmass_frdm
#
class nucmass_frdm
- parent nucmass_densmat
- double a1
- double J
- double K
- double a2
- double Q
- double a3
- double ca
- double W
- double ael
- double rp
- double r0
- double MH
- double Mn
- double e2
- double a
- double aden
- double rmac
- double h
- double L
- double C
- double gamma
- double amu
- double nn
- double np
- double Rn
- double Rp
#
# ------------------------------------------------------
# 
# Class nucmass_mnmsk
#
class nucmass_mnmsk
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_mnmsk_exp
#
class nucmass_mnmsk_exp
- parent nucmass_mnmsk
#
# ------------------------------------------------------
# 
# Class nucmass_gen
#
class nucmass_gen
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_dglg
#
class nucmass_dglg
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_hfb
#
class nucmass_hfb
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_hfb_sp
#
class nucmass_hfb_sp
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_ktuy
#
class nucmass_ktuy
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_sdnp
#
class nucmass_sdnp
- parent nucmass_table
#
# ------------------------------------------------------
# 
# Class nucmass_wlw
#
class nucmass_wlw
- parent nucmass_table
- function load
  - int
  - std::string model [""]
  - bool external [false]
#
# ------------------------------------------------------
# 
# Class nucmass_fit
#
class nucmass_fit
- int fit_method
- static const int rms_mass_excess
- static const int rms_binding_energy
- static const int chi_squared_me
- static const int chi_squared_be
- bool even_even
- int minZ
- int minN
- function fit
  - void
  - io nucmass_fit_base &n
  - out double &res
- function eval
  - void
  - io nucmass &n
  - out double &res
- function fit_covar
  - void
  - io nucmass_fit_base &n
  - out double &chi2
  - out ubmatrix &covar
#
# ------------------------------------------------------
#
# Class vector<nucleus>
#                              
class std::vector<nucleus>
- py_name std_vector_nucleus
- std_cc                             
- function resize
  - void
  - size_t n                             
- function size
  - size_t
- function operator[]
  - nucleus &
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
# HDF functions
#
# ------------------------------------------------------
#
function ame_load
- void
- nucmass_ame &ame
- std::string name ["20"]
- bool exp_only [False]
#
# ------------------------------------------------------
#
function ame_load_ext
- void
- nucmass_ame &ame
- std::string file_name
- std::string table_name
- bool exp_only [False]
#
# ------------------------------------------------------
#
function mnmsk_load
- void
- nucmass_mnmsk &mnmsk
- std::string model [""]
- std::string filename [""]
#
# ------------------------------------------------------
#
function hfb_load
- void
- nucmass_hfb &hfb
# We can't specify the default parameter for model because
# we cannot specify a default parameter for string objects.
- size_t model
- std::string filename
#
# ------------------------------------------------------
#
function hfb_sp_load
- void
- nucmass_hfb_sp &hfb
# We can't specify the default parameter for model because
# we cannot specify a default parameter for string objects.
- size_t model
- std::string filename
#
# ------------------------------------------------------
#
function nucdist_set
- void
- vector<nucleus> &dist
- nucmass &nm
- std::string expr ["1"]  
- int maxA [400]
- bool include_neutron [false]
#
# ------------------------------------------------------
#
function nucdist_pair_set
- void
- vector<nucleus> &dist
- nucmass &nm
- nucmass &nm2
- std::string expr ["1"]  
- int maxA [400]
- bool include_neutron [false]
#
# ------------------------------------------------------
#
function nucdist_set_ext
- void
- vector<nucleus> &dist
- vector<nucleus> &dist_ext
- nucmass &nm
- std::string expr ["1"]  
- int maxA [400]
- int n_chop [1]

