# Interface file for o2scl nuclei classes
# 
namespace o2scl
py_class_doc |
| Python interface for O\ :sub:`2`\ scl class ``%name%``,
| See
| https://neutronstars.utk.edu/code/o2scl-dev/part/html/class/%name%.html .
dll_name o2scl_part
rst_header |
| .. _nuclei:
|
| Nuclei and nuclear mass classes
| ===============================
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
py_header from o2sclpy.part import *
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
# Class nucmass_fit_base
#
class nucmass_fit_base abstract
- parent nucmass
- size_t nfit
#- function fit_fun
#  - int
#  - const ubvector &x
#- function guess_fun
#  - int
#  - ubvector &x
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
# Class nucmass_ame
#
class nucmass_ame
- parent nucmass_table
# 
# Class nucmass_dz_table
#
class nucmass_dz_table
- parent nucmass_table
# 
# Class nucmass_dz_fit
#
class nucmass_dz_fit
- parent nucmass_fit_base
# 
# Class nucmass_dz_fit_33
#
class nucmass_dz_fit_33
- parent nucmass_fit_base
# 
# Class nucmass_frdm
#
class nucmass_frdm
- parent nucmass_fit_base
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
# Class nucmass_mnmsk
#
class nucmass_mnmsk
- parent nucmass_table
# 
# Class nucmass_mnmsk_exp
#
class nucmass_mnmsk_exp
- parent nucmass_mnmsk
# 
# Class nucmass_gen
#
class nucmass_gen
- parent nucmass_table
# 
# Class nucmass_dglg
#
class nucmass_dglg
- parent nucmass_table
# 
# Class nucmass_hfb
#
class nucmass_hfb
- parent nucmass_table
# 
# Class nucmass_hfb_sp
#
class nucmass_hfb_sp
- parent nucmass_table
# 
# Class nucmass_ktuy
#
class nucmass_ktuy
- parent nucmass_table
# 
# Class nucmass_sdnp
#
class nucmass_sdnp
- parent nucmass_table
# 
# Class nucmass_wlw
#
class nucmass_wlw
- parent nucmass_table
# 
# HDF functions
#
function ame_load
- void
- nucmass_ame &ame
- std::string name
- bool exp_only
function ame_load_ext
- void
- nucmass_ame &ame
- std::string file_name
- std::string table_name
- bool exp_only
function mnmsk_load
- void
- nucmass_mnmsk &mnmsk
- std::string model
- std::string filename
function hfb_load
- void
- nucmass_hfb &hfb
- size_t model
- std::string filename
function hfb_sp_load
- void
- nucmass_hfb_sp &hfb
- size_t model
- std::string filename

