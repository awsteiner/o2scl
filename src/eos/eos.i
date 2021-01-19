# Interface file for o2scl eos classes
#
namespace o2scl
py_class_doc_pattern "Python interface for class :ref:`%name% <o2scle:%name%>`."
# 
# Include statements for C++ header file
# 
h_include <o2scl/eos_base.h>
h_include <o2scl/eos_had_base.h>
h_include <o2scl/eos_had_skyrme.h>
h_include <o2scl/eos_had_apr.h>
h_include <o2scl/eos_had_rmf.h>
h_include <o2scl/eos_quark_bag.h>
h_include <o2scl/eos_quark_njl.h>
h_include <o2scl/part.h>
h_include <o2scl/fermion_nonrel.h>
h_include <o2scl/fermion_deriv_nr.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/eos_python.h>
cpp_include <o2scl/part.h>
cpp_include <o2scl/fermion_nonrel.h>
cpp_include <o2scl/fermion_deriv_nr.h>
# 
# Namespace to use in C++ source code
# 
cpp_using std
cpp_using o2scl
#
# ------------------------------------------------------
# 
# Class eos_base
#
class eos_base
- o2scl::thermo def_thermo
# 1/17/2021 removing set_thermo and get_thermo for now, as cpppy does
# not yet handle these kind of functions.
# - function set_thermo
#  - void
#  - o2scl::thermo &th
# - function get_thermo
#   - const o2scl::thermo &
#
# ------------------------------------------------------
# 
# Class eos_had_base
#
class eos_had_base abstract
- parent eos_base
- double eoa
- double msom
- double comp
- double n0
- double esym
- double kprime
- bool err_nonconv
- o2scl::fermion def_neutron
- o2scl::fermion def_proton
- function calc_e
  - int
  - o2scl::fermion &n
  - o2scl::fermion &p
  - o2scl::thermo &th
# 
# Class eos_had_eden_base
#
class eos_had_eden_base abstract
- parent eos_had_base
# 
# Class eos_had_pres_base
#
class eos_had_pres_base abstract
- parent eos_had_base
# 
# Class eos_had_temp_base 
#
class eos_had_temp_base abstract
- parent eos_had_base
# 
# Class eos_had_temp_eden_base
#
class eos_had_temp_eden_base abstract
- parent eos_had_temp_base
# 
# Class eos_had_temp_pres_base
#
class eos_had_temp_pres_base abstract
- parent eos_had_temp_base
# 
# Class eos_had_skyrme
#
class eos_had_skyrme
- parent eos_had_temp_eden_base
- double t0
- double t1
- double t2
- double t3
- double x0
- double x1
- double x2
- double x3
- double alpha
- double a
- double b
- double W0
- double b4
- double b4p
- bool parent_method
- std::string reference
- o2scl::fermion_deriv_nr nrfd
# 
# Class eos_had_apr
#
class eos_had_apr
- parent eos_had_temp_eden_base
- int pion
- bool parent_method
# 
# Class eos_had_rmf
#
class eos_had_rmf
- parent eos_had_temp_pres_base
- size_t calc_e_steps
- bool calc_e_relative
- bool zm_mode
- int verbose
- bool err_nonconv
- double mnuc
- double ms
- double mw
- double mr
- double cs
- double cw
- double cr
- double b
- double c
- double zeta
- double xi
- double a1
- double a2
- double a3
- double a4
- double a5
- double a6
- double b1
- double b2
- double b3
# 
# Class eos_quark
#
class eos_quark
- parent eos_base
# 
# Class eos_quark_bag
#
class eos_quark_bag
- parent eos_quark
- double bag_constant
# 
# Class eos_quark_njl
#
class eos_quark_njl
- parent eos_quark
- double B0
- double L
- double G
- double K
- double limit
- bool fromqq
# - o2scl::quark def_up
# - o2scl::quark def_down
# - o2scl::quark def_strange

