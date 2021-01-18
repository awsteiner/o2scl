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


