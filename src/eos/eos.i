# Interface file for o2scl eos classes
#
namespace o2scl
py_class_doc_pattern "Python interface for class :ref:`%name% <o2sclp:%name%>`."
# 
# Include statements for C++ header file
# 
h_include <o2scl/eos_base.h>
h_include <o2scl/eos_had_base.h>
h_include <o2scl/part.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/eos_python.h>
cpp_include <o2scl/part.h>
# 
# Namespace to use in C++ source code
# 
cpp_using std
cpp_using o2scl
# 
# Class eos_base
#
class eos_base
- o2scl::thermo def_thermo
- function set_thermo
  - void
  - o2scl::thermo &th
- function get_thermo
  - const o2scl::thermo &
# 
# Class eos_had_base
#
class eos_had_base abstract
- parent eos_base
- double eoa
