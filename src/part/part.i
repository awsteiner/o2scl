# Interface file for o2scl part classes
namespace o2scl
# Include statements for header file
h_include <o2scl/part.h>
h_include <o2scl/fermion_rel.h>
h_include <o2scl/fermion_nonrel.h>
h_include <o2scl/fermion_deriv_nr.h>
h_include <o2scl/fermion_deriv_rel.h>
h_include <o2scl/boson_rel.h>
h_include <o2scl/classical.h>
h_include <o2scl/classical_deriv.h>
h_include <o2scl/fermion_mag_zerot.h>
# Include statement for C++ source code
cpp_include <o2scl/part_python.h>
# Namespace to use in C++ source code
cpp_using std
cpp_using o2scl
# 
# Classes
class part
- double g
- double m
- double ms
- double mu
- double nu
- double ed
- double pr
- double en
- bool inc_rest_mass
- bool non_interacting
class fermion
- parent part
- double kf
- double del
class fermion_rel
- function calc_density
  - int
  - fermion &f
  - double T
- function calc_mu
  - void
  - fermion &f
  - double T
class fermion_nonrel
class fermion_deriv_nr
class fermion_deriv_rel
class boson_rel
class classical_thermo
class classical_deriv_thermo
class fermion_mag_zerot



