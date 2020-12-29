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
class fermion_zerot
- function kf_from_density
  - void
  - fermion &f
- function energy_density_zerot
  - void
  - fermion &f
- function pressure_zerot
  - void
  - fermion &f
- function calc_mu_zerot
  - void
  - fermion &f
- function calc_density_zerot
  - void
  - fermion &f
# The class fermion_thermo is between fermion_rel and
# fermion_zerot, but it is abstract, so we just skip over it
# and add it's associated member functions to fermion_rel for now
class fermion_rel
- parent fermion_zerot
- function calc_mu_deg
  - bool
  - fermion &f
  - double T
  - double prec
- function calc_mu_ndeg
  - bool
  - fermion &f
  - double T
  - double prec
  - bool inc_antip
- function massless_calc_mu
  - void
  - fermion &f
  - double T
- function massless_pair_mu
  - void
  - fermion &f
  - double T
- function massless_calc_density
  - void
  - fermion &f
  - double T
- function massless_pair_density
  - void
  - fermion &f
  - double T
- function nu_from_n
  - int
  - fermion &f
  - double T
- function calc_density
  - int
  - fermion &f
  - double T
- function pair_density
  - int
  - fermion &f
  - double T
- function calc_mu
  - void
  - fermion &f
  - double T
- function pair_mu
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



