# Interface file for o2scl part classes
# 
namespace o2scl
py_class_doc |
| Python interface for O2scl class ``%name%``,
| See
| https://awsteiner.org/code/o2scl/html/class/%name%.html .
dll_name o2scl
rst_header |
| .. _particle:
|
| Particle classes
| ===============================
|
| :ref:`O2sclpy <o2sclpy>`
# 
# Include statements for C++ header file
# 
h_include <o2scl/part.h>
h_include <o2scl/fermion_rel.h>
h_include <o2scl/fermion_nonrel.h>
h_include <o2scl/fermion_deriv_nr.h>
h_include <o2scl/fermion_deriv_rel.h>
h_include <o2scl/boson_rel.h>
h_include <o2scl/classical.h>
h_include <o2scl/classical_deriv.h>
h_include <o2scl/fermion_mag_zerot.h>
h_include <o2scl/quark.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/part_python.h>
# 
# Namespace to use in C++ source code
# 
cpp_using std
cpp_using o2scl
#
# ------------------------------------------------------
# 
# Class thermo
#
class thermo
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- double ed
- double pr
- double en
#
# ------------------------------------------------------
# 
# Class part
#
class part
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- double g
- double m
- double ms
- double mu
- double nu
- double n
- double ed
- double pr
- double en
- bool inc_rest_mass
- bool non_interacting
- function init
  - void
  - double mass
  - double dof
- function anti
  - void
  - part &ax
#
# ------------------------------------------------------
# 
# Class fermion
#
class fermion
- py_class_doc |
| Python interface for O2scl class ``%name%``,
| which is a typedef of ``%name%_tl<double>``. See
| https://awsteiner.org/code/o2scl-dev/part/html/class/%name%_tl.html
| .
- parent part
- double kf
- double del
  - py_name delta
#
# ------------------------------------------------------
# 
# Class quark
#
class quark
- parent fermion
- double B
- double qq
#
# ------------------------------------------------------
# 
# Class fermion_zerot
#
class fermion_zerot
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
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
#
# ------------------------------------------------------
# 
# Class fermion_thermo
#
class fermion_thermo abstract
- py_class_doc |
| Python interface for :ref:`%name% <o2scl:%name%_tl>`.
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
#
# ------------------------------------------------------
# 
# Class fermion_rel
#   
class fermion_rel
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- parent fermion_thermo
- bool err_nonconv
- double min_psi
- double deg_limit
- double upper_limit_fac
- int verbose
- bool use_expansions
- double tol_expan
- bool verify_ti
- fermion unc
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
#
# ------------------------------------------------------
# 
# Class fermion_nonrel
#
class fermion_nonrel
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- parent fermion_zerot
- function calc_density
  - int
  - fermion &f
  - double T
- function calc_mu
  - void
  - fermion &f
  - double T
- function nu_from_n
  - void
  - fermion &f
  - double T
#
# ------------------------------------------------------
# 
# Class boson
#
class boson
- parent part
- double co
#
# ------------------------------------------------------
# 
# Class boson_rel
#
class boson_rel
- function calc_density
  - void
  - boson &b
  - double T
- function calc_mu
  - void
  - boson &b
  - double T
- function nu_from_n
  - void
  - boson &b
  - double T
- function pair_density
  - void
  - boson &b
  - double T
- function pair_mu
  - void
  - boson &b
  - double T
#
# ------------------------------------------------------
# 
# Class classical_thermo
#
class classical_thermo
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- function calc_density
  - void
  - part &p
  - double T
- function calc_mu
  - void
  - part &p
  - double T
#
# ------------------------------------------------------
# 
# Class thermo_np_deriv_press
#
class thermo_np_deriv_press
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- double dsdT
- double dnndT
- double dnpdT
- double dnndmun
- double dndmu_mixed
- double dnpdmup
#
# ------------------------------------------------------
# 
# Class thermo_np_deriv_helm
#
class thermo_np_deriv_helm
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- double dsdT
- double dmundT
- double dmupdT
- double dmundnn
- double dmudn_mixed
- double dmupdnp
#
# ------------------------------------------------------
# 
# Class part_deriv_press
#
class part_deriv_press
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- double dndmu
- double dndT
- double dsdT
- function deriv_f
  - void
  - out double &dmudn
  - out double &dmudT
  - out double &dsdT_n
#
# ------------------------------------------------------
# 
# Class part_deriv
#
class part_deriv
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- parent part
- parent part_deriv_press
#
# ------------------------------------------------------
# 
# Class fermion_deriv
#
class fermion_deriv
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- parent fermion
- parent part_deriv_press
#
# ------------------------------------------------------
# 
# Class deriv_thermo_base
#
class deriv_thermo_base
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- function heat_cap_ppart_const_vol
  - double
  - part_deriv &p
  - double T
- function heat_cap_ppart_const_press
  - double
  - part_deriv &p
  - double T
- function compress_adiabatic
  - double
  - part_deriv &p
  - double T
- function compress_const_tptr
  - double
  - part_deriv &p
  - double T
- function coeff_thermal_exp
  - double
  - part_deriv &p
  - double T
- function squared_sound_speed
  - double
  - part_deriv &p
  - double T
#
# ------------------------------------------------------
# 
# Class fermion_deriv_rel
#
# class fermion_deriv_rel
# - py_class_doc |
# | Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
# - double deg_limit
# - double upper_limit_fac
# - fermion_deriv unc
# # We cannot include fr until it has a real copy constructor
# # - fermion_rel fr
# - int method
# - int last_method
# - bool err_nonconv
# - function nu_from_n
#   - int
#   - fermion_deriv &f
#   - double T
# - function calc_density
#   - int
#   - fermion_deriv &f
#   - double T
# - function pair_density
#   - int
#   - fermion_deriv &f
#   - double T
# - function calc_mu
#   - int
#   - fermion_deriv &f
#   - double T
# - function pair_mu
#   - int
#   - fermion_deriv &f
#   - double T
#
# ------------------------------------------------------
# 
# Class fermion_deriv_nr
#
# class fermion_deriv_nr
# - py_class_doc |
# | Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
# - double flimit
# - fermion_deriv unc
# - function calc_density_zerot
#   - void
#   - fermion_deriv &f
# - function calc_mu_zerot
#   - void
#   - fermion_deriv &f
# - function nu_from_n
#   - int
#   - fermion_deriv &f
#   - double T
# - function calc_density
#   - int
#   - fermion_deriv &f
#   - double T
# - function calc_mu
#   - int
#   - fermion_deriv &f
#   - double T
#
# ------------------------------------------------------
# 
# Class classical_deriv_thermo
#
class classical_deriv_thermo
- py_class_doc |
| Python interface for class :ref:`%name% <o2scl:%name%_tl>`.
- function calc_density
  - void
  - part_deriv &p
  - double T
- function calc_mu
  - void
  - part_deriv &p
  - double T
#
# ------------------------------------------------------
# 
# Class fermion_mag_zerot
#
class fermion_mag_zerot
- int nmax_up
- int nmax_dn
- int sum_limit
- function calc_mu_zerot_mag
  - void
  - fermion &f
  - double qB
  - double kappa
- function calc_density_zerot_mag
  - void
  - fermion &f
  - double qB
  - double kappa
