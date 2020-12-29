# Interface file for o2scl part classes
# 
namespace o2scl
# 
# Include statements for header file
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
# Class thermo
#
class thermo
- double ed
- double pr
- double en
# 
# Class part
#
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
- function init
  - void
  - double mass
  - double dof
- function anti
  - void
  - part &ax
# 
# Class fermion
#
class fermion
- parent part
- double kf
- double del
# 
# Class fermion_zerot
#
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
# 
# Class fermion_rel
#   
# The class fermion_thermo is between fermion_rel and
# fermion_zerot, but it is abstract, so we just skip over it
# and add it's associated member functions to fermion_rel for now
class fermion_rel
- parent fermion_zerot
- bool err_nonconv
- double min_psi
- double deg_limit
- double exp_limit
- double upper_limit_fac
- double deg_entropy_fac
- int verbose
- bool use_expansions
- double tol_expan
- bool verify_ti
- double therm_ident
- fermion unc
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
# 
# Class fermion_nonrel
#
class fermion_nonrel
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
# Class boson
#
class boson
- parent part
- double co
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
# Class classical_thermo
#
class classical_thermo
- function calc_density
  - void
  - part &p
  - double T
- function calc_mu
  - void
  - part &p
  - double T
# 
# Class thermo_np_deriv_press
#
class thermo_np_deriv_press
- double dsdT
- double dnndT
- double dnpdT
- double dnndmun
- double dndmu_mixed
- double dnpdmup
# 
# Class thermo_np_deriv_helm
#
class thermo_np_deriv_helm
- double dsdT
- double dmundT
- double dmupdT
- double dmundnn
- double dmudn_mixed
- double dmupdnp
# 
# Class part_deriv_press
#
class part_deriv_press
- double dndmu
- double dndT
- double dsdT
- function deriv_f
  - void
  - double &dmudn
  - double &dmudT
  - double &dsdT_n
# 
# Class part_deriv
#
class part_deriv
- parent part
- parent part_deriv_press
# 
# Class fermion_deriv
#
class fermion_deriv
- parent fermion
- parent part_deriv_press
# 
# Class deriv_thermo_base
#
class deriv_thermo_base
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
# Class fermion_deriv_rel
#
class fermion_deriv_rel
- double exp_limit
- double deg_limit
- double upper_limit_fac
- fermion_deriv unc
- fermion_rel fr
- int method
- int last_method
- bool err_nonconv
- function nu_from_n
  - int
  - fermion_deriv &f
  - double T
- function calc_density
  - int
  - fermion_deriv &f
  - double T
- function pair_density
  - int
  - fermion_deriv &f
  - double T
- function calc_mu
  - int
  - fermion_deriv &f
  - double T
- function pair_mu
  - int
  - fermion_deriv &f
  - double T
# 
# Class fermion_deriv_nr
#
class fermion_deriv_nr
- double flimit
- fermion_deriv unc
- function calc_density_zerot
  - void
  - fermion_deriv &f
- function calc_mu_zerot
  - void
  - fermion_deriv &f
- function nu_from_n
  - int
  - fermion_deriv &f
  - double T
- function calc_density
  - int
  - fermion_deriv &f
  - double T
- function calc_mu
  - int
  - fermion_deriv &f
  - double T
# 
# Class classical_deriv_thermo
#
class classical_deriv_thermo
- function calc_density
  - void
  - part_deriv &p
  - double T
- function calc_mu
  - void
  - part_deriv &p
  - double T
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


