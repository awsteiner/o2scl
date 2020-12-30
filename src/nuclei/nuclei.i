# Interface file for o2scl nuclei classes
# 
namespace o2scl
# 
# Include statements for C++ header file
# 
h_include <o2scl/nucleus.h>
h_include <o2scl/nucmass.h>
# 
# Include statement for C++ source code
# 
cpp_include <o2scl/nuclei_python.h>
# 
# Namespace to use in C++ source code
# 
cpp_using std
cpp_using o2scl
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
  - int &Z
  - int &N
  - int &A
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
- function fit_fun
  - int
  - const ubvector &x
- function guess_fun
  - int
  - ubvector &x
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
